use std::collections::HashMap;

use spade::handles::FixedVertexHandle;
use spade::{
    AngleLimit, ConstrainedDelaunayTriangulation, Point2, RefinementParameters, Triangulation,
};

use crate::error::MeshError;
use crate::geometry::{EdgeRef, Panel, SharedEdgeDeclaration};
use crate::mesher::r#trait::Mesher;
use crate::model::{MeshConfig, MeshModel, MeshQuality, TriElement6};

// ── Key helpers ──────────────────────────────────────────────────────────

/// Lossless key for deduplicating float coordinates in a HashMap.
#[inline]
fn vertex_key(p: [f64; 2]) -> (u64, u64) {
    (p[0].to_bits(), p[1].to_bits())
}

fn midpoint(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5]
}

fn dist2(a: [f64; 2], b: [f64; 2]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    dx * dx + dy * dy
}

/// Signed area of a triangle (positive = CCW).
fn signed_area(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])
}

// ── Point-in-polygon (winding number) ───────────────────────────────────

/// Returns true if point `p` is strictly inside the polygon (not in any hole).
fn point_in_polygon(p: [f64; 2], outline: &[[f64; 2]], holes: &[Vec<[f64; 2]>]) -> bool {
    if !winding_inside(p, outline) {
        return false;
    }
    for hole in holes {
        if winding_inside(p, hole) {
            return false;
        }
    }
    true
}

/// Winding-number test: is `p` inside the polygon defined by `verts`?
fn winding_inside(p: [f64; 2], verts: &[[f64; 2]]) -> bool {
    let n = verts.len();
    let mut winding = 0i32;
    for i in 0..n {
        let a = verts[i];
        let b = verts[(i + 1) % n];
        if a[1] <= p[1] {
            if b[1] > p[1] && signed_area(a, b, p) > 0.0 {
                winding += 1;
            }
        } else if b[1] <= p[1] && signed_area(a, b, p) < 0.0 {
            winding -= 1;
        }
    }
    winding != 0
}

// ── SpadeMesher ──────────────────────────────────────────────────────────

pub struct SpadeMesher;

impl Mesher for SpadeMesher {
    fn triangulate(
        &self,
        panels: &[Panel],
        shared_edges: &[SharedEdgeDeclaration],
        config: &MeshConfig,
    ) -> Result<MeshModel, MeshError> {
        // ── 1. Collect and deduplicate vertices ──────────────────────────
        //
        // We build a global node list and a map from bit-cast coordinate pair
        // to global node ID. All subsequent operations reference node IDs.

        let mut nodes: Vec<[f64; 2]> = Vec::new();
        let mut coord_to_id: HashMap<(u64, u64), usize> = HashMap::new();

        let mut get_or_insert = |coord: [f64; 2]| -> usize {
            let key = vertex_key(coord);
            if let Some(&id) = coord_to_id.get(&key) {
                return id;
            }
            let id = nodes.len();
            nodes.push(coord);
            coord_to_id.insert(key, id);
            id
        };

        // Insert all outline and hole vertices in panel order so shared-edge
        // vertices get the same ID regardless of which panel inserts them first.
        for panel in panels {
            for &v in &panel.outline.vertices {
                get_or_insert(v);
            }
            for hole in &panel.holes {
                for &v in &hole.vertices {
                    get_or_insert(v);
                }
            }
        }

        // ── 2. Validate shared edges ────────────────────────────────────

        let bbox_diag = bounding_box_diagonal(&nodes);
        let tol = 1e-9_f64 * bbox_diag.max(1.0);

        for decl in shared_edges {
            let pa = panels
                .iter()
                .find(|p| p.id == decl.panel_a)
                .ok_or_else(|| MeshError::UnknownPanel(decl.panel_a.clone()))?;
            let pb = panels
                .iter()
                .find(|p| p.id == decl.panel_b)
                .ok_or_else(|| MeshError::UnknownPanel(decl.panel_b.clone()))?;

            let (a0, a1) = edge_endpoints(pa, decl.edge_a);
            let (b0, b1) = edge_endpoints(pb, decl.edge_b);

            // The two edges must coincide with opposite orientation (a0≈b1, a1≈b0).
            if dist2(a0, b1) > tol * tol || dist2(a1, b0) > tol * tol {
                return Err(MeshError::SharedEdgeMismatch {
                    panel_a: decl.panel_a.clone(),
                    edge_a: decl.edge_a,
                    panel_b: decl.panel_b.clone(),
                    edge_b: decl.edge_b,
                });
            }
        }

        // ── 3. Build CDT and insert constrained edges ───────────────────

        let mut cdt = ConstrainedDelaunayTriangulation::<Point2<f64>>::new();

        // Keep a map from bit-cast coord → Spade FixedVertexHandle.
        let mut spade_handles: HashMap<(u64, u64), FixedVertexHandle> = HashMap::new();

        for &coord in &nodes {
            let handle = cdt
                .insert(Point2::new(coord[0], coord[1]))
                .map_err(|e| MeshError::TriangulationFailed(e.to_string()))?;
            spade_handles.insert(vertex_key(coord), handle);
        }

        // Add constrained edges for every polygon boundary segment.
        for panel in panels {
            insert_polygon_constraints(&mut cdt, &panel.outline.vertices, &spade_handles);
            for hole in &panel.holes {
                insert_polygon_constraints(&mut cdt, &hole.vertices, &spade_handles);
            }
        }

        // ── 4. Ruppert refinement ────────────────────────────────────────

        // Note: exclude_outer_faces is intentionally omitted. For non-convex panels
        // Spade's flood-fill misidentifies concave interior regions as "outer", which
        // suppresses refinement there. Instead we refine the full convex hull and rely
        // on the centroid point-in-polygon test (step 5) to discard exterior triangles.
        //
        // Spade's default vertex cap is initial_vertices × 10, which is far too small
        // for fine meshes. We estimate an upper bound from the bounding-box area and
        // the requested element area, then add a generous safety factor of 4×.
        let bbox_area = {
            let (mut xmin, mut xmax, mut ymin, mut ymax) =
                (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY);
            for &[x, y] in &nodes {
                xmin = xmin.min(x);
                xmax = xmax.max(x);
                ymin = ymin.min(y);
                ymax = ymax.max(y);
            }
            (xmax - xmin).max(0.0) * (ymax - ymin).max(0.0)
        };
        let vertex_budget = ((bbox_area / config.max_element_area) * 4.0).ceil() as usize;
        let vertex_budget = vertex_budget.max(1_000);

        let params = RefinementParameters::<f64>::new()
            .with_max_allowed_area(config.max_element_area)
            .with_angle_limit(AngleLimit::from_deg(config.min_angle_deg))
            .with_max_additional_vertices(vertex_budget);

        cdt.refine(params);

        // ── 5. Extract 3-node triangles ──────────────────────────────────

        // After refinement new Steiner points were inserted into `cdt`. We need
        // to re-sync our `nodes` / `coord_to_id` maps with these new vertices.
        for vertex in cdt.vertices() {
            let p = vertex.position();
            let coord = [p.x, p.y];
            let key = vertex_key(coord);
            coord_to_id.entry(key).or_insert_with(|| {
                let id = nodes.len();
                nodes.push(coord);
                id
            });
        }

        // Collect panel outlines and holes for point-in-polygon queries.
        struct PanelGeom {
            id: String,
            outline: Vec<[f64; 2]>,
            holes: Vec<Vec<[f64; 2]>>,
        }
        #[allow(clippy::type_complexity)]
        let panel_outlines: Vec<PanelGeom> = panels
            .iter()
            .map(|p| PanelGeom {
                id: p.id.clone(),
                outline: p.outline.vertices.clone(),
                holes: p.holes.iter().map(|h| h.vertices.clone()).collect(),
            })
            .collect();

        let mut elements: Vec<TriElement6> = Vec::new();
        // Midpoint deduplication: same as vertices.
        let mut mid_id: HashMap<(u64, u64), usize> = coord_to_id.clone();

        let inner_faces: Vec<_> = cdt.inner_faces().collect();

        for face in inner_faces {
            let verts = face.vertices();
            let positions: [[f64; 2]; 3] = [
                [verts[0].position().x, verts[0].position().y],
                [verts[1].position().x, verts[1].position().y],
                [verts[2].position().x, verts[2].position().y],
            ];

            // Ensure CCW orientation (swap if needed).
            let [p0, p1, p2] = if signed_area(positions[0], positions[1], positions[2]) >= 0.0 {
                positions
            } else {
                [positions[0], positions[2], positions[1]]
            };

            let corners = [
                *coord_to_id.get(&vertex_key(p0)).unwrap(),
                *coord_to_id.get(&vertex_key(p1)).unwrap(),
                *coord_to_id.get(&vertex_key(p2)).unwrap(),
            ];

            // Assign panel_id via centroid test.
            let centroid = [(p0[0] + p1[0] + p2[0]) / 3.0, (p0[1] + p1[1] + p2[1]) / 3.0];
            let panel_id = panel_outlines
                .iter()
                .find(|pg| point_in_polygon(centroid, &pg.outline, &pg.holes))
                .map(|pg| pg.id.clone())
                .unwrap_or_else(|| "_outside".into());

            // Skip triangles that belong to no panel (outer region).
            if panel_id == "_outside" {
                continue;
            }

            // Compute midpoint node IDs.
            let corner_pts = [p0, p1, p2];
            let mut midpoints = [0usize; 3];
            for i in 0..3 {
                let m = midpoint(corner_pts[i], corner_pts[(i + 1) % 3]);
                let key = vertex_key(m);
                let mid = if let Some(&id) = mid_id.get(&key) {
                    id
                } else {
                    let id = nodes.len();
                    nodes.push(m);
                    mid_id.insert(key, id);
                    // Also update coord_to_id for the edge_nodes step later.
                    coord_to_id.insert(key, id);
                    id
                };
                midpoints[i] = mid;
            }

            elements.push(TriElement6 {
                id: elements.len(),
                corners,
                midpoints,
                panel_id,
            });
        }

        if elements.is_empty() {
            return Err(MeshError::TriangulationFailed(
                "No interior triangles produced".into(),
            ));
        }

        // ── 6. Build edge_nodes map ──────────────────────────────────────

        let mut edge_nodes: HashMap<EdgeRef, Vec<usize>> = HashMap::new();

        for panel in panels {
            let n = panel.outline.n_edges();
            let mut all_edge_ids: Vec<usize> = Vec::new();

            for i in 0..n {
                let v_start = panel.outline.vertices[i];
                let v_end = panel.outline.vertices[(i + 1) % n];
                let ordered = nodes_on_segment(&nodes, &coord_to_id, v_start, v_end);
                all_edge_ids.extend_from_slice(&ordered);
                edge_nodes.insert(
                    EdgeRef::PolygonEdge {
                        panel_id: panel.id.clone(),
                        edge_index: i,
                    },
                    ordered,
                );
            }
            // AllEdges: all outline edge node IDs in edge-index order, deduplicated.
            edge_nodes.insert(
                EdgeRef::AllEdges {
                    panel_id: panel.id.clone(),
                },
                dedup_preserve_order(all_edge_ids),
            );

            for (hole_idx, hole) in panel.holes.iter().enumerate() {
                let hn = hole.n_edges();
                let mut hole_all: Vec<usize> = Vec::new();
                for i in 0..hn {
                    let v_start = hole.vertices[i];
                    let v_end = hole.vertices[(i + 1) % hn];
                    let ordered = nodes_on_segment(&nodes, &coord_to_id, v_start, v_end);
                    hole_all.extend_from_slice(&ordered);
                    edge_nodes.insert(
                        EdgeRef::HoleEdge {
                            panel_id: panel.id.clone(),
                            hole_index: hole_idx,
                            edge_index: i,
                        },
                        ordered,
                    );
                }
                edge_nodes.insert(
                    EdgeRef::AllHoleEdges {
                        panel_id: panel.id.clone(),
                        hole_index: hole_idx,
                    },
                    dedup_preserve_order(hole_all),
                );
            }
        }

        // ── 7. Quality metrics ───────────────────────────────────────────

        let quality = compute_quality(&nodes, &elements);

        Ok(MeshModel {
            nodes,
            elements,
            edge_nodes,
            quality,
        })
    }
}

// ── Internal helpers ─────────────────────────────────────────────────────

fn edge_endpoints(panel: &Panel, edge_index: usize) -> ([f64; 2], [f64; 2]) {
    let n = panel.outline.vertices.len();
    (
        panel.outline.vertices[edge_index % n],
        panel.outline.vertices[(edge_index + 1) % n],
    )
}

fn bounding_box_diagonal(nodes: &[[f64; 2]]) -> f64 {
    if nodes.is_empty() {
        return 1.0;
    }
    let (mut xmin, mut xmax, mut ymin, mut ymax) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for &[x, y] in nodes {
        xmin = xmin.min(x);
        xmax = xmax.max(x);
        ymin = ymin.min(y);
        ymax = ymax.max(y);
    }
    let dx = xmax - xmin;
    let dy = ymax - ymin;
    (dx * dx + dy * dy).sqrt()
}

fn insert_polygon_constraints(
    cdt: &mut ConstrainedDelaunayTriangulation<Point2<f64>>,
    verts: &[[f64; 2]],
    handles: &HashMap<(u64, u64), FixedVertexHandle>,
) {
    let n = verts.len();
    for i in 0..n {
        let a = verts[i];
        let b = verts[(i + 1) % n];
        if let (Some(&ha), Some(&hb)) = (handles.get(&vertex_key(a)), handles.get(&vertex_key(b))) {
            cdt.add_constraint(ha, hb);
        }
    }
}

/// Collect all nodes that lie on the line segment from `start` to `end`,
/// sorted by distance from `start`.
fn nodes_on_segment(
    nodes: &[[f64; 2]],
    coord_to_id: &HashMap<(u64, u64), usize>,
    start: [f64; 2],
    end: [f64; 2],
) -> Vec<usize> {
    let dx = end[0] - start[0];
    let dy = end[1] - start[1];
    let seg_len2 = dx * dx + dy * dy;

    let mut on_seg: Vec<(f64, usize)> = Vec::new();

    for (id, &coord) in nodes.iter().enumerate() {
        // Parametric t along the segment.
        let t = if seg_len2 < 1e-30 {
            0.0
        } else {
            ((coord[0] - start[0]) * dx + (coord[1] - start[1]) * dy) / seg_len2
        };
        if !(-1e-9..=1.0 + 1e-9).contains(&t) {
            continue;
        }
        // Check that the node is actually on the line (perpendicular distance).
        let proj = [start[0] + t * dx, start[1] + t * dy];
        if dist2(coord, proj) < 1e-18 * seg_len2.max(1e-18) {
            on_seg.push((t, id));
        }
    }

    // Fallback: ensure start and end are always included.
    let start_id = *coord_to_id.get(&vertex_key(start)).unwrap();
    let end_id = *coord_to_id.get(&vertex_key(end)).unwrap();
    if !on_seg.iter().any(|(_, id)| *id == start_id) {
        on_seg.push((0.0, start_id));
    }
    if !on_seg.iter().any(|(_, id)| *id == end_id) {
        on_seg.push((1.0, end_id));
    }

    on_seg.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    on_seg.into_iter().map(|(_, id)| id).collect()
}

fn dedup_preserve_order(ids: Vec<usize>) -> Vec<usize> {
    let mut seen = std::collections::HashSet::new();
    ids.into_iter().filter(|id| seen.insert(*id)).collect()
}

fn compute_quality(nodes: &[[f64; 2]], elements: &[TriElement6]) -> MeshQuality {
    let mut min_angle = f64::INFINITY;
    let mut max_aspect = 0.0_f64;

    for el in elements {
        let [a, b, c] = [
            nodes[el.corners[0]],
            nodes[el.corners[1]],
            nodes[el.corners[2]],
        ];

        let lab = dist2(a, b).sqrt();
        let lbc = dist2(b, c).sqrt();
        let lca = dist2(c, a).sqrt();

        if lab < 1e-15 || lbc < 1e-15 || lca < 1e-15 {
            continue;
        }

        // Min angle via law of cosines.
        let cos_a = (lab * lab + lca * lca - lbc * lbc) / (2.0 * lab * lca);
        let cos_b = (lab * lab + lbc * lbc - lca * lca) / (2.0 * lab * lbc);
        let cos_c = (lbc * lbc + lca * lca - lab * lab) / (2.0 * lbc * lca);

        for &cos in &[cos_a, cos_b, cos_c] {
            let angle_deg = cos.clamp(-1.0, 1.0).acos().to_degrees();
            if angle_deg < min_angle {
                min_angle = angle_deg;
            }
        }

        // Aspect ratio: circumradius / (2 * inradius).
        let area = (signed_area(a, b, c) / 2.0).abs();
        let s = (lab + lbc + lca) / 2.0;
        let inradius = area / s;
        let circumradius = (lab * lbc * lca) / (4.0 * area);
        let aspect = circumradius / (2.0 * inradius);
        if aspect > max_aspect {
            max_aspect = aspect;
        }
    }

    MeshQuality {
        n_elements: elements.len(),
        n_nodes: nodes.len(),
        min_angle_deg: if min_angle.is_finite() {
            min_angle
        } else {
            0.0
        },
        max_aspect_ratio: max_aspect,
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{Panel, Polygon2D};

    fn square_panel(size: f64) -> Panel {
        Panel::new("p", Polygon2D::rectangle(0.0, 0.0, size, size))
    }

    fn default_config() -> MeshConfig {
        MeshConfig {
            max_element_area: 0.1,
            ..MeshConfig::default()
        }
    }

    #[test]
    fn single_panel_produces_elements() {
        let panels = vec![square_panel(2.0)];
        let mesh = SpadeMesher
            .triangulate(&panels, &[], &default_config())
            .unwrap();
        assert!(
            !mesh.elements.is_empty(),
            "should produce at least one element"
        );
    }

    #[test]
    fn elements_are_ccw() {
        let panels = vec![square_panel(2.0)];
        let mesh = SpadeMesher
            .triangulate(&panels, &[], &default_config())
            .unwrap();
        for el in &mesh.elements {
            let a = mesh.nodes[el.corners[0]];
            let b = mesh.nodes[el.corners[1]];
            let c = mesh.nodes[el.corners[2]];
            assert!(
                signed_area(a, b, c) > 0.0,
                "element {} has CW or degenerate corners",
                el.id
            );
        }
    }

    #[test]
    fn midpoint_is_midpoint_of_side() {
        let panels = vec![square_panel(2.0)];
        let mesh = SpadeMesher
            .triangulate(&panels, &[], &default_config())
            .unwrap();
        for el in &mesh.elements {
            for i in 0..3 {
                let a = mesh.nodes[el.corners[i]];
                let b = mesh.nodes[el.corners[(i + 1) % 3]];
                let m = mesh.nodes[el.midpoints[i]];
                let expected = midpoint(a, b);
                assert!(
                    dist2(m, expected) < 1e-20,
                    "midpoint[{i}] of element {} is wrong",
                    el.id
                );
            }
        }
    }

    #[test]
    fn edge_nodes_contain_endpoints() {
        let panels = vec![square_panel(2.0)];
        let mesh = SpadeMesher
            .triangulate(&panels, &[], &default_config())
            .unwrap();
        let key = EdgeRef::PolygonEdge {
            panel_id: "p".into(),
            edge_index: 0,
        };
        let ids = mesh.edge_nodes.get(&key).expect("edge 0 missing");
        assert!(!ids.is_empty());
        // First and last node should be at the corners of edge 0: (0,0) → (2,0).
        let first = mesh.nodes[*ids.first().unwrap()];
        let last = mesh.nodes[*ids.last().unwrap()];
        assert!(dist2(first, [0.0, 0.0]) < 1e-20 || dist2(first, [2.0, 0.0]) < 1e-20);
        assert!(dist2(last, [0.0, 0.0]) < 1e-20 || dist2(last, [2.0, 0.0]) < 1e-20);
    }

    #[test]
    fn two_panel_shared_edge_uses_same_node_ids() {
        // span_1: (0,0)→(2,0)→(2,2)→(0,2)  edge 1 = (2,0)→(2,2)
        // span_2: (2,0)→(4,0)→(4,2)→(2,2)  edge 3 = (2,2)→(2,0)  (opposite)
        let panels = vec![
            Panel::new(
                "s1",
                Polygon2D::new(vec![[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]]),
            ),
            Panel::new(
                "s2",
                Polygon2D::new(vec![[2.0, 0.0], [4.0, 0.0], [4.0, 2.0], [2.0, 2.0]]),
            ),
        ];
        let shared = vec![SharedEdgeDeclaration::between("s1", 1, "s2", 3)];
        let mesh = SpadeMesher
            .triangulate(&panels, &shared, &default_config())
            .unwrap();

        // Nodes on the shared edge from s1 perspective and s2 perspective should
        // be the same node IDs (same coordinates → same ID via deduplication).
        let s1_ids = mesh
            .edge_nodes
            .get(&EdgeRef::PolygonEdge {
                panel_id: "s1".into(),
                edge_index: 1,
            })
            .expect("s1 edge 1 missing");
        let s2_ids = mesh
            .edge_nodes
            .get(&EdgeRef::PolygonEdge {
                panel_id: "s2".into(),
                edge_index: 3,
            })
            .expect("s2 edge 3 missing");

        let mut s1_sorted = s1_ids.clone();
        let mut s2_sorted = s2_ids.clone();
        s1_sorted.sort();
        s2_sorted.sort();
        assert_eq!(
            s1_sorted, s2_sorted,
            "shared edge should have identical node IDs"
        );
    }
}
