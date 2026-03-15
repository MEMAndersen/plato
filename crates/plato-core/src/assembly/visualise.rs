//! SVG visualisation of the assembled DOF map — compiled only during tests.
//!
//! # Usage
//!
//! ```sh
//! cargo test -p plato-core dump_dof_map_svg -- --ignored --nocapture
//! ```
//!
//! Then open `target/dof_map_debug.svg` in a browser.
//!
//! **Visual encoding:**
//! - Blue **X**  — free w-DOF (displacement); label = global row index
//! - Gray **X**  — pinned w-DOF; label = "pin"
//! - Red  **O**  — active θ-DOF (rotation); label = global row index
//! - Gray **O**  — clamped θ-DOF; label = "clamp"
//!
//! θ markers are placed at 1/3 and 2/3 along each element side.
//! Labels are offset toward the element centroid so they sit inside the triangle.

use plato_mesh::model::MeshModel;

use super::DofMap;

// ── SVG canvas constants ───────────────────────────────────────────────────

const SVG_W: f64 = 900.0;
const SVG_H: f64 = 900.0;
const MARGIN: f64 = 80.0;

// ── Coordinate transform ───────────────────────────────────────────────────

/// World → SVG-pixel transform (uniform scale, Y-flipped).
struct Transform {
    x_min: f64,
    y_min: f64,
    scale: f64,
    x_off: f64,
    y_off: f64,
}

impl Transform {
    fn new(nodes: &[[f64; 2]]) -> Self {
        let inner_w = SVG_W - 2.0 * MARGIN;
        let inner_h = SVG_H - 2.0 * MARGIN;

        let x_min = nodes.iter().map(|n| n[0]).fold(f64::MAX, f64::min);
        let x_max = nodes.iter().map(|n| n[0]).fold(f64::MIN, f64::max);
        let y_min = nodes.iter().map(|n| n[1]).fold(f64::MAX, f64::min);
        let y_max = nodes.iter().map(|n| n[1]).fold(f64::MIN, f64::max);

        let wx = if x_max > x_min { x_max - x_min } else { 1.0 };
        let wy = if y_max > y_min { y_max - y_min } else { 1.0 };

        // Uniform scale — preserve aspect ratio, centre inside inner box.
        let scale = (inner_w / wx).min(inner_h / wy);
        let x_off = MARGIN + (inner_w - scale * wx) * 0.5;
        let y_off = MARGIN + (inner_h - scale * wy) * 0.5;

        Self {
            x_min,
            y_min,
            scale,
            x_off,
            y_off,
        }
    }

    /// Map world (x,y) → SVG pixel (sx, sy).  Y is flipped (SVG Y grows downward).
    #[inline]
    fn px(&self, x: f64, y: f64) -> (f64, f64) {
        let sx = self.x_off + (x - self.x_min) * self.scale;
        let sy = SVG_H - self.y_off - (y - self.y_min) * self.scale;
        (sx, sy)
    }
}

// ── Low-level SVG helpers ──────────────────────────────────────────────────

fn draw_x(svg: &mut String, cx: f64, cy: f64, r: f64, color: &str) {
    let d = r / std::f64::consts::SQRT_2;
    svg.push_str(&format!(
        "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" \
         stroke=\"{color}\" stroke-width=\"1.8\" stroke-linecap=\"round\"/>\n",
        cx - d,
        cy - d,
        cx + d,
        cy + d,
    ));
    svg.push_str(&format!(
        "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" \
         stroke=\"{color}\" stroke-width=\"1.8\" stroke-linecap=\"round\"/>\n",
        cx - d,
        cy + d,
        cx + d,
        cy - d,
    ));
}

fn draw_circle(svg: &mut String, cx: f64, cy: f64, r: f64, color: &str) {
    svg.push_str(&format!(
        "<circle cx=\"{:.2}\" cy=\"{:.2}\" r=\"{:.2}\" \
         fill=\"none\" stroke=\"{color}\" stroke-width=\"1.8\"/>\n",
        cx, cy, r,
    ));
}

fn draw_label(svg: &mut String, x: f64, y: f64, label: &str, color: &str) {
    svg.push_str(&format!(
        "<text x=\"{:.2}\" y=\"{:.2}\" font-size=\"8\" font-family=\"monospace\" \
         text-anchor=\"middle\" dominant-baseline=\"middle\" fill=\"{color}\">{label}</text>\n",
        x, y,
    ));
}

/// Offset `from` by `dist` pixels in the direction of `toward`.
fn offset_toward(from: (f64, f64), toward: (f64, f64), dist: f64) -> (f64, f64) {
    let dx = toward.0 - from.0;
    let dy = toward.1 - from.1;
    let len = (dx * dx + dy * dy).sqrt();
    if len < 1e-6 {
        return (from.0, from.1 - dist); // degenerate: move up
    }
    (from.0 + dx / len * dist, from.1 + dy / len * dist)
}

fn lerp(a: [f64; 2], b: [f64; 2], t: f64) -> [f64; 2] {
    [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]
}

/// Filled triangle arrowhead.  `(ux, uy)` is the unit direction the tip points toward.
fn draw_arrowhead(svg: &mut String, tip_x: f64, tip_y: f64, ux: f64, uy: f64, size: f64, color: &str) {
    let base_x = tip_x - ux * size;
    let base_y = tip_y - uy * size;
    let hw = size * 0.5; // half-width of base
    let lx = base_x - uy * hw;
    let ly = base_y + ux * hw;
    let rx = base_x + uy * hw;
    let ry = base_y - ux * hw;
    svg.push_str(&format!(
        "<polygon points=\"{:.2},{:.2} {:.2},{:.2} {:.2},{:.2}\" \
         fill=\"{color}\" opacity=\"0.85\"/>\n",
        tip_x, tip_y, lx, ly, rx, ry,
    ));
}

// ── Main render function ───────────────────────────────────────────────────

/// Render the mesh DOF layout to a self-contained SVG string.
pub(super) fn render_dof_map_svg(mesh: &MeshModel, dof_map: &DofMap) -> String {
    let tr = Transform::new(&mesh.nodes);
    let mut svg = String::with_capacity(128 * 1024);

    // Header
    svg.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" \
         width=\"{w}\" height=\"{h}\" viewBox=\"0 0 {w} {h}\">\n",
        w = SVG_W as u32,
        h = SVG_H as u32,
    ));
    svg.push_str("<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n");

    // ── 1. Triangle edges ─────────────────────────────────────────────────
    for elem in &mesh.elements {
        let c: Vec<(f64, f64)> = elem
            .corners
            .iter()
            .map(|&id| {
                let [x, y] = mesh.nodes[id];
                tr.px(x, y)
            })
            .collect();
        for (i, j) in [(0usize, 1usize), (1, 2), (2, 0)] {
            svg.push_str(&format!(
                "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" \
                 stroke=\"#ccc\" stroke-width=\"0.8\"/>\n",
                c[i].0, c[i].1, c[j].0, c[j].1,
            ));
        }
        // Element ID at centroid (small gray)
        let cx: f64 = elem
            .corners
            .iter()
            .map(|&id| mesh.nodes[id][0])
            .sum::<f64>()
            / 3.0;
        let cy: f64 = elem
            .corners
            .iter()
            .map(|&id| mesh.nodes[id][1])
            .sum::<f64>()
            / 3.0;
        let (sx, sy) = tr.px(cx, cy);
        svg.push_str(&format!(
            "<text x=\"{:.2}\" y=\"{:.2}\" font-size=\"7\" font-family=\"monospace\" \
             text-anchor=\"middle\" dominant-baseline=\"middle\" fill=\"#aaa\">e{}</text>\n",
            sx, sy, elem.id,
        ));
    }

    // ── 1b. Element orientation arrows ───────────────────────────────────────
    // For each edge (corner[i] → corner[(i+1)%3]), draw a filled arrowhead at
    // the edge midpoint, offset toward the element centroid.
    // Green = CCW (positive signed area in world coords).
    // Orange = CW (negative signed area).
    for elem in &mesh.elements {
        let corners: [[f64; 2]; 3] =
            std::array::from_fn(|i| mesh.nodes[elem.corners[i]]);
        let cx = (corners[0][0] + corners[1][0] + corners[2][0]) / 3.0;
        let cy = (corners[0][1] + corners[1][1] + corners[2][1]) / 3.0;
        let centroid_px = tr.px(cx, cy);

        let signed_area = 0.5
            * ((corners[1][0] - corners[0][0]) * (corners[2][1] - corners[0][1])
                - (corners[2][0] - corners[0][0]) * (corners[1][1] - corners[0][1]));
        let arrow_color = if signed_area > 0.0 { "#27ae60" } else { "#e67e22" };

        for i in 0..3usize {
            let a = corners[i];
            let b = corners[(i + 1) % 3];

            // Edge midpoint in screen space.
            let mid = [(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5];
            let (mx, my) = tr.px(mid[0], mid[1]);

            // Direction a→b in screen space (unit vector).
            let (ax, ay) = tr.px(a[0], a[1]);
            let (bx, by) = tr.px(b[0], b[1]);
            let dx = bx - ax;
            let dy = by - ay;
            let len = (dx * dx + dy * dy).sqrt();
            if len < 1e-6 {
                continue;
            }
            let (ux, uy) = (dx / len, dy / len);

            // Shift midpoint toward centroid before placing the arrowhead.
            let (tip_x, tip_y) = offset_toward((mx, my), centroid_px, 10.0);
            draw_arrowhead(&mut svg, tip_x, tip_y, ux, uy, 7.0, arrow_color);
        }
    }

    // ── 2. w-DOF markers (X) ──────────────────────────────────────────────
    // One marker per unique node ID; offset label toward the centroid of the
    // first element that claims this node.
    let mut seen_w = std::collections::HashSet::<usize>::new();
    for elem in &mesh.elements {
        let cx: f64 = elem
            .corners
            .iter()
            .map(|&id| mesh.nodes[id][0])
            .sum::<f64>()
            / 3.0;
        let cy: f64 = elem
            .corners
            .iter()
            .map(|&id| mesh.nodes[id][1])
            .sum::<f64>()
            / 3.0;
        let centroid_px = tr.px(cx, cy);

        let all_nodes: Vec<usize> = elem
            .corners
            .iter()
            .chain(elem.midpoints.iter())
            .copied()
            .collect();

        for node_id in all_nodes {
            if !seen_w.insert(node_id) {
                continue;
            }
            let [x, y] = mesh.nodes[node_id];
            let (sx, sy) = tr.px(x, y);

            let (color, label): (&str, String) = match dof_map.w_row(node_id) {
                Some(row) => ("#1a73e8", row.to_string()),
                None => ("#999", "pin".to_string()),
            };

            draw_x(&mut svg, sx, sy, 5.0, color);
            let (lx, ly) = offset_toward((sx, sy), centroid_px, 13.0);
            draw_label(&mut svg, lx, ly, &label, color);
        }
    }

    // ── 3. θ-DOF markers (O) ──────────────────────────────────────────────
    // Placed at 1/3 and 2/3 along each side.
    // The circle is deduplicated by screen position (draw once per physical
    // point), but the label is drawn once per element — offset toward that
    // element's centroid.  Interior edges therefore show two labels (one from
    // each adjacent element) on opposite sides of the shared circle.
    let mut seen_circle = std::collections::HashSet::<(i32, i32)>::new();

    // First pass: draw circles (once per physical position).
    for elem in &mesh.elements {
        let corners: [[f64; 2]; 3] = std::array::from_fn(|i| mesh.nodes[elem.corners[i]]);
        for side in 0..3usize {
            let c_a = corners[(side + 1) % 3];
            let c_b = corners[(side + 2) % 3];
            for sub in 0..2usize {
                let t = if sub == 0 { 1.0 / 3.0 } else { 2.0 / 3.0 };
                let pos = lerp(c_a, c_b, t);
                let (sx, sy) = tr.px(pos[0], pos[1]);
                let key = (sx.round() as i32, sy.round() as i32);
                if !seen_circle.insert(key) {
                    continue;
                }
                // Colour by the first element that claims this position.
                let e = elem.id;
                let color = match dof_map.theta_row_opt(e, side, sub) {
                    Some(_) => "#c0392b",
                    None => "#bbb",
                };
                draw_circle(&mut svg, sx, sy, 4.5, color);
            }
        }
    }

    // Second pass: draw one label per (element, side, sub) — each offset
    // toward its own element's centroid.
    for elem in &mesh.elements {
        let e = elem.id;
        let corners: [[f64; 2]; 3] = std::array::from_fn(|i| mesh.nodes[elem.corners[i]]);
        let cx = (corners[0][0] + corners[1][0] + corners[2][0]) / 3.0;
        let cy = (corners[0][1] + corners[1][1] + corners[2][1]) / 3.0;
        let centroid_px = tr.px(cx, cy);

        for side in 0..3usize {
            let c_a = corners[(side + 1) % 3];
            let c_b = corners[(side + 2) % 3];
            for sub in 0..2usize {
                let t = if sub == 0 { 1.0 / 3.0 } else { 2.0 / 3.0 };
                let pos = lerp(c_a, c_b, t);
                let (sx, sy) = tr.px(pos[0], pos[1]);

                let (color, label): (&str, String) = match dof_map.theta_row_opt(e, side, sub) {
                    Some(row) => ("#c0392b", row.to_string()),
                    None => ("#bbb", "clamp".to_string()),
                };

                let (lx, ly) = offset_toward((sx, sy), centroid_px, 15.0);
                draw_label(&mut svg, lx, ly, &label, color);
            }
        }
    }

    // ── 4. Legend ─────────────────────────────────────────────────────────
    let lx0 = 10.0;
    let mut ly = 14.0;
    let step = 18.0;

    for (symbol, color, desc) in [
        ("X", "#1a73e8", "w-DOF (free)"),
        ("X", "#999",    "w-DOF (pinned)"),
        ("O", "#c0392b", "θ-DOF (active)"),
        ("O", "#bbb",    "θ-DOF (clamped)"),
        ("A", "#27ae60", "CCW element (pos area)"),
        ("A", "#e67e22", "CW element (neg area)"),
    ] {
        let mx = lx0 + 7.0;
        let my = ly;
        match symbol {
            "X" => draw_x(&mut svg, mx, my, 5.0, color),
            "O" => draw_circle(&mut svg, mx, my, 4.5, color),
            _   => draw_arrowhead(&mut svg, mx + 4.0, my, 1.0, 0.0, 7.0, color),
        }
        svg.push_str(&format!(
            "<text x=\"{:.0}\" y=\"{:.1}\" font-size=\"10\" font-family=\"sans-serif\" \
             dominant-baseline=\"middle\" fill=\"#444\">{desc}</text>\n",
            lx0 + 18.0,
            ly,
        ));
        ly += step;
    }

    svg.push_str("</svg>\n");
    svg
}

// ── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use plato_mesh::geometry::{Panel, Polygon2D};
    use plato_mesh::mesher::{Mesher, SpadeMesher};
    use plato_mesh::model::MeshConfig;

    use super::super::dof_map::collect_all_displacement_nodes;
    use super::*;

    /// Write a DOF-map SVG for a boundary-pinned 2×2 m square mesh.
    ///
    /// Run with:
    /// ```sh
    /// cargo test -p plato-core dump_dof_map_svg -- --ignored --nocapture
    /// ```
    /// Then open `target/dof_map_debug.svg` in a browser.
    #[test]
    #[ignore]
    fn dump_dof_map_svg() {
        let panel = Panel::new("p", Polygon2D::rectangle(0.0, 0.0, 2.0, 2.0));
        let config = MeshConfig {
            max_element_area: 0.25,
            ..MeshConfig::default()
        };
        let mesh = SpadeMesher.triangulate(&[panel], &[], &config).unwrap();

        // Pin all boundary nodes (x=0, x=2, y=0, y=2).
        let tol = 1e-9;
        let pinned: HashSet<usize> = mesh
            .nodes
            .iter()
            .enumerate()
            .filter_map(|(id, &[x, y])| {
                if !(tol..=2.0 - tol).contains(&x) || !(tol..=2.0 - tol).contains(&y) {
                    Some(id)
                } else {
                    None
                }
            })
            .collect();

        let all_nodes = collect_all_displacement_nodes(&mesh);
        let dof_map = DofMap::new(&all_nodes, &pinned, &HashSet::new(), &mesh);

        let svg = render_dof_map_svg(&mesh, &dof_map);

        // Write to <workspace_root>/target/dof_map_debug.svg
        let out = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../target/dof_map_debug.svg");
        std::fs::create_dir_all(out.parent().unwrap()).ok();
        std::fs::write(&out, &svg).unwrap();

        println!(
            "SVG written → {}",
            out.canonicalize().unwrap_or(out).display()
        );
        println!(
            "  {} elements | {} nodes | n_free_w={} | n_rot={} | total_rows={}",
            mesh.elements.len(),
            mesh.nodes.len(),
            dof_map.n_free_w,
            dof_map.n_rot,
            dof_map.total_rows(),
        );
    }
}
