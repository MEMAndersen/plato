use std::collections::{HashMap, HashSet};

use plato_mesh::model::MeshModel;

/// State of a displacement DOF (w) at a corner or midpoint node.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DofState {
    /// Free DOF — carries its contiguous row index in the global B^T.
    Free(usize),
    /// Pinned support — row excluded from the global B^T.
    Constrained,
}

/// Maps corner/midpoint node IDs and element rotation DOFs to rows in the
/// assembled global B^T equilibrium matrix.
///
/// **Row layout:**
/// - Rows `0..n_free_w`: free displacement DOFs (w at corner + midpoint nodes)
/// - Rows `n_free_w..n_free_w + n_rot`: rotation DOFs (2 per unique edge)
///
/// **Theta DOF sharing:**
/// Each mesh edge has exactly 2 rotation DOF rows, one per endpoint of that edge.
/// For *interior* edges shared by two elements the two elements contribute to the
/// *same* rows with opposite signs — this enforces bending-moment continuity.
/// For *boundary* edges the single element's contribution enforces zero normal
/// moment (correct for simply-supported BCs).
///
/// Internally `theta_dof[elem.id][side][sub]` gives the global row index.
#[derive(Debug)]
pub struct DofMap {
    node_dof: HashMap<usize, DofState>,
    /// Number of free displacement DOF rows.
    pub n_free_w: usize,
    /// Number of rotation DOF rows (= 2 × number of unique mesh edges).
    pub n_rot: usize,
    /// theta_dof[elem_id][side 0..3][sub 0..2]  →  global row index
    theta_dof: HashMap<usize, [[usize; 2]; 3]>,
}

impl DofMap {
    /// Build a DOF map from the mesh.
    ///
    /// `all_displacement_nodes` is iterated in order; the first occurrence of
    /// each node ID determines its free row index. Pinned nodes get
    /// `DofState::Constrained`.
    ///
    /// For each mesh edge (identified by its shared midpoint node):
    /// - **Interior edge** (midpoint appears in two elements): the two elements
    ///   share the same pair of theta DOF rows, with swapped sub-indices so that
    ///   their combined equation enforces normal-moment continuity.
    /// - **Boundary edge** (midpoint appears in one element): two fresh theta DOF
    ///   rows are assigned to that element alone; their equation enforces zero
    ///   normal moment (simply-supported BC).
    ///
    /// **Side–midpoint convention (0-indexed):**
    /// Side `i` of an element is opposite corner `i` and spans corners
    /// `(i+1)%3` → `(i+2)%3`.  Its midpoint node is `elem.midpoints[(i+1)%3]`.
    pub fn new(
        all_displacement_nodes: &[usize],
        pinned: &HashSet<usize>,
        mesh: &MeshModel,
    ) -> Self {
        // ── Displacement DOF assignment ───────────────────────────────────
        let mut node_dof: HashMap<usize, DofState> = HashMap::new();
        let mut free_counter = 0usize;

        for &id in all_displacement_nodes {
            if node_dof.contains_key(&id) {
                continue;
            }
            let state = if pinned.contains(&id) {
                DofState::Constrained
            } else {
                let idx = free_counter;
                free_counter += 1;
                DofState::Free(idx)
            };
            node_dof.insert(id, state);
        }
        let n_free_w = free_counter;

        // ── Theta DOF assignment (shared between adjacent elements) ───────
        // mid_to_side[midpoint_node_id] = (elem_id, side_index) of the first
        // element that claimed this edge.  When a second element finds the same
        // midpoint it knows the edge is shared and inherits the already-assigned
        // theta rows with swapped sub-indices.
        let mut mid_to_side: HashMap<usize, (usize, usize)> = HashMap::new();
        let mut theta_dof: HashMap<usize, [[usize; 2]; 3]> =
            HashMap::with_capacity(mesh.elements.len());
        let mut rot_counter = n_free_w;

        for elem in &mesh.elements {
            let e = elem.id;
            let mut elem_theta = [[0usize; 2]; 3];

            #[allow(clippy::needless_range_loop)]
            for side in 0..3usize {
                // Side i's midpoint node is elem.midpoints[(i+1)%3].
                let mid_id = elem.midpoints[(side + 1) % 3];

                if let Some(&(adj_e, adj_side)) = mid_to_side.get(&mid_id) {
                    // ── Interior edge: inherit rows from the first element,
                    //    but with sub-indices SWAPPED.
                    //
                    //  A's (side, sub=0) ↔ B's (adj_side, sub=1) → same row
                    //  A's (side, sub=1) ↔ B's (adj_side, sub=0) → same row
                    //
                    // This gives A_sub0 · (+a) + B_sub1 · (-a) = a·(m_A - m_B) = 0.
                    let adj = theta_dof[&adj_e];
                    elem_theta[side][0] = adj[adj_side][1]; // e sub=0 ← adj sub=1
                    elem_theta[side][1] = adj[adj_side][0]; // e sub=1 ← adj sub=0
                } else {
                    // ── New (boundary or first-seen) edge: assign fresh rows.
                    elem_theta[side][0] = rot_counter;
                    rot_counter += 1;
                    elem_theta[side][1] = rot_counter;
                    rot_counter += 1;
                    mid_to_side.insert(mid_id, (e, side));
                }
            }

            theta_dof.insert(e, elem_theta);
        }

        let n_rot = rot_counter - n_free_w;

        Self {
            node_dof,
            n_free_w,
            n_rot,
            theta_dof,
        }
    }

    /// Global row index for the displacement DOF at `node_id`,
    /// or `None` if the node is constrained.
    pub fn w_row(&self, node_id: usize) -> Option<usize> {
        match self.node_dof.get(&node_id) {
            Some(DofState::Free(idx)) => Some(*idx),
            _ => None,
        }
    }

    /// Global row index for element-local rotation DOF
    /// (`element_id`, `side` ∈ 0..3, `sub` ∈ 0..2 for a/b).
    #[inline]
    pub fn theta_row(&self, element_id: usize, side: usize, sub: usize) -> usize {
        self.theta_dof[&element_id][side][sub]
    }

    /// Total number of rows in the assembled global B^T.
    #[inline]
    pub fn total_rows(&self) -> usize {
        self.n_free_w + self.n_rot
    }
}

/// Collect all displacement node IDs from a mesh in a deterministic,
/// first-seen order: for each element, corners[0..3] then midpoints[0..3].
pub fn collect_all_displacement_nodes(mesh: &MeshModel) -> Vec<usize> {
    let mut seen = HashSet::new();
    let mut ordered = Vec::new();
    for elem in &mesh.elements {
        for &id in elem.corners.iter().chain(elem.midpoints.iter()) {
            if seen.insert(id) {
                ordered.push(id);
            }
        }
    }
    ordered
}

// ── Tests ─────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use plato_mesh::geometry::{Panel, Polygon2D};
    use plato_mesh::mesher::{Mesher, SpadeMesher};
    use plato_mesh::model::MeshConfig;

    fn simple_mesh() -> MeshModel {
        let panel = Panel::new("p", Polygon2D::rectangle(0.0, 0.0, 2.0, 2.0));
        let config = MeshConfig {
            max_element_area: 0.5,
            ..MeshConfig::default()
        };
        SpadeMesher
            .triangulate(&[panel], &[], &config)
            .expect("mesh failed")
    }

    #[test]
    fn dofmap_all_free() {
        let mesh = simple_mesh();
        let all_nodes = collect_all_displacement_nodes(&mesh);
        let pinned: HashSet<usize> = HashSet::new();
        let dm = DofMap::new(&all_nodes, &pinned, &mesh);

        // All nodes are free, so n_free_w = number of unique nodes.
        assert_eq!(dm.n_free_w, all_nodes.len());
        for (i, &id) in all_nodes.iter().enumerate() {
            assert_eq!(dm.w_row(id), Some(i));
        }
        // n_rot = 2 × number of unique mesh edges.
        // Count unique edges via unique midpoint node IDs (each edge has exactly one).
        let n_edges: usize = {
            let unique: std::collections::HashSet<usize> = mesh
                .elements
                .iter()
                .flat_map(|e| e.midpoints.iter().copied())
                .collect();
            unique.len()
        };
        assert_eq!(dm.n_rot, 2 * n_edges, "n_rot should be 2 × n_edges");
    }

    #[test]
    fn dofmap_constrained_nodes() {
        let mesh = simple_mesh();
        let all_nodes = collect_all_displacement_nodes(&mesh);
        let pinned: HashSet<usize> = all_nodes.iter().take(3).copied().collect();
        let dm = DofMap::new(&all_nodes, &pinned, &mesh);

        assert_eq!(dm.n_free_w, all_nodes.len() - 3);
        for &id in pinned.iter() {
            assert_eq!(dm.w_row(id), None);
        }
    }

    #[test]
    fn theta_rows_shared_on_interior_edges() {
        let mesh = simple_mesh();
        let all_nodes = collect_all_displacement_nodes(&mesh);
        let dm = DofMap::new(&all_nodes, &HashSet::new(), &mesh);

        // For each interior edge (midpoint shared between 2 elements), confirm
        // that the theta rows of the two elements for that side are the SAME
        // but with swapped sub-indices (enforcing continuity, not double-zero).
        let mut mid_to_side: HashMap<usize, (usize, usize)> = HashMap::new();
        for elem in &mesh.elements {
            let e = elem.id;
            for side in 0..3usize {
                let mid_id = elem.midpoints[(side + 1) % 3];
                if let Some(&(adj_e, adj_side)) = mid_to_side.get(&mid_id) {
                    // Interior edge: rows should be swapped.
                    assert_eq!(
                        dm.theta_row(e, side, 0),
                        dm.theta_row(adj_e, adj_side, 1),
                        "e={e} side={side} sub=0 should equal adj_e={adj_e} adj_side={adj_side} sub=1"
                    );
                    assert_eq!(
                        dm.theta_row(e, side, 1),
                        dm.theta_row(adj_e, adj_side, 0),
                        "e={e} side={side} sub=1 should equal adj_e={adj_e} adj_side={adj_side} sub=0"
                    );
                } else {
                    mid_to_side.insert(mid_id, (e, side));
                }
            }
        }
    }
}
