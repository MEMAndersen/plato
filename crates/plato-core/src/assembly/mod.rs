use std::collections::{HashMap, HashSet};

use clarabel::algebra::CscMatrix;
use plato_mesh::model::MeshModel;

use crate::element::PlateElement;

// ── DofMap ────────────────────────────────────────────────────────────────

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
/// - Rows `n_free_w..n_free_w + n_rot`: rotation DOFs (always free;
///   6 per element in the order θ0a, θ0b, θ1a, θ1b, θ2a, θ2b)
#[derive(Debug)]
pub struct DofMap {
    node_dof: HashMap<usize, DofState>,
    /// Number of free displacement DOF rows.
    pub n_free_w: usize,
    /// Number of rotation DOF rows (= 6 × N_elements).
    pub n_rot: usize,
}

impl DofMap {
    /// Build a DOF map from an ordered list of all displacement node IDs and
    /// a set of pinned node IDs.
    ///
    /// `all_displacement_nodes` is iterated in order; the first occurrence of
    /// each node ID determines its row index assignment (subsequent duplicates
    /// are ignored). Pinned nodes get `DofState::Constrained`; all others get
    /// contiguous `DofState::Free(idx)` indices starting from 0.
    pub fn new(
        all_displacement_nodes: &[usize],
        pinned: &HashSet<usize>,
        n_elements: usize,
    ) -> Self {
        let mut node_dof: HashMap<usize, DofState> = HashMap::new();
        let mut counter = 0usize;

        for &id in all_displacement_nodes {
            if node_dof.contains_key(&id) {
                continue; // already assigned
            }
            let state = if pinned.contains(&id) {
                DofState::Constrained
            } else {
                let idx = counter;
                counter += 1;
                DofState::Free(idx)
            };
            node_dof.insert(id, state);
        }

        Self {
            node_dof,
            n_free_w: counter,
            n_rot: 6 * n_elements,
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
        self.n_free_w + 6 * element_id + 2 * side + sub
    }

    /// Total number of rows in the assembled global B^T.
    #[inline]
    pub fn total_rows(&self) -> usize {
        self.n_free_w + self.n_rot
    }
}

// ── Node collection helper ────────────────────────────────────────────────

/// Collect all displacement node IDs from a mesh in a deterministic,
/// first-seen order: for each element, corners[0..3] then midpoints[0..3].
///
/// The returned vec is suitable as the `all_displacement_nodes` argument to
/// `DofMap::new`.
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

// ── Global assembly ───────────────────────────────────────────────────────

/// Assemble the global equilibrium matrix **B^T** in CSC format.
///
/// The returned `CscMatrix<f64>` has shape `dof_map.total_rows() × 9·N_elements`.
/// Each element `e` contributes columns `9e..9e+8`, with rows mapped via
/// `dof_map`. Constrained displacement DOFs are excluded (their rows are absent).
/// Rotation DOFs are always included.
///
/// Shared nodes between adjacent elements have their B^T contributions summed
/// automatically (Clarabel's `new_from_triplets` handles duplicates).
pub fn assemble_bt(mesh: &MeshModel, dof_map: &DofMap) -> CscMatrix<f64> {
    let n_rows = dof_map.total_rows();
    let n_cols = 9 * mesh.elements.len();

    // Pre-allocate: empirically ~40-60 non-zeros per element.
    let cap = 60 * mesh.elements.len();
    let mut t_rows: Vec<usize> = Vec::with_capacity(cap);
    let mut t_cols: Vec<usize> = Vec::with_capacity(cap);
    let mut t_vals: Vec<f64> = Vec::with_capacity(cap);

    for elem in &mesh.elements {
        let e = elem.id;

        // Corner coordinates for element geometry.
        let corners: [[f64; 2]; 3] = std::array::from_fn(|i| mesh.nodes[elem.corners[i]]);

        let plate = PlateElement::new(elem, &corners);
        let bt_local = plate.local_bt(); // SMatrix<f64, 12, 9>

        // Map local rows 0-11 to global rows (None = constrained, skip).
        let global_rows: [Option<usize>; 12] = [
            dof_map.w_row(elem.corners[0]),
            dof_map.w_row(elem.corners[1]),
            dof_map.w_row(elem.corners[2]),
            dof_map.w_row(elem.midpoints[0]),
            dof_map.w_row(elem.midpoints[1]),
            dof_map.w_row(elem.midpoints[2]),
            Some(dof_map.theta_row(e, 0, 0)), // θ0a
            Some(dof_map.theta_row(e, 0, 1)), // θ0b
            Some(dof_map.theta_row(e, 1, 0)), // θ1a
            Some(dof_map.theta_row(e, 1, 1)), // θ1b
            Some(dof_map.theta_row(e, 2, 0)), // θ2a
            Some(dof_map.theta_row(e, 2, 1)), // θ2b
        ];

        for local_row in 0..12usize {
            let Some(global_row) = global_rows[local_row] else {
                continue;
            };
            for local_col in 0..9usize {
                let val = bt_local[(local_row, local_col)];
                if val.abs() > 16.0 * f64::EPSILON {
                    t_rows.push(global_row);
                    t_cols.push(9 * e + local_col);
                    t_vals.push(val);
                }
            }
        }
    }

    CscMatrix::new_from_triplets(n_rows, n_cols, t_rows, t_cols, t_vals)
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
        let all_nodes: Vec<usize> = (0..6).collect();
        let pinned: HashSet<usize> = HashSet::new();
        let dm = DofMap::new(&all_nodes, &pinned, 1);

        assert_eq!(dm.n_free_w, 6);
        assert_eq!(dm.n_rot, 6);
        assert_eq!(dm.total_rows(), 12);
        for i in 0..6 {
            assert_eq!(dm.w_row(i), Some(i));
        }
    }

    #[test]
    fn dofmap_constrained_nodes() {
        let all_nodes: Vec<usize> = (0..6).collect();
        let pinned: HashSet<usize> = [0usize, 2, 5].into_iter().collect();
        let dm = DofMap::new(&all_nodes, &pinned, 1);

        assert_eq!(dm.n_free_w, 3); // nodes 1, 3, 4 are free
        assert_eq!(dm.w_row(0), None);
        assert_eq!(dm.w_row(1), Some(0));
        assert_eq!(dm.w_row(2), None);
        assert_eq!(dm.w_row(3), Some(1));
        assert_eq!(dm.w_row(4), Some(2));
        assert_eq!(dm.w_row(5), None);
        // Rotation DOFs start right after free displacement DOFs.
        assert_eq!(dm.theta_row(0, 0, 0), 3);
    }

    #[test]
    fn assembly_shape() {
        let mesh = simple_mesh();
        let all_nodes = collect_all_displacement_nodes(&mesh);
        let dm = DofMap::new(&all_nodes, &HashSet::new(), mesh.elements.len());
        let bt = assemble_bt(&mesh, &dm);

        assert_eq!(bt.m, dm.total_rows(), "wrong row count");
        assert_eq!(bt.n, 9 * mesh.elements.len(), "wrong col count");
    }

    #[test]
    fn assembly_sparsity() {
        let mesh = simple_mesh();
        let all_nodes = collect_all_displacement_nodes(&mesh);
        let dm = DofMap::new(&all_nodes, &HashSet::new(), mesh.elements.len());
        let bt = assemble_bt(&mesh, &dm);

        // Each column belongs to exactly one element; at most 12 active rows.
        for j in 0..bt.n {
            let nnz = bt.colptr[j + 1] - bt.colptr[j];
            assert!(nnz <= 12, "col {j} has {nnz} non-zeros, expected ≤ 12");
        }
    }
}
