use clarabel::algebra::CscMatrix;
use plato_mesh::model::MeshModel;

use crate::element::PlateElement;

use super::dof_map::DofMap;

/// Assemble the global equilibrium matrix **B^T** in CSC format.
///
/// Shape: `dof_map.total_rows() × 9·N_elements`.
///
/// Each element `e` contributes columns `9·e..9·e+9`.  Shared nodes have their
/// contributions summed via `CscMatrix::new_from_triplets`.
///
/// Local row mapping:
/// - Local rows 0-2: displacement at corners → `dof_map.w_row(corner_id)`
/// - Local rows 3-5: displacement at midpoints → `dof_map.w_row(midpoint_id)`
/// - Local rows 6-11: rotation DOFs → `dof_map.theta_row(e, side, sub)`
///
/// For interior edges the theta rows receive contributions from **both** adjacent
/// elements (with the sign embedded in `B^T_θ`: sub=0 → +, sub=1 → −).  This
/// implements the moment-continuity condition automatically.
pub fn assemble_bt(mesh: &MeshModel, dof_map: &DofMap) -> CscMatrix<f64> {
    let n_rows = dof_map.total_rows();
    let n_cols = 9 * mesh.elements.len();

    let cap = 60 * mesh.elements.len();
    let mut t_rows: Vec<usize> = Vec::with_capacity(cap);
    let mut t_cols: Vec<usize> = Vec::with_capacity(cap);
    let mut t_vals: Vec<f64> = Vec::with_capacity(cap);

    for elem in &mesh.elements {
        let e = elem.id;

        let corners: [[f64; 2]; 3] = std::array::from_fn(|i| mesh.nodes[elem.corners[i]]);
        let plate = PlateElement::new(elem, &corners);
        let bt_local = plate.local_bt();

        // Map local rows 0-11 to global rows (None = constrained, skip).
        let global_rows: [Option<usize>; 12] = [
            dof_map.w_row(elem.corners[0]),
            dof_map.w_row(elem.corners[1]),
            dof_map.w_row(elem.corners[2]),
            dof_map.w_row(elem.midpoints[0]),
            dof_map.w_row(elem.midpoints[1]),
            dof_map.w_row(elem.midpoints[2]),
            dof_map.theta_row_opt(e, 0, 0), // θ0a
            dof_map.theta_row_opt(e, 0, 1), // θ0b
            dof_map.theta_row_opt(e, 1, 0), // θ1a
            dof_map.theta_row_opt(e, 1, 1), // θ1b
            dof_map.theta_row_opt(e, 2, 0), // θ2a
            dof_map.theta_row_opt(e, 2, 1), // θ2b
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
    use std::collections::HashSet;

    use plato_mesh::geometry::{Panel, Polygon2D};
    use plato_mesh::mesher::{Mesher, SpadeMesher};
    use plato_mesh::model::{MeshConfig, MeshModel, MeshQuality, TriElement6};

    use super::super::dof_map::collect_all_displacement_nodes;
    use super::*;

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
    fn assembly_shape() {
        let mesh = simple_mesh();
        let all_nodes = collect_all_displacement_nodes(&mesh);
        let dm = DofMap::new(&all_nodes, &HashSet::new(), &HashSet::new(), &mesh);
        let bt = assemble_bt(&mesh, &dm);

        assert_eq!(bt.m, dm.total_rows(), "wrong row count");
        assert_eq!(bt.n, 9 * mesh.elements.len(), "wrong col count");
    }

    #[test]
    fn assembly_sparsity() {
        let mesh = simple_mesh();
        let all_nodes = collect_all_displacement_nodes(&mesh);
        let dm = DofMap::new(&all_nodes, &HashSet::new(), &HashSet::new(), &mesh);
        let bt = assemble_bt(&mesh, &dm);

        // Each column belongs to one element; at most 12 active rows before
        // shared-edge contributions from adjacent elements add to the same rows.
        // With sharing, a rotation row can receive contributions from 2 elements,
        // so per-column nnz is still ≤ 12.
        for j in 0..bt.n {
            let nnz = bt.colptr[j + 1] - bt.colptr[j];
            assert!(nnz <= 12, "col {j} has {nnz} non-zeros, expected ≤ 12");
        }
    }

    /// Manually build a 2-element mesh (two right triangles covering the unit square)
    /// and verify that B^T · m_uniform = 0 at all displacement DOF rows.
    ///
    /// This is the global patch test: for a uniform moment field (all mx=my=mxy=C),
    /// the equilibrium residual at every free displacement node must be zero
    /// (uniform moments carry no distributed load).
    #[test]
    fn two_element_patch_test_uniform_moments() {
        let nodes: Vec<[f64; 2]> = vec![
            [0.0, 0.0], // 0
            [1.0, 0.0], // 1
            [0.0, 1.0], // 2
            [1.0, 1.0], // 3
            [0.5, 0.0], // 4  midpoint A[0]
            [0.5, 0.5], // 5  midpoint A[1] = B[2]  (shared)
            [0.0, 0.5], // 6  midpoint A[2]
            [1.0, 0.5], // 7  midpoint B[0]
            [0.5, 1.0], // 8  midpoint B[1]
        ];

        let elem_a = TriElement6 {
            id: 0,
            corners: [0, 1, 2],
            midpoints: [4, 5, 6],
            panel_id: "p".into(),
        };
        let elem_b = TriElement6 {
            id: 1,
            corners: [1, 3, 2],
            midpoints: [7, 8, 5],
            panel_id: "p".into(),
        };

        let mesh = MeshModel {
            nodes,
            elements: vec![elem_a, elem_b],
            edge_nodes: Default::default(),
            quality: MeshQuality {
                n_elements: 2,
                n_nodes: 9,
                min_angle_deg: 45.0,
                max_aspect_ratio: 2.0,
            },
        };

        // Pin all boundary nodes (x=0, x=1, y=0, y=1).
        let all_nodes = collect_all_displacement_nodes(&mesh);
        let pinned: HashSet<usize> = mesh
            .nodes
            .iter()
            .enumerate()
            .filter_map(|(id, &[x, y])| {
                if !(1e-9..=1.0 - 1e-9).contains(&x) || !(1e-9..=1.0 - 1e-9).contains(&y) {
                    Some(id)
                } else {
                    None
                }
            })
            .collect();

        assert_eq!(pinned.len(), 8, "8 boundary nodes should be pinned");

        let dm = DofMap::new(&all_nodes, &pinned, &HashSet::new(), &mesh);
        assert_eq!(dm.n_free_w, 1, "only the shared interior midpoint is free");

        let bt = assemble_bt(&mesh, &dm);

        assert_eq!(bt.m, dm.total_rows());
        assert_eq!(bt.n, 18);

        // Patch test: B^T · m_uniform at displacement DOF rows must be zero.
        let m_uniform = [1.0_f64; 18];
        let n_rows = bt.m;
        let mut result = vec![0.0f64; n_rows];
        for (col, &val) in m_uniform.iter().enumerate() {
            for ptr in bt.colptr[col]..bt.colptr[col + 1] {
                result[bt.rowval[ptr]] += bt.nzval[ptr] * val;
            }
        }

        for (r, &res) in result.iter().enumerate().take(dm.n_free_w) {
            assert!(
                res.abs() < 1e-10,
                "B^T · m_uniform at displacement row {r} = {:.2e} (expected 0)",
                res
            );
        }
    }
}
