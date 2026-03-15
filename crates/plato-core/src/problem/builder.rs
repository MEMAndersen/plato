use clarabel::algebra::CscMatrix;
use clarabel::solver::{
    DefaultSettings, DefaultSolver, IPSolver, NonnegativeConeT, SecondOrderConeT, SolverStatus,
    SupportedConeT, ZeroConeT,
};

use crate::criteria::JohansenCriterion;
use crate::result::{ElementMoments, SolveResult, SolveStatus};

// ── Variable layout ───────────────────────────────────────────────────────────

/// Maps logical variables to column indices in the Clarabel `x` vector.
///
/// Column layout:
/// - Col 0: λ (load factor)
/// - Cols 1 .. 1+9·N_e: moment components [mx, my, mxy] at 3 corners per element
/// - Cols 1+9·N_e .. 1+21·N_e: α auxiliaries (4 per corner × 3 corners per element)
pub struct VarLayout {
    pub n_e: usize,
}

impl VarLayout {
    pub fn n_vars(&self) -> usize {
        1 + 21 * self.n_e
    }

    /// Column for moment component `c` (0=mx, 1=my, 2=mxy) at corner `k` of element `e`.
    #[inline]
    pub fn moment_col(&self, e: usize, k: usize, c: usize) -> usize {
        1 + 9 * e + 3 * k + c
    }

    /// Column for auxiliary `a` (0=α₁, 1=α₂, 2=α₃, 3=α₄) at corner `k` of element `e`.
    #[inline]
    pub fn alpha_col(&self, e: usize, k: usize, a: usize) -> usize {
        1 + 9 * self.n_e + 12 * e + 4 * k + a
    }
}

// ── ClarabelProblem builder ───────────────────────────────────────────────────

/// Inputs to one limit-analysis SOCP solve.
pub struct ClarabelProblem<'a> {
    /// Global B^T matrix, shape n_dof × 9·n_e (CSC, column-major).
    pub bt: &'a CscMatrix<f64>,
    /// Variable load vector (scaled by λ), length n_dof.
    pub q_ext: &'a [f64],
    /// Permanent load vector (constant RHS), length n_dof or empty (treated as zero).
    pub q_dead: &'a [f64],
    /// Johansen yield criteria — one reference per element.
    /// Elements sharing a panel material point to the same criterion object.
    pub criteria: &'a [&'a JohansenCriterion],
    /// Number of elements.
    pub n_e: usize,
}

impl ClarabelProblem<'_> {
    /// Assemble and solve the Clarabel SOCP. Returns the plastic load factor λ.
    pub fn solve(&self) -> SolveResult {
        let n_e = self.n_e;
        let n_dof = self.q_ext.len();

        let layout = VarLayout { n_e };
        let n_vars = layout.n_vars();

        // ── Row offsets ───────────────────────────────────────────────────
        let n_eq = n_dof + 12 * n_e; // ZeroConeT rows
        let n_nonneg = 12 * n_e; // NonnegativeConeT rows
        let n_soc_rows = 18 * n_e; // 6 SOC(3) blocks per element

        let nonneg_start = n_eq;
        let soc_start = n_eq + n_nonneg;
        let n_rows = n_eq + n_nonneg + n_soc_rows;

        // ── Build A triplets ──────────────────────────────────────────────
        let cap = 60 * n_e + 3 * 12 * n_e + 12 * n_e + 3 * 18 * n_e;
        let mut t_rows: Vec<usize> = Vec::with_capacity(cap);
        let mut t_cols: Vec<usize> = Vec::with_capacity(cap);
        let mut t_vals: Vec<f64> = Vec::with_capacity(cap);

        // Initialise b; set dead-load values at equilibrium rows.
        let mut b = vec![0.0f64; n_rows];
        for (r, &dval) in self.q_dead.iter().enumerate() {
            b[r] = dval;
        }

        macro_rules! push {
            ($r:expr, $c:expr, $v:expr) => {{
                t_rows.push($r);
                t_cols.push($c);
                t_vals.push($v);
            }};
        }

        // ── 1. Equilibrium rows (ZeroConeT, rows 0..n_dof) ───────────────
        // A[r, 0] = -q_ext[r]  (λ column)
        for (r, &fval) in self.q_ext.iter().enumerate() {
            if fval.abs() > 0.0 {
                push!(r, 0, -fval);
            }
        }
        // A[r, 1+j] = B^T[r, j]
        for col_j in 0..self.bt.n {
            for ptr in self.bt.colptr[col_j]..self.bt.colptr[col_j + 1] {
                let row_r = self.bt.rowval[ptr];
                let val = self.bt.nzval[ptr];
                push!(row_r, 1 + col_j, val);
            }
        }

        // ── 2/3/4. Yield criterion rows (α defs, nonneg, SOC) ────────────
        // Delegated to JohansenCriterion::push_corner_blocks() per corner.
        for e in 0..n_e {
            for k in 0..3usize {
                self.criteria[e].push_corner_blocks(
                    &mut t_rows,
                    &mut t_cols,
                    &mut t_vals,
                    &mut b,
                    n_dof + 12 * e + 4 * k,
                    nonneg_start + 12 * e + 4 * k,
                    soc_start + (6 * e + 2 * k) * 3,
                    layout.moment_col(e, k, 0),
                    layout.moment_col(e, k, 1),
                    layout.moment_col(e, k, 2),
                    layout.alpha_col(e, k, 0),
                    layout.alpha_col(e, k, 1),
                    layout.alpha_col(e, k, 2),
                    layout.alpha_col(e, k, 3),
                );
            }
        }

        // ── Assemble CscMatrix A ──────────────────────────────────────────
        let a_mat = CscMatrix::new_from_triplets(n_rows, n_vars, t_rows, t_cols, t_vals);

        // ── Objective: minimise −λ ────────────────────────────────────────
        let p_mat = CscMatrix::zeros((n_vars, n_vars));
        let mut q = vec![0.0f64; n_vars];
        q[0] = -1.0;

        // ── Cone list ─────────────────────────────────────────────────────
        let mut cones: Vec<SupportedConeT<f64>> = Vec::with_capacity(2 + 6 * n_e);
        cones.push(ZeroConeT(n_eq));
        cones.push(NonnegativeConeT(n_nonneg));
        for _ in 0..6 * n_e {
            cones.push(SecondOrderConeT(3));
        }

        // ── Solve ─────────────────────────────────────────────────────────
        let settings = DefaultSettings::default();
        let mut solver = DefaultSolver::new(&p_mat, &q, &a_mat, &b, &cones, settings);
        solver.solve();

        let sol = &solver.solution;
        let status = match sol.status {
            SolverStatus::Solved => SolveStatus::Solved,
            SolverStatus::AlmostSolved => SolveStatus::AlmostSolved,
            SolverStatus::PrimalInfeasible | SolverStatus::AlmostPrimalInfeasible => {
                SolveStatus::Infeasible
            }
            _ => SolveStatus::Failed,
        };

        let load_factor = if matches!(status, SolveStatus::Solved | SolveStatus::AlmostSolved) {
            sol.x[0]
        } else {
            f64::NAN
        };

        // ── Extract element moments ───────────────────────────────────────
        let element_moments = if matches!(status, SolveStatus::Solved | SolveStatus::AlmostSolved) {
            (0..n_e)
                .map(|e| {
                    let mut corner_moments = [[0.0f64; 3]; 3];
                    for (k, corner) in corner_moments.iter_mut().enumerate() {
                        for (comp, val) in corner.iter_mut().enumerate() {
                            *val = sol.x[layout.moment_col(e, k, comp)];
                        }
                    }
                    let yield_utilisation = std::array::from_fn(|k| {
                        let [mx, my, mxy] = corner_moments[k];
                        self.criteria[e].yield_utilisation(mx, my, mxy)
                    });
                    ElementMoments {
                        element_id: e,
                        corner_moments,
                        yield_utilisation,
                    }
                })
                .collect()
        } else {
            Vec::new()
        };

        let duality_gap = (sol.obj_val - sol.obj_val_dual).abs();

        // ── Post-solve equilibrium diagnostic (test builds only) ──────────
        #[cfg(test)]
        if matches!(status, SolveStatus::Solved | SolveStatus::AlmostSolved) {
            let lf = load_factor;
            let n_dof = self.q_ext.len();

            // Compute B^T · m  at equilibrium rows
            let mut bt_m = vec![0.0f64; n_dof];
            for col_j in 0..self.bt.n {
                let m_val = sol.x[1 + col_j];
                for ptr in self.bt.colptr[col_j]..self.bt.colptr[col_j + 1] {
                    let row_r = self.bt.rowval[ptr];
                    if row_r < n_dof {
                        bt_m[row_r] += self.bt.nzval[ptr] * m_val;
                    }
                }
            }
            // Residual: B^T·m − λ·q_ext − q_dead
            let mut max_res = 0.0f64;
            let mut sum_res = 0.0f64;
            let mut n_eq = 0usize;
            for r in 0..n_dof {
                let q_val = self.q_ext[r];
                let d_val = self.q_dead.get(r).copied().unwrap_or(0.0);
                let res = (bt_m[r] - lf * q_val - d_val).abs();
                max_res = max_res.max(res);
                sum_res += res;
                n_eq += 1;
            }

            // Yield utilisation stats across all element corners
            let utils: Vec<f64> = element_moments
                .iter()
                .flat_map(|em| em.yield_utilisation.iter().copied())
                .collect();
            let util_max = utils.iter().cloned().fold(0.0f64, f64::max);
            let util_mean = utils.iter().sum::<f64>() / utils.len() as f64;

            // Moment magnitude stats
            let m_max = element_moments
                .iter()
                .flat_map(|em| em.corner_moments.iter())
                .flat_map(|m| m.iter().copied())
                .fold(0.0f64, |a, b| a.max(b.abs()));
            let m_mean = {
                let vals: Vec<f64> = element_moments
                    .iter()
                    .flat_map(|em| em.corner_moments.iter())
                    .flat_map(|m| m.iter().copied())
                    .collect();
                vals.iter().sum::<f64>() / vals.len() as f64
            };

            // q_ext stats
            let q_max = self.q_ext.iter().cloned().fold(0.0f64, f64::max);
            let q_sum: f64 = self.q_ext.iter().sum();

            eprintln!(
                "[eq-diag] λ={lf:.4} n_dof={n_eq} n_e={n_e} \
                 | eq_res max={max_res:.4e} mean={:.4e} \
                 | yield_util max={util_max:.4} mean={util_mean:.4} \
                 | |m| max={m_max:.4e} mean={m_mean:.4e} \
                 | q_ext max={q_max:.4e} sum={q_sum:.4e}",
                sum_res / n_eq as f64,
            );
        }

        SolveResult {
            status,
            load_factor,
            element_moments,
            solver_iterations: sol.iterations,
            duality_gap,
        }
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use plato_mesh::model::{MeshModel, MeshQuality, TriElement6};

    use crate::assembly::{DofMap, assemble_bt, collect_all_displacement_nodes};

    use super::*;

    /// Two right-angle triangles covering the unit square (1×1 m).
    /// All 8 boundary nodes pinned; 1 interior midpoint free.
    /// Uniform unit load (1 kN/m²). mp = 1 kNm/m (isotropic).
    ///
    /// Checks solver convergence and λ > 0. Does not assert a specific value
    /// (accurate λ requires proper BCs and load distribution — Phase 4 scope).
    #[test]
    fn two_element_square_socp_smoke() {
        let nodes: Vec<[f64; 2]> = vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0], // corners
            [0.5, 0.0],
            [0.5, 0.5],
            [0.0, 0.5], // A midpoints
            [1.0, 0.5],
            [0.5, 1.0], // B midpoints
        ];
        let mesh = MeshModel {
            nodes,
            elements: vec![
                TriElement6 {
                    id: 0,
                    corners: [0, 1, 2],
                    midpoints: [4, 5, 6],
                    panel_id: "p".into(),
                },
                TriElement6 {
                    id: 1,
                    corners: [1, 3, 2],
                    midpoints: [7, 8, 5],
                    panel_id: "p".into(),
                },
            ],
            edge_nodes: Default::default(),
            quality: MeshQuality {
                n_elements: 2,
                n_nodes: 9,
                min_angle_deg: 45.0,
                max_aspect_ratio: 2.0,
            },
        };

        // Pin all boundary nodes.
        let tol = 1e-9;
        let pinned: HashSet<usize> = mesh
            .nodes
            .iter()
            .enumerate()
            .filter_map(|(id, &[x, y])| {
                if !(tol..=1.0 - tol).contains(&x) || !(tol..=1.0 - tol).contains(&y) {
                    Some(id)
                } else {
                    None
                }
            })
            .collect();

        let all_nodes = collect_all_displacement_nodes(&mesh);
        let dof_map = DofMap::new(&all_nodes, &pinned, &HashSet::new(), &mesh);
        let n_dof = dof_map.total_rows();
        let bt = assemble_bt(&mesh, &dof_map);

        // Uniform unit load: each element midpoint gets area / 3.
        let mut f_ref = vec![0.0f64; n_dof];
        for elem in &mesh.elements {
            let p0 = mesh.nodes[elem.corners[0]];
            let p1 = mesh.nodes[elem.corners[1]];
            let p2 = mesh.nodes[elem.corners[2]];
            let area =
                0.5 * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1])).abs();
            let f_node = area / 3.0;
            for &mid_id in &elem.midpoints {
                if let Some(row) = dof_map.w_row(mid_id) {
                    f_ref[row] += f_node;
                }
            }
        }

        let criterion = JohansenCriterion::isotropic(1.0);
        let criteria = [&criterion; 2];
        let problem = ClarabelProblem {
            bt: &bt,
            q_ext: &f_ref,
            q_dead: &[],
            criteria: &criteria,
            n_e: 2,
        };
        let result = problem.solve();

        assert!(
            matches!(
                result.status,
                SolveStatus::Solved | SolveStatus::AlmostSolved
            ),
            "solver status: {:?}",
            result.status
        );
        let lambda = result.load_factor;
        println!("2-element smoke test λ = {lambda:.4}  (exact upper bound = 24.0)");
        assert!(
            lambda > 0.0 && lambda.is_finite(),
            "λ must be positive finite, got {lambda}"
        );
    }
}
