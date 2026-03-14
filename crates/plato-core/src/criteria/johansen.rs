/// Johansen yield criterion for a reinforced-concrete slab.
///
/// Supports the general orthotropic case with independent positive (sagging)
/// and negative (hogging) moment capacities in x and y.
///
/// The yield surface at each corner node is defined by two product inequalities:
///
/// Positive surface (sagging):
///     (m+_x − mx)(m+_y − my) ≥ mxy²
///
/// Negative surface (hogging):
///     (m-_x + mx)(m-_y + my) ≥ mxy²
///
/// Both surfaces must hold simultaneously.
#[derive(Debug, Clone, Copy)]
pub struct JohansenCriterion {
    /// Positive (sagging) moment capacity in x [kNm/m]. Must be > 0.
    pub mpx_pos: f64,
    /// Positive (sagging) moment capacity in y [kNm/m]. Must be > 0.
    pub mpy_pos: f64,
    /// Negative (hogging) moment capacity in x [kNm/m]. Must be > 0.
    pub mpx_neg: f64,
    /// Negative (hogging) moment capacity in y [kNm/m]. Must be > 0.
    pub mpy_neg: f64,
}

impl JohansenCriterion {
    /// General constructor. All capacities must be strictly positive.
    pub fn new(mpx_pos: f64, mpy_pos: f64, mpx_neg: f64, mpy_neg: f64) -> Self {
        assert!(
            mpx_pos > 0.0 && mpy_pos > 0.0 && mpx_neg > 0.0 && mpy_neg > 0.0,
            "all moment capacities must be positive"
        );
        Self {
            mpx_pos,
            mpy_pos,
            mpx_neg,
            mpy_neg,
        }
    }

    /// Isotropic slab: equal capacity in all directions, top and bottom.
    pub fn isotropic(mp: f64) -> Self {
        Self::new(mp, mp, mp, mp)
    }

    /// Orthotropic slab: different x/y capacities but equal top/bottom.
    pub fn orthotropic(mpx: f64, mpy: f64) -> Self {
        Self::new(mpx, mpy, mpx, mpy)
    }

    /// Returns `true` if `(mx, my, mxy)` satisfies both yield surfaces.
    pub fn is_admissible(&self, mx: f64, my: f64, mxy: f64) -> bool {
        let mxy2 = mxy * mxy;
        (self.mpx_pos - mx) * (self.mpy_pos - my) >= mxy2
            && (self.mpx_neg + mx) * (self.mpy_neg + my) >= mxy2
    }

    /// Push Clarabel constraint rows for a single yield-criterion evaluation point.
    ///
    /// Appends triplets to the caller's `(t_rows, t_cols, t_vals)` buffers and sets
    /// RHS values in `b`. The caller must pre-allocate `b` with at least
    /// `soc_row + 6` entries initialised to zero.
    ///
    /// **Row layout emitted (14 rows total):**
    ///
    /// | Rows | Cone | Constraint |
    /// |---|---|---|
    /// | `eq_row + 0..4`     | ZeroCone | α₁..α₄ definitions |
    /// | `nonneg_row + 0..4` | NonnegativeCone | α₁..α₄ ≥ 0 |
    /// | `soc_row + 0..3`    | SecondOrderCone(3) | positive surface |
    /// | `soc_row + 3..6`    | SecondOrderCone(3) | negative surface |
    ///
    /// **Local variable semantics (scaled by 1/√2 for O(1) matrix entries):**
    /// - α₁ = (1/√2)(mpx_pos − mx) ≥ 0,  α₂ = (1/√2)(mpy_pos − my) ≥ 0
    /// - α₃ = (1/√2)(mpx_neg + mx) ≥ 0,  α₄ = (1/√2)(mpy_neg + my) ≥ 0
    #[allow(clippy::too_many_arguments)]
    pub fn push_corner_blocks(
        &self,
        t_rows: &mut Vec<usize>,
        t_cols: &mut Vec<usize>,
        t_vals: &mut Vec<f64>,
        b: &mut [f64],
        eq_row: usize,
        nonneg_row: usize,
        soc_row: usize,
        mx_col: usize,
        my_col: usize,
        mxy_col: usize,
        a1_col: usize,
        a2_col: usize,
        a3_col: usize,
        a4_col: usize,
    ) {
        macro_rules! push {
            ($r:expr, $c:expr, $v:expr) => {{
                t_rows.push($r);
                t_cols.push($c);
                t_vals.push($v);
            }};
        }

        // 1/√2 — scales α to keep all matrix entries at O(1/√2) for
        // better Clarabel preconditioning (design doc §7.5).
        let s = std::f64::consts::FRAC_1_SQRT_2;

        // ── Equality rows: α definitions ──────────────────────────────────
        // α₁ + (1/√2)·mx = (1/√2)·m+_x   →   α₁ = (1/√2)(m+_x − mx)
        push!(eq_row, a1_col, 1.0);
        push!(eq_row, mx_col, s);
        b[eq_row] = s * self.mpx_pos;
        // α₂ + (1/√2)·my = (1/√2)·m+_y
        push!(eq_row + 1, a2_col, 1.0);
        push!(eq_row + 1, my_col, s);
        b[eq_row + 1] = s * self.mpy_pos;
        // α₃ − (1/√2)·mx = (1/√2)·m−_x   →   α₃ = (1/√2)(m−_x + mx)
        push!(eq_row + 2, a3_col, 1.0);
        push!(eq_row + 2, mx_col, -s);
        b[eq_row + 2] = s * self.mpx_neg;
        // α₄ − (1/√2)·my = (1/√2)·m−_y
        push!(eq_row + 3, a4_col, 1.0);
        push!(eq_row + 3, my_col, -s);
        b[eq_row + 3] = s * self.mpy_neg;

        // ── Nonneg rows: α ≥ 0 ────────────────────────────────────────────
        push!(nonneg_row, a1_col, -1.0);
        push!(nonneg_row + 1, a2_col, -1.0);
        push!(nonneg_row + 2, a3_col, -1.0);
        push!(nonneg_row + 3, a4_col, -1.0);

        // ── SOC rows: Johansen yield surfaces ─────────────────────────────
        // Positive surface: α₁α₂ ≥ (1/2)·mxy²
        //   t = (1/√2)·α₁ + (1/√2)·α₂,  s_val = (1/√2)·α₁ − (1/√2)·α₂,  w = mxy
        //   [t; s_val; w] ∈ SOC(3)  ↔  t² − s² = 2α₁α₂ ≥ mxy²
        push!(soc_row, a1_col, -s);
        push!(soc_row, a2_col, -s);
        push!(soc_row + 1, a1_col, -s);
        push!(soc_row + 1, a2_col, s);
        push!(soc_row + 2, mxy_col, -1.0);
        // Negative surface: α₃α₄ ≥ (1/2)·mxy²
        push!(soc_row + 3, a3_col, -s);
        push!(soc_row + 3, a4_col, -s);
        push!(soc_row + 4, a3_col, -s);
        push!(soc_row + 4, a4_col, s);
        push!(soc_row + 5, mxy_col, -1.0);
    }

    /// Approximate yield utilisation in [0, ∞). Values ≤ 1 are admissible.
    ///
    /// Returns `1 − min(pos_slack, neg_slack) / (mpx_pos · mpy_pos)` where
    /// `pos_slack` and `neg_slack` are the margins of the two surfaces.
    pub fn yield_utilisation(&self, mx: f64, my: f64, mxy: f64) -> f64 {
        let mxy2 = mxy * mxy;
        let pos = (self.mpx_pos - mx) * (self.mpy_pos - my) - mxy2;
        let neg = (self.mpx_neg + mx) * (self.mpy_neg + my) - mxy2;
        let denom = self.mpx_pos * self.mpy_pos;
        (1.0 - pos.min(neg) / denom).max(0.0)
    }
}

#[cfg(test)]
mod tests {
    use clarabel::algebra::CscMatrix;
    use clarabel::solver::{
        DefaultSettings, DefaultSolver, IPSolver, NonnegativeConeT, SecondOrderConeT,
        SupportedConeT, ZeroConeT,
    };

    use super::*;

    // ── push_corner_blocks SOCP tests ──────────────────────────────────────

    /// Build a tiny SOCP that fixes (mx, my, mxy) and lets the criterion
    /// constraints determine α₁..α₄. Checks the solved α values match
    /// the analytical values αᵢ = (capacity − moment) / √2.
    #[test]
    fn single_point_alpha_values() {
        // Variables: [mx=0, my=1, mxy=2, a1=3, a2=4, a3=5, a4=6]
        let n_vars = 7;
        let c = JohansenCriterion::isotropic(1.0);
        let mx_val = 0.5_f64;
        let my_val = 0.3_f64;
        let mxy_val = 0.0_f64;

        // Row layout:
        //   rows 0..3  ZeroCone: fix mx, my, mxy
        //   rows 3..7  ZeroCone: α definitions  (eq_row = 3)
        //   rows 7..11 NonnegativeCone          (nonneg_row = 7)
        //   rows 11..17 2×SOC(3)                (soc_row = 11)
        let n_rows = 17;
        let mut t_rows: Vec<usize> = Vec::new();
        let mut t_cols: Vec<usize> = Vec::new();
        let mut t_vals: Vec<f64> = Vec::new();
        let mut b = vec![0.0f64; n_rows];

        // Fix moment values.
        t_rows.push(0);
        t_cols.push(0);
        t_vals.push(1.0);
        b[0] = mx_val;
        t_rows.push(1);
        t_cols.push(1);
        t_vals.push(1.0);
        b[1] = my_val;
        t_rows.push(2);
        t_cols.push(2);
        t_vals.push(1.0);
        b[2] = mxy_val;

        c.push_corner_blocks(
            &mut t_rows,
            &mut t_cols,
            &mut t_vals,
            &mut b,
            3,
            7,
            11,
            0,
            1,
            2,
            3,
            4,
            5,
            6,
        );

        let a_mat = CscMatrix::new_from_triplets(n_rows, n_vars, t_rows, t_cols, t_vals);
        let p_mat = CscMatrix::zeros((n_vars, n_vars));
        let q = vec![0.0f64; n_vars];
        let cones: Vec<SupportedConeT<f64>> = vec![
            ZeroConeT(7),        // rows 0..7 (3 fix + 4 def)
            NonnegativeConeT(4), // rows 7..11
            SecondOrderConeT(3), // rows 11..14
            SecondOrderConeT(3), // rows 14..17
        ];
        let settings = DefaultSettings {
            verbose: false,
            ..DefaultSettings::default()
        };
        let mut solver = DefaultSolver::new(&p_mat, &q, &a_mat, &b, &cones, settings);
        solver.solve();

        let x = &solver.solution.x;
        let tol = 1e-6;
        let s = std::f64::consts::FRAC_1_SQRT_2;
        // α values are scaled by 1/√2: αᵢ = (capacity − moment) / √2
        let a1_exp = (1.0 - mx_val) * s;
        let a2_exp = (1.0 - my_val) * s;
        let a3_exp = (1.0 + mx_val) * s;
        let a4_exp = (1.0 + my_val) * s;
        assert!(
            (x[3] - a1_exp).abs() < tol,
            "a1 expected {a1_exp}, got {}",
            x[3]
        );
        assert!(
            (x[4] - a2_exp).abs() < tol,
            "a2 expected {a2_exp}, got {}",
            x[4]
        );
        assert!(
            (x[5] - a3_exp).abs() < tol,
            "a3 expected {a3_exp}, got {}",
            x[5]
        );
        assert!(
            (x[6] - a4_exp).abs() < tol,
            "a4 expected {a4_exp}, got {}",
            x[6]
        );
    }

    /// Maximise λ for moments proportional to (mx=λ, my=0, mxy=0) with
    /// an orthotropic criterion (mpx=2, mpy=3). Analytical λ_max = 2.0
    /// (positive surface hits capacity: (2−λ)·3 = 0 at λ=2).
    #[test]
    fn single_point_max_lambda() {
        // Variables: [λ=0, mx=1, my=2, mxy=3, a1=4, a2=5, a3=6, a4=7]
        let n_vars = 8;
        let c = JohansenCriterion::orthotropic(2.0, 3.0);

        // Row layout:
        //   rows 0..3  ZeroCone: mx=λ, my=0, mxy=0
        //   rows 3..7  ZeroCone: α definitions  (eq_row = 3)
        //   rows 7..11 NonnegativeCone          (nonneg_row = 7)
        //   rows 11..17 2×SOC(3)                (soc_row = 11)
        let n_rows = 17;
        let mut t_rows: Vec<usize> = Vec::new();
        let mut t_cols: Vec<usize> = Vec::new();
        let mut t_vals: Vec<f64> = Vec::new();
        let mut b = vec![0.0f64; n_rows];

        // mx − λ = 0  →  A[0,mx]=1, A[0,λ]=−1
        t_rows.push(0);
        t_cols.push(1);
        t_vals.push(1.0);
        t_rows.push(0);
        t_cols.push(0);
        t_vals.push(-1.0);
        // my = 0  →  A[1,my]=1
        t_rows.push(1);
        t_cols.push(2);
        t_vals.push(1.0);
        // mxy = 0  →  A[2,mxy]=1
        t_rows.push(2);
        t_cols.push(3);
        t_vals.push(1.0);

        c.push_corner_blocks(
            &mut t_rows,
            &mut t_cols,
            &mut t_vals,
            &mut b,
            3,
            7,
            11,
            1,
            2,
            3,
            4,
            5,
            6,
            7,
        );

        let a_mat = CscMatrix::new_from_triplets(n_rows, n_vars, t_rows, t_cols, t_vals);
        let p_mat = CscMatrix::zeros((n_vars, n_vars));
        let mut q = vec![0.0f64; n_vars];
        q[0] = -1.0; // minimise −λ
        let cones: Vec<SupportedConeT<f64>> = vec![
            ZeroConeT(7),
            NonnegativeConeT(4),
            SecondOrderConeT(3),
            SecondOrderConeT(3),
        ];
        let settings = DefaultSettings {
            verbose: false,
            ..DefaultSettings::default()
        };
        let mut solver = DefaultSolver::new(&p_mat, &q, &a_mat, &b, &cones, settings);
        solver.solve();

        let lambda = solver.solution.x[0];
        assert!(
            (lambda - 2.0).abs() < 1e-5,
            "λ_max expected 2.0, got {lambda:.6}"
        );
    }

    // ── Admissible cases ───────────────────────────────────────────────────

    #[test]
    fn isotropic_zero() {
        let c = JohansenCriterion::isotropic(1.0);
        assert!(c.is_admissible(0.0, 0.0, 0.0));
        assert_eq!(c.yield_utilisation(0.0, 0.0, 0.0), 0.0);
    }

    #[test]
    fn isotropic_on_pos_surface() {
        // mx = mp, my = 0, mxy = 0: pos = (1-1)(1-0) = 0 — on surface, admissible.
        let c = JohansenCriterion::isotropic(1.0);
        assert!(c.is_admissible(1.0, 0.0, 0.0));
        assert_eq!(c.yield_utilisation(1.0, 0.0, 0.0), 1.0);
    }

    #[test]
    fn isotropic_interior() {
        // pos = (0.5)(0.5) - 0.16 = 0.09 > 0, neg = (1.5)(1.5) - 0.16 = 2.09 > 0.
        let c = JohansenCriterion::isotropic(1.0);
        assert!(c.is_admissible(0.5, 0.5, 0.4));
    }

    #[test]
    fn ortho_sagging_interior() {
        // mpx_pos=1, mpy_pos=2. mx=0.5, my=1.5, mxy=0.3.
        // pos = (1-0.5)(2-1.5) - 0.09 = 0.25 - 0.09 = 0.16 > 0.
        // neg = (1+0.5)(2+1.5) - 0.09 = 5.25 - 0.09 > 0.
        let c = JohansenCriterion::orthotropic(1.0, 2.0);
        assert!(c.is_admissible(0.5, 1.5, 0.3));
    }

    #[test]
    fn asym_small_hogging() {
        // mpx_pos=1, mpy_pos=1, mpx_neg=0.5, mpy_neg=0.5.
        // mx=-0.4, my=-0.4, mxy=0.
        // neg = (0.5-0.4)(0.5-0.4) = 0.01 > 0.
        let c = JohansenCriterion::new(1.0, 1.0, 0.5, 0.5);
        assert!(c.is_admissible(-0.4, -0.4, 0.0));
    }

    // ── Violating cases ────────────────────────────────────────────────────

    #[test]
    fn isotropic_twist_violation() {
        // mx=0.5, my=0.5, mxy=0.6: pos = 0.25 - 0.36 = -0.11 < 0.
        let c = JohansenCriterion::isotropic(1.0);
        assert!(!c.is_admissible(0.5, 0.5, 0.6));
        assert!(c.yield_utilisation(0.5, 0.5, 0.6) > 1.0);
    }

    #[test]
    fn ortho_exceed_mpx_pos() {
        // mpx_pos=1, mpy_pos=2. mx=1.5 > mpx_pos.
        // pos = (1-1.5)(2-0) = -1.0 < 0.
        let c = JohansenCriterion::orthotropic(1.0, 2.0);
        assert!(!c.is_admissible(1.5, 0.0, 0.0));
    }

    #[test]
    fn asym_exceed_mpx_neg() {
        // mpx_neg=0.5, mpy_neg=0.5. mx=-0.6, mxy=0.3.
        // neg = (0.5-0.6)(0.5+0) - 0.09 = (-0.1)(0.5) - 0.09 = -0.05 - 0.09 < 0.
        let c = JohansenCriterion::new(1.0, 1.0, 0.5, 0.5);
        assert!(!c.is_admissible(-0.6, 0.0, 0.3));
    }
}
