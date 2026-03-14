use nalgebra::SMatrix;
use plato_mesh::model::TriElement6;

use super::geometry::ElementGeometry;

/// A 6-node triangular plate bending element with pre-computed geometry.
pub struct PlateElement<'a> {
    pub elem: &'a TriElement6,
    pub geom: ElementGeometry,
}

impl<'a> PlateElement<'a> {
    pub fn new(elem: &'a TriElement6, corners: &[[f64; 2]; 3]) -> Self {
        Self {
            elem,
            geom: ElementGeometry::from_corners(corners),
        }
    }

    /// Returns the 12×9 local equilibrium matrix **B^T**.
    ///
    /// **Row layout:**
    /// - Rows 0-2: displacement DOFs at corners 0, 1, 2
    /// - Rows 3-5: displacement DOFs at midpoints 0, 1, 2
    /// - Rows 6-11: rotation DOFs θ0a, θ0b, θ1a, θ1b, θ2a, θ2b
    ///
    /// **Column layout:** block `k` (cols `3k..3k+2`) = `[mx, my, mxy]` at corner `k`.
    ///
    /// B^T = B^T_r + B^T_q + B^T_t + B^T_θ  (design-doc Appendix A.2)
    pub fn local_bt(&self) -> SMatrix<f64, 12, 9> {
        let g = &self.geom;
        let mut bt = SMatrix::<f64, 12, 9>::zeros();

        // Add ±(s * v) into row r, moment block blk (cols 3*blk..3*blk+2).
        macro_rules! add {
            ($r:expr, $blk:expr, $v:expr, $s:expr) => {
                bt[($r, 3 * $blk)] += $s * $v[0];
                bt[($r, 3 * $blk + 1)] += $s * $v[1];
                bt[($r, 3 * $blk + 2)] += $s * $v[2];
            };
        }
        macro_rules! sub {
            ($r:expr, $blk:expr, $v:expr, $s:expr) => {
                bt[($r, 3 * $blk)] -= $s * $v[0];
                bt[($r, 3 * $blk + 1)] -= $s * $v[1];
                bt[($r, 3 * $blk + 2)] -= $s * $v[2];
            };
        }

        let b = &g.b;
        let a = &g.a;
        let c = &g.c;
        let inv12a = 1.0 / (12.0 * g.area);
        let inv6 = 1.0 / 6.0;

        // ── B^T_r: concentrated corner forces ─────────────────────────────
        for i in 0..3 {
            let j = (i + 1) % 3;
            let k = (i + 2) % 3;
            add!(i, i, sub3(&b[j], &b[k]), 1.0);
        }

        // ── B^T_q: Kirchhoff shear forces, prefactor 1/(12A) ─────────────
        for i in 0..3 {
            for j in 0..3 {
                sub!(i, j, c[j][i], inv12a);
                add!(3 + i, j, c[j][i], 4.0 * inv12a);
            }
        }

        // ── B^T_t: twisting moment gradient, prefactor 1/6 ───────────────
        // Corner rows
        add!(0, 0, sub3(&b[2], &b[1]), inv6);
        sub!(0, 1, b[2], inv6);
        add!(0, 2, b[1], inv6);

        add!(1, 0, b[2], inv6);
        add!(1, 1, sub3(&b[0], &b[2]), inv6);
        sub!(1, 2, b[0], inv6);

        sub!(2, 0, b[1], inv6);
        add!(2, 1, b[0], inv6);
        add!(2, 2, sub3(&b[1], &b[0]), inv6);

        // Midpoint rows
        add!(3, 1, b[0], 4.0 * inv6);
        sub!(3, 2, b[0], 4.0 * inv6);

        sub!(4, 0, b[1], 4.0 * inv6);
        add!(4, 2, b[1], 4.0 * inv6);

        add!(5, 0, b[2], 4.0 * inv6);
        sub!(5, 1, b[2], 4.0 * inv6);

        // ── B^T_θ: moment continuity at rotation nodes ────────────────────
        add!(6, 1, a[0][0], 1.0);
        sub!(7, 2, a[0][0], 1.0);

        add!(8, 2, a[1][1], 1.0);
        sub!(9, 0, a[1][1], 1.0);

        add!(10, 0, a[2][2], 1.0);
        sub!(11, 1, a[2][2], 1.0);

        bt
    }
}

// ── Private helpers ───────────────────────────────────────────────────────

#[inline]
fn sub3(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

// ── Tests ─────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    fn equilateral_corners() -> [[f64; 2]; 3] {
        let s3 = f64::sqrt(3.0);
        [[0.0, 0.0], [1.0, 0.0], [0.5, s3 / 2.0]]
    }

    fn dummy_elem() -> TriElement6 {
        TriElement6 {
            id: 0,
            corners: [0, 1, 2],
            midpoints: [3, 4, 5],
            panel_id: "p".into(),
        }
    }

    #[test]
    fn bt_rigid_body_null_uniform_w() {
        // Patch test: for uniform vertical displacement (w=1 at all 6 displacement
        // nodes), virtual work must vanish for all moment fields m.
        // Equivalently: column sums of B^T rows 0-5 must all be zero.
        let test_cases: &[[[f64; 2]; 3]] = &[
            [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
            [[0.0, 0.0], [3.0, 0.0], [1.0, 2.0]],
            equilateral_corners(),
        ];
        for corners in test_cases {
            let elem = dummy_elem();
            let plate = PlateElement::new(&elem, corners);
            let bt = plate.local_bt();
            for col in 0..9 {
                let col_sum: f64 = (0..6).map(|r| bt[(r, col)]).sum();
                assert!(
                    col_sum.abs() < 1e-10,
                    "col {col} sum = {col_sum:.2e} for corners {corners:?}"
                );
            }
        }
    }

    #[test]
    fn equilateral_bt_rotation_rows() {
        let s3 = f64::sqrt(3.0);
        let elem = dummy_elem();
        let plate = PlateElement::new(&elem, &equilateral_corners());
        let bt = plate.local_bt();

        // a[0][0] = [3/4, 1/4, s3/2]
        // Row 6 (θ0a): blk1 = +a[0][0], rest zero
        assert_abs_diff_eq!(bt[(6, 3)], 3.0 / 4.0, epsilon = 1e-10);
        assert_abs_diff_eq!(bt[(6, 4)], 1.0 / 4.0, epsilon = 1e-10);
        assert_abs_diff_eq!(bt[(6, 5)], s3 / 2.0, epsilon = 1e-10);
        for k in [0usize, 1, 2, 6, 7, 8] {
            assert_abs_diff_eq!(bt[(6, k)], 0.0, epsilon = 1e-10);
        }

        // Row 7 (θ0b): blk2 = -a[0][0], rest zero
        assert_abs_diff_eq!(bt[(7, 6)], -3.0 / 4.0, epsilon = 1e-10);
        assert_abs_diff_eq!(bt[(7, 7)], -1.0 / 4.0, epsilon = 1e-10);
        assert_abs_diff_eq!(bt[(7, 8)], -s3 / 2.0, epsilon = 1e-10);
        for k in [0usize, 1, 2, 3, 4, 5] {
            assert_abs_diff_eq!(bt[(7, k)], 0.0, epsilon = 1e-10);
        }

        // a[1][1]: n[1]=[-s3/2, 1/2], PT(n[1])@n[1] = [3/4, 1/4, -s3/2]
        //          (cross term: (1/2)(-s3/2) + (-s3/2)(1/2) = -s3/2)
        // Row 9 (θ1b): blk0 = -a[1][1] = [-3/4, -1/4, +s3/2]
        assert_abs_diff_eq!(bt[(9, 0)], -3.0 / 4.0, epsilon = 1e-10);
        assert_abs_diff_eq!(bt[(9, 1)], -1.0 / 4.0, epsilon = 1e-10);
        assert_abs_diff_eq!(bt[(9, 2)], s3 / 2.0, epsilon = 1e-10);

        // a[2][2]: n[2]=[0,-1], PT(n[2])@n[2] = [0, 1, 0]
        //          (cross term: (-1)(0) + (0)(-1) = 0)
        // Row 10 (θ2a): blk0 = +a[2][2] = [0, 1, 0]
        assert_abs_diff_eq!(bt[(10, 0)], 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(bt[(10, 1)], 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(bt[(10, 2)], 0.0, epsilon = 1e-10);
    }
}
