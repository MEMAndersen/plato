/// Geometric quantities for a 6-node triangular plate element.
///
/// Follows design-doc Appendix A.1 convention: **side `i` is opposite corner `i`**,
/// so side 0 spans corners 1→2, side 1 spans corners 2→0, side 2 spans corners 0→1.
///
/// All vectors are 3-component: the P(n)^T operator maps a 2-vector `v` to
/// `[nx·vx,  ny·vy,  ny·vx + nx·vy]` so that `a`, `b`, `c` align with the
/// moment column triplets `[mx, my, mxy]` in B^T.
#[derive(Debug, Clone, Copy)]
pub struct ElementGeometry {
    /// Element area (positive).
    pub area: f64,
    /// Side lengths: `l[i]` = length of side opposite corner `i`.
    pub l: [f64; 3],
    /// Outward unit normals: `n[i]` points away from corner `i`.
    pub n: [[f64; 2]; 3],
    /// Tangent vectors: `nhat[i]` = 90° CCW rotation of `n[i]`.
    pub nhat: [[f64; 2]; 3],
    /// `a[i][j] = P(n[i])^T · n[j]`  — 3-vector.
    pub a: [[[f64; 3]; 3]; 3],
    /// `b[i] = P(n[i])^T · nhat[i]`  — 3-vector.
    pub b: [[f64; 3]; 3],
    /// `c[i][j] = l[i] · l[j] · a[i][j]`  — 3-vector.
    pub c: [[[f64; 3]; 3]; 3],
}

impl ElementGeometry {
    /// Compute all geometric quantities from the three CCW corner coordinates.
    pub fn from_corners(p: &[[f64; 2]; 3]) -> Self {
        // ── Area ─────────────────────────────────────────────────────────
        let d10 = [p[1][0] - p[0][0], p[1][1] - p[0][1]];
        let d20 = [p[2][0] - p[0][0], p[2][1] - p[0][1]];
        let area = 0.5 * (d10[0] * d20[1] - d10[1] * d20[0]).abs();

        // ── Side lengths: side i spans corners (i+1)%3 → (i+2)%3 ────────
        let l = std::array::from_fn::<f64, 3, _>(|i| {
            let j = (i + 1) % 3;
            let k = (i + 2) % 3;
            let dx = p[k][0] - p[j][0];
            let dy = p[k][1] - p[j][1];
            (dx * dx + dy * dy).sqrt()
        });

        // ── Outward unit normals (design-doc exact formulas) ─────────────
        let n = [
            [(p[2][1] - p[1][1]) / l[0], -(p[2][0] - p[1][0]) / l[0]],
            [(p[0][1] - p[2][1]) / l[1], -(p[0][0] - p[2][0]) / l[1]],
            [(p[1][1] - p[0][1]) / l[2], -(p[1][0] - p[0][0]) / l[2]],
        ];

        // ── Tangents: 90° CCW from normal ────────────────────────────────
        let nhat = std::array::from_fn::<[f64; 2], 3, _>(|i| [-n[i][1], n[i][0]]);

        // ── P(n)^T applied to a 2-vector v → 3-vector ────────────────────
        // P(n)^T = [[nx, 0], [0, ny], [ny, nx]]
        let pt_apply = |ni: &[f64; 2], v: &[f64; 2]| -> [f64; 3] {
            [ni[0] * v[0], ni[1] * v[1], ni[1] * v[0] + ni[0] * v[1]]
        };

        // ── a[i][j] and b[i] ─────────────────────────────────────────────
        let a = std::array::from_fn::<[[f64; 3]; 3], 3, _>(|i| {
            std::array::from_fn::<[f64; 3], 3, _>(|j| pt_apply(&n[i], &n[j]))
        });
        let b = std::array::from_fn::<[f64; 3], 3, _>(|i| pt_apply(&n[i], &nhat[i]));

        // ── c[i][j] = l[i] * l[j] * a[i][j] ─────────────────────────────
        let c = std::array::from_fn::<[[f64; 3]; 3], 3, _>(|i| {
            std::array::from_fn::<[f64; 3], 3, _>(|j| {
                let s = l[i] * l[j];
                [s * a[i][j][0], s * a[i][j][1], s * a[i][j][2]]
            })
        });

        Self {
            area,
            l,
            n,
            nhat,
            a,
            b,
            c,
        }
    }
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

    #[test]
    fn equilateral_triangle_geom() {
        let s3 = f64::sqrt(3.0);
        let g = ElementGeometry::from_corners(&equilateral_corners());

        assert_abs_diff_eq!(g.area, s3 / 4.0, epsilon = 1e-12);
        for i in 0..3 {
            assert_abs_diff_eq!(g.l[i], 1.0, epsilon = 1e-12);
        }
        // n[0] = [s3/2, 1/2]
        assert_abs_diff_eq!(g.n[0][0], s3 / 2.0, epsilon = 1e-12);
        assert_abs_diff_eq!(g.n[0][1], 0.5, epsilon = 1e-12);
        // n[1] = [-s3/2, 1/2]
        assert_abs_diff_eq!(g.n[1][0], -s3 / 2.0, epsilon = 1e-12);
        assert_abs_diff_eq!(g.n[1][1], 0.5, epsilon = 1e-12);
        // n[2] = [0, -1]
        assert_abs_diff_eq!(g.n[2][0], 0.0, epsilon = 1e-12);
        assert_abs_diff_eq!(g.n[2][1], -1.0, epsilon = 1e-12);

        // b[0]: nhat[0]=[-1/2, s3/2], PT(n[0]) @ nhat[0]
        // = [(s3/2)(-1/2), (1/2)(s3/2), (1/2)(-1/2)+(s3/2)(s3/2)] = [-s3/4, s3/4, 1/2]
        assert_abs_diff_eq!(g.b[0][0], -s3 / 4.0, epsilon = 1e-12);
        assert_abs_diff_eq!(g.b[0][1], s3 / 4.0, epsilon = 1e-12);
        assert_abs_diff_eq!(g.b[0][2], 0.5, epsilon = 1e-12);

        // b[2]: n[2]=[0,-1], nhat[2]=[1,0] → [0, 0, -1]
        assert_abs_diff_eq!(g.b[2][0], 0.0, epsilon = 1e-12);
        assert_abs_diff_eq!(g.b[2][1], 0.0, epsilon = 1e-12);
        assert_abs_diff_eq!(g.b[2][2], -1.0, epsilon = 1e-12);

        // a[0][0] = PT(n[0]) @ n[0] = [3/4, 1/4, s3/2]
        assert_abs_diff_eq!(g.a[0][0][0], 3.0 / 4.0, epsilon = 1e-12);
        assert_abs_diff_eq!(g.a[0][0][1], 1.0 / 4.0, epsilon = 1e-12);
        assert_abs_diff_eq!(g.a[0][0][2], s3 / 2.0, epsilon = 1e-12);
    }
}
