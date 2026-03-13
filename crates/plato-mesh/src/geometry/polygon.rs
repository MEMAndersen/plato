use serde::{Deserialize, Serialize};

/// An ordered polygon defined by its vertices in counter-clockwise order.
///
/// Edge `i` runs from `vertices[i]` to `vertices[(i+1) % n]`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Polygon2D {
    /// Vertices in CCW order; minimum 3.
    pub vertices: Vec<[f64; 2]>,
}

impl Polygon2D {
    /// Construct from a vertex list. Panics if fewer than 3 vertices are given
    /// (invalid geometry is a programmer error at model-build time).
    pub fn new(vertices: Vec<[f64; 2]>) -> Self {
        assert!(
            vertices.len() >= 3,
            "Polygon2D requires at least 3 vertices"
        );
        Self { vertices }
    }

    /// Convenience constructor for an axis-aligned rectangle in CCW order:
    /// `(x_min,y_min) → (x_max,y_min) → (x_max,y_max) → (x_min,y_max)`.
    pub fn rectangle(x_min: f64, y_min: f64, x_max: f64, y_max: f64) -> Self {
        Self::new(vec![
            [x_min, y_min],
            [x_max, y_min],
            [x_max, y_max],
            [x_min, y_max],
        ])
    }

    /// Number of edges, which equals the number of vertices.
    pub fn n_edges(&self) -> usize {
        self.vertices.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn n_edges_equals_vertex_count() {
        let p = Polygon2D::new(vec![[0.0, 0.0], [1.0, 0.0], [0.5, 1.0]]);
        assert_eq!(p.n_edges(), 3);

        let q = Polygon2D::new(vec![[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [0.0, 1.0]]);
        assert_eq!(q.n_edges(), 4);
    }

    #[test]
    fn rectangle_has_four_ccw_vertices() {
        let r = Polygon2D::rectangle(0.0, 0.0, 3.0, 4.0);
        assert_eq!(r.vertices.len(), 4);
        // Check CCW winding via signed area (should be positive).
        let signed_area = {
            let v = &r.vertices;
            let n = v.len();
            (0..n)
                .map(|i| {
                    let j = (i + 1) % n;
                    v[i][0] * v[j][1] - v[j][0] * v[i][1]
                })
                .sum::<f64>()
                / 2.0
        };
        assert!(
            signed_area > 0.0,
            "rectangle should be CCW (positive signed area)"
        );
    }

    #[test]
    #[should_panic]
    fn panics_on_two_vertices() {
        Polygon2D::new(vec![[0.0, 0.0], [1.0, 0.0]]);
    }
}
