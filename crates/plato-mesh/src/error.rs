use thiserror::Error;

#[derive(Debug, Error)]
pub enum MeshError {
    #[error("Panel '{0}' not found")]
    UnknownPanel(String),
    #[error("Panel polygon is invalid: {0}")]
    InvalidPolygon(String),
    #[error("Hole polygon is outside or overlaps panel outline")]
    InvalidHole,
    #[error("Shared edge mismatch: '{panel_a}' edge {edge_a} vs '{panel_b}' edge {edge_b}")]
    SharedEdgeMismatch {
        panel_a: String,
        edge_a: usize,
        panel_b: String,
        edge_b: usize,
    },
    #[error("Triangulation failed: {0}")]
    TriangulationFailed(String),
}
