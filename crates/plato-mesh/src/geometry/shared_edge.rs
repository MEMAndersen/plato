use serde::{Deserialize, Serialize};

/// Declares that two panel outline edges are geometrically coincident
/// and should share nodes in the mesh.
///
/// The two edges must have opposite orientation (one runs A→B, the other B→A).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SharedEdgeDeclaration {
    pub panel_a: String,
    /// Zero-based edge index in `panel_a`'s outline.
    pub edge_a: usize,
    pub panel_b: String,
    /// Zero-based edge index in `panel_b`'s outline (opposite orientation to `edge_a`).
    pub edge_b: usize,
}

impl SharedEdgeDeclaration {
    pub fn between(
        panel_a: impl Into<String>,
        edge_a: usize,
        panel_b: impl Into<String>,
        edge_b: usize,
    ) -> Self {
        Self {
            panel_a: panel_a.into(),
            edge_a,
            panel_b: panel_b.into(),
            edge_b,
        }
    }
}
