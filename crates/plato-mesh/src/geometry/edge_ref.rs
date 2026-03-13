use serde::{Deserialize, Serialize};

/// A reference to one or more edges of a panel, used for applying boundary
/// conditions, loads, and resolving edge node lists.
///
/// Edge index `i` runs from `vertices[i]` to `vertices[(i+1) % n]`.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub enum EdgeRef {
    /// A single edge of the panel outline.
    PolygonEdge { panel_id: String, edge_index: usize },
    /// A single edge of a hole within the panel.
    HoleEdge {
        panel_id: String,
        hole_index: usize,
        edge_index: usize,
    },
    /// All outline edges of a panel (convenience shorthand).
    AllEdges { panel_id: String },
    /// All edges of a specific hole (convenience shorthand).
    AllHoleEdges { panel_id: String, hole_index: usize },
}
