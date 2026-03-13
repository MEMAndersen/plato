use serde::{Deserialize, Serialize};

use super::polygon::Polygon2D;

/// A slab panel defined by its geometry only (no material — material is added in `plato-api`).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Panel {
    /// Unique identifier.
    pub id: String,
    /// Outer boundary in CCW order.
    pub outline: Polygon2D,
    /// Interior voids. Each hole is a CCW polygon.
    /// `hole_index` is the position in this `Vec`.
    pub holes: Vec<Polygon2D>,
}

impl Panel {
    pub fn new(id: impl Into<String>, outline: Polygon2D) -> Self {
        Self {
            id: id.into(),
            outline,
            holes: vec![],
        }
    }
}
