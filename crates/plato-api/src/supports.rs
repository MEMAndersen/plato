use plato_mesh::geometry::EdgeRef;
use serde::{Deserialize, Serialize};

/// Boundary condition applied to an edge of a panel.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Support {
    /// w = 0 along the edge; normal rotation free.
    SimplySupportedEdge { edge: EdgeRef },
    /// w = 0 and normal rotation = 0 along the edge.
    ClampedEdge { edge: EdgeRef },
    /// Natural BC — no constraint. Explicitly stated for documentation clarity.
    FreeEdge { edge: EdgeRef },
    /// Symmetry line: mxy = 0 and odd-symmetric displacements removed.
    ///
    /// **Not yet implemented.** Returns [`AnalysisError::Validation`] if used.
    SymmetryLine { edge: EdgeRef, axis: SymmetryAxis },
}

/// The axis of symmetry for a [`Support::SymmetryLine`].
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum SymmetryAxis {
    X,
    Y,
}
