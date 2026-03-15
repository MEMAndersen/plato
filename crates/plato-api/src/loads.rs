use plato_mesh::geometry::EdgeRef;
use serde::{Deserialize, Serialize};

/// A load applied to the slab.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Load {
    /// Distributed load over a panel area.
    AreaLoad {
        panel_id: String,
        intensity: LoadIntensity,
        load_case: LoadCase,
    },
    /// Line load along a panel edge.
    ///
    /// `intensity` is force per unit length in the user's unit system.
    LineLoad {
        edge: EdgeRef,
        intensity: f64,
        load_case: LoadCase,
    },
}

/// Spatial distribution of an area load.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum LoadIntensity {
    /// Uniform pressure [force / area] over the entire panel.
    Uniform(f64),
}

/// Distinguishes permanent (dead) loads from variable (live) loads.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum LoadCase {
    /// Factored into the variable load vector q_ext; multiplied by λ.
    Variable,
    /// Factored into the permanent load vector q_dead; always present.
    Permanent,
}

/// Controls mesh density when building an `AnalysisModel` via `ModelBuilder`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum MeshDensity {
    /// Directly specify the maximum element area.
    MaxElementArea(f64),
    /// Target approximately `n` elements across the shortest bounding-box dimension.
    ElementsAcross(usize),
}
