use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::geometry::EdgeRef;

// ── Mesh configuration ────────────────────────────────────────────────────

/// Controls mesh density and quality for the Spade mesher.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MeshConfig {
    /// Maximum triangle area [length²]. Controls global mesh density.
    pub max_element_area: f64,
    /// Minimum interior angle [degrees] for Ruppert refinement.
    /// Default: 20.0. Values above 33.0 risk non-termination.
    pub min_angle_deg: f64,
    /// Local refinement zones (e.g. around re-entrant corners or openings).
    pub refinement_zones: Vec<RefinementZone>,
}

impl Default for MeshConfig {
    fn default() -> Self {
        Self {
            max_element_area: 0.05,
            min_angle_deg: 20.0,
            refinement_zones: vec![],
        }
    }
}

/// A circular region in which mesh density is locally increased.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RefinementZone {
    pub centre: [f64; 2],
    pub radius: f64,
    /// Must be smaller than `MeshConfig::max_element_area`.
    pub max_area: f64,
}

// ── MeshModel ────────────────────────────────────────────────────────────

/// Output of the mesher. Self-contained: carries node coordinates, element
/// connectivity, edge-to-node maps, and quality metrics.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MeshModel {
    /// Node coordinates `[x, y]`, indexed by global `node_id`.
    pub nodes: Vec<[f64; 2]>,
    /// All 6-node triangular elements.
    pub elements: Vec<TriElement6>,
    /// Ordered node IDs along each named edge (for BC and load application).
    /// `AllEdges` / `AllHoleEdges` keys are expanded at mesh time.
    pub edge_nodes: HashMap<EdgeRef, Vec<usize>>,
    /// Mesh quality diagnostics.
    pub quality: MeshQuality,
}

// ── TriElement6 ──────────────────────────────────────────────────────────

/// A 6-node triangular plate bending element.
///
/// Node ordering follows the design document convention:
/// - `corners[i]` for i in 0..3: corner nodes in CCW order
/// - `midpoints[i]`: midpoint node between `corners[i]` and `corners[(i+1) % 3]`
///
/// ```text
///   corners[0] --mid[0]-- corners[1]
///        \                    /
///       mid[2]            mid[1]
///           \              /
///           corners[2]
/// ```
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TriElement6 {
    pub id: usize,
    /// Corner node indices into `MeshModel::nodes`, in CCW order.
    pub corners: [usize; 3],
    /// Midpoint node indices.
    /// `midpoints[i]` is the midpoint of the side between `corners[i]` and `corners[(i+1)%3]`.
    pub midpoints: [usize; 3],
    /// ID of the panel this element belongs to.
    pub panel_id: String,
}

// ── MeshQuality ──────────────────────────────────────────────────────────

/// Summary quality metrics for the triangulation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MeshQuality {
    pub n_elements: usize,
    pub n_nodes: usize,
    pub min_angle_deg: f64,
    pub max_aspect_ratio: f64,
}

// ── Tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn minimal_mesh() -> MeshModel {
        MeshModel {
            nodes: vec![[0.0, 0.0], [1.0, 0.0], [0.5, 1.0], [0.5, 0.0]],
            elements: vec![TriElement6 {
                id: 0,
                corners: [0, 1, 2],
                midpoints: [3, 3, 3], // dummy — just for roundtrip
                panel_id: "p".into(),
            }],
            edge_nodes: HashMap::new(),
            quality: MeshQuality {
                n_elements: 1,
                n_nodes: 4,
                min_angle_deg: 45.0,
                max_aspect_ratio: 1.5,
            },
        }
    }

    #[test]
    fn mesh_serde_roundtrip() {
        let m = minimal_mesh();
        let json = serde_json::to_string(&m).expect("serialise");
        let m2: MeshModel = serde_json::from_str(&json).expect("deserialise");
        assert_eq!(m.nodes, m2.nodes);
        assert_eq!(m.elements.len(), m2.elements.len());
        assert_eq!(m.quality.n_elements, m2.quality.n_elements);
    }
}
