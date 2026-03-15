use plato_mesh::geometry::{Polygon2D, SharedEdgeDeclaration};
use plato_mesh::model::MeshConfig;
use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::loads::Load;
use crate::supports::Support;
use crate::units::UnitSystem;

// ── Material ──────────────────────────────────────────────────────────────

/// Reinforced concrete slab yield moments per unit width.
///
/// Positive moments produce sagging (tension at bottom).
/// All values are in the user's [`UnitSystem`] and converted to SI internally.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RCSlabMaterial {
    /// Sagging (positive) yield moment in the x-direction [moment/length].
    pub m_pos_x: f64,
    /// Sagging (positive) yield moment in the y-direction [moment/length].
    pub m_pos_y: f64,
    /// Hogging (negative) yield moment in the x-direction [moment/length].
    pub m_neg_x: f64,
    /// Hogging (negative) yield moment in the y-direction [moment/length].
    pub m_neg_y: f64,
}

impl RCSlabMaterial {
    /// Isotropic slab: same positive and negative capacity in both directions.
    pub fn isotropic(m_pos: f64, m_neg: f64) -> Self {
        Self {
            m_pos_x: m_pos,
            m_pos_y: m_pos,
            m_neg_x: m_neg,
            m_neg_y: m_neg,
        }
    }
}

// ── Panel (API) ───────────────────────────────────────────────────────────

/// A slab panel with geometry and material.
///
/// Wraps [`plato_mesh::geometry::Panel`] with an associated [`RCSlabMaterial`].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Panel {
    pub id: String,
    pub outline: Polygon2D,
    pub holes: Vec<Polygon2D>,
    pub material: RCSlabMaterial,
}

impl Panel {
    pub fn new(id: impl Into<String>, outline: Polygon2D, material: RCSlabMaterial) -> Self {
        Self {
            id: id.into(),
            outline,
            holes: vec![],
            material,
        }
    }

    /// Convert to the mesh-layer panel type (geometry only, no material).
    pub(crate) fn to_mesh_panel(&self) -> plato_mesh::geometry::Panel {
        plato_mesh::geometry::Panel {
            id: self.id.clone(),
            outline: self.outline.clone(),
            holes: self.holes.clone(),
        }
    }
}

// ── SolverConfig ──────────────────────────────────────────────────────────

/// Solver tolerances and limits passed to Clarabel.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SolverConfig {
    pub tol_gap_abs: f64,
    pub tol_gap_rel: f64,
    pub max_iter: u32,
}

impl Default for SolverConfig {
    fn default() -> Self {
        Self {
            tol_gap_abs: 1e-8,
            tol_gap_rel: 1e-8,
            max_iter: 200,
        }
    }
}

// ── AnalysisModel ─────────────────────────────────────────────────────────

/// Top-level serialisable input for a plastic limit-analysis run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnalysisModel {
    pub name: String,
    pub units: UnitSystem,
    pub panels: Vec<Panel>,
    pub shared_edges: Vec<SharedEdgeDeclaration>,
    pub supports: Vec<Support>,
    pub loads: Vec<Load>,
    pub mesh_config: MeshConfig,
    pub solver_config: SolverConfig,
}

// ── SolveResult ───────────────────────────────────────────────────────────

/// Status of a completed analysis run.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SolveStatus {
    /// Converged to the requested tolerance.
    Optimal,
    /// Terminated at iteration limit; result may be approximate.
    MaxIterationsReached,
    /// The problem is infeasible (e.g. fully unconstrained structure).
    Infeasible,
    /// Numerical failure inside the solver.
    NumericalError,
}

/// One triangular element in the solve result (for self-contained rendering).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SolveElement {
    pub id: usize,
    pub panel_id: String,
    /// `[corner0, corner1, corner2, mid01, mid12, mid20]` — indices into `SolveResult::nodes`.
    pub nodes: [usize; 6],
}

/// Output of a successful `run_analysis()` call.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SolveResult {
    pub status: SolveStatus,
    /// Plastic collapse load factor λ.  NaN if not solved.
    pub load_factor: f64,
    pub units: UnitSystem,
    /// Wall-clock solve time in milliseconds.
    pub solve_time_ms: u64,
    /// Node coordinates copied from the mesh.
    pub nodes: Vec<[f64; 2]>,
    /// Element connectivity copied from the mesh.
    pub elements: Vec<SolveElement>,
    /// Per-element moment fields (empty if not solved).
    pub element_moments: Vec<plato_core::result::ElementMoments>,
    pub n_elements: usize,
    pub n_nodes: usize,
    pub n_free_dofs: usize,
    pub n_variables: usize,
    pub solver_iterations: usize,
    pub duality_gap: f64,
}

// ── Errors ────────────────────────────────────────────────────────────────

/// Error returned by `run_analysis()`.
#[derive(Debug, Error)]
pub enum AnalysisError {
    #[error("mesh generation failed: {0}")]
    Mesh(#[from] plato_mesh::error::MeshError),
    #[error("assembly error: {0}")]
    Assembly(#[from] AssemblyError),
    #[error("validation error: {0}")]
    Validation(String),
    #[error("analysis was cancelled")]
    Cancelled,
}

/// Errors arising during DOF map / matrix assembly.
#[derive(Debug, Error)]
pub enum AssemblyError {
    #[error("all DOFs are constrained — structure is fully fixed")]
    FullyConstrained,
    #[error("element {0} is degenerate (zero area)")]
    DegenerateElement(usize),
}
