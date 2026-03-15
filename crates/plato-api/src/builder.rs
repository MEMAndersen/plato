use plato_mesh::geometry::{Polygon2D, SharedEdgeDeclaration};
use plato_mesh::model::MeshConfig;

use crate::loads::{Load, MeshDensity};
use crate::model::{AnalysisError, AnalysisModel, Panel, RCSlabMaterial, SolveResult, SolverConfig};
use crate::run::{ProgressCallback, run_analysis_with_progress};
use crate::supports::Support;
use crate::units::UnitSystem;

/// Fluent builder for [`AnalysisModel`].
///
/// All methods take `self` by value and return a new `ModelBuilder`,
/// allowing method chaining without `mut`.
pub struct ModelBuilder {
    name: String,
    units: UnitSystem,
    panels: Vec<Panel>,
    shared_edges: Vec<SharedEdgeDeclaration>,
    supports: Vec<Support>,
    loads: Vec<Load>,
    mesh_config: MeshConfig,
    solver_config: SolverConfig,
}

impl ModelBuilder {
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            units: UnitSystem::SI,
            panels: Vec::new(),
            shared_edges: Vec::new(),
            supports: Vec::new(),
            loads: Vec::new(),
            mesh_config: MeshConfig::default(),
            solver_config: SolverConfig::default(),
        }
    }

    pub fn with_units(mut self, units: UnitSystem) -> Self {
        self.units = units;
        self
    }

    pub fn set_mesh_config(mut self, config: MeshConfig) -> Self {
        self.mesh_config = config;
        self
    }

    /// Set mesh density from a [`MeshDensity`] hint.
    ///
    /// `ElementsAcross(n)` computes `max_element_area` from the bounding box
    /// of all panel outlines: `(shortest_dimension / n)²`.
    pub fn set_mesh_density(mut self, density: MeshDensity) -> Self {
        match density {
            MeshDensity::MaxElementArea(area) => {
                self.mesh_config.max_element_area = area;
            }
            MeshDensity::ElementsAcross(n) => {
                let area = self.compute_target_area(n);
                self.mesh_config.max_element_area = area;
            }
        }
        self
    }

    pub fn set_solver_config(mut self, config: SolverConfig) -> Self {
        self.solver_config = config;
        self
    }

    /// Add a panel with a placeholder zero-capacity material.
    ///
    /// Call [`set_material`](Self::set_material) afterwards to assign the real material.
    pub fn add_panel(mut self, id: impl Into<String>, outline: Polygon2D) -> Self {
        let id = id.into();
        self.panels.push(Panel {
            id,
            outline,
            holes: vec![],
            material: RCSlabMaterial::isotropic(0.0, 0.0),
        });
        self
    }

    /// Add a hole to an existing panel.
    ///
    /// # Panics
    /// Panics if `panel_id` is not found (programmer error).
    pub fn add_hole(mut self, panel_id: impl Into<String>, hole: Polygon2D) -> Self {
        let panel_id = panel_id.into();
        let panel = self
            .panels
            .iter_mut()
            .find(|p| p.id == panel_id)
            .unwrap_or_else(|| panic!("panel '{panel_id}' not found"));
        panel.holes.push(hole);
        self
    }

    pub fn declare_shared_edge(mut self, decl: SharedEdgeDeclaration) -> Self {
        self.shared_edges.push(decl);
        self
    }

    /// Set the material for a panel.
    ///
    /// # Panics
    /// Panics if `panel_id` is not found (programmer error).
    pub fn set_material(mut self, panel_id: impl Into<String>, mat: RCSlabMaterial) -> Self {
        let panel_id = panel_id.into();
        let panel = self
            .panels
            .iter_mut()
            .find(|p| p.id == panel_id)
            .unwrap_or_else(|| panic!("panel '{panel_id}' not found"));
        panel.material = mat;
        self
    }

    pub fn add_support(mut self, support: Support) -> Self {
        self.supports.push(support);
        self
    }

    pub fn add_load(mut self, load: Load) -> Self {
        self.loads.push(load);
        self
    }

    /// Validate and build the [`AnalysisModel`].
    pub fn build(self) -> Result<AnalysisModel, AnalysisError> {
        validate_builder(&self)?;
        Ok(AnalysisModel {
            name: self.name,
            units: self.units,
            panels: self.panels,
            shared_edges: self.shared_edges,
            supports: self.supports,
            loads: self.loads,
            mesh_config: self.mesh_config,
            solver_config: self.solver_config,
        })
    }

    /// Build, mesh, and solve without progress reporting.
    pub fn solve(self) -> Result<SolveResult, AnalysisError> {
        let model = self.build()?;
        crate::run::run_analysis(&model)
    }

    /// Build, mesh, and solve with a progress callback.
    pub fn solve_with_progress(self, p: &dyn ProgressCallback) -> Result<SolveResult, AnalysisError> {
        let model = self.build()?;
        run_analysis_with_progress(&model, p)
    }

    // ── Internal helpers ─────────────────────────────────────────────────

    fn compute_target_area(&self, n: usize) -> f64 {
        if self.panels.is_empty() || n == 0 {
            return MeshConfig::default().max_element_area;
        }
        let mut x_min = f64::MAX;
        let mut x_max = f64::MIN;
        let mut y_min = f64::MAX;
        let mut y_max = f64::MIN;
        for panel in &self.panels {
            for &[x, y] in &panel.outline.vertices {
                x_min = x_min.min(x);
                x_max = x_max.max(x);
                y_min = y_min.min(y);
                y_max = y_max.max(y);
            }
        }
        let w = x_max - x_min;
        let h = y_max - y_min;
        let shortest = w.min(h);
        let cell = shortest / n as f64;
        cell * cell
    }
}

fn validate_builder(b: &ModelBuilder) -> Result<(), AnalysisError> {
    if b.panels.is_empty() {
        return Err(AnalysisError::Validation("no panels defined".into()));
    }
    let has_variable_load = b.loads.iter().any(|l| {
        matches!(
            l,
            Load::AreaLoad { load_case: crate::loads::LoadCase::Variable, .. }
                | Load::LineLoad { load_case: crate::loads::LoadCase::Variable, .. }
        )
    });
    if !has_variable_load {
        return Err(AnalysisError::Validation(
            "at least one Variable load is required".into(),
        ));
    }
    Ok(())
}
