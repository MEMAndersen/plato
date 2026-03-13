use crate::model::AnalysisModel;
use std::ops::ControlFlow;

/// Events streamed from the solver to the caller.
/// Full definition in Phase 5.
pub enum ProgressEvent {
    MeshingStarted,
    MeshingDone { n_elements: usize, n_nodes: usize },
    AssemblyStarted,
    AssemblyDone { n_dofs: usize, n_variables: usize },
    SolverIteration { iteration: usize, gap: f64 },
    Done { load_factor: f64 },
    Cancelled,
}

/// Injectable progress / cancellation interface.
pub trait ProgressCallback: Send {
    /// Return `ControlFlow::Break(())` to request cancellation.
    fn on_event(&self, event: ProgressEvent) -> ControlFlow<()>;
}

/// Placeholder solve result. Full definition in Phase 5.
pub struct SolveResult;

/// Error type for analysis failures.
#[derive(Debug)]
pub enum AnalysisError {
    Cancelled,
    Other(String),
}

/// Run analysis without progress reporting.
pub fn run_analysis(_model: &AnalysisModel) -> Result<SolveResult, AnalysisError> {
    unimplemented!("run_analysis will be implemented in Phase 3–5")
}

/// Run analysis with a progress callback.
pub fn run_analysis_with_progress(
    _model: &AnalysisModel,
    _progress: &dyn ProgressCallback,
) -> Result<SolveResult, AnalysisError> {
    unimplemented!("run_analysis_with_progress will be implemented in Phase 3–5")
}
