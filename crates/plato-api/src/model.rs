use serde::{Deserialize, Serialize};

/// Top-level analysis input. Fully serialisable.
/// Full definition in Phase 4.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnalysisModel;
