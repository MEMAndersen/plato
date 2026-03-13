use crate::error::MeshError;
use crate::geometry::{Panel, SharedEdgeDeclaration};
use crate::model::{MeshConfig, MeshModel};

/// Backend-agnostic meshing interface.
pub trait Mesher: Send + Sync {
    fn triangulate(
        &self,
        panels: &[Panel],
        shared_edges: &[SharedEdgeDeclaration],
        config: &MeshConfig,
    ) -> Result<MeshModel, MeshError>;
}
