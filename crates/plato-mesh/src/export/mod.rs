use std::path::Path;

// Stub — will be implemented alongside MeshModel in Phase 1.
// Writes a flat SVG of the triangulated mesh for visual validation.
pub fn write_svg(_mesh: &crate::model::MeshModel, _path: &Path) -> std::io::Result<()> {
    unimplemented!("SVG export will be implemented in Phase 1 alongside MeshModel")
}
