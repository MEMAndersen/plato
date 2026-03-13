/// Quick visual check: mesh an L-shaped panel and write mesh.svg next to this file.
///
/// Run with:
///   cargo run -p plato-mesh --example mesh_svg
///
/// Then open `crates/plato-mesh/examples/mesh.svg` in a browser.
use std::path::Path;

use plato_mesh::export::write_svg;
use plato_mesh::geometry::{Panel, Polygon2D};
use plato_mesh::mesher::{Mesher, SpadeMesher};
use plato_mesh::model::MeshConfig;

fn main() {
    // An L-shaped panel (3 m × 3 m with a 1.5 m × 1.5 m corner cut out).
    //
    //   (0,3)──────(1.5,3)
    //     │            │
    //   (0,1.5)─(1.5,1.5)
    //     │            │
    //   (0,0)──────(3,0)──(3,3)
    //
    // Vertices in CCW order:
    let l_shape = Polygon2D::new(vec![
        [0.0, 0.0],
        [3.0, 0.0],
        [3.0, 3.0],
        [1.5, 3.0],
        [1.5, 1.5],
        [0.0, 1.5],
    ]);

    let panels = vec![Panel::new("L", l_shape)];

    let config = MeshConfig {
        max_element_area: 0.05,
        min_angle_deg: 20.0,
        ..MeshConfig::default()
    };

    let mesh = SpadeMesher
        .triangulate(&panels, &[], &config)
        .expect("triangulation failed");

    println!(
        "Meshed {} elements, {} nodes  (min angle {:.1}°, max aspect {:.2})",
        mesh.quality.n_elements,
        mesh.quality.n_nodes,
        mesh.quality.min_angle_deg,
        mesh.quality.max_aspect_ratio,
    );

    let out = Path::new("crates/plato-mesh/examples/mesh.svg");
    write_svg(&mesh, out).expect("failed to write SVG");
    println!("SVG written to {}", out.display());
}
