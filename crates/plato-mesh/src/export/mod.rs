use std::fmt::Write as FmtWrite;
use std::io::Write;
use std::path::Path;

use crate::model::MeshModel;

/// Write a flat SVG of the triangulated mesh for visual validation.
///
/// Draws wireframe triangles (corner nodes) with small dots at corner nodes
/// and midpoint nodes. No external crate is used — the SVG is built as a
/// plain string.
pub fn write_svg(mesh: &MeshModel, path: &Path) -> std::io::Result<()> {
    let svg = render_svg(mesh);
    let mut file = std::fs::File::create(path)?;
    file.write_all(svg.as_bytes())
}

fn render_svg(mesh: &MeshModel) -> String {
    const W: f64 = 800.0;
    const H: f64 = 600.0;
    const MARGIN: f64 = 0.05; // 5 % padding

    // Bounding box of all nodes.
    let (mut xmin, mut xmax, mut ymin, mut ymax) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for &[x, y] in &mesh.nodes {
        xmin = xmin.min(x);
        xmax = xmax.max(x);
        ymin = ymin.min(y);
        ymax = ymax.max(y);
    }
    if xmin == xmax {
        xmax = xmin + 1.0;
    }
    if ymin == ymax {
        ymax = ymin + 1.0;
    }

    let pad_x = (xmax - xmin) * MARGIN;
    let pad_y = (ymax - ymin) * MARGIN;
    xmin -= pad_x;
    xmax += pad_x;
    ymin -= pad_y;
    ymax += pad_y;

    let scale_x = W / (xmax - xmin);
    let scale_y = H / (ymax - ymin);
    let scale = scale_x.min(scale_y);

    // Map model coordinates to SVG coordinates (Y-axis flipped).
    let to_svg = |x: f64, y: f64| -> (f64, f64) { ((x - xmin) * scale, H - (y - ymin) * scale) };

    let mut buf = String::new();
    let _ = writeln!(
        buf,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{W}\" height=\"{H}\" \
         viewBox=\"0 0 {W} {H}\">"
    );
    let _ = writeln!(buf, "<rect width=\"{W}\" height=\"{H}\" fill=\"#1a1c22\"/>");

    // ── Wireframe triangles ───────────────────────────────────────────────
    buf.push_str("<g stroke=\"#4a9eff\" stroke-width=\"0.8\" fill=\"#1e2a3a\" fill-opacity=\"0.5\">\n");
    for el in &mesh.elements {
        let pts: String = el
            .corners
            .iter()
            .map(|&id| {
                let [x, y] = mesh.nodes[id];
                let (sx, sy) = to_svg(x, y);
                format!("{sx:.2},{sy:.2}")
            })
            .collect::<Vec<_>>()
            .join(" ");
        let _ = writeln!(buf, r#"  <polygon points="{pts}"/>"#);
    }
    let _ = writeln!(buf, "</g>");

    // ── Corner nodes ─────────────────────────────────────────────────────
    let mut corner_ids: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let mut midpoint_ids: std::collections::HashSet<usize> = std::collections::HashSet::new();
    for el in &mesh.elements {
        for &id in &el.corners {
            corner_ids.insert(id);
        }
        for &id in &el.midpoints {
            midpoint_ids.insert(id);
        }
    }

    buf.push_str("<g fill=\"#e8eaf6\">\n");
    for id in &corner_ids {
        let [x, y] = mesh.nodes[*id];
        let (sx, sy) = to_svg(x, y);
        let _ = writeln!(buf, r#"  <circle cx="{sx:.2}" cy="{sy:.2}" r="2.0"/>"#);
    }
    let _ = writeln!(buf, "</g>");

    // ── Midpoint nodes ────────────────────────────────────────────────────
    buf.push_str("<g fill=\"#ff8a65\">\n");
    for id in midpoint_ids.difference(&corner_ids) {
        let [x, y] = mesh.nodes[*id];
        let (sx, sy) = to_svg(x, y);
        let _ = writeln!(buf, r#"  <circle cx="{sx:.2}" cy="{sy:.2}" r="1.2"/>"#);
    }
    let _ = writeln!(buf, "</g>");

    // ── Quality legend ────────────────────────────────────────────────────
    let legend = format!(
        "  <text x=\"10\" y=\"20\" font-family=\"monospace\" font-size=\"12\" fill=\"#aaaaaa\">\
         {} elements  {} nodes  min angle {:.1} deg  max AR {:.2}</text>",
        mesh.quality.n_elements,
        mesh.quality.n_nodes,
        mesh.quality.min_angle_deg,
        mesh.quality.max_aspect_ratio,
    );
    let _ = writeln!(buf, "{legend}");

    let _ = writeln!(buf, "</svg>");
    buf
}
