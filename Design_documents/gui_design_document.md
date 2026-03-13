# GUI Design Document — PLATO (Plate Limit Analysis TOol)

**Version:** 0.2-draft  
**Status:** Proposal  
**Language:** Rust  
**Scope:** Native desktop GUI crate (`plato-gui`) integrated into the existing `limit-analysis` workspace

---

## 1. Overview

This document describes the architecture, library choices, layout, and implementation plan for the `plato-gui` crate. The GUI provides pre-processing (model definition and mesh preview), solver control with live progress, and post-processing (moment field contours and collapse mode visualisation) for the RC plate limit analysis engine.

The GUI crate is added directly to the existing `limit-analysis` Cargo workspace and links against `la-api` in-process — no IPC, no serialisation round-trip at runtime. Communication between GUI and engine happens via direct Rust function calls and the `ProgressCallback` trait already defined in `la-api`.

---

## 2. Workspace Integration

The new crate slots into the existing workspace without modifying any existing crates:

```
limit-analysis/
├── Cargo.toml                  # workspace — add "crates/plato-gui" to members
├── crates/
│   ├── la-core/                # unchanged
│   ├── la-mesh/                # unchanged
│   ├── la-api/                 # unchanged — plato-gui depends on this
│   ├── la-cli/                 # unchanged
│   └── plato-gui/                 # NEW
│       ├── Cargo.toml
│       └── src/
│           ├── main.rs         # eframe::run_native() entry point
│           ├── app.rs          # LimitAnalysisApp — top-level eframe::App impl
│           ├── model_state.rs  # Owned AnalysisModel + dirty tracking + undo stack
│           ├── solver_state.rs # Solver thread handle, ProgressEvent receiver
│           ├── panels/
│           │   ├── mod.rs
│           │   ├── model_panel.rs    # Left panel: geometry tree + material inputs
│           │   ├── viewer_panel.rs   # Central panel: 3D viewport
│           │   └── results_panel.rs  # Bottom panel: charts + scalar selector
│           ├── renderer/
│           │   ├── mod.rs
│           │   ├── mesh_renderer.rs  # wgpu mesh wireframe + face fill
│           │   ├── field_renderer.rs # Scalar field colouring (moments, collapse)
│           │   └── glyph_renderer.rs # BC and load glyphs (arrows, pins)
│           └── dialogs/
│               ├── material_dialog.rs
│               ├── mesh_settings_dialog.rs
│               └── export_dialog.rs
```

Dependency in `crates/plato-gui/Cargo.toml`:

```toml
[dependencies]
la-api  = { path = "../la-api" }
la-mesh = { path = "../la-mesh" }   # for MeshModel types in renderer
```

---

## 3. Technology Stack

### 3.1 GUI Framework — egui + eframe

**Choice: egui with the eframe native shell**

`egui` is an immediate-mode GUI library written in pure Rust. `eframe` is the official wrapper that handles the native OS window, event loop, and wgpu backend. This is the right choice for an engineering analysis tool:

- Immediate-mode maps naturally to analysis tools where UI state is always derived from the underlying model and result data — no complex retained widget trees to keep in sync
- Pure Rust, no C++ runtime, no system dependencies beyond the GPU driver
- First-class `wgpu` integration via `egui_wgpu`: custom 3D rendering in the same window as egui widgets, sharing the same GPU context
- Used in production engineering and scientific tools (Rerun, various CAD tooling)

```toml
egui      = "0.29"
eframe    = { version = "0.29", features = ["wgpu"] }
egui_wgpu = "0.29"
```

**Why not the alternatives:**

| Option | Verdict |
|--------|---------|
| **egui** | Immediate-mode, pure Rust, wgpu integration, active ecosystem. **Selected.** |
| Tauri | Web frontend (HTML/JS) — not Rust GUI |
| Slint | Declarative `.slint` language, good for traditional apps but less suited to engineering tools |
| iced | Elm-architecture, less mature for complex layouts and custom 3D |
| Bevy UI | Full game engine overhead, embedding it just for UI is heavy |

### 3.2 3D Rendering — wgpu via egui_wgpu

**Choice: wgpu with a custom renderer embedded via `egui_wgpu::Callback`**

`wgpu` is the standard Rust GPU API (WebGPU-based), portable across Vulkan, Metal, DX12, and WebGL. The 3D renderer is written directly against wgpu rather than pulling in a full scene-graph library. The mesh structures are simple flat triangles with per-face scalar values — a full 3D engine would be overkill and would add unnecessary build complexity.

The viewport is implemented as an `egui_wgpu::Callback`, which lets a custom wgpu render pass draw into an `egui::Rect` region of the window, composited seamlessly with the egui panels.

**What the custom wgpu renderer handles:**
- Flat-shaded triangulated mesh faces coloured by a per-face scalar field via colourmap texture
- Mesh wireframe overlay (edge lines)
- Load glyphs (instanced arrow meshes, scaled by magnitude)
- BC glyphs (pin/roller/clamp symbols as instanced meshes)
- Collapse mode: vertex positions displaced by the dual solution vector, CPU-side, scaled by a user-controlled warp factor

### 3.3 Plotting — egui_plot

Used for the live solver convergence chart (iteration vs. duality gap, log scale) and a utilisation histogram after solve. Integrates directly into the egui layout without a separate widget system.

```toml
egui_plot = "0.29"
```

### 3.4 Concurrency — std::thread + std::sync::mpsc

The solver runs on a background `std::thread`. Progress events stream back to the GUI via `std::sync::mpsc` channels, polled each frame in egui's `update()` loop. This maps directly onto the `ProgressCallback` trait already in `la-api` — no async runtime required.

### 3.5 Serialisation

Project save/load uses `serde_json` on `AnalysisModel`, which is already fully `Serialize + Deserialize` from the engine design. No additional serialisation work is needed in `plato-gui`.

### 3.6 Full Dependency Summary

```toml
[dependencies]
# Engine
la-api       = { path = "../la-api" }
la-mesh      = { path = "../la-mesh" }

# GUI
egui         = "0.29"
eframe       = { version = "0.29", features = ["wgpu"] }
egui_wgpu    = "0.29"
egui_plot    = "0.29"

# GPU
wgpu         = "22"
bytemuck     = { version = "1", features = ["derive"] }  # vertex buffer casting

# Serialisation
serde_json   = "1"
toml         = "0.8"

# Utilities
log          = "0.4"
env_logger   = "0.11"
rfd          = "0.14"   # native OS file open/save dialogs
```

---

## 4. Application State Design

The top-level `LimitAnalysisApp` struct owns all application state. In egui's immediate-mode model, this struct is mutated by UI code each frame and its fields drive what gets rendered.

```rust
pub struct LimitAnalysisApp {
    model_state:  ModelState,
    solver_state: SolverState,
    viewport:     Viewport,
    ui_state:     UiState,
}

pub struct ModelState {
    model:        AnalysisModel,       // live model being edited
    mesh_preview: Option<MeshModel>,   // available after meshing, before solve
    is_dirty:     bool,                // unsaved changes flag
    undo_stack:   Vec<AnalysisModel>,  // serialised snapshots
    redo_stack:   Vec<AnalysisModel>,
}

pub struct SolverState {
    status:        SolverStatus,
    thread_handle: Option<std::thread::JoinHandle<()>>,
    progress_rx:   Option<std::sync::mpsc::Receiver<ProgressEvent>>,
    result:        Option<SolveResult>,
    convergence:   Vec<(usize, f64)>,  // (iteration, duality_gap) for chart
}

pub enum SolverStatus {
    Idle,
    Meshing,
    Solving { iteration: usize, gap: f64 },
    Done,
    Failed(String),
}

pub struct UiState {
    selected_panel_id: Option<String>,
    scalar_field:      ScalarField,
    warp_factor:       f32,
    show_mesh_edges:   bool,
    show_load_glyphs:  bool,
    show_bc_glyphs:    bool,
    colormap:          Colormap,
}

pub enum ScalarField {
    None,
    MomentMx,
    MomentMy,
    MomentMxy,
    PrincipalM1,
    PrincipalM2,
    YieldUtilisation,
    CollapseDisplacement,
}
```

---

## 5. Layout

```
┌─────────────────────────────────────────────────────────────────────┐
│  TopPanel: Menu bar + Toolbar                                       │
│  [File] [Model] [Analysis] [View]   [New][Open][Save] | [▶ Run]    │
├──────────────┬──────────────────────────────────┬───────────────────┤
│              │                                  │                   │
│  SidePanel   │      CentralPanel                │  SidePanel        │
│  (Left)      │      3D Viewport                 │  (Right)          │
│              │      (wgpu via egui_wgpu)        │                   │
│  Model Tree  │                                  │  Properties       │
│  ─────────── │                                  │  ─────────────── │
│  📐 Panels   │                                  │  Selected item    │
│    span_1    │                                  │  attributes       │
│    span_2    │                                  │                   │
│  🧱 Materials│                                  │  Material:        │
│  🔗 Supports │                                  │  m⁺ₓ  [10.0 kNm] │
│  ⬇ Loads    │                                  │  m⁺ᵧ  [10.0 kNm] │
│  📊 Results  │                                  │  m⁻ₓ  [12.0 kNm] │
│    λ = 1.42  │                                  │  m⁻ᵧ  [12.0 kNm] │
│              │                                  │                   │
├──────────────┴──────────────────────────────────┴───────────────────┤
│  BottomPanel: Scalar selector + convergence chart + summary         │
│  Field: [m_x ▼]  Warp: [──●──]  |  [convergence chart]  | λ=1.423 │
├─────────────────────────────────────────────────────────────────────┤
│  StatusBar: [✓ Mesh: 3,564 el]  [Solver: 47 iter, gap 1×10⁻⁸]     │
└─────────────────────────────────────────────────────────────────────┘
```

All panels are resizable via egui splitters. The central viewport fills the remaining space.

### 5.1 Toolbar

```
[ New ] [ Open ] [ Save ]  ─  [ Preview Mesh ] [ ▶ Run ] [ ■ Stop ]  ─  [ Top ] [ ISO ] [ Reset ]
```

- **Preview Mesh** — calls `la-mesh` directly, displays `MeshModel` without solving
- **▶ Run** — spawns solver thread, transitions to `SolverStatus::Solving`
- **■ Stop** — sends cancellation signal to solver thread (see Open Questions §10.3)

### 5.2 Model Tree (Left Panel)

Built with `egui::CollapsingHeader`. Clicking a leaf node selects it (highlighted in the 3D viewport) and populates the right properties panel.

```
▼ Panels
    span_1
    span_2
▼ Materials
    RC Slab — m⁺=10 / m⁻=12 kNm/m
▼ Supports
    SS — span_1 edge 0 (south)
    SS — span_1 edge 2 (north)
▼ Loads
    UDL 1.0 kN/m² — span_1
    UDL 1.0 kN/m² — span_2
▼ Results  [visible only after successful solve]
    Load factor λ = 1.423
    3,564 elements · 1,847 nodes
    47 iterations · gap = 9.8×10⁻⁹
```

### 5.3 Viewport (Central Panel)

Camera controls follow standard CAD conventions:
- Left drag → orbit
- Scroll → zoom
- Middle drag / Shift+left drag → pan
- `R` → reset to isometric view
- `T` / `F` / `S` → top / front / side view

**Pre-solve rendering:** mesh faces filled with a flat colour per panel, wireframe overlay, load and BC glyphs.

**Post-solve rendering:** mesh faces coloured by the selected scalar field via colourmap texture, optional wireframe, optional warp by collapse displacement. Colourbar drawn as an egui widget overlaid on the viewport.

### 5.4 Properties Panel (Right Panel)

**Panel selected:**
```
Panel ID:   span_1
Vertices:   4
Area:       12.00 m²

Material
  m⁺ₓ  [10.00] kNm/m
  m⁺ᵧ  [10.00] kNm/m
  m⁻ₓ  [12.00] kNm/m
  m⁻ᵧ  [12.00] kNm/m

Mesh Density
  Max element area:  [0.05] m²
```

**Support selected:**
```
Type:   Simply Supported Edge
Panel:  span_1 / edge 0
        (0.0, 0.0) → (3.0, 0.0)
```

**Load selected:**
```
Type:       Area Load
Panel:      span_1
Intensity:  Uniform  [1.00] kN/m²
Load case:  Variable
```

### 5.5 Results / Charts Panel (Bottom Panel)

```
Display field: [ m_x ▼ ]   Warp (collapse): [────●─] ×12.4
Colormap:      [ RdBu  ▼ ] Edges: [✓]  Glyphs: [✓]
```

Followed by the `egui_plot` convergence chart (duality gap vs. iteration, log Y axis) updating live during solve. A summary table shows λ, element count, iteration count, gap, and wall-clock solve time.

---

## 6. Solver Integration

```rust
// In plato-gui/src/solver_state.rs
use std::sync::mpsc;
use la_api::{run_analysis, ProgressCallback, ProgressEvent, AnalysisModel};

struct ChannelProgress(mpsc::Sender<ProgressEvent>);

impl ProgressCallback for ChannelProgress {
    fn on_event(&self, event: ProgressEvent) {
        let _ = self.0.send(event);
    }
}

impl SolverState {
    pub fn launch(&mut self, model: AnalysisModel) {
        let (tx, rx) = mpsc::channel();
        self.progress_rx = Some(rx);
        self.status = SolverStatus::Meshing;

        self.thread_handle = Some(std::thread::spawn(move || {
            let cb = ChannelProgress(tx);
            // run_analysis is the existing la-api entry point
            let _ = run_analysis(&model, cb);
        }));
    }

    /// Called every frame from egui update() to drain the channel.
    pub fn poll(&mut self) {
        if let Some(rx) = &self.progress_rx {
            while let Ok(event) = rx.try_recv() {
                match event {
                    ProgressEvent::SolverIteration { iteration, gap } => {
                        self.status = SolverStatus::Solving { iteration, gap };
                        self.convergence.push((iteration, gap));
                    }
                    ProgressEvent::Done => { self.status = SolverStatus::Done; }
                    _ => {}
                }
            }
        }
    }
}
```

The `SolveResult` is sent back to the main thread via a second `mpsc` channel (or wrapped in `Arc<Mutex<Option<SolveResult>>>`), then stored in `SolverState::result` and uploaded to the GPU renderer.

---

## 7. wgpu Renderer Design

Three `wgpu::RenderPipeline` instances, each compiled once at startup:

**Face pipeline** — fills triangles. Per-vertex attributes: `position: vec3<f32>`, `scalar: f32`. The fragment shader samples a 256×1 RGBA colourmap texture using the scalar value as UV coordinate. For the collapse mode, vertex positions are computed CPU-side as `base_pos + warp_factor * displacement` and re-uploaded via `queue.write_buffer()` when the warp slider changes.

**Edge pipeline** — `LineList` topology. Per-vertex attributes: `position: vec3<f32>`, `color: vec4<f32>`. Used for wireframe overlay and boundary condition glyph outlines.

**Glyph pipeline** — instanced rendering of procedurally generated arrow and pin meshes. Per-instance attributes: `position: vec3<f32>`, `orientation: vec4<f32>` (quaternion), `scale: f32`, `color: vec4<f32>`.

Data is re-uploaded to GPU buffers only when the `MeshModel` or `SolveResult` changes (not every frame), keeping per-frame cost to a uniform + draw call.

---

## 8. File I/O

| Format | Direction | Notes |
|--------|-----------|-------|
| `.plato` (JSON) | Read/Write | `serde_json` on `AnalysisModel` — native project format |
| `.vtk` (ASCII legacy) | Write | Nodes + cells + scalar fields for ParaView |
| `.csv` | Write | Element ID, m_x, m_y, m_xy, utilisation per row |
| `.png` | Write | wgpu framebuffer capture via `Texture::copy_to_buffer` |

Native OS file dialogs via the `rfd` crate (pure Rust, no GTK/Qt dependency).

---

## 9. Visual Design and Theming

This section defines the egui `Style`, `Visuals`, and `Spacing` configuration for the application, and the colour conventions used in the wgpu renderer. The goal is a consistent, professional appearance suited to a structural engineering tool — dark-themed to reduce eye strain during long analysis sessions, with colour used purposefully to convey meaning rather than decoration.

### 9.1 Colour Palette

All colours are defined as named constants in `plato-gui/src/theme.rs` and referenced throughout both the egui panels and the wgpu renderer. Using a single source of truth ensures the UI widgets and the 3D viewport always agree on what a colour means.

```rust
// plato-gui/src/theme.rs
use egui::Color32;

// --- Base UI Colours ---
pub const BACKGROUND:        Color32 = Color32::from_rgb(28,  30,  36);  // near-black, main window bg
pub const PANEL_BG:          Color32 = Color32::from_rgb(36,  38,  46);  // slightly lighter, side panels
pub const WIDGET_BG:         Color32 = Color32::from_rgb(48,  51,  62);  // inputs, buttons
pub const BORDER:            Color32 = Color32::from_rgb(65,  68,  82);  // panel separators
pub const TEXT_PRIMARY:      Color32 = Color32::from_rgb(220, 222, 228); // main text
pub const TEXT_SECONDARY:    Color32 = Color32::from_rgb(140, 145, 160); // labels, hints
pub const TEXT_DISABLED:     Color32 = Color32::from_rgb(90,  94, 110);  // disabled items

// --- Semantic / Status Colours ---
pub const ACCENT:            Color32 = Color32::from_rgb(82,  139, 255); // primary interactive (blue)
pub const SUCCESS:           Color32 = Color32::from_rgb(72,  199, 142); // shared edge OK, solve done
pub const WARNING:           Color32 = Color32::from_rgb(245, 166,  35); // solver stall, numerical warning
pub const ERROR:             Color32 = Color32::from_rgb(240,  80,  80); // shared edge mismatch, solve fail

// --- Viewport / Renderer Colours ---
pub const VIEWPORT_BG:       Color32 = Color32::from_rgb(22,  24,  30);  // 3D canvas background
pub const MESH_FACE_DEFAULT: Color32 = Color32::from_rgb(80,  100, 140); // unselected panel face fill
pub const MESH_FACE_SELECTED:Color32 = Color32::from_rgb(120, 160, 220); // selected panel highlight
pub const MESH_EDGE:         Color32 = Color32::from_rgb(50,  55,  70);  // wireframe edges
pub const SHARED_EDGE_OK:    Color32 = Color32::from_rgb(72,  199, 142); // = SUCCESS — matched shared edge
pub const SHARED_EDGE_BAD:   Color32 = Color32::from_rgb(240,  80,  80); // = ERROR  — mismatched shared edge
pub const LOAD_GLYPH:        Color32 = Color32::from_rgb(245, 166,  35); // = WARNING — load arrows
pub const BC_GLYPH:          Color32 = Color32::from_rgb(160, 100, 220); // support pins/rollers (purple)
pub const COLLAPSE_GLYPH:    Color32 = Color32::from_rgb(220, 220, 100); // collapse mode warp overlay
```

The semantic colours (`SUCCESS`, `WARNING`, `ERROR`) are reused in both the egui log panel text and the wgpu renderer glyphs, so the same green that labels a "Solve succeeded" log line is the same green on a correctly matched shared edge. This consistency reinforces meaning across the UI.

### 9.2 Spacing System

egui's `Spacing` struct is configured once at startup and applied globally. All layout decisions use multiples of a 4px base unit, matching standard engineering UI conventions.

```rust
// Applied in app.rs during setup
fn apply_spacing(ctx: &egui::Context) {
    let mut style = (*ctx.style()).clone();
    style.spacing.item_spacing     = egui::vec2(8.0, 6.0);   // 2×base, 1.5×base
    style.spacing.window_margin    = egui::Margin::same(12.0);// 3×base
    style.spacing.button_padding   = egui::vec2(12.0, 6.0);
    style.spacing.indent           = 16.0;                    // 4×base — tree indentation
    style.spacing.slider_width     = 120.0;
    style.spacing.combo_width      = 140.0;
    ctx.set_style(style);
}
```

### 9.3 Typography

egui uses a single built-in proportional font by default. We configure two sizes:

| Role | Size | Usage |
|------|------|-------|
| `Body` | 14px (default) | All panel text, property labels, tree items |
| `Small` | 11px | Status bar, secondary labels, colourbar tick marks |
| `Heading` | 16px, bold | Panel section headers (e.g. "Material", "Mesh Density") |
| `Monospace` | 13px | Solver log output, numerical values in the results summary |

Section headers inside panels use a subtle horizontal rule beneath them to create visual grouping without requiring nested widgets.

### 9.4 egui Visuals Configuration

```rust
fn apply_visuals(ctx: &egui::Context) {
    let mut visuals = egui::Visuals::dark(); // start from the dark base
    
    visuals.panel_fill              = PANEL_BG;
    visuals.window_fill             = PANEL_BG;
    visuals.extreme_bg_color        = BACKGROUND;
    visuals.faint_bg_color          = WIDGET_BG;
    visuals.code_bg_color           = BACKGROUND;

    visuals.widgets.inactive.bg_fill   = WIDGET_BG;
    visuals.widgets.inactive.fg_stroke = egui::Stroke::new(1.0, TEXT_SECONDARY);
    visuals.widgets.hovered.bg_fill    = ACCENT.linear_multiply(0.3);
    visuals.widgets.active.bg_fill     = ACCENT.linear_multiply(0.5);
    visuals.selection.bg_fill          = ACCENT.linear_multiply(0.25);
    visuals.selection.stroke           = egui::Stroke::new(1.0, ACCENT);

    visuals.window_rounding   = egui::Rounding::same(6.0);
    visuals.menu_rounding     = egui::Rounding::same(4.0);
    visuals.widgets.inactive.rounding = egui::Rounding::same(3.0);

    ctx.set_visuals(visuals);
}
```

### 9.5 Colourmap — Moment Fields

The `RdBu` colourmap is stored as a 256×1 RGBA texture uploaded to the GPU once at startup. The texture is sampled by the face shader using the normalised scalar value as the U coordinate.

The range is always **symmetric about zero**: given the min and max values of the displayed field, the scale is set to `[-S, +S]` where `S = max(|min|, |max|)`. This ensures zero always maps to white (the neutral midpoint of `RdBu`), negative moments are always blue, and positive moments are always red, regardless of the absolute magnitudes.

```
Colourmap scale:
  -S ──────────── 0 ──────────── +S
  Blue          White            Red
  (hogging)                  (sagging)
```

A colourbar widget is drawn as an egui `Painter` element overlaid in the bottom-right corner of the viewport, showing the numeric range at each end and a zero tick at the midpoint.

### 9.6 Visual Hierarchy Summary

The design deliberately uses colour sparingly so that colour always carries meaning:

| Colour | Meaning |
|--------|---------|
| Blue (accent) | Interactive elements, selection highlight |
| Green | Correct / matched / solver success |
| Amber | Loads, warnings, solver stall |
| Red | Errors, mismatched shared edges, solver failure |
| Purple | Boundary conditions (distinct from loads) |
| Red–White–Blue gradient | Moment field magnitude and sign |
| Dim blue-grey | Default mesh face fill (neutral, recedes) |

---

## 10. Implementation Phases

**Phase 1 — Window and layout skeleton**  
`eframe` window with all five panel slots, placeholder content, camera orbit/zoom/pan, wgpu pipeline initialised and drawing a test triangle.

**Phase 2 — Mesh preview**  
Load a hardcoded `AnalysisModel`, call `la-mesh` for `MeshModel`, upload to GPU, render wireframe + filled faces, clickable face highlights parent panel.

**Phase 3 — Model tree and properties editor**  
Model tree populated from `AnalysisModel`, right panel properties editing for all entity types, undo/redo stack via `AnalysisModel` snapshots, `is_dirty` flag driving the title bar.

**Phase 4 — Solver integration**  
Run button, solver thread, `ChannelProgress`, live convergence chart updating via mpsc poll, `SolveResult` stored on completion.

**Phase 5 — Results rendering**  
Scalar field upload to GPU, colourmap texture, colourbar egui widget, warp slider driving CPU-side displacement, scalar field selector populated from available result arrays.

**Phase 6 — File I/O**  
Save/load `.plato`, VTK export, CSV export, PNG screenshot.

**Phase 7 — Polish**  
BC and load glyphs, material dialog with schematic cross-section preview, keyboard shortcuts, status bar diagnostics, high-DPI awareness via eframe's `NativeOptions::pixels_per_point`, apply full theme from Section 9 (`theme.rs` constants, `Visuals`, `Spacing`).

---

## 11. Resolved Design Decisions

1. **Geometry editing** — panel vertices are edited via the properties panel only. No direct viewport manipulation.

2. **Shared edge visualisation** — declared `SharedEdge` connections are rendered as a coloured line drawn between the two shared edge midpoints (or along the shared edge itself). If the edge pairing is geometrically inconsistent (midpoints do not coincide within a tolerance), the line is drawn in red as a warning; a correctly matched shared edge is drawn in green. This makes misconfigured continuity immediately apparent before the user runs the solver.

3. **Cancellation token** — `ProgressCallback::on_event()` will return a `ControlFlow` enum (`Continue` / `Cancel`). The solver checks the return value at each iteration boundary. This addition to `la-api` is required before Phase 4 of the GUI implementation.

4. **Colourmap** — `RdBu` (red = positive/sagging, blue = negative/hogging) with the range automatically centred on zero and symmetric (`max(|min|, |max|)` on both sides). This matches standard structural engineering convention for bending moment diagrams.
