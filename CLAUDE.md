# PLATO — Claude Code Context

PLATO (Plate Limit Analysis TOol) computes the plastic collapse load of reinforced concrete plates as a single convex SOCP — no incremental nonlinear FEM.

For full detail, see the design documents:
- [Design_documents/limit_analysis_engine_design.md](Design_documents/limit_analysis_engine_design.md) — solver theory, architecture, API, mesh generation, SOCP assembly
- [Design_documents/gui_design_document.md](Design_documents/gui_design_document.md) — GUI layout, rendering, state design, theming, implementation phases

---

## Workspace Structure

```
plato/
├── Cargo.toml
├── crates/
│   ├── plato-core/   # Elements, yield criterion, SOCP assembly
│   ├── plato-mesh/   # Spade-based Delaunay mesh generation
│   ├── plato-api/    # Public API: run_analysis(), ProgressCallback, AnalysisModel
│   ├── plato-cli/    # CLI frontend
│   └── plato-gui/    # egui/eframe/wgpu native GUI (planned)
└── examples/
```

## Key Technology Choices

| Concern | Choice |
|---------|--------|
| Language | Rust 2021 |
| Solver | Clarabel (native Rust SOCP) |
| Mesh | Spade (constrained Delaunay + Ruppert refinement) |
| GUI framework | egui + eframe |
| GPU rendering | wgpu via `egui_wgpu::Callback` |
| Plotting | egui_plot |
| Serialisation | serde_json (`.plato` project files) |

## Core Concepts

- **λ (load factor)** — the scalar maximised by the solver; the safe collapse load multiplier
- **AnalysisModel** — top-level serialisable input (panels, materials, supports, loads)
- **MeshModel** — output of `plato-mesh`; consumed by both solver and renderer
- **SolveResult** — contains λ, per-element moment fields, and the dual (collapse mode) vector
- **ProgressCallback** — trait injected into `run_analysis_with_progress()` to stream solver events; GUI uses mpsc channels

## GUI Architecture (plato-gui)

The GUI design document describes the planned `plato-gui` crate in detail. Key points:
- `PlatoApp` owns `ModelState`, `SolverState`, `Viewport`, `UiState`
- Solver runs on a background `std::thread`; progress polled each egui frame via `mpsc`
- Three wgpu render pipelines: face fill (scalar field), wireframe edges, instanced glyphs
- Theme constants in `theme.rs`; colourmap is `RdBu` symmetric about zero (red = sagging, blue = hogging)
- Geometry editing is via the properties panel only — no direct viewport manipulation

## Conventions

- Moment colourmap: red = positive (sagging), blue = negative (hogging), white = zero
- Colours carry semantic meaning: green = correct/success, amber = loads/warnings, red = errors, purple = BCs
- Shared edges shown green (geometrically matched) or red (mismatched) in the viewport
- Internal storage always in base SI units; `UnitSystem` controls input/output conversion
