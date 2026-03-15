# PLATO — Plate Limit Analysis TOol

PLATO computes the plastic collapse load of reinforced concrete plates as a single convex Second-Order Cone Program (SOCP). No incremental nonlinear FEM is required.

The solver is based on the **static (lower-bound) theorem of plasticity**: any stress field that satisfies equilibrium and does not violate the yield criterion is a safe lower bound to the true collapse load multiplier **λ**.

---

## Workspace layout

```
plato/
├── crates/
│   ├── plato-core/   # Elements, yield criterion, SOCP assembly
│   ├── plato-mesh/   # Spade-based Delaunay mesh generation
│   ├── plato-api/    # Public API: run_analysis(), AnalysisModel
│   └── plato-cli/    # CLI frontend (in development)
└── Design_documents/ # Solver theory, GUI design
```

---

## Prerequisites

- **Rust 1.85+** (2024 edition, set as MSRV in `clippy.toml`)

Install or update via [rustup](https://rustup.rs/):

```sh
rustup update stable
```

---

## Build

```sh
# debug
cargo build

# release (optimised — use this for actual analyses)
cargo build --release
```

---

## Tests

```sh
# run all tests across the workspace
cargo test

# run tests for a specific crate
cargo test -p plato-core
cargo test -p plato-mesh
cargo test -p plato-api

# run a single test by name (partial match)
cargo test johansen

# show stdout from tests
cargo test -- --nocapture
```

---

## Format

The workspace uses `rustfmt` with a 100-character line width (see `rustfmt.toml`).

```sh
# format all crates
cargo fmt

# check formatting without writing changes (useful in CI)
cargo fmt -- --check
```

---

## Lint

Clippy is configured with MSRV 1.85 (see `clippy.toml`).

```sh
# lint all crates
cargo clippy

# treat warnings as errors (matches CI strictness)
cargo clippy -- -D warnings
```

---

## Design documents

- [limit_analysis_engine_design.md](Design_documents/limit_analysis_engine_design.md) — solver theory, SOCP formulation, mesh generation, API
- [gui_design_document.md](Design_documents/gui_design_document.md) — GUI layout, rendering pipelines, theming

---

## Status

Alpha — active development. Breaking changes to data structures and APIs are expected.
