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

# run tests for a specific crate (cleaner output — less Clarabel noise)
cargo test -p plato-core   # 22 unit tests: elements, yield criterion, DOF map, assembly
cargo test -p plato-mesh   # mesh generation tests
cargo test -p plato-api    # end-to-end integration tests

# run a single test by name (partial match)
cargo test johansen

# suppress captured stdout (Clarabel solver prints to stderr; that always shows)
cargo test -q
```

> **Note:** The Clarabel solver prints diagnostic output to stderr during every solve,
> which makes the workspace test output noisy. Run per-crate to get a clean view.

> **Known failures (Phase 4 WIP):** Two integration tests in `plato-api` currently fail —
> `simply_supported_square` and `simply_supported_beam_8m`. This is a known equilibrium
> formulation issue under active investigation. All unit tests in `plato-core` pass.

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
