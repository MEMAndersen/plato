# PLATO — Plate Limit Analysis TOol: Design & Implementation Specification

**Scope:** Numerical lower-bound plastic limit analysis of reinforced concrete plates (pure bending, no in-plane action)
**Language:** Rust
**Solver:** Clarabel (native Rust SOCP)
**Mesh:** Spade (pure Rust constrained Delaunay triangulation)
**Version:** v1.1 — Draft, March 2026

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Theoretical Background](#2-theoretical-background)
3. [System Architecture](#3-system-architecture)
4. [Loads and Boundary Conditions](#4-loads-and-boundary-conditions)
5. [API Design](#5-api-design)
6. [Mesh Generation](#6-mesh-generation)
7. [Problem Assembly and Optimisation](#7-problem-assembly-and-optimisation)
8. [Implementation Plan](#8-implementation-plan)
9. [Dependencies](#9-dependencies)
10. [Open Questions](#10-open-questions)
11. [Appendix A — Mathematical Reference](#appendix-a--mathematical-reference)

---

## 1. Executive Summary

PLATO (Plate Limit Analysis TOol) computes the maximum load factor **λ** that a reinforced concrete plate can sustain before plastic collapse, without performing incremental nonlinear analysis. The theoretical basis is the **static (lower-bound) theorem of plasticity**, which guarantees that any stress field satisfying equilibrium and not violating the yield criterion is a safe lower bound to the true collapse load.

The problem is solved once as a **convex optimisation problem (SOCP)**, making it orders of magnitude faster than nonlinear FEM for the purpose of design checking, capacity assessment, and reinforcement optimisation.

### Scope (v1)

| Item | Decision |
|---|---|
| Structural type | RC plates — pure bending only, no in-plane (membrane) action |
| Element type | Triangular plate bending element (Kirchhoff, linear moment field) |
| Yield criterion | Johansen's criterion, enforced at all 3 corner nodes per element |
| Reinforcement | Constant per slab region (orthotropic, different top/bottom) |
| Geometry input | Programmatic Rust API, fully serialisable; native file format `.plato` |
| Multi-panel support | Multiple slab regions sharing edges (continuous slabs) |
| Unit system | User-selectable via `UnitSystem` enum |
| Loads | One variable load pattern per analysis, plus optional permanent loads |
| Boundary conditions | DOF removal (not penalty); no point supports in v1 |
| Solver | Clarabel — SOCP via interior-point method |
| Mesh | Spade CDT with Ruppert refinement |

---

## 2. Theoretical Background

### 2.1 Lower-Bound Limit Analysis

The static theorem states: any moment field **m** that satisfies equilibrium everywhere and satisfies the yield criterion everywhere constitutes a safe lower bound to the true collapse load. The engine maximises the load multiplier λ over all such admissible moment fields:

```
maximise    λ

subject to  B^T · m = λ · q_ext + q_dead    (equilibrium at all free DOFs)
            f_k(m_e) ≤ 0  ∀ elements e, corners k    (Johansen criterion)
            m free
```

Where:
- `B^T` is the global equilibrium matrix assembled from element contributions
- `m` is the vector of moment unknowns (`mx, my, mxy` at each element corner node)
- `q_ext` is the pattern load vector (scaled by λ)
- `q_dead` is the permanent load vector (not scaled by λ)

### 2.2 The Plate Bending Element

A **triangular element with a linear moment field** is used. The moment field is defined by three moment vectors, one at each corner node, giving 9 moment unknowns per element. The linear field is complete within the element and enables exact enforcement of internal equilibrium.

#### Element Topology

Each triangular element has:
- **3 corner nodes** — transverse displacement DOF `w`
- **3 midpoint nodes** (one per side) — transverse displacement DOF `w`; used to distribute Kirchhoff shear forces via the quadratic displacement interpolation
- **6 rotation nodes** (two per side, near corners) — rotation DOF `θ`; used to enforce bending moment continuity across element sides

The moment field at any point in the element is interpolated linearly from the corner values:

```
m(ξ,η) = N1(ξ,η)·m1 + N2(ξ,η)·m2 + N3(ξ,η)·m3

where  m_i = [mix, miy, mixy]^T    (3 components per corner, 9 total per element)
       N_i  are the standard linear triangular shape functions (barycentric coordinates)
```

#### Equilibrium Matrix

The element equilibrium matrix `B^T` maps the 9 moment unknowns to nodal force/moment contributions at the 12 DOF nodes (6 displacement + 6 rotation). It decomposes into:

```
B^T_displacement = B^T_r + B^T_q + B^T_t     (6 displacement node rows)
B^T_θ                                          (6 rotation node rows, moment continuity)
```

Full explicit matrices are given in [Appendix A](#appendix-a--mathematical-reference).

**B^T_r — Concentrated corner forces**

At each corner, the difference in twisting moments along the two meeting sides produces a concentrated nodal force. The twisting moment along side i is `m^i_t = b^T_i · m`.

**B^T_q — Kirchhoff shear forces**

The shear force along side j from the moment gradient:
```
q_j = -(1/2A) Σ_i  l_i · a^T_ij · m_i
```
Distributed to the 6 displacement nodes via quadratic interpolation (corner nodes weighted −1×, midpoint nodes 4×, then scaled by 1/(12A)).

**B^T_t — Twisting moment gradient**

The gradient of the twisting moment along side i:
```
t_i = (b^T_i / l_i) · (m_k - m_j)    where (i,j,k) cyclic permutations of (1,2,3)
```

**B^T_θ — Moment continuity**

The normal bending moment at side i is `m^i_n = a^T_ii · m`. This is enforced at the two rotation nodes near each corner of each side, giving 6 equations that ensure the moment field is continuous across shared element edges.

### 2.3 Johansen's Yield Criterion

Johansen's criterion is the standard yield criterion for RC slabs. It accounts for different positive (sagging) and negative (hogging) capacities in x and y, and couples them via the twisting moment.

For a moment vector `m = [mx, my, mxy]^T` and material parameters:
- `m+_x`, `m+_y` — positive (bottom reinforcement) moment capacities
- `m-_x`, `m-_y` — negative (top reinforcement) moment capacities (positive values)

The criterion consists of two yield surfaces:

```
Positive yield surface (sagging):
    f+(m) = -(m+_x - mx)(m+_y - my) + m²_xy ≤ 0

Negative yield surface (hogging):
    f-(m) = -(m-_x + mx)(m-_y + my) + m²_xy ≤ 0
```

Both surfaces are enforced at each of the 3 corner nodes of every element.

#### Reformulation as Second-Order Cone (SOC) Constraints

Clarabel supports the standard Second-Order Cone:

```
(t, z) ∈ SOC(n+1)  ↔  ||z||₂ ≤ t,   t ≥ 0,   z ∈ ℝⁿ
```

The Johansen yield surfaces are quadratic products `2·u·v ≥ w²` (a Rotated SOC shape), which must be converted to standard SOC form. The standard conversion uses the substitution:

```
Given  2·u·v ≥ w²,  u ≥ 0, v ≥ 0:

Let  t = (u + v) / √2      (sum, scaled)
     s = (u − v) / √2      (difference, scaled)

Then  2·u·v = t² − s²  ≥  w²
  ↔   t² ≥ s² + w²
  ↔   ||(s, w)||₂ ≤ t

i.e.  [t; s; w] ∈ SOC(3)
```

**Positive yield surface** at corner node k: introduce primary auxiliaries

```
α1_k = (1/√2)(m+_x − mx_k) ≥ 0
α2_k = (1/√2)(m+_y − my_k) ≥ 0
```

so that `2·α1_k·α2_k ≥ mxy_k²`. Apply the substitution:

```
t+_k = (α1_k + α2_k) / √2      (SOC height, added as new auxiliary)
s+_k = (α1_k − α2_k) / √2      (SOC lateral, added as new auxiliary)

Equality enforcing the definition:
    t+_k = (α1_k + α2_k) / √2   →   √2·t+_k − α1_k − α2_k = 0
    s+_k = (α1_k − α2_k) / √2   →   √2·s+_k − α1_k + α2_k = 0

SOC constraint (Clarabel SecondOrderCone(3)):
    [t+_k;  s+_k;  mxy_k] ∈ SOC(3)    ↔    ||(s+_k, mxy_k)||₂ ≤ t+_k
```

**Negative yield surface** at corner node k: introduce

```
α3_k = (1/√2)(m-_x + mx_k) ≥ 0
α4_k = (1/√2)(m-_y + my_k) ≥ 0

t-_k = (α3_k + α4_k) / √2
s-_k = (α3_k − α4_k) / √2

Equality definitions:
    √2·t-_k − α3_k − α4_k = 0
    √2·s-_k − α3_k + α4_k = 0

SOC constraint:
    [t-_k;  s-_k;  mxy_k] ∈ SOC(3)    ↔    ||(s-_k, mxy_k)||₂ ≤ t-_k
```

#### Auxiliary Variable Count

Per corner node k there are **6 auxiliary variables**: `α1_k, α2_k, α3_k, α4_k` (primary, linking moments to the geometry of the yield surface) and `t+_k, t-_k, s+_k, s-_k`... 

Actually, noting that `s+_k` and `s-_k` can be treated as free variables determined entirely by the SOC constraint and the equality definitions, the cleanest implementation **eliminates `s` and `t` as explicit variables** by directly expressing the SOC in terms of `α`:

```
Since  t+_k = (α1_k + α2_k)/√2  and  s+_k = (α1_k − α2_k)/√2,
the SOC  ||(s+_k, mxy_k)||₂ ≤ t+_k  becomes a linear affine cone constraint:

    || [ (α1_k − α2_k)/√2 ] ||      (α1_k + α2_k)
    || [      mxy_k        ] ||₂  ≤  ─────────────
                                           √2

In Clarabel's Ax + b ∈ SOC form, this is expressed as three affine rows:
    row 0 (t):   (1/√2)·α1_k + (1/√2)·α2_k               →  t+_k
    row 1 (s):   (1/√2)·α1_k − (1/√2)·α2_k               →  s+_k
    row 2 (w):   mxy_k                                      →  mxy_k
```

This avoids introducing `t, s` as explicit decision variables. The SOC block is expressed as an **affine SOC constraint** `A_soc · x + b_soc ∈ SOC(3)` where `A_soc` selects and linearly combines the `α` and `mxy` entries from `x`. Clarabel accepts constraints in this form directly.

```
Per corner node k:  α1_k, α2_k, α3_k, α4_k          (4 auxiliary variables)
Per element:        3 corners × 4 = 12 auxiliary variables
Cone constraints:   3 corners × 2 SOC(3) blocks = 6 SOC(3) blocks per element
```

#### Orthotropic and Isotropic Cases

For **orthotropic** reinforcement: `m+_x ≠ m+_y` and/or `m-_x ≠ m-_y`. The criterion takes the same SOC form.

For **isotropic** reinforcement with equal capacity `m_p`:
```
m+_x = m+_y = m-_x = m-_y = m_p
```

---

## 3. System Architecture

### 3.1 Crate Structure

```
plato/
├── Cargo.toml                        # workspace manifest
├── crates/
│   ├── plato-core/                      # Pure numerics: elements, assembly, optimisation
│   │   └── src/
│   │       ├── element/
│   │       │   ├── mod.rs
│   │       │   ├── geometry.rs       # ElementGeometry: area, normals, side lengths, aux vectors
│   │       │   └── plate.rs          # PlateElement: local B^T assembly
│   │       ├── criteria/
│   │       │   └── johansen.rs       # Johansen criterion → SOC cone blocks
│   │       ├── assembly/
│   │       │   ├── dof_map.rs        # Node → global DOF index; BC DOF removal
│   │       │   └── global.rs         # Sparse B^T assembly (triplet → CSC)
│   │       ├── problem/
│   │       │   └── builder.rs        # ClarabelProblem: P, q, A, b, cones
│   │       └── result/
│   │           └── mod.rs            # SolveResult, collapse mode extraction
│   │
│   ├── plato-mesh/                      # Geometry and triangulation
│   │   └── src/
│   │       ├── geometry/
│   │       │   ├── polygon.rs        # Polygon2D: ordered vertex list, numbered edges
│   │       │   └── panel.rs          # Panel: polygon + holes + panel_id
│   │       ├── mesher/
│   │       │   ├── trait.rs          # Mesher trait (backend-agnostic)
│   │       │   └── spade.rs          # SpadeMesher: CDT + midpoint insertion
│   │       └── model/
│   │           └── mod.rs            # MeshModel: nodes, elements, shared edges
│   │
│   ├── plato-api/                       # Public-facing API — GUI-ready
│   │   └── src/
│   │       ├── model.rs              # AnalysisModel (fully serialisable)
│   │       ├── builder.rs            # ModelBuilder fluent API
│   │       ├── units.rs              # UnitSystem, unit conversion
│   │       ├── loads.rs              # Load types
│   │       ├── supports.rs           # Support/BC types
│   │       └── run.rs                # run_analysis() + ProgressCallback
│   │
│   └── plato-cli/                       # Optional CLI runner
│
└── examples/
    ├── simply_supported_slab.rs
    ├── clamped_slab.rs
    ├── two_span_slab.rs
    └── slab_with_opening.rs
```

### 3.2 Data Flow

```
ModelBuilder
    │  user defines: panels, shared edges, materials, loads, supports
    ▼
AnalysisModel  ──serialize──▶  .plato  (GUI save/load/undo)
    │
    ▼  plato-mesh
MeshModel
    │  nodes, 6-node elements, shared edge node sets
    ▼  plato-core / assembly
DofMap  +  GlobalEquilibriumMatrix  (sparse CSC B^T)
    │
    ▼  plato-core / problem
ClarabelProblem  (P=0, q, A_eq, b_eq, cones)
    │
    ▼  clarabel
ClarabelSolution  (primal x*, dual y*)
    │
    ▼  plato-core / result
SolveResult  (λ, moment field, collapse mode displacements)
    │
    ▼
User / GUI
```

### 3.3 GUI-Readiness Principles

- **All public types** in `plato-api` implement `serde::Serialize + serde::Deserialize`
- **No rendering code** anywhere in the engine crates
- **`MeshModel` is exposed** before solving so a GUI can show the mesh for inspection
- **`SolveResult` is self-contained** — a GUI renders moment contours and collapse mode directly from it
- **Progress reporting** via an injectable `ProgressCallback` trait

```rust
/// Run without progress reporting.
pub fn run_analysis(
    model: &AnalysisModel,
) -> Result<SolveResult, AnalysisError>;

/// Run with a progress callback (progress events + cancellation support).
pub fn run_analysis_with_progress(
    model: &AnalysisModel,
    progress: &dyn ProgressCallback,
) -> Result<SolveResult, AnalysisError>;

pub enum ProgressEvent {
    MeshingStarted,
    MeshingDone     { n_elements: usize, n_nodes: usize },
    AssemblyStarted,
    AssemblyDone    { n_dofs: usize, n_variables: usize },
    /// Fired at each solver iteration. Return `ControlFlow::Cancel` to abort.
    SolverIteration { iteration: usize, gap: f64 },
    /// Fired on successful completion with the final load factor.
    Done            { load_factor: f64 },
    /// Confirmation that a previously returned `ControlFlow::Cancel` was received
    /// and acted upon. Fired on the callback thread before `AnalysisError::Cancelled`
    /// propagates back to the caller, so the GUI can reset progress indicators
    /// synchronously rather than waiting for the `Err` to cross any async boundary.
    Cancelled,
}
```

---

## 4. Loads and Boundary Conditions

### 4.1 How Loads Enter the Problem

In the lower-bound formulation, loads appear on the right-hand side of the equilibrium equations:

```
B^T · m = λ · q_ext + q_dead
```

The load vectors `q_ext` and `q_dead` have entries at **free displacement DOFs** only. The value at each DOF is the work-conjugate transverse force.

#### 4.1.1 Distributed Area Load

A uniformly distributed load `p` [force/length²] is the standard case for RC slabs. It is converted to equivalent nodal forces by integrating against the quadratic displacement shape functions over each element area.

For a **uniform load p** over an element of area A, the **consistent nodal load vector** for the 6-node quadratic triangle is:

```
// Node ordering: [n1, n2, n3, m12, m23, m31]
//   n_i  = corner nodes
//   m_ij = midpoint node between corners i and j

q_element = p · A · [0, 0, 0, 1/3, 1/3, 1/3]^T
```

Corner nodes carry zero; midpoint nodes each carry `p·A/3`. Total = `p·A`. ✓

For a **spatially varying load** `p(x,y)`, integrate numerically using 3-point Gaussian quadrature (exact for quadratic variation):

```
q_i = ∫∫_element  p(x,y) · N_i(x,y) dA
```

**Variable load** (scaled by λ):
```rust
Load::AreaLoad {
    panel_id: "slab_panel_1".into(),
    intensity: LoadIntensity::Uniform(5.0),   // e.g. kN/m²
    load_case: LoadCase::Variable,
}
```

**Permanent load** (not scaled by λ):
```rust
Load::AreaLoad {
    panel_id: "slab_panel_1".into(),
    intensity: LoadIntensity::Uniform(1.5),   // e.g. self-weight kN/m²
    load_case: LoadCase::Permanent,
}
```

#### 4.1.2 Line Load on an Edge

A distributed line load `p` [force/length] along a boundary edge. For a quadratic edge of length L (corner A — midpoint M — corner B), the consistent nodal distribution is:

```
F_A = p · L / 6
F_B = p · L / 6
F_M = p · L · 2/3
```

```rust
Load::LineLoad {
    edge: EdgeRef::PolygonEdge { panel_id: "slab".into(), edge_index: 2 },
    intensity: 10.0,   // kN/m
    load_case: LoadCase::Variable,
}
```

#### 4.1.3 Dead Load vs. Variable Load

λ scales **only** `LoadCase::Variable` loads. Permanent loads appear as a constant offset `q_dead` in the equilibrium equations, independent of λ.

### 4.2 Boundary Conditions

BCs are applied by **removing constrained DOFs** from the system. The global equilibrium matrix is assembled for all DOFs, then rows and columns for constrained DOFs are deleted before passing to Clarabel. This is numerically clean and avoids ill-conditioning from penalty terms.

#### 4.2.1 Simply Supported Edge

Transverse displacement `w = 0` along the edge. Rotation is **free**.

- **Removes:** displacement DOFs (`w`) for all corner and midpoint nodes on the edge
- **Keeps:** rotation DOFs (`θ`) — the slab rotates freely at the support
- **Moment:** no constraint on edge moment (the yield criterion alone limits it)

```rust
Support::SimplySupportedEdge {
    edge: EdgeRef::PolygonEdge { panel_id: "slab".into(), edge_index: 0 },
}
```

#### 4.2.2 Clamped Edge

Both `w = 0` and rotation `θ = 0` along the edge.

- **Removes:** displacement DOFs (`w`) and rotation DOFs (`θ`) for all nodes on the edge
- The slab can develop full negative (hogging) moments at the clamped edge up to `m-_x`, `m-_y`

```rust
Support::ClampedEdge {
    edge: EdgeRef::PolygonEdge { panel_id: "slab".into(), edge_index: 1 },
}
```

#### 4.2.3 Free Edge

No constraints. **Default** for any edge not referenced by a support. Zero moment and shear at the free edge are natural boundary conditions, automatically satisfied by the lower-bound formulation.

#### 4.2.4 Symmetry Line

Enforces zero rotation normal to the symmetry axis and `mxy = 0` along the line.

- **Removes:** rotation DOFs normal to the symmetry axis
- **Adds:** constraint `mxy = 0` at nodes on the line (enforced via equality constraint in the equilibrium system)

```rust
Support::SymmetryLine {
    edge: EdgeRef::PolygonEdge { panel_id: "slab".into(), edge_index: 2 },
    axis: SymmetryAxis::X,   // or Y
}
```

### 4.3 Multi-Panel Slabs and Shared Edges

For a continuous slab made up of multiple panels, panels share edges. Along a shared interior edge:
- **Displacement continuity:** nodes on the shared edge have a single shared `w` DOF — no duplication
- **Moment continuity:** the `B^T_θ` equations from elements on both sides of the edge enforce continuous normal moments automatically
- **No support:** a shared interior edge is not a support; it is a pure continuity condition

In the mesher, all panels are triangulated in a single CDT operation. Shared edges are declared explicitly so the mesher places a single set of nodes on the interface, used by elements from both panels.

```rust
ModelBuilder::new("two_span_slab")
    .add_panel("span_1", Polygon2D::new(vec![
        [0.0, 0.0], [3.0, 0.0], [3.0, 4.0], [0.0, 4.0],
    ]))
    .add_panel("span_2", Polygon2D::new(vec![
        [3.0, 0.0], [6.0, 0.0], [6.0, 4.0], [3.0, 4.0],
    ]))
    // edge_index 1 of span_1 = edge from (3,0)→(3,4)
    // edge_index 3 of span_2 = edge from (3,4)→(3,0)
    .declare_shared_edge(SharedEdgeDeclaration::between("span_1", 1, "span_2", 3))
```

### 4.4 Load and Support Summary

| Type | Effect on primal problem | Effect on dual (collapse mode) |
|---|---|---|
| Area load, variable | Adds to `q_ext` at displacement DOFs | Virtual work in objective |
| Area load, permanent | Adds to `q_dead` (constant RHS) | Fixed work, not in objective |
| Line load | Adds to `q_ext` or `q_dead` at edge DOFs | — |
| Simply supported edge | Removes `w` DOFs on edge | Zero virtual displacement there |
| Clamped edge | Removes `w` and `θ` DOFs on edge | Zero virtual displacement and rotation |
| Free edge | No change (natural BC) | Virtual displacement/rotation free |
| Symmetry line | Removes normal rotation DOF; adds `mxy=0` | — |
| Shared panel edge | Merges DOFs from both panels | Continuous virtual displacement field |

---

## 5. API Design

### 5.1 Unit System

All quantities are stored internally in **base SI units** (N, m, Pa, N·m/m). The `UnitSystem` enum controls input and output units. Conversion is applied at the API boundary.

```rust
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub enum UnitSystem {
    /// N, m, Pa, N·m/m
    SI,
    /// kN, m, kPa, kN·m/m  (standard for RC design)
    KiloNewtonMetre,
    /// kN, mm, MPa, kN·mm/mm
    KiloNewtonMillimetre,
}
```

The `UnitSystem` is stored in `AnalysisModel`. All values in `SolveResult` are returned in the same unit system as the input model.

### 5.2 Edge Reference System

Edges are referenced by **panel ID + zero-based edge index**. For a polygon with N vertices defined as `[v0, v1, ..., v_{N-1}]`:

```
edge_index i  =  edge from vertex[i]  to  vertex[(i+1) % N]
```

This is unambiguous for any polygon and maps directly to GUI behaviour (click on an edge → resolve its index from the vertex list). The same indexing scheme applies to **hole polygons** within a panel — holes are also `Polygon2D` values with vertices listed in order, so their edges are addressed with the same integer index.

```rust
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub enum EdgeRef {
    /// A specific edge of the panel outline.
    /// edge_index i = edge from outline.vertices[i] to outline.vertices[(i+1) % n]
    PolygonEdge { panel_id: String, edge_index: usize },

    /// A specific edge of a hole within a panel.
    /// hole_index is the position of the hole in Panel::holes.
    /// edge_index follows the same vertex-order convention as PolygonEdge.
    HoleEdge { panel_id: String, hole_index: usize, edge_index: usize },

    /// All outline edges of a panel (convenience shorthand).
    AllEdges { panel_id: String },

    /// All edges of a specific hole (convenience shorthand).
    AllHoleEdges { panel_id: String, hole_index: usize },
}
```

**Example — panel with rectangular hole:**

The diagram below shows a 6×4 m panel with a 2×2 m opening. Vertices are listed counter-clockwise (CCW) for both the outline and the hole. Edge index `i` always runs from `vertex[i]` to `vertex[(i+1) % n]`, so the numbering follows the arrow direction.

```
    (0,4)──────────────────────(6,4)
      │   v3        edge 2        v2  │
      │                              │
  e   │      (2,3)────────(4,3)      │  e
  d   │        │ h3  ←  h2 │        │  d
  g   │        │            │        │  g
  e   │   h    │ hole 0     │    e   │  e
      │        │            │        │
  3   │        │ h0  →  h1 │        │  1
      │      (2,1)────────(4,1)      │
      │                              │
      │   v0        edge 0        v1  │
    (0,0)──────────────────────(6,0)
                    ↑ edge 0
```

```
// Panel outline — 4 vertices, 4 edges:
Polygon2D::new(vec![[0.0,0.0], [6.0,0.0], [6.0,4.0], [0.0,4.0]])
  PolygonEdge { edge_index: 0 }  →  (0,0) → (6,0)   south
  PolygonEdge { edge_index: 1 }  →  (6,0) → (6,4)   east
  PolygonEdge { edge_index: 2 }  →  (6,4) → (0,4)   north
  PolygonEdge { edge_index: 3 }  →  (0,4) → (0,0)   west

// Hole 0 — 4 vertices, 4 edges (CCW, same convention):
Polygon2D::new(vec![[2.0,1.0], [4.0,1.0], [4.0,3.0], [2.0,3.0]])
  HoleEdge { hole_index: 0, edge_index: 0 }  →  (2,1) → (4,1)   south face of hole
  HoleEdge { hole_index: 0, edge_index: 1 }  →  (4,1) → (4,3)   east  face of hole
  HoleEdge { hole_index: 0, edge_index: 2 }  →  (4,3) → (2,3)   north face of hole
  HoleEdge { hole_index: 0, edge_index: 3 }  →  (2,3) → (2,1)   west  face of hole
```

> **Note on hole orientation:** The hole polygon must be wound **counter-clockwise** when viewed in the standard x-right, y-up coordinate frame. This is the same convention as the panel outline. The mesher treats hole edges as internal constraints that bound a void — the triangulation fills the annular region between the outline and the hole, leaving the interior of the hole empty.

Hole edges are **free by default** (no constraint — the slab has an open void at that location). To support the hole perimeter, apply a BC using `HoleEdge` or `AllHoleEdges` exactly as you would for an outline edge.


### 5.3 ModelBuilder — Fluent Interface

`ModelBuilder` is the primary entry point. Builder methods take `Self` by value and return `Self`, enabling a clean method-chain without binding a `mut` variable. Terminal methods (`solve`, `build`, `mesh`) consume the builder.

```rust
pub struct ModelBuilder {
    name: String,
    units: UnitSystem,
    panels: Vec<Panel>,
    shared_edges: Vec<SharedEdgeDeclaration>,
    supports: Vec<Support>,
    loads: Vec<Load>,
    mesh_config: MeshConfig,
    solver_config: SolverConfig,
}

impl ModelBuilder {
    pub fn new(name: impl Into<String>) -> Self;

    // ── Configuration ──────────────────────────────────────────────────────
    pub fn with_units(self, units: UnitSystem) -> Self;
    pub fn set_mesh_config(self, config: MeshConfig) -> Self;
    pub fn set_solver_config(self, config: SolverConfig) -> Self;

    // Convenience shorthand for the most common mesh control:
    pub fn set_mesh_density(self, density: MeshDensity) -> Self;

    // ── Geometry ───────────────────────────────────────────────────────────
    /// Add a slab panel. Outline vertices must be in CCW order.
    pub fn add_panel(self, id: impl Into<String>, outline: Polygon2D) -> Self;

    /// Add a hole to a panel. Hole vertices must be in CCW order.
    /// Appended to Panel::holes in order; hole_index = number of previous add_hole calls for this panel.
    pub fn add_hole(self, panel_id: impl Into<String>, hole: Polygon2D) -> Self;

    /// Declare that two panel edges are geometrically coincident (shared interior edge).
    pub fn declare_shared_edge(self, decl: SharedEdgeDeclaration) -> Self;

    // ── Material ───────────────────────────────────────────────────────────
    /// Panics at validation time if panel_id is not found.
    pub fn set_material(self, panel_id: impl Into<String>, material: RCSlabMaterial) -> Self;

    // ── Boundary conditions ────────────────────────────────────────────────
    pub fn add_support(self, support: Support) -> Self;

    // ── Loads ──────────────────────────────────────────────────────────────
    pub fn add_load(self, load: Load) -> Self;

    // ── Terminal methods (consume the builder) ─────────────────────────────
    /// Validate, mesh, assemble, and solve with no progress reporting.
    /// Equivalent to `solve_with_progress` with a no-op callback.
    pub fn solve(self) -> Result<SolveResult, AnalysisError>;

    /// Validate, mesh, assemble, and solve with a progress callback.
    /// The callback receives `ProgressEvent`s and can return `ControlFlow::Cancel`
    /// to abort the analysis cleanly. Use this variant when a GUI stop button or
    /// progress indicator is needed.
    pub fn solve_with_progress(
        self,
        progress: &dyn ProgressCallback,
    ) -> Result<SolveResult, AnalysisError>;

    /// Validate and build the AnalysisModel without solving (for serialisation / GUI preview).
    pub fn build(self) -> Result<AnalysisModel, AnalysisError>;

    /// Validate and mesh without solving (for mesh preview in GUI).
    pub fn mesh(self) -> Result<MeshModel, AnalysisError>;
}
```

**Two-span slab example:**

```rust
use plato_api::prelude::*;

let result = ModelBuilder::new("two_span_rc_slab")
    .with_units(UnitSystem::KiloNewtonMetre)

    .add_panel("span_1", Polygon2D::new(vec![
        [0.0, 0.0], [3.0, 0.0], [3.0, 4.0], [0.0, 4.0],
    ]))
    .add_panel("span_2", Polygon2D::new(vec![
        [3.0, 0.0], [6.0, 0.0], [6.0, 4.0], [3.0, 4.0],
    ]))
    .declare_shared_edge(SharedEdgeDeclaration::between("span_1", 1, "span_2", 3))

    .set_material("span_1", RCSlabMaterial {
        m_pos_x: 10.0, m_pos_y: 10.0,
        m_neg_x: 12.0, m_neg_y: 12.0,
    })
    .set_material("span_2", RCSlabMaterial {
        m_pos_x: 10.0, m_pos_y: 10.0,
        m_neg_x: 12.0, m_neg_y: 12.0,
    })

    .set_mesh_density(MeshDensity::MaxElementArea(0.05))

    .add_support(Support::SimplySupportedEdge {
        edge: EdgeRef::PolygonEdge { panel_id: "span_1".into(), edge_index: 0 },
    })
    .add_support(Support::SimplySupportedEdge {
        edge: EdgeRef::PolygonEdge { panel_id: "span_1".into(), edge_index: 2 },
    })
    .add_support(Support::SimplySupportedEdge {
        edge: EdgeRef::PolygonEdge { panel_id: "span_1".into(), edge_index: 3 },
    })
    .add_support(Support::SimplySupportedEdge {
        edge: EdgeRef::PolygonEdge { panel_id: "span_2".into(), edge_index: 0 },
    })
    .add_support(Support::SimplySupportedEdge {
        edge: EdgeRef::PolygonEdge { panel_id: "span_2".into(), edge_index: 1 },
    })
    .add_support(Support::SimplySupportedEdge {
        edge: EdgeRef::PolygonEdge { panel_id: "span_2".into(), edge_index: 2 },
    })

    .add_load(Load::AreaLoad {
        panel_id: "span_1".into(),
        intensity: LoadIntensity::Uniform(1.0),
        load_case: LoadCase::Variable,
    })
    .add_load(Load::AreaLoad {
        panel_id: "span_2".into(),
        intensity: LoadIntensity::Uniform(1.0),
        load_case: LoadCase::Variable,
    })

    .solve()?;

println!("Collapse load: {:.3} kN/m²", result.load_factor);
```

**Slab with hole example:**

```rust
// 6×4 m simply-supported slab with a 2×2 m rectangular opening.
// add_hole() appends to panel.holes — hole_index is 0 (first hole added).
let result = ModelBuilder::new("slab_with_opening")
    .with_units(UnitSystem::KiloNewtonMetre)

    .add_panel("slab",
        Polygon2D::new(vec![[0.0,0.0],[6.0,0.0],[6.0,4.0],[0.0,4.0]]),
    )
    .add_hole("slab",
        Polygon2D::new(vec![[2.0,1.0],[4.0,1.0],[4.0,3.0],[2.0,3.0]]),
    )
    .set_material("slab", RCSlabMaterial {
        m_pos_x: 10.0, m_pos_y: 10.0,
        m_neg_x: 10.0, m_neg_y: 10.0,
    })
    .set_mesh_density(MeshDensity::MaxElementArea(0.05))
    .add_support(Support::SimplySupportedEdge {
        edge: EdgeRef::AllEdges { panel_id: "slab".into() },
    })
    // Hole edges are free by default — no BC entry needed.
    // To support the hole perimeter:
    //   .add_support(Support::SimplySupportedEdge {
    //       edge: EdgeRef::AllHoleEdges { panel_id: "slab".into(), hole_index: 0 },
    //   })
    .add_load(Load::AreaLoad {
        panel_id: "slab".into(),
        intensity: LoadIntensity::Uniform(1.0),
        load_case: LoadCase::Variable,
    })
    .solve()?;
```

### 5.4 Core Types

```rust
// ── Top-level model ────────────────────────────────────────────────────────

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnalysisModel {
    pub name: String,
    pub units: UnitSystem,
    pub panels: Vec<Panel>,
    pub shared_edges: Vec<SharedEdgeDeclaration>,
    pub supports: Vec<Support>,
    pub loads: Vec<Load>,
    pub mesh_config: MeshConfig,
    pub solver_config: SolverConfig,
}

// ── Geometry ───────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Panel {
    pub id: String,
    /// Outer boundary, vertices in CCW order.
    pub outline: Polygon2D,
    /// Interior voids. Each hole is a CCW polygon. hole_index = position in this Vec.
    pub holes: Vec<Polygon2D>,
    pub material: RCSlabMaterial,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Polygon2D {
    /// Vertices in counter-clockwise order.
    /// Edge i runs from vertices[i] to vertices[(i+1) % len].
    pub vertices: Vec<[f64; 2]>,
}

impl Polygon2D {
    pub fn new(vertices: Vec<[f64; 2]>) -> Self;
    /// Convenience constructor for axis-aligned rectangles.
    pub fn rectangle(x_min: f64, y_min: f64, x_max: f64, y_max: f64) -> Self;
    pub fn n_edges(&self) -> usize;  // == vertices.len()
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SharedEdgeDeclaration {
    pub panel_a: String,
    pub edge_a: usize,   // edge_index in panel_a's outline
    pub panel_b: String,
    pub edge_b: usize,   // edge_index in panel_b's outline (opposite orientation)
}

impl SharedEdgeDeclaration {
    pub fn between(
        panel_a: impl Into<String>, edge_a: usize,
        panel_b: impl Into<String>, edge_b: usize,
    ) -> Self;
}

// ── Material ───────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RCSlabMaterial {
    /// Positive (sagging/bottom reinforcement) moment capacities [moment/length].
    pub m_pos_x: f64,
    pub m_pos_y: f64,
    /// Negative (hogging/top reinforcement) moment capacities [moment/length].
    /// Stored as positive values.
    pub m_neg_x: f64,
    pub m_neg_y: f64,
}

impl RCSlabMaterial {
    /// Isotropic convenience constructor (equal capacity in all directions).
    pub fn isotropic(m_pos: f64, m_neg: f64) -> Self;
}

// ── Supports and loads ─────────────────────────────────────────────────────

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Support {
    SimplySupportedEdge { edge: EdgeRef },
    ClampedEdge         { edge: EdgeRef },
    FreeEdge            { edge: EdgeRef },   // explicit form of the default
    SymmetryLine        { edge: EdgeRef, axis: SymmetryAxis },
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub enum SymmetryAxis { X, Y }

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Load {
    AreaLoad {
        panel_id: String,
        intensity: LoadIntensity,
        load_case: LoadCase,
    },
    LineLoad {
        edge: EdgeRef,
        intensity: f64,   // [force/length]
        load_case: LoadCase,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum LoadIntensity {
    Uniform(f64),   // [force/length²]
    // Future: Varying(Box<dyn Fn(f64, f64) -> f64>)
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub enum LoadCase {
    Variable,    // multiplied by λ
    Permanent,   // constant; not multiplied by λ
}

// ── Mesh configuration ─────────────────────────────────────────────────────
// Full definitions of MeshConfig and RefinementZone are in section 6.2 (plato-mesh).
// MeshDensity is a convenience enum used by ModelBuilder::set_mesh_density().

/// Convenience enum for the most common mesh control parameter.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum MeshDensity {
    /// Set maximum element area directly [length²].
    MaxElementArea(f64),
    /// Target approximately this many elements across the shortest panel dimension.
    ElementsAcross(usize),
}

// ── Solver configuration ───────────────────────────────────────────────────

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SolverConfig {
    /// Absolute duality gap tolerance. Default: 1e-6.
    pub tol_gap_abs: f64,
    /// Relative duality gap tolerance. Default: 1e-6.
    pub tol_gap_rel: f64,
    /// Maximum solver iterations. Default: 200.
    pub max_iter: usize,
    /// Enable diagonal scaling of constraint matrix. Default: true.
    pub scaling: bool,
}

impl Default for SolverConfig {
    fn default() -> Self {
        Self { tol_gap_abs: 1e-6, tol_gap_rel: 1e-6, max_iter: 200, scaling: true }
    }
}

// ── Solve result ───────────────────────────────────────────────────────────

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum SolveStatus {
    Optimal,
    /// Solver reached max iterations without converging. Result may be approximate.
    MaxIterationsReached,
    /// Problem is infeasible (no admissible moment field exists — check BCs and loads).
    Infeasible,
    /// Numerical failure inside the solver.
    NumericalError(String),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SolveResult {
    pub status: SolveStatus,
    pub load_factor: f64,
    pub units: UnitSystem,
    pub solve_time_ms: u64,

    // ── Embedded mesh geometry (self-contained for GUI rendering) ──────────
    /// Node coordinates [x, y], indexed by global node id.
    pub nodes: Vec<[f64; 2]>,
    /// Element connectivity. Indices reference `nodes`.
    pub elements: Vec<SolveElement>,

    // ── Primal solution ────────────────────────────────────────────────────
    pub element_moments: Vec<ElementMoments>,

    // ── Dual solution (collapse mode) ──────────────────────────────────────
    /// Transverse displacement w at each node (indexed by global node id).
    /// Magnitude is arbitrary; use for collapse shape visualisation only.
    pub nodal_displacements: Vec<f64>,

    // ── Diagnostics ────────────────────────────────────────────────────────
    pub n_elements: usize,
    pub n_nodes: usize,
    pub n_free_dofs: usize,
    /// Total Clarabel decision variables: 1 + 21·N_e
    pub n_variables: usize,
    pub solver_iterations: usize,
    pub duality_gap: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SolveElement {
    pub id: usize,
    pub panel_id: String,
    /// Indices into SolveResult::nodes: [corner_0, corner_1, corner_2, mid_01, mid_12, mid_20]
    pub nodes: [usize; 6],
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ElementMoments {
    pub element_id: usize,
    /// Moments at 3 corner nodes: [[mx, my, mxy], [mx, my, mxy], [mx, my, mxy]]
    pub corner_moments: [[f64; 3]; 3],
    /// Yield utilisation at each corner node: 0.0 = elastic, 1.0 = fully plastic.
    pub yield_utilisation: [f64; 3],
}

// ── Error types ────────────────────────────────────────────────────────────

#[derive(Debug, thiserror::Error)]
pub enum AnalysisError {
    #[error("Mesh error: {0}")]
    Mesh(#[from] MeshError),
    #[error("Assembly error: {0}")]
    Assembly(#[from] AssemblyError),
    #[error("Model validation failed: {0}")]
    Validation(String),
    #[error("Solver error: {0}")]
    Solver(String),
    /// Returned when the progress callback signals `ControlFlow::Cancel`.
    /// Distinct from `Solver` so the GUI can show "Cancelled" rather than "Error".
    #[error("Analysis cancelled by user")]
    Cancelled,
}

#[derive(Debug, thiserror::Error)]
pub enum MeshError {
    #[error("Panel '{0}' not found")]
    UnknownPanel(String),
    #[error("Panel polygon is not valid (e.g. self-intersecting or too few vertices): {0}")]
    InvalidPolygon(String),
    #[error("Hole polygon is outside or overlaps panel outline")]
    InvalidHole,
    #[error("Shared edge geometry mismatch between '{panel_a}' edge {edge_a} and '{panel_b}' edge {edge_b}")]
    SharedEdgeMismatch { panel_a: String, edge_a: usize, panel_b: String, edge_b: usize },
    #[error("Mesher failed to triangulate: {0}")]
    TriangulationFailed(String),
}

#[derive(Debug, thiserror::Error)]
pub enum AssemblyError {
    #[error("No free DOFs after applying boundary conditions — model is fully constrained")]
    FullyConstrained,
    #[error("Element {0} is degenerate (zero or near-zero area)")]
    DegenerateElement(usize),
}

// ── Progress reporting and cancellation ────────────────────────────────────

/// Returned by `ProgressCallback::on_event`. The solver checks this at every
/// iteration boundary and aborts with `AnalysisError::Cancelled` if `Cancel`
/// is returned.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ControlFlow {
    /// Continue the analysis normally.
    Continue,
    /// Abort the analysis cleanly at the next iteration boundary.
    Cancel,
}

/// Inject a progress callback into `run_analysis` or `ModelBuilder::solve_with_progress`
/// to receive structured events and optionally cancel the analysis.
///
/// `on_event` is called at every pipeline stage transition and at every solver
/// iteration. Returning `ControlFlow::Cancel` from **any** call causes the solver
/// to abort at the earliest safe opportunity and return `AnalysisError::Cancelled`.
///
/// "Earliest safe opportunity" means:
/// - During meshing or assembly: before the next pipeline stage begins.
/// - During solving: at the next iteration boundary.
///
/// In all cases the engine exits without corrupting internal state.
pub trait ProgressCallback: Send + Sync {
    fn on_event(&self, event: ProgressEvent) -> ControlFlow;
}

/// Blanket implementation so a plain closure can be used as a callback.
impl<F: Fn(ProgressEvent) -> ControlFlow + Send + Sync> ProgressCallback for F {
    fn on_event(&self, event: ProgressEvent) -> ControlFlow { self(event) }
}
```

---

## 6. Mesh Generation

### 6.1 Spade Mesher

Spade provides a Constrained Delaunay Triangulation (CDT) with Ruppert's refinement algorithm. The mesher operates in two passes:

1. **CDT pass:** Spade generates a 3-node triangulation of all panels together, respecting all boundary edges (outer polygon edges, hole edges, shared panel edges as constrained segments)
2. **Midpoint insertion:** `plato-mesh` inserts midpoint nodes on each triangle side to produce 6-node elements

For multi-panel models, the single CDT ensures shared edges produce one set of coincident nodes — no duplication or gaps at interfaces.

```rust
pub trait Mesher: Send + Sync {
    fn triangulate(
        &self,
        panels: &[Panel],
        shared_edges: &[SharedEdgeDeclaration],
        config: &MeshConfig,
    ) -> Result<MeshModel, MeshError>;
}

pub struct SpadeMesher;
// Future: pub struct TriangleMesher;  (C bindings, higher quality)
```

### 6.2 Mesh Configuration

```rust
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MeshConfig {
    /// Maximum triangle area [length²]. Controls global mesh density.
    pub max_element_area: f64,

    /// Minimum interior angle [degrees] for Ruppert refinement.
    /// Default: 20.0. Values above 33.0 risk non-termination.
    pub min_angle_deg: f64,

    /// Local refinement zones (e.g. around re-entrant corners).
    pub refinement_zones: Vec<RefinementZone>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RefinementZone {
    pub centre: [f64; 2],
    pub radius: f64,
    pub max_area: f64,   // must be smaller than global max_element_area
}
```

### 6.3 MeshModel

```rust
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MeshModel {
    /// Node coordinates, indexed by node_id
    pub nodes: Vec<[f64; 2]>,

    /// All 6-node triangular elements
    pub elements: Vec<TriElement6>,

    /// Ordered node ids along each named edge (for BC/load application).
    /// Keys can be PolygonEdge, HoleEdge, or AllEdges/AllHoleEdges — resolved at mesh time.
    pub edge_nodes: HashMap<EdgeRef, Vec<usize>>,

    pub quality: MeshQuality,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TriElement6 {
    pub id: usize,
    /// Corner nodes [n0, n1, n2] in counter-clockwise order.
    pub corners: [usize; 3],
    /// Midpoint nodes: midpoints[i] is the midpoint of the side between corners i
    /// and corner (i+1) % 3. i.e.:
    ///   midpoints[0] = midpoint of side between corners 0 and 1  (mid_01)
    ///   midpoints[1] = midpoint of side between corners 1 and 2  (mid_12)
    ///   midpoints[2] = midpoint of side between corners 2 and 0  (mid_20)
    /// This matches the SolveElement node ordering [c0, c1, c2, mid_01, mid_12, mid_20].
    pub midpoints: [usize; 3],
    pub panel_id: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MeshQuality {
    pub n_elements: usize,
    pub n_nodes: usize,
    pub min_angle_deg: f64,
    pub max_aspect_ratio: f64,
}
```

---

## 7. Problem Assembly and Optimisation

### 7.1 Degrees of Freedom

| DOF type | Count | Notes |
|---|---|---|
| Transverse displacement `w` | 1 per node | Constrained DOFs removed for BCs |
| Moment variables `m` | 9 per element | 3 corners × 3 components |
| Auxiliary variables `α` | 12 per element | 3 corners × 4 per corner (Johansen) |
| Load factor `λ` | 1 global | The objective variable |

For `N_n` total nodes (corner nodes + midpoint nodes; **not** counting the 6-node rotation nodes which have `θ` DOFs), `N_e` elements, `N_c` constrained displacement DOFs:

```
Displacement nodes:      N_n  (corner + midpoint nodes only; w DOFs)
Rotation nodes:          6 · N_e  (θ DOFs, internal to each element — shared across elements on shared edges)
Free displacement DOFs:  N_dof = N_n - N_c
Free rotation DOFs:      counted separately; also included in equilibrium rows via B^T_θ
Moment variables:        9 · N_e
Auxiliary variables:     12 · N_e
Total Clarabel variables: 1 + 21 · N_e   (λ + moments + auxiliaries only; θ DOFs are equilibrium rows, not decision variables)
```

> Note: rotation DOFs appear as **rows** in the global equilibrium matrix (via B^T_θ) but are **not** decision variables in the Clarabel problem. The moment continuity they enforce is expressed as equality constraints. This is why `N_dof` refers only to displacement (w) DOFs.

### 7.2 DOF Map and BC Application

The `DofMap` is constructed after meshing, before assembly:

1. Assign a provisional DOF index to every node's `w` DOF and every rotation node's `θ` DOF
2. For each `Support` in the model, mark the relevant DOFs as `Constrained`
3. Assign contiguous free DOF indices (0, 1, 2, ...) to all `Free` DOFs
4. Store: `(node_id, DofType) → FreeIndex | Constrained`

The equilibrium matrix is then assembled using free DOF indices only.

### 7.3 Global Equilibrium Matrix Assembly

The global `B^T` has shape `N_dof × 9·N_e`. Assembly is embarrassingly parallel:

```
for each element e (parallelised with rayon):
    1. Fetch corner coordinates from MeshModel
    2. Compute ElementGeometry (A, l_i, n_i, n̂_i, P_i, a_ij, b_i, c_ij)
    3. Compute B^T_r, B^T_q, B^T_t  →  B^T_disp  (6 rows × 9 cols)
    4. Compute B^T_θ                              (6 rows × 9 cols)
    5. Map local node ids → free global DOF indices via DofMap
       (skip rows where the local node is constrained)
    6. Emit (global_row, moment_col, value) triplets

Finalise:
    Collect all triplets → sort by column → build CSC matrix
```

### 7.4 Clarabel Problem Formulation

**Decision variable vector:**

```
x = [λ,  m_1, ..., m_{N_e},  α_1, ..., α_{N_e}]^T

m_e = [mx_1, my_1, mxy_1,  mx_2, my_2, mxy_2,  mx_3, my_3, mxy_3]^T   (9 per element)

α_e = [α1_1, α2_1, α3_1, α4_1,    ← corner 1
       α1_2, α2_2, α3_2, α4_2,    ← corner 2
       α1_3, α2_3, α3_3, α4_3]^T  ← corner 3  (12 per element)
```

**Objective:**

```
P = 0    (linear objective — no quadratic term)
q = [−1,  0, 0, ..., 0]^T    (minimise −λ  ≡  maximise λ)
```

**Equality constraints (equilibrium + auxiliary variable definitions):**

```
Equilibrium:
    [−q_ext | B^T_global | 0_aux] · x = q_dead
    shape: N_dof × (1 + 9·N_e + 12·N_e)
    (column blocks match x = [λ, m_1..m_{N_e}, α_1..α_{N_e}]^T)

Auxiliary definitions (per element e, per corner k):
    α1_k + (1/√2)·mx_k          =  (1/√2)·m+_x
    α2_k + (1/√2)·my_k          =  (1/√2)·m+_y
    α3_k − (1/√2)·mx_k          =  (1/√2)·m−_x
    α4_k − (1/√2)·my_k          =  (1/√2)·m−_y
    (4 equalities per corner × 3 corners × N_e elements)
```

**Cone constraints (per element e, per corner k ∈ {1,2,3}):**

```
// Non-negativity (NonnegativeCone, 4 constraints per corner):
α1_k ≥ 0,  α2_k ≥ 0,  α3_k ≥ 0,  α4_k ≥ 0

// Positive yield surface — SecondOrderCone(3), expressed as affine SOC:
// Clarabel form: ||z||₂ ≤ t, i.e. x[0] ≥ ||(x[1], x[2])||₂
//
//   t  =  (α1_k + α2_k) / √2        (affine combination of α variables)
//   s  =  (α1_k − α2_k) / √2        (affine combination of α variables)
//   w  =  mxy_k                      (moment variable, directly)
//
// SOC block rows (A_soc · x entries):
//   row 0:  (1/√2)·α1_k  +  (1/√2)·α2_k      →  t
//   row 1:  (1/√2)·α1_k  −  (1/√2)·α2_k      →  s
//   row 2:  mxy_k                               →  w
//
// i.e.  [t; s; w] ∈ SOC(3)    ↔    ||(s, mxy_k)||₂ ≤ t

// Negative yield surface — SecondOrderCone(3):
//   t  =  (α3_k + α4_k) / √2
//   s  =  (α3_k − α4_k) / √2
//   w  =  mxy_k
//
//   row 0:  (1/√2)·α3_k  +  (1/√2)·α4_k      →  t
//   row 1:  (1/√2)·α3_k  −  (1/√2)·α4_k      →  s
//   row 2:  mxy_k                               →  w
//
// i.e.  [t; s; w] ∈ SOC(3)

// Note: no extra t/s variables are added to x — the SOC rows are purely affine
//       expressions of the existing α and moment variables via A_cone · x.
// Note: mxy_k appears in both SOC blocks — this is correct and expected.
```

Per element totals: 12 nonneg + 6 SOC(3) cone constraints.

**Total problem size:**

```
Equality rows:     N_dof  +  12·N_e   (equilibrium + aux definitions)
Cone rows:         12·N_e  +  18·N_e  =  30·N_e   (nonneg + SOC, 4 nonneg + 6 SOC rows per corner × 3)
Variables:         1  +  21·N_e       (λ + moments + α auxiliaries; t/s are not explicit variables)
```

### 7.5 Numerical Considerations

- **Scaling:** Apply diagonal row/column scaling before passing to Clarabel. Target all non-zero entries in the constraint matrix to be O(1). Store scaling factors in `SolveResult` and invert to recover physical moments.
- **Sparsity:** Each element's 9 moment unknowns couple to at most 12 free displacement DOFs. Expected matrix density: < 1% for typical meshes.
- **Warm-starting:** Expose Clarabel's warm-start interface in `SolverConfig` for parametric studies.
- **Convergence warnings:** If `duality_gap > 1e-4` at termination, emit a warning — the result may be unreliable.

---

## 8. Implementation Plan

### Phase 0 — Workspace and Infrastructure
- Cargo workspace; crate stubs; CI (`cargo test`, `cargo clippy`, `cargo fmt --check`)
- Benchmark harness (criterion.rs)
- SVG mesh export for visual validation during development

### Phase 1 — Geometry and Mesh (`plato-mesh`)
- `Polygon2D`: vertex list, edge count, edge-index validation
- `SpadeMesher`: single-panel CDT → 3-node triangles → 6-node elements
- `MeshModel` serialisation roundtrip test
- Multi-panel CDT: shared edge detection, unified triangulation, node merging
- Edge → ordered node list resolution

### Phase 2 — Element Assembly (`plato-core`)
- `ElementGeometry`: A, `l_i`, `n_i`, `n̂_i`, `P_i`, `a_ij`, `b_i`, `c_ij`
- `PlateElement::local_bt()`: `B^T_r + B^T_q + B^T_t + B^T_θ`
- Unit test: hand-verified `B^T` for a unit equilateral triangle
- `DofMap`: provisional → constrained/free assignment; contiguous free index assignment
- Global CSC assembly (triplet collection → sort → CSC)

### Phase 3 — Johansen Criterion and Clarabel Integration (`plato-core`)
- `JohansenCriterion::cone_blocks()`: nonneg + affine SOC blocks per corner, per element
- `ClarabelProblem` builder: assemble `P, q, A_eq, b_eq, cones`
- First full solve: 2×2 m isotropic simply-supported square slab, uniform unit load
- Verify: `λ = 24·m_p / L² = 6.0 kN/m²` for `m_p = 1 kNm/m`, `L = 2 m`

### Phase 4 — Loads and Boundary Conditions (`plato-api`)
- `AreaLoad` consistent nodal distribution (uniform)
- `LineLoad` consistent nodal distribution along edges
- `SimplySupportedEdge`, `ClampedEdge`, `FreeEdge`: DOF removal
- `SymmetryLine`: rotation DOF removal + `mxy = 0` enforcement
- `UnitSystem` conversion layer
- Multi-panel shared edge DOF merging

### Phase 5 — Collapse Mode and Output
- Dual solution extraction from Clarabel output
- `nodal_displacements` assembly
- `yield_utilisation` per corner node per element
- `SolveResult` serialisation roundtrip test

### Phase 6 — Validation

| Test case | Target λ | Reference |
|---|---|---|
| Simply supported square slab (L×L), isotropic `m_p`, uniform | `24 m_p / L²` | Johansen yield-line (exact) |
| Clamped square slab (L×L), isotropic `m_p`, uniform | `≈ 42.4 m_p / L²` | Nielsen & Hoang (2010) |
| Two-span continuous slab, simply supported outer edges | Yield-line calculation | Hand calculation |
| Slab with rectangular hole, simply supported | Converges from below | Mesh convergence check |
| Orthotropic simply-supported square | `24 √(m+_x · m+_y) / L²` | Johansen equivalent |

### Phase 7 — Robustness and Performance
- Rayon parallel element assembly
- Sparse matrix benchmarks (assembly time vs. element count)
- Edge cases: near-degenerate triangles, very small holes, high-aspect panels
- Warning system: aspect ratio, near-singular geometry, solver non-convergence

---

## 9. Dependencies

| Crate | Role | Notes |
|---|---|---|
| `clarabel` | SOCP solver | Native Rust, no FFI. Standard SOC supported. CSC input. |
| `spade` ≥ 2.0 | CDT meshing | Pure Rust. Constrained edges for holes and shared panel boundaries. Ruppert refinement. |
| `nalgebra` | Dense linear algebra | Local element matrix ops (2×2 rotations, small dense B^T blocks). |
| `faer` | Sparse matrix construction | Fast CSC assembly. Fallback: `sprs`. |
| `serde` + `serde_json` + `toml` | Serialisation | Full roundtrip for all public types; native save format is `.plato` (JSON under the hood). |
| `rayon` | Parallelism | Data-parallel element assembly. |
| `thiserror` | Error handling | Typed error enums with context chains. |
| `tracing` | Internal logging and diagnostics | Structured spans per pipeline stage. Not the GUI progress mechanism — that is `ProgressCallback`. |
| `criterion` (dev) | Benchmarking | Solve-time tracking against analytical targets. |
| `proptest` (dev) | Property tests | Equilibrium identity for random admissible moment fields. |
| `approx` (dev) | Float assertions | Tolerance-based validation comparisons. |

---

## 10. Open Questions

### Resolved

| Question | Decision |
|---|---|
| Reinforcement variation | Constant per panel |
| Reinforcement orientation | x/y-aligned only (no skewed reinforcement in v1) |
| Number of load cases | One variable load pattern per analysis |
| Multi-panel slabs | Supported via `declare_shared_edge` |
| Unit system | `UnitSystem` enum; internal SI; output in user units |
| BC implementation | DOF removal |
| Yield criterion enforcement | All 3 corner nodes per element |
| Point supports | Not in v1 |
| Elastic line supports | Not in scope |
| In-plane action | Not in scope (pure bending) |
| `SolveResult` geometry embedding | Embed `nodes` + `elements` directly in `SolveResult` for self-contained GUI rendering |
| Hole boundary conditions | Supported — holes defined on `Panel`, addressed via `EdgeRef::HoleEdge { panel_id, hole_index, edge_index }` using the same vertex-order convention as outline edges |

### Still Open

None — all design questions are resolved for v1. The document is ready for implementation.

---

## Appendix A — Mathematical Reference

### A.1 Element Geometric Quantities

For a triangle with corners 1, 2, 3 (counter-clockwise) at coordinates `(x_i, y_i)`:

```
Area:
    A = (1/2) |(x2−x1)(y3−y1) − (x3−x1)(y2−y1)|

Side i is the side OPPOSITE node i:
    Side 1: between nodes 2 and 3
    Side 2: between nodes 3 and 1
    Side 3: between nodes 1 and 2

Side lengths:
    l_1 = sqrt((x3−x2)² + (y3−y2)²)
    l_2 = sqrt((x1−x3)² + (y1−y3)²)
    l_3 = sqrt((x2−x1)² + (y2−y1)²)

Outward unit normals (pointing away from opposite node):
    n_1 = (1/l_1) · [ (y3−y2), −(x3−x2) ]
    n_2 = (1/l_2) · [ (y1−y3), −(x1−x3) ]
    n_3 = (1/l_3) · [ (y2−y1), −(x2−x1) ]

    Verification: n_i · (x_i − x_j) < 0  for any j ≠ i on the opposite side.

Tangent along side i (90° CCW from outward normal):
    n̂_i = [−n_iy, n_ix]

Moment traction projection matrix (maps m=[mx,my,mxy] to traction on side i):
    P_i = [ [n_ix,  0,    n_iy],
            [0,     n_iy, n_ix] ]     shape: 2×3

Auxiliary vectors:
    a_ij = P_i · n_j               shape: (2,)
    b_i  = P_i · n̂_i              shape: (2,)
    c_ij = l_i · l_j · a_ij       shape: (2,)
```

### A.2 Local B^T Matrices (Explicit)

**Row ordering:** `[n1, n2, n3, m12, m23, m31]` for displacement rows; `[θ1a, θ1b, θ2a, θ2b, θ3a, θ3b]` for rotation rows.

**Column ordering:** Blocks for `[m1, m2, m3]`, each block a 3-vector `[mx, my, mxy]`.

All vectors `a_ij^T, b_i^T, c_ij^T` are row vectors of length 3.

```
B^T_r   (6 displacement rows):
    row n1:  [ (b2−b3)^T  |     0      |      0     ]
    row n2:  [     0      | (b3−b1)^T  |      0     ]
    row n3:  [     0      |     0      | (b1−b2)^T  ]
    row m12: [     0      |     0      |      0     ]
    row m23: [     0      |     0      |      0     ]
    row m31: [     0      |     0      |      0     ]

B^T_q   (6 displacement rows, prefactor 1/(12A)):
    row n1:  [ −c11^T | −c21^T | −c31^T ]
    row n2:  [ −c12^T | −c22^T | −c32^T ]
    row n3:  [ −c13^T | −c23^T | −c33^T ]
    row m12: [  4c11^T|  4c21^T|  4c31^T]
    row m23: [  4c12^T|  4c22^T|  4c32^T]
    row m31: [  4c13^T|  4c23^T|  4c33^T]

B^T_t   (6 displacement rows, prefactor 1/6):
    row n1:  [ (b3−b2)^T  | −b3^T      |  b2^T      ]
    row n2:  [  b3^T      | (b1−b3)^T  | −b1^T      ]
    row n3:  [ −b2^T      |  b1^T      | (b2−b1)^T  ]
    row m12: [    0       |  4b1^T     | −4b1^T     ]
    row m23: [ −4b2^T     |    0       |  4b2^T     ]
    row m31: [  4b3^T     | −4b3^T     |    0       ]

B^T_θ   (6 rotation rows):
    row θ1a: [    0    |  a11^T  |    0    ]
    row θ1b: [    0    |    0    | −a11^T  ]
    row θ2a: [    0    |    0    |  a22^T  ]
    row θ2b: [ −a22^T  |    0    |    0    ]
    row θ3a: [  a33^T  |    0    |    0    ]
    row θ3b: [    0    | −a33^T  |    0    ]
```

> **Sign verification:** Apply a uniform moment state `m = [m0, m0, 0]^T` at all 3 corners (constant field). The assembled `B^T · m_uniform` should equal the equivalent nodal forces for a uniform distributed moment, and moment continuity rows should give zero residual. Use this as the patch test.

### A.3 Johansen SOC Blocks for Clarabel

Full cone constraint specification for element `e`, corner `k`, material parameters `m+_x, m+_y, m-_x, m-_y`.

**Background: RSOC → SOC conversion**

The Johansen yield surfaces have the form `2·u·v ≥ w²` (a Rotated SOC). Since Clarabel only supports standard SOC, convert via:

```
2·u·v ≥ w²,  u ≥ 0, v ≥ 0
  ↔   t² ≥ s² + w²          where t = (u+v)/√2,  s = (u−v)/√2
  ↔   ||(s, w)||₂ ≤ t
  ↔   [t; s; w] ∈ SOC(3)
```

**Index layout in x:**

```
mx_k  at index: moment_offset_e + 3*(k-1) + 0
my_k  at index: moment_offset_e + 3*(k-1) + 1
mxy_k at index: moment_offset_e + 3*(k-1) + 2
α1_k  at index: aux_offset_e    + 4*(k-1) + 0
α2_k  at index: aux_offset_e    + 4*(k-1) + 1
α3_k  at index: aux_offset_e    + 4*(k-1) + 2
α4_k  at index: aux_offset_e    + 4*(k-1) + 3
```

**Auxiliary definition equalities (rows in A_eq, b_eq):**

```
α1_k + (1/√2)·mx_k  =  (1/√2)·m+_x
α2_k + (1/√2)·my_k  =  (1/√2)·m+_y
α3_k − (1/√2)·mx_k  =  (1/√2)·m-_x
α4_k − (1/√2)·my_k  =  (1/√2)·m-_y
```

**NonnegativeCone (4 per corner):**

```
[α1_k, α2_k, α3_k, α4_k] ∈ NonnegativeCone(4)
```

**Positive yield surface — SecondOrderCone(3):**

```
Affine rows selecting from x (these form the 3-vector fed into SOC(3)):
    row 0 [t]:  (1/√2)·α1_k + (1/√2)·α2_k
    row 1 [s]:  (1/√2)·α1_k − (1/√2)·α2_k
    row 2 [w]:  mxy_k

Constraint: [t; s; w] ∈ SOC(3)   ↔   ||(s, mxy_k)||₂ ≤ t

Interpretation: t = (α1_k+α2_k)/√2,  s = (α1_k−α2_k)/√2
                t² − s² = 2·α1_k·α2_k ≥ mxy_k²   ✓
```

**Negative yield surface — SecondOrderCone(3):**

```
    row 0 [t]:  (1/√2)·α3_k + (1/√2)·α4_k
    row 1 [s]:  (1/√2)·α3_k − (1/√2)·α4_k
    row 2 [w]:  mxy_k

Constraint: [t; s; w] ∈ SOC(3)   ↔   ||(s, mxy_k)||₂ ≤ t
```

**Key points:**
- `t` and `s` are **not** added to the decision vector `x` — they are purely affine expressions computed via the cone constraint matrix rows
- `mxy_k` appears in both SOC blocks as the same variable — this is correct; both yield surfaces are bounded by the same twisting moment
- The `NonnegativeCone` on `α1..α4` is still required; the SOC alone does not enforce `α ≥ 0`

### A.4 Analytical Validation Solutions

**Simply supported square slab, side L, isotropic capacity `m_p`, uniform load `p`:**

```
p_collapse = 24 · m_p / L²
```

The lower-bound FE solution converges from below: `p_FE(h) ≤ p_collapse → p_collapse` as `h → 0`.

**Clamped square slab, isotropic `m_p` (equal positive and negative), uniform load:**

```
p_collapse ≈ 42.4 · m_p / L²
```

**Orthotropic simply-supported square slab, uniform load:**

```
p_collapse = 24 · sqrt(m+_x · m+_y) / L²
```

---

*This document is the primary reference for Claude Code implementation.*
*Mathematical notation follows Krabbenhøft (2016) and Larsen (2017), DTU Civil Engineering.*
*Sign conventions in B^T_θ must be verified against the patch test before proceeding to integration tests.*
