use std::collections::{HashMap, HashSet};
use std::time::Instant;

use plato_core::assembly::{DofMap, assemble_bt, collect_all_displacement_nodes};
use plato_core::criteria::JohansenCriterion;
use plato_core::problem::ClarabelProblem;
use plato_core::result::SolveStatus as CoreStatus;
use plato_mesh::geometry::EdgeRef;
use plato_mesh::mesher::{Mesher, SpadeMesher};
use plato_mesh::model::MeshModel;

use crate::loads::{Load, LoadCase, LoadIntensity};
use crate::model::{
    AnalysisError, AnalysisModel, AssemblyError, SolveElement, SolveResult, SolveStatus,
};
use crate::supports::Support;

// ── ProgressCallback ──────────────────────────────────────────────────────

/// Control flow signal returned from a progress callback.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ControlFlow {
    Continue,
    Cancel,
}

/// Events emitted during a `run_analysis_with_progress()` call.
#[derive(Debug)]
pub enum ProgressEvent {
    MeshingStarted,
    MeshingDone { n_elements: usize, n_nodes: usize },
    AssemblyStarted,
    AssemblyDone { n_dofs: usize, n_variables: usize },
    Done { load_factor: f64 },
    Cancelled,
}

/// Injectable progress / cancellation interface.
pub trait ProgressCallback: Send + Sync {
    fn on_event(&self, event: ProgressEvent) -> ControlFlow;
}

/// Blanket impl so closures can be used as progress callbacks.
impl<F: Fn(ProgressEvent) -> ControlFlow + Send + Sync> ProgressCallback for F {
    fn on_event(&self, event: ProgressEvent) -> ControlFlow {
        self(event)
    }
}

// ── Public entry points ───────────────────────────────────────────────────

/// Run analysis without progress reporting.
pub fn run_analysis(model: &AnalysisModel) -> Result<SolveResult, AnalysisError> {
    run_analysis_with_progress(model, &|_| ControlFlow::Continue)
}

/// Run analysis with a progress callback.
pub fn run_analysis_with_progress(
    model: &AnalysisModel,
    progress: &dyn ProgressCallback,
) -> Result<SolveResult, AnalysisError> {
    let t0 = Instant::now();

    // ── 1. Validate ───────────────────────────────────────────────────────
    validate_model(model)?;

    // ── 2. Mesh ───────────────────────────────────────────────────────────
    if progress.on_event(ProgressEvent::MeshingStarted) == ControlFlow::Cancel {
        return Err(AnalysisError::Cancelled);
    }

    let mesh_panels: Vec<plato_mesh::geometry::Panel> =
        model.panels.iter().map(|p| p.to_mesh_panel()).collect();
    let mesh = SpadeMesher.triangulate(&mesh_panels, &model.shared_edges, &model.mesh_config)?;

    if progress.on_event(ProgressEvent::MeshingDone {
        n_elements: mesh.elements.len(),
        n_nodes: mesh.nodes.len(),
    }) == ControlFlow::Cancel
    {
        return Err(AnalysisError::Cancelled);
    }

    // ── 3. Resolve BCs ────────────────────────────────────────────────────
    progress.on_event(ProgressEvent::AssemblyStarted);
    let (pinned_w, clamped_theta_mids) = resolve_bcs(&model.supports, &mesh)?;

    // ── 4. DOF map ────────────────────────────────────────────────────────
    let all_nodes = collect_all_displacement_nodes(&mesh);
    let dof_map = DofMap::new(&all_nodes, &pinned_w, &clamped_theta_mids, &mesh);

    if dof_map.total_rows() == 0 {
        return Err(AnalysisError::Assembly(AssemblyError::FullyConstrained));
    }

    // ── 5. Assemble B^T ───────────────────────────────────────────────────
    let bt = assemble_bt(&mesh, &dof_map);
    let n_dof = dof_map.total_rows();

    // ── 6. Load vectors ───────────────────────────────────────────────────
    let (q_ext, q_dead) = build_load_vectors(&model.loads, &mesh, &dof_map, model)?;

    #[cfg(test)]
    {
        let q_sum: f64 = q_ext.iter().sum();
        let q_max: f64 = q_ext.iter().cloned().fold(0.0_f64, f64::max);
        let n_nonzero = q_ext.iter().filter(|&&v| v.abs() > 1e-30).count();
        eprintln!(
            "[diag] n_free_w={} n_rot={} n_dof={} | q_ext: sum={:.4e} max={:.4e} n_nonzero={}",
            dof_map.n_free_w,
            dof_map.n_rot,
            n_dof,
            q_sum,
            q_max,
            n_nonzero
        );
    }

    if q_ext.iter().all(|&v| v == 0.0) {
        return Err(AnalysisError::Validation(
            "no variable loads produce non-zero DOF forces".into(),
        ));
    }

    // ── 7. Per-element criteria ───────────────────────────────────────────
    let unit_factor = model.units.moment_per_length_factor();
    let panel_criteria: HashMap<&str, JohansenCriterion> = model
        .panels
        .iter()
        .map(|p| {
            let m = &p.material;
            let criterion = JohansenCriterion::new(
                m.m_pos_x * unit_factor, // mpx_pos
                m.m_pos_y * unit_factor, // mpy_pos
                m.m_neg_x * unit_factor, // mpx_neg
                m.m_neg_y * unit_factor,
            );
            (p.id.as_str(), criterion)
        })
        .collect();

    let criteria_refs: Vec<&JohansenCriterion> = mesh
        .elements
        .iter()
        .map(|elem| &panel_criteria[elem.panel_id.as_str()])
        .collect();

    // ── 8. Solve ──────────────────────────────────────────────────────────
    let n_e = mesh.elements.len();
    let layout = plato_core::problem::VarLayout { n_e };
    let n_variables = layout.n_vars();

    if progress.on_event(ProgressEvent::AssemblyDone {
        n_dofs: n_dof,
        n_variables,
    }) == ControlFlow::Cancel
    {
        return Err(AnalysisError::Cancelled);
    }

    let core_result = ClarabelProblem {
        bt: &bt,
        q_ext: &q_ext,
        q_dead: &q_dead,
        criteria: &criteria_refs,
        n_e,
    }
    .solve();

    progress.on_event(ProgressEvent::Done {
        load_factor: core_result.load_factor,
    });

    // ── 9. Build public SolveResult ───────────────────────────────────────
    let status = match core_result.status {
        CoreStatus::Solved => SolveStatus::Optimal,
        CoreStatus::AlmostSolved => SolveStatus::MaxIterationsReached,
        CoreStatus::Infeasible => SolveStatus::Infeasible,
        CoreStatus::Failed => SolveStatus::NumericalError,
    };

    let elements: Vec<SolveElement> = mesh
        .elements
        .iter()
        .map(|e| SolveElement {
            id: e.id,
            panel_id: e.panel_id.clone(),
            nodes: [
                e.corners[0],
                e.corners[1],
                e.corners[2],
                e.midpoints[0],
                e.midpoints[1],
                e.midpoints[2],
            ],
        })
        .collect();

    Ok(SolveResult {
        status,
        load_factor: core_result.load_factor,
        units: model.units,
        solve_time_ms: t0.elapsed().as_millis() as u64,
        nodes: mesh.nodes.clone(),
        elements,
        element_moments: core_result.element_moments,
        n_elements: mesh.elements.len(),
        n_nodes: mesh.nodes.len(),
        n_free_dofs: n_dof,
        n_variables,
        solver_iterations: core_result.solver_iterations as usize,
        duality_gap: core_result.duality_gap,
    })
}

// ── BC resolution ─────────────────────────────────────────────────────────

fn resolve_bcs(
    supports: &[Support],
    mesh: &MeshModel,
) -> Result<(HashSet<usize>, HashSet<usize>), AnalysisError> {
    let mut pinned_w: HashSet<usize> = HashSet::new();
    let mut clamped_theta_mids: HashSet<usize> = HashSet::new();

    for support in supports {
        match support {
            Support::SymmetryLine { .. } => {
                return Err(AnalysisError::Validation(
                    "SymmetryLine support is not yet implemented".into(),
                ));
            }
            Support::FreeEdge { .. } => {}
            Support::SimplySupportedEdge { edge } => {
                // Only pin corner nodes (even indices); midpoint nodes (odd indices)
                // remain free so their Kirchhoff shear DOFs carry support reactions.
                for (i, &node) in nodes_for_edge(edge, mesh)?.iter().enumerate() {
                    if i % 2 == 0 {
                        pinned_w.insert(node);
                    }
                }
            }
            Support::ClampedEdge { edge } => {
                let nodes = nodes_for_edge(edge, mesh)?;
                for (i, &node) in nodes.iter().enumerate() {
                    pinned_w.insert(node);
                    // Odd-indexed nodes are midpoints — remove their theta DOFs.
                    if i % 2 == 1 {
                        clamped_theta_mids.insert(node);
                    }
                }
            }
        }
    }

    Ok((pinned_w, clamped_theta_mids))
}

fn nodes_for_edge<'m>(
    edge: &EdgeRef,
    mesh: &'m MeshModel,
) -> Result<&'m Vec<usize>, AnalysisError> {
    mesh.edge_nodes.get(edge).ok_or_else(|| {
        AnalysisError::Validation(format!("edge {edge:?} not found in mesh edge_nodes"))
    })
}

// ── Load vector assembly ──────────────────────────────────────────────────

fn build_load_vectors(
    loads: &[Load],
    mesh: &MeshModel,
    dof_map: &DofMap,
    model: &AnalysisModel,
) -> Result<(Vec<f64>, Vec<f64>), AnalysisError> {
    let n_dof = dof_map.total_rows();
    let mut q_ext = vec![0.0f64; n_dof];
    let mut q_dead = vec![0.0f64; n_dof];

    for load in loads {
        match load {
            Load::AreaLoad {
                panel_id,
                intensity: LoadIntensity::Uniform(p_user),
                load_case,
            } => {
                let p = p_user * model.units.load_intensity_factor();
                let q = target_vec(load_case, &mut q_ext, &mut q_dead);

                for elem in mesh.elements.iter().filter(|e| &e.panel_id == panel_id) {
                    let p0 = mesh.nodes[elem.corners[0]];
                    let p1 = mesh.nodes[elem.corners[1]];
                    let p2 = mesh.nodes[elem.corners[2]];
                    let area = 0.5
                        * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]))
                            .abs();
                    // Nodal loads applied at corner nodes — p·A/3 per corner.
                    // Per the Kirchhoff plate equilibrium: the distributed load
                    // converts to three equal concentrated corner forces Q_D.
                    let f_node = p * area / 3.0;
                    for &corner in &elem.corners {
                        if let Some(row) = dof_map.w_row(corner) {
                            q[row] += f_node;
                        }
                    }
                }
            }
            Load::LineLoad {
                edge,
                intensity: p_user,
                load_case,
            } => {
                let p = p_user * model.units.line_load_factor();
                let nodes = nodes_for_edge(edge, mesh)?;
                let q = target_vec(load_case, &mut q_ext, &mut q_dead);

                // nodes = [c0, m0, c1, m1, c2, ...]; windows of 3 stepped by 2.
                for triplet in nodes.windows(3).step_by(2) {
                    let (ci, mi, ci1) = (triplet[0], triplet[1], triplet[2]);
                    let [x0, y0] = mesh.nodes[ci];
                    let [x1, y1] = mesh.nodes[ci1];
                    let length = ((x1 - x0).powi(2) + (y1 - y0).powi(2)).sqrt();
                    if let Some(row) = dof_map.w_row(ci) {
                        q[row] += p * length / 6.0;
                    }
                    if let Some(row) = dof_map.w_row(mi) {
                        q[row] += p * length * 2.0 / 3.0;
                    }
                    if let Some(row) = dof_map.w_row(ci1) {
                        q[row] += p * length / 6.0;
                    }
                }
            }
        }
    }

    Ok((q_ext, q_dead))
}

fn target_vec<'a>(
    load_case: &LoadCase,
    q_ext: &'a mut Vec<f64>,
    q_dead: &'a mut Vec<f64>,
) -> &'a mut Vec<f64> {
    match load_case {
        LoadCase::Variable => q_ext,
        LoadCase::Permanent => q_dead,
    }
}

// ── Model validation ──────────────────────────────────────────────────────

fn validate_model(model: &AnalysisModel) -> Result<(), AnalysisError> {
    if model.panels.is_empty() {
        return Err(AnalysisError::Validation("no panels defined".into()));
    }

    let panel_ids: HashSet<&str> = model.panels.iter().map(|p| p.id.as_str()).collect();

    for load in &model.loads {
        if let Load::AreaLoad { panel_id, .. } = load {
            if !panel_ids.contains(panel_id.as_str()) {
                return Err(AnalysisError::Validation(format!(
                    "load references unknown panel '{panel_id}'"
                )));
            }
        }
    }

    let has_variable = model.loads.iter().any(|l| {
        matches!(
            l,
            Load::AreaLoad {
                load_case: LoadCase::Variable,
                ..
            } | Load::LineLoad {
                load_case: LoadCase::Variable,
                ..
            }
        )
    });
    if !has_variable {
        return Err(AnalysisError::Validation(
            "at least one Variable load is required".into(),
        ));
    }

    Ok(())
}

// ── Tests ─────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use plato_mesh::geometry::{EdgeRef, Polygon2D};

    use crate::builder::ModelBuilder;
    use crate::loads::{LoadCase, LoadIntensity, MeshDensity};
    use crate::model::{RCSlabMaterial, SolveStatus};
    use crate::supports::Support;
    use crate::units::UnitSystem;

    use super::Load;

    /// Simply-supported 1×1 m isotropic slab, mp = 1 N·m/m, uniform q = 1 N/m² (SI).
    ///
    /// Exact yield-line solution: λ = 24·mp/L² = 24.0.
    /// The SOCP lower-bound converges toward 24.0 from below.
    #[test]
    fn simply_supported_square() {
        let result = ModelBuilder::new("ss_square")
            .with_units(UnitSystem::SI)
            .add_panel("p", Polygon2D::rectangle(0.0, 0.0, 1.0, 1.0))
            .set_material("p", RCSlabMaterial::isotropic(1.0, 1.0))
            .add_support(Support::SimplySupportedEdge {
                edge: EdgeRef::AllEdges {
                    panel_id: "p".into(),
                },
            })
            .add_load(Load::AreaLoad {
                panel_id: "p".into(),
                intensity: LoadIntensity::Uniform(1.0),
                load_case: LoadCase::Variable,
            })
            .set_mesh_density(MeshDensity::MaxElementArea(0.02))
            .solve()
            .expect("should converge");

        assert!(
            matches!(
                result.status,
                SolveStatus::Optimal | SolveStatus::MaxIterationsReached
            ),
            "status: {:?}",
            result.status
        );
        let lf = result.load_factor;
        println!("simply_supported_square λ = {lf:.4}  (exact = 24.0)");
        let util_max = result.element_moments.iter()
            .flat_map(|em| em.yield_utilisation.iter().copied())
            .fold(0.0f64, f64::max);
        let util_mean = {
            let v: Vec<f64> = result.element_moments.iter()
                .flat_map(|em| em.yield_utilisation.iter().copied())
                .collect();
            v.iter().sum::<f64>() / v.len() as f64
        };
        let m_max = result.element_moments.iter()
            .flat_map(|em| em.corner_moments.iter().flat_map(|m| m.iter().copied()))
            .fold(0.0f64, |a, b| a.max(b.abs()));
        println!(
            "  n_el={} n_dof={} | yield_util max={util_max:.4} mean={util_mean:.4} | |m|_max={m_max:.4e}",
            result.n_elements, result.n_free_dofs
        );
        assert!(
            lf > 0.0 && lf.is_finite(),
            "λ must be positive finite, got {lf}"
        );
        assert!(
            lf < 24.0 + 0.5,
            "lower bound must not exceed exact (24.0), got {lf}"
        );
        assert!(lf > 15.0, "coarse mesh should give at least 15, got {lf}");
    }

    /// All-clamped slab: λ must be higher than simply-supported.
    #[test]
    fn clamped_square_lambda_higher_than_ss() {
        let result = ModelBuilder::new("clamped_square")
            .with_units(UnitSystem::SI)
            .add_panel("p", Polygon2D::rectangle(0.0, 0.0, 1.0, 1.0))
            .set_material("p", RCSlabMaterial::isotropic(1.0, 1.0))
            .add_support(Support::ClampedEdge {
                edge: EdgeRef::AllEdges {
                    panel_id: "p".into(),
                },
            })
            .add_load(Load::AreaLoad {
                panel_id: "p".into(),
                intensity: LoadIntensity::Uniform(1.0),
                load_case: LoadCase::Variable,
            })
            .set_mesh_density(MeshDensity::MaxElementArea(0.02))
            .solve()
            .expect("should converge");

        let lf = result.load_factor;
        println!("clamped_square λ = {lf:.4}  (expected > 24.0)");
        assert!(
            lf > 0.0 && lf.is_finite(),
            "λ must be positive finite, got {lf}"
        );
        assert!(
            lf > 24.0,
            "clamped slab must have higher λ than SS (24.0), got {lf}"
        );
    }

    /// Simply-supported 1 m × 8 m beam (plate strip) — SI units.
    ///
    /// Geometry: rectangle (0,0)→(8,1).  CCW vertex order ⇒ edge indices:
    ///   0 : (0,0)→(8,0) — bottom long edge  [free]
    ///   1 : (8,0)→(8,1) — right end x = 8   [simply supported]
    ///   2 : (8,1)→(0,1) — top long edge      [free]
    ///   3 : (0,1)→(0,0) — left  end x = 0   [simply supported]
    ///
    /// Material: isotropic, Mpx+ = Mpx- = 1 N·m/m (SI).
    /// Load    : q = 1 N/m² (variable, uniform, SI).
    ///
    /// Exact yield-line collapse load factor:
    ///   λ_exact = 8·Mpx / (q·L²) = 8·1 / (1·64) = 0.125
    ///
    /// Moment profile at collapse (per unit width):
    ///   mx(x) = λ·q·x·(L−x)/2  — parabola, peak = Mpx at x = L/2.
    ///   my, mxy should be ≈ 0 (beam strip, unidirectional bending).
    #[test]
    fn simply_supported_beam_8m() {
        let l_span = 8.0_f64;
        let mpx = 1.0_f64;
        let q_applied = 1.0_f64;
        let lambda_exact = 8.0 * mpx / (q_applied * l_span * l_span); // 0.125

        let result = ModelBuilder::new("ss_beam_8m")
            .with_units(UnitSystem::SI)
            .add_panel("beam", Polygon2D::rectangle(0.0, 0.0, l_span, 1.0))
            .set_material("beam", RCSlabMaterial::isotropic(mpx, mpx))
            .add_support(Support::SimplySupportedEdge {
                edge: EdgeRef::PolygonEdge {
                    panel_id: "beam".into(),
                    edge_index: 1,
                }, // x = L
            })
            .add_support(Support::SimplySupportedEdge {
                edge: EdgeRef::PolygonEdge {
                    panel_id: "beam".into(),
                    edge_index: 3,
                }, // x = 0
            })
            .add_load(Load::AreaLoad {
                panel_id: "beam".into(),
                intensity: LoadIntensity::Uniform(q_applied),
                load_case: LoadCase::Variable,
            })
            .set_mesh_density(MeshDensity::MaxElementArea(0.1))
            .solve()
            .expect("ss_beam_8m should converge");

        let lambda = result.load_factor;

        println!("\n=== Simply-Supported Beam (1 m × 8 m) ===");
        println!("  Status     : {:?}", result.status);
        println!("  λ computed = {:.6}", lambda);
        println!("  λ exact    = {:.6}  (8·Mpx/(q·L²))", lambda_exact);
        println!(
            "  Error      = {:.2}%",
            (lambda - lambda_exact).abs() / lambda_exact * 100.0
        );
        println!(
            "  Mesh       : {} elements, {} nodes, {} free DOFs",
            result.n_elements, result.n_nodes, result.n_free_dofs
        );

        // ── Moment profile ─────────────────────────────────────────────────
        // For each element: x-centroid and average [mx, my, mxy] across 3 corners.
        let mut profile: Vec<(f64, f64, f64, f64)> = result
            .element_moments
            .iter()
            .enumerate()
            .map(|(i, em)| {
                let elem = &result.elements[i];
                // corners are nodes[0..3]; x-coordinate is node[*][0]
                let x_c = (0..3usize)
                    .map(|k| result.nodes[elem.nodes[k]][0])
                    .sum::<f64>()
                    / 3.0;
                // corner_moments[k] = [mx, my, mxy]
                let mx_avg = em.corner_moments.iter().map(|c| c[0]).sum::<f64>() / 3.0;
                let my_avg = em.corner_moments.iter().map(|c| c[1]).sum::<f64>() / 3.0;
                let mxy_avg = em.corner_moments.iter().map(|c| c[2]).sum::<f64>() / 3.0;
                (x_c, mx_avg, my_avg, mxy_avg)
            })
            .collect();
        profile.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        // Bin into 8 equal zones along the span.
        let n_bins: usize = 8;
        let bin_w = l_span / n_bins as f64;
        let mut bins_mx: Vec<Vec<f64>> = vec![vec![]; n_bins];
        let mut bins_my: Vec<Vec<f64>> = vec![vec![]; n_bins];
        let mut bins_mxy: Vec<Vec<f64>> = vec![vec![]; n_bins];
        for &(x, mx, my, mxy) in &profile {
            let bi = ((x / bin_w) as usize).min(n_bins - 1);
            bins_mx[bi].push(mx);
            bins_my[bi].push(my);
            bins_mxy[bi].push(mxy);
        }

        println!("\n  Moment profile along span (N·m/m, SI):");
        println!(
            "  {:>6}  {:>12}  {:>12}  {:>12}  {:>12}",
            "x [m]", "mx SOCP", "mx theory", "my SOCP", "mxy SOCP"
        );
        for bi in 0..n_bins {
            if bins_mx[bi].is_empty() {
                continue;
            }
            let x_mid = (bi as f64 + 0.5) * bin_w;
            let avg = |v: &Vec<f64>| v.iter().sum::<f64>() / v.len() as f64;
            let mx_socp = avg(&bins_mx[bi]);
            let my_socp = avg(&bins_my[bi]);
            let mxy_socp = avg(&bins_mxy[bi]);
            // Theoretical parabola using the *computed* λ so we compare shapes.
            let mx_theory = lambda * q_applied * x_mid * (l_span - x_mid) / 2.0;
            println!(
                "  {:>6.2}  {:>12.5}  {:>12.5}  {:>12.5}  {:>12.5}",
                x_mid, mx_socp, mx_theory, my_socp, mxy_socp
            );
        }

        // Peak |mx| at midspan.
        let midspan_mx_abs = bins_mx[n_bins / 2 - 1]
            .iter()
            .chain(bins_mx[n_bins / 2].iter())
            .map(|v| v.abs())
            .fold(0.0_f64, f64::max);
        println!("\n  Peak |mx| (midspan)  = {midspan_mx_abs:.5} N·m/m");
        println!("  Mpx capacity         = {mpx:.5} N·m/m");

        // ── Assertions ─────────────────────────────────────────────────────
        assert!(
            matches!(
                result.status,
                SolveStatus::Optimal | SolveStatus::MaxIterationsReached
            ),
            "status: {:?}",
            result.status
        );
        assert!(
            lambda > 0.0 && lambda.is_finite(),
            "λ must be positive finite, got {lambda}"
        );
        // Lower-bound SOCP must not exceed exact (within numerical noise).
        assert!(
            lambda <= lambda_exact + 1e-4,
            "lower bound exceeded exact: {lambda:.6} > {lambda_exact:.6}"
        );
        // Should recover at least 90% of exact for this mesh density.
        assert!(
            lambda >= lambda_exact * 0.90,
            "λ too far below exact: {lambda:.6} vs {:.6}",
            lambda_exact * 0.90
        );
    }

    /// Permanent dead load reduces the effective variable-load capacity.
    #[test]
    fn dead_load_reduces_lambda() {
        let base = ModelBuilder::new("base")
            .with_units(UnitSystem::SI)
            .add_panel("p", Polygon2D::rectangle(0.0, 0.0, 1.0, 1.0))
            .set_material("p", RCSlabMaterial::isotropic(1.0, 1.0))
            .add_support(Support::SimplySupportedEdge {
                edge: EdgeRef::AllEdges {
                    panel_id: "p".into(),
                },
            })
            .add_load(Load::AreaLoad {
                panel_id: "p".into(),
                intensity: LoadIntensity::Uniform(1.0),
                load_case: LoadCase::Variable,
            })
            .set_mesh_density(MeshDensity::MaxElementArea(0.05))
            .solve()
            .expect("base should converge");

        let with_dead = ModelBuilder::new("with_dead")
            .with_units(UnitSystem::SI)
            .add_panel("p", Polygon2D::rectangle(0.0, 0.0, 1.0, 1.0))
            .set_material("p", RCSlabMaterial::isotropic(1.0, 1.0))
            .add_support(Support::SimplySupportedEdge {
                edge: EdgeRef::AllEdges {
                    panel_id: "p".into(),
                },
            })
            .add_load(Load::AreaLoad {
                panel_id: "p".into(),
                intensity: LoadIntensity::Uniform(1.0),
                load_case: LoadCase::Variable,
            })
            .add_load(Load::AreaLoad {
                panel_id: "p".into(),
                intensity: LoadIntensity::Uniform(5.0),
                load_case: LoadCase::Permanent,
            })
            .set_mesh_density(MeshDensity::MaxElementArea(0.05))
            .solve()
            .expect("with_dead should converge");

        println!(
            "dead_load: base λ = {:.4}, with_dead λ = {:.4}",
            base.load_factor, with_dead.load_factor
        );
        assert!(
            with_dead.load_factor < base.load_factor,
            "permanent load must reduce effective λ: base={:.4}, with_dead={:.4}",
            base.load_factor,
            with_dead.load_factor
        );
    }
}
