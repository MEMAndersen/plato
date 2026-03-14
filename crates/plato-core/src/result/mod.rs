/// Status returned by the Clarabel solver.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolveStatus {
    /// Optimal solution found.
    Solved,
    /// Optimal solution found at reduced accuracy.
    AlmostSolved,
    /// Problem is infeasible.
    Infeasible,
    /// Solver failed (numerical issues, iteration limit, etc.).
    Failed,
}

/// Moment field and yield utilisation for one triangular element.
#[derive(Debug, Clone)]
pub struct ElementMoments {
    pub element_id: usize,
    /// Moments [mx, my, mxy] at each of the 3 corners (corner_k = 0..3).
    pub corner_moments: [[f64; 3]; 3],
    /// Johansen yield utilisation at each corner (0 = unyielded, 1 = yielded).
    pub yield_utilisation: [f64; 3],
}

/// Outcome of a limit-analysis solve.
#[derive(Debug)]
pub struct SolveResult {
    pub status: SolveStatus,
    /// Plastic load factor λ (multiplier on the applied variable load).
    pub load_factor: f64,
    pub element_moments: Vec<ElementMoments>,
}
