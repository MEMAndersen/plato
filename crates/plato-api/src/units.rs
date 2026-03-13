use serde::{Deserialize, Serialize};

/// Controls input/output unit conversion. Full implementation in Phase 4.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum UnitSystem {
    /// N, m, Pa, N·m/m
    SI,
    /// kN, m, kPa, kN·m/m  (standard for RC design)
    KiloNewtonMetre,
    /// kN, mm, MPa, kN·mm/mm
    KiloNewtonMillimetre,
}
