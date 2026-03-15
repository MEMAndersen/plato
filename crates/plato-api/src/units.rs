use serde::{Deserialize, Serialize};

/// Controls input/output unit conversion.
///
/// Internal storage is always in base SI units (N, m). These methods return
/// the factor to multiply user-facing values by to obtain SI values.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum UnitSystem {
    /// N, m, Pa, N·m/m
    SI,
    /// kN, m, kPa, kN·m/m  (standard for RC design)
    KiloNewtonMetre,
    /// kN, mm, MPa, kN·mm/mm
    KiloNewtonMillimetre,
}

impl UnitSystem {
    /// Factor to convert user moment-per-length to N·m/m.
    ///
    /// - SI: 1 N·m/m = 1 N·m/m           → ×1.0
    /// - kN·m/m = 1000 N·m/m              → ×1000
    /// - kN·mm/mm = kN (mm/mm cancel)     → ×1000 (= 1 kN = 1000 N·m/m)
    pub fn moment_per_length_factor(self) -> f64 {
        match self {
            UnitSystem::SI => 1.0,
            UnitSystem::KiloNewtonMetre => 1_000.0,
            UnitSystem::KiloNewtonMillimetre => 1_000.0,
        }
    }

    /// Factor to convert user load intensity to N/m².
    ///
    /// - SI: 1 N/m² → ×1.0
    /// - 1 kN/m² = 1000 N/m²              → ×1000
    /// - 1 kN/mm² = 1e3 N / (1e-3 m)²    → ×1e9
    pub fn load_intensity_factor(self) -> f64 {
        match self {
            UnitSystem::SI => 1.0,
            UnitSystem::KiloNewtonMetre => 1_000.0,
            UnitSystem::KiloNewtonMillimetre => 1_000.0 * 1_000.0 * 1_000.0, // 1e9
        }
    }

    /// Factor to convert user line-load intensity to N/m.
    ///
    /// - SI: 1 N/m → ×1.0
    /// - 1 kN/m = 1000 N/m                → ×1000
    /// - 1 kN/mm = 1e3 N / 1e-3 m        → ×1e6
    pub fn line_load_factor(self) -> f64 {
        match self {
            UnitSystem::SI => 1.0,
            UnitSystem::KiloNewtonMetre => 1_000.0,
            UnitSystem::KiloNewtonMillimetre => 1_000.0 * 1_000.0, // 1e6
        }
    }
}
