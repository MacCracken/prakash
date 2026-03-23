//! Wave optics — interference, diffraction, polarization.

use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

// ── Interference ────────────────────────────────────────────────────────────

/// Two-wave interference intensity.
///
/// Given amplitudes A1, A2 and phase difference δ (radians):
/// I = A1² + A2² + 2·A1·A2·cos(δ)
#[inline]
pub fn interference_intensity(a1: f64, a2: f64, phase_diff: f64) -> f64 {
    a1 * a1 + a2 * a2 + 2.0 * a1 * a2 * phase_diff.cos()
}

/// Constructive interference condition: path difference = m·λ.
#[inline]
pub fn is_constructive(path_diff: f64, wavelength: f64) -> bool {
    let m = path_diff / wavelength;
    (m - m.round()).abs() < 0.01
}

/// Destructive interference condition: path difference = (m + 0.5)·λ.
#[inline]
pub fn is_destructive(path_diff: f64, wavelength: f64) -> bool {
    let m = path_diff / wavelength - 0.5;
    (m - m.round()).abs() < 0.01
}

/// Phase difference from path difference and wavelength.
/// δ = 2π·Δ/λ
#[inline]
pub fn path_to_phase(path_diff: f64, wavelength: f64) -> f64 {
    2.0 * PI * path_diff / wavelength
}

// ── Thin Film Interference ──────────────────────────────────────────────────

/// Thin film interference: reflected intensity for a thin film of thickness `d`,
/// refractive index `n_film`, at near-normal incidence.
///
/// Returns intensity as a fraction of incident (0.0–1.0 approximate).
pub fn thin_film_reflectance(wavelength_nm: f64, thickness_nm: f64, n_film: f64) -> f64 {
    let path = 2.0 * n_film * thickness_nm;
    let phase = 2.0 * PI * path / wavelength_nm + PI; // +π for phase change on reflection
    // Simplified: R ∝ sin²(δ/2) for low-reflectance films
    let half_phase = phase / 2.0;
    half_phase.sin().powi(2)
}

// ── Diffraction ─────────────────────────────────────────────────────────────

/// Single-slit diffraction: intensity at angle θ.
///
/// I(θ) = I0 · (sin(β)/β)² where β = π·a·sin(θ)/λ
/// `slit_width` and `wavelength` in same units.
pub fn single_slit_intensity(slit_width: f64, wavelength: f64, angle: f64, i0: f64) -> f64 {
    let beta = PI * slit_width * angle.sin() / wavelength;
    if beta.abs() < 1e-10 {
        return i0; // central maximum
    }
    let sinc = beta.sin() / beta;
    i0 * sinc * sinc
}

/// Double-slit diffraction: intensity at angle θ.
///
/// Combines single-slit envelope with two-slit interference.
/// `slit_width`: width of each slit, `slit_spacing`: center-to-center distance.
pub fn double_slit_intensity(
    slit_width: f64,
    slit_spacing: f64,
    wavelength: f64,
    angle: f64,
    i0: f64,
) -> f64 {
    // Single-slit envelope
    let envelope = single_slit_intensity(slit_width, wavelength, angle, 1.0);
    // Two-slit interference
    let delta = PI * slit_spacing * angle.sin() / wavelength;
    let interference = delta.cos().powi(2);
    i0 * envelope * interference * 4.0 // 4x for coherent double slit
}

/// Diffraction grating: angular positions of maxima.
///
/// d·sin(θ) = m·λ → θ = asin(m·λ/d)
/// Returns angles for orders m = 0, ±1, ±2, ... up to `max_order`.
pub fn grating_maxima(grating_spacing: f64, wavelength: f64, max_order: u32) -> Vec<f64> {
    let mut angles = Vec::new();
    for m in 0..=max_order {
        let sin_theta = (m as f64) * wavelength / grating_spacing;
        if sin_theta.abs() <= 1.0 {
            let angle = sin_theta.asin();
            if m == 0 {
                angles.push(angle);
            } else {
                angles.push(angle);
                angles.push(-angle);
            }
        }
    }
    angles
}

// ── Polarization ────────────────────────────────────────────────────────────

/// Polarization state (Jones vector components).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Polarization {
    /// Amplitude of horizontal component.
    pub ex: f64,
    /// Amplitude of vertical component.
    pub ey: f64,
    /// Phase difference between components (radians).
    pub phase: f64,
}

impl Polarization {
    /// Horizontally polarized light.
    pub const HORIZONTAL: Self = Self { ex: 1.0, ey: 0.0, phase: 0.0 };
    /// Vertically polarized light.
    pub const VERTICAL: Self = Self { ex: 0.0, ey: 1.0, phase: 0.0 };
    /// Right circular polarization.
    pub fn circular_right() -> Self {
        Self { ex: 1.0 / 2.0f64.sqrt(), ey: 1.0 / 2.0f64.sqrt(), phase: -PI / 2.0 }
    }
    /// Left circular polarization.
    pub fn circular_left() -> Self {
        Self { ex: 1.0 / 2.0f64.sqrt(), ey: 1.0 / 2.0f64.sqrt(), phase: PI / 2.0 }
    }

    /// Intensity after passing through a linear polarizer at angle θ.
    /// Malus's law: I = I0 · cos²(θ - polarization_angle)
    #[inline]
    pub fn through_polarizer(&self, polarizer_angle: f64) -> f64 {
        let pol_angle = self.ey.atan2(self.ex);
        let diff = polarizer_angle - pol_angle;
        diff.cos().powi(2) * (self.ex * self.ex + self.ey * self.ey)
    }

    /// Total intensity.
    #[inline]
    pub fn intensity(&self) -> f64 {
        self.ex * self.ex + self.ey * self.ey
    }
}

/// Malus's law: intensity after a polarizer at angle θ relative to polarization.
///
/// I = I0 · cos²(θ)
#[inline]
pub fn malus_law(intensity: f64, angle: f64) -> f64 {
    intensity * angle.cos().powi(2)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-6;

    #[test]
    fn test_constructive_interference() {
        assert!(is_constructive(500.0, 500.0)); // m=1
        assert!(is_constructive(1000.0, 500.0)); // m=2
        assert!(!is_constructive(250.0, 500.0)); // m=0.5 → destructive
    }

    #[test]
    fn test_destructive_interference() {
        assert!(is_destructive(250.0, 500.0)); // m=0.5
        assert!(is_destructive(750.0, 500.0)); // m=1.5
        assert!(!is_destructive(500.0, 500.0)); // m=1 → constructive
    }

    #[test]
    fn test_interference_constructive_max() {
        // Same amplitude, zero phase diff → 4x single intensity
        let i = interference_intensity(1.0, 1.0, 0.0);
        assert!((i - 4.0).abs() < EPS);
    }

    #[test]
    fn test_interference_destructive_zero() {
        // Same amplitude, π phase diff → zero intensity
        let i = interference_intensity(1.0, 1.0, PI);
        assert!(i.abs() < EPS);
    }

    #[test]
    fn test_path_to_phase() {
        let phase = path_to_phase(500.0, 500.0);
        assert!((phase - 2.0 * PI).abs() < EPS);
    }

    #[test]
    fn test_single_slit_central_max() {
        let i = single_slit_intensity(1e-3, 500e-9, 0.0, 1.0);
        assert!((i - 1.0).abs() < EPS); // central maximum = I0
    }

    #[test]
    fn test_single_slit_decreases_off_axis() {
        let i_center = single_slit_intensity(1e-3, 500e-9, 0.0, 1.0);
        let i_off = single_slit_intensity(1e-3, 500e-9, 0.01, 1.0);
        assert!(i_off < i_center);
    }

    #[test]
    fn test_double_slit_central() {
        let i = double_slit_intensity(0.1e-3, 0.5e-3, 500e-9, 0.0, 1.0);
        assert!(i > 0.0); // central max
    }

    #[test]
    fn test_grating_maxima_zeroth_order() {
        let angles = grating_maxima(1e-6, 500e-9, 0);
        assert_eq!(angles.len(), 1); // just m=0
        assert!(angles[0].abs() < EPS);
    }

    #[test]
    fn test_grating_maxima_multiple_orders() {
        let angles = grating_maxima(1e-6, 500e-9, 2);
        assert!(angles.len() >= 3); // m=0, ±1, possibly ±2
    }

    #[test]
    fn test_thin_film() {
        let r = thin_film_reflectance(550.0, 100.0, 1.5);
        assert!(r >= 0.0 && r <= 1.0);
    }

    #[test]
    fn test_malus_law_aligned() {
        assert!((malus_law(1.0, 0.0) - 1.0).abs() < EPS);
    }

    #[test]
    fn test_malus_law_crossed() {
        assert!(malus_law(1.0, PI / 2.0).abs() < EPS);
    }

    #[test]
    fn test_malus_law_45() {
        assert!((malus_law(1.0, PI / 4.0) - 0.5).abs() < EPS);
    }

    #[test]
    fn test_polarization_horizontal() {
        let p = Polarization::HORIZONTAL;
        assert!((p.intensity() - 1.0).abs() < EPS);
    }

    #[test]
    fn test_polarization_through_aligned() {
        let p = Polarization::HORIZONTAL;
        let i = p.through_polarizer(0.0);
        assert!((i - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_polarization_circular() {
        let p = Polarization::circular_right();
        assert!((p.intensity() - 1.0).abs() < EPS);
    }
}
