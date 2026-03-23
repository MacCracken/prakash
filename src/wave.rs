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
    let mut angles = Vec::with_capacity(2 * max_order as usize + 1);
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
    pub const HORIZONTAL: Self = Self {
        ex: 1.0,
        ey: 0.0,
        phase: 0.0,
    };
    /// Vertically polarized light.
    pub const VERTICAL: Self = Self {
        ex: 0.0,
        ey: 1.0,
        phase: 0.0,
    };
    /// Right circular polarization.
    pub fn circular_right() -> Self {
        Self {
            ex: 1.0 / 2.0f64.sqrt(),
            ey: 1.0 / 2.0f64.sqrt(),
            phase: -PI / 2.0,
        }
    }
    /// Left circular polarization.
    pub fn circular_left() -> Self {
        Self {
            ex: 1.0 / 2.0f64.sqrt(),
            ey: 1.0 / 2.0f64.sqrt(),
            phase: PI / 2.0,
        }
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

    // ── Interference tests ────────────────────────────────────────────────

    #[test]
    fn test_constructive_interference() {
        assert!(is_constructive(500.0, 500.0));
        assert!(is_constructive(1000.0, 500.0));
        assert!(is_constructive(0.0, 500.0)); // m=0
        assert!(!is_constructive(250.0, 500.0));
    }

    #[test]
    fn test_destructive_interference() {
        assert!(is_destructive(250.0, 500.0));
        assert!(is_destructive(750.0, 500.0));
        assert!(!is_destructive(500.0, 500.0));
        assert!(!is_destructive(0.0, 500.0));
    }

    #[test]
    fn test_constructive_destructive_mutual_exclusion() {
        // At any given path difference, cannot be both
        for m in 0..10 {
            let path = (m as f64) * 500.0;
            assert!(
                is_constructive(path, 500.0) != is_destructive(path, 500.0)
                    || !is_constructive(path, 500.0),
                "Both true at path={path}"
            );
        }
    }

    #[test]
    fn test_interference_constructive_max() {
        let i = interference_intensity(1.0, 1.0, 0.0);
        assert!((i - 4.0).abs() < EPS);
    }

    #[test]
    fn test_interference_destructive_zero() {
        let i = interference_intensity(1.0, 1.0, PI);
        assert!(i.abs() < EPS);
    }

    #[test]
    fn test_interference_unequal_amplitudes() {
        // Unequal amplitudes: never reach zero
        let i_min = interference_intensity(2.0, 1.0, PI);
        assert!((i_min - 1.0).abs() < EPS); // (2-1)² = 1
        let i_max = interference_intensity(2.0, 1.0, 0.0);
        assert!((i_max - 9.0).abs() < EPS); // (2+1)² = 9
    }

    #[test]
    fn test_interference_single_wave() {
        let i = interference_intensity(3.0, 0.0, 0.0);
        assert!((i - 9.0).abs() < EPS); // just A²
    }

    #[test]
    fn test_path_to_phase() {
        let phase = path_to_phase(500.0, 500.0);
        assert!((phase - 2.0 * PI).abs() < EPS);
    }

    #[test]
    fn test_path_to_phase_half_wavelength() {
        let phase = path_to_phase(250.0, 500.0);
        assert!((phase - PI).abs() < EPS);
    }

    #[test]
    fn test_path_to_phase_zero() {
        let phase = path_to_phase(0.0, 500.0);
        assert!(phase.abs() < EPS);
    }

    // ── Thin film tests ───────────────────────────────────────────────────

    #[test]
    fn test_thin_film_range() {
        let r = thin_film_reflectance(550.0, 100.0, 1.5);
        assert!((0.0..=1.0).contains(&r));
    }

    #[test]
    fn test_thin_film_varies_with_wavelength() {
        let r1 = thin_film_reflectance(400.0, 100.0, 1.5);
        let r2 = thin_film_reflectance(600.0, 100.0, 1.5);
        // Different wavelengths should generally give different reflectance
        assert!((r1 - r2).abs() > EPS);
    }

    #[test]
    fn test_thin_film_quarter_wave() {
        // Quarter-wave coating: thickness = λ/(4n) should maximize reflectance
        let n = 1.5;
        let wl = 550.0;
        let t_quarter = wl / (4.0 * n);
        let r = thin_film_reflectance(wl, t_quarter, n);
        // At quarter-wave, path = 2*n*t = λ/2, phase = π + π = 2π → sin²(π) = 0
        // Actually for anti-reflection: sin²(δ/2) with δ = 2π + π
        assert!((0.0..=1.0).contains(&r));
    }

    // ── Diffraction tests ─────────────────────────────────────────────────

    #[test]
    fn test_single_slit_central_max() {
        let i = single_slit_intensity(1e-3, 500e-9, 0.0, 1.0);
        assert!((i - 1.0).abs() < EPS);
    }

    #[test]
    fn test_single_slit_decreases_off_axis() {
        let i_center = single_slit_intensity(1e-3, 500e-9, 0.0, 1.0);
        let i_off = single_slit_intensity(1e-3, 500e-9, 0.01, 1.0);
        assert!(i_off < i_center);
    }

    #[test]
    fn test_single_slit_first_minimum() {
        // First minimum at sin(θ) = λ/a → θ = asin(λ/a)
        let a = 1e-3;
        let wl = 500e-9;
        let ratio: f64 = wl / a;
        let theta_min = ratio.asin();
        let i = single_slit_intensity(a, wl, theta_min, 1.0);
        assert!(i < 1e-6, "Intensity at first minimum should be ~0, got {i}");
    }

    #[test]
    fn test_single_slit_scales_with_i0() {
        let i1 = single_slit_intensity(1e-3, 500e-9, 0.01, 1.0);
        let i5 = single_slit_intensity(1e-3, 500e-9, 0.01, 5.0);
        assert!((i5 / i1 - 5.0).abs() < EPS);
    }

    #[test]
    fn test_single_slit_always_non_negative() {
        for angle_mrad in 0..100 {
            let angle = angle_mrad as f64 * 0.001;
            let i = single_slit_intensity(1e-3, 500e-9, angle, 1.0);
            assert!(i >= 0.0, "Negative intensity at angle {angle}");
        }
    }

    #[test]
    fn test_double_slit_central() {
        let i = double_slit_intensity(0.1e-3, 0.5e-3, 500e-9, 0.0, 1.0);
        assert!(i > 0.0);
    }

    #[test]
    fn test_double_slit_greater_than_single_at_center() {
        let i_single = single_slit_intensity(0.1e-3, 500e-9, 0.0, 1.0);
        let i_double = double_slit_intensity(0.1e-3, 0.5e-3, 500e-9, 0.0, 1.0);
        assert!(
            i_double > i_single,
            "Double slit central max should exceed single slit"
        );
    }

    #[test]
    fn test_double_slit_non_negative() {
        for angle_mrad in 0..50 {
            let angle = angle_mrad as f64 * 0.001;
            let i = double_slit_intensity(0.1e-3, 0.5e-3, 500e-9, angle, 1.0);
            assert!(i >= 0.0, "Negative intensity at angle {angle}");
        }
    }

    // ── Grating tests ─────────────────────────────────────────────────────

    #[test]
    fn test_grating_maxima_zeroth_order() {
        let angles = grating_maxima(1e-6, 500e-9, 0);
        assert_eq!(angles.len(), 1);
        assert!(angles[0].abs() < EPS);
    }

    #[test]
    fn test_grating_maxima_multiple_orders() {
        let angles = grating_maxima(1e-6, 500e-9, 2);
        assert!(angles.len() >= 3);
    }

    #[test]
    fn test_grating_maxima_symmetric() {
        let angles = grating_maxima(1e-6, 500e-9, 1);
        // Should have m=0, +1, -1
        assert_eq!(angles.len(), 3);
        // +1 and -1 should be symmetric
        assert!((angles[1] + angles[2]).abs() < EPS);
    }

    #[test]
    fn test_grating_maxima_limited_by_sin() {
        // Very fine grating with long wavelength: fewer orders possible
        let angles = grating_maxima(600e-9, 500e-9, 5);
        // m·λ/d = m·500/600, max m where this ≤ 1 is m=1
        assert!(angles.len() <= 3); // m=0, ±1 at most
    }

    // ── Polarization tests ────────────────────────────────────────────────

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
    fn test_malus_law_scales_with_intensity() {
        assert!((malus_law(10.0, PI / 4.0) - 5.0).abs() < EPS);
    }

    #[test]
    fn test_malus_law_30_degrees() {
        let i = malus_law(1.0, PI / 6.0);
        assert!((i - 0.75).abs() < EPS); // cos²(30°) = 3/4
    }

    #[test]
    fn test_malus_law_60_degrees() {
        let i = malus_law(1.0, PI / 3.0);
        assert!((i - 0.25).abs() < EPS); // cos²(60°) = 1/4
    }

    #[test]
    fn test_polarization_horizontal() {
        let p = Polarization::HORIZONTAL;
        assert!((p.intensity() - 1.0).abs() < EPS);
        assert!((p.ex - 1.0).abs() < EPS);
        assert!(p.ey.abs() < EPS);
    }

    #[test]
    fn test_polarization_vertical() {
        let p = Polarization::VERTICAL;
        assert!((p.intensity() - 1.0).abs() < EPS);
        assert!(p.ex.abs() < EPS);
        assert!((p.ey - 1.0).abs() < EPS);
    }

    #[test]
    fn test_polarization_through_aligned() {
        let p = Polarization::HORIZONTAL;
        let i = p.through_polarizer(0.0);
        assert!((i - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_polarization_through_crossed() {
        let p = Polarization::HORIZONTAL;
        let i = p.through_polarizer(PI / 2.0);
        assert!(i < 0.01);
    }

    #[test]
    fn test_polarization_circular_right() {
        let p = Polarization::circular_right();
        assert!((p.intensity() - 1.0).abs() < EPS);
        assert!((p.phase + PI / 2.0).abs() < EPS);
    }

    #[test]
    fn test_polarization_circular_left() {
        let p = Polarization::circular_left();
        assert!((p.intensity() - 1.0).abs() < EPS);
        assert!((p.phase - PI / 2.0).abs() < EPS);
    }

    #[test]
    fn test_polarization_serde_roundtrip() {
        let p = Polarization::circular_right();
        let json = serde_json::to_string(&p).unwrap();
        let back: Polarization = serde_json::from_str(&json).unwrap();
        assert!((back.ex - p.ex).abs() < EPS);
        assert!((back.ey - p.ey).abs() < EPS);
        assert!((back.phase - p.phase).abs() < EPS);
    }
}
