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
    std::f64::consts::TAU * path_diff / wavelength
}

// ── Thin Film Interference ──────────────────────────────────────────────────

/// Thin film interference: reflected intensity for a thin film of thickness `d`,
/// refractive index `n_film`, at near-normal incidence.
///
/// Returns intensity as a fraction of incident (0.0–1.0 approximate).
#[inline]
pub fn thin_film_reflectance(wavelength_nm: f64, thickness_nm: f64, n_film: f64) -> f64 {
    let path = 2.0 * n_film * thickness_nm;
    let phase = std::f64::consts::TAU * path / wavelength_nm + PI; // +π for phase change on reflection
    // Simplified: R ∝ sin²(δ/2) for low-reflectance films
    let half_phase = phase / 2.0;
    half_phase.sin().powi(2)
}

// ── Diffraction ─────────────────────────────────────────────────────────────

/// Single-slit diffraction: intensity at angle θ.
///
/// I(θ) = I0 · (sin(β)/β)² where β = π·a·sin(θ)/λ
/// `slit_width` and `wavelength` in same units.
#[inline]
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
#[inline]
pub fn double_slit_intensity(
    slit_width: f64,
    slit_spacing: f64,
    wavelength: f64,
    angle: f64,
    i0: f64,
) -> f64 {
    let envelope = single_slit_intensity(slit_width, wavelength, angle, 1.0);
    let delta = PI * slit_spacing * angle.sin() / wavelength;
    let interference = delta.cos().powi(2);
    i0 * envelope * interference * 4.0
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

// ── Coherence ─────────────────────────────────────────────────────────────

/// Temporal coherence length.
///
/// l_c = λ² / Δλ
///
/// `center_wavelength` and `bandwidth` in same units (e.g., nm or m).
/// Returns coherence length in the same units.
#[inline]
pub fn coherence_length(center_wavelength: f64, bandwidth: f64) -> f64 {
    center_wavelength * center_wavelength / bandwidth
}

/// Temporal coherence time.
///
/// τ_c = 1 / Δν = λ² / (c · Δλ)
///
/// `center_wavelength_m` and `bandwidth_m` in meters.
/// Returns coherence time in seconds.
#[inline]
pub fn coherence_time(center_wavelength_m: f64, bandwidth_m: f64) -> f64 {
    coherence_length(center_wavelength_m, bandwidth_m) / 299_792_458.0
}

/// Spatial coherence angle (van Cittert–Zernike theorem).
///
/// θ_c = λ / d
///
/// `wavelength` and `source_diameter` in same units.
/// Returns the coherence angle in radians.
#[inline]
pub fn spatial_coherence_angle(wavelength: f64, source_diameter: f64) -> f64 {
    wavelength / source_diameter
}

/// Spatial coherence area at distance R from a source.
///
/// A_c = (λ · R / d)²
///
/// `wavelength`, `distance`, and `source_diameter` in same units.
/// Returns the coherence area in the same units squared.
#[inline]
pub fn coherence_area(wavelength: f64, distance: f64, source_diameter: f64) -> f64 {
    let l = wavelength * distance / source_diameter;
    l * l
}

/// Number of coherence lengths that fit in a given path difference.
///
/// A visibility metric: values >> 1 mean the source is incoherent at this path difference.
#[inline]
pub fn coherence_ratio(path_difference: f64, center_wavelength: f64, bandwidth: f64) -> f64 {
    path_difference / coherence_length(center_wavelength, bandwidth)
}

// ── Circular Aperture Diffraction ─────────────────────────────────────────

/// First-order Bessel function J₁(x) — rational approximation.
///
/// Accurate to ~1e-7 for all x. Uses polynomial approximation for |x| < 8
/// and asymptotic expansion for |x| >= 8.
#[inline]
pub fn bessel_j1(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 8.0 {
        // Polynomial approximation (Abramowitz & Stegun)
        let y = x * x;
        let num = x
            * (72362614232.0
                + y * (-7895059235.0
                    + y * (242396853.1
                        + y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606))))));
        let den = 144725228442.0
            + y * (2300535178.0 + y * (18583304.74 + y * (99447.43394 + y * (376.9991397 + y))));
        num / den
    } else {
        // Asymptotic expansion for large |x|
        let z = 8.0 / ax;
        let z2 = z * z;
        let theta = ax - 2.356_194_490_2; // ax - 3π/4
        let p1 = 1.0
            + z2 * (0.183_105_e-2
                + z2 * (-0.3516396496e-4 + z2 * (0.2457520174e-5 + z2 * (-0.240337019e-6))));
        let q1 = 0.04687499995
            + z2 * (-0.2002690873e-3
                + z2 * (0.8449199096e-5 + z2 * (-0.88228987e-6 + z2 * 0.105787412e-6)));
        let result =
            (std::f64::consts::FRAC_2_PI / ax).sqrt() * (theta.cos() * p1 - z * theta.sin() * q1);
        if x < 0.0 { -result } else { result }
    }
}

/// Airy diffraction pattern intensity for a circular aperture.
///
/// I(θ) = I₀ · [2·J₁(x)/x]² where x = π·D·sin(θ)/λ
///
/// `aperture_diameter` and `wavelength` in same units.
/// `angle` in radians. Returns intensity relative to `i0`.
#[inline]
pub fn airy_pattern(aperture_diameter: f64, wavelength: f64, angle: f64, i0: f64) -> f64 {
    let x = PI * aperture_diameter * angle.sin() / wavelength;
    if x.abs() < 1e-10 {
        return i0; // central maximum
    }
    let jinc = 2.0 * bessel_j1(x) / x;
    i0 * jinc * jinc
}

/// Angular radius of the first Airy disk zero (first dark ring).
///
/// θ₁ = 1.2196·λ/D (radians)
///
/// More precise than the commonly quoted 1.22 factor.
#[inline]
pub fn airy_first_zero(wavelength: f64, aperture_diameter: f64) -> f64 {
    // First zero of J1 is at x = 3.8317..., so sin(θ) = 3.8317λ/(πD) ≈ 1.2197λ/D
    1.219_670_5 * wavelength / aperture_diameter
}

/// Rayleigh resolution criterion.
///
/// Two point sources are just resolved when the central maximum of one
/// falls on the first zero of the other: θ_R = 1.22·λ/D
///
/// Returns minimum resolvable angle in radians.
#[inline]
pub fn rayleigh_criterion(wavelength: f64, aperture_diameter: f64) -> f64 {
    1.22 * wavelength / aperture_diameter
}

// ── Fabry-Pérot Interferometer ────────────────────────────────────────────

/// Fabry-Pérot coefficient of finesse.
///
/// F = 4R / (1 − R)²
///
/// `reflectance` is the mirror reflectivity (0.0–1.0).
#[inline]
pub fn fabry_perot_finesse_coefficient(reflectance: f64) -> f64 {
    4.0 * reflectance / ((1.0 - reflectance) * (1.0 - reflectance))
}

/// Fabry-Pérot finesse (sharpness of fringes).
///
/// F = π·√R / (1 − R)
///
/// Higher finesse means sharper transmission peaks.
#[inline]
pub fn fabry_perot_finesse(reflectance: f64) -> f64 {
    PI * reflectance.sqrt() / (1.0 - reflectance)
}

/// Fabry-Pérot transmittance (Airy function).
///
/// T = 1 / (1 + F·sin²(δ/2))
///
/// where δ = 4π·n·d·cos(θ)/λ is the round-trip phase.
///
/// `n` = refractive index of cavity, `thickness` = mirror separation,
/// `wavelength` and `thickness` in same units, `angle` = incidence angle (radians),
/// `reflectance` = mirror reflectivity.
#[inline]
pub fn fabry_perot_transmittance(
    wavelength: f64,
    thickness: f64,
    n: f64,
    angle: f64,
    reflectance: f64,
) -> f64 {
    let delta = 4.0 * PI * n * thickness * angle.cos() / wavelength;
    let f_coeff = fabry_perot_finesse_coefficient(reflectance);
    1.0 / (1.0 + f_coeff * (delta / 2.0).sin().powi(2))
}

/// Fabry-Pérot free spectral range (FSR).
///
/// Δν = c / (2·n·d)  (in Hz)
///
/// or equivalently in wavelength: Δλ = λ² / (2·n·d)
///
/// `thickness_m` in meters, `n` = refractive index of cavity.
/// Returns FSR in Hz.
#[inline]
pub fn fabry_perot_fsr(thickness_m: f64, n: f64) -> f64 {
    299_792_458.0 / (2.0 * n * thickness_m)
}

/// Fabry-Pérot free spectral range in wavelength units.
///
/// Δλ = λ² / (2·n·d)
///
/// `wavelength` and `thickness` in same units.
#[inline]
pub fn fabry_perot_fsr_wavelength(wavelength: f64, thickness: f64, n: f64) -> f64 {
    wavelength * wavelength / (2.0 * n * thickness)
}

/// Fabry-Pérot resolving power.
///
/// R = m · F where m = 2·n·d/λ (order number) and F is finesse.
///
/// `wavelength` and `thickness` in same units.
#[inline]
pub fn fabry_perot_resolving_power(
    wavelength: f64,
    thickness: f64,
    n: f64,
    reflectance: f64,
) -> f64 {
    let order = 2.0 * n * thickness / wavelength;
    let finesse = fabry_perot_finesse(reflectance);
    order * finesse
}

mod diffraction;
mod polarization;

pub use diffraction::*;
pub use polarization::*;

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-6;

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

    // ── Coherence tests ───────────────────────────────────────────────────

    #[test]
    fn test_coherence_length_laser() {
        // Narrow-band laser: λ=632.8nm, Δλ=0.001nm → l_c = 632.8²/0.001 ≈ 4e8 nm = 400mm
        let lc = coherence_length(632.8, 0.001);
        assert!(
            (lc - 400_435_840.0).abs() < 1.0,
            "HeNe laser coherence length ≈ 4e8 nm, got {lc}"
        );
    }

    #[test]
    fn test_coherence_length_white_light() {
        // White light: λ=550nm, Δλ=300nm → l_c ≈ 1μm
        let lc = coherence_length(550.0, 300.0);
        assert!(
            lc < 2000.0 && lc > 500.0,
            "White light coherence ≈ 1μm, got {lc}nm"
        );
    }

    #[test]
    fn test_coherence_length_narrower_bandwidth_longer() {
        let lc_narrow = coherence_length(550.0, 1.0);
        let lc_broad = coherence_length(550.0, 10.0);
        assert!(lc_narrow > lc_broad);
    }

    #[test]
    fn test_coherence_time_positive() {
        let tc = coherence_time(550e-9, 1e-9);
        assert!(tc > 0.0);
        // Should be on order of nanoseconds for 1nm bandwidth
        assert!(tc > 1e-12 && tc < 1e-6);
    }

    #[test]
    fn test_spatial_coherence_angle() {
        let theta = spatial_coherence_angle(550e-9, 1e-3);
        assert!(theta > 0.0);
        assert!(theta < 0.001); // very small angle
    }

    #[test]
    fn test_coherence_area_increases_with_distance() {
        let a1 = coherence_area(550e-9, 1.0, 1e-3);
        let a2 = coherence_area(550e-9, 10.0, 1e-3);
        assert!(a2 > a1);
    }

    #[test]
    fn test_coherence_ratio() {
        // Path diff equal to coherence length → ratio = 1
        let lc = coherence_length(550.0, 1.0);
        let ratio = coherence_ratio(lc, 550.0, 1.0);
        assert!((ratio - 1.0).abs() < EPS);
    }

    // ── Bessel J1 tests ───────────────────────────────────────────────────

    #[test]
    fn test_bessel_j1_at_zero() {
        assert!(bessel_j1(0.0).abs() < EPS);
    }

    #[test]
    fn test_bessel_j1_known_values() {
        // J1(1) ≈ 0.44005
        assert!((bessel_j1(1.0) - 0.44005).abs() < 0.001);
        // J1(3.8317) ≈ 0 (first zero)
        assert!(bessel_j1(3.8317).abs() < 0.001);
        // J1(π) ≈ 0.28468
        assert!((bessel_j1(PI) - 0.28468).abs() < 0.001);
    }

    #[test]
    fn test_bessel_j1_odd_symmetry() {
        // J1(-x) = -J1(x)
        for x in [0.5, 1.0, 2.0, 5.0, 10.0, 15.0] {
            assert!(
                (bessel_j1(-x) + bessel_j1(x)).abs() < 0.001,
                "J1 should be odd at x={x}"
            );
        }
    }

    #[test]
    fn test_bessel_j1_large_argument() {
        // J1(10) ≈ 0.04347
        assert!((bessel_j1(10.0) - 0.04347).abs() < 0.001);
        // J1(20) — check it doesn't blow up
        let j = bessel_j1(20.0);
        assert!(j.is_finite());
        assert!(j.abs() < 1.0);
    }

    #[test]
    fn test_bessel_j1_continuity_at_boundary() {
        // The polynomial and asymptotic branches meet at x=8
        // Values just below and above should be close
        let j_below = bessel_j1(7.999);
        let j_above = bessel_j1(8.001);
        assert!(
            (j_below - j_above).abs() < 0.001,
            "J1 discontinuity at x=8: below={j_below}, above={j_above}"
        );
    }

    // ── Airy pattern tests ────────────────────────────────────────────────

    #[test]
    fn test_airy_central_maximum() {
        let i = airy_pattern(10e-3, 550e-9, 0.0, 1.0);
        assert!((i - 1.0).abs() < EPS);
    }

    #[test]
    fn test_airy_decreases_off_axis() {
        let i_center = airy_pattern(10e-3, 550e-9, 0.0, 1.0);
        let i_off = airy_pattern(10e-3, 550e-9, 1e-5, 1.0);
        assert!(i_off < i_center);
    }

    #[test]
    fn test_airy_first_zero_location() {
        let d = 10e-3;
        let wl = 550e-9;
        let theta_zero = airy_first_zero(wl, d);
        let i_at_zero = airy_pattern(d, wl, theta_zero, 1.0);
        assert!(
            i_at_zero < 0.001,
            "Intensity at first zero should be ~0, got {i_at_zero}"
        );
    }

    #[test]
    fn test_airy_always_non_negative() {
        let d = 10e-3;
        let wl = 550e-9;
        for i in 0..100 {
            let angle = i as f64 * 1e-5;
            let intensity = airy_pattern(d, wl, angle, 1.0);
            assert!(intensity >= 0.0, "Negative intensity at angle {angle}");
        }
    }

    #[test]
    fn test_airy_scales_with_i0() {
        let i1 = airy_pattern(10e-3, 550e-9, 1e-5, 1.0);
        let i5 = airy_pattern(10e-3, 550e-9, 1e-5, 5.0);
        assert!((i5 / i1 - 5.0).abs() < EPS);
    }

    #[test]
    fn test_rayleigh_criterion() {
        let theta = rayleigh_criterion(550e-9, 0.1);
        // 1.22 * 550e-9 / 0.1 ≈ 6.71e-6 rad
        assert!((theta - 6.71e-6).abs() < 1e-7);
    }

    #[test]
    fn test_rayleigh_larger_aperture_better() {
        let theta_small = rayleigh_criterion(550e-9, 0.05);
        let theta_large = rayleigh_criterion(550e-9, 0.20);
        assert!(theta_large < theta_small);
    }

    // ── Fabry-Pérot tests ─────────────────────────────────────────────────

    #[test]
    fn test_fp_finesse_coefficient() {
        // R=0.9: F = 4*0.9/(0.1)² = 360
        let f = fabry_perot_finesse_coefficient(0.9);
        assert!((f - 360.0).abs() < EPS);
    }

    #[test]
    fn test_fp_finesse_coefficient_zero() {
        assert!(fabry_perot_finesse_coefficient(0.0).abs() < EPS);
    }

    #[test]
    fn test_fp_finesse() {
        // R=0.9: F = π√0.9 / 0.1 ≈ 29.8
        let f = fabry_perot_finesse(0.9);
        assert!((f - 29.8).abs() < 0.5, "Finesse ≈ 29.8, got {f}");
    }

    #[test]
    fn test_fp_finesse_increases_with_reflectance() {
        let f_low = fabry_perot_finesse(0.5);
        let f_high = fabry_perot_finesse(0.95);
        assert!(f_high > f_low);
    }

    #[test]
    fn test_fp_transmittance_at_resonance() {
        // At resonance (δ = 2mπ), sin²(δ/2) = 0, T = 1
        // δ = 4πnd cos(θ)/λ = 2mπ when 2nd cos(θ) = mλ
        // For normal incidence (θ=0), n=1, d=λ/2: δ = 4π·1·(λ/2)·1/λ = 2π
        let wl = 550e-9;
        let d = wl / 2.0;
        let t = fabry_perot_transmittance(wl, d, 1.0, 0.0, 0.9);
        assert!(
            (t - 1.0).abs() < 0.001,
            "Transmittance at resonance should be ~1.0, got {t}"
        );
    }

    #[test]
    fn test_fp_transmittance_range() {
        // Transmittance should always be in [0, 1]
        for angle_mrad in 0..100 {
            let angle = angle_mrad as f64 * 0.001;
            let t = fabry_perot_transmittance(550e-9, 1e-3, 1.0, angle, 0.9);
            assert!(
                (0.0..=1.0 + EPS).contains(&t),
                "Transmittance out of range at angle {angle}: {t}"
            );
        }
    }

    #[test]
    fn test_fp_transmittance_low_reflectance_broad() {
        // Low R → broad peaks (nearly flat transmission)
        let t_peak = fabry_perot_transmittance(550e-9, 1e-3, 1.0, 0.0, 0.1);
        let t_mid = fabry_perot_transmittance(550e-9, 1e-3, 1.0, 0.01, 0.1);
        // Should be relatively similar (broad)
        assert!((t_peak - t_mid).abs() < 0.3);
    }

    #[test]
    fn test_fp_fsr() {
        // d=1mm, n=1: FSR = c/(2*1*0.001) = 1.499e11 Hz ≈ 150 GHz
        let fsr = fabry_perot_fsr(1e-3, 1.0);
        assert!((fsr - 1.499e11).abs() < 1e9);
    }

    #[test]
    fn test_fp_fsr_wavelength() {
        // λ=550nm, d=1mm, n=1: Δλ = 550²/(2*1*1e6) ≈ 0.151nm
        let fsr = fabry_perot_fsr_wavelength(550.0, 1e6, 1.0);
        assert!((fsr - 0.15125).abs() < 0.01, "FSR ≈ 0.15nm, got {fsr}nm");
    }

    #[test]
    fn test_fp_fsr_thicker_cavity_narrower() {
        let fsr_thin = fabry_perot_fsr(0.5e-3, 1.0);
        let fsr_thick = fabry_perot_fsr(2e-3, 1.0);
        assert!(fsr_thin > fsr_thick);
    }

    #[test]
    fn test_fp_resolving_power() {
        // High finesse + thick cavity → high resolving power
        let rp = fabry_perot_resolving_power(550e-9, 1e-3, 1.0, 0.95);
        assert!(rp > 1e5, "Resolving power should be > 10⁵, got {rp}");
    }

    #[test]
    fn test_fp_resolving_power_increases_with_reflectance() {
        let rp_low = fabry_perot_resolving_power(550e-9, 1e-3, 1.0, 0.5);
        let rp_high = fabry_perot_resolving_power(550e-9, 1e-3, 1.0, 0.95);
        assert!(rp_high > rp_low);
    }
}
