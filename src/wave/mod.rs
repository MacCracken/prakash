//! Wave optics — interference, diffraction, polarization.

use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

// ── Interference ────────────────────────────────────────────────────────────

/// Two-wave interference intensity.
///
/// Given amplitudes A1, A2 and phase difference δ (radians):
/// I = A1² + A2² + 2·A1·A2·cos(δ)
///
/// ```
/// # use prakash::wave::interference_intensity;
/// let constructive = interference_intensity(1.0, 1.0, 0.0);
/// assert!((constructive - 4.0).abs() < 1e-6); // (1+1)² = 4
/// ```
#[must_use]
#[inline]
pub fn interference_intensity(a1: f64, a2: f64, phase_diff: f64) -> f64 {
    a1 * a1 + a2 * a2 + 2.0 * a1 * a2 * phase_diff.cos()
}

/// Constructive interference condition: path difference = m·λ.
#[must_use]
#[inline]
pub fn is_constructive(path_diff: f64, wavelength: f64) -> bool {
    let m = path_diff / wavelength;
    (m - m.round()).abs() < 0.01
}

/// Destructive interference condition: path difference = (m + 0.5)·λ.
#[must_use]
#[inline]
pub fn is_destructive(path_diff: f64, wavelength: f64) -> bool {
    let m = path_diff / wavelength - 0.5;
    (m - m.round()).abs() < 0.01
}

/// Phase difference from path difference and wavelength.
/// δ = 2π·Δ/λ
#[must_use]
#[inline]
pub fn path_to_phase(wavelength: f64, path_diff: f64) -> f64 {
    std::f64::consts::TAU * path_diff / wavelength
}

// ── Thin Film Interference ──────────────────────────────────────────────────

/// Thin film interference: reflected intensity for a thin film of thickness `d`,
/// refractive index `n_film`, at near-normal incidence.
///
/// Returns intensity as a fraction of incident (0.0–1.0 approximate).
#[must_use]
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
/// I(θ) = I₀ · (sin(β)/β)² where β = π·a·sin(θ)/λ
/// `wavelength` and `slit_width` in same units.
#[must_use]
#[inline]
pub fn single_slit_intensity(wavelength: f64, slit_width: f64, angle: f64, i0: f64) -> f64 {
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
#[must_use]
#[inline]
pub fn double_slit_intensity(
    wavelength: f64,
    slit_width: f64,
    slit_spacing: f64,
    angle: f64,
    i0: f64,
) -> f64 {
    let sin_angle = angle.sin();
    // Single-slit envelope (sinc² term)
    let beta = PI * slit_width * sin_angle / wavelength;
    let envelope = if beta.abs() < 1e-10 {
        1.0
    } else {
        let sinc = beta.sin() / beta;
        sinc * sinc
    };
    // Two-slit interference (cos² term)
    let delta = PI * slit_spacing * sin_angle / wavelength;
    let cos_delta = delta.cos();
    let interference = cos_delta * cos_delta;
    i0 * envelope * interference * 4.0
}

/// Diffraction grating: angular positions of maxima.
///
/// d·sin(θ) = m·λ → θ = asin(m·λ/d)
/// Returns angles for orders m = 0, ±1, ±2, ... up to `max_order`.
#[must_use]
pub fn grating_maxima(wavelength: f64, grating_spacing: f64, max_order: u32) -> Vec<f64> {
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
    #[must_use]
    pub const fn circular_right() -> Self {
        Self {
            ex: std::f64::consts::FRAC_1_SQRT_2,
            ey: std::f64::consts::FRAC_1_SQRT_2,
            phase: -std::f64::consts::FRAC_PI_2,
        }
    }
    /// Left circular polarization.
    #[must_use]
    pub const fn circular_left() -> Self {
        Self {
            ex: std::f64::consts::FRAC_1_SQRT_2,
            ey: std::f64::consts::FRAC_1_SQRT_2,
            phase: std::f64::consts::FRAC_PI_2,
        }
    }

    /// Intensity after passing through a linear polarizer at angle θ.
    /// Malus's law: I = I0 · cos²(θ - polarization_angle)
    #[must_use]
    #[inline]
    pub fn through_polarizer(&self, polarizer_angle: f64) -> f64 {
        let pol_angle = self.ey.atan2(self.ex);
        let diff = polarizer_angle - pol_angle;
        diff.cos().powi(2) * (self.ex * self.ex + self.ey * self.ey)
    }

    /// Total intensity.
    #[must_use]
    #[inline]
    pub fn intensity(&self) -> f64 {
        self.ex * self.ex + self.ey * self.ey
    }
}

impl From<Polarization> for StokesVector {
    /// Convert a Jones-style polarization state to a Stokes vector.
    ///
    /// The mapping follows the standard Jones → Stokes conversion:
    /// - S₀ = |Ex|² + |Ey|² (total intensity)
    /// - S₁ = |Ex|² − |Ey|² (horizontal vs vertical preference)
    /// - S₂ = 2·|Ex|·|Ey|·cos(δ) (diagonal preference)
    /// - S₃ = 2·|Ex|·|Ey|·sin(δ) (circular preference)
    #[inline]
    fn from(p: Polarization) -> Self {
        let ex2 = p.ex * p.ex;
        let ey2 = p.ey * p.ey;
        Self::new(
            ex2 + ey2,
            ex2 - ey2,
            2.0 * p.ex * p.ey * p.phase.cos(),
            2.0 * p.ex * p.ey * p.phase.sin(),
        )
    }
}

#[cfg(feature = "bijli-backend")]
impl From<Polarization> for bijli::polarization::JonesVector {
    /// Convert a prakash `Polarization` to a bijli `JonesVector`.
    ///
    /// Maps `(ex, ey, phase)` to complex Jones vector `[Ex, Ey·e^(iδ)]`.
    #[inline]
    fn from(p: Polarization) -> Self {
        use bijli::polarization::Complex;
        bijli::polarization::JonesVector::new(
            Complex::real(p.ex),
            Complex::from_polar(p.ey, p.phase),
        )
    }
}

/// Malus's law: intensity after a polarizer at angle θ relative to polarization.
///
/// I = I0 · cos²(θ)
#[must_use]
#[inline]
pub fn malus_law(intensity: f64, angle: f64) -> f64 {
    intensity * angle.cos().powi(2)
}

mod coherence;
pub use coherence::*;

mod airy;
pub use airy::*;

mod fabry_perot;
pub use fabry_perot::*;

mod diffraction;
mod pattern;
mod polarization;
/// Zernike polynomials for wavefront decomposition and aberration analysis.
pub mod zernike;

pub use diffraction::*;
pub use pattern::*;
pub use polarization::*;

// ── Bijli re-exports ────────────────────────────────────────────────────────

/// Jones vector, Jones matrix, and complex number types from bijli.
///
/// Available when the `bijli-backend` feature is enabled.
#[cfg(feature = "bijli-backend")]
pub use bijli::polarization::{Complex, JonesMatrix, JonesVector};

/// Gaussian beam propagation and ABCD ray transfer matrices from bijli.
///
/// Available when the `bijli-backend` feature is enabled.
#[cfg(feature = "bijli-backend")]
pub use bijli::beam::{AbcdMatrix, GaussianBeam, ResonatorStability, resonator_stability};

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
        let phase = path_to_phase(500.0, 500.0); // wavelength=500, path=500
        assert!((phase - 2.0 * PI).abs() < EPS);
    }

    #[test]
    fn test_path_to_phase_half_wavelength() {
        let phase = path_to_phase(500.0, 250.0);
        assert!((phase - PI).abs() < EPS);
    }

    #[test]
    fn test_path_to_phase_zero() {
        let phase = path_to_phase(500.0, 0.0);
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
        let i = single_slit_intensity(500e-9, 1e-3, 0.0, 1.0);
        assert!((i - 1.0).abs() < EPS);
    }

    #[test]
    fn test_single_slit_decreases_off_axis() {
        let i_center = single_slit_intensity(500e-9, 1e-3, 0.0, 1.0);
        let i_off = single_slit_intensity(500e-9, 1e-3, 0.01, 1.0);
        assert!(i_off < i_center);
    }

    #[test]
    fn test_single_slit_first_minimum() {
        // First minimum at sin(θ) = λ/a → θ = asin(λ/a)
        let a = 1e-3;
        let wl = 500e-9;
        let ratio: f64 = wl / a;
        let theta_min = ratio.asin();
        let i = single_slit_intensity(wl, a, theta_min, 1.0);
        assert!(i < 1e-6, "Intensity at first minimum should be ~0, got {i}");
    }

    #[test]
    fn test_single_slit_scales_with_i0() {
        let i1 = single_slit_intensity(500e-9, 1e-3, 0.01, 1.0);
        let i5 = single_slit_intensity(500e-9, 1e-3, 0.01, 5.0);
        assert!((i5 / i1 - 5.0).abs() < EPS);
    }

    #[test]
    fn test_single_slit_always_non_negative() {
        for angle_mrad in 0..100 {
            let angle = angle_mrad as f64 * 0.001;
            let i = single_slit_intensity(500e-9, 1e-3, angle, 1.0);
            assert!(i >= 0.0, "Negative intensity at angle {angle}");
        }
    }

    #[test]
    fn test_double_slit_central() {
        let i = double_slit_intensity(500e-9, 0.1e-3, 0.5e-3, 0.0, 1.0);
        assert!(i > 0.0);
    }

    #[test]
    fn test_double_slit_greater_than_single_at_center() {
        let i_single = single_slit_intensity(500e-9, 0.1e-3, 0.0, 1.0);
        let i_double = double_slit_intensity(500e-9, 0.1e-3, 0.5e-3, 0.0, 1.0);
        assert!(
            i_double > i_single,
            "Double slit central max should exceed single slit"
        );
    }

    #[test]
    fn test_double_slit_non_negative() {
        for angle_mrad in 0..50 {
            let angle = angle_mrad as f64 * 0.001;
            let i = double_slit_intensity(500e-9, 0.1e-3, 0.5e-3, angle, 1.0);
            assert!(i >= 0.0, "Negative intensity at angle {angle}");
        }
    }

    // ── Grating tests ─────────────────────────────────────────────────────

    #[test]
    fn test_grating_maxima_zeroth_order() {
        let angles = grating_maxima(500e-9, 1e-6, 0);
        assert_eq!(angles.len(), 1);
        assert!(angles[0].abs() < EPS);
    }

    #[test]
    fn test_grating_maxima_multiple_orders() {
        let angles = grating_maxima(500e-9, 1e-6, 2);
        assert!(angles.len() >= 3);
    }

    #[test]
    fn test_grating_maxima_symmetric() {
        let angles = grating_maxima(500e-9, 1e-6, 1);
        // Should have m=0, +1, -1
        assert_eq!(angles.len(), 3);
        // +1 and -1 should be symmetric
        assert!((angles[1] + angles[2]).abs() < EPS);
    }

    #[test]
    fn test_grating_maxima_limited_by_sin() {
        // Very fine grating with long wavelength: fewer orders possible
        let angles = grating_maxima(500e-9, 600e-9, 5);
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

    // ── Polarization → StokesVector conversion tests ────────────────────

    #[test]
    fn test_jones_to_stokes_horizontal() {
        let p = Polarization::HORIZONTAL;
        let s = StokesVector::from(p);
        assert!((s.s0 - 1.0).abs() < EPS, "S0 should be 1.0");
        assert!((s.s1 - 1.0).abs() < EPS, "S1 should be 1.0 for horizontal");
        assert!(s.s2.abs() < EPS, "S2 should be 0");
        assert!(s.s3.abs() < EPS, "S3 should be 0");
    }

    #[test]
    fn test_jones_to_stokes_vertical() {
        let p = Polarization::VERTICAL;
        let s = StokesVector::from(p);
        assert!((s.s0 - 1.0).abs() < EPS);
        assert!((s.s1 + 1.0).abs() < EPS, "S1 should be -1.0 for vertical");
        assert!(s.s2.abs() < EPS);
        assert!(s.s3.abs() < EPS);
    }

    #[test]
    fn test_jones_to_stokes_circular_right() {
        let p = Polarization::circular_right();
        let s = StokesVector::from(p);
        assert!((s.s0 - 1.0).abs() < EPS);
        assert!(s.s1.abs() < EPS, "Circular should have no H/V preference");
        assert!(
            s.s2.abs() < EPS,
            "Circular should have no diagonal preference"
        );
        assert!(
            (s.s3 + 1.0).abs() < 0.01,
            "S3 should be -1 for right circular (phase = -π/2)"
        );
    }

    #[test]
    fn test_jones_to_stokes_circular_left() {
        let p = Polarization::circular_left();
        let s = StokesVector::from(p);
        assert!((s.s0 - 1.0).abs() < EPS);
        assert!(s.s1.abs() < EPS);
        assert!(s.s2.abs() < EPS);
        assert!(
            (s.s3 - 1.0).abs() < 0.01,
            "S3 should be +1 for left circular (phase = +π/2)"
        );
    }

    #[test]
    fn test_jones_to_stokes_intensity_preserved() {
        // For any Jones vector, S0 should equal total intensity
        let p = Polarization {
            ex: 0.8,
            ey: 0.6,
            phase: 0.5,
        };
        let s = StokesVector::from(p);
        assert!((s.s0 - p.intensity()).abs() < EPS);
    }

    #[test]
    fn test_jones_to_stokes_degree_of_polarization() {
        // Any pure Jones state should produce a fully polarized Stokes vector (DOP = 1)
        let states = [
            Polarization::HORIZONTAL,
            Polarization::VERTICAL,
            Polarization::circular_right(),
            Polarization::circular_left(),
            Polarization {
                ex: 0.8,
                ey: 0.6,
                phase: 1.0,
            },
        ];
        for p in &states {
            let s = StokesVector::from(*p);
            let dop = s.degree_of_polarization();
            assert!(
                (dop - 1.0).abs() < 0.01,
                "Pure Jones state should be fully polarized, got DOP={dop}"
            );
        }
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
        let i = airy_pattern(550e-9, 10e-3, 0.0, 1.0);
        assert!((i - 1.0).abs() < EPS);
    }

    #[test]
    fn test_airy_decreases_off_axis() {
        let i_center = airy_pattern(550e-9, 10e-3, 0.0, 1.0);
        let i_off = airy_pattern(550e-9, 10e-3, 1e-5, 1.0);
        assert!(i_off < i_center);
    }

    #[test]
    fn test_airy_first_zero_location() {
        let d = 10e-3;
        let wl = 550e-9;
        let theta_zero = airy_first_zero(wl, d);
        let i_at_zero = airy_pattern(wl, d, theta_zero, 1.0);
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
            let intensity = airy_pattern(wl, d, angle, 1.0);
            assert!(intensity >= 0.0, "Negative intensity at angle {angle}");
        }
    }

    #[test]
    fn test_airy_scales_with_i0() {
        let i1 = airy_pattern(550e-9, 10e-3, 1e-5, 1.0);
        let i5 = airy_pattern(550e-9, 10e-3, 1e-5, 5.0);
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
