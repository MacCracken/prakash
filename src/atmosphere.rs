//! Atmospheric optics — Rayleigh/Mie scattering, sky color, air mass, optical depth.
//!
//! Models how light interacts with Earth's atmosphere: wavelength-dependent
//! scattering (blue sky, red sunsets), aerosol extinction, and the optical
//! path through curved atmospheric layers.

use std::f64::consts::PI;
use tracing::trace;

// ── Physical Constants ──────────────────────────────────────────────────────

/// Index of refraction of air at sea level (STP, ~550nm).
const N_AIR: f64 = 1.000_293;

/// Molecular number density at sea level (molecules/m³), standard atmosphere.
const N_S: f64 = 2.504e25;

/// Scale height of the atmosphere for Rayleigh scattering (meters).
/// Characteristic height at which density drops by 1/e.
pub const SCALE_HEIGHT_RAYLEIGH: f64 = 8500.0;

/// Scale height for Mie (aerosol) scattering (meters).
pub const SCALE_HEIGHT_MIE: f64 = 1200.0;

/// Typical sea-level Mie scattering coefficient (1/m) for clear atmosphere.
const BETA_M_SEA_LEVEL: f64 = 21.0e-6;

/// Representative wavelengths for RGB channels (meters): R=680nm, G=550nm, B=440nm.
const RGB_WAVELENGTHS: [f64; 3] = [680e-9, 550e-9, 440e-9];

/// Precomputed numerator for Rayleigh cross-section: (8π³/3) · (n²−1)² / N².
const RAYLEIGH_PREFACTOR: f64 = {
    // (8π³/3) · (N_AIR²−1)² / N_S²
    // Computed at compile time via const evaluation
    let n2m1_sq = (N_AIR * N_AIR - 1.0) * (N_AIR * N_AIR - 1.0);
    let pi3 = PI * PI * PI;
    (8.0 * pi3 / 3.0) * n2m1_sq / (N_S * N_S)
};

// ── Rayleigh Scattering ─────────────────────────────────────────────────────

/// Rayleigh scattering cross-section for a single molecule.
///
/// σ(λ) = (8π³/3) · (n²−1)² / (N²·λ⁴)
///
/// where n = refractive index of air, N = number density at STP.
///
/// `wavelength_m` in meters. Returns cross-section in m².
#[must_use]
#[inline]
pub fn rayleigh_cross_section(wavelength_m: f64) -> f64 {
    let lambda4 = wavelength_m * wavelength_m * wavelength_m * wavelength_m;
    RAYLEIGH_PREFACTOR / lambda4
}

/// Rayleigh scattering coefficient at sea level.
///
/// β(λ) = N · σ(λ) = (8π³/3) · (n²−1)² / (N·λ⁴)
///
/// `wavelength_m` in meters. Returns scattering coefficient in 1/m.
#[must_use]
#[inline]
pub fn rayleigh_scattering_coefficient(wavelength_m: f64) -> f64 {
    N_S * rayleigh_cross_section(wavelength_m)
}

/// Rayleigh scattering coefficient at a given altitude.
///
/// β(λ,h) = β₀(λ) · exp(-h / H_R)
///
/// `wavelength_m` in meters, `altitude_m` in meters above sea level.
/// Returns scattering coefficient in 1/m.
#[must_use]
#[inline]
pub fn rayleigh_scattering_at_altitude(wavelength_m: f64, altitude_m: f64) -> f64 {
    rayleigh_scattering_coefficient(wavelength_m) * (-altitude_m / SCALE_HEIGHT_RAYLEIGH).exp()
}

/// Rayleigh phase function.
///
/// p(θ) = 3(1 + cos²θ) / (16π)
///
/// `cos_theta` = cosine of scattering angle.
#[must_use]
#[inline]
pub fn rayleigh_phase(cos_theta: f64) -> f64 {
    3.0 * (1.0 + cos_theta * cos_theta) / (16.0 * PI)
}

// ── Mie Scattering ──────────────────────────────────────────────────────────

/// Mie scattering coefficient at sea level for typical aerosol conditions.
///
/// β_M(λ) = β_M_ref · (λ_ref / λ)^α
///
/// where α ≈ 0 for large particles (Mie regime). We use α = 0 (wavelength-independent)
/// which is the standard simplification for atmospheric aerosols.
///
/// `wavelength_m` is unused in this simplified model but kept for API consistency.
/// Returns scattering coefficient in 1/m.
#[must_use]
#[inline]
pub fn mie_scattering_coefficient(_wavelength_m: f64) -> f64 {
    BETA_M_SEA_LEVEL
}

/// Mie scattering coefficient at a given altitude.
///
/// β_M(h) = β_M₀ · exp(-h / H_M)
///
/// `altitude_m` in meters above sea level.
#[must_use]
#[inline]
pub fn mie_scattering_at_altitude(wavelength_m: f64, altitude_m: f64) -> f64 {
    mie_scattering_coefficient(wavelength_m) * (-altitude_m / SCALE_HEIGHT_MIE).exp()
}

/// Cornette-Shanks phase function for Mie scattering.
///
/// Improved Henyey-Greenstein that better models atmospheric aerosols:
///
/// p(θ) = (3/8π) · (1 − g²)(1 + cos²θ) / ((2 + g²)(1 + g² − 2g·cosθ)^(3/2))
///
/// `g` = asymmetry parameter (typically 0.76 for atmospheric aerosols).
#[must_use]
#[inline]
pub fn mie_phase_cornette_shanks(cos_theta: f64, g: f64) -> f64 {
    let g2 = g * g;
    let denom_base = 1.0 + g2 - 2.0 * g * cos_theta;
    let denom = (2.0 + g2) * denom_base * denom_base.sqrt();
    if denom < 1e-15 {
        return 0.0;
    }
    (3.0 / (8.0 * PI)) * (1.0 - g2) * (1.0 + cos_theta * cos_theta) / denom
}

/// Default Mie asymmetry parameter for Earth's atmosphere.
pub const MIE_G_DEFAULT: f64 = 0.76;

// ── Air Mass & Optical Depth ────────────────────────────────────────────────

/// Relative air mass using the Kasten & Young (1989) formula.
///
/// Accounts for atmospheric curvature. More accurate than sec(θ) at
/// high zenith angles (horizon).
///
/// `zenith_angle` in radians (0 = overhead, π/2 = horizon).
/// Returns the relative air mass (1.0 at zenith).
#[must_use]
#[inline]
pub fn air_mass(zenith_angle: f64) -> f64 {
    let zenith_deg = zenith_angle.to_degrees().min(90.0);
    // Kasten & Young (1989): m = 1 / (cos(θ) + 0.50572·(96.07995 − θ)^(-1.6364))
    let cos_z = zenith_angle.cos();
    let correction = 0.505_72 * (96.079_95 - zenith_deg).powf(-1.6364);
    1.0 / (cos_z + correction)
}

/// Rayleigh optical depth at sea level for a given wavelength.
///
/// τ_R(λ) = β_R(λ) · H_R
///
/// This is the total optical depth looking straight up (zenith) through
/// the entire atmosphere.
///
/// `wavelength_m` in meters.
#[must_use]
#[inline]
pub fn optical_depth_rayleigh(wavelength_m: f64) -> f64 {
    rayleigh_scattering_coefficient(wavelength_m) * SCALE_HEIGHT_RAYLEIGH
}

/// Mie optical depth at sea level.
///
/// τ_M = β_M · H_M
#[must_use]
#[inline]
pub fn optical_depth_mie() -> f64 {
    BETA_M_SEA_LEVEL * SCALE_HEIGHT_MIE
}

/// Total transmittance along a path through the atmosphere.
///
/// T(λ,θ) = exp(−(τ_R(λ) + τ_M) · m(θ))
///
/// `wavelength_m` in meters, `zenith_angle` in radians.
#[must_use]
#[inline]
pub fn atmospheric_transmittance(wavelength_m: f64, zenith_angle: f64) -> f64 {
    let m = air_mass(zenith_angle);
    let tau_r = optical_depth_rayleigh(wavelength_m);
    let tau_m = optical_depth_mie();
    (-(tau_r + tau_m) * m).exp()
}

// ── Sky Color ───────────────────────────────────────────────────────────────

/// Single-scattering sky radiance at a point in the sky.
///
/// Computes the spectral radiance (relative) at a given wavelength from a
/// single-scattering model. Integrates scattering along the view ray.
///
/// `wavelength_m` in meters. `sun_zenith` and `view_zenith` in radians.
/// `scattering_angle` = angle between sun direction and view direction (radians).
///
/// Returns relative radiance (arbitrary units, useful for ratios between wavelengths).
#[must_use]
#[inline]
pub fn sky_radiance_single_scatter(
    wavelength_m: f64,
    sun_zenith: f64,
    scattering_angle: f64,
) -> f64 {
    let cos_scatter = scattering_angle.cos();
    let beta_r = rayleigh_scattering_coefficient(wavelength_m);
    let beta_m = mie_scattering_coefficient(wavelength_m);

    // Rayleigh contribution
    let phase_r = rayleigh_phase(cos_scatter);
    let radiance_r = beta_r * phase_r;

    // Mie contribution
    let phase_m = mie_phase_cornette_shanks(cos_scatter, MIE_G_DEFAULT);
    let radiance_m = beta_m * phase_m;

    // Attenuate by transmittance along sun path
    let sun_transmittance = atmospheric_transmittance(wavelength_m, sun_zenith);

    (radiance_r + radiance_m) * sun_transmittance
}

/// Sky color as RGB for a given sun position and viewing direction.
///
/// Uses single-scattering Rayleigh + Mie model at three representative
/// wavelengths (R=680nm, G=550nm, B=440nm).
///
/// `sun_zenith` in radians (0 = overhead, π/2 = horizon).
/// `scattering_angle` in radians (angle between sun and view direction).
///
/// Returns `[r, g, b]` in relative units (normalize for display).
#[must_use]
pub fn sky_color_rgb(sun_zenith: f64, scattering_angle: f64) -> [f64; 3] {
    trace!(sun_zenith, scattering_angle, "sky_color_rgb");
    let mut rgb = [0.0; 3];
    for (i, &wl) in RGB_WAVELENGTHS.iter().enumerate() {
        rgb[i] = sky_radiance_single_scatter(wl, sun_zenith, scattering_angle);
    }
    rgb
}

// ── Sunset / Sunrise ────────────────────────────────────────────────────────

/// Color of direct sunlight after passing through the atmosphere.
///
/// At low sun angles, the long atmospheric path preferentially removes
/// short wavelengths (blue), leaving warm reds/oranges.
///
/// `sun_zenith` in radians (π/2 = horizon = deep sunset).
///
/// Returns `[r, g, b]` transmittance (multiply by sun's intrinsic white).
#[must_use]
pub fn sunlight_color(sun_zenith: f64) -> [f64; 3] {
    trace!(sun_zenith, "sunlight_color");
    let mut rgb = [0.0; 3];
    for (i, &wl) in RGB_WAVELENGTHS.iter().enumerate() {
        rgb[i] = atmospheric_transmittance(wl, sun_zenith);
    }
    rgb
}

/// Sunset/sunrise sky gradient model.
///
/// Returns the sky color `[r, g, b]` at a given elevation above the horizon,
/// with the sun at a given zenith angle. Combines direct sunlight scattering
/// near the sun with Rayleigh blue away from it.
///
/// `sun_zenith` in radians (values near π/2 give sunset conditions).
/// `view_elevation` in radians above the horizon (0 = horizon, π/2 = zenith).
/// `angular_distance_to_sun` in radians.
///
/// Returns `[r, g, b]` in relative units.
#[must_use]
pub fn sunset_gradient(
    sun_zenith: f64,
    view_elevation: f64,
    angular_distance_to_sun: f64,
) -> [f64; 3] {
    trace!(
        sun_zenith,
        view_elevation, angular_distance_to_sun, "sunset_gradient"
    );
    let view_zenith = (PI / 2.0 - view_elevation).clamp(0.0, PI / 2.0);

    // Lift angle-dependent computations out of per-wavelength loop
    let cos_scatter = angular_distance_to_sun.cos();
    let phase_r = rayleigh_phase(cos_scatter);
    let phase_m = mie_phase_cornette_shanks(cos_scatter, MIE_G_DEFAULT);
    let sun_air_mass = air_mass(sun_zenith);
    let view_air_mass = air_mass(view_zenith);
    let tau_m = optical_depth_mie();

    let mut rgb = [0.0; 3];
    for (i, &wl) in RGB_WAVELENGTHS.iter().enumerate() {
        let beta_r = rayleigh_scattering_coefficient(wl);
        let tau_r = beta_r * SCALE_HEIGHT_RAYLEIGH;
        let total_tau = tau_r + tau_m;

        // Sun and view path transmittance
        let sun_transmittance = (-total_tau * sun_air_mass).exp();
        let view_transmittance = (-total_tau * view_air_mass).exp();

        // Scattering along view ray (single-scatter approximation)
        let scatter = (beta_r * phase_r + BETA_M_SEA_LEVEL * phase_m) * sun_transmittance;
        rgb[i] = scatter * view_transmittance;
    }
    rgb
}

/// Compute the scattering angle between sun and view directions.
///
/// Given sun and view directions as zenith + azimuth angles:
/// cos(γ) = cos(θ_s)·cos(θ_v) + sin(θ_s)·sin(θ_v)·cos(φ_v − φ_s)
///
/// `sun_zenith`, `view_zenith` in radians. `azimuth_diff` = |φ_view − φ_sun| in radians.
#[must_use]
#[inline]
pub fn scattering_angle(sun_zenith: f64, view_zenith: f64, azimuth_diff: f64) -> f64 {
    let cos_gamma = sun_zenith.cos() * view_zenith.cos()
        + sun_zenith.sin() * view_zenith.sin() * azimuth_diff.cos();
    cos_gamma.clamp(-1.0, 1.0).acos()
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-6;

    // ── Rayleigh scattering ──────────────────────────────────────────────

    #[test]
    fn test_rayleigh_cross_section_positive() {
        let sigma = rayleigh_cross_section(550e-9);
        assert!(sigma > 0.0, "Cross-section must be positive");
    }

    #[test]
    fn test_rayleigh_cross_section_order_of_magnitude() {
        // At 550nm, σ ≈ 4.5e-31 m² (order of magnitude check)
        let sigma = rayleigh_cross_section(550e-9);
        assert!(
            sigma > 1e-32 && sigma < 1e-29,
            "Cross-section should be ~1e-31 m², got {sigma}"
        );
    }

    #[test]
    fn test_rayleigh_lambda4_dependence() {
        // Blue light should scatter much more than red (λ⁻⁴)
        let sigma_blue = rayleigh_cross_section(440e-9);
        let sigma_red = rayleigh_cross_section(680e-9);
        let ratio = sigma_blue / sigma_red;
        let expected = (680.0 / 440.0_f64).powi(4);
        assert!(
            (ratio - expected).abs() / expected < 0.01,
            "Should follow λ⁻⁴: ratio={ratio}, expected={expected}"
        );
    }

    #[test]
    fn test_rayleigh_scattering_coefficient_positive() {
        let beta = rayleigh_scattering_coefficient(550e-9);
        assert!(beta > 0.0);
        // At 550nm, β ≈ 1.1e-5 1/m
        assert!(
            beta > 1e-6 && beta < 1e-4,
            "β should be ~1e-5 1/m, got {beta}"
        );
    }

    #[test]
    fn test_rayleigh_blue_scatters_more_than_red() {
        let beta_blue = rayleigh_scattering_coefficient(440e-9);
        let beta_red = rayleigh_scattering_coefficient(680e-9);
        assert!(
            beta_blue > beta_red * 4.0,
            "Blue should scatter much more than red"
        );
    }

    #[test]
    fn test_rayleigh_at_altitude_decreases() {
        let beta_sea = rayleigh_scattering_at_altitude(550e-9, 0.0);
        let beta_high = rayleigh_scattering_at_altitude(550e-9, 8500.0);
        assert!(
            (beta_high / beta_sea - 1.0 / std::f64::consts::E).abs() < 0.01,
            "At scale height, density should be 1/e of sea level"
        );
    }

    #[test]
    fn test_rayleigh_at_altitude_zero_matches_sea_level() {
        let beta_0 = rayleigh_scattering_at_altitude(550e-9, 0.0);
        let beta_sea = rayleigh_scattering_coefficient(550e-9);
        assert!((beta_0 - beta_sea).abs() < EPS);
    }

    #[test]
    fn test_rayleigh_phase_symmetric() {
        let p_fwd = rayleigh_phase(1.0);
        let p_back = rayleigh_phase(-1.0);
        assert!((p_fwd - p_back).abs() < EPS);
    }

    #[test]
    fn test_rayleigh_phase_min_at_90() {
        let p_90 = rayleigh_phase(0.0);
        let p_0 = rayleigh_phase(1.0);
        assert!(p_0 > p_90);
    }

    #[test]
    fn test_rayleigh_phase_normalizes() {
        // Phase function should integrate to 1 over 4π steradians
        // Numerical check: ∫ p(θ) sin(θ) dθ dφ = 1 (over θ=0..π, φ=0..2π)
        let n = 1000;
        let mut integral = 0.0;
        for i in 0..n {
            let theta = PI * (i as f64 + 0.5) / n as f64;
            let cos_t = theta.cos();
            integral += rayleigh_phase(cos_t) * theta.sin() * 2.0 * PI * (PI / n as f64);
        }
        assert!(
            (integral - 1.0).abs() < 0.01,
            "Phase function integral ≈ 1, got {integral}"
        );
    }

    // ── Mie scattering ───────────────────────────────────────────────────

    #[test]
    fn test_mie_coefficient_positive() {
        assert!(mie_scattering_coefficient(550e-9) > 0.0);
    }

    #[test]
    fn test_mie_at_altitude_decreases() {
        let beta_sea = mie_scattering_at_altitude(550e-9, 0.0);
        let beta_high = mie_scattering_at_altitude(550e-9, 1200.0);
        assert!(
            (beta_high / beta_sea - 1.0 / std::f64::consts::E).abs() < 0.01,
            "At Mie scale height, should be 1/e of sea level"
        );
    }

    #[test]
    fn test_mie_phase_forward_peak() {
        let p_fwd = mie_phase_cornette_shanks(1.0, 0.76);
        let p_back = mie_phase_cornette_shanks(-1.0, 0.76);
        assert!(p_fwd > p_back * 10.0, "Mie should have strong forward peak");
    }

    #[test]
    fn test_mie_phase_positive() {
        for g in [0.0, 0.5, 0.76, 0.9] {
            for ct in [-1.0, -0.5, 0.0, 0.5, 1.0] {
                let p = mie_phase_cornette_shanks(ct, g);
                assert!(p >= 0.0, "Mie phase negative at g={g}, cos_theta={ct}");
            }
        }
    }

    #[test]
    fn test_mie_phase_isotropic_at_g0() {
        // g=0 should be nearly isotropic (but with cos² modulation from Cornette-Shanks)
        let p_fwd = mie_phase_cornette_shanks(1.0, 0.0);
        let p_90 = mie_phase_cornette_shanks(0.0, 0.0);
        // Cornette-Shanks at g=0: p = (3/16π)(1+cos²θ), not truly isotropic
        assert!(p_fwd > p_90, "Forward should still be > 90° even at g=0");
    }

    #[test]
    fn test_mie_phase_normalizes() {
        let n = 1000;
        let mut integral = 0.0;
        for i in 0..n {
            let theta = PI * (i as f64 + 0.5) / n as f64;
            let cos_t = theta.cos();
            integral +=
                mie_phase_cornette_shanks(cos_t, 0.76) * theta.sin() * 2.0 * PI * (PI / n as f64);
        }
        assert!(
            (integral - 1.0).abs() < 0.02,
            "Mie phase integral ≈ 1, got {integral}"
        );
    }

    // ── Air mass ─────────────────────────────────────────────────────────

    #[test]
    fn test_air_mass_zenith() {
        let m = air_mass(0.0);
        assert!((m - 1.0).abs() < 0.01, "Air mass at zenith should be ~1.0");
    }

    #[test]
    fn test_air_mass_60_degrees() {
        let m = air_mass(60.0_f64.to_radians());
        assert!(
            (m - 2.0).abs() < 0.05,
            "Air mass at 60° ≈ 2.0 (sec 60°), got {m}"
        );
    }

    #[test]
    fn test_air_mass_horizon_finite() {
        // At horizon, simple sec(90°) = ∞, but Kasten formula gives ~38
        let m = air_mass(PI / 2.0);
        assert!(m.is_finite(), "Air mass at horizon should be finite");
        assert!(m > 30.0 && m < 45.0, "Air mass at horizon ≈ 38, got {m}");
    }

    #[test]
    fn test_air_mass_increases_toward_horizon() {
        let m_30 = air_mass(30.0_f64.to_radians());
        let m_60 = air_mass(60.0_f64.to_radians());
        let m_80 = air_mass(80.0_f64.to_radians());
        assert!(m_30 < m_60);
        assert!(m_60 < m_80);
    }

    // ── Optical depth ────────────────────────────────────────────────────

    #[test]
    fn test_optical_depth_rayleigh_positive() {
        let tau = optical_depth_rayleigh(550e-9);
        assert!(tau > 0.0);
        // At 550nm, τ_R ≈ 0.09
        assert!(
            tau > 0.05 && tau < 0.2,
            "τ_R at 550nm should be ~0.1, got {tau}"
        );
    }

    #[test]
    fn test_optical_depth_blue_greater_than_red() {
        let tau_blue = optical_depth_rayleigh(440e-9);
        let tau_red = optical_depth_rayleigh(680e-9);
        assert!(tau_blue > tau_red);
    }

    #[test]
    fn test_optical_depth_mie_positive() {
        let tau = optical_depth_mie();
        assert!(tau > 0.0);
    }

    // ── Transmittance ────────────────────────────────────────────────────

    #[test]
    fn test_transmittance_range() {
        for angle_deg in [0.0_f64, 30.0, 60.0, 80.0] {
            let t = atmospheric_transmittance(550e-9, angle_deg.to_radians());
            assert!(
                (0.0..=1.0).contains(&t),
                "Transmittance out of range at {angle_deg}°: {t}"
            );
        }
    }

    #[test]
    fn test_transmittance_zenith_high() {
        let t = atmospheric_transmittance(550e-9, 0.0);
        assert!(
            t > 0.8,
            "Zenith transmittance at 550nm should be high, got {t}"
        );
    }

    #[test]
    fn test_transmittance_decreases_toward_horizon() {
        let t_0 = atmospheric_transmittance(550e-9, 0.0);
        let t_60 = atmospheric_transmittance(550e-9, 60.0_f64.to_radians());
        let t_85 = atmospheric_transmittance(550e-9, 85.0_f64.to_radians());
        assert!(t_0 > t_60);
        assert!(t_60 > t_85);
    }

    #[test]
    fn test_transmittance_red_higher_than_blue() {
        // Red light gets through better (less Rayleigh scattering)
        let t_red = atmospheric_transmittance(680e-9, 80.0_f64.to_radians());
        let t_blue = atmospheric_transmittance(440e-9, 80.0_f64.to_radians());
        assert!(
            t_red > t_blue,
            "Red should transmit better than blue at low angles"
        );
    }

    // ── Sky color ────────────────────────────────────────────────────────

    #[test]
    fn test_sky_color_overhead_is_blue() {
        // Sun overhead, looking at 90° scattering angle → blue sky
        let rgb = sky_color_rgb(0.0, PI / 2.0);
        assert!(
            rgb[2] > rgb[0],
            "Sky at 90° from overhead sun should be blue-dominant"
        );
    }

    #[test]
    fn test_sky_color_all_positive() {
        for sun_z in [0.0, 0.5, 1.0, 1.3] {
            for scatter in [0.1, 0.5, 1.0, PI / 2.0, PI] {
                let rgb = sky_color_rgb(sun_z, scatter);
                assert!(rgb[0] >= 0.0 && rgb[1] >= 0.0 && rgb[2] >= 0.0);
            }
        }
    }

    #[test]
    fn test_sky_radiance_blue_dominates_at_90deg() {
        let r_blue = sky_radiance_single_scatter(440e-9, 0.3, PI / 2.0);
        let r_red = sky_radiance_single_scatter(680e-9, 0.3, PI / 2.0);
        assert!(
            r_blue > r_red,
            "Blue should dominate Rayleigh scattering at 90°"
        );
    }

    // ── Sunlight color ───────────────────────────────────────────────────

    #[test]
    fn test_sunlight_overhead_near_white() {
        let rgb = sunlight_color(0.0);
        // All channels should be high (near 1.0)
        assert!(rgb[0] > 0.8, "Red transmittance at zenith: {}", rgb[0]);
        assert!(rgb[1] > 0.8, "Green transmittance at zenith: {}", rgb[1]);
        assert!(rgb[2] > 0.7, "Blue transmittance at zenith: {}", rgb[2]);
    }

    #[test]
    fn test_sunlight_sunset_is_red() {
        let rgb = sunlight_color(85.0_f64.to_radians());
        assert!(
            rgb[0] > rgb[1] && rgb[1] > rgb[2],
            "Sunset light should be R > G > B, got {:?}",
            rgb
        );
    }

    #[test]
    fn test_sunlight_deep_sunset_very_red() {
        let rgb = sunlight_color(89.0_f64.to_radians());
        // Blue should be nearly completely extinguished
        assert!(
            rgb[2] < 0.1,
            "Deep sunset should extinguish blue, got B={}",
            rgb[2]
        );
        assert!(rgb[0] > rgb[2] * 5.0, "Red should dominate at deep sunset");
    }

    #[test]
    fn test_sunlight_range() {
        for angle_deg in [0.0_f64, 30.0, 60.0, 80.0, 85.0, 89.0] {
            let rgb = sunlight_color(angle_deg.to_radians());
            for (i, &c) in rgb.iter().enumerate() {
                assert!(
                    (0.0..=1.0).contains(&c),
                    "Sunlight channel {i} out of range at {angle_deg}°: {c}"
                );
            }
        }
    }

    // ── Sunset gradient ──────────────────────────────────────────────────

    #[test]
    fn test_sunset_gradient_all_positive() {
        let rgb = sunset_gradient(85.0_f64.to_radians(), 10.0_f64.to_radians(), 0.3);
        assert!(rgb[0] >= 0.0 && rgb[1] >= 0.0 && rgb[2] >= 0.0);
    }

    #[test]
    fn test_sunset_gradient_near_sun_warm() {
        // Looking near the sun at sunset → warm colors
        let near = sunset_gradient(85.0_f64.to_radians(), 5.0_f64.to_radians(), 0.1);
        // Red should dominate near the sun at sunset
        assert!(
            near[0] > near[2],
            "Near sun at sunset should be warm: {:?}",
            near
        );
    }

    #[test]
    fn test_sunset_gradient_away_from_sun_bluer() {
        // Looking away from sun at sunset → more blue/neutral
        let near = sunset_gradient(85.0_f64.to_radians(), 20.0_f64.to_radians(), 0.2);
        let away = sunset_gradient(85.0_f64.to_radians(), 20.0_f64.to_radians(), PI / 2.0);
        // Ratio of blue/red should be higher when looking away
        let ratio_near = near[2] / (near[0] + 1e-15);
        let ratio_away = away[2] / (away[0] + 1e-15);
        assert!(
            ratio_away > ratio_near,
            "Away from sun should be bluer: near_ratio={ratio_near}, away_ratio={ratio_away}"
        );
    }

    // ── Scattering angle ─────────────────────────────────────────────────

    #[test]
    fn test_scattering_angle_same_direction() {
        let gamma = scattering_angle(0.5, 0.5, 0.0);
        assert!(gamma.abs() < EPS, "Same direction → 0 angle");
    }

    #[test]
    fn test_scattering_angle_opposite() {
        let gamma = scattering_angle(0.0, PI, 0.0);
        assert!((gamma - PI).abs() < 0.01, "Opposite directions → π");
    }

    #[test]
    fn test_scattering_angle_perpendicular() {
        // Sun at zenith (0), looking at horizon (π/2), any azimuth
        let gamma = scattering_angle(0.0, PI / 2.0, 0.0);
        assert!(
            (gamma - PI / 2.0).abs() < 0.01,
            "Zenith sun, horizon view → π/2"
        );
    }

    #[test]
    fn test_scattering_angle_range() {
        for sz in [0.0, 0.5, 1.0, PI / 2.0] {
            for vz in [0.0, 0.5, 1.0, PI / 2.0] {
                for az in [0.0, PI / 4.0, PI / 2.0, PI] {
                    let gamma = scattering_angle(sz, vz, az);
                    assert!(
                        (0.0..=PI + EPS).contains(&gamma),
                        "Angle out of range: {gamma}"
                    );
                }
            }
        }
    }
}
