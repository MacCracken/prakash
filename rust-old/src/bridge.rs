//! Cross-crate bridges — convert primitive values from other AGNOS science crates
//! into prakash optics parameters and vice versa.
//!
//! Always available — takes primitive values (f64), no science crate deps.
//!
//! # Architecture
//!
//! ```text
//! bijli (electromagnetism) ──┐
//! tara  (stellar physics)    ┼──> bridge ──> prakash optics parameters
//! badal (weather/atmosphere) ┘
//! ```

// ── Bijli bridges (electromagnetism) ───────────────────────────────────────

/// Convert EM frequency (Hz) to optical wavelength (nm).
///
/// λ = c / f, then metres → nanometres.
#[must_use]
#[inline]
pub fn em_frequency_to_wavelength_nm(frequency_hz: f64) -> f64 {
    if frequency_hz <= 0.0 {
        return 0.0;
    }
    // c = 299_792_458 m/s
    (299_792_458.0 / frequency_hz) * 1e9
}

/// Convert optical wavelength (nm) to EM frequency (Hz).
///
/// f = c / λ
#[must_use]
#[inline]
pub fn wavelength_nm_to_em_frequency(wavelength_nm: f64) -> f64 {
    if wavelength_nm <= 0.0 {
        return 0.0;
    }
    299_792_458.0 / (wavelength_nm * 1e-9)
}

/// Convert E-field amplitude (V/m) to optical intensity (W/m²).
///
/// I = ε₀ × c × E² / 2  (time-averaged Poynting vector magnitude).
#[must_use]
#[inline]
pub fn e_field_to_intensity(e_amplitude_v_per_m: f64) -> f64 {
    const EPSILON_0: f64 = 8.854_187_817e-12;
    const C: f64 = 299_792_458.0;
    0.5 * EPSILON_0 * C * e_amplitude_v_per_m * e_amplitude_v_per_m
}

/// Convert optical intensity (W/m²) to E-field amplitude (V/m).
///
/// E = sqrt(2 × I / (ε₀ × c))
#[must_use]
#[inline]
pub fn intensity_to_e_field(intensity_w_per_m2: f64) -> f64 {
    const EPSILON_0: f64 = 8.854_187_817e-12;
    const C: f64 = 299_792_458.0;
    if intensity_w_per_m2 <= 0.0 {
        return 0.0;
    }
    (2.0 * intensity_w_per_m2 / (EPSILON_0 * C)).sqrt()
}

/// Estimate refractive index from EM wavelength using Cauchy's equation
/// for common optical glass (BK7 approximation).
///
/// n(λ) ≈ A + B/λ² + C/λ⁴  where λ is in micrometres.
/// Returns refractive index (typically 1.45–1.55 for visible light).
#[must_use]
#[inline]
pub fn refractive_index_cauchy_bk7(wavelength_nm: f64) -> f64 {
    if wavelength_nm <= 0.0 {
        return 1.0;
    }
    let lambda_um = wavelength_nm * 1e-3;
    let l2 = lambda_um * lambda_um;
    // BK7 Cauchy coefficients (approximate)
    1.5046 + 0.00420 / l2 + 0.0000 / (l2 * l2)
}

// ── Tara bridges (stellar astrophysics) ────────────────────────────────────

/// Convert stellar effective temperature (K) to approximate RGB color.
///
/// Uses a piecewise blackbody-to-RGB mapping. Returns `[r, g, b]` in 0.0–1.0.
/// Valid range: ~1000 K – 40000 K.
#[must_use]
pub fn stellar_temperature_to_rgb(temperature_k: f64) -> [f64; 3] {
    // Simplified Planckian locus to RGB mapping
    let t = (temperature_k / 100.0).clamp(10.0, 400.0);

    let r = if t <= 66.0 {
        1.0
    } else {
        let x = t - 60.0;
        (329.698_727_446 * x.powf(-0.133_204_759_2) / 255.0).clamp(0.0, 1.0)
    };

    let g = if t <= 66.0 {
        let x = t;
        (99.470_802_586_1 * x.ln() - 161.119_568_166_1).clamp(0.0, 255.0) / 255.0
    } else {
        let x = t - 60.0;
        (288.122_169_528_3 * x.powf(-0.075_514_849_37) / 255.0).clamp(0.0, 1.0)
    };

    let b = if t >= 66.0 {
        1.0
    } else if t <= 19.0 {
        0.0
    } else {
        let x = t - 10.0;
        (138.517_731_223_1 * x.ln() - 305.044_792_730_7).clamp(0.0, 255.0) / 255.0
    };

    [r, g, b]
}

/// Convert stellar B-V color index to approximate temperature (K).
///
/// Uses Ballesteros (2012) formula: T = 4600 × (1/(0.92×(B-V)+1.7) + 1/(0.92×(B-V)+0.62))
#[must_use]
#[inline]
pub fn color_index_to_temperature(b_minus_v: f64) -> f64 {
    let bv = b_minus_v;
    4600.0 * (1.0 / (0.92 * bv + 1.7) + 1.0 / (0.92 * bv + 0.62))
}

/// Convert stellar spectral class letter to approximate temperature range midpoint (K).
///
/// Returns the midpoint of the temperature range for each spectral class.
#[must_use]
pub fn spectral_class_to_temperature(class: char) -> f64 {
    match class.to_ascii_uppercase() {
        'O' => 35_000.0,
        'B' => 17_500.0,
        'A' => 8_750.0,
        'F' => 6_750.0,
        'G' => 5_500.0,
        'K' => 4_250.0,
        'M' => 3_000.0,
        'L' => 1_750.0,
        'T' => 1_000.0,
        'Y' => 500.0,
        _ => 5_500.0, // default to solar
    }
}

/// Convert stellar absolute magnitude to luminosity relative to the Sun.
///
/// L/L_sun = 10^((M_sun - M) / 2.5) where M_sun = 4.83.
#[must_use]
#[inline]
pub fn absolute_magnitude_to_luminosity(absolute_magnitude: f64) -> f64 {
    const M_SUN: f64 = 4.83;
    10.0_f64.powf((M_SUN - absolute_magnitude) / 2.5)
}

// ── Badal bridges (weather/atmosphere) ──────────────────────────────────────

/// Convert atmospheric density at altitude to Rayleigh scattering coefficient
/// scale factor.
///
/// Rayleigh scattering scales linearly with number density (∝ ρ/ρ₀).
/// `density_kg_m3`: air density at altitude. Sea level ≈ 1.225 kg/m³.
/// Returns scale factor relative to sea level (1.0 at sea level).
#[must_use]
#[inline]
pub fn density_to_rayleigh_scale(density_kg_m3: f64) -> f64 {
    const SEA_LEVEL_DENSITY: f64 = 1.225;
    if density_kg_m3 <= 0.0 {
        return 0.0;
    }
    density_kg_m3 / SEA_LEVEL_DENSITY
}

/// Convert relative humidity (%) to Mie scattering coefficient scale factor.
///
/// Higher humidity increases aerosol particle size, boosting Mie scattering.
/// Returns a scale factor (1.0 at 50% RH baseline).
#[must_use]
#[inline]
pub fn humidity_to_mie_scale(relative_humidity_percent: f64) -> f64 {
    let rh = relative_humidity_percent.clamp(0.0, 100.0);
    // Mie scattering roughly doubles from 50% to 100% RH
    // due to hygroscopic aerosol growth
    if rh < 50.0 {
        0.5 + rh / 100.0
    } else {
        1.0 + (rh - 50.0) / 50.0
    }
}

/// Convert cloud cover fraction (0.0–1.0) to diffuse/direct light ratio.
///
/// Returns the fraction of incoming solar radiation that is diffuse (scattered).
/// 0.0 cloud cover → mostly direct. 1.0 cloud cover → almost all diffuse.
#[must_use]
#[inline]
pub fn cloud_cover_to_diffuse_fraction(cloud_fraction: f64) -> f64 {
    let cf = cloud_fraction.clamp(0.0, 1.0);
    // Empirical model: diffuse fraction increases roughly linearly
    // with a floor of ~0.1 (clear sky Rayleigh) and ceiling of ~0.95
    0.1 + 0.85 * cf
}

/// Convert visibility (km) to an approximate Mie extinction coefficient (1/m).
///
/// Koschmieder relation: β_ext ≈ 3.912 / V  where V is visibility in metres.
#[must_use]
#[inline]
pub fn visibility_to_extinction(visibility_km: f64) -> f64 {
    let v_m = visibility_km * 1000.0;
    if v_m <= 0.0 {
        return 0.0;
    }
    3.912 / v_m
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── Bijli bridges ──────────────────────────────────────────────────

    #[test]
    fn em_frequency_to_wavelength_visible() {
        // Green light ~550 nm → ~545 THz
        let nm = em_frequency_to_wavelength_nm(545e12);
        assert!((nm - 550.0).abs() < 5.0);
    }

    #[test]
    fn wavelength_to_frequency_roundtrip() {
        let freq = wavelength_nm_to_em_frequency(550.0);
        let nm = em_frequency_to_wavelength_nm(freq);
        assert!((nm - 550.0).abs() < 0.01);
    }

    #[test]
    fn em_frequency_zero() {
        assert_eq!(em_frequency_to_wavelength_nm(0.0), 0.0);
        assert_eq!(wavelength_nm_to_em_frequency(0.0), 0.0);
    }

    #[test]
    fn e_field_to_intensity_and_back() {
        let e = 100.0; // V/m
        let i = e_field_to_intensity(e);
        let e_back = intensity_to_e_field(i);
        assert!((e - e_back).abs() < 0.01);
    }

    #[test]
    fn intensity_zero() {
        assert_eq!(intensity_to_e_field(0.0), 0.0);
    }

    #[test]
    fn refractive_index_bk7_visible() {
        let n = refractive_index_cauchy_bk7(550.0);
        assert!(n > 1.5 && n < 1.55);
    }

    // ── Tara bridges ───────────────────────────────────────────────────

    #[test]
    fn sun_temperature_yellowish() {
        let rgb = stellar_temperature_to_rgb(5778.0);
        // Sun should be warmish white
        assert!(rgb[0] > 0.9);
        assert!(rgb[1] > 0.8);
        assert!(rgb[2] > 0.7);
    }

    #[test]
    fn hot_star_bluish() {
        let rgb = stellar_temperature_to_rgb(30000.0);
        assert!(rgb[2] > rgb[0]); // blue > red
    }

    #[test]
    fn cool_star_reddish() {
        let rgb = stellar_temperature_to_rgb(3000.0);
        assert!(rgb[0] > rgb[2]); // red > blue
    }

    #[test]
    fn color_index_sun() {
        // Sun B-V ≈ 0.65
        let t = color_index_to_temperature(0.65);
        assert!((t - 5800.0).abs() < 300.0);
    }

    #[test]
    fn spectral_class_ordering() {
        // O > B > A > F > G > K > M
        let o = spectral_class_to_temperature('O');
        let m = spectral_class_to_temperature('M');
        assert!(o > m);
    }

    #[test]
    fn spectral_class_sun_g() {
        let t = spectral_class_to_temperature('G');
        assert!((t - 5500.0).abs() < 500.0);
    }

    #[test]
    fn absolute_magnitude_sun() {
        let l = absolute_magnitude_to_luminosity(4.83);
        assert!((l - 1.0).abs() < 0.01);
    }

    #[test]
    fn absolute_magnitude_brighter() {
        // Magnitude 0 → much more luminous than Sun
        let l = absolute_magnitude_to_luminosity(0.0);
        assert!(l > 50.0);
    }

    // ── Badal bridges ──────────────────────────────────────────────────

    #[test]
    fn density_sea_level_unity() {
        let scale = density_to_rayleigh_scale(1.225);
        assert!((scale - 1.0).abs() < 0.001);
    }

    #[test]
    fn density_half_altitude() {
        let scale = density_to_rayleigh_scale(0.6125);
        assert!((scale - 0.5).abs() < 0.001);
    }

    #[test]
    fn density_zero() {
        assert_eq!(density_to_rayleigh_scale(0.0), 0.0);
    }

    #[test]
    fn humidity_baseline() {
        let scale = humidity_to_mie_scale(50.0);
        assert!((scale - 1.0).abs() < 0.01);
    }

    #[test]
    fn humidity_high_increases() {
        let high = humidity_to_mie_scale(100.0);
        assert!(high > 1.5);
    }

    #[test]
    fn cloud_cover_clear() {
        let d = cloud_cover_to_diffuse_fraction(0.0);
        assert!((d - 0.1).abs() < 0.01);
    }

    #[test]
    fn cloud_cover_overcast() {
        let d = cloud_cover_to_diffuse_fraction(1.0);
        assert!((d - 0.95).abs() < 0.01);
    }

    #[test]
    fn visibility_to_extinction_clear() {
        // 20 km visibility → low extinction
        let ext = visibility_to_extinction(20.0);
        assert!(ext > 0.0 && ext < 0.001);
    }

    #[test]
    fn visibility_zero() {
        assert_eq!(visibility_to_extinction(0.0), 0.0);
    }
}
