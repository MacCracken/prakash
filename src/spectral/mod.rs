//! Spectral math — wavelength ↔ RGB, blackbody radiation, color temperature.

use serde::{Deserialize, Serialize};
use tracing::trace;

use crate::error::{PrakashError, Result};

/// Visible light range in nanometers.
pub const VISIBLE_MIN_NM: f64 = 380.0;
pub const VISIBLE_MAX_NM: f64 = 780.0;

/// Speed of light in m/s.
pub const SPEED_OF_LIGHT: f64 = 299_792_458.0;
/// Planck's constant in J·s.
pub const PLANCK_H: f64 = 6.626_070_15e-34;
/// Boltzmann constant in J/K.
pub const BOLTZMANN_K: f64 = 1.380_649e-23;

/// Representative wavelengths for RGB channels (nanometers): R=650nm, G=550nm, B=450nm.
pub const RGB_WAVELENGTHS_NM: [f64; 3] = [650.0, 550.0, 450.0];

/// Representative wavelengths for RGB channels (meters): R=650nm, G=550nm, B=450nm.
pub const RGB_WAVELENGTHS_M: [f64; 3] = [650e-9, 550e-9, 450e-9];

/// RGB color with floating-point components (0.0–1.0).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Rgb {
    pub r: f64,
    pub g: f64,
    pub b: f64,
}

impl Rgb {
    #[must_use]
    pub const fn new(r: f64, g: f64, b: f64) -> Self {
        Self { r, g, b }
    }

    pub const BLACK: Rgb = Rgb::new(0.0, 0.0, 0.0);
    pub const WHITE: Rgb = Rgb::new(1.0, 1.0, 1.0);

    /// Clamp all components to 0.0–1.0.
    #[must_use]
    #[inline]
    pub fn clamp(self) -> Self {
        Self {
            r: self.r.clamp(0.0, 1.0),
            g: self.g.clamp(0.0, 1.0),
            b: self.b.clamp(0.0, 1.0),
        }
    }

    /// Convert to 8-bit sRGB.
    #[must_use]
    #[inline]
    pub fn to_u8(self) -> [u8; 3] {
        let c = self.clamp();
        [
            (c.r * 255.0 + 0.5) as u8,
            (c.g * 255.0 + 0.5) as u8,
            (c.b * 255.0 + 0.5) as u8,
        ]
    }

    /// Luminance (perceived brightness, Rec. 709).
    #[must_use]
    #[inline]
    pub fn luminance(self) -> f64 {
        0.2126 * self.r + 0.7152 * self.g + 0.0722 * self.b
    }
}

/// Convert a wavelength (nm) to an approximate RGB color.
///
/// Uses a piecewise linear approximation of the CIE 1931 color matching functions.
/// Only valid for visible range (380–780 nm).
#[must_use = "returns the computed RGB color"]
#[inline]
pub fn wavelength_to_rgb(wavelength_nm: f64) -> Result<Rgb> {
    trace!(wavelength_nm, "wavelength_to_rgb");
    if !(VISIBLE_MIN_NM..=VISIBLE_MAX_NM).contains(&wavelength_nm) {
        return Err(PrakashError::WavelengthOutOfRange { nm: wavelength_nm });
    }

    let w = wavelength_nm;
    let (r, g, b) = if w < 440.0 {
        (-(w - 440.0) / (440.0 - 380.0), 0.0, 1.0)
    } else if w < 490.0 {
        (0.0, (w - 440.0) / (490.0 - 440.0), 1.0)
    } else if w < 510.0 {
        (0.0, 1.0, -(w - 510.0) / (510.0 - 490.0))
    } else if w < 580.0 {
        ((w - 510.0) / (580.0 - 510.0), 1.0, 0.0)
    } else if w < 645.0 {
        (1.0, -(w - 645.0) / (645.0 - 580.0), 0.0)
    } else {
        (1.0, 0.0, 0.0)
    };

    // Intensity falloff at edges of visible spectrum
    let factor = if w < 420.0 {
        0.3 + 0.7 * (w - 380.0) / (420.0 - 380.0)
    } else if w > 700.0 {
        0.3 + 0.7 * (780.0 - w) / (780.0 - 700.0)
    } else {
        1.0
    };

    Ok(Rgb::new(r * factor, g * factor, b * factor))
}

/// First radiation constant: 2·h·c²
const PLANCK_C1: f64 = 2.0 * PLANCK_H * SPEED_OF_LIGHT * SPEED_OF_LIGHT;
/// Second radiation constant: h·c/k_B
const PLANCK_C2: f64 = PLANCK_H * SPEED_OF_LIGHT / BOLTZMANN_K;

/// Planck's law: spectral radiance of a blackbody at temperature T (Kelvin)
/// and wavelength λ (meters).
///
/// Returns spectral radiance in W·sr⁻¹·m⁻³.
#[must_use]
#[inline]
pub fn planck_radiance(wavelength_m: f64, temperature_k: f64) -> f64 {
    let lambda5 = wavelength_m.powi(5);
    PLANCK_C1 / (lambda5 * ((PLANCK_C2 / (wavelength_m * temperature_k)).exp() - 1.0))
}

/// Wien's displacement law: peak wavelength (meters) for a blackbody at temperature T (Kelvin).
///
/// λ_max = b / T, where b ≈ 2.898 × 10⁻³ m·K
#[must_use]
#[inline]
pub fn wien_peak(temperature_k: f64) -> f64 {
    2.897_771_955e-3 / temperature_k
}

/// Approximate color temperature (Kelvin) from correlated color temperature.
///
/// Converts a color temperature to an RGB color using Tanner Helland's algorithm.
#[must_use]
#[inline]
pub fn color_temperature_to_rgb(kelvin: f64) -> Rgb {
    trace!(kelvin, "color_temperature_to_rgb");
    let temp = (kelvin / 100.0).clamp(10.0, 400.0);

    let r = if temp <= 66.0 {
        1.0
    } else {
        let x = temp - 60.0;
        (329.698_727_446 * x.powf(-0.133_204_759_2) / 255.0).clamp(0.0, 1.0)
    };

    let g = if temp <= 66.0 {
        let x = temp;
        (99.470_802_586_1 * x.ln() - 161.119_568_166_1).clamp(0.0, 255.0) / 255.0
    } else {
        let x = temp - 60.0;
        (288.122_169_528_3 * x.powf(-0.075_514_849_37) / 255.0).clamp(0.0, 1.0)
    };

    let b = if temp >= 66.0 {
        1.0
    } else if temp <= 19.0 {
        0.0
    } else {
        let x = temp - 10.0;
        (138.517_731_223_1 * x.ln() - 305.044_792_730_7).clamp(0.0, 255.0) / 255.0
    };

    Rgb::new(r, g, b)
}

/// Wavelength (nm) to frequency (Hz).
#[must_use]
#[inline]
pub fn wavelength_to_frequency(wavelength_nm: f64) -> f64 {
    SPEED_OF_LIGHT / (wavelength_nm * 1e-9)
}

/// Frequency (Hz) to wavelength (nm).
#[must_use]
#[inline]
pub fn frequency_to_wavelength(frequency_hz: f64) -> f64 {
    (SPEED_OF_LIGHT / frequency_hz) * 1e9
}

/// Photon energy (Joules) for a given wavelength (nm).
#[must_use]
#[inline]
pub fn photon_energy(wavelength_nm: f64) -> f64 {
    PLANCK_H * wavelength_to_frequency(wavelength_nm)
}

/// Photon energy in electron-volts for a given wavelength (nm).
#[must_use]
#[inline]
pub fn photon_energy_ev(wavelength_nm: f64) -> f64 {
    photon_energy(wavelength_nm) / 1.602_176_634e-19
}

mod cie;
pub use cie::*;

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-6;

    #[test]
    fn test_rgb_luminance() {
        assert!((Rgb::WHITE.luminance() - 1.0).abs() < EPS);
        assert!(Rgb::BLACK.luminance().abs() < EPS);
        // Pure green has highest luminance weight
        let green = Rgb::new(0.0, 1.0, 0.0);
        let red = Rgb::new(1.0, 0.0, 0.0);
        assert!(green.luminance() > red.luminance());
    }

    #[test]
    fn test_rgb_luminance_coefficients_sum_to_one() {
        // Rec. 709: 0.2126 + 0.7152 + 0.0722 = 1.0
        assert!((Rgb::WHITE.luminance() - 1.0).abs() < EPS);
    }

    #[test]
    fn test_rgb_to_u8() {
        assert_eq!(Rgb::WHITE.to_u8(), [255, 255, 255]);
        assert_eq!(Rgb::BLACK.to_u8(), [0, 0, 0]);
        assert_eq!(Rgb::new(0.5, 0.5, 0.5).to_u8(), [128, 128, 128]);
    }

    #[test]
    fn test_rgb_to_u8_clamps_overrange() {
        let c = Rgb::new(2.0, -1.0, 0.5);
        let u = c.to_u8();
        assert_eq!(u[0], 255);
        assert_eq!(u[1], 0);
        assert_eq!(u[2], 128);
    }

    #[test]
    fn test_rgb_clamp() {
        let c = Rgb::new(1.5, -0.5, 0.5).clamp();
        assert!((c.r - 1.0).abs() < EPS);
        assert!(c.g.abs() < EPS);
        assert!((c.b - 0.5).abs() < EPS);
    }

    #[test]
    fn test_rgb_clamp_already_valid() {
        let c = Rgb::new(0.3, 0.6, 0.9);
        let clamped = c.clamp();
        assert!((c.r - clamped.r).abs() < EPS);
        assert!((c.g - clamped.g).abs() < EPS);
        assert!((c.b - clamped.b).abs() < EPS);
    }

    #[test]
    fn test_rgb_serde_roundtrip() {
        let c = Rgb::new(0.1, 0.2, 0.3);
        let json = serde_json::to_string(&c).unwrap();
        let back: Rgb = serde_json::from_str(&json).unwrap();
        assert!((back.r - c.r).abs() < EPS);
        assert!((back.g - c.g).abs() < EPS);
        assert!((back.b - c.b).abs() < EPS);
    }

    // ── Wavelength → RGB tests ────────────────────────────────────────────

    #[test]
    fn test_wavelength_to_rgb_red() {
        let rgb = wavelength_to_rgb(650.0).unwrap();
        assert!(rgb.r > 0.8);
        assert!(rgb.g < 0.2);
        assert!(rgb.b < 0.1);
    }

    #[test]
    fn test_wavelength_to_rgb_green() {
        let rgb = wavelength_to_rgb(530.0).unwrap();
        assert!(rgb.g > 0.8);
    }

    #[test]
    fn test_wavelength_to_rgb_blue() {
        let rgb = wavelength_to_rgb(450.0).unwrap();
        assert!(rgb.b > 0.8);
    }

    #[test]
    fn test_wavelength_to_rgb_yellow() {
        let rgb = wavelength_to_rgb(580.0).unwrap();
        assert!(rgb.r > 0.5);
        assert!(rgb.g > 0.5);
        assert!(rgb.b < 0.1);
    }

    #[test]
    fn test_wavelength_to_rgb_cyan() {
        let rgb = wavelength_to_rgb(490.0).unwrap();
        assert!(rgb.g > 0.5);
        assert!(rgb.b > 0.5);
    }

    #[test]
    fn test_wavelength_out_of_range() {
        assert!(wavelength_to_rgb(300.0).is_err());
        assert!(wavelength_to_rgb(800.0).is_err());
        assert!(wavelength_to_rgb(0.0).is_err());
        assert!(wavelength_to_rgb(-100.0).is_err());
        assert!(wavelength_to_rgb(f64::NAN).is_err());
    }

    #[test]
    fn test_wavelength_boundary() {
        assert!(wavelength_to_rgb(380.0).is_ok());
        assert!(wavelength_to_rgb(780.0).is_ok());
        assert!(wavelength_to_rgb(379.9).is_err());
        assert!(wavelength_to_rgb(780.1).is_err());
    }

    #[test]
    fn test_wavelength_to_rgb_all_components_valid() {
        // Sweep visible spectrum — all components should be in [0, 1]
        let mut nm = 380.0;
        while nm <= 780.0 {
            let rgb = wavelength_to_rgb(nm).unwrap();
            assert!(
                (0.0..=1.0).contains(&rgb.r),
                "r out of range at {nm}nm: {}",
                rgb.r
            );
            assert!(
                (0.0..=1.0).contains(&rgb.g),
                "g out of range at {nm}nm: {}",
                rgb.g
            );
            assert!(
                (0.0..=1.0).contains(&rgb.b),
                "b out of range at {nm}nm: {}",
                rgb.b
            );
            nm += 5.0;
        }
    }

    #[test]
    fn test_wavelength_edge_falloff() {
        // At spectrum edges, intensity should be reduced
        let edge = wavelength_to_rgb(390.0).unwrap();
        let center = wavelength_to_rgb(500.0).unwrap();
        let edge_lum = edge.r + edge.g + edge.b;
        let center_lum = center.r + center.g + center.b;
        assert!(edge_lum < center_lum, "Edge should have lower intensity");
    }

    // ── Blackbody & Wien tests ────────────────────────────────────────────

    #[test]
    fn test_planck_radiance_positive() {
        let r = planck_radiance(500e-9, 5778.0);
        assert!(r > 0.0);
    }

    #[test]
    fn test_planck_radiance_peaks_at_wien() {
        let temp = 5778.0;
        let peak_wl = wien_peak(temp);
        let r_peak = planck_radiance(peak_wl, temp);
        let r_shorter = planck_radiance(peak_wl * 0.5, temp);
        let r_longer = planck_radiance(peak_wl * 2.0, temp);
        assert!(r_peak > r_shorter, "Peak should be > shorter wavelength");
        assert!(r_peak > r_longer, "Peak should be > longer wavelength");
    }

    #[test]
    fn test_planck_hotter_means_more_radiance() {
        let wl = 500e-9;
        assert!(planck_radiance(wl, 6000.0) > planck_radiance(wl, 3000.0));
    }

    #[test]
    fn test_wien_peak_sun() {
        let peak = wien_peak(5778.0);
        assert!((peak * 1e9 - 502.0).abs() < 5.0);
    }

    #[test]
    fn test_wien_peak_hot() {
        assert!(wien_peak(10000.0) < wien_peak(5000.0));
    }

    #[test]
    fn test_wien_inverse_proportionality() {
        // λ_max ∝ 1/T
        let p1 = wien_peak(3000.0);
        let p2 = wien_peak(6000.0);
        assert!((p1 / p2 - 2.0).abs() < EPS);
    }

    // ── Color temperature tests ───────────────────────────────────────────

    #[test]
    fn test_color_temp_warm() {
        let rgb = color_temperature_to_rgb(2700.0);
        assert!(rgb.r > rgb.b);
    }

    #[test]
    fn test_color_temp_cool() {
        let rgb = color_temperature_to_rgb(10000.0);
        assert!(rgb.b > rgb.r * 0.5);
    }

    #[test]
    fn test_color_temp_daylight() {
        let rgb = color_temperature_to_rgb(6500.0);
        assert!(rgb.r > 0.8);
        assert!(rgb.g > 0.8);
        assert!(rgb.b > 0.8);
    }

    #[test]
    fn test_color_temp_candle() {
        let rgb = color_temperature_to_rgb(1800.0);
        assert!(rgb.r > 0.8);
        assert!(rgb.b < 0.3);
    }

    #[test]
    fn test_color_temp_range_produces_valid_rgb() {
        for k in (1000..=20000).step_by(500) {
            let rgb = color_temperature_to_rgb(k as f64);
            assert!((0.0..=1.0).contains(&rgb.r), "r invalid at {k}K");
            assert!((0.0..=1.0).contains(&rgb.g), "g invalid at {k}K");
            assert!((0.0..=1.0).contains(&rgb.b), "b invalid at {k}K");
        }
    }

    // ── Frequency & energy tests ──────────────────────────────────────────

    #[test]
    fn test_wavelength_frequency_roundtrip() {
        for nm in [380.0, 500.0, 550.0, 650.0, 780.0] {
            let freq = wavelength_to_frequency(nm);
            let back = frequency_to_wavelength(freq);
            assert!((back - nm).abs() < EPS, "Roundtrip failed at {nm}nm");
        }
    }

    #[test]
    fn test_shorter_wavelength_higher_frequency() {
        let f_blue = wavelength_to_frequency(450.0);
        let f_red = wavelength_to_frequency(700.0);
        assert!(f_blue > f_red);
    }

    #[test]
    fn test_photon_energy_visible() {
        let e = photon_energy_ev(550.0);
        assert!((e - 2.25).abs() < 0.1);
    }

    #[test]
    fn test_photon_energy_blue_greater_than_red() {
        let e_blue = photon_energy(450.0);
        let e_red = photon_energy(700.0);
        assert!(e_blue > e_red);
    }

    #[test]
    fn test_photon_energy_ev_positive() {
        let e = photon_energy_ev(550.0);
        assert!(e > 0.0);
    }

    #[test]
    fn test_photon_energy_joules_tiny() {
        let e = photon_energy(550.0);
        assert!(e > 0.0);
        assert!(e < 1e-17); // should be on order of 10^-19 J
    }

    // ── Physical constants sanity ─────────────────────────────────────────

    #[test]
    fn test_speed_of_light_value() {
        assert!((SPEED_OF_LIGHT - 2.998e8).abs() < 1e5);
    }

    #[test]
    fn test_planck_constant_value() {
        assert!((PLANCK_H - 6.626e-34).abs() < 1e-36);
    }

    #[test]
    fn test_boltzmann_constant_value() {
        assert!((BOLTZMANN_K - 1.381e-23).abs() < 1e-25);
    }
}
