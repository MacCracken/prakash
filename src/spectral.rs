//! Spectral math — wavelength ↔ RGB, blackbody radiation, color temperature.

use serde::{Deserialize, Serialize};

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

/// RGB color with floating-point components (0.0–1.0).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Rgb {
    pub r: f64,
    pub g: f64,
    pub b: f64,
}

impl Rgb {
    pub const fn new(r: f64, g: f64, b: f64) -> Self {
        Self { r, g, b }
    }

    pub const BLACK: Rgb = Rgb::new(0.0, 0.0, 0.0);
    pub const WHITE: Rgb = Rgb::new(1.0, 1.0, 1.0);

    /// Clamp all components to 0.0–1.0.
    #[inline]
    pub fn clamp(self) -> Self {
        Self {
            r: self.r.clamp(0.0, 1.0),
            g: self.g.clamp(0.0, 1.0),
            b: self.b.clamp(0.0, 1.0),
        }
    }

    /// Convert to 8-bit sRGB.
    #[inline]
    pub fn to_u8(self) -> [u8; 3] {
        let c = self.clamp();
        [
            (c.r * 255.0) as u8,
            (c.g * 255.0) as u8,
            (c.b * 255.0) as u8,
        ]
    }

    /// Luminance (perceived brightness, Rec. 709).
    #[inline]
    pub fn luminance(self) -> f64 {
        0.2126 * self.r + 0.7152 * self.g + 0.0722 * self.b
    }
}

/// Convert a wavelength (nm) to an approximate RGB color.
///
/// Uses a piecewise linear approximation of the CIE 1931 color matching functions.
/// Only valid for visible range (380–780 nm).
pub fn wavelength_to_rgb(wavelength_nm: f64) -> Result<Rgb> {
    if wavelength_nm < VISIBLE_MIN_NM || wavelength_nm > VISIBLE_MAX_NM {
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

/// Planck's law: spectral radiance of a blackbody at temperature T (Kelvin)
/// and wavelength λ (meters).
///
/// Returns spectral radiance in W·sr⁻¹·m⁻³.
#[inline]
pub fn planck_radiance(wavelength_m: f64, temperature_k: f64) -> f64 {
    let c1 = 2.0 * PLANCK_H * SPEED_OF_LIGHT * SPEED_OF_LIGHT;
    let c2 = PLANCK_H * SPEED_OF_LIGHT / (BOLTZMANN_K * temperature_k);
    let lambda5 = wavelength_m.powi(5);
    c1 / (lambda5 * ((c2 / wavelength_m).exp() - 1.0))
}

/// Wien's displacement law: peak wavelength (meters) for a blackbody at temperature T (Kelvin).
///
/// λ_max = b / T, where b ≈ 2.898 × 10⁻³ m·K
#[inline]
pub fn wien_peak(temperature_k: f64) -> f64 {
    2.897_771_955e-3 / temperature_k
}

/// Approximate color temperature (Kelvin) from correlated color temperature.
///
/// Converts a color temperature to an RGB color using Tanner Helland's algorithm.
pub fn color_temperature_to_rgb(kelvin: f64) -> Rgb {
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
#[inline]
pub fn wavelength_to_frequency(wavelength_nm: f64) -> f64 {
    SPEED_OF_LIGHT / (wavelength_nm * 1e-9)
}

/// Frequency (Hz) to wavelength (nm).
#[inline]
pub fn frequency_to_wavelength(frequency_hz: f64) -> f64 {
    (SPEED_OF_LIGHT / frequency_hz) * 1e9
}

/// Photon energy (Joules) for a given wavelength (nm).
#[inline]
pub fn photon_energy(wavelength_nm: f64) -> f64 {
    PLANCK_H * wavelength_to_frequency(wavelength_nm)
}

/// Photon energy in electron-volts for a given wavelength (nm).
#[inline]
pub fn photon_energy_ev(wavelength_nm: f64) -> f64 {
    photon_energy(wavelength_nm) / 1.602_176_634e-19
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-6;

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
    fn test_wavelength_out_of_range() {
        assert!(wavelength_to_rgb(300.0).is_err());
        assert!(wavelength_to_rgb(800.0).is_err());
    }

    #[test]
    fn test_wavelength_boundary() {
        assert!(wavelength_to_rgb(380.0).is_ok());
        assert!(wavelength_to_rgb(780.0).is_ok());
    }

    #[test]
    fn test_planck_radiance_positive() {
        let r = planck_radiance(500e-9, 5778.0); // Sun's peak
        assert!(r > 0.0);
    }

    #[test]
    fn test_wien_peak_sun() {
        let peak = wien_peak(5778.0);
        // Sun's peak ≈ 502 nm
        assert!((peak * 1e9 - 502.0).abs() < 5.0);
    }

    #[test]
    fn test_wien_peak_hot() {
        // Hotter → shorter wavelength
        assert!(wien_peak(10000.0) < wien_peak(5000.0));
    }

    #[test]
    fn test_color_temp_warm() {
        let rgb = color_temperature_to_rgb(2700.0); // warm white
        assert!(rgb.r > rgb.b); // warm = more red than blue
    }

    #[test]
    fn test_color_temp_cool() {
        let rgb = color_temperature_to_rgb(10000.0); // cool/blue
        assert!(rgb.b > rgb.r * 0.5); // blue-ish
    }

    #[test]
    fn test_color_temp_daylight() {
        let rgb = color_temperature_to_rgb(6500.0); // D65 daylight
        // Should be roughly white
        assert!(rgb.r > 0.8);
        assert!(rgb.g > 0.8);
        assert!(rgb.b > 0.8);
    }

    #[test]
    fn test_wavelength_frequency_roundtrip() {
        let nm = 550.0;
        let freq = wavelength_to_frequency(nm);
        let back = frequency_to_wavelength(freq);
        assert!((back - nm).abs() < EPS);
    }

    #[test]
    fn test_photon_energy_visible() {
        let e = photon_energy_ev(550.0);
        // Green light ≈ 2.25 eV
        assert!((e - 2.25).abs() < 0.1);
    }

    #[test]
    fn test_rgb_luminance() {
        assert!((Rgb::WHITE.luminance() - 1.0).abs() < EPS);
        assert!(Rgb::BLACK.luminance().abs() < EPS);
    }

    #[test]
    fn test_rgb_to_u8() {
        assert_eq!(Rgb::WHITE.to_u8(), [255, 255, 255]);
        assert_eq!(Rgb::BLACK.to_u8(), [0, 0, 0]);
        assert_eq!(Rgb::new(0.5, 0.5, 0.5).to_u8(), [127, 127, 127]);
    }

    #[test]
    fn test_rgb_clamp() {
        let c = Rgb::new(1.5, -0.5, 0.5).clamp();
        assert!((c.r - 1.0).abs() < EPS);
        assert!(c.g.abs() < EPS);
        assert!((c.b - 0.5).abs() < EPS);
    }
}
