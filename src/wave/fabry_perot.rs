//! Fabry-Perot interferometer: finesse, transmittance, free spectral range, resolving power.

use std::f64::consts::PI;

/// Fabry-Perot coefficient of finesse.
///
/// F = 4R / (1 - R)²
///
/// `reflectance` is the mirror reflectivity (0.0-1.0).
#[must_use]
#[inline]
pub fn fabry_perot_finesse_coefficient(reflectance: f64) -> f64 {
    4.0 * reflectance / ((1.0 - reflectance) * (1.0 - reflectance))
}

/// Fabry-Perot finesse (sharpness of fringes).
///
/// F = π·√R / (1 - R)
///
/// Higher finesse means sharper transmission peaks.
#[must_use]
#[inline]
pub fn fabry_perot_finesse(reflectance: f64) -> f64 {
    PI * reflectance.sqrt() / (1.0 - reflectance)
}

/// Fabry-Perot transmittance (Airy function).
///
/// T = 1 / (1 + F·sin²(δ/2))
///
/// where δ = 4π·n·d·cos(θ)/λ is the round-trip phase.
///
/// `n` = refractive index of cavity, `thickness` = mirror separation,
/// `wavelength` and `thickness` in same units, `angle` = incidence angle (radians),
/// `reflectance` = mirror reflectivity.
#[must_use]
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

/// Fabry-Perot free spectral range (FSR).
///
/// Δν = c / (2·n·d)  (in Hz)
///
/// `thickness_m` in meters, `n` = refractive index of cavity.
/// Returns FSR in Hz.
#[must_use]
#[inline]
pub fn fabry_perot_fsr(thickness_m: f64, n: f64) -> f64 {
    299_792_458.0 / (2.0 * n * thickness_m)
}

/// Fabry-Perot free spectral range in wavelength units.
///
/// Δλ = λ² / (2·n·d)
///
/// `wavelength` and `thickness` in same units.
#[must_use]
#[inline]
pub fn fabry_perot_fsr_wavelength(wavelength: f64, thickness: f64, n: f64) -> f64 {
    wavelength * wavelength / (2.0 * n * thickness)
}

/// Fabry-Perot resolving power.
///
/// R = m · F where m = 2·n·d/λ (order number) and F is finesse.
///
/// `wavelength` and `thickness` in same units.
#[must_use]
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
