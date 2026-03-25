//! Circular aperture diffraction: Bessel J₁, Airy pattern, Rayleigh criterion.

use std::f64::consts::PI;

/// First-order Bessel function J₁(x) — rational approximation.
///
/// Accurate to ~1e-7 for all x. Uses polynomial approximation for |x| < 8
/// and asymptotic expansion for |x| >= 8.
#[must_use]
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
/// `wavelength` and `aperture_diameter` in same units.
/// `angle` in radians. Returns intensity relative to `i0`.
#[must_use]
#[inline]
pub fn airy_pattern(wavelength: f64, aperture_diameter: f64, angle: f64, i0: f64) -> f64 {
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
#[must_use]
#[inline]
pub fn airy_first_zero(wavelength: f64, aperture_diameter: f64) -> f64 {
    1.219_670_5 * wavelength / aperture_diameter
}

/// Rayleigh resolution criterion.
///
/// Two point sources are just resolved when the central maximum of one
/// falls on the first zero of the other: θ_R = 1.22·λ/D
///
/// Returns minimum resolvable angle in radians.
///
/// This is the same formula as [`lens::diffraction_limit`](crate::lens::diffraction_limit),
/// provided here for wave-optics contexts.
#[must_use]
#[inline]
pub fn rayleigh_criterion(wavelength: f64, aperture_diameter: f64) -> f64 {
    1.22 * wavelength / aperture_diameter
}
