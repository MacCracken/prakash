//! Fresnel equations (real and complex), Brewster's angle, and Beer-Lambert attenuation.

use serde::{Deserialize, Serialize};

use super::critical_angle;
use crate::error::{PrakashError, Result};

// ── Fresnel Equations (real refractive index) ────────────────────────────────

/// Fresnel reflectance for s-polarized light (TE mode).
#[must_use]
#[inline]
pub fn fresnel_s(n1: f64, n2: f64, cos_i: f64, cos_t: f64) -> f64 {
    let num = n1 * cos_i - n2 * cos_t;
    let den = n1 * cos_i + n2 * cos_t;
    if den.abs() < 1e-15 {
        return 1.0;
    }
    let r = num / den;
    r * r
}

/// Fresnel reflectance for p-polarized light (TM mode).
#[must_use]
#[inline]
pub fn fresnel_p(n1: f64, n2: f64, cos_i: f64, cos_t: f64) -> f64 {
    let num = n2 * cos_i - n1 * cos_t;
    let den = n2 * cos_i + n1 * cos_t;
    if den.abs() < 1e-15 {
        return 1.0;
    }
    let r = num / den;
    r * r
}

/// Average Fresnel reflectance for unpolarized light.
///
/// Returns the fraction of light reflected (0.0 = none, 1.0 = total reflection).
///
#[must_use = "returns the Fresnel reflectance"]
#[inline]
pub fn fresnel_unpolarized(n1: f64, n2: f64, incident_angle: f64) -> Result<f64> {
    fresnel_unpolarized_impl(n1, n2, incident_angle)
}

#[inline]
fn fresnel_unpolarized_impl(n1: f64, n2: f64, incident_angle: f64) -> Result<f64> {
    let (sin_i, cos_i) = incident_angle.sin_cos();
    let sin_t = (n1 / n2) * sin_i;
    if sin_t.abs() > 1.0 {
        let critical = critical_angle(n1, n2)?;
        return Err(PrakashError::TotalInternalReflection {
            angle_deg: incident_angle.to_degrees(),
            critical_deg: critical.to_degrees(),
            n1,
            n2,
        });
    }
    let cos_t = (1.0 - sin_t * sin_t).sqrt();
    Ok(0.5 * (fresnel_s(n1, n2, cos_i, cos_t) + fresnel_p(n1, n2, cos_i, cos_t)))
}

/// Fresnel reflectance at normal incidence (θ = 0).
///
/// Simplified: R = ((n1 - n2) / (n1 + n2))²
///
/// ```
/// # use prakash::ray::fresnel_normal;
/// let r = fresnel_normal(1.0, 1.52); // air→glass ≈ 4%
/// assert!((r - 0.04).abs() < 0.01);
/// ```
#[must_use]
#[inline]
pub fn fresnel_normal(n1: f64, n2: f64) -> f64 {
    fresnel_normal_impl(n1, n2)
}

#[inline]
fn fresnel_normal_impl(n1: f64, n2: f64) -> f64 {
    let r = (n1 - n2) / (n1 + n2);
    r * r
}

// ── Brewster's Angle ────────────────────────────────────────────────────────

/// Brewster's angle — the angle at which reflected light is fully polarized.
///
/// At this angle, p-polarized reflectance is zero.
///
#[must_use]
#[inline]
pub fn brewster_angle(n1: f64, n2: f64) -> f64 {
    brewster_angle_impl(n1, n2)
}

#[inline]
fn brewster_angle_impl(n1: f64, n2: f64) -> f64 {
    (n2 / n1).atan()
}

// ── Complex Fresnel (absorbing media) ────────────────────────────────────────

/// Medium with complex refractive index for absorbing materials (metals, semiconductors).
///
/// The complex refractive index is **n~ = n + ik** where:
/// - `n` = real part (phase velocity ratio)
/// - `k` = extinction coefficient (absorption)
///
/// For dielectrics, k = 0. For metals, k is typically 1–10.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ComplexMedium {
    /// Real part of the refractive index.
    pub n: f64,
    /// Extinction coefficient (imaginary part of the refractive index).
    pub k: f64,
    /// Human-readable name.
    pub name: &'static str,
}

impl ComplexMedium {
    /// Gold at 550nm (visible).
    pub const GOLD_550NM: Self = Self {
        n: 0.43,
        k: 2.46,
        name: "gold (550nm)",
    };
    /// Silver at 550nm (visible).
    pub const SILVER_550NM: Self = Self {
        n: 0.059,
        k: 3.33,
        name: "silver (550nm)",
    };
    /// Copper at 550nm (visible).
    pub const COPPER_550NM: Self = Self {
        n: 1.13,
        k: 2.60,
        name: "copper (550nm)",
    };
    /// Aluminum at 550nm (visible).
    pub const ALUMINUM_550NM: Self = Self {
        n: 0.96,
        k: 6.69,
        name: "aluminum (550nm)",
    };

    /// Create a dielectric (non-absorbing) complex medium from a real refractive index.
    #[must_use]
    #[inline]
    pub const fn dielectric(n: f64, name: &'static str) -> Self {
        Self { n, k: 0.0, name }
    }
}

/// Fresnel reflectance at normal incidence for absorbing media.
///
/// R = ((n₁ − n₂)² + k₂²) / ((n₁ + n₂)² + k₂²)
///
/// `n1` = incident medium (real, e.g., air = 1.0),
/// `medium` = absorbing medium with complex refractive index.
#[must_use]
#[inline]
pub fn fresnel_normal_complex(n1: f64, medium: &ComplexMedium) -> f64 {
    let dn = n1 - medium.n;
    let dp = n1 + medium.n;
    let k2 = medium.k * medium.k;
    (dn * dn + k2) / (dp * dp + k2)
}

/// Fresnel reflectance for s-polarized light on an absorbing medium.
///
/// Uses the full complex Snell's law to compute the transmitted complex angle,
/// then evaluates the Fresnel coefficient for s-polarization.
///
/// Returns `R_s = |r_s|²`.
#[must_use]
#[inline]
pub fn fresnel_s_complex(n1: f64, medium: &ComplexMedium, incident_angle: f64) -> f64 {
    let (sin_i, cos_i) = incident_angle.sin_cos();
    let sin2_i = sin_i * sin_i;

    // Complex cos(theta_t) via Snell's law: cos²(θt) = 1 - (n1/n~)² sin²(θi)
    let n2_sq = medium.n * medium.n - medium.k * medium.k; // Re(n~²)
    let n2_im = 2.0 * medium.n * medium.k; // Im(n~²)
    let denom = n2_sq * n2_sq + n2_im * n2_im; // |n~²|²

    let n1_sin2 = n1 * n1 * sin2_i;
    let u = 1.0 - n1_sin2 * n2_sq / denom; // Re(cos²θt)
    let v = n1_sin2 * n2_im / denom; // Im(cos²θt)

    // cos(θt) = sqrt(u + iv) — principal square root
    let mag = (u * u + v * v).sqrt();
    let cos_t_re = ((mag + u) / 2.0).sqrt();
    let cos_t_im = if v >= 0.0 {
        ((mag - u) / 2.0).sqrt()
    } else {
        -((mag - u) / 2.0).sqrt()
    };

    // r_s = (n1 cos θi - n~ cos θt) / (n1 cos θi + n~ cos θt)
    let nt_re = medium.n * cos_t_re - medium.k * cos_t_im;
    let nt_im = medium.n * cos_t_im + medium.k * cos_t_re;

    let num_re = n1 * cos_i - nt_re;
    let num_im = -nt_im;
    let den_re = n1 * cos_i + nt_re;
    let den_im = nt_im;

    (num_re * num_re + num_im * num_im) / (den_re * den_re + den_im * den_im + 1e-30)
}

/// Fresnel reflectance for p-polarized light on an absorbing medium.
///
/// Uses the full complex Snell's law.
/// Returns `R_p = |r_p|²`.
#[must_use]
#[inline]
pub fn fresnel_p_complex(n1: f64, medium: &ComplexMedium, incident_angle: f64) -> f64 {
    let (sin_i, cos_i) = incident_angle.sin_cos();
    let sin2_i = sin_i * sin_i;

    let n2_sq = medium.n * medium.n - medium.k * medium.k;
    let n2_im = 2.0 * medium.n * medium.k;
    let denom = n2_sq * n2_sq + n2_im * n2_im;

    let n1_sin2 = n1 * n1 * sin2_i;
    let u = 1.0 - n1_sin2 * n2_sq / denom;
    let v = n1_sin2 * n2_im / denom;

    let mag = (u * u + v * v).sqrt();
    let cos_t_re = ((mag + u) / 2.0).sqrt();
    let cos_t_im = if v >= 0.0 {
        ((mag - u) / 2.0).sqrt()
    } else {
        -((mag - u) / 2.0).sqrt()
    };

    let nci_re = medium.n * cos_i;
    let nci_im = medium.k * cos_i;

    let num_re = nci_re - n1 * cos_t_re;
    let num_im = nci_im - n1 * cos_t_im;
    let den_re = nci_re + n1 * cos_t_re;
    let den_im = nci_im + n1 * cos_t_im;

    (num_re * num_re + num_im * num_im) / (den_re * den_re + den_im * den_im + 1e-30)
}

/// Average Fresnel reflectance for unpolarized light on an absorbing medium.
///
/// R = (R_s + R_p) / 2
#[must_use]
#[inline]
pub fn fresnel_unpolarized_complex(n1: f64, medium: &ComplexMedium, incident_angle: f64) -> f64 {
    0.5 * (fresnel_s_complex(n1, medium, incident_angle)
        + fresnel_p_complex(n1, medium, incident_angle))
}

// ── Attenuation ─────────────────────────────────────────────────────────────

/// Beer-Lambert law: intensity after traveling distance `d` through a medium
/// with absorption coefficient `alpha`.
///
/// I = I₀ · exp(−α · d)
///
/// Related: [`pbr::volume_transmittance`](crate::pbr::volume_transmittance) returns
/// the transmittance factor alone (equivalent to `beer_lambert(1.0, σt, d)`).
#[must_use]
#[inline]
pub fn beer_lambert(intensity: f64, alpha: f64, distance: f64) -> f64 {
    intensity * (-alpha * distance).exp()
}
