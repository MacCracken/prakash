//! Geometric (ray) optics — reflection, refraction, Snell's law, Fresnel equations.
//!
//! All angles in radians unless suffixed with `_deg`.

use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

use crate::error::{PrakashError, Result};

// ── Refractive indices for common materials ─────────────────────────────────

/// Common material refractive indices at ~589 nm (sodium D line).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Medium {
    /// Refractive index (n >= 1.0).
    pub n: f64,
    /// Human-readable name.
    pub name: &'static str,
}

impl Medium {
    pub const VACUUM: Medium = Medium { n: 1.0, name: "vacuum" };
    pub const AIR: Medium = Medium { n: 1.000293, name: "air" };
    pub const WATER: Medium = Medium { n: 1.333, name: "water" };
    pub const GLASS: Medium = Medium { n: 1.52, name: "glass" };
    pub const CROWN_GLASS: Medium = Medium { n: 1.523, name: "crown glass" };
    pub const FLINT_GLASS: Medium = Medium { n: 1.62, name: "flint glass" };
    pub const DIAMOND: Medium = Medium { n: 2.417, name: "diamond" };
    pub const ICE: Medium = Medium { n: 1.31, name: "ice" };
    pub const QUARTZ: Medium = Medium { n: 1.544, name: "quartz" };
    pub const SAPPHIRE: Medium = Medium { n: 1.77, name: "sapphire" };
    pub const ACRYLIC: Medium = Medium { n: 1.49, name: "acrylic" };
    pub const POLYCARBONATE: Medium = Medium { n: 1.585, name: "polycarbonate" };

    /// Custom medium with a given refractive index.
    pub fn custom(n: f64, name: &'static str) -> Result<Self> {
        if n < 1.0 {
            return Err(PrakashError::InvalidRefractiveIndex { n });
        }
        Ok(Self { n, name })
    }
}

// ── Snell's Law ─────────────────────────────────────────────────────────────

/// Apply Snell's law: n1 * sin(θ1) = n2 * sin(θ2).
///
/// Returns the refracted angle in radians, or `TotalInternalReflection` error
/// if the angle exceeds the critical angle.
#[inline]
pub fn snell(n1: f64, n2: f64, incident_angle: f64) -> Result<f64> {
    let sin_t = (n1 / n2) * incident_angle.sin();
    if sin_t.abs() > 1.0 {
        let critical = critical_angle(n1, n2)?;
        return Err(PrakashError::TotalInternalReflection {
            angle_deg: incident_angle.to_degrees(),
            critical_deg: critical.to_degrees(),
            n1,
            n2,
        });
    }
    Ok(sin_t.asin())
}

/// Critical angle for total internal reflection (n1 > n2).
///
/// Returns the angle in radians beyond which all light is reflected.
#[inline]
pub fn critical_angle(n1: f64, n2: f64) -> Result<f64> {
    if n1 <= n2 {
        return Err(PrakashError::InvalidParameter {
            reason: format!("critical angle requires n1 > n2, got n1={n1}, n2={n2}"),
        });
    }
    Ok((n2 / n1).asin())
}

// ── Reflection ──────────────────────────────────────────────────────────────

/// Angle of reflection equals angle of incidence.
#[inline]
pub fn reflect_angle(incident_angle: f64) -> f64 {
    incident_angle
}

/// Reflect a 2D direction vector about a surface normal.
///
/// Both vectors should be normalized.
#[inline]
pub fn reflect_2d(direction: [f64; 2], normal: [f64; 2]) -> [f64; 2] {
    let dot = direction[0] * normal[0] + direction[1] * normal[1];
    [
        direction[0] - 2.0 * dot * normal[0],
        direction[1] - 2.0 * dot * normal[1],
    ]
}

/// Reflect a 3D direction vector about a surface normal.
#[inline]
pub fn reflect_3d(direction: [f64; 3], normal: [f64; 3]) -> [f64; 3] {
    let dot = direction[0] * normal[0] + direction[1] * normal[1] + direction[2] * normal[2];
    [
        direction[0] - 2.0 * dot * normal[0],
        direction[1] - 2.0 * dot * normal[1],
        direction[2] - 2.0 * dot * normal[2],
    ]
}

// ── Fresnel Equations ───────────────────────────────────────────────────────

/// Fresnel reflectance for s-polarized light (TE mode).
#[inline]
pub fn fresnel_s(n1: f64, n2: f64, cos_i: f64, cos_t: f64) -> f64 {
    let num = n1 * cos_i - n2 * cos_t;
    let den = n1 * cos_i + n2 * cos_t;
    if den.abs() < 1e-15 {
        return 1.0;
    }
    (num / den).powi(2)
}

/// Fresnel reflectance for p-polarized light (TM mode).
#[inline]
pub fn fresnel_p(n1: f64, n2: f64, cos_i: f64, cos_t: f64) -> f64 {
    let num = n2 * cos_i - n1 * cos_t;
    let den = n2 * cos_i + n1 * cos_t;
    if den.abs() < 1e-15 {
        return 1.0;
    }
    (num / den).powi(2)
}

/// Average Fresnel reflectance for unpolarized light.
///
/// Returns the fraction of light reflected (0.0 = none, 1.0 = total reflection).
#[inline]
pub fn fresnel_unpolarized(n1: f64, n2: f64, incident_angle: f64) -> Result<f64> {
    let cos_i = incident_angle.cos();
    let refracted = snell(n1, n2, incident_angle)?;
    let cos_t = refracted.cos();
    Ok(0.5 * (fresnel_s(n1, n2, cos_i, cos_t) + fresnel_p(n1, n2, cos_i, cos_t)))
}

/// Fresnel reflectance at normal incidence (θ = 0).
///
/// Simplified: R = ((n1 - n2) / (n1 + n2))²
#[inline]
pub fn fresnel_normal(n1: f64, n2: f64) -> f64 {
    let r = (n1 - n2) / (n1 + n2);
    r * r
}

// ── Brewster's Angle ────────────────────────────────────────────────────────

/// Brewster's angle — the angle at which reflected light is fully polarized.
///
/// At this angle, p-polarized reflectance is zero.
#[inline]
pub fn brewster_angle(n1: f64, n2: f64) -> f64 {
    (n2 / n1).atan()
}

// ── Attenuation ─────────────────────────────────────────────────────────────

/// Beer-Lambert law: intensity after traveling distance `d` through a medium
/// with absorption coefficient `alpha`.
///
/// I = I0 * exp(-alpha * d)
#[inline]
pub fn beer_lambert(intensity: f64, alpha: f64, distance: f64) -> f64 {
    intensity * (-alpha * distance).exp()
}

// ── Utility ─────────────────────────────────────────────────────────────────

/// Convert degrees to radians.
#[inline]
pub fn deg_to_rad(deg: f64) -> f64 {
    deg * PI / 180.0
}

/// Convert radians to degrees.
#[inline]
pub fn rad_to_deg(rad: f64) -> f64 {
    rad * 180.0 / PI
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-6;

    #[test]
    fn test_snell_air_to_glass() {
        let angle_i = deg_to_rad(30.0);
        let angle_t = snell(Medium::AIR.n, Medium::GLASS.n, angle_i).unwrap();
        // sin(30°) = 0.5, sin(θt) = 0.5 * 1.000293 / 1.52 ≈ 0.329
        assert!((angle_t.to_degrees() - 19.2).abs() < 0.5);
    }

    #[test]
    fn test_snell_normal_incidence() {
        let angle_t = snell(1.0, 1.5, 0.0).unwrap();
        assert!(angle_t.abs() < EPS); // normal incidence → no bending
    }

    #[test]
    fn test_snell_tir() {
        let angle_i = deg_to_rad(45.0);
        let result = snell(Medium::GLASS.n, Medium::AIR.n, angle_i);
        assert!(result.is_err());
    }

    #[test]
    fn test_critical_angle_glass_air() {
        let ca = critical_angle(Medium::GLASS.n, Medium::AIR.n).unwrap();
        // Critical angle for glass→air ≈ 41.1°
        assert!((ca.to_degrees() - 41.1).abs() < 0.5);
    }

    #[test]
    fn test_critical_angle_requires_n1_gt_n2() {
        assert!(critical_angle(1.0, 1.5).is_err());
    }

    #[test]
    fn test_reflect_angle() {
        assert!((reflect_angle(0.5) - 0.5).abs() < EPS);
    }

    #[test]
    fn test_reflect_2d() {
        // Light hitting a horizontal surface from above-left
        let dir = [0.707, -0.707]; // 45° downward
        let normal = [0.0, 1.0]; // surface pointing up
        let r = reflect_2d(dir, normal);
        assert!((r[0] - 0.707).abs() < 0.01);
        assert!((r[1] - 0.707).abs() < 0.01); // reflected upward
    }

    #[test]
    fn test_reflect_3d() {
        let dir = [0.0, -1.0, 0.0]; // straight down
        let normal = [0.0, 1.0, 0.0]; // surface pointing up
        let r = reflect_3d(dir, normal);
        assert!((r[0]).abs() < EPS);
        assert!((r[1] - 1.0).abs() < EPS); // reflected straight up
        assert!((r[2]).abs() < EPS);
    }

    #[test]
    fn test_fresnel_normal_air_glass() {
        let r = fresnel_normal(Medium::AIR.n, Medium::GLASS.n);
        // ~4% reflection at normal incidence for glass
        assert!((r - 0.04).abs() < 0.01);
    }

    #[test]
    fn test_fresnel_normal_symmetric() {
        // Same result regardless of direction
        let r1 = fresnel_normal(1.0, 1.5);
        let r2 = fresnel_normal(1.5, 1.0);
        assert!((r1 - r2).abs() < EPS);
    }

    #[test]
    fn test_fresnel_unpolarized_low_angle() {
        let r = fresnel_unpolarized(1.0, 1.5, deg_to_rad(10.0)).unwrap();
        // Near-normal: should be close to fresnel_normal
        assert!((r - fresnel_normal(1.0, 1.5)).abs() < 0.01);
    }

    #[test]
    fn test_fresnel_unpolarized_high_angle() {
        let r = fresnel_unpolarized(1.0, 1.5, deg_to_rad(80.0)).unwrap();
        // Grazing angle: high reflectance
        assert!(r > 0.3);
    }

    #[test]
    fn test_brewster_angle() {
        let ba = brewster_angle(1.0, 1.5);
        // Brewster's angle for air→glass ≈ 56.3°
        assert!((ba.to_degrees() - 56.3).abs() < 0.5);
    }

    #[test]
    fn test_beer_lambert() {
        // No absorption → no change
        assert!((beer_lambert(1.0, 0.0, 10.0) - 1.0).abs() < EPS);
        // Some absorption
        let result = beer_lambert(1.0, 0.1, 10.0);
        assert!(result < 1.0);
        assert!(result > 0.0);
        // exp(-1) ≈ 0.368
        assert!((beer_lambert(1.0, 1.0, 1.0) - (-1.0f64).exp()).abs() < EPS);
    }

    #[test]
    fn test_medium_custom_valid() {
        let m = Medium::custom(1.8, "custom").unwrap();
        assert!((m.n - 1.8).abs() < EPS);
    }

    #[test]
    fn test_medium_custom_invalid() {
        assert!(Medium::custom(0.5, "invalid").is_err());
    }

    #[test]
    fn test_deg_rad_roundtrip() {
        let deg = 45.0;
        assert!((rad_to_deg(deg_to_rad(deg)) - deg).abs() < EPS);
    }

    #[test]
    fn test_diamond_high_reflectance() {
        let r = fresnel_normal(Medium::AIR.n, Medium::DIAMOND.n);
        // Diamond reflects ~17% at normal incidence
        assert!((r - 0.17).abs() < 0.02);
    }

    #[test]
    fn test_snell_symmetry() {
        let angle_i = deg_to_rad(30.0);
        let angle_t = snell(1.0, 1.5, angle_i).unwrap();
        // Going back: snell(1.5, 1.0, angle_t) should give angle_i
        let angle_back = snell(1.5, 1.0, angle_t).unwrap();
        assert!((angle_back - angle_i).abs() < EPS);
    }
}
