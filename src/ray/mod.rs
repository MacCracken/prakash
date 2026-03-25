//! Geometric (ray) optics — reflection, refraction, Snell's law, Fresnel equations.
//!
//! All angles in radians unless suffixed with `_deg`.

use serde::{Deserialize, Serialize};
use tracing::trace;

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
    /// Vacuum (n = 1.0).
    pub const VACUUM: Medium = Medium {
        n: 1.0,
        name: "vacuum",
    };
    /// Air at STP (n = 1.000293).
    pub const AIR: Medium = Medium {
        n: 1.000293,
        name: "air",
    };
    /// Water (n = 1.333).
    pub const WATER: Medium = Medium {
        n: 1.333,
        name: "water",
    };
    /// Generic glass (n = 1.52).
    pub const GLASS: Medium = Medium {
        n: 1.52,
        name: "glass",
    };
    /// Crown glass (n = 1.523).
    pub const CROWN_GLASS: Medium = Medium {
        n: 1.523,
        name: "crown glass",
    };
    /// Flint glass (n = 1.62).
    pub const FLINT_GLASS: Medium = Medium {
        n: 1.62,
        name: "flint glass",
    };
    /// Diamond (n = 2.417).
    pub const DIAMOND: Medium = Medium {
        n: 2.417,
        name: "diamond",
    };
    /// Ice (n = 1.31).
    pub const ICE: Medium = Medium {
        n: 1.31,
        name: "ice",
    };
    /// Quartz (n = 1.544).
    pub const QUARTZ: Medium = Medium {
        n: 1.544,
        name: "quartz",
    };
    /// Sapphire (n = 1.77).
    pub const SAPPHIRE: Medium = Medium {
        n: 1.77,
        name: "sapphire",
    };
    /// Acrylic (n = 1.49).
    pub const ACRYLIC: Medium = Medium {
        n: 1.49,
        name: "acrylic",
    };
    /// Polycarbonate (n = 1.585).
    pub const POLYCARBONATE: Medium = Medium {
        n: 1.585,
        name: "polycarbonate",
    };

    /// Custom medium with a given refractive index.
    #[must_use = "returns the constructed Medium"]
    pub fn custom(n: f64, name: &'static str) -> Result<Self> {
        if n < 1.0 {
            return Err(PrakashError::InvalidRefractiveIndex { n });
        }
        Ok(Self { n, name })
    }

    /// Relative dielectric permittivity (ε_r = n²).
    ///
    /// Bridges prakash's optical refractive index to bijli's electromagnetic
    /// permittivity representation. For non-magnetic materials, n = √ε_r.
    #[must_use]
    #[inline]
    pub fn permittivity(&self) -> f64 {
        self.n * self.n
    }
}

// ── Snell's Law ─────────────────────────────────────────────────────────────

/// Apply Snell's law: n1 * sin(θ1) = n2 * sin(θ2).
///
/// Returns the refracted angle in radians, or `TotalInternalReflection` error
/// if the angle exceeds the critical angle.
///
/// Delegates to [`bijli::wave::snell_refraction_angle`] when the `bijli-backend` feature is enabled.
#[must_use = "returns the refracted angle"]
#[inline]
pub fn snell(n1: f64, n2: f64, incident_angle: f64) -> Result<f64> {
    trace!(n1, n2, incident_angle, "snell");
    snell_impl(n1, n2, incident_angle)
}

#[cfg(feature = "bijli-backend")]
#[inline]
fn snell_impl(n1: f64, n2: f64, incident_angle: f64) -> Result<f64> {
    match bijli::wave::snell_refraction_angle(n1, n2, incident_angle)? {
        Some(angle) => Ok(angle),
        None => {
            let ca = critical_angle(n1, n2).unwrap_or(0.0);
            Err(PrakashError::TotalInternalReflection {
                angle_deg: incident_angle.to_degrees(),
                critical_deg: ca.to_degrees(),
                n1,
                n2,
            })
        }
    }
}

#[cfg(not(feature = "bijli-backend"))]
#[inline]
fn snell_impl(n1: f64, n2: f64, incident_angle: f64) -> Result<f64> {
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
///
/// Delegates to [`bijli::wave::critical_angle`] when the `bijli-backend` feature is enabled.
#[must_use = "returns the critical angle"]
#[inline]
pub fn critical_angle(n1: f64, n2: f64) -> Result<f64> {
    critical_angle_impl(n1, n2)
}

#[cfg(feature = "bijli-backend")]
#[inline]
fn critical_angle_impl(n1: f64, n2: f64) -> Result<f64> {
    match bijli::wave::critical_angle(n1, n2)? {
        Some(angle) => Ok(angle),
        None => Err(PrakashError::InvalidParameter {
            reason: "critical angle requires n1 > n2".into(),
        }),
    }
}

#[cfg(not(feature = "bijli-backend"))]
#[inline]
fn critical_angle_impl(n1: f64, n2: f64) -> Result<f64> {
    if n1 <= n2 {
        return Err(PrakashError::InvalidParameter {
            reason: "critical angle requires n1 > n2".into(),
        });
    }
    Ok((n2 / n1).asin())
}

// ── Reflection ──────────────────────────────────────────────────────────────

/// Angle of reflection equals angle of incidence.
#[must_use]
#[inline]
pub fn reflect_angle(incident_angle: f64) -> f64 {
    incident_angle
}

/// Reflect a 2D direction vector about a surface normal.
///
/// Both vectors should be normalized.
#[must_use]
#[inline]
pub fn reflect_2d(direction: [f64; 2], normal: [f64; 2]) -> [f64; 2] {
    let dot = direction[0] * normal[0] + direction[1] * normal[1];
    [
        direction[0] - 2.0 * dot * normal[0],
        direction[1] - 2.0 * dot * normal[1],
    ]
}

/// Reflect a 3D direction vector about a surface normal.
#[must_use]
#[inline]
pub fn reflect_3d(direction: [f64; 3], normal: [f64; 3]) -> [f64; 3] {
    let dot = direction[0] * normal[0] + direction[1] * normal[1] + direction[2] * normal[2];
    [
        direction[0] - 2.0 * dot * normal[0],
        direction[1] - 2.0 * dot * normal[1],
        direction[2] - 2.0 * dot * normal[2],
    ]
}

/// Refract a 3D direction vector through a surface.
///
/// Given an incident direction and outward surface normal (both should be normalized),
/// computes the refracted direction using vector form of Snell's law:
///
///   t = (n1/n2)(d + (cos_θi − cos_θt)·n)
///
/// where cos_θi = −d·n and cos_θt = √(1 − (n1/n2)²(1 − cos²_θi)).
///
/// Returns `TotalInternalReflection` if the discriminant is negative.
#[must_use = "returns the refracted direction"]
#[inline]
pub fn refract_3d(direction: [f64; 3], normal: [f64; 3], n1: f64, n2: f64) -> Result<[f64; 3]> {
    let ratio = n1 / n2;

    // cos_i = -dot(d, n), ensuring the normal faces the incident ray
    let cos_i = -(direction[0] * normal[0] + direction[1] * normal[1] + direction[2] * normal[2]);

    let sin2_t = ratio * ratio * (1.0 - cos_i * cos_i);
    if sin2_t > 1.0 {
        // Compute critical angle for the error message
        let ca = if n1 > n2 {
            (n2 / n1).asin()
        } else {
            std::f64::consts::FRAC_PI_2
        };
        return Err(PrakashError::TotalInternalReflection {
            angle_deg: cos_i.acos().to_degrees(),
            critical_deg: ca.to_degrees(),
            n1,
            n2,
        });
    }

    let cos_t = (1.0 - sin2_t).sqrt();
    let factor = ratio * cos_i - cos_t;
    Ok([
        ratio * direction[0] + factor * normal[0],
        ratio * direction[1] + factor * normal[1],
        ratio * direction[2] + factor * normal[2],
    ])
}

/// 3D Snell's law: compute refracted direction with Fresnel reflectance.
///
/// Returns `(refracted_direction, reflectance)` where reflectance is the
/// unpolarized Fresnel reflectance (average of s and p).
///
/// Both `direction` and `normal` should be normalized. The normal should
/// point outward (toward the incident side).
#[must_use = "returns the refracted direction and reflectance"]
#[inline]
pub fn snell_3d(
    direction: [f64; 3],
    normal: [f64; 3],
    n1: f64,
    n2: f64,
) -> Result<([f64; 3], f64)> {
    let cos_i = -(direction[0] * normal[0] + direction[1] * normal[1] + direction[2] * normal[2]);
    let ratio = n1 / n2;
    let sin2_t = ratio * ratio * (1.0 - cos_i * cos_i);

    if sin2_t > 1.0 {
        let ca = if n1 > n2 {
            (n2 / n1).asin()
        } else {
            std::f64::consts::FRAC_PI_2
        };
        return Err(PrakashError::TotalInternalReflection {
            angle_deg: cos_i.acos().to_degrees(),
            critical_deg: ca.to_degrees(),
            n1,
            n2,
        });
    }

    let cos_t = (1.0 - sin2_t).sqrt();
    let factor = ratio * cos_i - cos_t;
    let refracted = [
        ratio * direction[0] + factor * normal[0],
        ratio * direction[1] + factor * normal[1],
        ratio * direction[2] + factor * normal[2],
    ];

    let rs = fresnel_s(n1, n2, cos_i, cos_t);
    let rp = fresnel_p(n1, n2, cos_i, cos_t);
    let reflectance = 0.5 * (rs + rp);

    Ok((refracted, reflectance))
}

mod fresnel;
pub use fresnel::*;

// ── Utility ─────────────────────────────────────────────────────────────────

/// Convert degrees to radians.
#[must_use]
#[inline]
pub fn deg_to_rad(deg: f64) -> f64 {
    deg.to_radians()
}

/// Convert radians to degrees.
#[must_use]
#[inline]
pub fn rad_to_deg(rad: f64) -> f64 {
    rad.to_degrees()
}

mod dispersion;
mod simulate;
mod system;
mod trace;

pub use dispersion::*;
pub use simulate::*;
pub use system::*;
pub use trace::*;

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    const EPS: f64 = 1e-6;

    #[test]
    fn test_medium_presets_valid() {
        let presets = [
            Medium::VACUUM,
            Medium::AIR,
            Medium::WATER,
            Medium::GLASS,
            Medium::CROWN_GLASS,
            Medium::FLINT_GLASS,
            Medium::DIAMOND,
            Medium::ICE,
            Medium::QUARTZ,
            Medium::SAPPHIRE,
            Medium::ACRYLIC,
            Medium::POLYCARBONATE,
        ];
        for m in &presets {
            assert!(m.n >= 1.0, "{} has n={} < 1.0", m.name, m.n);
        }
    }

    #[test]
    fn test_medium_custom_valid() {
        let m = Medium::custom(1.8, "custom").unwrap();
        assert!((m.n - 1.8).abs() < EPS);
        assert_eq!(m.name, "custom");
    }

    #[test]
    fn test_medium_custom_exactly_one() {
        let m = Medium::custom(1.0, "edge").unwrap();
        assert!((m.n - 1.0).abs() < EPS);
    }

    #[test]
    fn test_medium_custom_invalid() {
        assert!(Medium::custom(0.5, "invalid").is_err());
        assert!(Medium::custom(0.999, "invalid").is_err());
        assert!(Medium::custom(-1.0, "negative").is_err());
    }

    #[test]
    fn test_medium_ordering() {
        // Physical ordering: vacuum < air < water < glass < diamond
        let ns = [
            Medium::VACUUM.n,
            Medium::AIR.n,
            Medium::WATER.n,
            Medium::GLASS.n,
            Medium::DIAMOND.n,
        ];
        for w in ns.windows(2) {
            assert!(w[0] < w[1], "{} should be < {}", w[0], w[1]);
        }
    }

    #[test]
    fn test_medium_serializes() {
        let m = Medium::GLASS;
        let json = serde_json::to_string(&m).unwrap();
        assert!(json.contains("1.52"));
        assert!(json.contains("glass"));
    }

    // ── Snell's law tests ─────────────────────────────────────────────────

    #[test]
    fn test_snell_air_to_glass() {
        let angle_i = deg_to_rad(30.0);
        let angle_t = snell(Medium::AIR.n, Medium::GLASS.n, angle_i).unwrap();
        assert!((angle_t.to_degrees() - 19.2).abs() < 0.5);
    }

    #[test]
    fn test_snell_normal_incidence() {
        let angle_t = snell(1.0, 1.5, 0.0).unwrap();
        assert!(angle_t.abs() < EPS);
    }

    #[test]
    fn test_snell_tir() {
        let angle_i = deg_to_rad(45.0);
        let result = snell(Medium::GLASS.n, Medium::AIR.n, angle_i);
        assert!(result.is_err());
    }

    #[test]
    fn test_snell_at_critical_angle() {
        // At exactly the critical angle, sin_t = 1.0 (grazing)
        let ca = critical_angle(Medium::GLASS.n, Medium::AIR.n).unwrap();
        let result = snell(Medium::GLASS.n, Medium::AIR.n, ca);
        // Should succeed — refracted ray goes along surface (90°)
        assert!(result.is_ok());
        let angle_t = result.unwrap();
        assert!((angle_t - PI / 2.0).abs() < 0.01);
    }

    #[test]
    fn test_snell_symmetry() {
        let angle_i = deg_to_rad(30.0);
        let angle_t = snell(1.0, 1.5, angle_i).unwrap();
        let angle_back = snell(1.5, 1.0, angle_t).unwrap();
        assert!((angle_back - angle_i).abs() < EPS);
    }

    #[test]
    fn test_snell_bends_toward_normal_entering_denser() {
        for angle_deg in [10.0, 20.0, 30.0, 40.0, 50.0] {
            let angle_i = deg_to_rad(angle_deg);
            let angle_t = snell(1.0, 1.5, angle_i).unwrap();
            assert!(
                angle_t < angle_i,
                "Should bend toward normal at {angle_deg}°"
            );
        }
    }

    #[test]
    fn test_snell_bends_away_from_normal_leaving_denser() {
        for angle_deg in [10.0, 15.0, 20.0, 25.0] {
            let angle_i = deg_to_rad(angle_deg);
            let angle_t = snell(1.5, 1.0, angle_i).unwrap();
            assert!(angle_t > angle_i, "Should bend away at {angle_deg}°");
        }
    }

    #[test]
    fn test_snell_same_medium_no_bending() {
        let angle_i = deg_to_rad(45.0);
        let angle_t = snell(1.5, 1.5, angle_i).unwrap();
        assert!((angle_t - angle_i).abs() < EPS);
    }

    #[test]
    fn test_snell_multiple_materials() {
        let pairs = [
            (Medium::AIR, Medium::WATER),
            (Medium::AIR, Medium::DIAMOND),
            (Medium::WATER, Medium::GLASS),
        ];
        for (m1, m2) in &pairs {
            let result = snell(m1.n, m2.n, deg_to_rad(30.0));
            assert!(result.is_ok(), "Snell failed for {} → {}", m1.name, m2.name);
        }
    }

    // ── Critical angle tests ──────────────────────────────────────────────

    #[test]
    fn test_critical_angle_glass_air() {
        let ca = critical_angle(Medium::GLASS.n, Medium::AIR.n).unwrap();
        assert!((ca.to_degrees() - 41.1).abs() < 0.5);
    }

    #[test]
    fn test_critical_angle_requires_n1_gt_n2() {
        assert!(critical_angle(1.0, 1.5).is_err());
        assert!(critical_angle(1.5, 1.5).is_err());
    }

    #[test]
    fn test_critical_angle_diamond_has_small_critical_angle() {
        let ca_diamond = critical_angle(Medium::DIAMOND.n, Medium::AIR.n).unwrap();
        let ca_glass = critical_angle(Medium::GLASS.n, Medium::AIR.n).unwrap();
        assert!(
            ca_diamond < ca_glass,
            "Diamond should have smaller critical angle (more TIR)"
        );
    }

    // ── Reflection tests ──────────────────────────────────────────────────

    #[test]
    fn test_reflect_angle() {
        assert!((reflect_angle(0.5) - 0.5).abs() < EPS);
        assert!((reflect_angle(0.0) - 0.0).abs() < EPS);
        assert!((reflect_angle(PI / 4.0) - PI / 4.0).abs() < EPS);
    }

    #[test]
    fn test_reflect_2d() {
        let dir = [0.707, -0.707];
        let normal = [0.0, 1.0];
        let r = reflect_2d(dir, normal);
        assert!((r[0] - 0.707).abs() < 0.01);
        assert!((r[1] - 0.707).abs() < 0.01);
    }

    #[test]
    fn test_reflect_2d_normal_incidence() {
        let dir = [0.0, -1.0];
        let normal = [0.0, 1.0];
        let r = reflect_2d(dir, normal);
        assert!(r[0].abs() < EPS);
        assert!((r[1] - 1.0).abs() < EPS);
    }

    #[test]
    fn test_reflect_3d() {
        let dir = [0.0, -1.0, 0.0];
        let normal = [0.0, 1.0, 0.0];
        let r = reflect_3d(dir, normal);
        assert!((r[0]).abs() < EPS);
        assert!((r[1] - 1.0).abs() < EPS);
        assert!((r[2]).abs() < EPS);
    }

    #[test]
    fn test_reflect_3d_preserves_magnitude() {
        let dir = [0.577, -0.577, 0.577]; // approximately normalized
        let normal = [0.0, 1.0, 0.0];
        let r = reflect_3d(dir, normal);
        let mag_in = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
        let mag_out = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
        assert!((mag_in - mag_out).abs() < 0.01);
    }

    #[test]
    fn test_reflect_3d_diagonal() {
        // 45° angle onto XY plane
        let s = 1.0 / 2.0f64.sqrt();
        let dir = [s, -s, 0.0];
        let normal = [0.0, 1.0, 0.0];
        let r = reflect_3d(dir, normal);
        assert!((r[0] - s).abs() < EPS);
        assert!((r[1] - s).abs() < EPS);
        assert!(r[2].abs() < EPS);
    }

    // ── Fresnel tests ─────────────────────────────────────────────────────

    #[test]
    fn test_fresnel_normal_air_glass() {
        let r = fresnel_normal(Medium::AIR.n, Medium::GLASS.n);
        assert!((r - 0.04).abs() < 0.01);
    }

    #[test]
    fn test_fresnel_normal_symmetric() {
        let r1 = fresnel_normal(1.0, 1.5);
        let r2 = fresnel_normal(1.5, 1.0);
        assert!((r1 - r2).abs() < EPS);
    }

    #[test]
    fn test_fresnel_normal_same_medium() {
        assert!(fresnel_normal(1.5, 1.5).abs() < EPS);
    }

    #[test]
    fn test_fresnel_normal_range_0_to_1() {
        let pairs = [(1.0, 1.5), (1.0, 2.417), (1.333, 1.52), (1.0, 4.0)];
        for (n1, n2) in pairs {
            let r = fresnel_normal(n1, n2);
            assert!((0.0..=1.0).contains(&r), "fresnel_normal({n1}, {n2}) = {r}");
        }
    }

    #[test]
    fn test_fresnel_s_and_p_at_brewster() {
        let ba = brewster_angle(1.0, 1.5);
        let cos_i = ba.cos();
        let cos_t = snell(1.0, 1.5, ba).unwrap().cos();
        let rp = fresnel_p(1.0, 1.5, cos_i, cos_t);
        assert!(
            rp < 0.001,
            "p-polarized reflectance should be ~0 at Brewster's angle, got {rp}"
        );
        let rs = fresnel_s(1.0, 1.5, cos_i, cos_t);
        assert!(
            rs > rp,
            "s-polarized should be higher than p at Brewster's angle"
        );
    }

    #[test]
    fn test_fresnel_unpolarized_low_angle() {
        let r = fresnel_unpolarized(1.0, 1.5, deg_to_rad(10.0)).unwrap();
        assert!((r - fresnel_normal(1.0, 1.5)).abs() < 0.01);
    }

    #[test]
    fn test_fresnel_unpolarized_high_angle() {
        let r = fresnel_unpolarized(1.0, 1.5, deg_to_rad(80.0)).unwrap();
        assert!(r > 0.3);
    }

    #[test]
    fn test_fresnel_unpolarized_monotonic_increase() {
        // Reflectance should increase monotonically with angle (air→glass)
        let mut prev = 0.0;
        for deg in (0..=85).step_by(5) {
            let r = fresnel_unpolarized(1.0, 1.5, deg_to_rad(deg as f64)).unwrap();
            assert!(
                r >= prev - EPS,
                "Reflectance should increase: at {deg}° got {r} < {prev}"
            );
            prev = r;
        }
    }

    // ── Brewster & Beer-Lambert ───────────────────────────────────────────

    #[test]
    fn test_brewster_angle() {
        let ba = brewster_angle(1.0, 1.5);
        assert!((ba.to_degrees() - 56.3).abs() < 0.5);
    }

    #[test]
    fn test_brewster_complementary() {
        // Brewster angles from both sides should sum to 90°
        let b1 = brewster_angle(1.0, 1.5);
        let b2 = brewster_angle(1.5, 1.0);
        assert!((b1 + b2 - PI / 2.0).abs() < EPS);
    }

    #[test]
    fn test_beer_lambert_no_absorption() {
        assert!((beer_lambert(1.0, 0.0, 10.0) - 1.0).abs() < EPS);
    }

    #[test]
    fn test_beer_lambert_exact() {
        assert!((beer_lambert(1.0, 1.0, 1.0) - (-1.0f64).exp()).abs() < EPS);
    }

    #[test]
    fn test_beer_lambert_doubles_distance_squares_attenuation() {
        let i1 = beer_lambert(1.0, 0.5, 1.0);
        let i2 = beer_lambert(1.0, 0.5, 2.0);
        // exp(-0.5*2) = exp(-1) = exp(-0.5)^2
        assert!((i2 - i1 * i1).abs() < EPS);
    }

    #[test]
    fn test_beer_lambert_zero_distance() {
        assert!((beer_lambert(5.0, 100.0, 0.0) - 5.0).abs() < EPS);
    }

    #[test]
    fn test_beer_lambert_always_positive() {
        let i = beer_lambert(1.0, 10.0, 100.0);
        assert!(i >= 0.0);
    }

    // ── Degree/Radian ─────────────────────────────────────────────────────

    #[test]
    fn test_deg_rad_roundtrip() {
        for deg in [0.0, 30.0, 45.0, 60.0, 90.0, 180.0, 360.0] {
            assert!((rad_to_deg(deg_to_rad(deg)) - deg).abs() < EPS);
        }
    }

    #[test]
    fn test_deg_to_rad_known_values() {
        assert!((deg_to_rad(0.0)).abs() < EPS);
        assert!((deg_to_rad(90.0) - PI / 2.0).abs() < EPS);
        assert!((deg_to_rad(180.0) - PI).abs() < EPS);
        assert!((deg_to_rad(360.0) - 2.0 * PI).abs() < EPS);
    }

    // ── Cross-material tests ──────────────────────────────────────────────

    #[test]
    fn test_diamond_high_reflectance() {
        let r = fresnel_normal(Medium::AIR.n, Medium::DIAMOND.n);
        assert!((r - 0.17).abs() < 0.02);
    }

    #[test]
    fn test_higher_index_contrast_means_more_reflection() {
        let r_water = fresnel_normal(1.0, Medium::WATER.n);
        let r_glass = fresnel_normal(1.0, Medium::GLASS.n);
        let r_diamond = fresnel_normal(1.0, Medium::DIAMOND.n);
        assert!(r_water < r_glass);
        assert!(r_glass < r_diamond);
    }

    // ── 3D Refraction tests ───────────────────────────────────────────────

    fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    }

    fn len3(v: [f64; 3]) -> f64 {
        dot3(v, v).sqrt()
    }

    fn normalize3(v: [f64; 3]) -> [f64; 3] {
        let l = len3(v);
        [v[0] / l, v[1] / l, v[2] / l]
    }

    #[test]
    fn test_refract_3d_normal_incidence() {
        // Ray straight down onto horizontal surface → no bending
        let dir = [0.0, 0.0, -1.0]; // going -z
        let normal = [0.0, 0.0, 1.0]; // surface pointing +z
        let refracted = refract_3d(dir, normal, 1.0, 1.5).unwrap();
        assert!((refracted[0]).abs() < EPS);
        assert!((refracted[1]).abs() < EPS);
        assert!(refracted[2] < 0.0); // still going -z
        assert!((len3(refracted) - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_refract_3d_bends_toward_normal() {
        // 45° incidence, air → glass
        let s = 1.0 / 2.0f64.sqrt();
        let dir = [s, 0.0, -s]; // 45° in xz-plane
        let normal = [0.0, 0.0, 1.0];
        let refracted = refract_3d(dir, normal, 1.0, 1.5).unwrap();

        // Refracted angle should be less than 45°
        let cos_t = -refracted[2]; // dot with -normal = cos of refracted angle
        let angle_t = cos_t.acos();
        assert!(
            angle_t < PI / 4.0,
            "Should bend toward normal entering denser medium"
        );
    }

    #[test]
    fn test_refract_3d_tir() {
        // Glass → air at steep angle → TIR
        let dir = normalize3([0.8, 0.0, -0.6]);
        let normal = [0.0, 0.0, 1.0];
        let result = refract_3d(dir, normal, 1.5, 1.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_refract_3d_snell_law_holds() {
        // Verify n1*sin(θi) = n2*sin(θt) for the 3D refraction
        let n1 = 1.0;
        let n2 = 1.5;
        let dir = normalize3([0.5, 0.0, -0.866]); // ~30° incidence
        let normal = [0.0, 0.0, 1.0];

        let cos_i = -dot3(dir, normal);
        let sin_i = (1.0 - cos_i * cos_i).sqrt();

        let refracted = refract_3d(dir, normal, n1, n2).unwrap();
        let cos_t = -dot3(refracted, normal);
        let sin_t = (1.0 - cos_t * cos_t).sqrt();

        // n1 * sin(θi) should equal n2 * sin(θt)
        assert!(
            (n1 * sin_i - n2 * sin_t).abs() < 0.01,
            "Snell's law violated: {} != {}",
            n1 * sin_i,
            n2 * sin_t
        );
    }

    #[test]
    fn test_refract_3d_reversible() {
        let dir = normalize3([0.3, 0.2, -0.9]);
        let normal = [0.0, 0.0, 1.0];
        let refracted = refract_3d(dir, normal, 1.0, 1.5).unwrap();

        // Reverse: negate direction, negate normal, swap n1/n2
        let neg_refracted = [-refracted[0], -refracted[1], -refracted[2]];
        let neg_normal = [0.0, 0.0, -1.0];
        let back = refract_3d(neg_refracted, neg_normal, 1.5, 1.0).unwrap();

        // Should recover original direction (negated)
        let neg_dir = [-dir[0], -dir[1], -dir[2]];
        assert!((back[0] - neg_dir[0]).abs() < 0.01);
        assert!((back[1] - neg_dir[1]).abs() < 0.01);
        assert!((back[2] - neg_dir[2]).abs() < 0.01);
    }

    #[test]
    fn test_refract_3d_same_medium() {
        let dir = normalize3([0.5, 0.3, -0.8]);
        let normal = [0.0, 0.0, 1.0];
        let refracted = refract_3d(dir, normal, 1.5, 1.5).unwrap();
        // Same medium → no bending
        assert!((refracted[0] - dir[0]).abs() < 0.01);
        assert!((refracted[1] - dir[1]).abs() < 0.01);
        assert!((refracted[2] - dir[2]).abs() < 0.01);
    }

    #[test]
    fn test_refract_3d_preserves_tangent_plane() {
        // Refracted ray should stay in the plane of incidence (incident + normal).
        // Check: tangent components of incident and refracted are parallel.
        let dir = normalize3([0.4, 0.3, -0.8]);
        let normal = [0.0, 0.0, 1.0];
        let refracted = refract_3d(dir, normal, 1.0, 1.5).unwrap();

        // For normal=[0,0,1], tangent components are [x,y]
        // They should be parallel (2D cross product ≈ 0)
        let cross_2d = dir[0] * refracted[1] - dir[1] * refracted[0];
        assert!(
            cross_2d.abs() < 0.01,
            "Refracted ray should stay in plane of incidence"
        );
    }

    // ── 3D Snell's law tests ──────────────────────────────────────────────

    #[test]
    fn test_snell_3d_returns_reflectance() {
        let dir = normalize3([0.5, 0.0, -0.866]);
        let normal = [0.0, 0.0, 1.0];
        let (refracted, reflectance) = snell_3d(dir, normal, 1.0, 1.5).unwrap();
        assert!((0.0..=1.0).contains(&reflectance));
        assert!(refracted[2] < 0.0); // still propagating forward
    }

    #[test]
    fn test_snell_3d_normal_incidence_matches_fresnel_normal() {
        let dir = [0.0, 0.0, -1.0];
        let normal = [0.0, 0.0, 1.0];
        let (_, reflectance) = snell_3d(dir, normal, 1.0, 1.5).unwrap();
        let expected = fresnel_normal(1.0, 1.5);
        assert!(
            (reflectance - expected).abs() < 0.01,
            "3D Snell reflectance at normal incidence should match fresnel_normal"
        );
    }

    #[test]
    fn test_snell_3d_tir() {
        let dir = normalize3([0.8, 0.0, -0.6]);
        let normal = [0.0, 0.0, 1.0];
        assert!(snell_3d(dir, normal, 1.5, 1.0).is_err());
    }

    #[test]
    fn test_snell_3d_grazing_angle_high_reflectance() {
        let dir = normalize3([0.98, 0.0, -0.2]); // very grazing
        let normal = [0.0, 0.0, 1.0];
        let (_, reflectance) = snell_3d(dir, normal, 1.0, 1.5).unwrap();
        assert!(
            reflectance > 0.3,
            "Grazing angle should have high reflectance, got {reflectance}"
        );
    }

    // ── Complex Fresnel tests ────────────────────────────────────────────

    #[test]
    fn test_complex_medium_dielectric() {
        let m = ComplexMedium::dielectric(1.52, "glass");
        assert!((m.n - 1.52).abs() < EPS);
        assert!(m.k.abs() < EPS);
    }

    #[test]
    fn test_fresnel_normal_complex_dielectric_matches_real() {
        // For k=0, complex Fresnel should match real Fresnel
        let m = ComplexMedium::dielectric(1.52, "glass");
        let r_complex = fresnel_normal_complex(1.0, &m);
        let r_real = fresnel_normal(1.0, 1.52);
        assert!(
            (r_complex - r_real).abs() < 0.001,
            "Complex Fresnel for dielectric should match real: {r_complex} vs {r_real}"
        );
    }

    #[test]
    fn test_fresnel_normal_complex_gold() {
        // Gold at 550nm: n=0.43, k=2.46 → R ≈ 0.82
        let r = fresnel_normal_complex(1.0, &ComplexMedium::GOLD_550NM);
        assert!(
            (r - 0.82).abs() < 0.05,
            "Gold normal reflectance ≈ 0.82, got {r}"
        );
    }

    #[test]
    fn test_fresnel_normal_complex_silver() {
        // Silver at 550nm: very high reflectance
        let r = fresnel_normal_complex(1.0, &ComplexMedium::SILVER_550NM);
        assert!(r > 0.9, "Silver should be highly reflective, got {r}");
    }

    #[test]
    fn test_fresnel_normal_complex_aluminum() {
        // Aluminum at 550nm: very high reflectance
        let r = fresnel_normal_complex(1.0, &ComplexMedium::ALUMINUM_550NM);
        assert!(r > 0.9, "Aluminum should be highly reflective, got {r}");
    }

    #[test]
    fn test_fresnel_s_complex_normal_incidence_matches() {
        // At normal incidence (θ=0), s and p should equal fresnel_normal_complex
        let m = ComplexMedium::GOLD_550NM;
        let r_normal = fresnel_normal_complex(1.0, &m);
        let r_s = fresnel_s_complex(1.0, &m, 0.0);
        let r_p = fresnel_p_complex(1.0, &m, 0.0);
        assert!(
            (r_s - r_normal).abs() < 0.01,
            "s at θ=0 should match normal: {r_s} vs {r_normal}"
        );
        assert!(
            (r_p - r_normal).abs() < 0.01,
            "p at θ=0 should match normal: {r_p} vs {r_normal}"
        );
    }

    #[test]
    fn test_fresnel_complex_dielectric_matches_real_at_angle() {
        // For k=0, complex s/p should match real s/p
        let m = ComplexMedium::dielectric(1.52, "glass");
        let angle = deg_to_rad(30.0);
        let r_unpol_complex = fresnel_unpolarized_complex(1.0, &m, angle);
        let r_unpol_real = fresnel_unpolarized(1.0, 1.52, angle).unwrap();
        assert!(
            (r_unpol_complex - r_unpol_real).abs() < 0.01,
            "Complex unpolarized should match real for dielectric: {r_unpol_complex} vs {r_unpol_real}"
        );
    }

    #[test]
    fn test_fresnel_complex_grazing_high_reflectance() {
        // At grazing incidence, all materials should be highly reflective
        let m = ComplexMedium::dielectric(1.52, "glass");
        let r = fresnel_unpolarized_complex(1.0, &m, deg_to_rad(85.0));
        assert!(r > 0.5, "Grazing should have high reflectance, got {r}");
    }

    #[test]
    fn test_fresnel_complex_range_valid() {
        // Reflectance should always be in [0, 1]
        let metals = [
            ComplexMedium::GOLD_550NM,
            ComplexMedium::SILVER_550NM,
            ComplexMedium::COPPER_550NM,
            ComplexMedium::ALUMINUM_550NM,
        ];
        for m in &metals {
            for deg in (0..=85).step_by(5) {
                let angle = deg_to_rad(deg as f64);
                let r_s = fresnel_s_complex(1.0, m, angle);
                let r_p = fresnel_p_complex(1.0, m, angle);
                assert!(
                    (0.0..=1.0 + 0.001).contains(&r_s),
                    "R_s out of range for {} at {deg}°: {r_s}",
                    m.name
                );
                assert!(
                    (0.0..=1.0 + 0.001).contains(&r_p),
                    "R_p out of range for {} at {deg}°: {r_p}",
                    m.name
                );
            }
        }
    }

    #[test]
    fn test_fresnel_complex_metal_always_high() {
        // Metals should have high reflectance at all angles
        let m = ComplexMedium::ALUMINUM_550NM;
        for deg in (0..=80).step_by(10) {
            let r = fresnel_unpolarized_complex(1.0, &m, deg_to_rad(deg as f64));
            assert!(
                r > 0.8,
                "Aluminum reflectance should be >0.8, got {r} at {deg}°"
            );
        }
    }
}
