//! Geometric (ray) optics — reflection, refraction, Snell's law, Fresnel equations.
//!
//! All angles in radians unless suffixed with `_deg`.

use serde::{Deserialize, Serialize};

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
    pub const VACUUM: Medium = Medium {
        n: 1.0,
        name: "vacuum",
    };
    pub const AIR: Medium = Medium {
        n: 1.000293,
        name: "air",
    };
    pub const WATER: Medium = Medium {
        n: 1.333,
        name: "water",
    };
    pub const GLASS: Medium = Medium {
        n: 1.52,
        name: "glass",
    };
    pub const CROWN_GLASS: Medium = Medium {
        n: 1.523,
        name: "crown glass",
    };
    pub const FLINT_GLASS: Medium = Medium {
        n: 1.62,
        name: "flint glass",
    };
    pub const DIAMOND: Medium = Medium {
        n: 2.417,
        name: "diamond",
    };
    pub const ICE: Medium = Medium {
        n: 1.31,
        name: "ice",
    };
    pub const QUARTZ: Medium = Medium {
        n: 1.544,
        name: "quartz",
    };
    pub const SAPPHIRE: Medium = Medium {
        n: 1.77,
        name: "sapphire",
    };
    pub const ACRYLIC: Medium = Medium {
        n: 1.49,
        name: "acrylic",
    };
    pub const POLYCARBONATE: Medium = Medium {
        n: 1.585,
        name: "polycarbonate",
    };

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

    let refracted = refract_3d(direction, normal, n1, n2)?;

    let cos_t = (1.0 - sin2_t).sqrt();
    let rs = fresnel_s(n1, n2, cos_i, cos_t);
    let rp = fresnel_p(n1, n2, cos_i, cos_t);
    let reflectance = 0.5 * (rs + rp);

    Ok((refracted, reflectance))
}

// ── Dispersion ────────────────────────────────────────────────────────────

/// Cauchy dispersion coefficients.
///
/// Models wavelength-dependent refractive index as:
///   n(λ) = b + c/λ²
///
/// Wavelength λ is in micrometers.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CauchyCoefficients {
    pub b: f64,
    pub c: f64,
}

impl CauchyCoefficients {
    /// Refractive index at a given wavelength (micrometers).
    ///
    /// Wavelength must be positive. Returns `b` coefficient for very large wavelengths.
    #[inline]
    pub fn n_at(&self, wavelength_um: f64) -> f64 {
        debug_assert!(wavelength_um > 0.0, "wavelength must be positive");
        self.b + self.c / (wavelength_um * wavelength_um)
    }

    /// Common glass (approximate BK7).
    pub const BK7: Self = Self {
        b: 1.5046,
        c: 0.004_20,
    };

    /// Fused silica (approximate).
    pub const FUSED_SILICA: Self = Self {
        b: 1.4580,
        c: 0.003_54,
    };
}

/// Sellmeier dispersion coefficients (3-term).
///
/// Models wavelength-dependent refractive index as:
///   n²(λ) = 1 + B₁λ²/(λ²−C₁) + B₂λ²/(λ²−C₂) + B₃λ²/(λ²−C₃)
///
/// Wavelength λ is in micrometers. Cᵢ values are in μm².
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SellmeierCoefficients {
    pub b1: f64,
    pub c1: f64,
    pub b2: f64,
    pub c2: f64,
    pub b3: f64,
    pub c3: f64,
}

impl SellmeierCoefficients {
    /// Refractive index at a given wavelength (micrometers).
    ///
    /// Returns n >= 1.0 for valid wavelengths away from resonance poles.
    /// Near resonances (λ² ≈ Cᵢ), results may be unphysical.
    #[inline]
    pub fn n_at(&self, wavelength_um: f64) -> f64 {
        let l2 = wavelength_um * wavelength_um;
        // Guard denominators against resonance poles (l² ≈ Cᵢ)
        let d1 = l2 - self.c1;
        let d2 = l2 - self.c2;
        let d3 = l2 - self.c3;
        let t1 = if d1.abs() > 1e-15 { self.b1 * l2 / d1 } else { 0.0 };
        let t2 = if d2.abs() > 1e-15 { self.b2 * l2 / d2 } else { 0.0 };
        let t3 = if d3.abs() > 1e-15 { self.b3 * l2 / d3 } else { 0.0 };
        let n2 = 1.0 + t1 + t2 + t3;
        // Guard against negative n² (unphysical region between resonances)
        if n2 < 1.0 {
            return 1.0;
        }
        n2.sqrt()
    }

    /// Schott N-BK7 borosilicate crown glass.
    pub const BK7: Self = Self {
        b1: 1.039_612_12,
        c1: 0.006_000_698_67,
        b2: 0.231_792_344,
        c2: 0.020_017_914_4,
        b3: 1.010_469_45,
        c3: 103.560_653,
    };

    /// Schott N-SF11 dense flint glass.
    pub const SF11: Self = Self {
        b1: 1.737_596_26,
        c1: 0.013_188_707_0,
        b2: 0.313_747_346,
        c2: 0.062_306_814_2,
        b3: 1.898_781_01,
        c3: 155.236_29,
    };

    /// Fused silica (SiO₂).
    pub const FUSED_SILICA: Self = Self {
        b1: 0.696_166_3,
        c1: 0.004_679_148_2,
        b2: 0.407_942_6,
        c2: 0.013_512_063,
        b3: 0.897_479_4,
        c3: 97.934_002_5,
    };

    /// Sapphire (Al₂O₃, ordinary ray).
    pub const SAPPHIRE: Self = Self {
        b1: 1.431_349_3,
        c1: 0.005_279_25,
        b2: 0.650_547_13,
        c2: 0.014_218_26,
        b3: 5.341_482_2,
        c3: 325.017_83,
    };

    /// Water at 25°C (Daimon & Masumura, 2007).
    pub const WATER: Self = Self {
        b1: 0.567_019_72,
        c1: 0.005_085_50,
        b2: 0.172_629_26,
        c2: 0.018_180_00,
        b3: 0.020_624_60,
        c3: 0.026_250_00,
    };

    /// Diamond (C).
    pub const DIAMOND: Self = Self {
        b1: 4.335_8,
        c1: 0.010_6,
        b2: 0.306_0,
        c2: 0.017_5,
        b3: 0.0,
        c3: 1.0, // two-term fit
    };
}

/// Fraunhofer spectral lines used for Abbe number calculation (in micrometers).
pub const FRAUNHOFER_D: f64 = 0.587_56; // Helium d-line (yellow)
pub const FRAUNHOFER_F: f64 = 0.486_13; // Hydrogen F-line (blue)
pub const FRAUNHOFER_C: f64 = 0.656_27; // Hydrogen C-line (red)

/// Abbe number (constringence) — a measure of a material's dispersion.
///
/// V = (n_d − 1) / (n_F − n_C)
///
/// Higher V means lower dispersion. Crown glasses typically V > 50,
/// flint glasses V < 50.
#[inline]
pub fn abbe_number(sellmeier: &SellmeierCoefficients) -> f64 {
    let n_d = sellmeier.n_at(FRAUNHOFER_D);
    let n_f = sellmeier.n_at(FRAUNHOFER_F);
    let n_c = sellmeier.n_at(FRAUNHOFER_C);
    (n_d - 1.0) / (n_f - n_c)
}

// ── Sequential Ray Trace ──────────────────────────────────────────────────

/// An optical surface for sequential ray tracing.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum SurfaceShape {
    /// Spherical surface with given radius of curvature.
    /// Positive radius means center of curvature is to the right (along +z).
    Sphere { radius: f64 },
    /// Flat surface (infinite radius of curvature).
    Plane,
}

/// An interface between two media at a given position along the optical axis.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct OpticalSurface {
    /// Shape of the surface.
    pub shape: SurfaceShape,
    /// Position along the optical axis (z-coordinate of surface vertex).
    pub z_position: f64,
    /// Refractive index of the medium after this surface.
    pub n_after: f64,
    /// Clear aperture radius (limits which rays pass through).
    pub aperture_radius: f64,
}

/// A ray state during sequential tracing.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TraceRay {
    /// Current position [x, y, z].
    pub position: [f64; 3],
    /// Current direction (normalized) [dx, dy, dz].
    pub direction: [f64; 3],
    /// Current medium refractive index.
    pub n: f64,
}

/// Result of tracing a ray through one surface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TraceHit {
    /// Intersection point on the surface.
    pub hit_point: [f64; 3],
    /// Surface normal at the hit point (outward).
    pub normal: [f64; 3],
    /// Refracted ray after the surface.
    pub ray_after: TraceRay,
    /// Fresnel reflectance at this surface.
    pub reflectance: f64,
}

/// Trace a ray through a single optical surface.
///
/// The surface is centered on the optical axis (z-axis) at the given z position.
/// For spherical surfaces, the center of curvature is at z_position + radius.
#[inline]
pub fn trace_surface(ray: &TraceRay, surface: &OpticalSurface) -> Result<TraceHit> {
    // Find intersection point
    let (hit_point, normal) = match surface.shape {
        SurfaceShape::Plane => {
            // Intersection with z = z_position plane
            let dz = ray.direction[2];
            if dz.abs() < 1e-15 {
                return Err(PrakashError::InvalidParameter {
                    reason: "ray parallel to plane surface".into(),
                });
            }
            let t = (surface.z_position - ray.position[2]) / dz;
            if t < 0.0 {
                return Err(PrakashError::InvalidParameter {
                    reason: "surface is behind the ray".into(),
                });
            }
            let hit = [
                ray.position[0] + t * ray.direction[0],
                ray.position[1] + t * ray.direction[1],
                ray.position[2] + t * ray.direction[2],
            ];
            // Normal points back toward the incoming ray (against propagation)
            let n = if dz > 0.0 {
                [0.0, 0.0, -1.0]
            } else {
                [0.0, 0.0, 1.0]
            };
            (hit, n)
        }
        SurfaceShape::Sphere { radius } => {
            // Center of curvature
            let center = [0.0, 0.0, surface.z_position + radius];

            // Ray-sphere intersection: solve |P + t*D - C|² = R²
            // Direction is assumed normalized, so a = dot(D,D) = 1
            let oc = [
                ray.position[0] - center[0],
                ray.position[1] - center[1],
                ray.position[2] - center[2],
            ];
            let half_b =
                oc[0] * ray.direction[0] + oc[1] * ray.direction[1] + oc[2] * ray.direction[2];
            let c = oc[0] * oc[0] + oc[1] * oc[1] + oc[2] * oc[2] - radius * radius;
            let discriminant = half_b * half_b - c;

            if discriminant < 0.0 {
                return Err(PrakashError::InvalidParameter {
                    reason: "ray misses spherical surface".into(),
                });
            }

            let sqrt_d = discriminant.sqrt();
            let t1 = -half_b - sqrt_d;
            let t2 = -half_b + sqrt_d;

            // Pick the nearest positive intersection
            let t = if t1 > 1e-10 {
                t1
            } else if t2 > 1e-10 {
                t2
            } else {
                return Err(PrakashError::InvalidParameter {
                    reason: "surface is behind the ray".into(),
                });
            };

            let hit = [
                ray.position[0] + t * ray.direction[0],
                ray.position[1] + t * ray.direction[1],
                ray.position[2] + t * ray.direction[2],
            ];

            // Normal = (hit - center) / |hit - center|, pointing outward
            let mut n = [hit[0] - center[0], hit[1] - center[1], hit[2] - center[2]];
            let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            n[0] /= len;
            n[1] /= len;
            n[2] /= len;

            // Ensure normal faces the incoming ray
            let dot_dn =
                ray.direction[0] * n[0] + ray.direction[1] * n[1] + ray.direction[2] * n[2];
            if dot_dn > 0.0 {
                n[0] = -n[0];
                n[1] = -n[1];
                n[2] = -n[2];
            }

            (hit, n)
        }
    };

    // Check aperture
    let r2 = hit_point[0] * hit_point[0] + hit_point[1] * hit_point[1];
    if r2 > surface.aperture_radius * surface.aperture_radius {
        return Err(PrakashError::InvalidParameter {
            reason: "ray outside aperture".into(),
        });
    }

    // Refract through the surface
    let (refracted_dir, reflectance) = snell_3d(ray.direction, normal, ray.n, surface.n_after)?;

    Ok(TraceHit {
        hit_point,
        normal,
        ray_after: TraceRay {
            position: hit_point,
            direction: refracted_dir,
            n: surface.n_after,
        },
        reflectance,
    })
}

/// Trace a ray through a sequence of optical surfaces.
///
/// Returns the list of hits at each surface. The ray propagates through
/// surfaces in order, refracting at each interface.
pub fn trace_sequential(
    initial_ray: &TraceRay,
    surfaces: &[OpticalSurface],
) -> Result<Vec<TraceHit>> {
    let mut hits = Vec::with_capacity(surfaces.len());
    let mut current_ray = *initial_ray;

    for surface in surfaces {
        let hit = trace_surface(&current_ray, surface)?;
        current_ray = hit.ray_after;
        hits.push(hit);
    }

    Ok(hits)
}

// ── Prism ─────────────────────────────────────────────────────────────────

/// Minimum deviation angle for a prism.
///
/// For a prism with apex angle A and refractive index n:
///   δ_min = 2·arcsin(n·sin(A/2)) − A
///
/// `apex_angle` is in radians.
#[inline]
pub fn prism_deviation(apex_angle: f64, n: f64) -> Result<f64> {
    let sin_half_a = (apex_angle / 2.0).sin();
    let arg = n * sin_half_a;
    if arg.abs() > 1.0 {
        return Err(PrakashError::TotalInternalReflection {
            angle_deg: (apex_angle / 2.0).to_degrees(),
            critical_deg: (1.0 / n).asin().to_degrees(),
            n1: n,
            n2: 1.0,
        });
    }
    Ok(2.0 * arg.asin() - apex_angle)
}

/// Angular dispersion of a prism for a given wavelength.
///
/// Returns the deviation angle using the Sellmeier equation to determine
/// the refractive index at the specified wavelength.
///
/// `apex_angle` is in radians, `wavelength_um` is in micrometers.
#[inline]
pub fn prism_dispersion(
    apex_angle: f64,
    sellmeier: &SellmeierCoefficients,
    wavelength_um: f64,
) -> Result<f64> {
    let n = sellmeier.n_at(wavelength_um);
    prism_deviation(apex_angle, n)
}

/// Angular spread of a prism across a wavelength range.
///
/// Returns (deviation_short, deviation_long, angular_spread) in radians.
/// The angular spread is the difference in deviation between the shortest
/// and longest wavelengths.
pub fn prism_angular_spread(
    apex_angle: f64,
    sellmeier: &SellmeierCoefficients,
    wavelength_short_um: f64,
    wavelength_long_um: f64,
) -> Result<(f64, f64, f64)> {
    let dev_short = prism_dispersion(apex_angle, sellmeier, wavelength_short_um)?;
    let dev_long = prism_dispersion(apex_angle, sellmeier, wavelength_long_um)?;
    Ok((dev_short, dev_long, (dev_short - dev_long).abs()))
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
    deg.to_radians()
}

/// Convert radians to degrees.
#[inline]
pub fn rad_to_deg(rad: f64) -> f64 {
    rad.to_degrees()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    const EPS: f64 = 1e-6;

    // ── Medium tests ──────────────────────────────────────────────────────

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

    // ── Cauchy dispersion tests ───────────────────────────────────────────

    #[test]
    fn test_cauchy_bk7_visible_range() {
        let n_blue = CauchyCoefficients::BK7.n_at(0.486);
        let n_red = CauchyCoefficients::BK7.n_at(0.656);
        let n_d = CauchyCoefficients::BK7.n_at(0.5876);

        assert!(n_blue > n_d, "Blue should have higher n than yellow");
        assert!(n_d > n_red, "Yellow should have higher n than red");
        assert!((n_d - 1.517).abs() < 0.02, "BK7 n_d ≈ 1.517");
    }

    #[test]
    fn test_cauchy_normal_dispersion() {
        // Normal dispersion: dn/dλ < 0 (n decreases with wavelength)
        let c = CauchyCoefficients::BK7;
        let n1 = c.n_at(0.4);
        let n2 = c.n_at(0.5);
        let n3 = c.n_at(0.6);
        let n4 = c.n_at(0.7);
        assert!(n1 > n2);
        assert!(n2 > n3);
        assert!(n3 > n4);
    }

    #[test]
    fn test_cauchy_fused_silica() {
        let n = CauchyCoefficients::FUSED_SILICA.n_at(0.5876);
        assert!((n - 1.458).abs() < 0.02);
    }

    // ── Sellmeier dispersion tests ────────────────────────────────────────

    #[test]
    fn test_sellmeier_bk7_nd() {
        let n = SellmeierCoefficients::BK7.n_at(FRAUNHOFER_D);
        assert!(
            (n - 1.5168).abs() < 0.001,
            "BK7 n_d should be ≈1.5168, got {n}"
        );
    }

    #[test]
    fn test_sellmeier_bk7_normal_dispersion() {
        let s = SellmeierCoefficients::BK7;
        let n_f = s.n_at(FRAUNHOFER_F);
        let n_d = s.n_at(FRAUNHOFER_D);
        let n_c = s.n_at(FRAUNHOFER_C);
        assert!(n_f > n_d, "n_F > n_d (blue > yellow)");
        assert!(n_d > n_c, "n_d > n_C (yellow > red)");
    }

    #[test]
    fn test_sellmeier_sf11() {
        let n = SellmeierCoefficients::SF11.n_at(FRAUNHOFER_D);
        assert!(
            (n - 1.7847).abs() < 0.001,
            "SF11 n_d should be ≈1.7847, got {n}"
        );
    }

    #[test]
    fn test_sellmeier_fused_silica() {
        let n = SellmeierCoefficients::FUSED_SILICA.n_at(FRAUNHOFER_D);
        assert!(
            (n - 1.4585).abs() < 0.001,
            "Fused silica n_d should be ≈1.4585, got {n}"
        );
    }

    #[test]
    fn test_sellmeier_water() {
        let n = SellmeierCoefficients::WATER.n_at(FRAUNHOFER_D);
        assert!(
            (n - 1.333).abs() < 0.005,
            "Water n_d should be ≈1.333, got {n}"
        );
    }

    #[test]
    fn test_sellmeier_diamond() {
        let n = SellmeierCoefficients::DIAMOND.n_at(FRAUNHOFER_D);
        assert!(
            (n - 2.417).abs() < 0.02,
            "Diamond n_d should be ≈2.417, got {n}"
        );
    }

    #[test]
    fn test_sellmeier_all_presets_reasonable() {
        let presets = [
            ("BK7", SellmeierCoefficients::BK7),
            ("SF11", SellmeierCoefficients::SF11),
            ("Fused Silica", SellmeierCoefficients::FUSED_SILICA),
            ("Sapphire", SellmeierCoefficients::SAPPHIRE),
            ("Water", SellmeierCoefficients::WATER),
            ("Diamond", SellmeierCoefficients::DIAMOND),
        ];
        for (name, s) in &presets {
            let n = s.n_at(FRAUNHOFER_D);
            assert!(
                (1.0..3.0).contains(&n),
                "{name}: n_d = {n} out of reasonable range"
            );
        }
    }

    #[test]
    fn test_sellmeier_serde_roundtrip() {
        let s = SellmeierCoefficients::BK7;
        let json = serde_json::to_string(&s).unwrap();
        let back: SellmeierCoefficients = serde_json::from_str(&json).unwrap();
        assert!((back.b1 - s.b1).abs() < EPS);
        assert!((back.c1 - s.c1).abs() < EPS);
    }

    // ── Abbe number tests ─────────────────────────────────────────────────

    #[test]
    fn test_abbe_number_bk7() {
        let v = abbe_number(&SellmeierCoefficients::BK7);
        // BK7 Abbe number ≈ 64.17
        assert!(
            (v - 64.17).abs() < 1.0,
            "BK7 Abbe number should be ≈64.17, got {v}"
        );
    }

    #[test]
    fn test_abbe_number_sf11() {
        let v = abbe_number(&SellmeierCoefficients::SF11);
        // SF11 Abbe number ≈ 25.76 (high dispersion flint)
        assert!(
            (v - 25.76).abs() < 1.0,
            "SF11 Abbe number should be ≈25.76, got {v}"
        );
    }

    #[test]
    fn test_abbe_crown_vs_flint() {
        let v_bk7 = abbe_number(&SellmeierCoefficients::BK7);
        let v_sf11 = abbe_number(&SellmeierCoefficients::SF11);
        assert!(
            v_bk7 > v_sf11,
            "Crown glass should have higher Abbe number than flint"
        );
        assert!(v_bk7 > 50.0, "Crown glass V > 50");
        assert!(v_sf11 < 50.0, "Flint glass V < 50");
    }

    #[test]
    fn test_sellmeier_resonance_guard() {
        // At a resonance pole (l² = c1), n_at should not return NaN or panic
        let s = SellmeierCoefficients::BK7;
        let resonance_wl = s.c1.sqrt(); // wavelength where l² = c1
        let n = s.n_at(resonance_wl);
        assert!(n.is_finite(), "n_at resonance should be finite, got {n}");
        assert!(n >= 1.0, "n_at resonance should be >= 1.0, got {n}");
    }

    #[test]
    fn test_refract_3d_output_normalized() {
        // Verify refracted direction has unit length
        let dir = normalize3([0.4, 0.3, -0.8]);
        let normal = [0.0, 0.0, 1.0];
        let refracted = refract_3d(dir, normal, 1.0, 1.5).unwrap();
        let len = len3(refracted);
        assert!(
            (len - 1.0).abs() < 0.01,
            "Refracted direction should be normalized, got length {len}"
        );
    }

    #[test]
    fn test_snell_3d_consistency_with_scalar() {
        // snell_3d at a known angle should agree with scalar snell()
        let angle_i = deg_to_rad(30.0);
        let n1 = 1.0;
        let n2 = 1.5;

        // Scalar
        let angle_t_scalar = snell(n1, n2, angle_i).unwrap();

        // 3D: construct direction at 30° in xz plane
        let dir = [angle_i.sin(), 0.0, -angle_i.cos()];
        let normal = [0.0, 0.0, 1.0];
        let (refracted, _) = snell_3d(dir, normal, n1, n2).unwrap();
        let cos_t_3d = -refracted[2]; // dot with -normal
        let angle_t_3d = cos_t_3d.acos();

        assert!(
            (angle_t_3d - angle_t_scalar).abs() < 0.001,
            "3D ({:.4}°) should match scalar ({:.4}°)",
            angle_t_3d.to_degrees(),
            angle_t_scalar.to_degrees()
        );
    }

    #[test]
    fn test_trace_reflectance_at_each_surface() {
        let surfaces = [
            OpticalSurface {
                shape: SurfaceShape::Plane,
                z_position: 10.0,
                n_after: 1.5,
                aperture_radius: 25.0,
            },
            OpticalSurface {
                shape: SurfaceShape::Plane,
                z_position: 15.0,
                n_after: 1.0,
                aperture_radius: 25.0,
            },
        ];
        let ray = TraceRay {
            position: [0.0, 0.0, 0.0],
            direction: [0.0, 0.0, 1.0],
            n: 1.0,
        };
        let hits = trace_sequential(&ray, &surfaces).unwrap();
        for (i, hit) in hits.iter().enumerate() {
            assert!(
                (0.0..=1.0).contains(&hit.reflectance),
                "Surface {i}: reflectance {:.4} out of range",
                hit.reflectance
            );
        }
    }

    #[test]
    fn test_prism_zero_apex() {
        // Zero apex angle → zero deviation
        let dev = prism_deviation(0.0, 1.5).unwrap();
        assert!(dev.abs() < EPS, "Zero apex should give zero deviation");
    }

    // ── Sequential ray trace tests ────────────────────────────────────────

    #[test]
    fn test_trace_plane_surface() {
        let ray = TraceRay {
            position: [0.0, 0.0, 0.0],
            direction: [0.0, 0.0, 1.0],
            n: 1.0,
        };
        let surface = OpticalSurface {
            shape: SurfaceShape::Plane,
            z_position: 10.0,
            n_after: 1.5,
            aperture_radius: 5.0,
        };
        let hit = trace_surface(&ray, &surface).unwrap();
        assert!((hit.hit_point[2] - 10.0).abs() < EPS);
        assert!((hit.ray_after.direction[2] - 1.0).abs() < 0.01); // normal incidence, no bending
    }

    #[test]
    fn test_trace_sphere_surface() {
        let ray = TraceRay {
            position: [0.0, 0.0, 0.0],
            direction: [0.0, 0.0, 1.0],
            n: 1.0,
        };
        let surface = OpticalSurface {
            shape: SurfaceShape::Sphere { radius: 50.0 },
            z_position: 10.0,
            n_after: 1.5,
            aperture_radius: 25.0,
        };
        let hit = trace_surface(&ray, &surface).unwrap();
        // On-axis ray should hit near z=10 and continue on-axis
        assert!((hit.hit_point[0]).abs() < EPS);
        assert!((hit.hit_point[1]).abs() < EPS);
        assert!(hit.hit_point[2] > 9.0 && hit.hit_point[2] < 11.0);
    }

    #[test]
    fn test_trace_aperture_blocks_ray() {
        let ray = TraceRay {
            position: [10.0, 0.0, 0.0], // far off-axis
            direction: [0.0, 0.0, 1.0],
            n: 1.0,
        };
        let surface = OpticalSurface {
            shape: SurfaceShape::Plane,
            z_position: 10.0,
            n_after: 1.5,
            aperture_radius: 5.0,
        };
        assert!(trace_surface(&ray, &surface).is_err());
    }

    #[test]
    fn test_trace_sequential_biconvex_lens() {
        // Simple biconvex lens: two spherical surfaces
        let surfaces = [
            OpticalSurface {
                shape: SurfaceShape::Sphere { radius: 100.0 },
                z_position: 0.0,
                n_after: 1.5,
                aperture_radius: 25.0,
            },
            OpticalSurface {
                shape: SurfaceShape::Sphere { radius: -100.0 },
                z_position: 5.0, // 5mm thick lens
                n_after: 1.0,
                aperture_radius: 25.0,
            },
        ];

        let ray = TraceRay {
            position: [0.0, 0.0, -50.0],
            direction: [0.0, 0.0, 1.0],
            n: 1.0,
        };

        let hits = trace_sequential(&ray, &surfaces).unwrap();
        assert_eq!(hits.len(), 2);
        // On-axis ray through biconvex should continue on-axis
        assert!((hits[1].ray_after.direction[0]).abs() < 0.01);
    }

    #[test]
    fn test_trace_off_axis_converges() {
        // Off-axis ray through converging lens should bend toward axis
        let surfaces = [
            OpticalSurface {
                shape: SurfaceShape::Sphere { radius: 50.0 },
                z_position: 0.0,
                n_after: 1.5,
                aperture_radius: 25.0,
            },
            OpticalSurface {
                shape: SurfaceShape::Sphere { radius: -50.0 },
                z_position: 5.0,
                n_after: 1.0,
                aperture_radius: 25.0,
            },
        ];

        let ray = TraceRay {
            position: [5.0, 0.0, -50.0],
            direction: [0.0, 0.0, 1.0], // parallel to axis but offset
            n: 1.0,
        };

        let hits = trace_sequential(&ray, &surfaces).unwrap();
        // After the lens, ray should be bending toward axis (negative x-direction)
        assert!(
            hits[1].ray_after.direction[0] < 0.0,
            "Off-axis ray through converging lens should bend toward axis"
        );
    }

    #[test]
    fn test_trace_sequential_empty() {
        let ray = TraceRay {
            position: [0.0, 0.0, 0.0],
            direction: [0.0, 0.0, 1.0],
            n: 1.0,
        };
        let hits = trace_sequential(&ray, &[]).unwrap();
        assert!(hits.is_empty());
    }

    // ── Prism tests ───────────────────────────────────────────────────────

    #[test]
    fn test_prism_deviation_60_degree() {
        // 60° equilateral prism with BK7 glass
        let apex = deg_to_rad(60.0);
        let n = SellmeierCoefficients::BK7.n_at(FRAUNHOFER_D);
        let dev = prism_deviation(apex, n).unwrap();
        // Minimum deviation for BK7 60° prism ≈ 40°
        let dev_deg = dev.to_degrees();
        assert!(
            dev_deg > 35.0 && dev_deg < 50.0,
            "60° BK7 prism deviation ≈ 40°, got {dev_deg}"
        );
    }

    #[test]
    fn test_prism_deviation_increases_with_n() {
        let apex = deg_to_rad(60.0);
        let dev_low = prism_deviation(apex, 1.4).unwrap();
        let dev_high = prism_deviation(apex, 1.7).unwrap();
        assert!(dev_high > dev_low, "Higher n should give larger deviation");
    }

    #[test]
    fn test_prism_deviation_tir() {
        // Very high n with large apex → TIR possible
        let apex = deg_to_rad(120.0);
        let result = prism_deviation(apex, 2.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_prism_dispersion_blue_more_than_red() {
        let apex = deg_to_rad(60.0);
        let s = SellmeierCoefficients::BK7;
        let dev_blue = prism_dispersion(apex, &s, FRAUNHOFER_F).unwrap();
        let dev_red = prism_dispersion(apex, &s, FRAUNHOFER_C).unwrap();
        assert!(
            dev_blue > dev_red,
            "Blue light should be deviated more than red"
        );
    }

    #[test]
    fn test_prism_angular_spread() {
        let apex = deg_to_rad(60.0);
        let s = SellmeierCoefficients::BK7;
        let (dev_short, dev_long, spread) =
            prism_angular_spread(apex, &s, FRAUNHOFER_F, FRAUNHOFER_C).unwrap();
        assert!(dev_short > dev_long);
        assert!(spread > 0.0);
        assert!(
            spread.to_degrees() < 5.0,
            "BK7 angular spread should be modest"
        );
    }

    #[test]
    fn test_prism_flint_more_dispersion_than_crown() {
        let apex = deg_to_rad(60.0);
        let (_, _, spread_crown) = prism_angular_spread(
            apex,
            &SellmeierCoefficients::BK7,
            FRAUNHOFER_F,
            FRAUNHOFER_C,
        )
        .unwrap();
        let (_, _, spread_flint) = prism_angular_spread(
            apex,
            &SellmeierCoefficients::SF11,
            FRAUNHOFER_F,
            FRAUNHOFER_C,
        )
        .unwrap();
        assert!(
            spread_flint > spread_crown,
            "Flint glass should have more dispersion than crown"
        );
    }
}
