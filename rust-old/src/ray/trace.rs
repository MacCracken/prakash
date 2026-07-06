//! Sequential ray tracing through optical surfaces.

use serde::{Deserialize, Serialize};

use tracing::trace;

use super::snell_3d;
use crate::error::{PrakashError, Result};

// ── Sequential Ray Trace ──────────────────────────────────────────────────

/// An optical surface for sequential ray tracing.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
pub enum SurfaceShape {
    /// Spherical surface with given radius of curvature.
    /// Positive radius means center of curvature is to the right (along +z).
    Sphere {
        /// Radius of curvature (positive = center to the right).
        radius: f64,
    },
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
#[must_use = "returns the trace hit result"]
#[inline]
pub fn trace_surface(ray: &TraceRay, surface: &OpticalSurface) -> Result<TraceHit> {
    trace!(
        z = surface.z_position,
        n_after = surface.n_after,
        "trace_surface"
    );
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
#[must_use = "returns the trace hits"]
#[inline]
pub fn trace_sequential(
    initial_ray: &TraceRay,
    surfaces: &[OpticalSurface],
) -> Result<Vec<TraceHit>> {
    trace!(num_surfaces = surfaces.len(), "trace_sequential");
    let mut hits = Vec::with_capacity(surfaces.len());
    let mut current_ray = *initial_ray;

    for surface in surfaces {
        let hit = trace_surface(&current_ray, surface)?;
        current_ray = hit.ray_after;
        hits.push(hit);
    }

    Ok(hits)
}

// ── Polarization Ray Tracing ─────────────────────────────────────────────

/// Result of polarization-aware ray tracing through a surface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PolarizedTraceHit {
    /// The geometric trace hit (position, normal, refracted ray, reflectance).
    pub hit: TraceHit,
    /// Cumulative Fresnel reflectance ratio R_p / R_s at this surface.
    ///
    /// Values > 1 mean p-polarization reflects more, < 1 means s reflects more.
    /// At Brewster's angle, R_p ≈ 0 while R_s > 0.
    pub rp_over_rs: f64,
    /// Cumulative s-polarization transmittance through all surfaces so far.
    pub transmittance_s: f64,
    /// Cumulative p-polarization transmittance through all surfaces so far.
    pub transmittance_p: f64,
}

/// Trace a ray through a sequence of optical surfaces with polarization tracking.
///
/// At each surface, computes the s- and p-polarized Fresnel transmittances
/// separately, accumulating through the system. This enables polarization-dependent
/// analysis (e.g., finding Brewster windows, polarization-sensitive coatings).
///
/// Returns the geometric trace hits augmented with polarization data.
#[must_use = "returns the polarized trace hits"]
pub fn trace_sequential_polarized(
    initial_ray: &TraceRay,
    surfaces: &[OpticalSurface],
) -> Result<Vec<PolarizedTraceHit>> {
    trace!(num_surfaces = surfaces.len(), "trace_sequential_polarized");
    let mut results = Vec::with_capacity(surfaces.len());
    let mut current_ray = *initial_ray;
    let mut cum_ts = 1.0;
    let mut cum_tp = 1.0;

    for surface in surfaces {
        let hit = trace_surface(&current_ray, surface)?;

        // Compute s and p reflectance from the dot products
        let cos_i = {
            let d = &current_ray.direction;
            let n = &hit.normal;
            -(d[0] * n[0] + d[1] * n[1] + d[2] * n[2]).abs()
        };
        let cos_t = {
            let d = &hit.ray_after.direction;
            let n = &hit.normal;
            -(d[0] * n[0] + d[1] * n[1] + d[2] * n[2]).abs()
        };

        let n1 = current_ray.n;
        let n2 = surface.n_after;

        let rs = super::fresnel_s(n1, n2, cos_i, cos_t);
        let rp = super::fresnel_p(n1, n2, cos_i, cos_t);

        let ts = 1.0 - rs;
        let tp = 1.0 - rp;
        cum_ts *= ts;
        cum_tp *= tp;

        let rp_over_rs = if rs > 1e-15 { rp / rs } else { 0.0 };

        results.push(PolarizedTraceHit {
            hit,
            rp_over_rs,
            transmittance_s: cum_ts,
            transmittance_p: cum_tp,
        });

        current_ray = hit.ray_after;
    }

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-6;

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

    // ── Polarization ray tracing tests ───────────────────────────────────

    #[test]
    fn test_polarized_trace_normal_incidence() {
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
        let results = trace_sequential_polarized(&ray, &[surface]).unwrap();
        assert_eq!(results.len(), 1);
        // At normal incidence, s and p should be equal
        assert!(
            (results[0].transmittance_s - results[0].transmittance_p).abs() < 0.01,
            "s and p should match at normal incidence"
        );
        // Transmittance should be close to 1 - R_normal ≈ 0.96
        assert!(results[0].transmittance_s > 0.9);
    }

    #[test]
    fn test_polarized_trace_oblique() {
        // Off-axis ray — s and p should differ
        let s2 = 1.0 / 2.0f64.sqrt();
        let ray = TraceRay {
            position: [-5.0, 0.0, 0.0],
            direction: [s2, 0.0, s2], // 45° incidence
            n: 1.0,
        };
        let surface = OpticalSurface {
            shape: SurfaceShape::Plane,
            z_position: 5.0,
            n_after: 1.5,
            aperture_radius: 10.0,
        };
        let results = trace_sequential_polarized(&ray, &[surface]).unwrap();
        // At 45°, s-reflectance should be higher than p-reflectance
        assert!(
            results[0].transmittance_s < results[0].transmittance_p,
            "s should transmit less than p at 45°: Ts={}, Tp={}",
            results[0].transmittance_s,
            results[0].transmittance_p
        );
    }

    #[test]
    fn test_polarized_trace_multi_surface() {
        let ray = TraceRay {
            position: [0.0, 0.0, 0.0],
            direction: [0.0, 0.0, 1.0],
            n: 1.0,
        };
        let surfaces = [
            OpticalSurface {
                shape: SurfaceShape::Plane,
                z_position: 10.0,
                n_after: 1.5,
                aperture_radius: 5.0,
            },
            OpticalSurface {
                shape: SurfaceShape::Plane,
                z_position: 20.0,
                n_after: 1.0,
                aperture_radius: 5.0,
            },
        ];
        let results = trace_sequential_polarized(&ray, &surfaces).unwrap();
        assert_eq!(results.len(), 2);
        // After two surfaces, transmittance should be reduced
        assert!(results[1].transmittance_s < results[0].transmittance_s);
        // Cumulative: T_total = T_surface1 × T_surface2
        let t1 = results[0].transmittance_s;
        let t2_single = 1.0 - crate::ray::fresnel_normal(1.5, 1.0);
        assert!(
            (results[1].transmittance_s - t1 * t2_single).abs() < 0.01,
            "Cumulative transmittance should multiply"
        );
    }

    #[test]
    fn test_polarized_trace_transmittance_range() {
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
        let results = trace_sequential_polarized(&ray, &[surface]).unwrap();
        assert!(
            (0.0..=1.0).contains(&results[0].transmittance_s),
            "Ts out of range"
        );
        assert!(
            (0.0..=1.0).contains(&results[0].transmittance_p),
            "Tp out of range"
        );
    }
}
