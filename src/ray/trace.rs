//! Sequential ray tracing through optical surfaces.

use serde::{Deserialize, Serialize};

use super::snell_3d;
use crate::error::{PrakashError, Result};

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
#[inline]
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
}
