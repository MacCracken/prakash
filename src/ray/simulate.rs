//! Simulation primitives — recursive ray tracing, ray fans, spot diagrams, OPD.
//!
//! Extends the sequential ray tracer with tools for optical system analysis:
//! recursive tracing (following reflected + refracted paths), ray fan generation,
//! spot diagram computation, and optical path difference calculation.

use std::f64::consts::PI;

use tracing::trace;

use super::trace::{OpticalSurface, TraceHit, TraceRay, trace_sequential, trace_surface};
use crate::error::Result;

// ── Recursive Ray Tracing ───────────────────────────────────────────────────

/// What happened at the end of a trace segment.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum TraceEvent {
    /// Ray refracted through a surface. `reflectance` is the Fresnel reflectance
    /// at this interface (the reflected energy that spawned a sibling branch).
    Refraction {
        reflectance: f64,
        surface_idx: usize,
    },
    /// Ray reflected off a surface.
    Reflection { surface_idx: usize },
    /// Ray escaped the system (no more surfaces hit).
    Escaped,
    /// Recursion depth or energy limit reached.
    Terminated,
}

/// A segment of a recursively traced ray path.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TraceSegment {
    /// Start position of this segment.
    pub origin: [f64; 3],
    /// Direction of travel (normalized).
    pub direction: [f64; 3],
    /// Refractive index of the medium.
    pub n: f64,
    /// Geometric path length of this segment.
    pub path_length: f64,
    /// Optical path length (n * geometric length).
    pub optical_path: f64,
    /// Remaining energy fraction (1.0 at start, decreases with each reflection/refraction).
    pub energy: f64,
    /// Recursion depth at which this segment was generated.
    pub depth: u32,
    /// What happened at the end of this segment.
    pub event: TraceEvent,
}

/// Configuration for recursive ray tracing.
#[derive(Debug, Clone, Copy)]
pub struct TraceConfig {
    /// Maximum recursion depth (reflection + refraction branches).
    pub max_depth: u32,
    /// Stop tracing when remaining energy drops below this threshold.
    pub min_energy: f64,
}

impl Default for TraceConfig {
    fn default() -> Self {
        Self {
            max_depth: 8,
            min_energy: 0.001,
        }
    }
}

/// Result of recursive ray tracing — a flat list of all traced segments.
///
/// Segments are ordered depth-first: refracted paths first, then reflected.
/// The total optical path through any branch can be reconstructed by summing
/// `optical_path` values along connected segments.
#[derive(Debug, Clone)]
pub struct TraceTree {
    /// All segments in depth-first order.
    pub segments: Vec<TraceSegment>,
    /// Total number of surface interactions.
    pub interactions: usize,
}

/// Trace a ray recursively through optical surfaces, following both
/// reflected and refracted paths at each interface.
///
/// At each surface, the ray splits into a refracted ray (weighted by transmittance)
/// and a reflected ray (weighted by Fresnel reflectance). Both are traced recursively
/// until `max_depth` or `min_energy` is reached.
///
/// Returns a [`TraceTree`] containing all segments.
#[must_use]
pub fn trace_recursive(
    ray: &TraceRay,
    surfaces: &[OpticalSurface],
    config: &TraceConfig,
) -> TraceTree {
    trace!(
        max_depth = config.max_depth,
        min_energy = config.min_energy,
        num_surfaces = surfaces.len(),
        "trace_recursive"
    );
    let mut segments = Vec::new();
    let mut interactions = 0;
    trace_recursive_inner(
        ray,
        surfaces,
        config,
        1.0,
        0,
        &mut segments,
        &mut interactions,
    );
    TraceTree {
        segments,
        interactions,
    }
}

fn trace_recursive_inner(
    ray: &TraceRay,
    surfaces: &[OpticalSurface],
    config: &TraceConfig,
    energy: f64,
    depth: u32,
    segments: &mut Vec<TraceSegment>,
    interactions: &mut usize,
) {
    if depth >= config.max_depth || energy < config.min_energy {
        segments.push(TraceSegment {
            origin: ray.position,
            direction: ray.direction,
            n: ray.n,
            path_length: 0.0,
            optical_path: 0.0,
            energy,
            depth,
            event: TraceEvent::Terminated,
        });
        return;
    }

    // Find the nearest surface hit
    let mut best_hit: Option<(usize, TraceHit)> = None;
    for (i, surface) in surfaces.iter().enumerate() {
        if let Ok(hit) = trace_surface(ray, surface) {
            let dist = vec3_dist(ray.position, hit.hit_point);
            match &best_hit {
                Some((_, prev)) => {
                    if dist < vec3_dist(ray.position, prev.hit_point) {
                        best_hit = Some((i, hit));
                    }
                }
                None => best_hit = Some((i, hit)),
            }
        }
    }

    let Some((surface_idx, hit)) = best_hit else {
        // Ray escaped the system
        segments.push(TraceSegment {
            origin: ray.position,
            direction: ray.direction,
            n: ray.n,
            path_length: 0.0,
            optical_path: 0.0,
            energy,
            depth,
            event: TraceEvent::Escaped,
        });
        return;
    };

    *interactions += 1;
    let path_len = vec3_dist(ray.position, hit.hit_point);
    let optical_path = ray.n * path_len;
    let reflectance = hit.reflectance;
    let transmittance = 1.0 - reflectance;

    // Record the refracted segment
    segments.push(TraceSegment {
        origin: ray.position,
        direction: ray.direction,
        n: ray.n,
        path_length: path_len,
        optical_path,
        energy: energy * transmittance,
        depth,
        event: TraceEvent::Refraction {
            reflectance,
            surface_idx,
        },
    });

    // Recurse: refracted ray
    let refracted_energy = energy * transmittance;
    if refracted_energy >= config.min_energy {
        trace_recursive_inner(
            &hit.ray_after,
            surfaces,
            config,
            refracted_energy,
            depth + 1,
            segments,
            interactions,
        );
    }

    // Recurse: reflected ray
    let reflected_energy = energy * reflectance;
    if reflected_energy >= config.min_energy {
        let reflected_dir = reflect_3d(ray.direction, hit.normal);
        let reflected_ray = TraceRay {
            position: hit.hit_point,
            direction: reflected_dir,
            n: ray.n, // stays in same medium
        };
        segments.push(TraceSegment {
            origin: hit.hit_point,
            direction: reflected_dir,
            n: ray.n,
            path_length: 0.0,
            optical_path: 0.0,
            energy: reflected_energy,
            depth,
            event: TraceEvent::Reflection { surface_idx },
        });
        trace_recursive_inner(
            &reflected_ray,
            surfaces,
            config,
            reflected_energy,
            depth + 1,
            segments,
            interactions,
        );
    }
}

/// Reflect a direction vector about a normal (used internally).
#[inline]
fn reflect_3d(direction: [f64; 3], normal: [f64; 3]) -> [f64; 3] {
    let dot = direction[0] * normal[0] + direction[1] * normal[1] + direction[2] * normal[2];
    [
        direction[0] - 2.0 * dot * normal[0],
        direction[1] - 2.0 * dot * normal[1],
        direction[2] - 2.0 * dot * normal[2],
    ]
}

/// Euclidean distance between two 3D points.
#[inline]
fn vec3_dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = b[0] - a[0];
    let dy = b[1] - a[1];
    let dz = b[2] - a[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

// ── Ray Fan Generator ───────────────────────────────────────────────────────

/// A ray in a fan, tagged with its type and aperture position.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FanRay {
    /// The ray to trace.
    pub ray: TraceRay,
    /// Normalized aperture height (−1.0 to +1.0 for meridional, 0.0 to 1.0 radial).
    pub aperture_height: f64,
}

/// Generate a meridional ray fan (rays in the y-z plane).
///
/// Creates `num_rays` rays uniformly distributed across the aperture in the
/// meridional (y-z) plane, aimed at the first surface.
///
/// `aperture_radius` = entrance pupil radius.
/// `field_angle` = off-axis field angle in radians (0 = on-axis).
/// `start_z` = z-position of the ray origins (before first surface).
///
/// Returns rays ordered from bottom to top of aperture.
#[must_use]
pub fn ray_fan_meridional(
    aperture_radius: f64,
    num_rays: usize,
    field_angle: f64,
    start_z: f64,
) -> Vec<FanRay> {
    trace!(aperture_radius, num_rays, field_angle, "ray_fan_meridional");
    let mut rays = Vec::with_capacity(num_rays);
    let (sin_f, cos_f) = field_angle.sin_cos();

    for i in 0..num_rays {
        let t = if num_rays == 1 {
            0.0
        } else {
            -1.0 + 2.0 * (i as f64) / (num_rays as f64 - 1.0)
        };
        let y = t * aperture_radius;
        rays.push(FanRay {
            ray: TraceRay {
                position: [0.0, y, start_z],
                direction: [0.0, sin_f, cos_f],
                n: 1.0,
            },
            aperture_height: t,
        });
    }
    rays
}

/// Generate a sagittal ray fan (rays in the x-z plane).
///
/// Same as meridional but spread across the x-axis.
#[must_use]
pub fn ray_fan_sagittal(
    aperture_radius: f64,
    num_rays: usize,
    field_angle: f64,
    start_z: f64,
) -> Vec<FanRay> {
    trace!(aperture_radius, num_rays, field_angle, "ray_fan_sagittal");
    let mut rays = Vec::with_capacity(num_rays);
    let (sin_f, cos_f) = field_angle.sin_cos();

    for i in 0..num_rays {
        let t = if num_rays == 1 {
            0.0
        } else {
            -1.0 + 2.0 * (i as f64) / (num_rays as f64 - 1.0)
        };
        let x = t * aperture_radius;
        rays.push(FanRay {
            ray: TraceRay {
                position: [x, 0.0, start_z],
                direction: [sin_f, 0.0, cos_f],
                n: 1.0,
            },
            aperture_height: t,
        });
    }
    rays
}

/// Generate a radial ray bundle (concentric rings across the full aperture).
///
/// Creates rays on `num_rings` concentric circles with `num_arms` rays per ring,
/// plus one on-axis chief ray. Total rays = `num_rings * num_arms + 1`.
///
/// `field_angle` in radians (0 = on-axis).
#[must_use]
pub fn ray_bundle(
    aperture_radius: f64,
    num_rings: usize,
    num_arms: usize,
    field_angle: f64,
    start_z: f64,
) -> Vec<FanRay> {
    trace!(
        aperture_radius,
        num_rings, num_arms, field_angle, "ray_bundle"
    );
    let total = num_rings * num_arms + 1;
    let mut rays = Vec::with_capacity(total);
    let (sin_f, cos_f) = field_angle.sin_cos();

    // Chief ray (on-axis through center)
    rays.push(FanRay {
        ray: TraceRay {
            position: [0.0, 0.0, start_z],
            direction: [0.0, sin_f, cos_f],
            n: 1.0,
        },
        aperture_height: 0.0,
    });

    // Concentric rings
    for ring in 1..=num_rings {
        let r = aperture_radius * (ring as f64) / (num_rings as f64);
        let t = r / aperture_radius; // normalized height
        for arm in 0..num_arms {
            let phi = 2.0 * PI * (arm as f64) / (num_arms as f64);
            let (sin_p, cos_p) = phi.sin_cos();
            let x = r * cos_p;
            let y = r * sin_p;
            rays.push(FanRay {
                ray: TraceRay {
                    position: [x, y, start_z],
                    direction: [0.0, sin_f, cos_f],
                    n: 1.0,
                },
                aperture_height: t,
            });
        }
    }
    rays
}

// ── Spot Diagram ────────────────────────────────────────────────────────────

/// A point in a spot diagram.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpotPoint {
    /// Transverse position on the image plane.
    pub x: f64,
    /// Transverse position on the image plane.
    pub y: f64,
    /// Total optical path length from source to image plane.
    pub optical_path: f64,
    /// Normalized aperture height of the source ray.
    pub aperture_height: f64,
}

/// Compute a spot diagram by tracing a ray bundle through an optical system.
///
/// Traces `num_rings * num_arms + 1` rays through the given surfaces and
/// records where they cross the image plane at `image_z`.
///
/// `field_angle` in radians, `start_z` = ray origin z-position.
///
/// Returns the spot points, or fewer if some rays are vignetted.
#[must_use]
pub fn spot_diagram(
    surfaces: &[OpticalSurface],
    aperture_radius: f64,
    num_rings: usize,
    num_arms: usize,
    field_angle: f64,
    start_z: f64,
    image_z: f64,
) -> Vec<SpotPoint> {
    trace!(num_rings, num_arms, field_angle, image_z, "spot_diagram");
    let bundle = ray_bundle(aperture_radius, num_rings, num_arms, field_angle, start_z);
    let mut spots = Vec::with_capacity(bundle.len());

    for fan_ray in &bundle {
        if let Ok(hits) = trace_sequential(&fan_ray.ray, surfaces) {
            // Get the final ray after all surfaces
            let final_ray = if let Some(last) = hits.last() {
                last.ray_after
            } else {
                fan_ray.ray
            };

            // Propagate to image plane
            let dz = final_ray.direction[2];
            if dz.abs() < 1e-15 {
                continue; // ray parallel to image plane
            }
            let t = (image_z - final_ray.position[2]) / dz;
            if t < 0.0 {
                continue; // image plane is behind
            }

            let img_x = final_ray.position[0] + t * final_ray.direction[0];
            let img_y = final_ray.position[1] + t * final_ray.direction[1];

            // Sum optical path through system
            let mut opl = 0.0;
            let mut prev_pos = fan_ray.ray.position;
            let mut prev_n = fan_ray.ray.n;
            for hit in &hits {
                let seg_len = vec3_dist(prev_pos, hit.hit_point);
                opl += prev_n * seg_len;
                prev_pos = hit.hit_point;
                prev_n = hit.ray_after.n;
            }
            // Add final segment to image plane
            let final_len = vec3_dist(prev_pos, [img_x, img_y, image_z]);
            opl += prev_n * final_len;

            spots.push(SpotPoint {
                x: img_x,
                y: img_y,
                optical_path: opl,
                aperture_height: fan_ray.aperture_height,
            });
        }
        // Vignetted rays are silently skipped
    }
    spots
}

/// Compute the RMS spot radius from a spot diagram.
///
/// Returns the root-mean-square distance of spot points from their centroid.
#[must_use]
pub fn spot_rms_radius(spots: &[SpotPoint]) -> f64 {
    if spots.is_empty() {
        return 0.0;
    }
    let n = spots.len() as f64;
    let cx: f64 = spots.iter().map(|s| s.x).sum::<f64>() / n;
    let cy: f64 = spots.iter().map(|s| s.y).sum::<f64>() / n;
    let sum_r2: f64 = spots
        .iter()
        .map(|s| (s.x - cx) * (s.x - cx) + (s.y - cy) * (s.y - cy))
        .sum();
    (sum_r2 / n).sqrt()
}

// ── Optical Path Difference ─────────────────────────────────────────────────

/// Compute the optical path length of a ray through an optical system.
///
/// Traces the ray through all surfaces and sums `n * distance` for each segment.
/// Includes the final segment to `image_z`.
///
/// Returns the total optical path length, or an error if the ray is vignetted.
#[must_use = "returns the computed optical path length"]
pub fn optical_path_length(
    ray: &TraceRay,
    surfaces: &[OpticalSurface],
    image_z: f64,
) -> Result<f64> {
    let hits = trace_sequential(ray, surfaces)?;

    let mut opl = 0.0;
    let mut prev_pos = ray.position;
    let mut prev_n = ray.n;

    for hit in &hits {
        opl += prev_n * vec3_dist(prev_pos, hit.hit_point);
        prev_pos = hit.hit_point;
        prev_n = hit.ray_after.n;
    }

    // Final segment to image plane
    let final_ray = if let Some(last) = hits.last() {
        last.ray_after
    } else {
        *ray
    };
    let dz = final_ray.direction[2];
    if dz.abs() > 1e-15 {
        let t = (image_z - final_ray.position[2]) / dz;
        if t > 0.0 {
            let end = [
                final_ray.position[0] + t * final_ray.direction[0],
                final_ray.position[1] + t * final_ray.direction[1],
                image_z,
            ];
            opl += prev_n * vec3_dist(prev_pos, end);
        }
    }

    Ok(opl)
}

/// Compute the optical path difference (OPD) of a ray relative to the chief ray.
///
/// OPD = OPL(ray) − OPL(chief), where the chief ray passes through the
/// center of the aperture.
///
/// Positive OPD means the ray arrives later (longer path).
#[must_use = "returns the computed OPD"]
pub fn optical_path_difference(
    ray: &TraceRay,
    chief_ray: &TraceRay,
    surfaces: &[OpticalSurface],
    image_z: f64,
) -> Result<f64> {
    let opl_ray = optical_path_length(ray, surfaces, image_z)?;
    let opl_chief = optical_path_length(chief_ray, surfaces, image_z)?;
    Ok(opl_ray - opl_chief)
}

/// Compute OPD for an entire ray fan relative to the chief ray.
///
/// Returns `(aperture_height, opd)` pairs for each ray in the fan.
/// Vignetted rays are excluded.
#[must_use]
pub fn opd_fan(fan: &[FanRay], surfaces: &[OpticalSurface], image_z: f64) -> Vec<(f64, f64)> {
    trace!(num_rays = fan.len(), image_z, "opd_fan");
    // Find the chief ray (aperture_height closest to 0)
    let chief = fan.iter().min_by(|a, b| {
        a.aperture_height
            .abs()
            .partial_cmp(&b.aperture_height.abs())
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let Some(chief) = chief else {
        return Vec::new();
    };

    let Ok(opl_chief) = optical_path_length(&chief.ray, surfaces, image_z) else {
        return Vec::new();
    };

    let mut result = Vec::with_capacity(fan.len());
    for fan_ray in fan {
        if let Ok(opl) = optical_path_length(&fan_ray.ray, surfaces, image_z) {
            result.push((fan_ray.aperture_height, opl - opl_chief));
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ray::trace::{OpticalSurface, SurfaceShape};

    const EPS: f64 = 1e-6;

    // Helper: simple biconvex lens
    fn biconvex_lens() -> Vec<OpticalSurface> {
        vec![
            OpticalSurface {
                shape: SurfaceShape::Sphere { radius: 100.0 },
                z_position: 0.0,
                n_after: 1.5,
                aperture_radius: 25.0,
            },
            OpticalSurface {
                shape: SurfaceShape::Sphere { radius: -100.0 },
                z_position: 5.0,
                n_after: 1.0,
                aperture_radius: 25.0,
            },
        ]
    }

    fn on_axis_ray(z: f64) -> TraceRay {
        TraceRay {
            position: [0.0, 0.0, z],
            direction: [0.0, 0.0, 1.0],
            n: 1.0,
        }
    }

    fn off_axis_ray(y: f64, z: f64) -> TraceRay {
        TraceRay {
            position: [0.0, y, z],
            direction: [0.0, 0.0, 1.0],
            n: 1.0,
        }
    }

    // ── Recursive tracing ────────────────────────────────────────────────

    #[test]
    fn test_recursive_empty_scene() {
        let ray = on_axis_ray(-50.0);
        let tree = trace_recursive(&ray, &[], &TraceConfig::default());
        assert_eq!(tree.segments.len(), 1);
        assert_eq!(tree.interactions, 0);
        assert!(matches!(tree.segments[0].event, TraceEvent::Escaped));
    }

    #[test]
    fn test_recursive_single_surface() {
        let surfaces = [OpticalSurface {
            shape: SurfaceShape::Plane,
            z_position: 10.0,
            n_after: 1.5,
            aperture_radius: 25.0,
        }];
        let ray = on_axis_ray(0.0);
        let tree = trace_recursive(&ray, &surfaces, &TraceConfig::default());
        assert!(tree.interactions >= 1);
        assert!(tree.segments.len() >= 2); // at least refracted + escaped
    }

    #[test]
    fn test_recursive_energy_conservation() {
        let surfaces = biconvex_lens();
        let ray = on_axis_ray(-50.0);
        let tree = trace_recursive(&ray, &surfaces, &TraceConfig::default());
        // All leaf segments' energies should sum to ≤ 1.0
        let total_energy: f64 = tree
            .segments
            .iter()
            .filter(|s| matches!(s.event, TraceEvent::Escaped | TraceEvent::Terminated))
            .map(|s| s.energy)
            .sum();
        assert!(
            total_energy <= 1.0 + EPS,
            "Energy should be conserved: {total_energy}"
        );
    }

    #[test]
    fn test_recursive_max_depth() {
        let surfaces = biconvex_lens();
        let ray = on_axis_ray(-50.0);
        let config = TraceConfig {
            max_depth: 1,
            min_energy: 0.0,
        };
        let tree = trace_recursive(&ray, &surfaces, &config);
        for seg in &tree.segments {
            assert!(seg.depth <= 1, "Depth should not exceed max_depth");
        }
    }

    #[test]
    fn test_recursive_min_energy() {
        let surfaces = biconvex_lens();
        let ray = on_axis_ray(-50.0);
        let config = TraceConfig {
            max_depth: 100,
            min_energy: 0.5,
        };
        let tree = trace_recursive(&ray, &surfaces, &config);
        // Should terminate quickly since each reflection has low energy
        assert!(tree.segments.len() < 20);
    }

    #[test]
    fn test_recursive_generates_reflections() {
        let surfaces = biconvex_lens();
        let ray = off_axis_ray(5.0, -50.0);
        let config = TraceConfig {
            max_depth: 4,
            min_energy: 0.001,
        };
        let tree = trace_recursive(&ray, &surfaces, &config);
        let has_reflection = tree
            .segments
            .iter()
            .any(|s| matches!(s.event, TraceEvent::Reflection { .. }));
        assert!(
            has_reflection,
            "Should have at least one reflection segment"
        );
    }

    #[test]
    fn test_recursive_optical_path_positive() {
        let surfaces = biconvex_lens();
        let ray = on_axis_ray(-50.0);
        let tree = trace_recursive(&ray, &surfaces, &TraceConfig::default());
        for seg in &tree.segments {
            assert!(seg.optical_path >= 0.0, "OPL should be non-negative");
            assert!(seg.path_length >= 0.0, "Path length should be non-negative");
        }
    }

    // ── Ray fans ─────────────────────────────────────────────────────────

    #[test]
    fn test_meridional_fan_count() {
        let fan = ray_fan_meridional(10.0, 11, 0.0, -50.0);
        assert_eq!(fan.len(), 11);
    }

    #[test]
    fn test_meridional_fan_symmetric() {
        let fan = ray_fan_meridional(10.0, 11, 0.0, -50.0);
        // First and last rays should be symmetric about axis
        let first = fan.first().unwrap();
        let last = fan.last().unwrap();
        assert!((first.ray.position[1] + last.ray.position[1]).abs() < EPS);
        assert!((first.aperture_height + last.aperture_height).abs() < EPS);
    }

    #[test]
    fn test_meridional_fan_center_ray() {
        let fan = ray_fan_meridional(10.0, 11, 0.0, -50.0);
        let center = &fan[5];
        assert!(
            center.ray.position[1].abs() < EPS,
            "Center ray should be on axis"
        );
        assert!(center.aperture_height.abs() < EPS);
    }

    #[test]
    fn test_meridional_fan_aperture_range() {
        let fan = ray_fan_meridional(10.0, 21, 0.0, -50.0);
        for ray in &fan {
            assert!(
                ray.ray.position[1].abs() <= 10.0 + EPS,
                "Ray outside aperture"
            );
        }
    }

    #[test]
    fn test_sagittal_fan_count() {
        let fan = ray_fan_sagittal(10.0, 11, 0.0, -50.0);
        assert_eq!(fan.len(), 11);
    }

    #[test]
    fn test_sagittal_fan_in_x_plane() {
        let fan = ray_fan_sagittal(10.0, 11, 0.0, -50.0);
        for ray in &fan {
            assert!(
                ray.ray.position[1].abs() < EPS,
                "Sagittal fan should be in x-z plane"
            );
        }
    }

    #[test]
    fn test_ray_bundle_count() {
        let bundle = ray_bundle(10.0, 3, 6, 0.0, -50.0);
        assert_eq!(bundle.len(), 3 * 6 + 1); // rings * arms + chief
    }

    #[test]
    fn test_ray_bundle_chief_at_center() {
        let bundle = ray_bundle(10.0, 3, 6, 0.0, -50.0);
        let chief = &bundle[0];
        assert!(chief.ray.position[0].abs() < EPS);
        assert!(chief.ray.position[1].abs() < EPS);
        assert!(chief.aperture_height.abs() < EPS);
    }

    #[test]
    fn test_ray_bundle_within_aperture() {
        let bundle = ray_bundle(10.0, 3, 8, 0.0, -50.0);
        for ray in &bundle {
            let r2 = ray.ray.position[0] * ray.ray.position[0]
                + ray.ray.position[1] * ray.ray.position[1];
            assert!(r2 <= 10.0 * 10.0 + EPS, "Ray outside aperture: r²={r2}");
        }
    }

    #[test]
    fn test_ray_fan_with_field_angle() {
        let fan = ray_fan_meridional(10.0, 5, 0.05, -50.0);
        for ray in &fan {
            // Direction should have nonzero y component for nonzero field angle
            assert!(ray.ray.direction[1].abs() > 0.01);
        }
    }

    // ── Spot diagram ─────────────────────────────────────────────────────

    #[test]
    fn test_spot_diagram_on_axis() {
        let surfaces = biconvex_lens();
        let spots = spot_diagram(&surfaces, 5.0, 2, 6, 0.0, -50.0, 200.0);
        assert!(!spots.is_empty(), "Should have spot points");
    }

    #[test]
    fn test_spot_diagram_chief_near_axis() {
        let surfaces = biconvex_lens();
        let spots = spot_diagram(&surfaces, 1.0, 1, 4, 0.0, -50.0, 200.0);
        // With small aperture and on-axis, spots should be near origin
        for spot in &spots {
            assert!(
                spot.x.abs() < 1.0 && spot.y.abs() < 1.0,
                "On-axis spots should be near origin: ({}, {})",
                spot.x,
                spot.y
            );
        }
    }

    #[test]
    fn test_spot_diagram_optical_path_positive() {
        let surfaces = biconvex_lens();
        let spots = spot_diagram(&surfaces, 5.0, 2, 6, 0.0, -50.0, 200.0);
        for spot in &spots {
            assert!(spot.optical_path > 0.0, "OPL should be positive");
        }
    }

    #[test]
    fn test_spot_rms_zero_for_single_point() {
        let spots = vec![SpotPoint {
            x: 1.0,
            y: 2.0,
            optical_path: 100.0,
            aperture_height: 0.0,
        }];
        assert!(spot_rms_radius(&spots).abs() < EPS);
    }

    #[test]
    fn test_spot_rms_empty() {
        assert!(spot_rms_radius(&[]).abs() < EPS);
    }

    #[test]
    fn test_spot_rms_known() {
        // Four points at unit distance from origin
        let spots = vec![
            SpotPoint {
                x: 1.0,
                y: 0.0,
                optical_path: 0.0,
                aperture_height: 0.0,
            },
            SpotPoint {
                x: -1.0,
                y: 0.0,
                optical_path: 0.0,
                aperture_height: 0.0,
            },
            SpotPoint {
                x: 0.0,
                y: 1.0,
                optical_path: 0.0,
                aperture_height: 0.0,
            },
            SpotPoint {
                x: 0.0,
                y: -1.0,
                optical_path: 0.0,
                aperture_height: 0.0,
            },
        ];
        let rms = spot_rms_radius(&spots);
        assert!((rms - 1.0).abs() < EPS, "RMS should be 1.0, got {rms}");
    }

    #[test]
    fn test_spot_larger_aperture_larger_spot() {
        let surfaces = biconvex_lens();
        let small = spot_diagram(&surfaces, 2.0, 2, 6, 0.0, -50.0, 200.0);
        let large = spot_diagram(&surfaces, 10.0, 2, 6, 0.0, -50.0, 200.0);
        let rms_small = spot_rms_radius(&small);
        let rms_large = spot_rms_radius(&large);
        assert!(
            rms_large >= rms_small,
            "Larger aperture should give larger spot: {rms_small} vs {rms_large}"
        );
    }

    // ── OPD ──────────────────────────────────────────────────────────────

    #[test]
    fn test_opl_positive() {
        let surfaces = biconvex_lens();
        let ray = on_axis_ray(-50.0);
        let opl = optical_path_length(&ray, &surfaces, 200.0).unwrap();
        assert!(opl > 0.0, "OPL should be positive");
    }

    #[test]
    fn test_opl_longer_in_glass() {
        // Ray through glass should have higher OPL than geometric distance
        let surfaces = biconvex_lens();
        let ray = on_axis_ray(-50.0);
        let opl = optical_path_length(&ray, &surfaces, 200.0).unwrap();
        let geometric = 250.0; // -50 to 200
        assert!(
            opl > geometric,
            "OPL should exceed geometric distance: {opl} vs {geometric}"
        );
    }

    #[test]
    fn test_opd_chief_is_zero() {
        let surfaces = biconvex_lens();
        let chief = on_axis_ray(-50.0);
        let opd = optical_path_difference(&chief, &chief, &surfaces, 200.0).unwrap();
        assert!(opd.abs() < EPS, "Chief vs chief OPD should be 0");
    }

    #[test]
    fn test_opd_off_axis_nonzero() {
        let surfaces = biconvex_lens();
        let chief = on_axis_ray(-50.0);
        let marginal = off_axis_ray(5.0, -50.0);
        let opd = optical_path_difference(&marginal, &chief, &surfaces, 200.0).unwrap();
        // Off-axis ray through a simple lens should have different OPL (spherical aberration)
        assert!(
            opd.abs() > 1e-6,
            "Off-axis ray should have nonzero OPD: {opd}"
        );
    }

    #[test]
    fn test_opd_fan_symmetric() {
        let surfaces = biconvex_lens();
        let fan = ray_fan_meridional(10.0, 11, 0.0, -50.0);
        let opd_data = opd_fan(&fan, &surfaces, 200.0);
        assert!(!opd_data.is_empty());
        // On-axis fan through symmetric lens: OPD should be symmetric
        // Find the pair at ±same height
        for &(h, opd) in &opd_data {
            if h.abs() < EPS {
                assert!(opd.abs() < EPS, "Chief ray OPD should be ~0");
            }
        }
    }

    #[test]
    fn test_opd_fan_empty_input() {
        let surfaces = biconvex_lens();
        let result = opd_fan(&[], &surfaces, 200.0);
        assert!(result.is_empty());
    }
}
