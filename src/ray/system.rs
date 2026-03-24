//! Optical bench — system builder, paraxial ray trace, cardinal point finder, prescriptions.
//!
//! Higher-level tools for defining and analyzing complete optical systems.
//! Builds on the sequential ray tracer with first-order (paraxial) analysis
//! and a builder API for constructing multi-element systems.

use serde::{Deserialize, Serialize};
use tracing::trace;

use super::trace::{OpticalSurface, SurfaceShape};
use crate::error::{PrakashError, Result};

// ── Paraxial Ray Trace (y-nu method) ────────────────────────────────────────

/// A paraxial ray state: height and optical direction cosine.
///
/// In the paraxial approximation, a ray is described by:
/// - `y` — ray height at a surface
/// - `nu` — n·u where n is refractive index and u is paraxial angle (slope)
///
/// This pair propagates linearly through refraction and transfer operations.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ParaxialRay {
    /// Ray height at the current surface.
    pub y: f64,
    /// Optical direction cosine: n · u (refractive index × paraxial angle).
    pub nu: f64,
}

impl ParaxialRay {
    #[must_use]
    #[inline]
    pub const fn new(y: f64, nu: f64) -> Self {
        Self { y, nu }
    }

    /// Marginal ray: enters at the edge of the aperture, parallel to axis.
    #[must_use]
    #[inline]
    pub const fn marginal(aperture_height: f64) -> Self {
        Self {
            y: aperture_height,
            nu: 0.0,
        }
    }

    /// Chief ray: enters on-axis from a field point at angle u.
    #[must_use]
    #[inline]
    pub fn chief(field_angle: f64, n: f64) -> Self {
        Self {
            y: 0.0,
            nu: n * field_angle,
        }
    }
}

/// Refract a paraxial ray at a surface with optical power φ.
///
/// nu' = nu − y·φ
///
/// `power` = φ = (n' − n) / R for a single refracting surface.
#[must_use]
#[inline]
pub fn paraxial_refract(ray: &ParaxialRay, power: f64) -> ParaxialRay {
    ParaxialRay {
        y: ray.y,
        nu: ray.nu - ray.y * power,
    }
}

/// Transfer a paraxial ray across a gap of reduced thickness t/n.
///
/// y' = y + (t/n)·nu
///
/// `reduced_thickness` = physical thickness / refractive index of the gap.
#[must_use]
#[inline]
pub fn paraxial_transfer(ray: &ParaxialRay, reduced_thickness: f64) -> ParaxialRay {
    ParaxialRay {
        y: ray.y + reduced_thickness * ray.nu,
        nu: ray.nu,
    }
}

/// Trace a paraxial ray through an optical system defined by surfaces and spacings.
///
/// `surfaces` is a list of `(power, reduced_thickness_after)` pairs where:
/// - `power` = φ = (n_after − n_before) / R
/// - `reduced_thickness_after` = distance_to_next_surface / n_after
///
/// The last surface's thickness is ignored (image space).
///
/// Returns the ray state after each surface (including the input surface).
#[must_use]
pub fn paraxial_trace(ray: &ParaxialRay, surfaces: &[(f64, f64)]) -> Vec<ParaxialRay> {
    trace!(num_surfaces = surfaces.len(), "paraxial_trace");
    let mut states = Vec::with_capacity(surfaces.len() + 1);
    let mut current = *ray;
    states.push(current);

    for (i, &(power, reduced_t)) in surfaces.iter().enumerate() {
        current = paraxial_refract(&current, power);
        states.push(current);
        if i < surfaces.len() - 1 {
            current = paraxial_transfer(&current, reduced_t);
        }
    }
    states
}

// ── Optical System Builder ──────────────────────────────────────────────────

/// An element in an optical system prescription.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PrescriptionSurface {
    /// Radius of curvature (positive = center right). `f64::INFINITY` for flat.
    pub radius: f64,
    /// Distance to the next surface along the axis.
    pub thickness: f64,
    /// Refractive index of the medium after this surface.
    pub n_after: f64,
    /// Clear aperture radius.
    pub aperture_radius: f64,
}

/// A complete optical system prescription (serializable).
///
/// Defines an ordered sequence of surfaces with their radii, thicknesses,
/// refractive indices, and aperture sizes. This is the standard way to
/// specify an optical system in lens design.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Prescription {
    /// Human-readable name for the system.
    pub name: String,
    /// Ordered list of surfaces from object side to image side.
    pub surfaces: Vec<PrescriptionSurface>,
    /// Refractive index of the medium before the first surface (typically 1.0 for air).
    pub n_initial: f64,
}

impl Prescription {
    /// Create a new empty prescription.
    #[must_use]
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            surfaces: Vec::new(),
            n_initial: 1.0,
        }
    }

    /// Set the initial medium refractive index (before the first surface).
    #[must_use]
    pub fn with_initial_medium(mut self, n: f64) -> Self {
        self.n_initial = n;
        self
    }

    /// Add a surface to the system.
    #[must_use]
    pub fn add_surface(
        mut self,
        radius: f64,
        thickness: f64,
        n_after: f64,
        aperture_radius: f64,
    ) -> Self {
        self.surfaces.push(PrescriptionSurface {
            radius,
            thickness,
            n_after,
            aperture_radius,
        });
        self
    }

    /// Convert to a list of `OpticalSurface` for ray tracing.
    ///
    /// Computes absolute z-positions from the cumulative thicknesses.
    #[must_use]
    pub fn to_trace_surfaces(&self) -> Vec<OpticalSurface> {
        let mut surfaces = Vec::with_capacity(self.surfaces.len());
        let mut z = 0.0;

        for ps in &self.surfaces {
            let shape = if ps.radius.is_infinite() {
                SurfaceShape::Plane
            } else {
                SurfaceShape::Sphere { radius: ps.radius }
            };
            surfaces.push(OpticalSurface {
                shape,
                z_position: z,
                n_after: ps.n_after,
                aperture_radius: ps.aperture_radius,
            });
            z += ps.thickness;
        }
        surfaces
    }

    /// Convert to paraxial surface parameters: `(power, reduced_thickness)` pairs.
    #[must_use]
    pub fn to_paraxial_surfaces(&self) -> Vec<(f64, f64)> {
        let mut result = Vec::with_capacity(self.surfaces.len());
        let mut n_before = self.n_initial;

        for ps in &self.surfaces {
            let power = if ps.radius.is_infinite() {
                0.0
            } else {
                (ps.n_after - n_before) / ps.radius
            };
            let reduced_t = if ps.n_after.abs() > 1e-15 {
                ps.thickness / ps.n_after
            } else {
                ps.thickness
            };
            result.push((power, reduced_t));
            n_before = ps.n_after;
        }
        result
    }

    /// Number of surfaces in the prescription.
    #[must_use]
    #[inline]
    pub fn len(&self) -> usize {
        self.surfaces.len()
    }

    /// Whether the prescription has no surfaces.
    #[must_use]
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.surfaces.is_empty()
    }
}

// ── System Cardinal Points ──────────────────────────────────────────────────

/// Cardinal points and properties of a complete optical system.
///
/// Found by tracing marginal and chief paraxial rays through the system.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SystemProperties {
    /// Effective focal length (EFL).
    pub focal_length: f64,
    /// Back focal distance (BFD) from the last surface to the rear focal point.
    pub bfd: f64,
    /// Front focal distance (FFD) from the first surface to the front focal point.
    pub ffd: f64,
    /// Total optical power (1/EFL).
    pub power: f64,
    /// Entrance pupil position (from first surface, positive = right).
    pub entrance_pupil: f64,
    /// Exit pupil position (from last surface, positive = right).
    pub exit_pupil: f64,
}

/// Find the cardinal points of an optical system from its prescription.
///
/// Traces a marginal paraxial ray (y=1, nu=0) through the system to determine
/// the effective focal length and back focal distance. Also traces a reverse
/// marginal ray for the front focal distance.
///
/// Returns `Err` if the system has zero power (afocal).
#[must_use = "returns the computed system properties"]
pub fn find_system_properties(prescription: &Prescription) -> Result<SystemProperties> {
    trace!(
        name = %prescription.name,
        num_surfaces = prescription.surfaces.len(),
        "find_system_properties"
    );

    if prescription.is_empty() {
        return Err(PrakashError::InvalidParameter {
            reason: "empty prescription has no optical power".into(),
        });
    }

    let paraxial_surfaces = prescription.to_paraxial_surfaces();

    // Forward marginal ray: y=1, nu=0 (collimated input)
    let marginal = ParaxialRay::marginal(1.0);
    let forward_states = paraxial_trace(&marginal, &paraxial_surfaces);

    // The final ray state gives us EFL and BFD
    let final_ray = forward_states.last().copied().unwrap_or(marginal);

    // EFL = -y_initial / nu_final (from a collimated input ray)
    if final_ray.nu.abs() < 1e-15 {
        return Err(PrakashError::DivisionByZero {
            context: "system is afocal (zero power)".into(),
        });
    }

    let efl = -marginal.y / final_ray.nu;
    let power = 1.0 / efl;

    // BFD = -y_final / nu_final
    let bfd = -final_ray.y / final_ray.nu;

    // For FFD, trace a reverse marginal ray
    // Reverse the system: reverse surface order, negate radii, swap n_before/n_after
    let mut reverse_surfaces = Vec::with_capacity(paraxial_surfaces.len());
    let mut n_values: Vec<f64> = vec![prescription.n_initial];
    for ps in &prescription.surfaces {
        n_values.push(ps.n_after);
    }

    for i in (0..prescription.surfaces.len()).rev() {
        let ps = &prescription.surfaces[i];
        let n_before_reverse = n_values[i + 1];
        let n_after_reverse = n_values[i];
        let power_reverse = if ps.radius.is_infinite() {
            0.0
        } else {
            // Reversed surface: radius negated, media swapped
            (n_after_reverse - n_before_reverse) / (-ps.radius)
        };
        let thickness = if i > 0 {
            prescription.surfaces[i - 1].thickness
        } else {
            0.0
        };
        let reduced_t = if n_after_reverse.abs() > 1e-15 {
            thickness / n_after_reverse
        } else {
            thickness
        };
        reverse_surfaces.push((power_reverse, reduced_t));
    }

    let reverse_marginal = ParaxialRay::marginal(1.0);
    let reverse_states = paraxial_trace(&reverse_marginal, &reverse_surfaces);
    let reverse_final = reverse_states.last().copied().unwrap_or(reverse_marginal);

    let ffd = if reverse_final.nu.abs() > 1e-15 {
        reverse_final.y / reverse_final.nu
    } else {
        f64::INFINITY
    };

    // Entrance pupil at front surface for now (aperture stop = first surface)
    let entrance_pupil = 0.0;
    // Exit pupil = BFD relative position
    let exit_pupil = bfd - efl;

    Ok(SystemProperties {
        focal_length: efl,
        bfd,
        ffd,
        power,
        entrance_pupil,
        exit_pupil,
    })
}

// ── Common Prescriptions ────────────────────────────────────────────────────

/// Prescription for a symmetric biconvex thin lens.
///
/// `focal_length` = desired focal length, `n` = glass refractive index,
/// `aperture` = clear aperture radius, `thickness` = center thickness.
#[must_use]
pub fn prescription_biconvex(
    focal_length: f64,
    n: f64,
    aperture: f64,
    thickness: f64,
) -> Prescription {
    // From lensmaker's: 1/f = (n-1)(2/R) for symmetric biconvex
    let r = 2.0 * focal_length * (n - 1.0);
    Prescription::new("biconvex")
        .add_surface(r, thickness, n, aperture)
        .add_surface(-r, 0.0, 1.0, aperture)
}

/// Prescription for a plano-convex lens (flat side first).
#[must_use]
pub fn prescription_planoconvex(
    focal_length: f64,
    n: f64,
    aperture: f64,
    thickness: f64,
) -> Prescription {
    // 1/f = (n-1)/R for plano-convex with flat first surface
    let r = -focal_length * (n - 1.0);
    Prescription::new("plano-convex")
        .add_surface(f64::INFINITY, thickness, n, aperture)
        .add_surface(r, 0.0, 1.0, aperture)
}

/// Prescription for a cemented doublet (achromatic).
///
/// `f` = desired focal length, `n1`/`n2` = crown/flint glass indices,
/// `v1`/`v2` = Abbe numbers for chromatic correction.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn prescription_doublet(
    focal_length: f64,
    n1: f64,
    v1: f64,
    n2: f64,
    v2: f64,
    aperture: f64,
    t1: f64,
    t2: f64,
) -> Prescription {
    // Achromatic condition: φ1/v1 + φ2/v2 = 0, φ1 + φ2 = 1/f
    let phi = 1.0 / focal_length;
    let phi1 = phi * v1 / (v1 - v2);
    let phi2 = phi - phi1;

    // R1 from φ1 = (n1-1)(1/R1 - 1/R2), using equi-convex for element 1
    let r1 = 2.0 * (n1 - 1.0) / phi1;
    // Cemented interface: R2 = R3 (shared radius)
    let r2 = (n1 - 1.0) / (phi1 / 2.0); // from the second surface of element 1
    // R4 from φ2 = (n2-1)(1/R3 - 1/R4)
    let r4 = if phi2.abs() > 1e-15 {
        -r2 * (n2 - 1.0) / (r2 * phi2 + (n2 - 1.0))
    } else {
        f64::INFINITY
    };

    Prescription::new("achromatic doublet")
        .add_surface(r1, t1, n1, aperture)
        .add_surface(r2, t2, n2, aperture)
        .add_surface(r4, 0.0, 1.0, aperture)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ray::trace::{TraceRay, trace_sequential};

    // ── Paraxial ray trace ───────────────────────────────────────────────

    #[test]
    fn test_paraxial_refract_zero_power() {
        let ray = ParaxialRay::new(5.0, 0.1);
        let refracted = paraxial_refract(&ray, 0.0);
        assert!((refracted.y - ray.y).abs() < 1e-15);
        assert!((refracted.nu - ray.nu).abs() < 1e-15);
    }

    #[test]
    fn test_paraxial_refract_converging() {
        let ray = ParaxialRay::marginal(10.0);
        let refracted = paraxial_refract(&ray, 0.01); // positive power
        assert!(
            refracted.nu < 0.0,
            "Converging surface should deflect ray downward"
        );
    }

    #[test]
    fn test_paraxial_transfer_free_space() {
        let ray = ParaxialRay::new(5.0, 0.1);
        let transferred = paraxial_transfer(&ray, 10.0);
        assert!((transferred.y - 6.0).abs() < 1e-15); // 5 + 10*0.1
        assert!((transferred.nu - 0.1).abs() < 1e-15); // unchanged
    }

    #[test]
    fn test_paraxial_trace_single_lens() {
        // Single thin lens: power = 1/f = 0.01 (f=100mm)
        let ray = ParaxialRay::marginal(10.0);
        let surfaces = [(0.01, 0.0)]; // power=0.01, no thickness
        let states = paraxial_trace(&ray, &surfaces);
        assert_eq!(states.len(), 2); // input + after refraction
        let final_ray = states[1];
        // nu' = 0 - 10 * 0.01 = -0.1
        assert!((final_ray.nu - (-0.1)).abs() < 1e-10);
    }

    #[test]
    fn test_paraxial_trace_two_surfaces() {
        let ray = ParaxialRay::marginal(10.0);
        let surfaces = [(0.005, 1.0), (0.005, 0.0)];
        let states = paraxial_trace(&ray, &surfaces);
        assert_eq!(states.len(), 3);
    }

    #[test]
    fn test_paraxial_trace_empty() {
        let ray = ParaxialRay::marginal(10.0);
        let states = paraxial_trace(&ray, &[]);
        assert_eq!(states.len(), 1);
        assert!((states[0].y - 10.0).abs() < 1e-15);
    }

    #[test]
    fn test_paraxial_marginal_constructor() {
        let ray = ParaxialRay::marginal(5.0);
        assert!((ray.y - 5.0).abs() < 1e-15);
        assert!(ray.nu.abs() < 1e-15);
    }

    #[test]
    fn test_paraxial_chief_constructor() {
        let ray = ParaxialRay::chief(0.1, 1.5);
        assert!(ray.y.abs() < 1e-15);
        assert!((ray.nu - 0.15).abs() < 1e-10);
    }

    // ── Prescription builder ─────────────────────────────────────────────

    #[test]
    fn test_prescription_new() {
        let p = Prescription::new("test");
        assert_eq!(p.name, "test");
        assert!(p.is_empty());
        assert_eq!(p.len(), 0);
    }

    #[test]
    fn test_prescription_add_surface() {
        let p = Prescription::new("test")
            .add_surface(100.0, 5.0, 1.5, 25.0)
            .add_surface(-100.0, 0.0, 1.0, 25.0);
        assert_eq!(p.len(), 2);
        assert!(!p.is_empty());
    }

    #[test]
    fn test_prescription_to_trace_surfaces() {
        let p = Prescription::new("test")
            .add_surface(100.0, 5.0, 1.5, 25.0)
            .add_surface(-100.0, 10.0, 1.0, 25.0);
        let surfaces = p.to_trace_surfaces();
        assert_eq!(surfaces.len(), 2);
        assert!((surfaces[0].z_position - 0.0).abs() < 1e-15);
        assert!((surfaces[1].z_position - 5.0).abs() < 1e-15);
    }

    #[test]
    fn test_prescription_flat_surface() {
        let p = Prescription::new("window")
            .add_surface(f64::INFINITY, 5.0, 1.5, 25.0)
            .add_surface(f64::INFINITY, 0.0, 1.0, 25.0);
        let surfaces = p.to_trace_surfaces();
        assert!(matches!(surfaces[0].shape, SurfaceShape::Plane));
        assert!(matches!(surfaces[1].shape, SurfaceShape::Plane));
    }

    #[test]
    fn test_prescription_to_paraxial() {
        let p = Prescription::new("test")
            .add_surface(100.0, 5.0, 1.5, 25.0)
            .add_surface(-100.0, 0.0, 1.0, 25.0);
        let paraxial = p.to_paraxial_surfaces();
        assert_eq!(paraxial.len(), 2);
        // First surface: power = (1.5-1.0)/100 = 0.005
        assert!((paraxial[0].0 - 0.005).abs() < 1e-10);
        // Second surface: power = (1.0-1.5)/(-100) = 0.005
        assert!((paraxial[1].0 - 0.005).abs() < 1e-10);
    }

    #[test]
    fn test_prescription_serde_roundtrip() {
        let p = Prescription::new("test lens")
            .add_surface(50.0, 3.0, 1.52, 12.5)
            .add_surface(-50.0, 0.0, 1.0, 12.5);
        let json = serde_json::to_string(&p).unwrap();
        let back: Prescription = serde_json::from_str(&json).unwrap();
        assert_eq!(back.name, p.name);
        assert_eq!(back.surfaces.len(), p.surfaces.len());
        assert!((back.surfaces[0].radius - 50.0).abs() < 1e-10);
    }

    #[test]
    fn test_prescription_initial_medium() {
        let p = Prescription::new("underwater").with_initial_medium(1.333);
        assert!((p.n_initial - 1.333).abs() < 1e-10);
    }

    // ── System properties ────────────────────────────────────────────────

    #[test]
    fn test_system_properties_single_thin_lens() {
        // Symmetric biconvex, R=100mm, n=1.5, thin → f≈100mm
        let p = Prescription::new("biconvex")
            .add_surface(100.0, 0.1, 1.5, 25.0)
            .add_surface(-100.0, 0.0, 1.0, 25.0);
        let props = find_system_properties(&p).unwrap();
        // EFL ≈ 100mm for thin biconvex with R=100, n=1.5
        assert!(
            (props.focal_length - 100.0).abs() < 2.0,
            "EFL should be ~100mm, got {}",
            props.focal_length
        );
    }

    #[test]
    fn test_system_properties_bfd_positive() {
        let p = Prescription::new("test")
            .add_surface(100.0, 5.0, 1.5, 25.0)
            .add_surface(-100.0, 0.0, 1.0, 25.0);
        let props = find_system_properties(&p).unwrap();
        assert!(
            props.bfd > 0.0,
            "BFD should be positive for converging lens"
        );
    }

    #[test]
    fn test_system_properties_power() {
        let p = Prescription::new("test")
            .add_surface(100.0, 5.0, 1.5, 25.0)
            .add_surface(-100.0, 0.0, 1.0, 25.0);
        let props = find_system_properties(&p).unwrap();
        assert!(
            (props.power - 1.0 / props.focal_length).abs() < 1e-10,
            "Power should be 1/EFL"
        );
    }

    #[test]
    fn test_system_properties_empty() {
        let p = Prescription::new("empty");
        assert!(find_system_properties(&p).is_err());
    }

    #[test]
    fn test_system_properties_flat_window() {
        // Two flat surfaces (window) — zero power, should error
        let p = Prescription::new("window")
            .add_surface(f64::INFINITY, 5.0, 1.5, 25.0)
            .add_surface(f64::INFINITY, 0.0, 1.0, 25.0);
        assert!(find_system_properties(&p).is_err());
    }

    #[test]
    fn test_system_properties_diverging() {
        // Biconcave lens → negative focal length
        let p = Prescription::new("biconcave")
            .add_surface(-100.0, 0.1, 1.5, 25.0)
            .add_surface(100.0, 0.0, 1.0, 25.0);
        let props = find_system_properties(&p).unwrap();
        assert!(
            props.focal_length < 0.0,
            "Biconcave should have negative EFL"
        );
    }

    // ── Common prescriptions ─────────────────────────────────────────────

    #[test]
    fn test_prescription_biconvex() {
        let p = prescription_biconvex(100.0, 1.5, 25.0, 5.0);
        assert_eq!(p.len(), 2);
        assert_eq!(p.name, "biconvex");
        // Radii should be equal and opposite
        assert!(
            (p.surfaces[0].radius + p.surfaces[1].radius).abs() < 1e-10,
            "Symmetric biconvex should have R1 = -R2"
        );
        // Check focal length is approximately correct
        let props = find_system_properties(&p).unwrap();
        assert!(
            (props.focal_length - 100.0).abs() < 5.0,
            "Focal length should be ~100mm, got {}",
            props.focal_length
        );
    }

    #[test]
    fn test_prescription_planoconvex() {
        let p = prescription_planoconvex(100.0, 1.5, 25.0, 3.0);
        assert_eq!(p.len(), 2);
        assert!(
            p.surfaces[0].radius.is_infinite(),
            "First surface should be flat"
        );
    }

    #[test]
    fn test_prescription_doublet() {
        // BK7 crown (n=1.517, V=64) + SF2 flint (n=1.648, V=34)
        let p = prescription_doublet(200.0, 1.517, 64.0, 1.648, 34.0, 25.0, 4.0, 2.0);
        assert_eq!(p.len(), 3);
        assert_eq!(p.name, "achromatic doublet");
    }

    // ── Paraxial + full trace consistency ─────────────────────────────────

    #[test]
    fn test_paraxial_matches_full_trace_on_axis() {
        // For small aperture on-axis, paraxial and full trace should agree
        let p = Prescription::new("test")
            .add_surface(100.0, 5.0, 1.5, 25.0)
            .add_surface(-100.0, 0.0, 1.0, 25.0);

        // Paraxial: marginal at y=0.1 (very paraxial)
        let paraxial = p.to_paraxial_surfaces();
        let marginal = ParaxialRay::marginal(0.1);
        let p_states = paraxial_trace(&marginal, &paraxial);
        let p_final = p_states.last().unwrap();

        // Full trace: ray at height 0.1
        let surfaces = p.to_trace_surfaces();
        let ray = TraceRay {
            position: [0.0, 0.1, -100.0],
            direction: [0.0, 0.0, 1.0],
            n: 1.0,
        };
        let hits = trace_sequential(&ray, &surfaces).unwrap();
        let final_ray = hits.last().unwrap().ray_after;

        // Compare exit angle: paraxial nu/n vs full trace direction[1]/direction[2]
        let paraxial_angle = p_final.nu; // nu = n*u, n=1 in air
        let full_angle = final_ray.direction[1] / final_ray.direction[2];

        assert!(
            (paraxial_angle - full_angle).abs() < 0.01,
            "Paraxial ({paraxial_angle:.6}) should match full trace ({full_angle:.6}) for small rays"
        );
    }
}
