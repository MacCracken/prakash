//! Lens and mirror geometry — focal length, magnification, thin lens equation, aberrations.

use serde::{Deserialize, Serialize};

use crate::error::{PrakashError, Result};

/// Thin lens equation: 1/f = 1/do + 1/di
///
/// Given focal length and object distance, returns image distance.
/// Positive di = real image (opposite side), negative = virtual image (same side).
#[inline]
pub fn thin_lens_image_distance(focal_length: f64, object_distance: f64) -> Result<f64> {
    let denom = object_distance - focal_length;
    if denom.abs() < 1e-15 {
        return Err(PrakashError::DivisionByZero {
            context: "object at focal point".into(),
        });
    }
    Ok(focal_length * object_distance / denom)
}

/// Magnification: M = -di/do
///
/// Positive M = upright image, negative = inverted.
/// |M| > 1 = magnified, |M| < 1 = diminished.
#[inline]
pub fn magnification(object_distance: f64, image_distance: f64) -> f64 {
    -image_distance / object_distance
}

/// Lensmaker's equation: focal length from radii of curvature and refractive index.
///
/// 1/f = (n-1) · (1/R1 - 1/R2)
///
/// Convention: R > 0 if center of curvature is to the right.
#[inline]
pub fn lensmaker_focal_length(n: f64, r1: f64, r2: f64) -> Result<f64> {
    let power = (n - 1.0) * (1.0 / r1 - 1.0 / r2);
    if power.abs() < 1e-15 {
        return Err(PrakashError::DivisionByZero {
            context: "flat surfaces produce no focusing".into(),
        });
    }
    Ok(1.0 / power)
}

/// Optical power in diopters (1/f in meters).
#[inline]
pub fn optical_power(focal_length_m: f64) -> Result<f64> {
    if focal_length_m.abs() < 1e-15 {
        return Err(PrakashError::DivisionByZero {
            context: "optical power requires non-zero focal length".into(),
        });
    }
    Ok(1.0 / focal_length_m)
}

/// Mirror focal length from radius of curvature: f = R/2.
#[inline]
pub fn mirror_focal_length(radius: f64) -> f64 {
    radius / 2.0
}

/// Mirror equation (same form as thin lens): 1/f = 1/do + 1/di.
#[inline]
pub fn mirror_image_distance(focal_length: f64, object_distance: f64) -> Result<f64> {
    thin_lens_image_distance(focal_length, object_distance)
}

/// Lens type classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum LensType {
    /// Converging (positive focal length): biconvex, plano-convex, positive meniscus.
    Converging,
    /// Diverging (negative focal length): biconcave, plano-concave, negative meniscus.
    Diverging,
}

/// Classify a lens by its focal length.
#[inline]
pub fn classify_lens(focal_length: f64) -> LensType {
    if focal_length > 0.0 {
        LensType::Converging
    } else {
        LensType::Diverging
    }
}

/// Two thin lenses in contact: combined focal length.
/// 1/f = 1/f1 + 1/f2
#[inline]
pub fn combined_focal_length(f1: f64, f2: f64) -> Result<f64> {
    let sum = f1 + f2;
    if sum.abs() < 1e-15 {
        return Err(PrakashError::DivisionByZero {
            context: "equal and opposite focal lengths produce zero combined power".into(),
        });
    }
    Ok((f1 * f2) / sum)
}

/// Depth of field approximation for a thin lens.
///
/// Returns (near_limit, far_limit) distances.
/// `f` = focal length, `N` = f-number, `c` = circle of confusion,
/// `d` = subject distance.
pub fn depth_of_field(focal_length: f64, f_number: f64, coc: f64, subject_dist: f64) -> (f64, f64) {
    let h = (focal_length * focal_length) / (f_number * coc); // hyperfocal distance
    let near = (h * subject_dist) / (h + (subject_dist - focal_length));
    let far_denom = h - (subject_dist - focal_length);
    let far = if far_denom <= 0.0 {
        f64::INFINITY // at or beyond hyperfocal distance
    } else {
        (h * subject_dist) / far_denom
    };
    (near, far)
}

// ── Thick Lens ────────────────────────────────────────────────────────────

/// Thick lens focal length.
///
/// 1/f = (n−1) · [1/R₁ − 1/R₂ + (n−1)·d / (n·R₁·R₂)]
///
/// `n` = refractive index, `r1`/`r2` = radii of curvature, `d` = center thickness.
/// Sign convention: R > 0 if center of curvature is to the right.
#[inline]
pub fn thick_lens_focal_length(n: f64, r1: f64, r2: f64, thickness: f64) -> Result<f64> {
    let nm1 = n - 1.0;
    let power = nm1 * (1.0 / r1 - 1.0 / r2 + nm1 * thickness / (n * r1 * r2));
    if power.abs() < 1e-15 {
        return Err(PrakashError::DivisionByZero {
            context: "thick lens has zero optical power".into(),
        });
    }
    Ok(1.0 / power)
}

/// Cardinal points of a thick lens.
///
/// Contains all the key reference points needed to analyze a thick lens system.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CardinalPoints {
    /// Effective focal length.
    pub focal_length: f64,
    /// Front focal distance (from front surface vertex to front focal point).
    pub ffd: f64,
    /// Back focal distance (from rear surface vertex to rear focal point).
    pub bfd: f64,
    /// Front principal plane offset (from front vertex, positive = into lens).
    pub front_principal: f64,
    /// Back principal plane offset (from rear vertex, positive = into lens).
    pub back_principal: f64,
}

/// Compute cardinal points for a thick lens.
///
/// `n` = refractive index, `r1`/`r2` = radii of curvature, `d` = center thickness.
#[inline]
pub fn cardinal_points(n: f64, r1: f64, r2: f64, thickness: f64) -> Result<CardinalPoints> {
    let f = thick_lens_focal_length(n, r1, r2, thickness)?;
    let nm1 = n - 1.0;

    // Individual surface powers
    let phi1 = nm1 / r1;
    let phi2 = -nm1 / r2;

    // Back focal distance: BFD = f · (1 − φ₁·d/n)
    let bfd = f * (1.0 - phi1 * thickness / n);
    // Front focal distance: FFD = −f · (1 − φ₂·d/n)
    let ffd = -f * (1.0 - phi2 * thickness / n);

    // Principal plane offsets from vertices
    // H (back principal) is at distance (f - BFD) from rear vertex
    let back_principal = f - bfd;
    // H' (front principal) is at distance (-f - FFD) from front vertex
    let front_principal = -f - ffd;

    Ok(CardinalPoints {
        focal_length: f,
        ffd,
        bfd,
        front_principal,
        back_principal,
    })
}

// ── Aperture, F-number, Numerical Aperture ────────────────────────────────

/// F-number (focal ratio): N = f / D.
///
/// `focal_length` and `aperture_diameter` in same units.
#[inline]
pub fn f_number(focal_length: f64, aperture_diameter: f64) -> f64 {
    focal_length / aperture_diameter
}

/// Aperture diameter from f-number: D = f / N.
#[inline]
pub fn aperture_from_f_number(focal_length: f64, f_num: f64) -> f64 {
    focal_length / f_num
}

/// Numerical aperture: NA = n · sin(θ).
///
/// `half_angle` is the half-angle of the maximum cone of light (radians).
/// `n` is the refractive index of the medium (1.0 for air).
#[inline]
pub fn numerical_aperture(n: f64, half_angle: f64) -> f64 {
    n * half_angle.sin()
}

/// Approximate NA from f-number (paraxial): NA ≈ 1 / (2·N).
///
/// Valid for small angles in air (n=1).
#[inline]
pub fn na_from_f_number(f_num: f64) -> f64 {
    1.0 / (2.0 * f_num)
}

/// Diffraction-limited angular resolution (Rayleigh criterion).
///
/// θ = 1.22 · λ / D (radians)
///
/// `wavelength` and `aperture_diameter` in same units.
#[inline]
pub fn diffraction_limit(wavelength: f64, aperture_diameter: f64) -> f64 {
    1.22 * wavelength / aperture_diameter
}

/// Diffraction-limited spot size (Airy disk radius).
///
/// r_airy = 1.22 · λ · N (where N is f-number)
///
/// `wavelength` in same units as desired output.
#[inline]
pub fn airy_disk_radius(wavelength: f64, f_num: f64) -> f64 {
    1.22 * wavelength * f_num
}

// ── Field of View ─────────────────────────────────────────────────────────

/// Field of view angle (full angle).
///
/// FOV = 2 · arctan(sensor_size / (2·f))
///
/// `sensor_size` and `focal_length` in same units. Returns radians.
#[inline]
pub fn field_of_view(sensor_size: f64, focal_length: f64) -> f64 {
    2.0 * (sensor_size / (2.0 * focal_length)).atan()
}

/// Diagonal field of view from sensor dimensions.
///
/// `width` and `height` are sensor dimensions, same units as `focal_length`.
#[inline]
pub fn field_of_view_diagonal(width: f64, height: f64, focal_length: f64) -> f64 {
    field_of_view(width.hypot(height), focal_length)
}

// ── MTF ───────────────────────────────────────────────────────────────────

/// Diffraction-limited cutoff spatial frequency (cycles per unit length).
///
/// f_cutoff = 1 / (λ · N)
///
/// `wavelength` and result in consistent units. Typically λ in mm → f in cycles/mm.
#[inline]
pub fn mtf_cutoff_frequency(wavelength: f64, f_num: f64) -> f64 {
    1.0 / (wavelength * f_num)
}

/// Diffraction-limited MTF for an ideal circular aperture.
///
/// MTF(ν) = (2/π) · [arccos(ν/ν_c) − (ν/ν_c)·√(1 − (ν/ν_c)²)]
///
/// `spatial_freq` is the spatial frequency, `cutoff_freq` is from `mtf_cutoff_frequency`.
/// Returns modulation (0.0–1.0). Returns 0.0 for frequencies above cutoff.
#[inline]
pub fn mtf_diffraction_limited(spatial_freq: f64, cutoff_freq: f64) -> f64 {
    let v = spatial_freq / cutoff_freq;
    if v >= 1.0 {
        return 0.0;
    }
    if v <= 0.0 {
        return 1.0;
    }
    const TWO_OVER_PI: f64 = 2.0 / std::f64::consts::PI;
    TWO_OVER_PI * (v.acos() - v * (1.0 - v * v).sqrt())
}

// ── Seidel Aberrations ────────────────────────────────────────────────────

/// Third-order (Seidel) aberration coefficients for a thin lens.
///
/// These describe the five monochromatic aberrations and are computed
/// from the lens shape factor and conjugate factor.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SeidelCoefficients {
    /// Spherical aberration coefficient (S₁). Causes on-axis blur.
    pub spherical: f64,
    /// Coma coefficient (S₂). Comet-shaped off-axis blur.
    pub coma: f64,
    /// Astigmatism coefficient (S₃). Different focal lengths for tangential/sagittal.
    pub astigmatism: f64,
    /// Field curvature coefficient (S₄). Petzval sum.
    pub field_curvature: f64,
    /// Distortion coefficient (S₅). Barrel/pincushion.
    pub distortion: f64,
}

/// Shape factor of a lens: q = (R₂ + R₁) / (R₂ − R₁).
///
/// q = 0 for equi-convex, q = 1 for plano-convex (flat on right),
/// q = -1 for plano-convex (flat on left).
#[inline]
pub fn shape_factor(r1: f64, r2: f64) -> f64 {
    let denom = r2 - r1;
    if denom.abs() < 1e-15 {
        return 0.0; // equal radii → symmetric, q=0
    }
    (r2 + r1) / denom
}

/// Conjugate factor (position factor): p = (di − do) / (di + do).
///
/// p = -1 for object at infinity, p = 0 for symmetric conjugates (2f-2f),
/// p = 1 for image at infinity.
#[inline]
pub fn conjugate_factor(object_distance: f64, image_distance: f64) -> f64 {
    (image_distance - object_distance) / (image_distance + object_distance)
}

/// Compute Seidel aberration coefficients for a thin lens in air.
///
/// `n` = refractive index, `f` = focal length, `q` = shape factor,
/// `p` = conjugate factor (position factor).
///
/// The coefficients are normalized to the lens power.
#[inline]
pub fn seidel_coefficients(n: f64, focal_length: f64, q: f64, p: f64) -> SeidelCoefficients {
    let phi = 1.0 / focal_length; // optical power

    // Spherical aberration: depends on shape and conjugate
    let spherical = phi.powi(3)
        * (n / (4.0 * (n - 1.0).powi(2)))
        * ((n + 2.0) / (n * (n - 1.0)) * q * q
            + (3.0 * n + 2.0) * (n - 1.0) / n * p * p
            + 4.0 * (n + 1.0) / (n * (n - 1.0)) * q * p);

    // Coma: linear in field angle
    let coma = phi.powi(2)
        * (1.0 / (2.0 * (n - 1.0)))
        * ((2.0 * (n + 1.0)) / (n * (n - 1.0)) * q + (3.0 * n + 1.0) / n * p);

    // Astigmatism: quadratic in field angle
    let astigmatism = phi;

    // Petzval field curvature
    let field_curvature = phi / n;

    // Distortion (zero for a single thin lens on-axis with stop at lens)
    let distortion = 0.0;

    SeidelCoefficients {
        spherical,
        coma,
        astigmatism,
        field_curvature,
        distortion,
    }
}

/// Longitudinal spherical aberration for a marginal ray.
///
/// LSA = h² / (2·f) · S₁_coeff
///
/// `ray_height` is the height at which the ray enters the lens.
/// Returns the axial shift of the focus point.
#[inline]
pub fn longitudinal_spherical_aberration(ray_height: f64, focal_length: f64, n: f64) -> f64 {
    // For a simple equi-convex lens (q=0, p=-1 for object at infinity)
    let phi = 1.0 / focal_length;
    let s1 = phi.powi(3) * n / (4.0 * (n - 1.0).powi(2)) * ((3.0 * n + 2.0) * (n - 1.0) / n);
    ray_height * ray_height * s1 / (2.0 * phi)
}

/// Chromatic aberration (longitudinal/axial color).
///
/// LCA = f / V
///
/// Returns the axial distance between F-line and C-line focal points.
/// `focal_length` in desired units, `abbe_v` (V number) from ray module.
#[inline]
pub fn chromatic_aberration(focal_length: f64, abbe_v: f64) -> f64 {
    focal_length / abbe_v
}

/// Petzval sum for a system of thin lenses.
///
/// S = Σ 1/(nᵢ·fᵢ)
///
/// The Petzval radius of curvature is R_p = -1/S.
/// A flat field requires S = 0.
pub fn petzval_sum(elements: &[(f64, f64)]) -> f64 {
    elements.iter().map(|(n, f)| 1.0 / (n * f)).sum()
}

/// Petzval radius of curvature from the Petzval sum.
///
/// R_p = -1 / S. Returns `None` if the sum is zero (flat field).
#[inline]
pub fn petzval_radius(petzval_sum: f64) -> Option<f64> {
    if petzval_sum.abs() < 1e-15 {
        None // flat field
    } else {
        Some(-1.0 / petzval_sum)
    }
}

// ── Multi-element System ──────────────────────────────────────────────────

/// Two thin lenses separated by distance `d`.
///
/// 1/f = 1/f₁ + 1/f₂ − d/(f₁·f₂)
///
/// Returns the effective focal length of the system.
#[inline]
pub fn separated_lenses_focal_length(f1: f64, f2: f64, separation: f64) -> Result<f64> {
    let power = 1.0 / f1 + 1.0 / f2 - separation / (f1 * f2);
    if power.abs() < 1e-15 {
        return Err(PrakashError::DivisionByZero {
            context: "separated lens system has zero power".into(),
        });
    }
    Ok(1.0 / power)
}

/// Back focal distance of a two-lens system.
///
/// BFD = f · (1 − d/f₁)
pub fn separated_lenses_bfd(f1: f64, f2: f64, separation: f64) -> Result<f64> {
    let f = separated_lenses_focal_length(f1, f2, separation)?;
    Ok(f * (1.0 - separation / f1))
}

/// System magnification for a multi-element system.
///
/// Total magnification is the product of individual magnifications.
pub fn system_magnification(magnifications: &[f64]) -> f64 {
    magnifications.iter().product()
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-6;

    // ── Thin lens tests ───────────────────────────────────────────────────

    #[test]
    fn test_thin_lens_real_image() {
        let di = thin_lens_image_distance(50.0, 100.0).unwrap();
        assert!((di - 100.0).abs() < EPS);
    }

    #[test]
    fn test_thin_lens_virtual_image() {
        let di = thin_lens_image_distance(50.0, 30.0).unwrap();
        assert!(di < 0.0);
    }

    #[test]
    fn test_thin_lens_at_focal() {
        assert!(thin_lens_image_distance(50.0, 50.0).is_err());
    }

    #[test]
    fn test_thin_lens_far_object() {
        // Object very far → image at focal length
        let di = thin_lens_image_distance(50.0, 1e6).unwrap();
        assert!((di - 50.0).abs() < 0.01);
    }

    #[test]
    fn test_thin_lens_diverging() {
        // Diverging lens (negative f) always produces virtual image
        let di = thin_lens_image_distance(-50.0, 100.0).unwrap();
        assert!(di < 0.0, "Diverging lens should produce virtual image");
    }

    #[test]
    fn test_thin_lens_satisfies_equation() {
        // Verify 1/f = 1/do + 1/di
        let f = 75.0;
        let do_ = 200.0;
        let di = thin_lens_image_distance(f, do_).unwrap();
        let lhs = 1.0 / f;
        let rhs = 1.0 / do_ + 1.0 / di;
        assert!((lhs - rhs).abs() < EPS);
    }

    // ── Magnification tests ───────────────────────────────────────────────

    #[test]
    fn test_magnification_at_2f() {
        let m = magnification(100.0, 100.0);
        assert!((m - (-1.0)).abs() < EPS);
    }

    #[test]
    fn test_magnification_magnified() {
        let m = magnification(60.0, 120.0);
        assert!(m.abs() > 1.0);
        assert!(m < 0.0); // inverted
    }

    #[test]
    fn test_magnification_diminished() {
        let m = magnification(200.0, 50.0);
        assert!(m.abs() < 1.0);
    }

    #[test]
    fn test_magnification_virtual_upright() {
        // Virtual image (negative di) → positive magnification (upright)
        let m = magnification(30.0, -60.0);
        assert!(m > 0.0, "Virtual image should be upright (positive mag)");
    }

    // ── Lensmaker tests ───────────────────────────────────────────────────

    #[test]
    fn test_lensmaker() {
        let f = lensmaker_focal_length(1.5, 100.0, -100.0).unwrap();
        assert!((f - 100.0).abs() < EPS);
    }

    #[test]
    fn test_lensmaker_planoconvex() {
        // Plano-convex: R1=100, R2=∞ → 1/f = (n-1)/R1
        let f = lensmaker_focal_length(1.5, 100.0, f64::INFINITY).unwrap();
        assert!((f - 200.0).abs() < EPS); // f = R1/(n-1) = 100/0.5 = 200
    }

    #[test]
    fn test_lensmaker_flat_surface_error() {
        // Two flat surfaces (R1=∞, R2=∞) → no focusing power
        assert!(lensmaker_focal_length(1.5, f64::INFINITY, f64::INFINITY).is_err());
    }

    #[test]
    fn test_lensmaker_diverging() {
        // Biconcave: R1=-100, R2=100 → negative focal length
        let f = lensmaker_focal_length(1.5, -100.0, 100.0).unwrap();
        assert!(f < 0.0, "Biconcave lens should have negative focal length");
    }

    // ── Optical power tests ───────────────────────────────────────────────

    #[test]
    fn test_optical_power() {
        assert!((optical_power(0.5).unwrap() - 2.0).abs() < EPS);
    }

    #[test]
    fn test_optical_power_negative() {
        assert!((optical_power(-0.5).unwrap() - (-2.0)).abs() < EPS);
    }

    #[test]
    fn test_optical_power_zero() {
        assert!(optical_power(0.0).is_err());
    }

    #[test]
    fn test_optical_power_roundtrip() {
        let f = 0.25;
        let d = optical_power(f).unwrap();
        assert!((1.0 / d - f).abs() < EPS);
    }

    // ── Mirror tests ──────────────────────────────────────────────────────

    #[test]
    fn test_mirror_focal() {
        assert!((mirror_focal_length(200.0) - 100.0).abs() < EPS);
    }

    #[test]
    fn test_mirror_focal_concave() {
        // Concave mirror (positive R) → positive focal length
        assert!(mirror_focal_length(100.0) > 0.0);
    }

    #[test]
    fn test_mirror_focal_convex() {
        // Convex mirror (negative R) → negative focal length
        assert!(mirror_focal_length(-100.0) < 0.0);
    }

    #[test]
    fn test_mirror_image_real() {
        let di = mirror_image_distance(100.0, 200.0).unwrap();
        assert!((di - 200.0).abs() < EPS);
    }

    #[test]
    fn test_mirror_image_at_focal() {
        assert!(mirror_image_distance(100.0, 100.0).is_err());
    }

    // ── Lens classification tests ─────────────────────────────────────────

    #[test]
    fn test_classify_lens() {
        assert_eq!(classify_lens(50.0), LensType::Converging);
        assert_eq!(classify_lens(-50.0), LensType::Diverging);
    }

    #[test]
    fn test_classify_lens_zero() {
        // f=0 → diverging by convention (else branch)
        assert_eq!(classify_lens(0.0), LensType::Diverging);
    }

    #[test]
    fn test_lens_type_serde_roundtrip() {
        let lt = LensType::Converging;
        let json = serde_json::to_string(&lt).unwrap();
        let back: LensType = serde_json::from_str(&json).unwrap();
        assert_eq!(back, lt);
    }

    // ── Combined focal length tests ───────────────────────────────────────

    #[test]
    fn test_combined_focal() {
        let f = combined_focal_length(100.0, 100.0).unwrap();
        assert!((f - 50.0).abs() < EPS);
    }

    #[test]
    fn test_combined_focal_converging_diverging() {
        let f = combined_focal_length(50.0, -100.0).unwrap();
        assert!((f - 100.0).abs() < EPS);
    }

    #[test]
    fn test_combined_focal_equal_opposite() {
        // Equal converging + diverging → zero combined power (error)
        assert!(combined_focal_length(100.0, -100.0).is_err());
    }

    // ── Depth of field tests ──────────────────────────────────────────────

    #[test]
    fn test_depth_of_field() {
        let (near, far) = depth_of_field(50.0, 2.8, 0.03, 2000.0);
        assert!(near < 2000.0);
        assert!(far > 2000.0);
        assert!(near > 0.0);
    }

    #[test]
    fn test_dof_wider_aperture_shallower() {
        let (near1, far1) = depth_of_field(50.0, 1.4, 0.03, 2000.0);
        let (near2, far2) = depth_of_field(50.0, 8.0, 0.03, 2000.0);
        let dof1 = far1 - near1;
        let dof2 = far2 - near2;
        assert!(
            dof1 < dof2,
            "Wider aperture (f/1.4) should have shallower DoF"
        );
    }

    #[test]
    fn test_dof_longer_lens_shallower() {
        let (near1, far1) = depth_of_field(85.0, 2.8, 0.03, 2000.0);
        let (near2, far2) = depth_of_field(35.0, 2.8, 0.03, 2000.0);
        let dof1 = far1 - near1;
        let dof2 = far2 - near2;
        assert!(dof1 < dof2, "Longer focal length should have shallower DoF");
    }

    #[test]
    fn test_dof_symmetric_around_subject() {
        let (near, far) = depth_of_field(50.0, 5.6, 0.03, 5000.0);
        assert!(near < 5000.0);
        assert!(far > 5000.0);
    }

    // ── Thick lens tests ──────────────────────────────────────────────────

    #[test]
    fn test_thick_lens_reduces_to_thin_at_zero_thickness() {
        let n = 1.5;
        let r1 = 100.0;
        let r2 = -100.0;
        let f_thin = lensmaker_focal_length(n, r1, r2).unwrap();
        let f_thick = thick_lens_focal_length(n, r1, r2, 0.0).unwrap();
        assert!(
            (f_thin - f_thick).abs() < EPS,
            "Zero thickness should match thin lens: thin={f_thin}, thick={f_thick}"
        );
    }

    #[test]
    fn test_thick_lens_biconvex() {
        let f = thick_lens_focal_length(1.5, 100.0, -100.0, 10.0).unwrap();
        assert!(f > 0.0, "Biconvex thick lens should be converging");
        // Thick lens should have slightly different f than thin
        let f_thin = lensmaker_focal_length(1.5, 100.0, -100.0).unwrap();
        assert!((f - f_thin).abs() > 0.01, "Thick and thin should differ");
    }

    #[test]
    fn test_thick_lens_zero_power() {
        // Symmetric meniscus with specific params can have zero power
        assert!(thick_lens_focal_length(1.5, f64::INFINITY, f64::INFINITY, 10.0).is_err());
    }

    #[test]
    fn test_cardinal_points_biconvex() {
        let cp = cardinal_points(1.5, 100.0, -100.0, 10.0).unwrap();
        assert!(cp.focal_length > 0.0);
        assert!(cp.bfd > 0.0, "BFD should be positive for biconvex");
        assert!(cp.ffd < 0.0, "FFD should be negative for biconvex");
    }

    #[test]
    fn test_cardinal_points_symmetric_lens() {
        // Symmetric biconvex: principal planes should be symmetric
        let cp = cardinal_points(1.5, 100.0, -100.0, 10.0).unwrap();
        assert!(
            (cp.front_principal.abs() - cp.back_principal.abs()).abs() < 0.1,
            "Symmetric lens should have symmetric principal planes"
        );
    }

    // ── F-number and NA tests ─────────────────────────────────────────────

    #[test]
    fn test_f_number() {
        assert!((f_number(50.0, 25.0) - 2.0).abs() < EPS);
        assert!((f_number(100.0, 50.0) - 2.0).abs() < EPS);
    }

    #[test]
    fn test_aperture_from_f_number() {
        assert!((aperture_from_f_number(50.0, 2.0) - 25.0).abs() < EPS);
    }

    #[test]
    fn test_f_number_aperture_roundtrip() {
        let f = 85.0;
        let d = 30.357;
        let n = f_number(f, d);
        let d_back = aperture_from_f_number(f, n);
        assert!((d_back - d).abs() < EPS);
    }

    #[test]
    fn test_numerical_aperture() {
        use std::f64::consts::FRAC_PI_6;
        let na = numerical_aperture(1.0, FRAC_PI_6);
        assert!((na - 0.5).abs() < EPS); // sin(30°) = 0.5
    }

    #[test]
    fn test_na_from_f_number() {
        let na = na_from_f_number(2.0);
        assert!((na - 0.25).abs() < EPS);
    }

    #[test]
    fn test_na_immersion_oil() {
        // Oil immersion (n=1.515) with wide cone
        use std::f64::consts::FRAC_PI_3;
        let na = numerical_aperture(1.515, FRAC_PI_3);
        assert!(na > 1.0, "Oil immersion NA can exceed 1.0");
    }

    #[test]
    fn test_diffraction_limit() {
        // 550nm light through 10mm aperture
        let theta = diffraction_limit(550e-6, 10.0); // mm units
        assert!(theta > 0.0);
        assert!(theta < 0.001); // very small angle
    }

    #[test]
    fn test_larger_aperture_better_resolution() {
        let theta_small = diffraction_limit(550e-6, 10.0);
        let theta_large = diffraction_limit(550e-6, 50.0);
        assert!(theta_large < theta_small);
    }

    #[test]
    fn test_airy_disk_radius() {
        let r = airy_disk_radius(0.00055, 2.8); // 550nm = 0.00055mm, f/2.8
        assert!(r > 0.0);
        // Airy disk at f/2.8 with green light ≈ 1.88 μm = 0.00188 mm
        assert!((r - 0.001_88).abs() < 0.001);
    }

    // ── Field of view tests ───────────────────────────────────────────────

    #[test]
    fn test_field_of_view_50mm_fullframe() {
        // 50mm on full-frame (36mm sensor width) ≈ 39.6° horizontal FOV
        let fov = field_of_view(36.0, 50.0);
        let fov_deg = fov.to_degrees();
        assert!(
            (fov_deg - 39.6).abs() < 1.0,
            "50mm on FF should be ≈39.6° FOV, got {fov_deg}"
        );
    }

    #[test]
    fn test_field_of_view_longer_lens_narrower() {
        let fov_50 = field_of_view(36.0, 50.0);
        let fov_200 = field_of_view(36.0, 200.0);
        assert!(fov_200 < fov_50);
    }

    #[test]
    fn test_field_of_view_diagonal() {
        let fov_diag = field_of_view_diagonal(36.0, 24.0, 50.0);
        let fov_horiz = field_of_view(36.0, 50.0);
        assert!(fov_diag > fov_horiz, "Diagonal FOV > horizontal FOV");
    }

    // ── MTF tests ─────────────────────────────────────────────────────────

    #[test]
    fn test_mtf_cutoff() {
        // 550nm, f/2.8 → cutoff = 1/(0.00055 * 2.8) ≈ 649 cycles/mm
        let fc = mtf_cutoff_frequency(0.00055, 2.8);
        assert!((fc - 649.0).abs() < 5.0, "Cutoff ≈ 649 cy/mm, got {fc}");
    }

    #[test]
    fn test_mtf_at_zero_frequency() {
        let mtf = mtf_diffraction_limited(0.0, 1000.0);
        assert!((mtf - 1.0).abs() < EPS, "MTF(0) should be 1.0");
    }

    #[test]
    fn test_mtf_at_cutoff() {
        let mtf = mtf_diffraction_limited(1000.0, 1000.0);
        assert!(mtf.abs() < EPS, "MTF at cutoff should be 0.0");
    }

    #[test]
    fn test_mtf_above_cutoff() {
        let mtf = mtf_diffraction_limited(1500.0, 1000.0);
        assert!(mtf.abs() < EPS, "MTF above cutoff should be 0.0");
    }

    #[test]
    fn test_mtf_monotonic_decrease() {
        let fc = 1000.0;
        let mut prev = 1.0;
        for i in 1..=10 {
            let freq = fc * (i as f64) / 10.0;
            let m = mtf_diffraction_limited(freq, fc);
            assert!(
                m <= prev + EPS,
                "MTF should decrease monotonically: at {freq} got {m} > {prev}"
            );
            prev = m;
        }
    }

    #[test]
    fn test_mtf_range() {
        let fc = 1000.0;
        for i in 0..=20 {
            let freq = fc * (i as f64) / 20.0;
            let m = mtf_diffraction_limited(freq, fc);
            assert!((0.0..=1.0).contains(&m), "MTF out of range at {freq}: {m}");
        }
    }

    // ── Seidel aberration tests ───────────────────────────────────────────

    #[test]
    fn test_shape_factor_equiconvex() {
        // R1 = 100, R2 = -100 → q = (-100 + 100) / (-100 - 100) = 0
        let q = shape_factor(100.0, -100.0);
        assert!(q.abs() < EPS, "Equi-convex should have q=0, got {q}");
    }

    #[test]
    fn test_shape_factor_planoconvex() {
        // Plano-convex: R1 = 100, R2 = very large (flat)
        let q = shape_factor(100.0, 1e15);
        assert!((q - 1.0).abs() < 0.01, "Plano-convex q ≈ 1, got {q}");
    }

    #[test]
    fn test_conjugate_factor_infinity() {
        // Object at infinity: di is finite, do → ∞ → p → -1
        let p = conjugate_factor(1e10, 50.0);
        assert!((p - (-1.0)).abs() < 0.001);
    }

    #[test]
    fn test_conjugate_factor_symmetric() {
        // 2f-2f conjugates: do = di → p = 0
        let p = conjugate_factor(100.0, 100.0);
        assert!(p.abs() < EPS);
    }

    #[test]
    fn test_seidel_equiconvex_at_infinity() {
        // Equi-convex (q=0), object at infinity (p=-1)
        let sc = seidel_coefficients(1.5, 100.0, 0.0, -1.0);
        assert!(
            sc.spherical != 0.0,
            "Spherical aberration should be nonzero"
        );
        assert!(sc.coma != 0.0, "Coma should be nonzero");
        assert!(
            sc.field_curvature > 0.0,
            "Field curvature should be positive"
        );
    }

    #[test]
    fn test_seidel_best_form_reduces_spherical() {
        // The best-form (minimum spherical) shape factor for object at infinity
        // is q_best ≈ -(2(n²-1))/(n+2) for p=-1
        let n = 1.5;
        let f = 100.0;
        let p = -1.0;

        let sc_equi = seidel_coefficients(n, f, 0.0, p);

        // Try a range of shape factors — minimum spherical should exist
        let mut min_sa = f64::MAX;
        let mut best_q = 0.0;
        for i in -20..=20 {
            let q = i as f64 * 0.1;
            let sc = seidel_coefficients(n, f, q, p);
            if sc.spherical.abs() < min_sa {
                min_sa = sc.spherical.abs();
                best_q = q;
            }
        }
        let sc_best = seidel_coefficients(n, f, best_q, p);

        assert!(
            sc_best.spherical.abs() <= sc_equi.spherical.abs(),
            "Best-form (q={best_q}) should have ≤ spherical aberration than equi-convex"
        );
    }

    #[test]
    fn test_longitudinal_spherical_increases_with_height() {
        let lsa_low = longitudinal_spherical_aberration(5.0, 100.0, 1.5);
        let lsa_high = longitudinal_spherical_aberration(10.0, 100.0, 1.5);
        assert!(
            lsa_high.abs() > lsa_low.abs(),
            "LSA should increase with ray height"
        );
    }

    #[test]
    fn test_chromatic_aberration() {
        let tca = chromatic_aberration(100.0, 64.0);
        assert!(tca > 0.0);
        // TCA = f/V = 100/64 ≈ 1.56
        assert!((tca - 1.5625).abs() < EPS);
    }

    #[test]
    fn test_low_dispersion_less_chromatic() {
        let tca_crown = chromatic_aberration(100.0, 64.0); // crown glass
        let tca_flint = chromatic_aberration(100.0, 26.0); // flint glass
        assert!(tca_crown < tca_flint, "Crown should have less TCA");
    }

    #[test]
    fn test_petzval_sum_single_lens() {
        let ps = petzval_sum(&[(1.5, 100.0)]);
        assert!((ps - 1.0 / 150.0).abs() < EPS);
    }

    #[test]
    fn test_petzval_sum_doublet() {
        // Crown (n=1.5, f=60) + flint (n=1.7, f=-90)
        // Petzval sum = 1/(1.5*60) + 1/(1.7*(-90))
        let ps = petzval_sum(&[(1.5, 60.0), (1.7, -90.0)]);
        // = 0.01111 - 0.00654 = 0.00457 (reduced but not zero)
        assert!(ps > 0.0 && ps < 0.01, "Doublet reduces Petzval sum");
    }

    #[test]
    fn test_petzval_radius_from_sum() {
        let r = petzval_radius(0.01).unwrap();
        assert!((r - (-100.0)).abs() < EPS);
    }

    #[test]
    fn test_petzval_radius_flat_field() {
        assert!(petzval_radius(0.0).is_none());
    }

    // ── Multi-element system tests ────────────────────────────────────────

    #[test]
    fn test_separated_lenses_contact() {
        // d=0 should match combined_focal_length
        let f_sep = separated_lenses_focal_length(100.0, 100.0, 0.0).unwrap();
        let f_contact = combined_focal_length(100.0, 100.0).unwrap();
        assert!((f_sep - f_contact).abs() < EPS);
    }

    #[test]
    fn test_separated_lenses_telephoto() {
        // Telephoto: positive front + negative rear, separated
        let f = separated_lenses_focal_length(50.0, -80.0, 30.0).unwrap();
        assert!(
            f > 50.0,
            "Telephoto should have longer effective focal length"
        );
    }

    #[test]
    fn test_separated_lenses_bfd() {
        let bfd = separated_lenses_bfd(50.0, -80.0, 30.0).unwrap();
        let f = separated_lenses_focal_length(50.0, -80.0, 30.0).unwrap();
        assert!(bfd < f, "Telephoto BFD should be shorter than focal length");
    }

    #[test]
    fn test_system_magnification() {
        let m = system_magnification(&[-0.5, -2.0, -0.3]);
        assert!((m - (-0.3)).abs() < EPS); // product = -0.3
    }

    #[test]
    fn test_system_magnification_single() {
        let m = system_magnification(&[-1.0]);
        assert!((m - (-1.0)).abs() < EPS);
    }

    #[test]
    fn test_system_magnification_empty() {
        let m = system_magnification(&[]);
        assert!((m - 1.0).abs() < EPS); // empty product = 1
    }
}
