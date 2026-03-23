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
    let far = (h * subject_dist) / (h - (subject_dist - focal_length));
    (near, far)
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
        // At moderate distances, DoF should bracket the subject
        let (near, far) = depth_of_field(50.0, 5.6, 0.03, 5000.0);
        assert!(near < 5000.0);
        assert!(far > 5000.0);
    }
}
