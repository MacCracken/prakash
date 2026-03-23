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
pub fn optical_power(focal_length_m: f64) -> f64 {
    1.0 / focal_length_m
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
pub fn combined_focal_length(f1: f64, f2: f64) -> f64 {
    (f1 * f2) / (f1 + f2)
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

    #[test]
    fn test_thin_lens_real_image() {
        // f=50mm, object at 100mm → image at 100mm (2f → 2f)
        let di = thin_lens_image_distance(50.0, 100.0).unwrap();
        assert!((di - 100.0).abs() < EPS);
    }

    #[test]
    fn test_thin_lens_virtual_image() {
        // f=50mm, object at 30mm (inside focal) → virtual image
        let di = thin_lens_image_distance(50.0, 30.0).unwrap();
        assert!(di < 0.0); // virtual image
    }

    #[test]
    fn test_thin_lens_at_focal() {
        // Object at focal point → image at infinity (error)
        assert!(thin_lens_image_distance(50.0, 50.0).is_err());
    }

    #[test]
    fn test_magnification_at_2f() {
        let m = magnification(100.0, 100.0);
        assert!((m - (-1.0)).abs() < EPS); // inverted, same size
    }

    #[test]
    fn test_magnification_magnified() {
        let m = magnification(60.0, 120.0);
        assert!(m.abs() > 1.0); // magnified
    }

    #[test]
    fn test_lensmaker() {
        // Biconvex lens, n=1.5, R1=100mm, R2=-100mm
        let f = lensmaker_focal_length(1.5, 100.0, -100.0).unwrap();
        // 1/f = 0.5 * (1/100 + 1/100) = 0.5 * 0.02 = 0.01 → f = 100mm
        assert!((f - 100.0).abs() < EPS);
    }

    #[test]
    fn test_optical_power() {
        assert!((optical_power(0.5) - 2.0).abs() < EPS); // 500mm → 2 diopters
    }

    #[test]
    fn test_mirror_focal() {
        assert!((mirror_focal_length(200.0) - 100.0).abs() < EPS);
    }

    #[test]
    fn test_classify_lens() {
        assert_eq!(classify_lens(50.0), LensType::Converging);
        assert_eq!(classify_lens(-50.0), LensType::Diverging);
    }

    #[test]
    fn test_combined_focal() {
        // Two 100mm lenses → 50mm combined
        let f = combined_focal_length(100.0, 100.0);
        assert!((f - 50.0).abs() < EPS);
    }

    #[test]
    fn test_depth_of_field() {
        let (near, far) = depth_of_field(50.0, 2.8, 0.03, 2000.0);
        assert!(near < 2000.0);
        assert!(far > 2000.0);
        assert!(near > 0.0);
    }

    #[test]
    fn test_mirror_image_real() {
        let di = mirror_image_distance(100.0, 200.0).unwrap();
        assert!((di - 200.0).abs() < EPS);
    }
}
