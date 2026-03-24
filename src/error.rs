//! Error types for prakash.

use std::borrow::Cow;

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum PrakashError {
    #[error(
        "total internal reflection: incident angle {angle_deg:.1}° exceeds critical angle {critical_deg:.1}° (n1={n1}, n2={n2})"
    )]
    TotalInternalReflection {
        angle_deg: f64,
        critical_deg: f64,
        n1: f64,
        n2: f64,
    },

    #[error("invalid refractive index: {n} (must be >= 1.0)")]
    InvalidRefractiveIndex { n: f64 },

    #[error("wavelength out of visible range: {nm} nm (visible: 380-780 nm)")]
    WavelengthOutOfRange { nm: f64 },

    #[error("invalid angle: {degrees}° (must be 0-90 for surface interactions)")]
    InvalidAngle { degrees: f64 },

    #[error("negative focal length not supported for this lens type: {focal_mm} mm")]
    InvalidFocalLength { focal_mm: f64 },

    #[error("division by zero in optical calculation: {context}")]
    DivisionByZero { context: Cow<'static, str> },

    #[error("invalid parameter: {reason}")]
    InvalidParameter { reason: Cow<'static, str> },
}

pub type Result<T> = std::result::Result<T, PrakashError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tir_error() {
        let e = PrakashError::TotalInternalReflection {
            angle_deg: 45.0,
            critical_deg: 41.8,
            n1: 1.5,
            n2: 1.0,
        };
        let msg = e.to_string();
        assert!(msg.contains("45.0"));
        assert!(msg.contains("41.8"));
        assert!(msg.contains("n1=1.5"));
        assert!(msg.contains("n2=1"));
        assert!(msg.contains("total internal reflection"));
    }

    #[test]
    fn test_invalid_index() {
        let e = PrakashError::InvalidRefractiveIndex { n: 0.5 };
        let msg = e.to_string();
        assert!(msg.contains("0.5"));
        assert!(msg.contains("must be >= 1.0"));
    }

    #[test]
    fn test_wavelength_range() {
        let e = PrakashError::WavelengthOutOfRange { nm: 200.0 };
        let msg = e.to_string();
        assert!(msg.contains("200"));
        assert!(msg.contains("380-780"));
    }

    #[test]
    fn test_invalid_focal() {
        let e = PrakashError::InvalidFocalLength { focal_mm: -50.0 };
        assert!(e.to_string().contains("-50"));
    }

    #[test]
    fn test_invalid_angle() {
        let e = PrakashError::InvalidAngle { degrees: 95.0 };
        let msg = e.to_string();
        assert!(msg.contains("95"));
        assert!(msg.contains("0-90"));
    }

    #[test]
    fn test_division_by_zero() {
        let e = PrakashError::DivisionByZero {
            context: "test context".into(),
        };
        assert!(e.to_string().contains("test context"));
    }

    #[test]
    fn test_invalid_parameter() {
        let e = PrakashError::InvalidParameter {
            reason: "bad input".into(),
        };
        assert!(e.to_string().contains("bad input"));
    }

    #[test]
    fn test_result_alias_ok() {
        let ok: Result<f64> = Ok(1.5);
        match ok {
            Ok(v) => assert_eq!(v, 1.5),
            Err(_) => panic!("expected Ok"),
        }
    }

    #[test]
    fn test_result_alias_err() {
        let err: Result<f64> = Err(PrakashError::InvalidRefractiveIndex { n: 0.0 });
        assert!(err.is_err());
    }

    #[test]
    fn test_error_is_debug() {
        let e = PrakashError::InvalidRefractiveIndex { n: 0.5 };
        let debug = format!("{:?}", e);
        assert!(debug.contains("InvalidRefractiveIndex"));
    }
}
