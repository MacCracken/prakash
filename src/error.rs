//! Error types for prakash.

use std::borrow::Cow;

/// Errors that can occur during optical computations.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum PrakashError {
    /// Total internal reflection — incident angle exceeds the critical angle.
    #[error(
        "total internal reflection: incident angle {angle_deg:.1}° exceeds critical angle {critical_deg:.1}° (n1={n1}, n2={n2})"
    )]
    TotalInternalReflection {
        /// Incident angle in degrees.
        angle_deg: f64,
        /// Critical angle in degrees.
        critical_deg: f64,
        /// Refractive index of the incident medium.
        n1: f64,
        /// Refractive index of the transmitting medium.
        n2: f64,
    },

    /// Refractive index is below the physical minimum of 1.0.
    #[error("invalid refractive index: {n} (must be >= 1.0)")]
    InvalidRefractiveIndex {
        /// The invalid refractive index value.
        n: f64,
    },

    /// Wavelength is outside the visible range (380–780 nm).
    #[error("wavelength out of visible range: {nm} nm (visible: 380-780 nm)")]
    WavelengthOutOfRange {
        /// The out-of-range wavelength in nanometers.
        nm: f64,
    },

    /// Angle is outside the valid range for the operation.
    #[error("invalid angle: {degrees}° (must be 0-90 for surface interactions)")]
    InvalidAngle {
        /// The invalid angle in degrees.
        degrees: f64,
    },

    /// Negative focal length where a positive one is required.
    #[error("negative focal length not supported for this lens type: {focal_mm} mm")]
    InvalidFocalLength {
        /// The invalid focal length in millimeters.
        focal_mm: f64,
    },

    /// Division by zero in an optical calculation.
    #[error("division by zero in optical calculation: {context}")]
    DivisionByZero {
        /// Description of the calculation context.
        context: Cow<'static, str>,
    },

    /// A parameter is outside the valid range for the operation.
    #[error("invalid parameter: {reason}")]
    InvalidParameter {
        /// Explanation of why the parameter is invalid.
        reason: Cow<'static, str>,
    },
}

/// Convenience alias for `Result<T, PrakashError>`.
pub type Result<T> = std::result::Result<T, PrakashError>;

#[cfg(feature = "bijli-backend")]
impl From<bijli::BijliError> for PrakashError {
    fn from(e: bijli::BijliError) -> Self {
        match e {
            bijli::BijliError::DivisionByZero { context } => PrakashError::DivisionByZero {
                context: context.into(),
            },
            bijli::BijliError::InvalidParameter { reason } => PrakashError::InvalidParameter {
                reason: reason.into(),
            },
            other => PrakashError::InvalidParameter {
                reason: other.to_string().into(),
            },
        }
    }
}

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
