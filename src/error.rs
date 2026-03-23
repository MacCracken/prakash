//! Error types for prakash.

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum PrakashError {
    #[error("total internal reflection: incident angle {angle_deg:.1}° exceeds critical angle {critical_deg:.1}° (n1={n1}, n2={n2})")]
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
    DivisionByZero { context: String },

    #[error("invalid parameter: {reason}")]
    InvalidParameter { reason: String },
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
        assert!(e.to_string().contains("45.0"));
        assert!(e.to_string().contains("41.8"));
    }

    #[test]
    fn test_invalid_index() {
        let e = PrakashError::InvalidRefractiveIndex { n: 0.5 };
        assert!(e.to_string().contains("0.5"));
    }

    #[test]
    fn test_wavelength_range() {
        let e = PrakashError::WavelengthOutOfRange { nm: 200.0 };
        assert!(e.to_string().contains("200"));
    }

    #[test]
    fn test_invalid_focal() {
        let e = PrakashError::InvalidFocalLength { focal_mm: -50.0 };
        assert!(e.to_string().contains("-50"));
    }

    #[test]
    fn test_result_alias() {
        let ok: Result<f64> = Ok(1.5);
        assert!(ok.is_ok());
    }
}
