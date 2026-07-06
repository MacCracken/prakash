//! Fiber optics — numerical aperture, modes, coupling efficiency.

use std::f64::consts::PI;

use crate::error::{PrakashError, Result};

/// Numerical aperture of an optical fiber.
///
/// NA = √(n_core² − n_clad²)
///
/// `n_core` and `n_clad` are the refractive indices of the core and cladding.
#[must_use = "returns the computed NA"]
#[inline]
pub fn fiber_na(n_core: f64, n_clad: f64) -> Result<f64> {
    let na2 = n_core * n_core - n_clad * n_clad;
    if na2 < 0.0 {
        return Err(PrakashError::InvalidParameter {
            reason: "core index must be greater than cladding index".into(),
        });
    }
    Ok(na2.sqrt())
}

/// V-number (normalized frequency) of an optical fiber.
///
/// V = 2π · a · NA / λ
///
/// where `a` is the core radius. Single-mode cutoff at V < 2.405.
///
/// `core_radius` and `wavelength` in same units.
#[must_use = "returns the V-number"]
#[inline]
pub fn v_number(core_radius: f64, na: f64, wavelength: f64) -> f64 {
    2.0 * PI * core_radius * na / wavelength
}

/// Approximate number of guided modes in a step-index multimode fiber.
///
/// N ≈ V² / 2
///
/// Only valid for V >> 2.405 (many modes). For single-mode fibers (V < 2.405),
/// there is exactly 1 guided mode (LP01 / HE11).
#[must_use]
#[inline]
pub fn num_modes(v: f64) -> u32 {
    if v < 2.405 {
        1
    } else {
        ((v * v / 2.0) as u32).max(1)
    }
}

/// Whether the fiber supports only a single mode at this V-number.
///
/// Single-mode cutoff: V < 2.405 (first zero of J₀ Bessel function).
#[must_use]
#[inline]
pub fn is_single_mode(v: f64) -> bool {
    v < 2.405
}

/// Mode field diameter (MFD) of a single-mode fiber (Marcuse approximation).
///
/// MFD ≈ 2a · (0.65 + 1.619·V^{−3/2} + 2.879·V^{−6})
///
/// `core_radius` in same units as desired output. Valid for 0.8 < V < 2.5.
#[must_use]
#[inline]
pub fn mode_field_diameter(core_radius: f64, v: f64) -> f64 {
    2.0 * core_radius * (0.65 + 1.619 * v.powf(-1.5) + 2.879 * v.powi(-6))
}

/// Coupling efficiency between a Gaussian beam and a single-mode fiber.
///
/// η = (2·w_beam·w_fiber / (w_beam² + w_fiber²))²
///
/// where `w_beam` and `w_fiber` are the 1/e² beam radii (half of beam diameter
/// and half of MFD respectively). Perfect coupling (η=1) when w_beam = w_fiber.
///
/// Does not account for angular or positional misalignment.
#[must_use]
#[inline]
pub fn coupling_efficiency_gaussian(w_beam: f64, w_fiber: f64) -> f64 {
    let num = 2.0 * w_beam * w_fiber;
    let den = w_beam * w_beam + w_fiber * w_fiber;
    let ratio = num / den;
    ratio * ratio
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-6;

    #[test]
    fn test_fiber_na_smf28() {
        // SMF-28: n_core ≈ 1.4682, n_clad ≈ 1.4629 → NA ≈ 0.12
        let na = fiber_na(1.4682, 1.4629).unwrap();
        assert!((na - 0.12).abs() < 0.02, "SMF-28 NA ≈ 0.12, got {na}");
    }

    #[test]
    fn test_fiber_na_invalid() {
        assert!(fiber_na(1.45, 1.47).is_err()); // core < clad
    }

    #[test]
    fn test_v_number_single_mode() {
        // SMF-28: core_radius=4.1μm, NA=0.12, λ=1.55μm → V ≈ 1.99
        let na = 0.12;
        let v = v_number(4.1e-6, na, 1.55e-6);
        assert!(v < 2.405, "SMF-28 at 1550nm should be single-mode, V={v}");
    }

    #[test]
    fn test_v_number_multimode() {
        // Large core fiber: radius=25μm, NA=0.22, λ=0.85μm
        let v = v_number(25e-6, 0.22, 0.85e-6);
        assert!(v > 2.405, "25μm core should be multimode, V={v}");
    }

    #[test]
    fn test_num_modes_single() {
        assert_eq!(num_modes(2.0), 1);
    }

    #[test]
    fn test_num_modes_multi() {
        let v = v_number(25e-6, 0.22, 0.85e-6);
        let n = num_modes(v);
        assert!(n > 100, "Many modes expected, got {n}");
    }

    #[test]
    fn test_is_single_mode() {
        assert!(is_single_mode(2.0));
        assert!(!is_single_mode(3.0));
    }

    #[test]
    fn test_mfd_typical() {
        // SMF-28 at 1550nm: MFD ≈ 10.4μm
        let v = v_number(4.1e-6, 0.12, 1.55e-6);
        let mfd = mode_field_diameter(4.1e-6, v);
        assert!(
            (mfd * 1e6 - 10.4).abs() < 2.0,
            "MFD ≈ 10.4μm, got {:.1}μm",
            mfd * 1e6
        );
    }

    #[test]
    fn test_coupling_perfect() {
        // Equal beam and fiber radii → perfect coupling
        let eta = coupling_efficiency_gaussian(5e-6, 5e-6);
        assert!((eta - 1.0).abs() < EPS);
    }

    #[test]
    fn test_coupling_mismatch() {
        // Mismatched radii → reduced coupling
        let eta = coupling_efficiency_gaussian(5e-6, 10e-6);
        assert!(eta < 1.0);
        assert!(eta > 0.5);
    }

    #[test]
    fn test_coupling_symmetric() {
        let e1 = coupling_efficiency_gaussian(3e-6, 5e-6);
        let e2 = coupling_efficiency_gaussian(5e-6, 3e-6);
        assert!((e1 - e2).abs() < EPS, "Coupling should be symmetric");
    }
}
