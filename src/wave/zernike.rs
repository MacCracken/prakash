//! Zernike polynomials for wavefront decomposition and aberration analysis.
//!
//! Zernike polynomials form an orthogonal basis on the unit disk, making them
//! the standard representation for wavefront aberrations in optical systems.
//!
//! # Ordering conventions
//!
//! This module uses **Noll sequential ordering** (J. Opt. Soc. Am. 66, 207, 1976)
//! as the primary indexing scheme. Functions accept `(n, m)` radial/azimuthal
//! indices directly.
//!
//! | Noll j | (n, m) | Name |
//! |--------|--------|------|
//! | 1 | (0, 0) | Piston |
//! | 2 | (1, 1) | Tilt X (cos) |
//! | 3 | (1, −1) | Tilt Y (sin) |
//! | 4 | (2, 0) | Defocus |
//! | 5 | (2, −2) | Astigmatism 45° |
//! | 6 | (2, 2) | Astigmatism 0° |
//! | 7 | (3, −1) | Coma Y |
//! | 8 | (3, 1) | Coma X |
//! | 9 | (3, −3) | Trefoil Y |
//! | 10 | (3, 3) | Trefoil X |
//! | 11 | (4, 0) | Spherical |

use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

// ── Radial polynomial ───────────────────────────────────────────────────────

/// Evaluate the radial Zernike polynomial R_n^m(ρ).
///
/// R_n^m(ρ) = Σ_{s=0}^{(n-|m|)/2} (-1)^s · (n-s)! / (s! · ((n+|m|)/2 - s)! · ((n-|m|)/2 - s)!) · ρ^(n-2s)
///
/// `rho` must be in [0, 1]. `n` and `m` must satisfy n ≥ |m| and (n - |m|) even.
#[must_use]
#[inline]
pub fn radial_polynomial(n: u32, m_abs: u32, rho: f64) -> f64 {
    if (n < m_abs) || !(n - m_abs).is_multiple_of(2) {
        return 0.0;
    }
    let half = ((n - m_abs) / 2) as usize;
    let mut sum = 0.0;
    for s in 0..=half {
        let sign = if s.is_multiple_of(2) { 1.0 } else { -1.0 };
        let num = factorial((n as usize) - s);
        let d1 = factorial(s);
        let d2 = factorial(((n + m_abs) / 2) as usize - s);
        let d3 = factorial(half - s);
        sum +=
            sign * (num as f64) / ((d1 * d2 * d3) as f64) * rho.powi((n as i32) - 2 * (s as i32));
    }
    sum
}

/// Factorial for small arguments (sufficient for Zernike up to order ~20).
#[inline]
fn factorial(n: usize) -> u64 {
    match n {
        0 | 1 => 1,
        2 => 2,
        3 => 6,
        4 => 24,
        5 => 120,
        6 => 720,
        7 => 5040,
        8 => 40320,
        9 => 362880,
        10 => 3628800,
        11 => 39916800,
        12 => 479001600,
        13 => 6227020800,
        14 => 87178291200,
        15 => 1307674368000,
        16 => 20922789888000,
        17 => 355687428096000,
        18 => 6402373705728000,
        19 => 121645100408832000,
        20 => 2432902008176640000,
        _ => (2..=n as u64).product(),
    }
}

// ── Full Zernike polynomial ─────────────────────────────────────────────────

/// Evaluate a single Zernike polynomial Z_n^m(ρ, θ) at a point on the unit disk.
///
/// The normalization follows the convention where:
/// - m > 0: Z = N · R_n^m(ρ) · cos(mθ)
/// - m < 0: Z = N · R_n^{|m|}(ρ) · sin(|m|θ)
/// - m = 0: Z = N · R_n^0(ρ)
///
/// where N = √(2(n+1)) for m ≠ 0, N = √(n+1) for m = 0.
///
/// `rho` in [0, 1], `theta` in radians, `n` = radial order, `m` = azimuthal frequency.
/// Returns 0.0 for points outside the unit disk (ρ > 1).
#[must_use]
#[inline]
pub fn zernike(n: u32, m: i32, rho: f64, theta: f64) -> f64 {
    if rho > 1.0 {
        return 0.0;
    }
    let m_abs = m.unsigned_abs();
    let r = radial_polynomial(n, m_abs, rho);
    let norm = if m == 0 {
        ((n + 1) as f64).sqrt()
    } else {
        (2.0 * (n + 1) as f64).sqrt()
    };
    let angular = if m > 0 {
        (m_abs as f64 * theta).cos()
    } else if m < 0 {
        (m_abs as f64 * theta).sin()
    } else {
        1.0
    };
    norm * r * angular
}

// ── Noll index conversion ───────────────────────────────────────────────────

/// Convert Noll sequential index j (1-based) to (n, m) indices.
///
/// Follows Noll 1976 convention (J. Opt. Soc. Am. 66, 207).
///
/// The Noll ordering enumerates Zernike polynomials sequentially:
/// j=1→piston, j=2→tilt Y, j=3→tilt X, j=4→defocus, etc.
#[must_use]
pub fn noll_to_nm(j: u32) -> (u32, i32) {
    if j == 0 {
        return (0, 0);
    }
    // Find radial order n: n(n+1)/2 < j <= (n+1)(n+2)/2
    let mut n = 0u32;
    while (n + 1) * (n + 2) / 2 < j {
        n += 1;
    }
    // Remainder within this order (0-based)
    let k = j - n * (n + 1) / 2 - 1;

    // |m| from remainder: for n with parity p, valid |m| are p, p+2, p+4, ..., n
    // Each |m| > 0 gets two slots (sin, cos). |m|=0 gets one slot.
    let m_abs = if n.is_multiple_of(2) {
        // even n: slot 0 → m=0, slots (1,2) → m=2, slots (3,4) → m=4, ...
        if k == 0 { 0 } else { k.div_ceil(2) * 2 }
    } else {
        // odd n: slots (0,1) → m=1, slots (2,3) → m=3, ...
        (k / 2) * 2 + 1
    };

    // Sign: even j → positive m (cosine), odd j → negative m (sine)
    let m = if m_abs == 0 {
        0i32
    } else if j.is_multiple_of(2) {
        m_abs as i32
    } else {
        -(m_abs as i32)
    };
    (n, m)
}

/// Convert (n, m) indices to Noll sequential index j (1-based).
#[must_use]
pub fn nm_to_noll(n: u32, m: i32) -> u32 {
    let m_abs = m.unsigned_abs();

    // Base of this radial order (j of first term)
    let base = n * (n + 1) / 2 + 1;

    // Find 0-based slot for |m|
    let slot_start = if n.is_multiple_of(2) {
        if m_abs == 0 { 0 } else { 2 * (m_abs / 2) - 1 }
    } else {
        2 * ((m_abs - 1) / 2)
    };

    if m_abs == 0 {
        return base;
    }

    // Two slots per |m|: odd j gets -m, even j gets +m
    let j_first = base + slot_start;
    if m > 0 {
        // Positive m → even j
        if j_first.is_multiple_of(2) {
            j_first
        } else {
            j_first + 1
        }
    } else {
        // Negative m → odd j
        if j_first % 2 == 1 {
            j_first
        } else {
            j_first + 1
        }
    }
}

// ── Wavefront from Zernike coefficients ─────────────────────────────────────

/// Zernike wavefront expansion — a set of coefficients defining a wavefront.
///
/// Coefficients are stored in Noll order (j = 1, 2, 3, ...).
/// The wavefront at any point is the sum of Zernike terms:
///
/// W(ρ, θ) = Σ_j c_j · Z_j(ρ, θ)
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ZernikeWavefront {
    /// Coefficients in Noll order (index 0 = j=1 = piston).
    pub coefficients: Vec<f64>,
}

impl ZernikeWavefront {
    /// Create a wavefront from Noll-ordered coefficients.
    ///
    /// `coeffs[0]` = piston (j=1), `coeffs[1]` = tilt Y (j=2), etc.
    #[must_use]
    pub fn new(coefficients: Vec<f64>) -> Self {
        Self { coefficients }
    }

    /// Create a flat (zero-aberration) wavefront.
    #[must_use]
    pub fn flat() -> Self {
        Self {
            coefficients: Vec::new(),
        }
    }

    /// Evaluate the wavefront at a point on the unit pupil.
    ///
    /// `rho` in [0, 1], `theta` in radians. Returns 0.0 outside the pupil.
    #[must_use]
    #[inline]
    pub fn evaluate(&self, rho: f64, theta: f64) -> f64 {
        if rho > 1.0 {
            return 0.0;
        }
        let mut w = 0.0;
        for (idx, &coeff) in self.coefficients.iter().enumerate() {
            if coeff.abs() < 1e-15 {
                continue; // skip zero coefficients
            }
            let j = (idx + 1) as u32;
            let (n, m) = noll_to_nm(j);
            w += coeff * zernike(n, m, rho, theta);
        }
        w
    }

    /// Compute the wavefront on a grid of points within the unit pupil.
    ///
    /// Returns a flat array of `grid_size × grid_size` values. Points outside
    /// the unit disk are set to 0.0.
    #[must_use]
    pub fn to_grid(&self, grid_size: usize) -> Vec<f64> {
        let mut grid = vec![0.0; grid_size * grid_size];
        let center = (grid_size as f64 - 1.0) / 2.0;
        for iy in 0..grid_size {
            for ix in 0..grid_size {
                let x = (ix as f64 - center) / center;
                let y = (iy as f64 - center) / center;
                let rho = (x * x + y * y).sqrt();
                if rho <= 1.0 {
                    let theta = y.atan2(x);
                    grid[iy * grid_size + ix] = self.evaluate(rho, theta);
                }
            }
        }
        grid
    }

    /// RMS wavefront error (excluding piston, tilt X, and tilt Y).
    ///
    /// Due to Zernike orthogonality: RMS² = Σ c_j² (for j > 3).
    #[must_use]
    pub fn rms_error(&self) -> f64 {
        // Skip j=1 (piston), j=2 (tilt Y), j=3 (tilt X)
        self.coefficients
            .iter()
            .skip(3)
            .map(|c| c * c)
            .sum::<f64>()
            .sqrt()
    }

    /// Peak-to-valley wavefront error on the unit pupil.
    ///
    /// Sampled on a grid. Higher `samples` gives better accuracy.
    #[must_use]
    pub fn peak_to_valley(&self, samples: usize) -> f64 {
        let mut min = f64::INFINITY;
        let mut max = f64::NEG_INFINITY;
        let center = (samples as f64 - 1.0) / 2.0;
        for iy in 0..samples {
            for ix in 0..samples {
                let x = (ix as f64 - center) / center;
                let y = (iy as f64 - center) / center;
                let rho = (x * x + y * y).sqrt();
                if rho <= 1.0 {
                    let theta = y.atan2(x);
                    let w = self.evaluate(rho, theta);
                    if w < min {
                        min = w;
                    }
                    if w > max {
                        max = w;
                    }
                }
            }
        }
        if min.is_finite() && max.is_finite() {
            max - min
        } else {
            0.0
        }
    }

    /// Strehl ratio estimate from RMS wavefront error (Marechal approximation).
    ///
    /// S ≈ exp(−(2π·σ)²) where σ is the RMS wavefront error in waves.
    ///
    /// `rms_waves` = RMS error expressed in units of wavelength.
    /// Valid for σ < 0.1 waves (S > 0.6).
    #[must_use]
    #[inline]
    pub fn strehl_ratio(&self, wavelength: f64) -> f64 {
        let sigma = self.rms_error() / wavelength;
        let x = 2.0 * PI * sigma;
        (-x * x).exp()
    }
}

// ── Named aberrations ───────────────────────────────────────────────────────

/// Create a Zernike wavefront with a single defocus term.
///
/// `amount` is the coefficient in length units (same as wavelength for Strehl).
#[must_use]
pub fn defocus(amount: f64) -> ZernikeWavefront {
    let mut coeffs = vec![0.0; 4];
    coeffs[3] = amount; // j=4 = defocus
    ZernikeWavefront::new(coeffs)
}

/// Create a Zernike wavefront with primary spherical aberration.
///
/// Noll j=11, (n=4, m=0).
#[must_use]
pub fn spherical(amount: f64) -> ZernikeWavefront {
    let mut coeffs = vec![0.0; 11];
    coeffs[10] = amount; // j=11
    ZernikeWavefront::new(coeffs)
}

/// Create a Zernike wavefront with primary coma (X).
///
/// Noll j=8, (n=3, m=1).
#[must_use]
pub fn coma_x(amount: f64) -> ZernikeWavefront {
    let mut coeffs = vec![0.0; 8];
    coeffs[7] = amount; // j=8
    ZernikeWavefront::new(coeffs)
}

/// Create a Zernike wavefront with primary coma (Y).
///
/// Noll j=7, (n=3, m=-1).
#[must_use]
pub fn coma_y(amount: f64) -> ZernikeWavefront {
    let mut coeffs = vec![0.0; 7];
    coeffs[6] = amount; // j=7
    ZernikeWavefront::new(coeffs)
}

/// Create a Zernike wavefront with astigmatism (0°).
///
/// Noll j=6, (n=2, m=2).
#[must_use]
pub fn astigmatism_0(amount: f64) -> ZernikeWavefront {
    let mut coeffs = vec![0.0; 6];
    coeffs[5] = amount; // j=6
    ZernikeWavefront::new(coeffs)
}

/// Create a Zernike wavefront with astigmatism (45°).
///
/// Noll j=5, (n=2, m=-2).
#[must_use]
pub fn astigmatism_45(amount: f64) -> ZernikeWavefront {
    let mut coeffs = vec![0.0; 5];
    coeffs[4] = amount; // j=5
    ZernikeWavefront::new(coeffs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_PI_4;

    const EPS: f64 = 1e-6;

    // ── Radial polynomial tests ──────────────────────────────────────────

    #[test]
    fn test_radial_piston() {
        // R_0^0(ρ) = 1 for all ρ
        assert!((radial_polynomial(0, 0, 0.0) - 1.0).abs() < EPS);
        assert!((radial_polynomial(0, 0, 0.5) - 1.0).abs() < EPS);
        assert!((radial_polynomial(0, 0, 1.0) - 1.0).abs() < EPS);
    }

    #[test]
    fn test_radial_tilt() {
        // R_1^1(ρ) = ρ
        assert!((radial_polynomial(1, 1, 0.0)).abs() < EPS);
        assert!((radial_polynomial(1, 1, 0.5) - 0.5).abs() < EPS);
        assert!((radial_polynomial(1, 1, 1.0) - 1.0).abs() < EPS);
    }

    #[test]
    fn test_radial_defocus() {
        // R_2^0(ρ) = 2ρ² - 1
        assert!((radial_polynomial(2, 0, 0.0) + 1.0).abs() < EPS);
        assert!((radial_polynomial(2, 0, 1.0) - 1.0).abs() < EPS);
        let mid = 1.0 / 2.0f64.sqrt(); // ρ where 2ρ²-1 = 0
        assert!((radial_polynomial(2, 0, mid)).abs() < 0.01);
    }

    #[test]
    fn test_radial_spherical() {
        // R_4^0(ρ) = 6ρ⁴ - 6ρ² + 1
        assert!((radial_polynomial(4, 0, 0.0) - 1.0).abs() < EPS);
        assert!((radial_polynomial(4, 0, 1.0) - 1.0).abs() < EPS);
    }

    #[test]
    fn test_radial_invalid_nm() {
        // n < m should return 0
        assert!((radial_polynomial(1, 2, 0.5)).abs() < EPS);
        // n-m odd should return 0
        assert!((radial_polynomial(3, 2, 0.5)).abs() < EPS);
    }

    // ── Full Zernike polynomial tests ────────────────────────────────────

    #[test]
    fn test_zernike_piston() {
        // Z(0,0) = √1 · 1 = 1
        let z = zernike(0, 0, 0.5, 0.0);
        assert!((z - 1.0).abs() < EPS);
    }

    #[test]
    fn test_zernike_tilt_x() {
        // Z(1,1) = √4 · ρ · cos(θ) = 2ρ cos(θ)
        let z = zernike(1, 1, 1.0, 0.0);
        assert!((z - 2.0).abs() < EPS); // ρ=1, cos(0)=1, norm=√4=2
    }

    #[test]
    fn test_zernike_tilt_y() {
        // Z(1,-1) = √4 · ρ · sin(θ) = 2ρ sin(θ)
        let z = zernike(1, -1, 1.0, PI / 2.0);
        assert!((z - 2.0).abs() < EPS); // ρ=1, sin(π/2)=1
    }

    #[test]
    fn test_zernike_defocus() {
        // Z(2,0) = √3 · (2ρ²-1)
        let z_center = zernike(2, 0, 0.0, 0.0);
        assert!((z_center + 3.0f64.sqrt()).abs() < EPS); // -√3
        let z_edge = zernike(2, 0, 1.0, 0.0);
        assert!((z_edge - 3.0f64.sqrt()).abs() < EPS); // +√3
    }

    #[test]
    fn test_zernike_outside_pupil() {
        assert!((zernike(2, 0, 1.5, 0.0)).abs() < EPS);
    }

    #[test]
    fn test_zernike_orthogonality_low_order() {
        // Numerical test: integral of Z_j * Z_k over unit disk ≈ π·δ_jk
        // Use grid integration for j=1 (piston) and j=4 (defocus)
        let n = 100;
        let mut sum_11 = 0.0;
        let mut sum_14 = 0.0;
        let mut sum_44 = 0.0;
        let da = 4.0 / (n * n) as f64;
        for iy in 0..n {
            for ix in 0..n {
                let x = (ix as f64 + 0.5) / (n as f64 / 2.0) - 1.0;
                let y = (iy as f64 + 0.5) / (n as f64 / 2.0) - 1.0;
                let rho = (x * x + y * y).sqrt();
                if rho <= 1.0 {
                    let theta = y.atan2(x);
                    let z1 = zernike(0, 0, rho, theta);
                    let z4 = zernike(2, 0, rho, theta);
                    sum_11 += z1 * z1 * da;
                    sum_14 += z1 * z4 * da;
                    sum_44 += z4 * z4 * da;
                }
            }
        }
        // Integral of Z_j^2 over unit disk = π for normalized Zernike
        assert!(
            (sum_11 - PI).abs() < 0.1,
            "Z1·Z1 should integrate to π, got {sum_11}"
        );
        assert!(
            sum_14.abs() < 0.1,
            "Z1·Z4 should integrate to 0, got {sum_14}"
        );
        assert!(
            (sum_44 - PI).abs() < 0.1,
            "Z4·Z4 should integrate to π, got {sum_44}"
        );
    }

    // ── Noll index tests ─────────────────────────────────────────────────

    #[test]
    fn test_noll_first_terms() {
        // Noll 1976: even j → cos (positive m), odd j → sin (negative m)
        assert_eq!(noll_to_nm(1), (0, 0)); // piston
        assert_eq!(noll_to_nm(2), (1, 1)); // tilt X (cos)
        assert_eq!(noll_to_nm(3), (1, -1)); // tilt Y (sin)
        assert_eq!(noll_to_nm(4), (2, 0)); // defocus
        assert_eq!(noll_to_nm(5), (2, -2)); // astigmatism 45°
        assert_eq!(noll_to_nm(6), (2, 2)); // astigmatism 0°
        assert_eq!(noll_to_nm(7), (3, -1)); // coma Y
        assert_eq!(noll_to_nm(8), (3, 1)); // coma X
        assert_eq!(noll_to_nm(11), (4, 0)); // spherical
    }

    #[test]
    fn test_noll_roundtrip() {
        for j in 1..=11 {
            let (n, m) = noll_to_nm(j);
            let back = nm_to_noll(n, m);
            assert_eq!(back, j, "Noll roundtrip failed: j={j} → ({n},{m}) → {back}");
        }
    }

    // ── Wavefront tests ──────────────────────────────────────────────────

    #[test]
    fn test_wavefront_flat() {
        let w = ZernikeWavefront::flat();
        assert!((w.evaluate(0.5, 0.0)).abs() < EPS);
        assert!((w.rms_error()).abs() < EPS);
    }

    #[test]
    fn test_wavefront_defocus() {
        let w = defocus(0.5);
        // Defocus: Z(2,0) = √3·(2ρ²-1)
        let center = w.evaluate(0.0, 0.0);
        let edge = w.evaluate(1.0, 0.0);
        assert!(center < 0.0, "Defocus should be negative at center");
        assert!(edge > 0.0, "Defocus should be positive at edge");
    }

    #[test]
    fn test_wavefront_rms() {
        // RMS of a single Zernike term is just |coefficient| (orthogonality)
        let w = defocus(0.1);
        let rms = w.rms_error();
        assert!(
            (rms - 0.1).abs() < EPS,
            "RMS of defocus 0.1 should be 0.1, got {rms}"
        );
    }

    #[test]
    fn test_wavefront_rms_excludes_piston_tilt() {
        // Piston and tilt should not contribute to RMS
        let mut coeffs = vec![1.0, 0.5, 0.5, 0.0]; // piston=1, tilt_y=0.5, tilt_x=0.5, defocus=0
        let w = ZernikeWavefront::new(coeffs.clone());
        assert!(
            w.rms_error().abs() < EPS,
            "Piston+tilt should have zero RMS"
        );
        // Add defocus
        coeffs.push(0.0); // j=5 = astigmatism
        coeffs[3] = 0.2; // j=4 = defocus
        let w2 = ZernikeWavefront::new(coeffs);
        assert!((w2.rms_error() - 0.2).abs() < EPS);
    }

    #[test]
    fn test_wavefront_strehl_diffraction_limited() {
        let w = ZernikeWavefront::flat();
        let s = w.strehl_ratio(550e-9);
        assert!((s - 1.0).abs() < EPS, "Flat wavefront should have Strehl=1");
    }

    #[test]
    fn test_wavefront_strehl_marechal() {
        // Marechal criterion: λ/14 RMS → Strehl ≈ 0.8
        let wl = 550e-9;
        let rms = wl / 14.0;
        let mut coeffs = vec![0.0; 11];
        coeffs[10] = rms; // spherical
        let w = ZernikeWavefront::new(coeffs);
        let s = w.strehl_ratio(wl);
        assert!(
            (s - 0.8).abs() < 0.05,
            "Marechal criterion: Strehl ≈ 0.8, got {s}"
        );
    }

    #[test]
    fn test_wavefront_to_grid() {
        let w = defocus(0.1);
        let grid = w.to_grid(32);
        assert_eq!(grid.len(), 32 * 32);
        // Center should be negative (defocus), corners should be zero (outside pupil)
        let center = grid[16 * 32 + 16];
        let corner = grid[0]; // outside unit disk
        assert!(center != 0.0, "Center should have nonzero wavefront");
        assert!(corner.abs() < EPS, "Corner should be zero (outside pupil)");
    }

    #[test]
    fn test_wavefront_pv() {
        let w = defocus(0.5);
        let pv = w.peak_to_valley(64);
        assert!(pv > 0.0, "P-V should be positive");
        // For defocus: PV ≈ 2·|coeff|·√3 (range of Z4 from -√3 to +√3)
        let expected = 2.0 * 0.5 * 3.0f64.sqrt();
        assert!(
            (pv - expected).abs() < 0.1,
            "Defocus PV ≈ {expected}, got {pv}"
        );
    }

    #[test]
    fn test_named_aberrations() {
        // Each named constructor should produce a wavefront with exactly one nonzero coefficient
        let aberrations = [
            defocus(1.0),
            spherical(1.0),
            coma_x(1.0),
            coma_y(1.0),
            astigmatism_0(1.0),
            astigmatism_45(1.0),
        ];
        for w in &aberrations {
            let nonzero: usize = w.coefficients.iter().filter(|c| c.abs() > EPS).count();
            assert_eq!(
                nonzero, 1,
                "Named aberration should have exactly 1 nonzero coefficient"
            );
        }
    }

    #[test]
    fn test_coma_symmetry() {
        // Coma X should be symmetric about x-axis
        let w = coma_x(0.5);
        let v1 = w.evaluate(0.7, FRAC_PI_4);
        let v2 = w.evaluate(0.7, -FRAC_PI_4);
        assert!((v1 - v2).abs() < EPS, "Coma X should be x-axis symmetric");
    }
}
