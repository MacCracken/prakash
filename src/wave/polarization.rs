//! Stokes parameters, Mueller matrices, and birefringence.

use serde::{Deserialize, Serialize};

/// Stokes vector representing a polarization state.
///
/// - S0: total intensity
/// - S1: preference for horizontal (>0) vs vertical (<0) polarization
/// - S2: preference for +45° (>0) vs −45° (<0) polarization
/// - S3: preference for right circular (>0) vs left circular (<0)
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct StokesVector {
    pub s0: f64,
    pub s1: f64,
    pub s2: f64,
    pub s3: f64,
}

impl StokesVector {
    #[must_use]
    #[inline]
    pub const fn new(s0: f64, s1: f64, s2: f64, s3: f64) -> Self {
        Self { s0, s1, s2, s3 }
    }

    /// Unpolarized light of given intensity.
    #[must_use]
    #[inline]
    pub const fn unpolarized(intensity: f64) -> Self {
        Self::new(intensity, 0.0, 0.0, 0.0)
    }

    /// Horizontally polarized light.
    #[must_use]
    #[inline]
    pub const fn horizontal(intensity: f64) -> Self {
        Self::new(intensity, intensity, 0.0, 0.0)
    }

    /// Vertically polarized light.
    #[must_use]
    #[inline]
    pub const fn vertical(intensity: f64) -> Self {
        Self::new(intensity, -intensity, 0.0, 0.0)
    }

    /// +45° linearly polarized light.
    #[must_use]
    #[inline]
    pub const fn diagonal_plus(intensity: f64) -> Self {
        Self::new(intensity, 0.0, intensity, 0.0)
    }

    /// −45° linearly polarized light.
    #[must_use]
    #[inline]
    pub const fn diagonal_minus(intensity: f64) -> Self {
        Self::new(intensity, 0.0, -intensity, 0.0)
    }

    /// Right circular polarization.
    #[must_use]
    #[inline]
    pub const fn circular_right(intensity: f64) -> Self {
        Self::new(intensity, 0.0, 0.0, intensity)
    }

    /// Left circular polarization.
    #[must_use]
    #[inline]
    pub const fn circular_left(intensity: f64) -> Self {
        Self::new(intensity, 0.0, 0.0, -intensity)
    }

    /// Degree of polarization (0.0 = unpolarized, 1.0 = fully polarized).
    #[must_use]
    #[inline]
    pub fn degree_of_polarization(&self) -> f64 {
        if self.s0.abs() < 1e-15 {
            return 0.0;
        }
        (self.s1 * self.s1 + self.s2 * self.s2 + self.s3 * self.s3).sqrt() / self.s0
    }

    /// Total intensity.
    #[must_use]
    #[inline]
    pub fn intensity(&self) -> f64 {
        self.s0
    }

    /// Ellipticity angle χ: tan(2χ) = S3 / √(S1² + S2²).
    /// χ = 0 for linear, ±π/4 for circular.
    #[must_use]
    #[inline]
    pub fn ellipticity_angle(&self) -> f64 {
        let linear_part = (self.s1 * self.s1 + self.s2 * self.s2).sqrt();
        if linear_part < 1e-15 && self.s3.abs() < 1e-15 {
            return 0.0;
        }
        0.5 * self.s3.atan2(linear_part)
    }

    /// Orientation angle ψ of the polarization ellipse: tan(2ψ) = S2/S1.
    /// Range: [−π/2, π/2].
    #[must_use]
    #[inline]
    pub fn orientation_angle(&self) -> f64 {
        0.5 * self.s2.atan2(self.s1)
    }
}

// ── Mueller Matrices ──────────────────────────────────────────────────────

/// A 4×4 Mueller matrix for transforming Stokes vectors.
///
/// Stored in row-major order: `m[row][col]`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MuellerMatrix {
    pub m: [[f64; 4]; 4],
}

impl MuellerMatrix {
    /// Create a Mueller matrix from a 4×4 array (row-major).
    #[must_use]
    #[inline]
    pub const fn new(m: [[f64; 4]; 4]) -> Self {
        Self { m }
    }

    /// Identity matrix (no effect on polarization).
    pub const IDENTITY: Self = Self::new([
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]);

    /// Ideal horizontal linear polarizer.
    pub const POLARIZER_HORIZONTAL: Self = Self::new([
        [0.5, 0.5, 0.0, 0.0],
        [0.5, 0.5, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
    ]);

    /// Ideal vertical linear polarizer.
    pub const POLARIZER_VERTICAL: Self = Self::new([
        [0.5, -0.5, 0.0, 0.0],
        [-0.5, 0.5, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
    ]);

    /// Linear polarizer at arbitrary angle θ (radians).
    #[must_use]
    #[inline]
    pub fn polarizer(angle: f64) -> Self {
        let (s2, c2) = (2.0 * angle).sin_cos();
        let cs = c2 * s2;
        Self::new([
            [0.5, 0.5 * c2, 0.5 * s2, 0.0],
            [0.5 * c2, 0.5 * c2 * c2, 0.5 * cs, 0.0],
            [0.5 * s2, 0.5 * cs, 0.5 * s2 * s2, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ])
    }

    /// Quarter-wave plate with fast axis horizontal.
    pub const QUARTER_WAVE_HORIZONTAL: Self = Self::new([
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0, 0.0],
    ]);

    /// Half-wave plate with fast axis horizontal.
    pub const HALF_WAVE_HORIZONTAL: Self = Self::new([
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, -1.0, 0.0],
        [0.0, 0.0, 0.0, -1.0],
    ]);

    /// General linear retarder with fast axis horizontal and retardance δ (radians).
    #[must_use]
    #[inline]
    pub fn retarder(retardance: f64) -> Self {
        let c = retardance.cos();
        let s = retardance.sin();
        Self::new([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, c, s],
            [0.0, 0.0, -s, c],
        ])
    }

    /// Rotation matrix: rotates the reference frame by angle θ.
    /// Used to orient optical elements at arbitrary angles.
    #[must_use]
    #[inline]
    pub fn rotation(angle: f64) -> Self {
        let (s2, c2) = (2.0 * angle).sin_cos();
        Self::new([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, c2, s2, 0.0],
            [0.0, -s2, c2, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ])
    }

    /// Apply this Mueller matrix to a Stokes vector.
    #[must_use]
    #[inline]
    pub fn apply(&self, s: &StokesVector) -> StokesVector {
        let sv = [s.s0, s.s1, s.s2, s.s3];
        StokesVector::new(
            self.m[0][0] * sv[0]
                + self.m[0][1] * sv[1]
                + self.m[0][2] * sv[2]
                + self.m[0][3] * sv[3],
            self.m[1][0] * sv[0]
                + self.m[1][1] * sv[1]
                + self.m[1][2] * sv[2]
                + self.m[1][3] * sv[3],
            self.m[2][0] * sv[0]
                + self.m[2][1] * sv[1]
                + self.m[2][2] * sv[2]
                + self.m[2][3] * sv[3],
            self.m[3][0] * sv[0]
                + self.m[3][1] * sv[1]
                + self.m[3][2] * sv[2]
                + self.m[3][3] * sv[3],
        )
    }

    /// Multiply two Mueller matrices: self · other.
    /// Result represents applying `other` first, then `self`.
    #[must_use]
    #[inline]
    pub fn multiply(&self, other: &MuellerMatrix) -> MuellerMatrix {
        let mut result = [[0.0; 4]; 4];
        for (i, row) in result.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                *cell = self.m[i][0] * other.m[0][j]
                    + self.m[i][1] * other.m[1][j]
                    + self.m[i][2] * other.m[2][j]
                    + self.m[i][3] * other.m[3][j];
            }
        }
        MuellerMatrix::new(result)
    }
}

/// Apply a chain of Mueller matrices to a Stokes vector (left to right = first to last).
#[must_use]
#[inline]
pub fn mueller_chain(stokes: &StokesVector, elements: &[MuellerMatrix]) -> StokesVector {
    let mut s = *stokes;
    for m in elements {
        s = m.apply(&s);
    }
    s
}

// ── Birefringence ─────────────────────────────────────────────────────────

/// Birefringent material properties.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BirefringentMaterial {
    /// Ordinary refractive index.
    pub n_o: f64,
    /// Extraordinary refractive index.
    pub n_e: f64,
    /// Human-readable name.
    pub name: &'static str,
}

impl BirefringentMaterial {
    /// Birefringence: Δn = n_e − n_o.
    /// Positive = positive uniaxial, negative = negative uniaxial.
    #[must_use]
    #[inline]
    pub fn birefringence(&self) -> f64 {
        self.n_e - self.n_o
    }

    /// Phase retardation for a given thickness and wavelength.
    ///
    /// δ = 2π · d · |n_e − n_o| / λ (radians)
    ///
    /// `thickness` and `wavelength` in same units.
    #[must_use]
    #[inline]
    pub fn retardation(&self, thickness: f64, wavelength: f64) -> f64 {
        std::f64::consts::TAU * thickness * self.birefringence().abs() / wavelength
    }

    /// Thickness needed for a quarter-wave plate at a given wavelength.
    ///
    /// d = λ / (4·|Δn|)
    #[must_use]
    #[inline]
    pub fn quarter_wave_thickness(&self, wavelength: f64) -> f64 {
        wavelength / (4.0 * self.birefringence().abs())
    }

    /// Thickness needed for a half-wave plate at a given wavelength.
    ///
    /// d = λ / (2·|Δn|)
    #[must_use]
    #[inline]
    pub fn half_wave_thickness(&self, wavelength: f64) -> f64 {
        wavelength / (2.0 * self.birefringence().abs())
    }

    /// Mueller matrix for this birefringent material at a given thickness and wavelength.
    /// Fast axis horizontal.
    #[must_use]
    #[inline]
    pub fn to_mueller(&self, thickness: f64, wavelength: f64) -> MuellerMatrix {
        MuellerMatrix::retarder(self.retardation(thickness, wavelength))
    }

    /// Calcite (CaCO₃) — strong negative birefringence.
    pub const CALCITE: Self = Self {
        n_o: 1.6584,
        n_e: 1.4864,
        name: "calcite",
    };

    /// Quartz (SiO₂) — weak positive birefringence.
    pub const QUARTZ: Self = Self {
        n_o: 1.5443,
        n_e: 1.5534,
        name: "quartz",
    };

    /// Rutile (TiO₂) — strong positive birefringence.
    pub const RUTILE: Self = Self {
        n_o: 2.616,
        n_e: 2.903,
        name: "rutile",
    };

    /// Mica (muscovite) — moderate birefringence.
    pub const MICA: Self = Self {
        n_o: 1.5601,
        n_e: 1.5936,
        name: "mica",
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    const EPS: f64 = 1e-6;

    #[test]
    fn test_stokes_unpolarized() {
        let s = StokesVector::unpolarized(1.0);
        assert!((s.s0 - 1.0).abs() < EPS);
        assert!(s.s1.abs() < EPS);
        assert!(s.s2.abs() < EPS);
        assert!(s.s3.abs() < EPS);
        assert!(s.degree_of_polarization() < EPS);
    }

    #[test]
    fn test_stokes_horizontal() {
        let s = StokesVector::horizontal(1.0);
        assert!((s.degree_of_polarization() - 1.0).abs() < EPS);
        assert!((s.s1 - 1.0).abs() < EPS);
    }

    #[test]
    fn test_stokes_vertical() {
        let s = StokesVector::vertical(1.0);
        assert!((s.degree_of_polarization() - 1.0).abs() < EPS);
        assert!((s.s1 + 1.0).abs() < EPS);
    }

    #[test]
    fn test_stokes_circular() {
        let r = StokesVector::circular_right(1.0);
        assert!((r.degree_of_polarization() - 1.0).abs() < EPS);
        assert!((r.s3 - 1.0).abs() < EPS);

        let l = StokesVector::circular_left(1.0);
        assert!((l.s3 + 1.0).abs() < EPS);
    }

    #[test]
    fn test_stokes_dop_partial() {
        let s = StokesVector::new(1.0, 0.5, 0.0, 0.0);
        let dop = s.degree_of_polarization();
        assert!((dop - 0.5).abs() < EPS);
    }

    #[test]
    fn test_stokes_ellipticity_linear() {
        let s = StokesVector::horizontal(1.0);
        assert!(s.ellipticity_angle().abs() < EPS);
    }

    #[test]
    fn test_stokes_ellipticity_circular() {
        let s = StokesVector::circular_right(1.0);
        assert!((s.ellipticity_angle() - PI / 4.0).abs() < 0.01);
    }

    #[test]
    fn test_stokes_orientation_horizontal() {
        let s = StokesVector::horizontal(1.0);
        assert!(s.orientation_angle().abs() < EPS);
    }

    #[test]
    fn test_stokes_orientation_45() {
        let s = StokesVector::diagonal_plus(1.0);
        assert!((s.orientation_angle() - PI / 4.0).abs() < 0.01);
    }

    #[test]
    fn test_stokes_zero_intensity() {
        let s = StokesVector::new(0.0, 0.0, 0.0, 0.0);
        assert!((s.degree_of_polarization()).abs() < EPS);
    }

    // ── Mueller matrix tests ──────────────────────────────────────────────

    #[test]
    fn test_mueller_identity() {
        let s = StokesVector::horizontal(1.0);
        let result = MuellerMatrix::IDENTITY.apply(&s);
        assert!((result.s0 - s.s0).abs() < EPS);
        assert!((result.s1 - s.s1).abs() < EPS);
        assert!((result.s2 - s.s2).abs() < EPS);
        assert!((result.s3 - s.s3).abs() < EPS);
    }

    #[test]
    fn test_mueller_horizontal_polarizer_passes_horizontal() {
        let s = StokesVector::horizontal(1.0);
        let result = MuellerMatrix::POLARIZER_HORIZONTAL.apply(&s);
        assert!((result.s0 - 1.0).abs() < EPS);
    }

    #[test]
    fn test_mueller_horizontal_polarizer_blocks_vertical() {
        let s = StokesVector::vertical(1.0);
        let result = MuellerMatrix::POLARIZER_HORIZONTAL.apply(&s);
        assert!(result.s0.abs() < EPS);
    }

    #[test]
    fn test_mueller_horizontal_polarizer_halves_unpolarized() {
        let s = StokesVector::unpolarized(1.0);
        let result = MuellerMatrix::POLARIZER_HORIZONTAL.apply(&s);
        assert!((result.s0 - 0.5).abs() < EPS);
    }

    #[test]
    fn test_mueller_crossed_polarizers_zero() {
        let s = StokesVector::unpolarized(1.0);
        let after_h = MuellerMatrix::POLARIZER_HORIZONTAL.apply(&s);
        let after_v = MuellerMatrix::POLARIZER_VERTICAL.apply(&after_h);
        assert!(after_v.s0.abs() < EPS);
    }

    #[test]
    fn test_mueller_polarizer_arbitrary_angle() {
        // Polarizer at 0° should match horizontal
        let p0 = MuellerMatrix::polarizer(0.0);
        let ph = MuellerMatrix::POLARIZER_HORIZONTAL;
        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    (p0.m[i][j] - ph.m[i][j]).abs() < EPS,
                    "Mismatch at [{i}][{j}]"
                );
            }
        }
    }

    #[test]
    fn test_mueller_quarter_wave_h_to_circular() {
        // Horizontal through QWP → right circular
        let s = StokesVector::diagonal_plus(1.0);
        let result = MuellerMatrix::QUARTER_WAVE_HORIZONTAL.apply(&s);
        assert!(
            result.s3.abs() > 0.9,
            "Should produce circular: s3={}",
            result.s3
        );
    }

    #[test]
    fn test_mueller_half_wave_flips_s2_s3() {
        let s = StokesVector::new(1.0, 0.0, 0.5, 0.5);
        let result = MuellerMatrix::HALF_WAVE_HORIZONTAL.apply(&s);
        assert!((result.s0 - 1.0).abs() < EPS);
        assert!((result.s2 + 0.5).abs() < EPS); // flipped
        assert!((result.s3 + 0.5).abs() < EPS); // flipped
    }

    #[test]
    fn test_mueller_retarder_zero_is_identity() {
        let r = MuellerMatrix::retarder(0.0);
        for i in 0..4 {
            for j in 0..4 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (r.m[i][j] - expected).abs() < EPS,
                    "Retarder(0) should be identity at [{i}][{j}]"
                );
            }
        }
    }

    #[test]
    fn test_mueller_retarder_pi_is_half_wave() {
        let r = MuellerMatrix::retarder(PI);
        let hw = MuellerMatrix::HALF_WAVE_HORIZONTAL;
        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    (r.m[i][j] - hw.m[i][j]).abs() < EPS,
                    "Retarder(π) should match HWP at [{i}][{j}]"
                );
            }
        }
    }

    #[test]
    fn test_mueller_retarder_halfpi_is_quarter_wave() {
        let r = MuellerMatrix::retarder(PI / 2.0);
        let qw = MuellerMatrix::QUARTER_WAVE_HORIZONTAL;
        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    (r.m[i][j] - qw.m[i][j]).abs() < EPS,
                    "Retarder(π/2) should match QWP at [{i}][{j}]"
                );
            }
        }
    }

    #[test]
    fn test_mueller_multiply_identity() {
        let m = MuellerMatrix::POLARIZER_HORIZONTAL;
        let result = MuellerMatrix::IDENTITY.multiply(&m);
        for i in 0..4 {
            for j in 0..4 {
                assert!((result.m[i][j] - m.m[i][j]).abs() < EPS);
            }
        }
    }

    #[test]
    fn test_mueller_chain_two_polarizers() {
        let s = StokesVector::unpolarized(1.0);
        let result = mueller_chain(
            &s,
            &[
                MuellerMatrix::POLARIZER_HORIZONTAL,
                MuellerMatrix::POLARIZER_VERTICAL,
            ],
        );
        assert!(
            result.s0.abs() < EPS,
            "Crossed polarizers should block all light"
        );
    }

    #[test]
    fn test_mueller_rotation_preserves_intensity() {
        let s = StokesVector::horizontal(1.0);
        let r = MuellerMatrix::rotation(PI / 6.0);
        let result = r.apply(&s);
        assert!(
            (result.s0 - 1.0).abs() < EPS,
            "Rotation should preserve intensity"
        );
    }

    // ── Birefringence tests ───────────────────────────────────────────────

    #[test]
    fn test_calcite_negative_birefringence() {
        assert!(
            BirefringentMaterial::CALCITE.birefringence() < 0.0,
            "Calcite is negative uniaxial"
        );
    }

    #[test]
    fn test_quartz_positive_birefringence() {
        assert!(
            BirefringentMaterial::QUARTZ.birefringence() > 0.0,
            "Quartz is positive uniaxial"
        );
    }

    #[test]
    fn test_birefringence_magnitude() {
        // Calcite has strong birefringence |Δn| ≈ 0.172
        let dn = BirefringentMaterial::CALCITE.birefringence().abs();
        assert!((dn - 0.172).abs() < 0.001);
    }

    #[test]
    fn test_retardation_quarter_wave() {
        let mat = BirefringentMaterial::QUARTZ;
        let wl = 550.0;
        let d = mat.quarter_wave_thickness(wl);
        let ret = mat.retardation(d, wl);
        assert!(
            (ret - PI / 2.0).abs() < 0.01,
            "QWP retardation should be π/2, got {ret}"
        );
    }

    #[test]
    fn test_retardation_half_wave() {
        let mat = BirefringentMaterial::QUARTZ;
        let wl = 550.0;
        let d = mat.half_wave_thickness(wl);
        let ret = mat.retardation(d, wl);
        assert!(
            (ret - PI).abs() < 0.01,
            "HWP retardation should be π, got {ret}"
        );
    }

    #[test]
    fn test_quarter_wave_half_wave_relationship() {
        let mat = BirefringentMaterial::CALCITE;
        let wl = 632.8;
        let d_qwp = mat.quarter_wave_thickness(wl);
        let d_hwp = mat.half_wave_thickness(wl);
        assert!(
            (d_hwp / d_qwp - 2.0).abs() < EPS,
            "HWP should be 2× QWP thickness"
        );
    }

    #[test]
    fn test_birefringent_mueller() {
        let mat = BirefringentMaterial::QUARTZ;
        let wl = 550.0;
        let d = mat.quarter_wave_thickness(wl);
        let m = mat.to_mueller(d, wl);
        // Should be approximately a quarter-wave plate
        let qw = MuellerMatrix::QUARTER_WAVE_HORIZONTAL;
        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    (m.m[i][j] - qw.m[i][j]).abs() < 0.01,
                    "QWP Mueller mismatch at [{i}][{j}]: {} vs {}",
                    m.m[i][j],
                    qw.m[i][j]
                );
            }
        }
    }

    #[test]
    fn test_all_presets_valid() {
        let presets = [
            BirefringentMaterial::CALCITE,
            BirefringentMaterial::QUARTZ,
            BirefringentMaterial::RUTILE,
            BirefringentMaterial::MICA,
        ];
        for mat in &presets {
            assert!(mat.n_o >= 1.0, "{} n_o invalid", mat.name);
            assert!(mat.n_e >= 1.0, "{} n_e invalid", mat.name);
            assert!(
                mat.birefringence().abs() > 0.001,
                "{} should have measurable birefringence",
                mat.name
            );
        }
    }
}
