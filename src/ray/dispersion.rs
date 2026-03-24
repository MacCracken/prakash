//! Dispersion models (Cauchy, Sellmeier) and prism optics.

use serde::{Deserialize, Serialize};

use crate::error::{PrakashError, Result};

/// Cauchy dispersion coefficients.
///
/// Models wavelength-dependent refractive index as:
///   n(λ) = b + c/λ²
///
/// Wavelength λ is in micrometers.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CauchyCoefficients {
    pub b: f64,
    pub c: f64,
}

impl CauchyCoefficients {
    /// Refractive index at a given wavelength (micrometers).
    ///
    /// Wavelength must be positive. Returns `b` coefficient for very large wavelengths.
    #[inline]
    pub fn n_at(&self, wavelength_um: f64) -> f64 {
        debug_assert!(wavelength_um > 0.0, "wavelength must be positive");
        self.b + self.c / (wavelength_um * wavelength_um)
    }

    /// Common glass (approximate BK7).
    pub const BK7: Self = Self {
        b: 1.5046,
        c: 0.004_20,
    };

    /// Fused silica (approximate).
    pub const FUSED_SILICA: Self = Self {
        b: 1.4580,
        c: 0.003_54,
    };
}

/// Sellmeier dispersion coefficients (3-term).
///
/// Models wavelength-dependent refractive index as:
///   n²(λ) = 1 + B₁λ²/(λ²−C₁) + B₂λ²/(λ²−C₂) + B₃λ²/(λ²−C₃)
///
/// Wavelength λ is in micrometers. Cᵢ values are in μm².
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SellmeierCoefficients {
    pub b1: f64,
    pub c1: f64,
    pub b2: f64,
    pub c2: f64,
    pub b3: f64,
    pub c3: f64,
}

impl SellmeierCoefficients {
    /// Refractive index at a given wavelength (micrometers).
    ///
    /// Returns n >= 1.0 for valid wavelengths away from resonance poles.
    /// Near resonances (λ² ≈ Cᵢ), results may be unphysical.
    #[inline]
    pub fn n_at(&self, wavelength_um: f64) -> f64 {
        let l2 = wavelength_um * wavelength_um;
        // Guard denominators against resonance poles (l² ≈ Cᵢ)
        let d1 = l2 - self.c1;
        let d2 = l2 - self.c2;
        let d3 = l2 - self.c3;
        let t1 = if d1.abs() > 1e-15 {
            self.b1 * l2 / d1
        } else {
            0.0
        };
        let t2 = if d2.abs() > 1e-15 {
            self.b2 * l2 / d2
        } else {
            0.0
        };
        let t3 = if d3.abs() > 1e-15 {
            self.b3 * l2 / d3
        } else {
            0.0
        };
        let n2 = 1.0 + t1 + t2 + t3;
        // Guard against negative n² (unphysical region between resonances)
        if n2 < 1.0 {
            return 1.0;
        }
        n2.sqrt()
    }

    /// Schott N-BK7 borosilicate crown glass.
    pub const BK7: Self = Self {
        b1: 1.039_612_12,
        c1: 0.006_000_698_67,
        b2: 0.231_792_344,
        c2: 0.020_017_914_4,
        b3: 1.010_469_45,
        c3: 103.560_653,
    };

    /// Schott N-SF11 dense flint glass.
    pub const SF11: Self = Self {
        b1: 1.737_596_26,
        c1: 0.013_188_707_0,
        b2: 0.313_747_346,
        c2: 0.062_306_814_2,
        b3: 1.898_781_01,
        c3: 155.236_29,
    };

    /// Fused silica (SiO₂).
    pub const FUSED_SILICA: Self = Self {
        b1: 0.696_166_3,
        c1: 0.004_679_148_2,
        b2: 0.407_942_6,
        c2: 0.013_512_063,
        b3: 0.897_479_4,
        c3: 97.934_002_5,
    };

    /// Sapphire (Al₂O₃, ordinary ray).
    pub const SAPPHIRE: Self = Self {
        b1: 1.431_349_3,
        c1: 0.005_279_25,
        b2: 0.650_547_13,
        c2: 0.014_218_26,
        b3: 5.341_482_2,
        c3: 325.017_83,
    };

    /// Water at 25°C (Daimon & Masumura, 2007).
    pub const WATER: Self = Self {
        b1: 0.567_019_72,
        c1: 0.005_085_50,
        b2: 0.172_629_26,
        c2: 0.018_180_00,
        b3: 0.020_624_60,
        c3: 0.026_250_00,
    };

    /// Diamond (C).
    pub const DIAMOND: Self = Self {
        b1: 4.335_8,
        c1: 0.010_6,
        b2: 0.306_0,
        c2: 0.017_5,
        b3: 0.0,
        c3: 1.0, // two-term fit
    };
}

/// Fraunhofer spectral lines used for Abbe number calculation (in micrometers).
pub const FRAUNHOFER_D: f64 = 0.587_56; // Helium d-line (yellow)
pub const FRAUNHOFER_F: f64 = 0.486_13; // Hydrogen F-line (blue)
pub const FRAUNHOFER_C: f64 = 0.656_27; // Hydrogen C-line (red)

/// Abbe number (constringence) — a measure of a material's dispersion.
///
/// V = (n_d − 1) / (n_F − n_C)
///
/// Higher V means lower dispersion. Crown glasses typically V > 50,
/// flint glasses V < 50.
#[inline]
pub fn abbe_number(sellmeier: &SellmeierCoefficients) -> f64 {
    let n_d = sellmeier.n_at(FRAUNHOFER_D);
    let n_f = sellmeier.n_at(FRAUNHOFER_F);
    let n_c = sellmeier.n_at(FRAUNHOFER_C);
    (n_d - 1.0) / (n_f - n_c)
}

// ── Prism ─────────────────────────────────────────────────────────────────

/// Minimum deviation angle for a prism.
///
/// For a prism with apex angle A and refractive index n:
///   δ_min = 2·arcsin(n·sin(A/2)) − A
///
/// `apex_angle` is in radians.
#[inline]
pub fn prism_deviation(apex_angle: f64, n: f64) -> Result<f64> {
    let sin_half_a = (apex_angle / 2.0).sin();
    let arg = n * sin_half_a;
    if arg.abs() > 1.0 {
        return Err(PrakashError::TotalInternalReflection {
            angle_deg: (apex_angle / 2.0).to_degrees(),
            critical_deg: (1.0 / n).asin().to_degrees(),
            n1: n,
            n2: 1.0,
        });
    }
    Ok(2.0 * arg.asin() - apex_angle)
}

/// Angular dispersion of a prism for a given wavelength.
///
/// Returns the deviation angle using the Sellmeier equation to determine
/// the refractive index at the specified wavelength.
///
/// `apex_angle` is in radians, `wavelength_um` is in micrometers.
#[inline]
pub fn prism_dispersion(
    apex_angle: f64,
    sellmeier: &SellmeierCoefficients,
    wavelength_um: f64,
) -> Result<f64> {
    let n = sellmeier.n_at(wavelength_um);
    prism_deviation(apex_angle, n)
}

/// Angular spread of a prism across a wavelength range.
///
/// Returns (deviation_short, deviation_long, angular_spread) in radians.
/// The angular spread is the difference in deviation between the shortest
/// and longest wavelengths.
pub fn prism_angular_spread(
    apex_angle: f64,
    sellmeier: &SellmeierCoefficients,
    wavelength_short_um: f64,
    wavelength_long_um: f64,
) -> Result<(f64, f64, f64)> {
    let dev_short = prism_dispersion(apex_angle, sellmeier, wavelength_short_um)?;
    let dev_long = prism_dispersion(apex_angle, sellmeier, wavelength_long_um)?;
    Ok((dev_short, dev_long, (dev_short - dev_long).abs()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ray::*;

    const EPS: f64 = 1e-6;

    fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    }
    fn len3(v: [f64; 3]) -> f64 {
        dot3(v, v).sqrt()
    }
    fn normalize3(v: [f64; 3]) -> [f64; 3] {
        let l = len3(v);
        [v[0] / l, v[1] / l, v[2] / l]
    }

    #[test]
    fn test_cauchy_bk7_visible_range() {
        let n_blue = CauchyCoefficients::BK7.n_at(0.486);
        let n_red = CauchyCoefficients::BK7.n_at(0.656);
        let n_d = CauchyCoefficients::BK7.n_at(0.5876);

        assert!(n_blue > n_d, "Blue should have higher n than yellow");
        assert!(n_d > n_red, "Yellow should have higher n than red");
        assert!((n_d - 1.517).abs() < 0.02, "BK7 n_d ≈ 1.517");
    }

    #[test]
    fn test_cauchy_normal_dispersion() {
        // Normal dispersion: dn/dλ < 0 (n decreases with wavelength)
        let c = CauchyCoefficients::BK7;
        let n1 = c.n_at(0.4);
        let n2 = c.n_at(0.5);
        let n3 = c.n_at(0.6);
        let n4 = c.n_at(0.7);
        assert!(n1 > n2);
        assert!(n2 > n3);
        assert!(n3 > n4);
    }

    #[test]
    fn test_cauchy_fused_silica() {
        let n = CauchyCoefficients::FUSED_SILICA.n_at(0.5876);
        assert!((n - 1.458).abs() < 0.02);
    }

    // ── Sellmeier dispersion tests ────────────────────────────────────────

    #[test]
    fn test_sellmeier_bk7_nd() {
        let n = SellmeierCoefficients::BK7.n_at(FRAUNHOFER_D);
        assert!(
            (n - 1.5168).abs() < 0.001,
            "BK7 n_d should be ≈1.5168, got {n}"
        );
    }

    #[test]
    fn test_sellmeier_bk7_normal_dispersion() {
        let s = SellmeierCoefficients::BK7;
        let n_f = s.n_at(FRAUNHOFER_F);
        let n_d = s.n_at(FRAUNHOFER_D);
        let n_c = s.n_at(FRAUNHOFER_C);
        assert!(n_f > n_d, "n_F > n_d (blue > yellow)");
        assert!(n_d > n_c, "n_d > n_C (yellow > red)");
    }

    #[test]
    fn test_sellmeier_sf11() {
        let n = SellmeierCoefficients::SF11.n_at(FRAUNHOFER_D);
        assert!(
            (n - 1.7847).abs() < 0.001,
            "SF11 n_d should be ≈1.7847, got {n}"
        );
    }

    #[test]
    fn test_sellmeier_fused_silica() {
        let n = SellmeierCoefficients::FUSED_SILICA.n_at(FRAUNHOFER_D);
        assert!(
            (n - 1.4585).abs() < 0.001,
            "Fused silica n_d should be ≈1.4585, got {n}"
        );
    }

    #[test]
    fn test_sellmeier_water() {
        let n = SellmeierCoefficients::WATER.n_at(FRAUNHOFER_D);
        assert!(
            (n - 1.333).abs() < 0.005,
            "Water n_d should be ≈1.333, got {n}"
        );
    }

    #[test]
    fn test_sellmeier_diamond() {
        let n = SellmeierCoefficients::DIAMOND.n_at(FRAUNHOFER_D);
        assert!(
            (n - 2.417).abs() < 0.02,
            "Diamond n_d should be ≈2.417, got {n}"
        );
    }

    #[test]
    fn test_sellmeier_all_presets_reasonable() {
        let presets = [
            ("BK7", SellmeierCoefficients::BK7),
            ("SF11", SellmeierCoefficients::SF11),
            ("Fused Silica", SellmeierCoefficients::FUSED_SILICA),
            ("Sapphire", SellmeierCoefficients::SAPPHIRE),
            ("Water", SellmeierCoefficients::WATER),
            ("Diamond", SellmeierCoefficients::DIAMOND),
        ];
        for (name, s) in &presets {
            let n = s.n_at(FRAUNHOFER_D);
            assert!(
                (1.0..3.0).contains(&n),
                "{name}: n_d = {n} out of reasonable range"
            );
        }
    }

    #[test]
    fn test_sellmeier_serde_roundtrip() {
        let s = SellmeierCoefficients::BK7;
        let json = serde_json::to_string(&s).unwrap();
        let back: SellmeierCoefficients = serde_json::from_str(&json).unwrap();
        assert!((back.b1 - s.b1).abs() < EPS);
        assert!((back.c1 - s.c1).abs() < EPS);
    }

    // ── Abbe number tests ─────────────────────────────────────────────────

    #[test]
    fn test_abbe_number_bk7() {
        let v = abbe_number(&SellmeierCoefficients::BK7);
        // BK7 Abbe number ≈ 64.17
        assert!(
            (v - 64.17).abs() < 1.0,
            "BK7 Abbe number should be ≈64.17, got {v}"
        );
    }

    #[test]
    fn test_abbe_number_sf11() {
        let v = abbe_number(&SellmeierCoefficients::SF11);
        // SF11 Abbe number ≈ 25.76 (high dispersion flint)
        assert!(
            (v - 25.76).abs() < 1.0,
            "SF11 Abbe number should be ≈25.76, got {v}"
        );
    }

    #[test]
    fn test_abbe_crown_vs_flint() {
        let v_bk7 = abbe_number(&SellmeierCoefficients::BK7);
        let v_sf11 = abbe_number(&SellmeierCoefficients::SF11);
        assert!(
            v_bk7 > v_sf11,
            "Crown glass should have higher Abbe number than flint"
        );
        assert!(v_bk7 > 50.0, "Crown glass V > 50");
        assert!(v_sf11 < 50.0, "Flint glass V < 50");
    }

    #[test]
    fn test_sellmeier_resonance_guard() {
        // At a resonance pole (l² = c1), n_at should not return NaN or panic
        let s = SellmeierCoefficients::BK7;
        let resonance_wl = s.c1.sqrt(); // wavelength where l² = c1
        let n = s.n_at(resonance_wl);
        assert!(n.is_finite(), "n_at resonance should be finite, got {n}");
        assert!(n >= 1.0, "n_at resonance should be >= 1.0, got {n}");
    }

    #[test]
    fn test_refract_3d_output_normalized() {
        // Verify refracted direction has unit length
        let dir = normalize3([0.4, 0.3, -0.8]);
        let normal = [0.0, 0.0, 1.0];
        let refracted = refract_3d(dir, normal, 1.0, 1.5).unwrap();
        let len = len3(refracted);
        assert!(
            (len - 1.0).abs() < 0.01,
            "Refracted direction should be normalized, got length {len}"
        );
    }

    #[test]
    fn test_snell_3d_consistency_with_scalar() {
        // snell_3d at a known angle should agree with scalar snell()
        let angle_i = deg_to_rad(30.0);
        let n1 = 1.0;
        let n2 = 1.5;

        // Scalar
        let angle_t_scalar = snell(n1, n2, angle_i).unwrap();

        // 3D: construct direction at 30° in xz plane
        let dir = [angle_i.sin(), 0.0, -angle_i.cos()];
        let normal = [0.0, 0.0, 1.0];
        let (refracted, _) = snell_3d(dir, normal, n1, n2).unwrap();
        let cos_t_3d = -refracted[2]; // dot with -normal
        let angle_t_3d = cos_t_3d.acos();

        assert!(
            (angle_t_3d - angle_t_scalar).abs() < 0.001,
            "3D ({:.4}°) should match scalar ({:.4}°)",
            angle_t_3d.to_degrees(),
            angle_t_scalar.to_degrees()
        );
    }

    #[test]
    fn test_trace_reflectance_at_each_surface() {
        let surfaces = [
            OpticalSurface {
                shape: SurfaceShape::Plane,
                z_position: 10.0,
                n_after: 1.5,
                aperture_radius: 25.0,
            },
            OpticalSurface {
                shape: SurfaceShape::Plane,
                z_position: 15.0,
                n_after: 1.0,
                aperture_radius: 25.0,
            },
        ];
        let ray = TraceRay {
            position: [0.0, 0.0, 0.0],
            direction: [0.0, 0.0, 1.0],
            n: 1.0,
        };
        let hits = trace_sequential(&ray, &surfaces).unwrap();
        for (i, hit) in hits.iter().enumerate() {
            assert!(
                (0.0..=1.0).contains(&hit.reflectance),
                "Surface {i}: reflectance {:.4} out of range",
                hit.reflectance
            );
        }
    }

    #[test]
    fn test_prism_zero_apex() {
        // Zero apex angle → zero deviation
        let dev = prism_deviation(0.0, 1.5).unwrap();
        assert!(dev.abs() < EPS, "Zero apex should give zero deviation");
    }

    // ── Prism tests ───────────────────────────────────────────────────────

    #[test]
    fn test_prism_deviation_60_degree() {
        // 60° equilateral prism with BK7 glass
        let apex = deg_to_rad(60.0);
        let n = SellmeierCoefficients::BK7.n_at(FRAUNHOFER_D);
        let dev = prism_deviation(apex, n).unwrap();
        // Minimum deviation for BK7 60° prism ≈ 40°
        let dev_deg = dev.to_degrees();
        assert!(
            dev_deg > 35.0 && dev_deg < 50.0,
            "60° BK7 prism deviation ≈ 40°, got {dev_deg}"
        );
    }

    #[test]
    fn test_prism_deviation_increases_with_n() {
        let apex = deg_to_rad(60.0);
        let dev_low = prism_deviation(apex, 1.4).unwrap();
        let dev_high = prism_deviation(apex, 1.7).unwrap();
        assert!(dev_high > dev_low, "Higher n should give larger deviation");
    }

    #[test]
    fn test_prism_deviation_tir() {
        // Very high n with large apex → TIR possible
        let apex = deg_to_rad(120.0);
        let result = prism_deviation(apex, 2.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_prism_dispersion_blue_more_than_red() {
        let apex = deg_to_rad(60.0);
        let s = SellmeierCoefficients::BK7;
        let dev_blue = prism_dispersion(apex, &s, FRAUNHOFER_F).unwrap();
        let dev_red = prism_dispersion(apex, &s, FRAUNHOFER_C).unwrap();
        assert!(
            dev_blue > dev_red,
            "Blue light should be deviated more than red"
        );
    }

    #[test]
    fn test_prism_angular_spread() {
        let apex = deg_to_rad(60.0);
        let s = SellmeierCoefficients::BK7;
        let (dev_short, dev_long, spread) =
            prism_angular_spread(apex, &s, FRAUNHOFER_F, FRAUNHOFER_C).unwrap();
        assert!(dev_short > dev_long);
        assert!(spread > 0.0);
        assert!(
            spread.to_degrees() < 5.0,
            "BK7 angular spread should be modest"
        );
    }

    #[test]
    fn test_prism_flint_more_dispersion_than_crown() {
        let apex = deg_to_rad(60.0);
        let (_, _, spread_crown) = prism_angular_spread(
            apex,
            &SellmeierCoefficients::BK7,
            FRAUNHOFER_F,
            FRAUNHOFER_C,
        )
        .unwrap();
        let (_, _, spread_flint) = prism_angular_spread(
            apex,
            &SellmeierCoefficients::SF11,
            FRAUNHOFER_F,
            FRAUNHOFER_C,
        )
        .unwrap();
        assert!(
            spread_flint > spread_crown,
            "Flint glass should have more dispersion than crown"
        );
    }
}
