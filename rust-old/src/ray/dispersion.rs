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
    /// Constant term of the Cauchy equation.
    pub b: f64,
    /// Dispersive coefficient (μm²).
    pub c: f64,
}

impl CauchyCoefficients {
    /// Refractive index at a given wavelength (micrometers).
    ///
    /// Wavelength must be positive. Returns `b` coefficient for very large wavelengths.
    #[must_use]
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
    /// First oscillator strength.
    pub b1: f64,
    /// First resonance wavelength squared (μm²).
    pub c1: f64,
    /// Second oscillator strength.
    pub b2: f64,
    /// Second resonance wavelength squared (μm²).
    pub c2: f64,
    /// Third oscillator strength.
    pub b3: f64,
    /// Third resonance wavelength squared (μm²).
    pub c3: f64,
}

impl SellmeierCoefficients {
    /// Refractive index at a given wavelength (micrometers).
    ///
    /// Returns n >= 1.0 for valid wavelengths away from resonance poles.
    /// Near resonances (λ² ≈ Cᵢ), results may be unphysical.
    #[must_use]
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

/// Fraunhofer d-line (helium yellow) wavelength in micrometers.
pub const FRAUNHOFER_D: f64 = 0.587_56;
/// Fraunhofer F-line (hydrogen blue) wavelength in micrometers.
pub const FRAUNHOFER_F: f64 = 0.486_13;
/// Fraunhofer C-line (hydrogen red) wavelength in micrometers.
pub const FRAUNHOFER_C: f64 = 0.656_27;

/// Abbe number (constringence) — a measure of a material's dispersion.
///
/// V = (n_d − 1) / (n_F − n_C)
///
/// Higher V means lower dispersion. Crown glasses typically V > 50,
/// flint glasses V < 50.
#[must_use]
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
#[must_use = "returns the minimum deviation angle"]
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
#[must_use = "returns the deviation angle"]
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
#[must_use = "returns the angular spread"]
#[inline]
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

// ── Additional Dispersion Models ─────────────────────────────────────────────

/// Herzberger dispersion coefficients (5-term).
///
/// n(λ) = A + B·L + C·L² + D·λ² + E·λ⁴
///
/// where L = 1/(λ² − 0.028) and λ is in micrometers.
/// Preferred for IR materials and wide spectral ranges.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct HerzbergerCoefficients {
    /// Constant term.
    pub a: f64,
    /// First L coefficient.
    pub b: f64,
    /// Second L coefficient.
    pub c: f64,
    /// λ² coefficient.
    pub d: f64,
    /// λ⁴ coefficient.
    pub e: f64,
}

impl HerzbergerCoefficients {
    /// Refractive index at a given wavelength (micrometers).
    #[must_use]
    #[inline]
    pub fn n_at(&self, wavelength_um: f64) -> f64 {
        let l2 = wavelength_um * wavelength_um;
        let l = 1.0 / (l2 - 0.028);
        self.a + self.b * l + self.c * l * l + self.d * l2 + self.e * l2 * l2
    }
}

/// Schott dispersion coefficients (6-term, legacy).
///
/// n²(λ) = a₀ + a₁·λ² + a₂·λ⁻² + a₃·λ⁻⁴ + a₄·λ⁻⁶ + a₅·λ⁻⁸
///
/// λ in micrometers. Used by legacy Schott glass catalogs (pre-1992).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SchottCoefficients {
    /// Constant term a₀.
    pub a0: f64,
    /// λ² coefficient a₁.
    pub a1: f64,
    /// λ⁻² coefficient a₂.
    pub a2: f64,
    /// λ⁻⁴ coefficient a₃.
    pub a3: f64,
    /// λ⁻⁶ coefficient a₄.
    pub a4: f64,
    /// λ⁻⁸ coefficient a₅.
    pub a5: f64,
}

impl SchottCoefficients {
    /// Refractive index at a given wavelength (micrometers).
    #[must_use]
    #[inline]
    pub fn n_at(&self, wavelength_um: f64) -> f64 {
        let l2 = wavelength_um * wavelength_um;
        let l2_inv = 1.0 / l2;
        let n2 = self.a0
            + self.a1 * l2
            + self.a2 * l2_inv
            + self.a3 * l2_inv * l2_inv
            + self.a4 * l2_inv * l2_inv * l2_inv
            + self.a5 * l2_inv * l2_inv * l2_inv * l2_inv;
        if n2 > 0.0 { n2.sqrt() } else { 1.0 }
    }

    /// N-BK7 in Schott format (legacy).
    pub const BK7: Self = Self {
        a0: 2.2718929,
        a1: -0.010108077,
        a2: 0.010592509,
        a3: 0.000200816,
        a4: -7.6472538e-6,
        a5: 4.9240991e-7,
    };
}

/// Conrady dispersion coefficients (3-term).
///
/// n(λ) = n₀ + A/λ + B/λ^3.5
///
/// λ in micrometers. Quick fit for visible-range glass data.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ConradyCoefficients {
    /// Constant (reference index).
    pub n0: f64,
    /// A coefficient (1/λ term).
    pub a: f64,
    /// B coefficient (1/λ^3.5 term).
    pub b: f64,
}

impl ConradyCoefficients {
    /// Refractive index at a given wavelength (micrometers).
    #[must_use]
    #[inline]
    pub fn n_at(&self, wavelength_um: f64) -> f64 {
        self.n0 + self.a / wavelength_um + self.b / wavelength_um.powf(3.5)
    }
}

// ── Partial Dispersion & Chromatic Aberration ────────────────────────────────

/// Partial dispersion ratio P_g,F.
///
/// P = (n_g − n_F) / (n_F − n_C)
///
/// Used for secondary spectrum analysis. For a doublet to be achromatic,
/// both elements must have equal P values. Normal glasses have P ≈ 0.53.
#[must_use]
#[inline]
pub fn partial_dispersion(sellmeier: &SellmeierCoefficients) -> f64 {
    let n_g = sellmeier.n_at(0.435_84); // Fraunhofer g-line
    let n_f = sellmeier.n_at(FRAUNHOFER_F);
    let n_c = sellmeier.n_at(FRAUNHOFER_C);
    let dn_fc = n_f - n_c;
    if dn_fc.abs() < 1e-15 {
        return 0.0;
    }
    (n_g - n_f) / dn_fc
}

/// Longitudinal chromatic aberration (axial color).
///
/// δf = f / V
///
/// where f is the focal length and V is the Abbe number.
/// Returns the axial focus shift between F and C Fraunhofer lines.
#[must_use]
#[inline]
pub fn longitudinal_chromatic_aberration(focal_length: f64, abbe_num: f64) -> f64 {
    focal_length / abbe_num
}

/// Lateral (transverse) chromatic aberration.
///
/// TCA = h · δf / f = h / V
///
/// where h is the image height. Returns the difference in image height
/// between F and C lines.
#[must_use]
#[inline]
pub fn lateral_chromatic_aberration(image_height: f64, abbe_num: f64) -> f64 {
    image_height / abbe_num
}

/// Secondary spectrum — residual chromatic aberration after achromatic correction.
///
/// δf_secondary = f · ΔP / V
///
/// where ΔP = P₁ - P₂ is the partial dispersion difference between the
/// two glass types of the doublet.
///
/// For an achromat (V₁ ≠ V₂ but same P), the secondary spectrum is zero.
/// For normal glasses, secondary spectrum ≈ f/2000.
#[must_use]
#[inline]
pub fn secondary_spectrum(focal_length: f64, delta_p: f64, abbe_num: f64) -> f64 {
    focal_length * delta_p / abbe_num
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

    // ── Additional dispersion model tests ────────────────────────────────

    #[test]
    fn test_schott_bk7_matches_sellmeier() {
        let n_schott = SchottCoefficients::BK7.n_at(FRAUNHOFER_D);
        let n_sellmeier = SellmeierCoefficients::BK7.n_at(FRAUNHOFER_D);
        assert!(
            (n_schott - n_sellmeier).abs() < 0.002,
            "Schott BK7 n_d={n_schott} vs Sellmeier n_d={n_sellmeier}"
        );
    }

    #[test]
    fn test_schott_reasonable_range() {
        for wl in [0.4, 0.5, 0.6, 0.7] {
            let n = SchottCoefficients::BK7.n_at(wl);
            assert!(n > 1.4 && n < 1.6, "Schott BK7 n={n} at {wl}μm");
        }
    }

    #[test]
    fn test_conrady_reasonable() {
        let c = ConradyCoefficients {
            n0: 1.51,
            a: 0.005,
            b: 0.0001,
        };
        let n = c.n_at(0.5876);
        assert!(n > 1.5 && n < 1.55, "Conrady n={n} at d-line");
    }

    #[test]
    fn test_herzberger_reasonable() {
        let h = HerzbergerCoefficients {
            a: 1.51,
            b: 0.001,
            c: 0.0,
            d: -0.001,
            e: 0.0,
        };
        let n = h.n_at(0.5876);
        assert!(n > 1.4 && n < 1.6, "Herzberger n={n}");
    }

    // ── Chromatic aberration tests ───────────────────────────────────────

    #[test]
    fn test_partial_dispersion_bk7() {
        let p = partial_dispersion(&SellmeierCoefficients::BK7);
        // BK7 partial dispersion P_gF ≈ 0.53
        assert!(
            (p - 0.53).abs() < 0.05,
            "BK7 partial dispersion ≈ 0.53, got {p}"
        );
    }

    #[test]
    fn test_longitudinal_ca() {
        let v = abbe_number(&SellmeierCoefficients::BK7);
        let lca = longitudinal_chromatic_aberration(100.0, v);
        // BK7 V ≈ 64, so LCA ≈ 100/64 ≈ 1.56mm
        assert!(lca > 1.0 && lca < 2.0, "LCA ≈ 1.5mm, got {lca}");
    }

    #[test]
    fn test_lateral_ca() {
        let v = abbe_number(&SellmeierCoefficients::BK7);
        let tca = lateral_chromatic_aberration(10.0, v);
        assert!(tca > 0.0 && tca < 1.0);
    }

    #[test]
    fn test_secondary_spectrum_zero_for_matched_p() {
        let ss = secondary_spectrum(100.0, 0.0, 64.0);
        assert!(
            ss.abs() < 1e-15,
            "Zero ΔP should give zero secondary spectrum"
        );
    }

    #[test]
    fn test_secondary_spectrum_positive() {
        let ss = secondary_spectrum(100.0, 0.01, 64.0);
        assert!(
            ss > 0.0,
            "Positive ΔP should give positive secondary spectrum"
        );
    }
}
