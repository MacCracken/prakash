//! Physically-based rendering primitives — BRDFs, Fresnel-Schlick, Cook-Torrance.
//!
//! These are the math functions behind PBR shading models. They produce
//! reflectance values that feed into rendering pipelines (aethersafta/kiran).

use std::f64::consts::PI;

// ── Fresnel-Schlick ─────────────────────────────────────────────────────────

/// Fresnel-Schlick approximation.
///
/// F(θ) = F0 + (1 - F0) · (1 - cos(θ))⁵
///
/// `f0` = reflectance at normal incidence (e.g. 0.04 for dielectrics, 0.95 for metals).
/// `cos_theta` = dot(view, half-vector), clamped to [0, 1].
#[inline]
pub fn fresnel_schlick(f0: f64, cos_theta: f64) -> f64 {
    let ct = cos_theta.clamp(0.0, 1.0);
    f0 + (1.0 - f0) * (1.0 - ct).powi(5)
}

/// Fresnel-Schlick for RGB reflectance (per-channel F0).
#[inline]
pub fn fresnel_schlick_rgb(f0: [f64; 3], cos_theta: f64) -> [f64; 3] {
    let ct = cos_theta.clamp(0.0, 1.0);
    let factor = (1.0 - ct).powi(5);
    [
        f0[0] + (1.0 - f0[0]) * factor,
        f0[1] + (1.0 - f0[1]) * factor,
        f0[2] + (1.0 - f0[2]) * factor,
    ]
}

// ── Normal Distribution Functions (NDF) ─────────────────────────────────────

/// GGX/Trowbridge-Reitz normal distribution function.
///
/// D(h) = α² / (π · ((n·h)²·(α²-1) + 1)²)
///
/// `n_dot_h` = dot(normal, half-vector).
/// `roughness` = material roughness (0 = mirror, 1 = fully rough).
#[inline]
pub fn distribution_ggx(n_dot_h: f64, roughness: f64) -> f64 {
    let a = roughness * roughness;
    let a2 = a * a;
    let ndh = n_dot_h.clamp(0.0, 1.0);
    let denom = ndh * ndh * (a2 - 1.0) + 1.0;
    a2 / (PI * denom * denom + 1e-15)
}

/// Beckmann normal distribution function.
///
/// Alternative to GGX, sharper specular highlights.
#[inline]
pub fn distribution_beckmann(n_dot_h: f64, roughness: f64) -> f64 {
    let a = roughness * roughness;
    let ndh = n_dot_h.clamp(0.0, 1.0);
    let cos2 = ndh * ndh;
    let tan2 = (1.0 - cos2) / (cos2 + 1e-15);
    let exp_term = (-tan2 / (a * a)).exp();
    exp_term / (PI * a * a * cos2 * cos2 + 1e-15)
}

// ── Geometry Functions ──────────────────────────────────────────────────────

/// Schlick-GGX geometry function (single direction).
///
/// G1(v) = n·v / (n·v · (1 - k) + k)
///
/// `k` depends on whether this is direct or IBL lighting.
#[inline]
pub fn geometry_schlick_ggx(n_dot_v: f64, roughness: f64) -> f64 {
    let r = roughness + 1.0;
    let k = (r * r) / 8.0; // k for direct lighting
    let ndv = n_dot_v.clamp(0.0, 1.0);
    ndv / (ndv * (1.0 - k) + k + 1e-15)
}

/// Smith's geometry function — combines view and light directions.
///
/// G(l, v) = G1(n·v) · G1(n·l)
#[inline]
pub fn geometry_smith(n_dot_v: f64, n_dot_l: f64, roughness: f64) -> f64 {
    geometry_schlick_ggx(n_dot_v, roughness) * geometry_schlick_ggx(n_dot_l, roughness)
}

// ── Cook-Torrance BRDF ──────────────────────────────────────────────────────

/// Cook-Torrance specular BRDF.
///
/// f_spec = D · F · G / (4 · (n·v) · (n·l))
///
/// Returns the specular reflectance multiplier.
#[inline]
pub fn cook_torrance(n_dot_h: f64, n_dot_v: f64, n_dot_l: f64, roughness: f64, f0: f64) -> f64 {
    let h_dot_v = n_dot_h; // approximation: h·v ≈ n·h for visualization
    let d = distribution_ggx(n_dot_h, roughness);
    let f = fresnel_schlick(f0, h_dot_v);
    let g = geometry_smith(n_dot_v, n_dot_l, roughness);
    let denom = 4.0 * n_dot_v.clamp(0.001, 1.0) * n_dot_l.clamp(0.001, 1.0);
    d * f * g / denom
}

// ── Lambert Diffuse ─────────────────────────────────────────────────────────

/// Lambertian diffuse BRDF: f_diff = albedo / π.
#[inline]
pub fn lambert_diffuse(albedo: f64) -> f64 {
    albedo / PI
}

/// Lambertian diffuse for RGB albedo.
#[inline]
pub fn lambert_diffuse_rgb(albedo: [f64; 3]) -> [f64; 3] {
    [albedo[0] / PI, albedo[1] / PI, albedo[2] / PI]
}

// ── IOR to F0 ───────────────────────────────────────────────────────────────

/// Convert index of refraction to F0 (reflectance at normal incidence).
///
/// F0 = ((n1 - n2) / (n1 + n2))²
/// For most materials, n1 = 1.0 (air).
#[inline]
pub fn ior_to_f0(ior: f64) -> f64 {
    let r = (1.0 - ior) / (1.0 + ior);
    r * r
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-6;

    // ── Fresnel-Schlick tests ─────────────────────────────────────────────

    #[test]
    fn test_fresnel_schlick_normal() {
        assert!((fresnel_schlick(0.04, 1.0) - 0.04).abs() < EPS);
    }

    #[test]
    fn test_fresnel_schlick_grazing() {
        assert!((fresnel_schlick(0.04, 0.0) - 1.0).abs() < EPS);
    }

    #[test]
    fn test_fresnel_schlick_metal() {
        assert!((fresnel_schlick(0.95, 1.0) - 0.95).abs() < EPS);
    }

    #[test]
    fn test_fresnel_schlick_monotonic() {
        // As cos_theta decreases (more grazing), reflectance should increase
        let mut prev = fresnel_schlick(0.04, 1.0);
        for i in (0..10).rev() {
            let ct = i as f64 / 10.0;
            let f = fresnel_schlick(0.04, ct);
            assert!(
                f >= prev - EPS,
                "Fresnel should increase as angle increases"
            );
            prev = f;
        }
    }

    #[test]
    fn test_fresnel_schlick_range() {
        for ct in [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0] {
            let f = fresnel_schlick(0.04, ct);
            assert!(
                (0.04 - EPS..=1.0 + EPS).contains(&f),
                "Fresnel out of range at cos_theta={ct}: {f}"
            );
        }
    }

    #[test]
    fn test_fresnel_schlick_clamps_input() {
        // Negative cos_theta should be clamped to 0
        let f = fresnel_schlick(0.04, -0.5);
        assert!((f - 1.0).abs() < EPS);
        // cos_theta > 1 should be clamped to 1
        let f = fresnel_schlick(0.04, 1.5);
        assert!((f - 0.04).abs() < EPS);
    }

    #[test]
    fn test_fresnel_schlick_rgb_normal() {
        let f = fresnel_schlick_rgb([0.04, 0.04, 0.04], 1.0);
        assert!((f[0] - 0.04).abs() < EPS);
        assert!((f[1] - 0.04).abs() < EPS);
        assert!((f[2] - 0.04).abs() < EPS);
    }

    #[test]
    fn test_fresnel_schlick_rgb_grazing() {
        let f = fresnel_schlick_rgb([0.04, 0.5, 0.95], 0.0);
        assert!((f[0] - 1.0).abs() < EPS);
        assert!((f[1] - 1.0).abs() < EPS);
        assert!((f[2] - 1.0).abs() < EPS);
    }

    #[test]
    fn test_fresnel_schlick_rgb_matches_scalar() {
        for ct in [0.0, 0.3, 0.7, 1.0] {
            let scalar = fresnel_schlick(0.04, ct);
            let rgb = fresnel_schlick_rgb([0.04, 0.04, 0.04], ct);
            assert!((rgb[0] - scalar).abs() < EPS);
        }
    }

    // ── NDF tests ─────────────────────────────────────────────────────────

    #[test]
    fn test_ggx_smooth() {
        let d_peak = distribution_ggx(1.0, 0.01);
        let d_off = distribution_ggx(0.5, 0.01);
        assert!(d_peak > d_off * 100.0);
    }

    #[test]
    fn test_ggx_rough() {
        let d_peak = distribution_ggx(1.0, 0.9);
        let d_off = distribution_ggx(0.5, 0.9);
        assert!(d_peak > d_off);
        assert!(d_peak < d_off * 100.0);
    }

    #[test]
    fn test_ggx_always_positive() {
        for ndh in [0.0, 0.1, 0.5, 0.9, 1.0] {
            for rough in [0.01, 0.1, 0.5, 0.9, 1.0] {
                let d = distribution_ggx(ndh, rough);
                assert!(d >= 0.0, "GGX negative at ndh={ndh}, rough={rough}: {d}");
            }
        }
    }

    #[test]
    fn test_ggx_peak_at_one() {
        // For any roughness, peak should be at n·h = 1
        for rough in [0.1, 0.3, 0.5, 0.7, 0.9] {
            let d_peak = distribution_ggx(1.0, rough);
            let d_half = distribution_ggx(0.5, rough);
            assert!(
                d_peak >= d_half,
                "GGX should peak at n·h=1 for roughness={rough}"
            );
        }
    }

    #[test]
    fn test_beckmann_positive() {
        let d = distribution_beckmann(0.8, 0.5);
        assert!(d > 0.0);
    }

    #[test]
    fn test_beckmann_peak_at_normal() {
        let d_normal = distribution_beckmann(1.0, 0.3);
        let d_off = distribution_beckmann(0.7, 0.3);
        assert!(d_normal > d_off);
    }

    #[test]
    fn test_beckmann_always_positive() {
        for ndh in [0.1, 0.3, 0.5, 0.7, 0.9, 1.0] {
            for rough in [0.1, 0.3, 0.5, 0.7, 0.9] {
                let d = distribution_beckmann(ndh, rough);
                assert!(d >= 0.0, "Beckmann negative at ndh={ndh}, rough={rough}");
            }
        }
    }

    // ── Geometry function tests ───────────────────────────────────────────

    #[test]
    fn test_geometry_schlick_ggx_normal() {
        let g = geometry_schlick_ggx(1.0, 0.5);
        assert!(g > 0.9, "Should have minimal occlusion head-on");
    }

    #[test]
    fn test_geometry_schlick_ggx_grazing() {
        let g = geometry_schlick_ggx(0.01, 0.5);
        assert!(g < 0.1, "Should have heavy occlusion at grazing");
    }

    #[test]
    fn test_geometry_schlick_ggx_range() {
        for ndv in [0.0, 0.1, 0.5, 0.9, 1.0] {
            for rough in [0.0, 0.3, 0.5, 0.7, 1.0] {
                let g = geometry_schlick_ggx(ndv, rough);
                assert!(
                    (0.0..=1.0 + EPS).contains(&g),
                    "G out of range at ndv={ndv}, rough={rough}: {g}"
                );
            }
        }
    }

    #[test]
    fn test_geometry_smith_no_occlusion() {
        let g = geometry_smith(1.0, 1.0, 0.1);
        assert!(g > 0.9);
    }

    #[test]
    fn test_geometry_smith_grazing() {
        let g = geometry_smith(0.01, 0.01, 0.5);
        assert!(g < 0.1);
    }

    #[test]
    fn test_geometry_smith_symmetric() {
        let g1 = geometry_smith(0.8, 0.6, 0.3);
        let g2 = geometry_smith(0.6, 0.8, 0.3);
        assert!((g1 - g2).abs() < EPS, "Smith geometry should be symmetric");
    }

    // ── Cook-Torrance tests ───────────────────────────────────────────────

    #[test]
    fn test_cook_torrance_positive() {
        let spec = cook_torrance(0.9, 0.8, 0.7, 0.3, 0.04);
        assert!(spec > 0.0);
    }

    #[test]
    fn test_cook_torrance_rougher_is_broader() {
        // At same geometry, rougher surface should have lower peak
        let smooth = cook_torrance(1.0, 0.8, 0.7, 0.1, 0.04);
        let rough = cook_torrance(1.0, 0.8, 0.7, 0.9, 0.04);
        assert!(
            smooth > rough,
            "Smooth surface should have higher specular peak"
        );
    }

    #[test]
    fn test_cook_torrance_metal_higher() {
        let dielectric = cook_torrance(0.9, 0.8, 0.7, 0.3, 0.04);
        let metal = cook_torrance(0.9, 0.8, 0.7, 0.3, 0.95);
        assert!(
            metal > dielectric,
            "Metal should have higher specular reflectance"
        );
    }

    // ── Lambert diffuse tests ─────────────────────────────────────────────

    #[test]
    fn test_lambert_diffuse() {
        let d = lambert_diffuse(1.0);
        assert!((d - 1.0 / PI).abs() < EPS);
    }

    #[test]
    fn test_lambert_diffuse_zero() {
        assert!(lambert_diffuse(0.0).abs() < EPS);
    }

    #[test]
    fn test_lambert_diffuse_rgb() {
        let d = lambert_diffuse_rgb([0.5, 0.3, 0.1]);
        assert!((d[0] - 0.5 / PI).abs() < EPS);
        assert!((d[1] - 0.3 / PI).abs() < EPS);
        assert!((d[2] - 0.1 / PI).abs() < EPS);
    }

    #[test]
    fn test_lambert_diffuse_rgb_matches_scalar() {
        let albedo = 0.7;
        let scalar = lambert_diffuse(albedo);
        let rgb = lambert_diffuse_rgb([albedo, albedo, albedo]);
        assert!((rgb[0] - scalar).abs() < EPS);
        assert!((rgb[1] - scalar).abs() < EPS);
        assert!((rgb[2] - scalar).abs() < EPS);
    }

    // ── IOR to F0 tests ──────────────────────────────────────────────────

    #[test]
    fn test_ior_to_f0_glass() {
        let f0 = ior_to_f0(1.5);
        assert!((f0 - 0.04).abs() < 0.001);
    }

    #[test]
    fn test_ior_to_f0_water() {
        let f0 = ior_to_f0(1.333);
        assert!((f0 - 0.02).abs() < 0.005);
    }

    #[test]
    fn test_ior_to_f0_diamond() {
        let f0 = ior_to_f0(2.417);
        assert!(f0 > 0.15);
    }

    #[test]
    fn test_ior_to_f0_air() {
        // n=1 → F0=0 (no reflection at identical media)
        assert!(ior_to_f0(1.0).abs() < EPS);
    }

    #[test]
    fn test_ior_to_f0_increases_with_ior() {
        let f0_water = ior_to_f0(1.333);
        let f0_glass = ior_to_f0(1.5);
        let f0_diamond = ior_to_f0(2.417);
        assert!(f0_water < f0_glass);
        assert!(f0_glass < f0_diamond);
    }

    #[test]
    fn test_ior_to_f0_range() {
        for n in [1.0, 1.1, 1.5, 2.0, 2.5, 3.0, 4.0] {
            let f0 = ior_to_f0(n);
            assert!((0.0..1.0).contains(&f0), "F0 out of range for n={n}: {f0}");
        }
    }
}
