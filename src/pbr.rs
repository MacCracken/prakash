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
    a2 / (PI * denom * denom).max(1e-15)
}

/// Beckmann normal distribution function.
///
/// Alternative to GGX, sharper specular highlights.
#[inline]
pub fn distribution_beckmann(n_dot_h: f64, roughness: f64) -> f64 {
    let a = roughness * roughness;
    let ndh = n_dot_h.clamp(0.0, 1.0);
    let cos2 = ndh * ndh;
    let tan2 = (1.0 - cos2) / cos2.max(1e-15);
    let exp_term = (-tan2 / (a * a)).exp();
    exp_term / (PI * a * a * cos2 * cos2).max(1e-15)
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
    ndv / (ndv * (1.0 - k) + k).max(1e-15)
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
pub fn cook_torrance(
    n_dot_h: f64,
    n_dot_v: f64,
    n_dot_l: f64,
    roughness: f64,
    f0: f64,
) -> f64 {
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

    #[test]
    fn test_fresnel_schlick_normal() {
        // At normal incidence, F = F0
        assert!((fresnel_schlick(0.04, 1.0) - 0.04).abs() < EPS);
    }

    #[test]
    fn test_fresnel_schlick_grazing() {
        // At grazing angle (cos=0), F → 1.0
        assert!((fresnel_schlick(0.04, 0.0) - 1.0).abs() < EPS);
    }

    #[test]
    fn test_fresnel_schlick_metal() {
        // Metal at normal: F0 = 0.95
        assert!((fresnel_schlick(0.95, 1.0) - 0.95).abs() < EPS);
    }

    #[test]
    fn test_fresnel_schlick_rgb() {
        let f = fresnel_schlick_rgb([0.04, 0.04, 0.04], 1.0);
        assert!((f[0] - 0.04).abs() < EPS);
        assert!((f[1] - 0.04).abs() < EPS);
        assert!((f[2] - 0.04).abs() < EPS);
    }

    #[test]
    fn test_ggx_smooth() {
        // Very smooth surface (roughness near 0): sharp peak at n·h = 1
        let d_peak = distribution_ggx(1.0, 0.01);
        let d_off = distribution_ggx(0.5, 0.01);
        assert!(d_peak > d_off * 100.0); // much sharper
    }

    #[test]
    fn test_ggx_rough() {
        // Rough surface: broader distribution
        let d_peak = distribution_ggx(1.0, 0.9);
        let d_off = distribution_ggx(0.5, 0.9);
        assert!(d_peak > d_off); // still peaks at 1
        assert!(d_peak < d_off * 100.0); // but much broader
    }

    #[test]
    fn test_beckmann_positive() {
        let d = distribution_beckmann(0.8, 0.5);
        assert!(d > 0.0);
    }

    #[test]
    fn test_geometry_smith_no_occlusion() {
        // Head-on view + light: minimal occlusion
        let g = geometry_smith(1.0, 1.0, 0.1);
        assert!(g > 0.9);
    }

    #[test]
    fn test_geometry_smith_grazing() {
        // Grazing angle: heavy occlusion
        let g = geometry_smith(0.01, 0.01, 0.5);
        assert!(g < 0.1);
    }

    #[test]
    fn test_cook_torrance_positive() {
        let spec = cook_torrance(0.9, 0.8, 0.7, 0.3, 0.04);
        assert!(spec > 0.0);
    }

    #[test]
    fn test_lambert_diffuse() {
        let d = lambert_diffuse(1.0);
        assert!((d - 1.0 / PI).abs() < EPS);
    }

    #[test]
    fn test_lambert_diffuse_rgb() {
        let d = lambert_diffuse_rgb([0.5, 0.3, 0.1]);
        assert!((d[0] - 0.5 / PI).abs() < EPS);
    }

    #[test]
    fn test_ior_to_f0_glass() {
        let f0 = ior_to_f0(1.5);
        // ((1 - 1.5) / (1 + 1.5))² = 0.04
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
        assert!(f0 > 0.15); // diamond is highly reflective
    }
}
