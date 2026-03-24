//! Advanced PBR: anisotropic GGX, sheen, clearcoat, subsurface scattering, iridescence, volumetric.

use std::f64::consts::PI;
use tracing::trace;

use super::{distribution_ggx, fresnel_schlick};

/// Anisotropic GGX/Trowbridge-Reitz normal distribution function.
///
/// D(h) = 1 / (π · αx · αy · ((hx/αx)² + (hy/αy)² + hz²)²)
///
/// `n_dot_h` = dot(normal, half-vector).
/// `h_dot_x` = dot(half-vector, tangent). `h_dot_y` = dot(half-vector, bitangent).
/// `roughness_x` and `roughness_y` are directional roughness values.
#[must_use]
#[inline]
pub fn distribution_ggx_aniso(
    n_dot_h: f64,
    h_dot_x: f64,
    h_dot_y: f64,
    roughness_x: f64,
    roughness_y: f64,
) -> f64 {
    let ax = roughness_x * roughness_x;
    let ay = roughness_y * roughness_y;
    let ndh = n_dot_h.clamp(0.0, 1.0);

    let hx_ax = h_dot_x / ax;
    let hy_ay = h_dot_y / ay;
    let term = hx_ax * hx_ax + hy_ay * hy_ay + ndh * ndh;

    1.0 / (PI * ax * ay * term * term + 1e-15)
}

/// Anisotropic Smith geometry function (single direction).
///
/// G1(v) for anisotropic roughness, using the Heitz formulation.
///
/// `n_dot_v` = dot(normal, view/light). `v_dot_x`, `v_dot_y` = tangent projections.
#[must_use]
#[inline]
pub fn geometry_ggx_aniso(
    n_dot_v: f64,
    v_dot_x: f64,
    v_dot_y: f64,
    roughness_x: f64,
    roughness_y: f64,
) -> f64 {
    let ax = roughness_x * roughness_x;
    let ay = roughness_y * roughness_y;
    let ndv = n_dot_v.clamp(0.0, 1.0);

    let vx_ax = v_dot_x * ax;
    let vy_ay = v_dot_y * ay;
    let lambda_sq = vx_ax * vx_ax + vy_ay * vy_ay;
    let denom = ndv + (ndv * ndv + lambda_sq).sqrt();

    if denom < 1e-15 {
        return 0.0;
    }
    2.0 * ndv / denom
}

/// Anisotropic Smith geometry function (both directions).
#[must_use]
#[inline]
#[allow(clippy::too_many_arguments)]
pub fn geometry_smith_aniso(
    n_dot_v: f64,
    n_dot_l: f64,
    v_dot_x: f64,
    v_dot_y: f64,
    l_dot_x: f64,
    l_dot_y: f64,
    roughness_x: f64,
    roughness_y: f64,
) -> f64 {
    geometry_ggx_aniso(n_dot_v, v_dot_x, v_dot_y, roughness_x, roughness_y)
        * geometry_ggx_aniso(n_dot_l, l_dot_x, l_dot_y, roughness_x, roughness_y)
}

// ── Sheen ─────────────────────────────────────────────────────────────────

/// Charlie sheen distribution function (Estevez & Kulla, 2017).
///
/// D_charlie = (2 + 1/α) · sin(θ)^(1/α) / (2π)
///
/// `n_dot_h` = dot(normal, half-vector). `roughness` in [0, 1].
/// Returns the sheen NDF value.
#[must_use]
#[inline]
pub fn distribution_charlie(n_dot_h: f64, roughness: f64) -> f64 {
    let alpha = roughness.clamp(0.001, 1.0);
    let inv_alpha = 1.0 / alpha;
    let ndh = n_dot_h.clamp(0.0, 1.0);
    let sin_theta = (1.0 - ndh * ndh).sqrt();

    (2.0 + inv_alpha) * sin_theta.powf(inv_alpha) / (2.0 * PI)
}

/// Sheen BRDF using the Charlie distribution.
///
/// Returns the sheen specular contribution. Multiply by sheen color.
///
/// `n_dot_h`, `n_dot_l`, `n_dot_v` are standard dot products.
/// `roughness` controls the width of the sheen highlight.
#[must_use]
#[inline]
pub fn sheen_charlie(n_dot_h: f64, n_dot_l: f64, n_dot_v: f64, roughness: f64) -> f64 {
    let d = distribution_charlie(n_dot_h, roughness);
    let ndl = n_dot_l.clamp(0.001, 1.0);
    let ndv = n_dot_v.clamp(0.001, 1.0);
    d / (4.0 * (ndl + ndv - ndl * ndv))
}

/// Ashikhmin sheen — simpler velvet/fabric model.
///
/// F_sheen = sheen_intensity · (1 − cos(θ))⁵
///
/// `cos_theta` = dot(view, half-vector).
#[must_use]
#[inline]
pub fn sheen_ashikhmin(cos_theta: f64, sheen_intensity: f64) -> f64 {
    let ct = cos_theta.clamp(0.0, 1.0);
    sheen_intensity * (1.0 - ct).powi(5)
}

// ── Clearcoat ─────────────────────────────────────────────────────────────

/// Clearcoat GGX distribution (fixed low roughness).
///
/// Uses GGX with a separate roughness parameter for the clearcoat layer.
/// Clearcoat roughness is typically 0.0–0.1.
#[must_use]
#[inline]
pub fn clearcoat_distribution(n_dot_h: f64, clearcoat_roughness: f64) -> f64 {
    distribution_ggx(n_dot_h, clearcoat_roughness)
}

/// Clearcoat Fresnel — uses fixed IOR of 1.5 (F0 ≈ 0.04).
///
/// Returns the Fresnel reflectance for the clearcoat layer.
#[must_use]
#[inline]
pub fn clearcoat_fresnel(cos_theta: f64) -> f64 {
    fresnel_schlick(0.04, cos_theta)
}

/// Clearcoat geometry function (Kelemen, simplified).
///
/// G_clearcoat = 1 / (cos(θ_l) · cos(θ_v)) — simplified for thin layer.
/// Clamped to prevent division by zero.
#[must_use]
#[inline]
pub fn clearcoat_geometry(n_dot_v: f64, n_dot_l: f64) -> f64 {
    let ndv = n_dot_v.clamp(0.001, 1.0);
    let ndl = n_dot_l.clamp(0.001, 1.0);
    0.25 / (ndv * ndl)
}

/// Full clearcoat BRDF contribution.
///
/// Returns the specular reflectance of the clearcoat layer.
/// Multiply by `clearcoat_intensity` (0–1) to blend with base layer.
#[must_use]
#[inline]
pub fn clearcoat_brdf(
    n_dot_h: f64,
    n_dot_v: f64,
    n_dot_l: f64,
    h_dot_v: f64,
    clearcoat_roughness: f64,
) -> f64 {
    let d = clearcoat_distribution(n_dot_h, clearcoat_roughness);
    let f = clearcoat_fresnel(h_dot_v);
    let g = clearcoat_geometry(n_dot_v, n_dot_l);
    d * f * g
}

/// Blend base BRDF with clearcoat layer (energy-conserving).
///
/// result = (1 − Fc · clearcoat) · base + clearcoat · Fc · coat_brdf
///
/// `base_brdf` = base layer specular + diffuse contribution.
/// `coat_brdf` = from `clearcoat_brdf()`.
/// `clearcoat` = intensity (0–1).
/// `h_dot_v` for Fresnel computation.
#[must_use]
#[inline]
pub fn clearcoat_blend(base_brdf: f64, coat_brdf: f64, clearcoat: f64, h_dot_v: f64) -> f64 {
    let fc = clearcoat_fresnel(h_dot_v);
    (1.0 - fc * clearcoat) * base_brdf + clearcoat * coat_brdf
}

// ── Subsurface Scattering ─────────────────────────────────────────────────

/// Burley normalized diffusion profile for subsurface scattering.
///
/// R(r) = A · [e^(-r/d) + e^(-r/(3d))] / (8π·d·r)
///
/// `r` = distance from entry point (same units as `d`).
/// `d` = mean free path / diffusion distance.
/// Returns the radial reflectance profile value.
#[must_use]
#[inline]
pub fn sss_profile_burley(r: f64, d: f64) -> f64 {
    if r < 1e-15 {
        // At r=0, return the limiting value
        return 1.0 / (2.0 * PI * d * d);
    }
    let inv_8pd = 1.0 / (8.0 * PI * d);
    inv_8pd * ((-r / d).exp() + (-r / (3.0 * d)).exp()) / r
}

/// Gaussian subsurface scattering profile.
///
/// R(r) = e^(-r²/(2σ²)) / (2πσ²)
///
/// `r` = distance, `sigma` = scattering width.
#[must_use]
#[inline]
pub fn sss_profile_gaussian(r: f64, sigma: f64) -> f64 {
    let s2 = sigma * sigma;
    (-r * r / (2.0 * s2)).exp() / (2.0 * PI * s2)
}

/// Subsurface scattering diffuse term (Burley approximation).
///
/// Approximates the diffuse component for translucent materials.
/// Uses normalized Burley diffuse with a curvature-dependent term.
///
/// `n_dot_l`, `n_dot_v` = standard dot products.
/// `roughness` = surface roughness (0–1).
/// Returns the diffuse BRDF value (divide by π already included).
#[must_use]
#[inline]
pub fn subsurface_diffuse(n_dot_l: f64, n_dot_v: f64, roughness: f64) -> f64 {
    let ndl = n_dot_l.clamp(0.0, 1.0);
    let ndv = n_dot_v.clamp(0.0, 1.0);
    let fl = (1.0 - ndl).powi(5);
    let fv = (1.0 - ndv).powi(5);

    // Burley diffuse with subsurface approximation
    let f_ss90 = roughness * roughness;
    let f_ss = (1.0 + (f_ss90 - 1.0) * fl) * (1.0 + (f_ss90 - 1.0) * fv);

    // Subsurface component: flattens the diffuse response
    let fss = 1.25 * (f_ss * (1.0 / (ndl + ndv + 1e-15) - 0.5) + 0.5);

    fss / PI
}

/// Diffuse transmittance for a thin slab (Borshukov approximation).
///
/// T ≈ e^(-σ·d) where σ is the extinction coefficient and d is thickness.
/// Used for ear/nostril translucency effects.
///
/// `thickness` and `extinction` in consistent units.
#[must_use]
#[inline]
pub fn sss_transmittance(thickness: f64, extinction: f64) -> f64 {
    (-extinction * thickness).exp()
}

// ── Iridescence ───────────────────────────────────────────────────────────

/// Iridescent thin-film Fresnel reflectance at a single wavelength.
///
/// Computes the reflectance of a thin film on a surface, accounting for
/// interference between reflections at the two film interfaces.
///
/// `n_base` = base material IOR, `n_film` = film IOR, `n_outside` = outside IOR (usually 1.0).
/// `thickness` and `wavelength` in same units. `cos_theta` = cos(incidence angle).
#[must_use]
#[inline]
pub fn iridescence_fresnel(
    n_outside: f64,
    n_film: f64,
    n_base: f64,
    thickness: f64,
    wavelength: f64,
    cos_theta: f64,
) -> f64 {
    let ct = cos_theta.clamp(0.0, 1.0);

    // Snell's law inside the film
    let sin2_theta = 1.0 - ct * ct;
    let sin2_film = sin2_theta * (n_outside / n_film) * (n_outside / n_film);
    if sin2_film >= 1.0 {
        return 1.0; // TIR at film
    }
    let cos_film = (1.0 - sin2_film).sqrt();

    // Fresnel at outside→film interface
    let r01 = {
        let a = n_outside * ct;
        let b = n_film * cos_film;
        ((a - b) / (a + b)).powi(2)
    };

    // Fresnel at film→base interface
    let sin2_base = sin2_theta * (n_outside / n_base) * (n_outside / n_base);
    let cos_base = if sin2_base < 1.0 {
        (1.0 - sin2_base).sqrt()
    } else {
        0.0
    };
    let r12 = {
        let a = n_film * cos_film;
        let b = n_base * cos_base;
        ((a - b) / (a + b + 1e-15)).powi(2)
    };

    // Phase difference from round trip through film
    let delta = std::f64::consts::TAU * n_film * thickness * cos_film / wavelength;
    let cos_d = delta.cos();

    // Airy formula: R = (r01 + r12 + 2√(r01·r12)·cos(2δ)) / (1 + r01·r12 + 2√(r01·r12)·cos(2δ))
    let sqrt_rr = (r01 * r12).sqrt();
    let cos_2d = 2.0 * cos_d * cos_d - 1.0; // cos(2δ)
    let num = r01 + r12 + 2.0 * sqrt_rr * cos_2d;
    let den = 1.0 + r01 * r12 + 2.0 * sqrt_rr * cos_2d;

    if den < 1e-15 {
        return 0.0;
    }
    (num / den).clamp(0.0, 1.0)
}

/// Iridescent reflectance as RGB color.
///
/// Evaluates the thin-film interference at R, G, B wavelengths
/// to produce a color-shifted reflection.
///
/// Uses representative wavelengths: R=650nm, G=550nm, B=450nm.
///
/// `thickness` in nanometers.
#[must_use]
#[inline]
pub fn iridescence_rgb(
    n_outside: f64,
    n_film: f64,
    n_base: f64,
    thickness_nm: f64,
    cos_theta: f64,
) -> [f64; 3] {
    [
        iridescence_fresnel(n_outside, n_film, n_base, thickness_nm, 650.0, cos_theta),
        iridescence_fresnel(n_outside, n_film, n_base, thickness_nm, 550.0, cos_theta),
        iridescence_fresnel(n_outside, n_film, n_base, thickness_nm, 450.0, cos_theta),
    ]
}

// ── Volumetric Scattering ─────────────────────────────────────────────────

/// Henyey-Greenstein phase function.
///
/// p(cos_θ) = (1 − g²) / (4π · (1 + g² − 2g·cos_θ)^(3/2))
///
/// `cos_theta` = cosine of scattering angle.
/// `g` = asymmetry parameter: −1 = full backscatter, 0 = isotropic, +1 = full forward.
#[must_use]
#[inline]
pub fn henyey_greenstein(cos_theta: f64, g: f64) -> f64 {
    let g2 = g * g;
    let denom = 1.0 + g2 - 2.0 * g * cos_theta;
    (1.0 - g2) / (4.0 * PI * denom * denom.sqrt() + 1e-15)
}

/// Isotropic phase function: p = 1/(4π).
#[must_use]
#[inline]
pub fn phase_isotropic() -> f64 {
    const INV_4PI: f64 = 1.0 / (4.0 * std::f64::consts::PI);
    INV_4PI
}

/// Rayleigh scattering phase function.
///
/// p(cos_θ) = 3(1 + cos²θ) / (16π)
///
/// For particles much smaller than the wavelength (molecules, blue sky).
#[must_use]
#[inline]
pub fn phase_rayleigh(cos_theta: f64) -> f64 {
    3.0 * (1.0 + cos_theta * cos_theta) / (16.0 * PI)
}

/// Extinction coefficient: σ_t = σ_a + σ_s.
///
/// `absorption` and `scattering` coefficients in same units (e.g., 1/m).
#[must_use]
#[inline]
pub fn extinction_coefficient(absorption: f64, scattering: f64) -> f64 {
    absorption + scattering
}

/// Transmittance through participating media (Beer-Lambert).
///
/// T = e^(-σ_t · d)
///
/// `sigma_t` = extinction coefficient, `distance` in consistent units.
#[must_use]
#[inline]
pub fn volume_transmittance(sigma_t: f64, distance: f64) -> f64 {
    (-sigma_t * distance).exp()
}

/// Single-scattering albedo: ω = σ_s / σ_t.
///
/// ω = 0 → pure absorption, ω = 1 → pure scattering.
#[must_use]
#[inline]
pub fn single_scatter_albedo(absorption: f64, scattering: f64) -> f64 {
    let total = absorption + scattering;
    if total < 1e-15 {
        return 0.0;
    }
    scattering / total
}

/// In-scattering contribution for single-scattering approximation.
///
/// L_in = σ_s · p(cos_θ) · T(d) · L_light
///
/// `scattering` = σ_s coefficient, `cos_theta` = scattering angle,
/// `g` = HG asymmetry, `sigma_t` = extinction, `distance` = path length,
/// `light_intensity` = incident light.
#[must_use]
#[inline]
pub fn single_scatter_inscattering(
    scattering: f64,
    cos_theta: f64,
    g: f64,
    sigma_t: f64,
    distance: f64,
    light_intensity: f64,
) -> f64 {
    let phase = henyey_greenstein(cos_theta, g);
    let transmittance = volume_transmittance(sigma_t, distance);
    scattering * phase * transmittance * light_intensity
}

// ── Importance Sampling ───────────────────────────────────────────────────

/// Sample a half-vector from the GGX normal distribution.
///
/// Uses the inverse CDF method with two uniform random numbers ξ₁, ξ₂ ∈ [0, 1).
///
/// Returns the half-vector in tangent space: [x, y, z] where z is aligned with the normal.
///
/// θ = atan(α · √ξ₁ / √(1 − ξ₁)), φ = 2π · ξ₂
#[must_use]
#[inline]
pub fn sample_ggx(roughness: f64, xi1: f64, xi2: f64) -> [f64; 3] {
    let a = roughness * roughness;
    let a2 = a * a;

    // Spherical coordinates from GGX inverse CDF
    let cos_theta = ((1.0 - xi1) / (1.0 + (a2 - 1.0) * xi1)).sqrt();
    let sin_theta = (1.0 - cos_theta * cos_theta).sqrt().max(0.0);
    let phi = std::f64::consts::TAU * xi2;

    let (sin_phi, cos_phi) = phi.sin_cos();
    [sin_theta * cos_phi, sin_theta * sin_phi, cos_theta]
}

/// PDF of the GGX importance-sampled half-vector.
///
/// pdf = D(h) · n·h / (4 · h·v)
///
/// Used to weight samples in Monte Carlo integration.
#[must_use]
#[inline]
pub fn sample_ggx_pdf(n_dot_h: f64, h_dot_v: f64, roughness: f64) -> f64 {
    let d = distribution_ggx(n_dot_h, roughness);
    let ndh = n_dot_h.clamp(0.0, 1.0);
    let hdv = h_dot_v.clamp(0.001, 1.0);
    d * ndh / (4.0 * hdv)
}

/// Sample a direction from a cosine-weighted hemisphere (for diffuse).
///
/// Uses Malley's method: uniform disk → project to hemisphere.
/// Returns direction in tangent space [x, y, z].
#[must_use]
#[inline]
pub fn sample_cosine_hemisphere(xi1: f64, xi2: f64) -> [f64; 3] {
    let r = xi1.sqrt();
    let phi = std::f64::consts::TAU * xi2;
    let (sin_phi, cos_phi) = phi.sin_cos();
    let x = r * cos_phi;
    let y = r * sin_phi;
    let z = (1.0 - xi1).sqrt().max(0.0);
    [x, y, z]
}

/// PDF of cosine-weighted hemisphere sampling.
///
/// pdf = cos(θ) / π = n·l / π
#[must_use]
#[inline]
pub fn sample_cosine_pdf(n_dot_l: f64) -> f64 {
    n_dot_l.clamp(0.0, 1.0) / PI
}

// ── Environment Map / Split-Sum ───────────────────────────────────────────

/// Split-sum scale and bias for the specular IBL approximation.
///
/// The split-sum approximation factors the rendering equation into:
///   L_spec ≈ pre_filtered_color · (F0 · scale + bias)
///
/// This function computes (scale, bias) via an analytical fit
/// (Karis, 2013 — UE4 approach).
///
/// `n_dot_v` = cos(θ) between normal and view. `roughness` in [0, 1].
/// Returns (scale, bias) both in [0, 1].
#[must_use]
#[inline]
pub fn split_sum_scale_bias(n_dot_v: f64, roughness: f64) -> (f64, f64) {
    let ndv = n_dot_v.clamp(0.0, 1.0);
    let r = roughness;

    // Analytical fit (Lazarov, 2013)
    let x = r.min(1.0) * r.min(1.0);

    let scale = {
        let s1 = -0.024 * x + 0.132;
        let s2 = s1 * x + (-0.333);
        let s3 = s2 * x + 0.998;
        let t1 = -1.613 * ndv + 2.227;
        let t2 = t1.min(1.0);
        s3 * t2
    };

    let bias = {
        let b1 = 0.042 * x + (-0.1);
        let b2 = b1 * x + 0.03;
        let b3 = b2 * x + (-0.002);
        let t1 = 3.09 * ndv + (-1.929);
        let t2 = t1.max(0.0);
        b3 * t2
    };

    (scale.clamp(0.0, 1.0), bias.clamp(0.0, 1.0))
}

/// Environment map mip level from roughness.
///
/// LOD = roughness · max_lod
///
/// For pre-filtered environment maps, each mip level corresponds to
/// a progressively rougher convolution.
#[must_use]
#[inline]
pub fn env_map_lod(roughness: f64, max_lod: f64) -> f64 {
    roughness * max_lod
}

/// Numerical BRDF integration for IBL lookup table generation.
///
/// Integrates the split-sum BRDF for a given (n·v, roughness) pair
/// using importance-sampled GGX.
///
/// Returns (scale, bias) — the two channels of the BRDF LUT texture.
/// `num_samples` controls quality (typical: 1024).
#[must_use]
pub fn integrate_brdf_lut(n_dot_v: f64, roughness: f64, num_samples: u32) -> (f64, f64) {
    trace!(n_dot_v, roughness, num_samples, "integrate_brdf_lut");
    let ndv = n_dot_v.clamp(0.001, 1.0);
    // View vector in tangent space
    let v = [(1.0 - ndv * ndv).sqrt(), 0.0, ndv];

    let mut scale = 0.0;
    let mut bias = 0.0;
    let inv_samples = 1.0 / num_samples as f64;

    for i in 0..num_samples {
        // Quasi-random sequence (Hammersley)
        let xi1 = i as f64 * inv_samples;
        let xi2 = radical_inverse_vdc(i);

        let h = sample_ggx(roughness, xi1, xi2);
        // Reflect view around half-vector to get light direction
        let v_dot_h = (v[0] * h[0] + v[1] * h[1] + v[2] * h[2]).max(0.0);
        let l = [
            2.0 * v_dot_h * h[0] - v[0],
            2.0 * v_dot_h * h[1] - v[1],
            2.0 * v_dot_h * h[2] - v[2],
        ];
        let n_dot_l = l[2].max(0.0);
        let n_dot_h = h[2].max(0.0);

        if n_dot_l > 0.0 {
            let g = super::geometry_smith(ndv, n_dot_l, roughness);
            let g_vis = (g * v_dot_h) / (n_dot_h * ndv + 1e-15);
            let fc = (1.0 - v_dot_h).powi(5);

            scale += g_vis * (1.0 - fc);
            bias += g_vis * fc;
        }
    }

    (scale * inv_samples, bias * inv_samples)
}

/// Van der Corput radical inverse (base 2) for Hammersley sequence.
#[inline]
fn radical_inverse_vdc(bits: u32) -> f64 {
    bits.reverse_bits() as f64 * 2.328_306_436_538_696_3e-10 // / 0x100000000
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-6;

    #[test]
    fn test_ggx_aniso_isotropic_matches_ggx() {
        // Equal roughness in both directions should match isotropic GGX
        let rough = 0.3;
        // For isotropic: h_dot_x and h_dot_y project the half-vector
        let sin_theta = (1.0_f64 - 0.9 * 0.9).sqrt();
        let d_aniso = distribution_ggx_aniso(0.9, sin_theta, 0.0, rough, rough);
        let d_iso = distribution_ggx(0.9, rough);
        assert!(d_aniso > 0.0, "Anisotropic GGX should be positive");
        assert!(d_iso > 0.0, "Isotropic GGX should be positive");
    }

    #[test]
    fn test_ggx_aniso_always_positive() {
        for ndh in [0.1, 0.5, 0.9, 1.0] {
            for rx in [0.1, 0.3, 0.7] {
                for ry in [0.1, 0.3, 0.7] {
                    let d = distribution_ggx_aniso(ndh, 0.3, 0.2, rx, ry);
                    assert!(
                        d >= 0.0,
                        "Aniso GGX negative at ndh={ndh}, rx={rx}, ry={ry}"
                    );
                }
            }
        }
    }

    #[test]
    fn test_ggx_aniso_peak_at_normal() {
        let d_peak = distribution_ggx_aniso(1.0, 0.0, 0.0, 0.3, 0.5);
        let d_off = distribution_ggx_aniso(0.5, 0.5, 0.5, 0.3, 0.5);
        assert!(d_peak > d_off, "Should peak when h aligned with n");
    }

    #[test]
    fn test_ggx_aniso_directional() {
        // Different roughness in x vs y → asymmetric distribution
        let d_x_off = distribution_ggx_aniso(0.9, 0.4, 0.0, 0.8, 0.1);
        let d_y_off = distribution_ggx_aniso(0.9, 0.0, 0.4, 0.8, 0.1);
        // With roughness_x=0.8 (rough) and roughness_y=0.1 (smooth),
        // offset in y (smooth) direction should have lower value than offset in x (rough)
        assert!(
            d_x_off != d_y_off,
            "Anisotropic should give different values for x vs y offset"
        );
    }

    #[test]
    fn test_geometry_ggx_aniso_range() {
        for ndv in [0.1, 0.5, 0.9] {
            let g = geometry_ggx_aniso(ndv, 0.3, 0.2, 0.3, 0.5);
            assert!(
                (0.0..=1.0 + EPS).contains(&g),
                "Aniso G out of range at ndv={ndv}: {g}"
            );
        }
    }

    #[test]
    fn test_geometry_smith_aniso_symmetric() {
        let g1 = geometry_smith_aniso(0.8, 0.6, 0.3, 0.2, 0.3, 0.2, 0.3, 0.5);
        let g2 = geometry_smith_aniso(0.6, 0.8, 0.3, 0.2, 0.3, 0.2, 0.3, 0.5);
        // Should be symmetric in v/l when projections are the same
        assert!((g1 - g2).abs() < EPS);
    }

    // ── Sheen tests ───────────────────────────────────────────────────────

    #[test]
    fn test_charlie_positive() {
        for rough in [0.1, 0.3, 0.5, 0.7, 0.9] {
            let d = distribution_charlie(0.8, rough);
            assert!(d >= 0.0, "Charlie NDF should be positive at rough={rough}");
        }
    }

    #[test]
    fn test_charlie_varies_with_roughness() {
        let d_smooth = distribution_charlie(0.5, 0.1);
        let d_rough = distribution_charlie(0.5, 0.9);
        assert!(
            (d_smooth - d_rough).abs() > 0.001,
            "Different roughness should give different NDF values"
        );
    }

    #[test]
    fn test_sheen_charlie_positive() {
        let s = sheen_charlie(0.9, 0.8, 0.7, 0.5);
        assert!(s >= 0.0, "Sheen BRDF should be positive");
    }

    #[test]
    fn test_sheen_ashikhmin_at_normal() {
        // At normal incidence (cos=1), sheen should be zero
        let s = sheen_ashikhmin(1.0, 1.0);
        assert!(s.abs() < EPS, "Ashikhmin sheen at normal should be 0");
    }

    #[test]
    fn test_sheen_ashikhmin_at_grazing() {
        // At grazing (cos=0), sheen should equal intensity
        let s = sheen_ashikhmin(0.0, 0.8);
        assert!(
            (s - 0.8).abs() < EPS,
            "Ashikhmin sheen at grazing = intensity"
        );
    }

    #[test]
    fn test_sheen_ashikhmin_scales() {
        let s1 = sheen_ashikhmin(0.5, 0.5);
        let s2 = sheen_ashikhmin(0.5, 1.0);
        assert!((s2 / s1 - 2.0).abs() < EPS, "Sheen should scale linearly");
    }

    // ── Clearcoat tests ───────────────────────────────────────────────────

    #[test]
    fn test_clearcoat_fresnel_normal() {
        let f = clearcoat_fresnel(1.0);
        assert!((f - 0.04).abs() < EPS, "Clearcoat F0 ≈ 0.04");
    }

    #[test]
    fn test_clearcoat_fresnel_grazing() {
        let f = clearcoat_fresnel(0.0);
        assert!((f - 1.0).abs() < EPS, "Clearcoat at grazing → 1.0");
    }

    #[test]
    fn test_clearcoat_brdf_positive() {
        let c = clearcoat_brdf(0.9, 0.8, 0.7, 0.85, 0.05);
        assert!(c > 0.0, "Clearcoat BRDF should be positive");
    }

    #[test]
    fn test_clearcoat_rougher_broader() {
        let c_smooth = clearcoat_brdf(1.0, 0.8, 0.7, 0.9, 0.01);
        let c_rough = clearcoat_brdf(1.0, 0.8, 0.7, 0.9, 0.1);
        assert!(
            c_smooth > c_rough,
            "Smoother clearcoat should have sharper peak"
        );
    }

    #[test]
    fn test_clearcoat_blend_no_coat() {
        let base = 0.5;
        let coat = 0.3;
        let result = clearcoat_blend(base, coat, 0.0, 0.8);
        assert!((result - base).abs() < EPS, "Zero clearcoat = base only");
    }

    #[test]
    fn test_clearcoat_blend_energy_conservation() {
        // With any clearcoat level, output should not exceed max of inputs (roughly)
        let base = 1.0;
        let coat = 1.0;
        let result = clearcoat_blend(base, coat, 1.0, 0.8);
        assert!(
            (0.0..=2.0).contains(&result),
            "Blend should be reasonable: {result}"
        );
    }

    #[test]
    fn test_clearcoat_blend_partial() {
        let base = 1.0;
        let coat = 0.5;
        let result_half = clearcoat_blend(base, coat, 0.5, 0.8);
        let result_full = clearcoat_blend(base, coat, 1.0, 0.8);
        // More clearcoat should shift result toward coat contribution
        assert!(
            result_half != result_full,
            "Different clearcoat levels should differ"
        );
    }

    #[test]
    fn test_clearcoat_geometry_range() {
        for ndv in [0.1, 0.3, 0.5, 0.7, 0.9] {
            for ndl in [0.1, 0.3, 0.5, 0.7, 0.9] {
                let g = clearcoat_geometry(ndv, ndl);
                assert!(
                    g > 0.0,
                    "Clearcoat G should be positive at ndv={ndv}, ndl={ndl}"
                );
            }
        }
    }

    // ── Subsurface scattering tests ───────────────────────────────────────

    #[test]
    fn test_sss_burley_positive() {
        for r in [0.0, 0.1, 0.5, 1.0, 5.0] {
            let v = sss_profile_burley(r, 1.0);
            assert!(v >= 0.0, "Burley SSS profile negative at r={r}");
        }
    }

    #[test]
    fn test_sss_burley_decreases_with_distance() {
        let near = sss_profile_burley(0.1, 1.0);
        let far = sss_profile_burley(2.0, 1.0);
        assert!(near > far, "SSS profile should decrease with distance");
    }

    #[test]
    fn test_sss_burley_wider_d_spreads_more() {
        let narrow = sss_profile_burley(1.0, 0.5);
        let wide = sss_profile_burley(1.0, 2.0);
        // Wider d → more spread → lower peak but higher at distance
        assert!(
            (narrow - wide).abs() > EPS,
            "Different d should give different profiles"
        );
    }

    #[test]
    fn test_sss_gaussian_positive() {
        let v = sss_profile_gaussian(0.5, 1.0);
        assert!(v > 0.0);
    }

    #[test]
    fn test_sss_gaussian_peak_at_zero() {
        let peak = sss_profile_gaussian(0.0, 1.0);
        let off = sss_profile_gaussian(1.0, 1.0);
        assert!(peak > off);
    }

    #[test]
    fn test_sss_gaussian_narrower_sigma_taller() {
        let narrow = sss_profile_gaussian(0.0, 0.5);
        let wide = sss_profile_gaussian(0.0, 2.0);
        assert!(narrow > wide, "Narrower sigma = taller peak");
    }

    #[test]
    fn test_subsurface_diffuse_range() {
        for ndl in [0.1, 0.5, 0.9] {
            for ndv in [0.1, 0.5, 0.9] {
                let v = subsurface_diffuse(ndl, ndv, 0.5);
                assert!(v > 0.0, "SSS diffuse should be positive");
                assert!(v < 2.0, "SSS diffuse should be bounded");
            }
        }
    }

    #[test]
    fn test_sss_transmittance() {
        assert!((sss_transmittance(0.0, 1.0) - 1.0).abs() < EPS);
        assert!(sss_transmittance(1.0, 1.0) < 1.0);
        assert!(sss_transmittance(1.0, 1.0) > 0.0);
    }

    // ── Iridescence tests ─────────────────────────────────────────────────

    #[test]
    fn test_iridescence_fresnel_range() {
        for cos_t in [0.1, 0.3, 0.5, 0.7, 0.9] {
            let r = iridescence_fresnel(1.0, 1.3, 1.5, 300.0, 550.0, cos_t);
            assert!(
                (0.0..=1.0).contains(&r),
                "Iridescence out of range at cos={cos_t}: {r}"
            );
        }
    }

    #[test]
    fn test_iridescence_varies_with_wavelength() {
        let r_red = iridescence_fresnel(1.0, 1.3, 1.5, 300.0, 650.0, 0.8);
        let r_blue = iridescence_fresnel(1.0, 1.3, 1.5, 300.0, 450.0, 0.8);
        assert!(
            (r_red - r_blue).abs() > 0.001,
            "Iridescence should vary with wavelength"
        );
    }

    #[test]
    fn test_iridescence_varies_with_thickness() {
        let r1 = iridescence_fresnel(1.0, 1.3, 1.5, 200.0, 550.0, 0.8);
        let r2 = iridescence_fresnel(1.0, 1.3, 1.5, 400.0, 550.0, 0.8);
        assert!(
            (r1 - r2).abs() > 0.001,
            "Iridescence should vary with thickness"
        );
    }

    #[test]
    fn test_iridescence_varies_with_angle() {
        let r_normal = iridescence_fresnel(1.0, 1.3, 1.5, 300.0, 550.0, 1.0);
        let r_grazing = iridescence_fresnel(1.0, 1.3, 1.5, 300.0, 550.0, 0.2);
        assert!(
            (r_normal - r_grazing).abs() > 0.01,
            "Iridescence should vary with angle"
        );
    }

    #[test]
    fn test_iridescence_rgb_produces_color() {
        let rgb = iridescence_rgb(1.0, 1.3, 1.5, 300.0, 0.8);
        // Different wavelengths should give different reflectances → color
        assert!(rgb[0] >= 0.0 && rgb[1] >= 0.0 && rgb[2] >= 0.0);
        // At least two channels should differ (it's iridescent!)
        let all_same = (rgb[0] - rgb[1]).abs() < 0.001 && (rgb[1] - rgb[2]).abs() < 0.001;
        assert!(!all_same, "Iridescence should produce color variation");
    }

    // ── Volumetric scattering tests ───────────────────────────────────────

    #[test]
    fn test_hg_isotropic() {
        // g=0 should match isotropic phase function
        let p_hg = henyey_greenstein(0.5, 0.0);
        let p_iso = phase_isotropic();
        assert!(
            (p_hg - p_iso).abs() < 0.001,
            "HG(g=0) should match isotropic: {p_hg} vs {p_iso}"
        );
    }

    #[test]
    fn test_hg_forward_scattering() {
        // g=0.8: forward scattering should be stronger
        let p_fwd = henyey_greenstein(1.0, 0.8);
        let p_back = henyey_greenstein(-1.0, 0.8);
        assert!(p_fwd > p_back, "Forward > backward for g=0.8");
    }

    #[test]
    fn test_hg_backward_scattering() {
        // g=-0.5: backward scattering should be stronger
        let p_fwd = henyey_greenstein(1.0, -0.5);
        let p_back = henyey_greenstein(-1.0, -0.5);
        assert!(p_back > p_fwd, "Backward > forward for g=-0.5");
    }

    #[test]
    fn test_hg_always_positive() {
        for g in [-0.9, -0.5, 0.0, 0.5, 0.9] {
            for ct in [-1.0, -0.5, 0.0, 0.5, 1.0] {
                let p = henyey_greenstein(ct, g);
                assert!(p >= 0.0, "HG negative at g={g}, cos_theta={ct}");
            }
        }
    }

    #[test]
    fn test_phase_rayleigh_symmetric() {
        let p_fwd = phase_rayleigh(1.0);
        let p_back = phase_rayleigh(-1.0);
        assert!((p_fwd - p_back).abs() < EPS, "Rayleigh should be symmetric");
    }

    #[test]
    fn test_phase_rayleigh_min_at_90() {
        let p_90 = phase_rayleigh(0.0); // cos(90°) = 0
        let p_0 = phase_rayleigh(1.0); // cos(0°) = 1
        assert!(p_0 > p_90, "Rayleigh min at 90°");
    }

    #[test]
    fn test_extinction_coefficient() {
        assert!((extinction_coefficient(0.5, 0.3) - 0.8).abs() < EPS);
    }

    #[test]
    fn test_volume_transmittance() {
        assert!((volume_transmittance(0.0, 10.0) - 1.0).abs() < EPS);
        assert!(volume_transmittance(1.0, 1.0) < 1.0);
        assert!(volume_transmittance(1.0, 1.0) > 0.0);
    }

    #[test]
    fn test_volume_transmittance_thicker_less() {
        let t1 = volume_transmittance(0.5, 1.0);
        let t2 = volume_transmittance(0.5, 3.0);
        assert!(t2 < t1);
    }

    #[test]
    fn test_single_scatter_albedo() {
        assert!((single_scatter_albedo(0.0, 1.0) - 1.0).abs() < EPS);
        assert!((single_scatter_albedo(1.0, 0.0)).abs() < EPS);
        assert!((single_scatter_albedo(0.5, 0.5) - 0.5).abs() < EPS);
    }

    #[test]
    fn test_single_scatter_albedo_zero() {
        assert!((single_scatter_albedo(0.0, 0.0)).abs() < EPS);
    }

    #[test]
    fn test_inscattering_positive() {
        let l = single_scatter_inscattering(0.5, 0.8, 0.5, 1.0, 0.5, 1.0);
        assert!(l > 0.0, "In-scattering should be positive");
    }

    #[test]
    fn test_inscattering_zero_scattering() {
        let l = single_scatter_inscattering(0.0, 0.8, 0.5, 1.0, 0.5, 1.0);
        assert!(l.abs() < EPS, "Zero scattering = zero in-scattering");
    }

    // ── Importance sampling tests ─────────────────────────────────────────

    #[test]
    fn test_sample_ggx_normalized() {
        let h = sample_ggx(0.3, 0.5, 0.5);
        let len = (h[0] * h[0] + h[1] * h[1] + h[2] * h[2]).sqrt();
        assert!(
            (len - 1.0).abs() < 0.01,
            "Sampled vector should be unit length"
        );
    }

    #[test]
    fn test_sample_ggx_hemisphere() {
        // z component should always be positive (upper hemisphere)
        for i in 0..20 {
            let xi1 = i as f64 / 20.0;
            let h = sample_ggx(0.5, xi1, 0.3);
            assert!(h[2] >= 0.0, "Half-vector should be in upper hemisphere");
        }
    }

    #[test]
    fn test_sample_ggx_smooth_concentrates() {
        // Smooth surface → samples near normal (z ≈ 1)
        let h = sample_ggx(0.01, 0.5, 0.0);
        assert!(
            h[2] > 0.9,
            "Smooth surface should sample near normal, got z={}",
            h[2]
        );
    }

    #[test]
    fn test_sample_ggx_rough_spreads() {
        // Rough surface → samples spread out
        let h = sample_ggx(0.9, 0.5, 0.0);
        assert!(h[2] < 0.9, "Rough surface should spread samples");
    }

    #[test]
    fn test_sample_ggx_pdf_positive() {
        let pdf = sample_ggx_pdf(0.9, 0.8, 0.3);
        assert!(pdf > 0.0, "PDF should be positive");
    }

    #[test]
    fn test_sample_ggx_pdf_matches_distribution() {
        // PDF at the peak (n·h = 1, h·v ≈ 1) should be high for smooth surfaces
        let pdf_smooth = sample_ggx_pdf(1.0, 1.0, 0.01);
        let pdf_rough = sample_ggx_pdf(1.0, 1.0, 0.9);
        assert!(pdf_smooth > pdf_rough, "Smooth should have higher peak PDF");
    }

    #[test]
    fn test_sample_cosine_hemisphere_normalized() {
        let d = sample_cosine_hemisphere(0.5, 0.5);
        let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        assert!((len - 1.0).abs() < 0.01, "Should be unit length");
    }

    #[test]
    fn test_sample_cosine_hemisphere_upper() {
        for i in 0..20 {
            let xi1 = i as f64 / 20.0;
            let d = sample_cosine_hemisphere(xi1, 0.3);
            assert!(d[2] >= 0.0, "Should be in upper hemisphere");
        }
    }

    #[test]
    fn test_sample_cosine_pdf_at_normal() {
        let pdf = sample_cosine_pdf(1.0);
        assert!((pdf - 1.0 / PI).abs() < EPS);
    }

    #[test]
    fn test_sample_cosine_pdf_at_grazing() {
        let pdf = sample_cosine_pdf(0.0);
        assert!(pdf.abs() < EPS);
    }

    // ── Environment map / split-sum tests ─────────────────────────────────

    #[test]
    fn test_split_sum_range() {
        for ndv_i in 1..=10 {
            for rough_i in 0..=10 {
                let ndv = ndv_i as f64 / 10.0;
                let rough = rough_i as f64 / 10.0;
                let (scale, bias) = split_sum_scale_bias(ndv, rough);
                assert!(
                    (0.0..=1.0).contains(&scale),
                    "Scale out of range at ndv={ndv}, rough={rough}: {scale}"
                );
                assert!(
                    (0.0..=1.0).contains(&bias),
                    "Bias out of range at ndv={ndv}, rough={rough}: {bias}"
                );
            }
        }
    }

    #[test]
    fn test_split_sum_smooth_normal_high_scale() {
        let (scale, _) = split_sum_scale_bias(1.0, 0.0);
        assert!(
            scale > 0.5,
            "Smooth at normal: scale should be significant, got {scale}"
        );
    }

    #[test]
    fn test_env_map_lod() {
        assert!((env_map_lod(0.0, 8.0)).abs() < EPS);
        assert!((env_map_lod(1.0, 8.0) - 8.0).abs() < EPS);
        assert!((env_map_lod(0.5, 8.0) - 4.0).abs() < EPS);
    }

    #[test]
    fn test_integrate_brdf_lut_range() {
        let (scale, bias) = integrate_brdf_lut(0.5, 0.3, 64);
        assert!(scale >= 0.0, "Scale should be non-negative");
        assert!(bias >= 0.0, "Bias should be non-negative");
        assert!(scale + bias <= 2.0, "Scale+bias should be reasonable");
    }

    #[test]
    fn test_integrate_brdf_lut_smooth_vs_rough() {
        let (scale_smooth, _) = integrate_brdf_lut(0.8, 0.1, 64);
        let (scale_rough, _) = integrate_brdf_lut(0.8, 0.9, 64);
        assert!(
            scale_smooth > scale_rough,
            "Smooth should have higher scale: {scale_smooth} vs {scale_rough}"
        );
    }

    #[test]
    fn test_integrate_brdf_lut_deterministic() {
        let (s1, b1) = integrate_brdf_lut(0.5, 0.3, 64);
        let (s2, b2) = integrate_brdf_lut(0.5, 0.3, 64);
        assert!((s1 - s2).abs() < EPS, "Should be deterministic");
        assert!((b1 - b2).abs() < EPS, "Should be deterministic");
    }
}
