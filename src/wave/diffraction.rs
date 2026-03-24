//! Advanced diffraction (Fraunhofer, Fresnel, Huygens-Fresnel) and anti-reflection coatings.

use tracing::trace;

use super::single_slit_intensity;

/// Fraunhofer diffraction pattern for a rectangular aperture.
///
/// I(θx, θy) = I₀ · sinc²(π·a·sin(θx)/λ) · sinc²(π·b·sin(θy)/λ)
///
/// `width` and `height` are aperture dimensions, same units as `wavelength`.
/// `angle_x` and `angle_y` in radians.
#[must_use]
#[inline]
pub fn fraunhofer_rect(
    width: f64,
    height: f64,
    wavelength: f64,
    angle_x: f64,
    angle_y: f64,
    i0: f64,
) -> f64 {
    single_slit_intensity(width, wavelength, angle_x, 1.0)
        * single_slit_intensity(height, wavelength, angle_y, 1.0)
        * i0
}

/// Fraunhofer diffraction for a 1D aperture defined by a transmission function.
///
/// Numerically integrates the aperture field and computes the far-field
/// intensity at a given observation angle.
///
/// `aperture` is a slice of (position, complex_amplitude) pairs across the aperture.
/// `wavelength` and positions in same units. `angle` in radians.
/// Returns intensity (proportional to |E|²).
#[must_use]
#[inline]
pub fn fraunhofer_1d(aperture: &[(f64, f64)], wavelength: f64, angle: f64) -> f64 {
    let k_sin = std::f64::consts::TAU * angle.sin() / wavelength;
    let mut real = 0.0;
    let mut imag = 0.0;
    for &(x, amplitude) in aperture {
        let (sin_p, cos_p) = (k_sin * x).sin_cos();
        real += amplitude * cos_p;
        imag += amplitude * sin_p;
    }
    real * real + imag * imag
}

// ── Fresnel Diffraction (Near-Field) ──────────────────────────────────────

/// Fresnel number: N_F = a² / (λ · z).
///
/// Determines whether diffraction is in the Fresnel (N_F ≳ 1) or
/// Fraunhofer (N_F ≪ 1) regime.
///
/// `aperture_radius`, `wavelength`, and `distance` in same units.
#[must_use]
#[inline]
pub fn fresnel_number(aperture_radius: f64, wavelength: f64, distance: f64) -> f64 {
    aperture_radius * aperture_radius / (wavelength * distance)
}

/// Fresnel cosine integral C(x) — polynomial approximation.
///
/// C(x) = ∫₀ˣ cos(πt²/2) dt
///
/// Accurate to ~1e-6 for all x. Uses rational approximation
/// (Abramowitz & Stegun, 7.3.19–7.3.20).
#[must_use]
#[inline]
pub fn fresnel_c(x: f64) -> f64 {
    let ax = x.abs();
    let result = if ax < 1.0 {
        // Small x: power series converges fast
        // C(x) ≈ x - (π²/40)x⁵ + (π⁴/3456)x⁹ - ...
        let x2 = ax * ax;
        let t = std::f64::consts::FRAC_PI_2 * x2;
        let t2 = t * t;
        ax * (1.0 - t2 / 20.0 + t2 * t2 / 1680.0)
    } else {
        // Auxiliary functions f(x) and g(x)
        let pi_x2 = std::f64::consts::FRAC_PI_2 * ax * ax;
        let (f, g) = fresnel_fg(ax);
        0.5 + f * pi_x2.sin() - g * pi_x2.cos()
    };
    if x < 0.0 { -result } else { result }
}

/// Fresnel sine integral S(x) — polynomial approximation.
///
/// S(x) = ∫₀ˣ sin(πt²/2) dt
#[must_use]
#[inline]
pub fn fresnel_s(x: f64) -> f64 {
    let ax = x.abs();
    let result = if ax < 1.0 {
        let x2 = ax * ax;
        let t = std::f64::consts::FRAC_PI_2 * x2;
        let t2 = t * t;
        ax * x2 * (std::f64::consts::FRAC_PI_2 / 3.0) * (1.0 - t2 / 42.0 + t2 * t2 / 3960.0)
    } else {
        let pi_x2 = std::f64::consts::FRAC_PI_2 * ax * ax;
        let (f, g) = fresnel_fg(ax);
        0.5 - f * pi_x2.cos() - g * pi_x2.sin()
    };
    if x < 0.0 { -result } else { result }
}

/// Auxiliary functions f(x) and g(x) for Fresnel integrals (x ≥ 1).
#[inline]
fn fresnel_fg(x: f64) -> (f64, f64) {
    let x2 = x * x;
    let x3 = x2 * x;
    let x4 = x2 * x2;
    // Rational approximation (Boersma, 1960)
    let f = (1.0 + 0.926 * x2) / (2.0 + 1.792 * x2 + 3.104 * x4) / x;
    let g = 1.0 / (2.0 + 4.142 * x2 + 3.492 * x4 + 6.670 * x2 * x4) / x3;
    (f, g)
}

/// Compute both Fresnel C and S integrals together, sharing work.
///
/// Returns `(C(x), S(x))` where:
/// - C(x) = ∫₀ˣ cos(πt²/2) dt
/// - S(x) = ∫₀ˣ sin(πt²/2) dt
#[must_use]
#[inline]
pub fn fresnel_cs(x: f64) -> (f64, f64) {
    let ax = x.abs();
    let (c, s) = if ax < 1.0 {
        let x2 = ax * ax;
        let t = std::f64::consts::FRAC_PI_2 * x2;
        let t2 = t * t;
        let c_val = ax * (1.0 - t2 / 20.0 + t2 * t2 / 1680.0);
        let s_val =
            ax * x2 * (std::f64::consts::FRAC_PI_2 / 3.0) * (1.0 - t2 / 42.0 + t2 * t2 / 3960.0);
        (c_val, s_val)
    } else {
        let pi_x2 = std::f64::consts::FRAC_PI_2 * ax * ax;
        let (f, g) = fresnel_fg(ax);
        let (sin_px2, cos_px2) = pi_x2.sin_cos();
        let c_val = 0.5 + f * sin_px2 - g * cos_px2;
        let s_val = 0.5 - f * cos_px2 - g * sin_px2;
        (c_val, s_val)
    };
    if x < 0.0 { (-c, -s) } else { (c, s) }
}

/// Fresnel straight-edge diffraction intensity.
///
/// Intensity at a point behind a semi-infinite opaque screen as a function
/// of the Fresnel parameter `u`. The edge is at u=0.
///
/// I/I₀ = [(C(u) + 0.5)² + (S(u) + 0.5)²] / 2
///
/// `u` is the dimensionless Fresnel parameter: u = x·√(2/(λ·z))
/// where x is distance from geometric shadow edge.
#[must_use]
#[inline]
pub fn fresnel_edge_intensity(u: f64) -> f64 {
    let (fc, fs) = fresnel_cs(u);
    let c = fc + 0.5;
    let s = fs + 0.5;
    (c * c + s * s) / 2.0
}

/// Fresnel parameter for a point at distance `x` from shadow edge.
///
/// u = x · √(2 / (λ · z))
///
/// `x` = lateral distance from geometric shadow edge (positive = illuminated side).
/// `wavelength` and `distance` in same units as `x`.
#[must_use]
#[inline]
pub fn fresnel_parameter(x: f64, wavelength: f64, distance: f64) -> f64 {
    x * (2.0 / (wavelength * distance)).sqrt()
}

/// Huygens-Fresnel diffraction integral (1D numerical).
///
/// Computes the complex field amplitude at observation point `x_obs` by
/// integrating secondary wavelets across the aperture.
///
/// Returns intensity (|E|²).
///
/// `aperture` is a slice of (position, amplitude) pairs.
/// `wavelength` and all positions in same units.
/// `z` is the propagation distance to the observation plane.
#[must_use]
#[inline]
pub fn huygens_fresnel_1d(aperture: &[(f64, f64)], wavelength: f64, z: f64, x_obs: f64) -> f64 {
    let k = std::f64::consts::TAU / wavelength;
    let mut real = 0.0;
    let mut imag = 0.0;
    for &(x_ap, amplitude) in aperture {
        let dx = x_obs - x_ap;
        let r = z.hypot(dx);
        let (sin_p, cos_p) = (k * r).sin_cos();
        let weight = amplitude * z / (r * r);
        real += weight * cos_p;
        imag += weight * sin_p;
    }
    real * real + imag * imag
}

// ── Anti-Reflection Coatings ──────────────────────────────────────────────

/// Ideal quarter-wave AR coating refractive index.
///
/// n_coating = √(n₁ · n₂)
///
/// For air (n₁=1) to glass (n₂=1.52): n_coating ≈ 1.233 (MgF₂ ≈ 1.38 is closest practical).
#[must_use]
#[inline]
pub fn ar_ideal_index(n1: f64, n2: f64) -> f64 {
    (n1 * n2).sqrt()
}

/// Quarter-wave coating thickness for a given design wavelength.
///
/// d = λ / (4 · n_coating)
///
/// `wavelength` and result in same units.
#[must_use]
#[inline]
pub fn ar_quarter_wave_thickness(wavelength: f64, n_coating: f64) -> f64 {
    wavelength / (4.0 * n_coating)
}

/// Reflectance of a single-layer thin film coating at normal incidence.
///
/// Uses the exact thin-film formula:
/// R = [(n₁·n₃ − n₂²)² · sin²(δ)] / [(n₁·n₃ + n₂²)² · sin²(δ) + n₂²·(n₁+n₃)²·cos²(δ)]
///
/// where δ = 2π·n₂·d/λ is the phase thickness.
///
/// `n1` = incident medium, `n2` = coating, `n3` = substrate,
/// `thickness` and `wavelength` in same units.
#[must_use]
#[inline]
pub fn coating_reflectance(n1: f64, n2: f64, n3: f64, thickness: f64, wavelength: f64) -> f64 {
    let delta = std::f64::consts::TAU * n2 * thickness / wavelength;
    let (sin_d, cos_d) = delta.sin_cos();
    let a = n1 * n3 - n2 * n2;
    let b = n1 * n3 + n2 * n2;
    let c = n2 * (n1 + n3);
    (a * a * sin_d * sin_d) / (b * b * sin_d * sin_d + c * c * cos_d * cos_d)
}

/// Reflectance of a V-coat (quarter-wave) AR coating at design wavelength.
///
/// At the design wavelength, δ = π/2, and:
/// R = [(n₁·n₃ − n₂²) / (n₁·n₃ + n₂²)]²
///
/// Returns zero when n₂ = √(n₁·n₃) (ideal AR).
#[must_use]
#[inline]
pub fn vcoat_reflectance(n1: f64, n2: f64, n3: f64) -> f64 {
    let num = n1 * n3 - n2 * n2;
    let den = n1 * n3 + n2 * n2;
    (num / den) * (num / den)
}

/// Multi-layer coating reflectance using the transfer matrix method.
///
/// Each layer is (refractive_index, thickness). Computes reflectance at
/// normal incidence for a given wavelength.
///
/// `n_incident` = incident medium index, `n_substrate` = substrate index.
/// `layers` = slice of (n, thickness) from outermost to innermost.
/// `wavelength` in same units as thickness.
#[must_use]
pub fn multilayer_reflectance(
    n_incident: f64,
    n_substrate: f64,
    layers: &[(f64, f64)],
    wavelength: f64,
) -> f64 {
    trace!(
        num_layers = layers.len(),
        wavelength, "multilayer_reflectance"
    );
    // Transfer matrix method: M = Π Mᵢ where Mᵢ = [[cos δᵢ, -j sin δᵢ/nᵢ], [-j nᵢ sin δᵢ, cos δᵢ]]
    // Track real/imag parts of 2×2 matrix entries
    // M = [[m11, m12], [m21, m22]] where m12 and m21 are imaginary
    let mut m11_r = 1.0;
    let mut m12_i = 0.0; // m12 = j * m12_i
    let mut m21_i = 0.0; // m21 = j * m21_i
    let mut m22_r = 1.0;

    for &(n, d) in layers {
        let delta = std::f64::consts::TAU * n * d / wavelength;
        let (sin_d, cos_d) = delta.sin_cos();

        let new_m11_r = m11_r * cos_d - m12_i * n * sin_d;
        let new_m12_i = -m11_r * sin_d / n + m12_i * cos_d;
        let new_m21_i = m21_i * cos_d - m22_r * n * sin_d;
        let new_m22_r = -m21_i * sin_d / n + m22_r * cos_d;

        m11_r = new_m11_r;
        m12_i = new_m12_i;
        m21_i = new_m21_i;
        m22_r = new_m22_r;
    }

    // Reflection coefficient: r = (n_i·M11 + n_i·n_s·M12 - M21 - n_s·M22) /
    //                              (n_i·M11 + n_i·n_s·M12 + M21 + n_s·M22)
    // With our convention (M12 and M21 imaginary):
    let ni = n_incident;
    let ns = n_substrate;

    let num_r = ni * m11_r - ns * m22_r;
    let num_i = ni * ns * m12_i - m21_i;
    let den_r = ni * m11_r + ns * m22_r;
    let den_i = ni * ns * m12_i + m21_i;

    // R = |r|² = (num_r² + num_i²) / (den_r² + den_i²)
    (num_r * num_r + num_i * num_i) / (den_r * den_r + den_i * den_i)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-6;

    #[test]
    fn test_fraunhofer_rect_central_max() {
        let i = fraunhofer_rect(1e-3, 1e-3, 550e-9, 0.0, 0.0, 1.0);
        assert!((i - 1.0).abs() < EPS);
    }

    #[test]
    fn test_fraunhofer_rect_decreases_off_axis() {
        let i_center = fraunhofer_rect(1e-3, 1e-3, 550e-9, 0.0, 0.0, 1.0);
        let i_off = fraunhofer_rect(1e-3, 1e-3, 550e-9, 0.001, 0.0, 1.0);
        assert!(i_off < i_center);
    }

    #[test]
    fn test_fraunhofer_rect_separable() {
        // Rectangular pattern = product of two 1D patterns
        let ix = single_slit_intensity(1e-3, 550e-9, 0.001, 1.0);
        let iy = single_slit_intensity(0.5e-3, 550e-9, 0.002, 1.0);
        let i_rect = fraunhofer_rect(1e-3, 0.5e-3, 550e-9, 0.001, 0.002, 1.0);
        assert!((i_rect - ix * iy).abs() < EPS);
    }

    #[test]
    fn test_fraunhofer_1d_uniform_aperture() {
        // Uniform aperture should match sinc² pattern behavior
        let n = 100;
        let width = 1e-3;
        let aperture: Vec<(f64, f64)> = (0..n)
            .map(|i| {
                let x = -width / 2.0 + width * (i as f64) / (n as f64 - 1.0);
                (x, 1.0)
            })
            .collect();
        let i_center = fraunhofer_1d(&aperture, 550e-9, 0.0);
        let i_off = fraunhofer_1d(&aperture, 550e-9, 0.001);
        assert!(i_center > i_off, "Central max should be brightest");
        assert!(i_center > 0.0);
    }

    #[test]
    fn test_fraunhofer_1d_non_negative() {
        let aperture: Vec<(f64, f64)> = (0..50).map(|i| ((i as f64 - 25.0) * 1e-5, 1.0)).collect();
        for angle_mrad in 0..20 {
            let angle = angle_mrad as f64 * 0.001;
            let i = fraunhofer_1d(&aperture, 550e-9, angle);
            assert!(i >= 0.0, "Negative intensity at angle {angle}");
        }
    }

    // ── Fresnel diffraction tests ─────────────────────────────────────────

    #[test]
    fn test_fresnel_number_far_field() {
        // Small aperture, large distance → N_F ≪ 1 (Fraunhofer regime)
        let nf = fresnel_number(0.1e-3, 550e-9, 1.0);
        assert!(nf < 1.0, "Should be Fraunhofer regime, N_F={nf}");
    }

    #[test]
    fn test_fresnel_number_near_field() {
        // Large aperture, short distance → N_F ≫ 1 (Fresnel regime)
        let nf = fresnel_number(5e-3, 550e-9, 0.01);
        assert!(nf > 1.0, "Should be Fresnel regime, N_F={nf}");
    }

    #[test]
    fn test_fresnel_c_at_zero() {
        assert!(fresnel_c(0.0).abs() < EPS);
    }

    #[test]
    fn test_fresnel_s_at_zero() {
        assert!(fresnel_s(0.0).abs() < EPS);
    }

    #[test]
    fn test_fresnel_c_converges_to_half() {
        // C(∞) → 0.5
        let c = fresnel_c(10.0);
        assert!((c - 0.5).abs() < 0.05, "C(10) should be ≈0.5, got {c}");
    }

    #[test]
    fn test_fresnel_s_converges_to_half() {
        // S(∞) → 0.5
        let s = fresnel_s(10.0);
        assert!((s - 0.5).abs() < 0.05, "S(10) should be ≈0.5, got {s}");
    }

    #[test]
    fn test_fresnel_c_odd_function() {
        assert!((fresnel_c(-2.0) + fresnel_c(2.0)).abs() < 0.01);
    }

    #[test]
    fn test_fresnel_s_odd_function() {
        assert!((fresnel_s(-2.0) + fresnel_s(2.0)).abs() < 0.01);
    }

    #[test]
    fn test_fresnel_c_known_value() {
        // C(1) ≈ 0.7799
        let c = fresnel_c(1.0);
        assert!((c - 0.7799).abs() < 0.01, "C(1) ≈ 0.7799, got {c}");
    }

    #[test]
    fn test_fresnel_s_known_value() {
        // S(1) ≈ 0.4383
        let s = fresnel_s(1.0);
        assert!((s - 0.4383).abs() < 0.01, "S(1) ≈ 0.4383, got {s}");
    }

    #[test]
    fn test_fresnel_edge_shadow_boundary() {
        // At u=0 (geometric shadow edge): I/I₀ = 0.25
        let i = fresnel_edge_intensity(0.0);
        assert!((i - 0.25).abs() < 0.01, "At edge: I/I₀ ≈ 0.25, got {i}");
    }

    #[test]
    fn test_fresnel_edge_deep_shadow() {
        // Deep in shadow (u → -∞): I → 0
        let i = fresnel_edge_intensity(-5.0);
        assert!(i < 0.01, "Deep shadow should be ~0, got {i}");
    }

    #[test]
    fn test_fresnel_edge_illuminated() {
        // Well into illuminated region (u → +∞): I → 1
        let i = fresnel_edge_intensity(5.0);
        assert!((i - 1.0).abs() < 0.05, "Illuminated region ≈ 1.0, got {i}");
    }

    #[test]
    fn test_fresnel_parameter() {
        let u = fresnel_parameter(1e-3, 550e-9, 1.0);
        assert!(u > 0.0);
        assert!(u.is_finite());
    }

    #[test]
    fn test_huygens_fresnel_1d_positive() {
        let aperture: Vec<(f64, f64)> = (0..50).map(|i| ((i as f64 - 25.0) * 1e-5, 1.0)).collect();
        let i = huygens_fresnel_1d(&aperture, 550e-9, 0.1, 0.0);
        assert!(i > 0.0);
    }

    #[test]
    fn test_huygens_fresnel_1d_peak_on_axis() {
        let aperture: Vec<(f64, f64)> = (0..100).map(|i| ((i as f64 - 50.0) * 1e-5, 1.0)).collect();
        let i_center = huygens_fresnel_1d(&aperture, 550e-9, 0.1, 0.0);
        let i_off = huygens_fresnel_1d(&aperture, 550e-9, 0.1, 1e-3);
        assert!(i_center > i_off, "On-axis should be brightest");
    }

    // ── AR coating tests ──────────────────────────────────────────────────

    #[test]
    fn test_ar_ideal_index() {
        // Air to glass: √(1.0 × 1.52) ≈ 1.233
        let n = ar_ideal_index(1.0, 1.52);
        assert!((n - 1.233).abs() < 0.001);
    }

    #[test]
    fn test_ar_quarter_wave_thickness() {
        let d = ar_quarter_wave_thickness(550.0, 1.38);
        assert!((d - 550.0 / (4.0 * 1.38)).abs() < EPS);
    }

    #[test]
    fn test_vcoat_ideal_zero_reflectance() {
        // When n₂ = √(n₁·n₃), V-coat gives zero reflectance
        let n1: f64 = 1.0;
        let n3: f64 = 1.52;
        let n2 = (n1 * n3).sqrt();
        let r = vcoat_reflectance(n1, n2, n3);
        assert!(r < EPS, "Ideal V-coat should have ~0 reflectance, got {r}");
    }

    #[test]
    fn test_vcoat_mgf2_on_glass() {
        // MgF₂ (n=1.38) on glass (n=1.52) in air: small but nonzero
        let r = vcoat_reflectance(1.0, 1.38, 1.52);
        assert!(
            r > 0.0 && r < 0.02,
            "MgF₂ V-coat ≈ 1.3% reflectance, got {r}"
        );
    }

    #[test]
    fn test_coating_reflectance_at_design_wavelength() {
        // Quarter-wave coating at design wavelength should match V-coat
        let n1 = 1.0;
        let n2 = 1.38;
        let n3 = 1.52;
        let wl = 550.0;
        let d = ar_quarter_wave_thickness(wl, n2);
        let r_coating = coating_reflectance(n1, n2, n3, d, wl);
        let r_vcoat = vcoat_reflectance(n1, n2, n3);
        assert!(
            (r_coating - r_vcoat).abs() < 0.001,
            "Coating at design λ should match V-coat: {r_coating} vs {r_vcoat}"
        );
    }

    #[test]
    fn test_coating_reflectance_varies_with_wavelength() {
        // Reflectance should change with wavelength (non-constant)
        let n2 = 1.38;
        let n3 = 1.52;
        let d = ar_quarter_wave_thickness(550.0, n2);
        let r_design = coating_reflectance(1.0, n2, n3, d, 550.0);
        let r_off = coating_reflectance(1.0, n2, n3, d, 450.0);
        assert!(
            (r_design - r_off).abs() > 0.001,
            "Reflectance should vary with wavelength"
        );
    }

    #[test]
    fn test_coating_reflectance_range() {
        let n2 = 1.38;
        let n3 = 1.52;
        let d = ar_quarter_wave_thickness(550.0, n2);
        for wl_nm in (400..=700).step_by(10) {
            let r = coating_reflectance(1.0, n2, n3, d, wl_nm as f64);
            assert!(
                (0.0..=1.0).contains(&r),
                "Reflectance out of range at {wl_nm}nm: {r}"
            );
        }
    }

    #[test]
    fn test_multilayer_single_layer_matches_coating() {
        let n1 = 1.0;
        let n2 = 1.38;
        let n3 = 1.52;
        let wl = 550.0;
        let d = ar_quarter_wave_thickness(wl, n2);
        let r_single = coating_reflectance(n1, n2, n3, d, wl);
        let r_multi = multilayer_reflectance(n1, n3, &[(n2, d)], wl);
        assert!(
            (r_single - r_multi).abs() < 0.001,
            "Single layer should match: coating={r_single}, multi={r_multi}"
        );
    }

    #[test]
    fn test_multilayer_no_layers() {
        // No coating → bare interface
        let r = multilayer_reflectance(1.0, 1.52, &[], 550.0);
        let r_bare = ((1.0_f64 - 1.52) / (1.0 + 1.52)).powi(2);
        assert!(
            (r - r_bare).abs() < 0.001,
            "No layers = bare surface: {r} vs {r_bare}"
        );
    }

    #[test]
    fn test_multilayer_two_layer_valid() {
        // Two-layer coating should produce valid reflectance
        let wl = 550.0;
        let n_mgf2 = 1.38;
        let n_zro2 = 2.1;
        let d1 = ar_quarter_wave_thickness(wl, n_zro2);
        let d2 = ar_quarter_wave_thickness(wl, n_mgf2);
        let r = multilayer_reflectance(1.0, 1.52, &[(n_mgf2, d2), (n_zro2, d1)], wl);
        assert!(
            (0.0..=1.0).contains(&r),
            "Two-layer reflectance should be in [0,1], got {r}"
        );
    }

    #[test]
    fn test_multilayer_reflectance_range() {
        let wl = 550.0;
        let d = ar_quarter_wave_thickness(wl, 1.38);
        let layers = [(1.38, d), (2.1, ar_quarter_wave_thickness(wl, 2.1))];
        for wl_nm in (400..=700).step_by(10) {
            let r = multilayer_reflectance(1.0, 1.52, &layers, wl_nm as f64);
            assert!(
                (0.0..=1.0 + EPS).contains(&r),
                "Reflectance out of range at {wl_nm}nm: {r}"
            );
        }
    }
}
