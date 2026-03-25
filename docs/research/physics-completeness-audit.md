# Prakash Physics Completeness Audit

> Deep external research report -- 2026-03-25
>
> Goal: identify every gap between prakash's current implementation and
> what a world-class, 100%-accurate optics library must contain.
>
> **Implementation status** (updated 2026-03-25):
> - Section 1 (Complex Fresnel): **IMPLEMENTED** -- `ComplexMedium`, `fresnel_s_complex`, `fresnel_p_complex`, metal presets
> - Section 10 (Planck stability): **IMPLEMENTED** -- `exp_m1()`, Wien branch, T<=0 guard
> - Section 13 (Polarization ray tracing): **IMPLEMENTED** -- `trace_sequential_polarized` with s/p tracking
> - Section 15 (Zernike polynomials): **IMPLEMENTED** -- `wave::zernike` module, Noll indexing, wavefront analysis

---

## Table of Contents

1. [Fresnel Equations: Complex Refractive Index](#1-fresnel-equations-complex-refractive-index)
2. [Dispersion Models](#2-dispersion-models)
3. [Thin Film / Multilayer Coatings](#3-thin-film--multilayer-coatings)
4. [Aberrations: Chromatic and Higher-Order](#4-aberrations-chromatic-and-higher-order)
5. [MTF Completeness](#5-mtf-completeness)
6. [Gaussian Beam & Higher-Order Modes](#6-gaussian-beam--higher-order-modes)
7. [Fiber Optics](#7-fiber-optics)
8. [Photometry vs Radiometry](#8-photometry-vs-radiometry)
9. [CIE Observer Accuracy](#9-cie-observer-accuracy)
10. [Planck Function Numerical Stability](#10-planck-function-numerical-stability)
11. [Rayleigh Scattering: King Correction Factor](#11-rayleigh-scattering-king-correction-factor)
12. [Mie Scattering: Cornette-Shanks Limitations](#12-mie-scattering-cornette-shanks-limitations)
13. [Polarization Ray Tracing](#13-polarization-ray-tracing)
14. [Vectorial Diffraction (Richards-Wolf)](#14-vectorial-diffraction-richards-wolf)
15. [Zernike Polynomials](#15-zernike-polynomials)
16. [GRIN Optics](#16-grin-optics)
17. [Diffractive Optical Elements](#17-diffractive-optical-elements)
18. [Fluorescence](#18-fluorescence)
19. [Nonlinear Optics](#19-nonlinear-optics)
20. [Orbital Angular Momentum](#20-orbital-angular-momentum)
21. [Competitive Analysis](#21-competitive-analysis)
22. [P2 Roadmap Item Assessment](#22-p2-roadmap-item-assessment)
23. [Prioritized Recommendations](#23-prioritized-recommendations)

---

## 1. Fresnel Equations: Complex Refractive Index

### Current state

`ray::fresnel_s` and `ray::fresnel_p` accept real `n1`, `n2` only. The
`pbr::fresnel_schlick` is an approximation for real-time rendering. Neither
handles absorbing media (metals, semiconductors).

### What's missing

For metals and semiconductors, the refractive index is complex:
**n~ = n + ik** where k is the extinction coefficient. The full Fresnel
equations for absorbing media are:

```
r_s = (n1 cos(theta_i) - n2~ cos(theta_t)) / (n1 cos(theta_i) + n2~ cos(theta_t))
r_p = (n2~ cos(theta_i) - n1 cos(theta_t)) / (n2~ cos(theta_i) + n1 cos(theta_t))
```

where `n2~ = n2 + ik2` and `cos(theta_t)` is computed from the complex
Snell's law:

```
cos(theta_t) = sqrt(1 - (n1/n2~)^2 * sin^2(theta_i))
```

The reflectances `R_s = |r_s|^2` and `R_p = |r_p|^2` are then real.

### Impact

- **HIGH**. Without complex n, prakash cannot accurately model: gold, silver,
  copper, aluminum mirrors, semiconductor wafers, metal-dielectric interfaces
  in thin films, or any PBR metallic material with correct physics (as opposed
  to Schlick approximation).

### Recommendation

Add `fresnel_s_complex(n1: f64, n2_real: f64, n2_imag: f64, angle: f64) -> f64`
and the p-polarized counterpart. This is the single most important accuracy
gap. Consider a `ComplexMedium { n: f64, k: f64 }` type.

### References

- [RP Photonics: Fresnel equations](https://www.rp-photonics.com/fresnel_equations.html)
- [Wikipedia: Fresnel equations](https://en.wikipedia.org/wiki/Fresnel_equations)
- [Mitsuba3 complex Fresnel issue #1113](https://github.com/mitsuba-renderer/mitsuba3/issues/1113)

---

## 2. Dispersion Models

### Current state

Prakash implements:
- **Cauchy** (2-term: `n = b + c/lambda^2`)
- **Sellmeier** (3-term: standard form)

### What's missing

Industry-standard optical design uses at least 11 dispersion formulas.
The important ones prakash lacks:

| Model | Formula | Use case |
|-------|---------|----------|
| **Schott** (6-term) | `n^2 = a0 + a1*l^2 + a2*l^-2 + a3*l^-4 + a4*l^-6 + a5*l^-8` | Legacy glass catalogs (pre-1992) |
| **Herzberger** (5-term) | `n = A + B*L + C*L^2 + D*l^2 + E*l^4` where `L = 1/(l^2 - 0.028)` | IR materials, wide spectral range |
| **Conrady** (3-term) | `n = n0 + A/l + B/l^3.5` | Quick visible-range fit |
| **Extended Cauchy** (6-term) | Adds `d/l^4`, `e/l^6`, etc. | Better than 2-term Cauchy |
| **Laurent** | `n^2 = A + B/l^2 + C/l^4 + ...` | Alternative to Sellmeier |

**Schott abandoned their eponymous formula in 1992 in favor of Sellmeier.**
However, legacy glass catalogs and refractiveindex.info still use Schott
format. The Herzberger formula is important for IR materials.

### Impact

- **MEDIUM**. Sellmeier covers 95% of use cases for modern glass catalogs.
  But for IR materials, legacy data, and maximum accuracy, additional models
  matter.

### Recommendation

1. Add `HerzbergerCoefficients` for IR material support.
2. Add `SchottCoefficients` for legacy glass catalog compatibility.
3. Add `ConradyCoefficients` for quick fits.
4. Consider a generic `DispersionModel` enum that wraps all of them.

### References

- [Arizona Optics: Dispersion Equations](https://wp.optics.arizona.edu/jpalmer/wp-content/uploads/sites/65/2017/03/dispeqns.pdf)
- [Zemax: Fitting index data](https://support.zemax.com/hc/en-us/articles/1500005487941)
- [Sellmeier equation -- Wikipedia](https://en.wikipedia.org/wiki/Sellmeier_equation)

---

## 3. Thin Film / Multilayer Coatings

### Current state

Prakash has:
- `thin_film_reflectance()` -- simplified sin^2 model, normal incidence only
- `coating_reflectance()` -- single-layer exact formula, normal incidence
- `multilayer_reflectance()` -- transfer matrix method, normal incidence,
  real refractive indices only

### What's missing

A world-class thin-film model requires:

1. **Oblique incidence** -- the transfer matrix changes for non-zero angle;
   the phase thickness becomes `delta = 2*pi*n*d*cos(theta)/lambda` and the
   effective admittance differs for s- and p-polarization.

2. **Complex refractive indices** in the TMM -- metals and absorbing layers
   require complex n (ties back to Section 1). The current implementation
   only stores `(f64, f64)` for layers -- needs `(Complex, f64)`.

3. **Separate s/p polarization results** -- the current TMM returns a single
   reflectance, but for polarized light the s and p reflectances differ.

4. **Transmittance** -- the current TMM only returns R, not T. For coatings
   analysis, both R and T (and absorptance A = 1 - R - T) are needed.

5. **Electric field distribution** through the stack -- for damage threshold
   analysis, knowing where the E-field peaks inside the multilayer is critical.

### Impact

- **HIGH**. The current TMM is correct for normal incidence on dielectric
  stacks, but cannot handle metallic layers, oblique light, or polarization.
  These are standard requirements for coating design.

### Recommendation

Extend `multilayer_reflectance` to:
```rust
pub fn multilayer_rt(
    wavelength: f64,
    angle: f64,          // incidence angle
    n_incident: Complex,
    n_substrate: Complex,
    layers: &[(Complex, f64)],  // (complex n, thickness)
) -> ThinFilmResult { R_s, R_p, T_s, T_p, A_s, A_p }
```

### References

- [Transfer-matrix method (optics) -- Wikipedia](https://en.wikipedia.org/wiki/Transfer-matrix_method_(optics))
- [Matrix method for thin film optics (arXiv)](https://arxiv.org/pdf/1809.07708)

---

## 4. Aberrations: Chromatic and Higher-Order

### Current state

Prakash has:
- **Seidel (3rd-order) monochromatic**: spherical, coma, astigmatism, field
  curvature, distortion
- **Longitudinal chromatic aberration**: `LCA = f/V`
- **Petzval sum/radius**

### What's missing

1. **Lateral (transverse) chromatic aberration** -- the difference in image
   height for different wavelengths. Formula:
   ```
   TCA = h * (1/V) * (y_bar / y)
   ```
   where h is image height, y_bar/y is the ratio of chief ray to marginal ray
   height. This is distinct from LCA and critical for wide-field imaging.

2. **Secondary spectrum** -- residual chromatic aberration after achromatic
   correction. An achromat brings F and C lines together, but the d-line
   focus remains offset. The secondary spectrum is:
   ```
   delta_f = f * P / V
   ```
   where P is the partial dispersion ratio `P = (n_d - n_C) / (n_F - n_C)`.
   For BK7, P ~ 0.3.

3. **Higher-order (5th, 7th) monochromatic aberrations** -- the Seidel sums
   are 3rd-order only. For high-performance systems, 5th-order (Buchdahl)
   and even 7th-order terms matter.

4. **Wavefront aberration coefficients** -- the Seidel sums give transverse
   ray aberration. Converting to wavefront (OPD) coefficients requires the
   standard mapping (W_040 for spherical, W_131 for coma, etc.).

### Impact

- **MEDIUM-HIGH**. Missing TCA and secondary spectrum means prakash can't
  fully analyze achromatic/apochromatic lens design. Higher-order aberrations
  are advanced but expected in a world-class library.

### Recommendation

1. Add `lateral_chromatic_aberration()`.
2. Add `secondary_spectrum()` using partial dispersion.
3. Add `partial_dispersion(sellmeier, lambda_short, lambda_long, lambda_ref)`.
4. Consider Buchdahl 5th-order coefficients as a stretch goal.

### References

- [Chromatic Aberrations -- Arizona](https://wp.optics.arizona.edu/jsasian/wp-content/uploads/sites/33/2016/03/L9_OPTI517_Chromatic_Aberrations.pdf)
- [Chromatic aberration -- Wikipedia](https://en.wikipedia.org/wiki/Chromatic_aberration)
- [Edmund Optics: Chromatic and Monochromatic Aberrations](https://www.edmundoptics.com/knowledge-center/application-notes/optics/chromatic-and-monochromatic-optical-aberrations/)

---

## 5. MTF Completeness

### Current state

- `mtf_cutoff_frequency()` -- monochromatic
- `mtf_diffraction_limited()` -- monochromatic, circular aperture

### What's missing

1. **Polychromatic MTF** -- computed by weighted vector sum of monochromatic
   MTFs at multiple wavelengths:
   ```
   MTF_poly(v) = |sum_i w_i * MTF_i(v) * exp(j * Phase_i(v))| / sum_i w_i
   ```
   where w_i are spectral weights (e.g., from V(lambda) or detector response).

2. **Through-focus MTF** -- MTF as a function of defocus at a fixed spatial
   frequency. This is the standard tool for analyzing depth of focus and
   chromatic focal shift.

3. **MTF with aberrations** -- the current model is diffraction-limited only.
   Aberrated MTF requires computing the OTF from the generalized pupil
   function via autocorrelation.

4. **Sagittal/Tangential MTF** -- off-axis, MTF differs in the sagittal and
   tangential planes due to astigmatism/coma.

### Impact

- **MEDIUM**. The current monochromatic diffraction-limited MTF is a useful
  baseline but insufficient for real lens analysis.

### Recommendation

1. Add `mtf_polychromatic(wavelengths, weights, f_num, spatial_freq)`.
2. Add `mtf_through_focus(wavelength, f_num, spatial_freq, defocus_range)`.
3. Aberrated MTF can wait for Zernike polynomial support.

### References

- [Zemax: Methods for analyzing MTF](https://support.zemax.com/hc/en-us/articles/1500005575102)
- [Eckhardt Optics: Through Focus MTF](https://www.eckop.com/resources/optics/mtf-modulation-transfer-function/through-focus-mtf/)

---

## 6. Gaussian Beam & Higher-Order Modes

### Current state

Via bijli backend: Gaussian beam propagation with ABCD matrix re-exported.

### What's missing

1. **Hermite-Gaussian modes (HG_mn)** -- field amplitude:
   ```
   E_mn(x,y,z) = E_0 * H_m(sqrt(2)*x/w) * H_n(sqrt(2)*y/w)
                  * exp(-r^2/w^2) * exp(-jkz - jk*r^2/(2R) + j*(m+n+1)*zeta)
   ```
   where H_m are Hermite polynomials, w(z) is the beam width, R(z) is the
   radius of curvature, and zeta(z) is the Gouy phase.

2. **Laguerre-Gaussian modes (LG_pl)** -- field amplitude in cylindrical
   coordinates:
   ```
   E_pl(r,phi,z) = E_0 * (r*sqrt(2)/w)^|l| * L_p^|l|(2r^2/w^2)
                    * exp(-r^2/w^2) * exp(-jlphi) * exp(...)
   ```
   where L_p^l are generalized Laguerre polynomials. These carry OAM.

3. **Higher-order mode ABCD propagation** -- the ABCD law transforms q
   identically for all modes, but the Gouy phase gains a factor of
   `(m+n+1)` or `(2p+|l|+1)`.

4. **M^2 beam quality factor** -- relates a real beam to a fundamental
   Gaussian: `M^2 = pi * w0 * theta / lambda`. For HG_mn: `M_x^2 = 2m+1`,
   `M_y^2 = 2n+1`.

### Impact

- **MEDIUM**. The fundamental Gaussian is sufficient for many applications,
  but laser physics consumers (joshua?) need higher-order modes.

### Recommendation

1. Add `hermite_gaussian_field(m, n, x, y, z, w0, wavelength)`.
2. Add `laguerre_gaussian_field(p, l, r, phi, z, w0, wavelength)`.
3. Add `beam_quality_m2(divergence_half_angle, beam_waist, wavelength)`.

### References

- [Gaussian, Hermite-Gaussian, and Laguerre-Gaussian beams: A primer (arXiv)](https://arxiv.org/pdf/physics/0410021)
- [Gaussian beam -- Wikipedia](https://en.wikipedia.org/wiki/Gaussian_beam)

---

## 7. Fiber Optics

### Current state

Not implemented. `lens::numerical_aperture()` exists for free-space optics.

### What's missing

Fiber optics is a common consumer need. Key formulas:

1. **Fiber NA**: `NA = sqrt(n_core^2 - n_clad^2)`
2. **V-number (normalized frequency)**: `V = 2*pi*a*NA / lambda`
   where a is fiber core radius. Single-mode cutoff: `V < 2.405`.
3. **Number of modes**: `N_modes ~ V^2 / 2` (step-index multimode)
4. **Mode field diameter** (Gaussian approximation for single-mode):
   `MFD ~ 2*a*(0.65 + 1.619*V^{-3/2} + 2.879*V^{-6})`
   (Marcuse approximation)
5. **Coupling efficiency**: overlap integral of input and fiber modes:
   `eta = |integral(E_in * E_fiber* dA)|^2 / (integral(|E_in|^2) * integral(|E_fiber|^2))`
6. **Fiber dispersion**: material + waveguide dispersion per unit length.

### Impact

- **MEDIUM**. Not critical for current consumers (PBR, imaging), but
  important for a comprehensive optics library and the joshua simulation
  consumer.

### Recommendation

Add a `fiber` submodule with `fiber_na`, `v_number`, `num_modes`,
`mode_field_diameter`, `coupling_efficiency_gaussian`. This is
self-contained and high-value.

### References

- [Newport: Fiber Optic Basics](https://www.newport.com/t/fiber-optic-basics)
- [RP Photonics: Numerical Aperture](https://www.rp-photonics.com/numerical_aperture.html)

---

## 8. Photometry vs Radiometry

### Current state

Prakash has `Rgb::luminance()` using Rec. 709 coefficients, and spectral
radiance via `planck_radiance()`. But there is no explicit photometric
pathway.

### What's missing

The radiometry-to-photometry bridge requires the **CIE spectral luminous
efficiency function V(lambda)**:

```
Luminous flux (lm) = 683 * integral(Phi_e(lambda) * V(lambda) d_lambda)
```

where `Phi_e` is spectral radiant flux (W/nm) and 683 lm/W is the maximum
luminous efficacy at 555 nm.

Specific missing items:
1. **V(lambda) table** -- tabulated at 1nm or 5nm, 380-780nm
2. **Scotopic V'(lambda)** -- night vision, peak at 507nm
3. **Luminous flux** from SPD: `Phi_v = 683 * sum(SPD * V * delta_lambda)`
4. **Luminous intensity** (candela), **illuminance** (lux), **luminance**
   (cd/m^2) from spectral data
5. **Luminous efficacy** of a source: `eta = Phi_v / Phi_e`
6. **Color Rendering Index** (CRI) already exists -- good
7. **CCT from chromaticity** -- Robertson's method or Ohno's method for
   computing correlated color temperature from (u',v') coordinates

### Impact

- **MEDIUM-HIGH**. Lighting applications (kiran, aethersafta, tazama) need
  photometric quantities. The V(lambda) function is foundational and the
  data is public domain (CIE).

### Recommendation

1. Add `V_LAMBDA_PHOTOPIC: [(f64, f64); 81]` (wavelength, value) at 5nm.
2. Add `luminous_flux(spd)`, `luminous_efficacy(spd)`.
3. Add `scotopic_v_lambda` for completeness.
4. Consider adding `cct_from_xy(x, y)` using Robertson's or Ohno's method.

### References

- [CIE spectral luminous efficiency](https://cie.co.at/datatable/cie-spectral-luminous-efficiency-photopic-vision)
- [Luminous efficiency function -- Wikipedia](https://en.wikipedia.org/wiki/Luminous_efficiency_function)
- [Radiometry and Photometry -- Ossila](https://www.ossila.com/pages/radiometry-and-photometry)

---

## 9. CIE Observer Accuracy

### Current state

`CIE_1931_2DEG` at 5nm intervals, 380-780nm. Used for XYZ tristimulus and
SPD integration.

### What's missing

1. **CIE 2015 cone-fundamental-based CMFs** -- the 2015 CMFs are
   physiologically derived from actual cone sensitivities rather than
   1920s color-matching experiments. They provide improved matching accuracy
   (confirmed by 54-observer studies). Available for both 2-deg and 10-deg
   field sizes.

2. **CIE 1964 10-degree observer** -- already widely used for large-field
   colorimetry (e.g., uniform color patches in displays).

3. **1nm resolution tables** -- 5nm is adequate for most applications but
   1nm tables exist and eliminate interpolation error for precise spectral
   integration.

4. **Age-dependent observer** -- the CIE 2006 model allows computing CMFs
   for observers of different ages (relevant for vision science).

### Current accuracy of CIE 1931 2-deg values

The values in prakash match published data. However, the original 1931
CMFs are known to have a systematic error in the blue region (below 460nm)
due to measurement limitations in the original Wright/Guild experiments.
The 2015 CMFs correct this.

### Impact

- **MEDIUM**. For 95% of rendering/display applications, CIE 1931 is
  standard and expected. Adding CIE 2015 positions prakash as genuinely
  world-class and forward-looking.

### Recommendation

1. Add `CIE_2015_2DEG` and `CIE_2015_10DEG` tables.
2. Add `CIE_1964_10DEG` table (widely used for paints, displays).
3. Provide functions that accept a generic observer type.

### References

- [CIE 2015 Cone-Fundamental-Based Colorimetry (PDF)](https://files.cie.co.at/LpR%2060%20CIE%20RESEARCH%20Special%20-%20Cone-Fundamental-Based%20CIE%20Colorimetry.pdf)
- [CIE 1931 color space -- Wikipedia](https://en.wikipedia.org/wiki/CIE_1931_color_space)
- [colour-science Python: CMFs dataset](https://colour.readthedocs.io/en/v0.3.7/colour.colorimetry.dataset.cmfs.html)

---

## 10. Planck Function Numerical Stability

### Current state

```rust
pub fn planck_radiance(wavelength_m: f64, temperature_k: f64) -> f64 {
    let lambda5 = wavelength_m.powi(5);
    PLANCK_C1 / (lambda5 * ((PLANCK_C2 / (wavelength_m * temperature_k)).exp() - 1.0))
}
```

### Issues

1. **Overflow at low temperatures or short wavelengths**: when
   `hc/(lambda*kT)` is large (e.g., UV at room temperature), `exp(x)` for
   x > ~709 overflows to `f64::INFINITY`. The function then returns
   `PLANCK_C1 / (lambda5 * INFINITY)` = 0.0, which is physically correct
   (negligible radiation) but loses precision.

2. **Underflow at high temperatures or long wavelengths**: when
   `hc/(lambda*kT)` is very small, `exp(x) - 1 ~ x` and floating-point
   subtraction loses significance. Should use `expm1()` instead of
   `exp() - 1.0` for this regime.

3. **Edge case: T = 0** -- division by zero in the exponent.

### Fix

```rust
pub fn planck_radiance(wavelength_m: f64, temperature_k: f64) -> f64 {
    if temperature_k <= 0.0 { return 0.0; }
    let x = PLANCK_C2 / (wavelength_m * temperature_k);
    let lambda5 = wavelength_m.powi(5);
    if x > 500.0 {
        // Wien approximation regime -- exp(x) >> 1
        PLANCK_C1 / (lambda5 * x.exp())
    } else {
        PLANCK_C1 / (lambda5 * x.exp_m1())  // expm1 for numerical stability
    }
}
```

### Impact

- **HIGH for correctness**. Using `exp() - 1.0` instead of `expm1()` is
  a known numerical issue. The current code works for typical visible-range
  temperatures but silently loses precision in edge cases.

### Recommendation

Replace `exp() - 1.0` with `expm1()` (Rust's `f64::exp_m1()`). Add Wien
approximation branch for very large exponent arguments. Add `T <= 0` guard.

---

## 11. Rayleigh Scattering: King Correction Factor

### Current state

```rust
const RAYLEIGH_PREFACTOR: f64 = {
    let n2m1_sq = (N_AIR * N_AIR - 1.0) * (N_AIR * N_AIR - 1.0);
    let pi3 = PI * PI * PI;
    (8.0 * pi3 / 3.0) * n2m1_sq / (N_S * N_S)
};
```

This is the simple `lambda^{-4}` formula without the King correction factor.

### What's missing

The full Rayleigh scattering cross-section is:

```
sigma(lambda) = (24*pi^3 / (lambda^4 * N_s^2)) * ((n^2-1)/(n^2+2))^2 * F_K(lambda)
```

where `F_K(lambda) = (6 + 3*rho_n) / (6 - 7*rho_n)` is the **King
correction factor** and `rho_n` is the depolarization ratio.

For air:
- `F_K` ~ 1.048 (about 4.8% correction)
- `rho_n` varies with wavelength (0.0279 at 500nm to 0.0301 at 400nm)

Research (Bodhaine et al. 1999, Bucholtz 1995) shows that using a constant
depolarization factor introduces errors of **1.3% to 10%** compared to
wavelength-dependent values, with the largest errors in the UV.

### Impact

- **MEDIUM for atmospheric accuracy**. The 4.8% correction is systematic
  and matters for quantitative atmospheric modeling. For rendering (sky
  color), it's less critical.

### Recommendation

1. Add wavelength-dependent King correction factor `king_factor(wavelength_m)`.
2. Provide corrected `rayleigh_cross_section_corrected(wavelength_m)`.
3. Keep the simple version as the default for performance-critical paths.

### References

- [Bodhaine et al. 1999: Improved Rayleigh-scattering optical depth](https://opg.optica.org/ao/abstract.cfm?uri=ao-44-16-3320)
- [Bucholtz 1995: Rayleigh-scattering calculations](https://opg.optica.org/ao/abstract.cfm?uri=ao-34-15-2765)
- [Rayleigh scattering -- Wikipedia](https://en.wikipedia.org/wiki/Rayleigh_scattering)

---

## 12. Mie Scattering: Cornette-Shanks Limitations

### Current state

`mie_phase_cornette_shanks()` with g=0.76 default. The roadmap notes "Full
Mie theory" as a P2 item. bijli provides `bijli::scattering::mie()` when
the bijli-backend feature is enabled.

### When Cornette-Shanks breaks down

Research (NVIDIA 2023, Wiscombe) shows:

1. **Sharp forward/backward peaks** -- CS cannot reproduce the narrow
   forward lobe of large particles (size parameter x = 2*pi*r/lambda > 10).
   Real Mie has glory peaks and rainbow features that CS completely misses.

2. **Fog, clouds, tissue** -- for water droplets (x ~ 100-1000), CS is
   qualitatively wrong. The 2023 NVIDIA paper shows that a blend of
   Henyey-Greenstein + Draine phase function matches 95% of exact Mie.

3. **Wavelength dependence** -- CS uses a single g parameter that doesn't
   capture how the phase function shape changes with wavelength.

### Impact

- **MEDIUM-HIGH for simulation accuracy**. For sky rendering, CS is
  adequate. For fog, clouds, biomedical optics, or particle characterization,
  it's insufficient.

### Recommendation

1. Promote the bijli full Mie solver to first-class status with ergonomic
   wrappers.
2. Add a `DrainePhaseFunction` as a better analytical approximation for
   intermediate cases.
3. Consider the NVIDIA approximate Mie function for real-time fog/cloud.

### References

- [NVIDIA: Approximate Mie Scattering Function for Fog and Cloud](https://research.nvidia.com/labs/rtr/approximate-mie/)
- [Cornette & Shanks phase function (ResearchGate)](https://www.researchgate.net/publication/45721602)
- [Mie scattering -- Wikipedia](https://en.wikipedia.org/wiki/Mie_scattering)

---

## 13. Polarization Ray Tracing

### Current state

Prakash has Stokes vectors, Mueller matrices, Jones calculus for discrete
elements, birefringent materials. The roadmap lists "polarization ray
tracing" as P2.

### State of the art

Modern polarization ray tracing (PRT) uses a **3x3 polarization ray-tracing
matrix** (Chipman, 2010) per surface, which generalizes Jones calculus to
3D. At each surface:

1. Compute the local s/p coordinate frame from ray direction and surface
   normal.
2. Compute Fresnel coefficients (complex, for s and p separately).
3. Construct the 3x3 PRT matrix for that surface.
4. Transform the polarization state through the surface.

The product of all surface PRT matrices gives the total system Jones matrix
for that ray.

Zemax 25R2 (2025) uses Rodrigues' rotation formula to handle off-axis
incidence correctly in their Jones matrix surface implementation.

### Impact

- **HIGH for a world-class library**. This is the bridge between prakash's
  existing polarization calculus and its ray tracer.

### Recommendation

1. Add `PolarizationRayTraceMatrix` -- 3x3 complex matrix per surface.
2. Extend `TraceRay` to optionally carry polarization state.
3. At each surface in `trace_sequential`, compute local s/p frame and
   construct the PRT matrix from the complex Fresnel coefficients.

### References

- [Chipman 2010: Three-dimensional polarization ray-tracing calculus](https://opg.optica.org/ao/abstract.cfm?uri=ao-50-18-2855)
- [Zemax: Jones Matrix surface](https://optics.ansys.com/hc/en-us/articles/43071140222099)

---

## 14. Vectorial Diffraction (Richards-Wolf)

### Current state

Prakash has scalar Fraunhofer and Fresnel diffraction, Huygens-Fresnel
integral (1D). All scalar (no polarization).

### What's missing

The **Richards-Wolf vectorial diffraction integral** (1959) computes the
full 3D electromagnetic field near the focus of a high-NA lens:

```
E(r_p, phi_p, z_p) = -jf/(2*lambda) * integral_0^alpha integral_0^{2pi}
    P(theta, phi) * a(theta) * exp(jk * r_p * sin(theta) * cos(phi - phi_p))
    * exp(jk * z_p * cos(theta)) * sin(theta) * d_phi * d_theta
```

where alpha = arcsin(NA/n) and P is the pupil polarization distribution.

This is essential for:
- Tight focusing (NA > 0.7) where scalar diffraction fails
- Focusing of radially/azimuthally polarized beams
- Near-field optical microscopy simulation

### Impact

- **LOW-MEDIUM for current consumers**. This is advanced research-grade
  physics. Not needed for PBR or standard imaging, but positions prakash
  at the cutting edge.

### Recommendation

This is a P3 item. If implemented, use the FFT-based fast algorithm rather
than direct double integration (4-7 orders of magnitude faster).

### References

- [Richards & Wolf 1959 (ADS)](https://ui.adsabs.harvard.edu/abs/1959RSPSA.253..358R)
- [Simple program for tightly focused fields (arXiv)](https://arxiv.org/pdf/2211.06725)

---

## 15. Zernike Polynomials

### Current state

Listed as P2 roadmap item. Not implemented.

### Why this matters

Zernike polynomials are the standard basis for describing wavefront
aberrations over a circular pupil:

```
Z_n^m(rho, theta) = N_n^m * R_n^|m|(rho) * { cos(m*theta) for m >= 0
                                              { sin(|m|*theta) for m < 0
```

where R_n^m is the radial polynomial and N_n^m is the normalization.

**Every major optical design tool** (Zemax, Code V, OSLO, FRED) uses
Zernike decomposition for:
- Wavefront error specification
- PSF computation from wavefront
- Tolerancing
- Adaptive optics correction

The standard ordering schemes are: Noll (1976), ANSI/OSA, Fringe (Zygo).

### Impact

- **HIGH**. This is arguably the most important missing feature for
  serious optical system analysis. It bridges aberration theory, wavefront
  sensing, and PSF/MTF computation.

### Recommendation

1. Add `zernike(n, m, rho, theta) -> f64` using standard Noll ordering.
2. Add `zernike_radial(n, m, rho) -> f64`.
3. Add `ZernikeCoefficients` struct with named fields for common terms
   (defocus, astigmatism, coma, spherical, trefoil, etc.).
4. Add `wavefront_from_zernike(coeffs, rho, theta) -> f64`.
5. Add `psf_from_zernike(coeffs, wavelength, pupil_radius) -> 2D array`.

### References

- [Zernike polynomials -- Wikipedia](https://en.wikipedia.org/wiki/Zernike_polynomials)
- [Wyant: Zernike Polynomials (Arizona)](https://wp.optics.arizona.edu/jcwyant/wp-content/uploads/sites/13/2016/08/Zernike_Polynomials_For_The_Web.pdf)
- [2025: Conversion of Zernike to power series (ScienceDirect)](https://www.sciencedirect.com/science/article/pii/S2211379725003249)

---

## 16. GRIN Optics

### Current state

Listed as P2 roadmap item. Not implemented.

### State of the art

GRIN ray tracing requires solving the eikonal equation numerically:

```
d/ds (n(r) * dr/ds) = grad(n)
```

This is an ODE system solved by:
1. **Runge-Kutta methods** (RK4) -- most common, good accuracy
2. **Optical path length discretization** -- shown to be consistently
   superior in accuracy (Optica 2026)
3. **Medium discretization** -- divide into homogeneous layers, apply Snell
   at each interface (simplest but least accurate)

Common GRIN profiles: axial (`n(z)`), radial (`n(r) = n0 - n2*r^2`),
spherical.

### Impact

- **LOW-MEDIUM**. GRIN lenses are used in endoscopes, gradient-index fiber,
  and some camera lenses, but it's a niche application.

### Recommendation

Implement basic radial GRIN ray tracing using RK4. This is conceptually
clean and self-contained.

### References

- [Ray trace algorithms for GRIN media](https://opg.optica.org/ao/abstract.cfm?uri=ao-26-15-2943)
- [Comparative analysis of GRIN ray-tracing methods](https://opg.optica.org/ao/abstract.cfm?uri=ao-65-1-10)

---

## 17. Diffractive Optical Elements

### Current state

Listed as P2 roadmap item. Not implemented.

### Key formulas

1. **Grating equation**: `n2*sin(theta_m) - n1*sin(theta_i) = m*lambda/d`
   where d is grating period, m is diffraction order.

2. **Diffractive lens phase function**: `phi(r) = -pi*r^2 / (lambda*f)`
   (quadratic phase, equivalent to thin lens).

3. **Diffraction efficiency** (scalar): for a blazed grating with blaze
   angle optimized for order m at design wavelength:
   `eta_m(lambda) = sinc^2(m - lambda_design/lambda)`

4. **Kinoform efficiency**: for a multi-level DOE with N phase levels:
   `eta = sinc^2(1/N)`

### Impact

- **LOW-MEDIUM**. DOEs are used in AR/VR optics, laser beam shaping, and
  spectrometers. The grating equation alone is high-value and trivial to add.

### Recommendation

1. Add `grating_diffraction_angle(n1, n2, theta_i, order, wavelength, period)`.
2. Add `grating_blaze_efficiency(order, wavelength, design_wavelength)`.
3. Defer full DOE simulation (RCWA) to a future release.

---

## 18. Fluorescence

### Current state

Listed as P2 roadmap item. Not implemented.

### Key physics

1. **Stokes shift**: `delta_nu = nu_absorption - nu_emission` (always > 0)
2. **Quantum yield**: `QY = photons_emitted / photons_absorbed`
3. **Fluorescence intensity**: `I_f = I_0 * epsilon * c * l * QY * phi`
   where epsilon = molar absorptivity, c = concentration, l = path length
4. **Excitation/emission spectra**: typically represented as normalized
   spectral distributions (similar to SPD)
5. **Kasha's rule**: emission spectrum is independent of excitation wavelength
6. **Forster resonance energy transfer (FRET)**: `E = 1/(1 + (r/R0)^6)`

### Impact

- **LOW**. Fluorescence is specialized. It matters for biological imaging,
  forensics simulation, and certain display technologies, but is not core
  to prakash's primary consumers.

### Recommendation

Defer unless joshua specifically needs it. If implemented, start with
Stokes shift representation and quantum yield, which are simple additions
to the spectral module.

---

## 19. Nonlinear Optics

### Current state

Listed as P2 roadmap item. Not implemented.

### Key formulas

1. **Second-harmonic generation** (chi-2):
   Coupled wave equations in the slowly-varying envelope approximation:
   ```
   dA_2w/dz = j * kappa * A_w^2 * exp(j*Delta_k*z)
   ```
   where `kappa = omega * d_eff / (n_2w * c)` and `Delta_k = k_2w - 2*k_w`
   is the phase mismatch.

2. **Kerr effect** (chi-3): intensity-dependent refractive index:
   ```
   n(I) = n_0 + n_2 * I
   ```
   where n_2 is the nonlinear refractive index coefficient.

3. **Self-phase modulation**: `phi_NL = n_2 * k * I * L`

4. **Phase matching**: `Delta_k = 0` condition for efficient conversion.

### Impact

- **LOW**. Nonlinear optics is a specialized field. Unless joshua needs
  laser simulation with frequency doubling, this is low priority.

### Recommendation

Defer. If needed, the Kerr effect (`n_effective(n0, n2, intensity)`) is a
trivial one-liner to add.

---

## 20. Orbital Angular Momentum

### Current state

Not in roadmap.

### Key physics

Light beams can carry OAM of l*hbar per photon, distinct from spin angular
momentum (polarization). Implemented via Laguerre-Gaussian modes (Section 6).

### Impact

- **LOW**. This is primarily a quantum optics / telecommunications topic.
  If LG modes are implemented, OAM comes free.

---

## 21. Competitive Analysis

### Rust ecosystem

There is no comparable Rust optics library. The closest are:
- `rustgeotrace` -- basic geometric ray tracing only
- `optics` crate -- functional programming optics (lenses/prisms), not physics
- Various ray tracers -- rendering focused, not optical engineering

Prakash is **unique in the Rust ecosystem** as a physics-complete optics library.

### Python: colour-science

The `colour` Python library is the gold standard for color science. It has:
- CIE 1931, 1964, 2006, 2015 observers
- 200+ illuminants
- Full photometric/radiometric conversions
- CCT computation (Robertson, Ohno)
- Chromatic adaptation transforms (Bradford, Von Kries, CAT02, CAT16)

Prakash is missing: CIE 2015 observers, photometric functions, CCT from
chromaticity, chromatic adaptation transforms.

### C++: Mitsuba3, PBRT

- Full complex Fresnel equations
- Spectral rendering with arbitrary basis functions
- Polarization ray tracing (Mitsuba3)

### Commercial: Zemax, Code V, FRED

Core math that these use which prakash lacks:
- Zernike polynomials (essential)
- Polarization ray tracing (3x3 PRT matrices)
- 11+ dispersion formulas
- GRIN ray tracing
- Full Mie theory
- Complex Fresnel equations
- Polychromatic/through-focus MTF
- Fiber coupling analysis

---

## 22. P2 Roadmap Item Assessment

| P2 Item | Worth it? | Priority | Notes |
|---------|-----------|----------|-------|
| Polarization ray tracing | YES | **P1** | Essential for world-class. Builds on existing Mueller/Jones. |
| GRIN optics | Maybe | P2 | Niche but clean to implement. RK4 integration. |
| Zernike polynomials | YES | **P1** | Most important single missing feature. |
| Diffractive optical elements | YES | P2 | Grating equation is trivial. Full DOE is complex. |
| Fluorescence | No (for now) | P3 | Wait for consumer demand. |
| Nonlinear optics | No (for now) | P3 | Wait for joshua needs. |
| GPU-friendly f32 API | Maybe | P2 | Depends on aethersafta/kiran needs. |
| Full Mie theory | YES | P2 | bijli already has it. Need ergonomic wrappers. |
| CIE 2015 10-deg observer | YES | **P1** | Trivial to add, high completeness value. |

---

## 23. Prioritized Recommendations

### Tier 1: Critical for 100% Accuracy (do first)

1. **Complex Fresnel equations** -- `fresnel_s_complex`, `fresnel_p_complex`
   with complex refractive index. This unlocks metals, semiconductors, and
   correct thin-film calculations.

2. **Planck function fix** -- replace `exp() - 1.0` with `expm1()`, add
   Wien approximation branch, add T=0 guard.

3. **Zernike polynomials** -- `zernike(n, m, rho, theta)`, standard Noll
   ordering, wavefront reconstruction, PSF from Zernike coefficients.

4. **Transfer matrix method upgrade** -- oblique incidence, complex n,
   separate s/p results, transmittance output.

### Tier 2: Important for Completeness (do next)

5. **Photometric functions** -- V(lambda) table, luminous flux, luminous
   efficacy, scotopic V'(lambda).

6. **CIE 2015 + 1964 observers** -- add the tables, parameterize functions
   on observer choice.

7. **Chromatic aberration expansion** -- lateral CA, secondary spectrum,
   partial dispersion.

8. **Polychromatic MTF** -- weighted sum of monochromatic MTFs.

9. **Polarization ray tracing** -- 3x3 PRT matrix per surface in the
   sequential tracer.

10. **King correction factor** -- wavelength-dependent Rayleigh scattering
    correction.

### Tier 3: Nice to Have (world-class polish)

11. **Herzberger/Schott/Conrady dispersion** -- legacy compatibility.

12. **Fiber optics module** -- NA, V-number, modes, coupling.

13. **Hermite-Gaussian / Laguerre-Gaussian modes** -- higher-order beam
    profiles.

14. **Full Mie wrappers** -- ergonomic API around bijli's solver.

15. **Through-focus MTF** -- MTF vs defocus.

16. **Grating equation + blaze efficiency** -- basic DOE support.

### Tier 4: Research-Grade (future)

17. Richards-Wolf vectorial diffraction
18. GRIN ray tracing
19. Nonlinear optics primitives
20. Fluorescence / Stokes shift
21. Orbital angular momentum (comes with LG modes)

---

## Summary

Prakash has an impressively solid foundation. The core ray optics, wave
optics, spectral, and PBR modules are well-implemented with correct physics.
The main gaps fall into three categories:

1. **Accuracy gaps**: complex Fresnel equations, Planck expm1, King factor.
   These are cases where the current code gives approximately right answers
   but not 100% correct ones.

2. **Missing standard tools**: Zernike polynomials, photometric functions,
   polychromatic MTF, CIE 2015 observers. These are expected by any serious
   optics practitioner.

3. **Advanced physics**: polarization ray tracing, vectorial diffraction,
   GRIN, nonlinear optics. These separate a good library from a world-class
   one.

The Tier 1 items (complex Fresnel, Planck fix, Zernike, TMM upgrade) would
close the most critical gaps. The Tier 2 items (photometry, CIE 2015,
chromatic aberrations, polarization ray tracing) would make prakash
genuinely comprehensive. Together, these 10 items would position prakash
as the most complete open-source optics library in any language.
