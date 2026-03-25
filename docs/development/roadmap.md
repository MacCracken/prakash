# Prakash Roadmap

> **Prakash** is the optics/light simulation crate. Math foundations come from [hisab](https://github.com/MacCracken/hisab). Rendering integration is via [aethersafta](https://github.com/MacCracken/aethersafta) and [kiran](https://github.com/MacCracken/kiran).

## Scope

Prakash owns the **physics of light**: how light travels, bends, reflects, interferes, diffracts, polarizes, and interacts with materials. It provides the math; consumers decide what to do with it (render pixels, simulate experiments, process images).

Prakash does NOT own:
- **Rendering pipeline** → aethersafta/kiran (they consume prakash for lighting math)
- **Image processing** → ranga (pixel operations, color spaces, filters)
- **Math primitives** → hisab (vectors, geometry, calculus)
- **Color science beyond spectral** → ranga (ICC profiles, gamut mapping)

## Completed

| Phase | Release | Summary |
|-------|---------|---------|
| V0.1 | Foundation | Ray optics (Snell, Fresnel, Beer-Lambert), wave basics (interference, diffraction, polarization), spectral (Planck, Wien, CIE), lens (thin lens, DoF), PBR (Cook-Torrance, GGX, Lambert) |
| V0.2 | Optical Systems | 3D refraction, sequential ray trace, Cauchy/Sellmeier dispersion, thick lens, cardinal points, Seidel aberrations, MTF, CIE 1931, SPD, illuminants, CRI |
| V0.3 | Wave Optics Expansion | Coherence, Airy/Bessel, Fabry-Pérot, Stokes/Mueller, birefringence, Fraunhofer/Fresnel diffraction, Huygens-Fresnel, AR coatings |
| V0.4 | Advanced PBR | Anisotropic GGX, sheen, clearcoat, SSS, iridescence, volumetric scattering, importance sampling, environment maps |
| V0.25 | Atmospheric Optics | Rayleigh/Mie scattering, sky color, air mass (Kasten & Young), optical depth, sunset gradient |
| V0.26 | Simulation Primitives | Recursive ray tracer, ray fans (meridional/sagittal/bundle), spot diagrams + RMS, OPD |
| V0.27 | Optical Bench | Paraxial y-nu trace, Prescription builder, system cardinal point finder, common presets (biconvex, planoconvex, doublet) |
| V0.28 | Pattern Computation | 2D FFT diffraction, N-source interference, spectrum strips, PSF from wavefront |
| V0.29 | Examples & Documentation | Rainbow simulation, camera lens simulator, PBR material preview, module docs |
| V1.0 | Stable Release | API review (naming, params), `#[must_use]`/`#[non_exhaustive]` audit, feature gate audit, doc coverage, README, architecture docs |

## Engineering Backlog

### P1 — Post-V1.0 Audit (complete)

- [x] Consolidate duplicate `Polarization` (Jones) and `StokesVector` — added `From<Polarization> for StokesVector` conversion
- [x] Consolidate `diffraction_limit` (lens) and `rayleigh_criterion` (wave) — cross-referenced docs, both kept for domain clarity
- [x] Consolidate `beer_lambert` (ray) and `volume_transmittance` (pbr) — cross-referenced docs showing equivalence
- [x] Unit-suffix audit: documented convention (`_nm`/`_m`/`_um`/unsuffixed) in lib.rs module docs
- [x] `Spd` type: changed to `Cow<'static, [f64]>` — illuminant functions now zero-alloc with `from_static()`
- [x] `#![warn(missing_docs)]` enabled in lib.rs — all 55 gaps fixed
- [x] Wavefront parameter ordering: wavelength-first convention across all wave/diffraction functions

### P1 — Consumer Integration Gaps

- [x] bijli EM backend: scalar Snell/Fresnel/Brewster delegate to bijli, polarization bridge (StokesVector/MuellerMatrix/JonesVector), Gaussian beam + ABCD re-exports, Mie scattering re-export, Medium::permittivity()
- [ ] hisab geometry bridge: adapter types between prakash `TraceRay`/`[f64;3]` and hisab `Ray`/`Vec3` (f32↔f64)
- [ ] kiran: verify Cook-Torrance `h_dot_v` parameter migration (breaking change from pre-V1.0)
- [ ] ranga: chromatic aberration filter needs `SellmeierCoefficients` — verify import path works
- [ ] joshua: atmospheric scattering API review with simulation team — does `sunset_gradient` cover their use cases?

### P1 — Accuracy & Completeness (promoted from research audit)

- [x] **Planck `expm1()` fix**: `exp_m1()` for numerical stability, Wien approximation for x>500, T≤0 guard
- [x] **Complex Fresnel equations**: `ComplexMedium`, `fresnel_s_complex`/`fresnel_p_complex`/`fresnel_unpolarized_complex`/`fresnel_normal_complex` with presets for gold, silver, copper, aluminum
- [x] **Zernike polynomials**: `wave::zernike` module — radial polynomials, full Z_n^m, Noll index conversion, `ZernikeWavefront` (evaluate, to_grid, RMS, P-V, Strehl), named aberrations (defocus, spherical, coma, astigmatism)
- [x] **Polarization ray tracing**: `trace_sequential_polarized` — tracks s/p transmittance per surface, cumulative through system, `PolarizedTraceHit` output

### P2 — Important Gaps (from research audit)

- [x] **Photometric functions**: V(λ) photopic + scotopic tables, luminous flux, luminous efficacy, scotopic flux — `spectral::photometry` module
- [ ] **CIE 2015 cone-fundamental observers**: 2-deg and 10-deg, plus CIE 1964 10-deg for legacy (data table entry needed)
- [x] **Transfer matrix method upgrade**: `multilayer_rt()` with oblique incidence, separate s/p reflectance + transmittance, `ThinFilmResult`
- [x] **Chromatic aberration expansion**: `lateral_chromatic_aberration()`, `secondary_spectrum()`, `partial_dispersion()`
- [x] **Polychromatic MTF**: `mtf_polychromatic()` with spectral weights, `mtf_through_focus()` with defocus range
- [x] **King correction factor**: `king_factor(wavelength_m)`, `rayleigh_cross_section_corrected()` (~4.8% correction)
- [x] **Fiber optics**: `ray::fiber` module — `fiber_na`, `v_number`, `num_modes`, `is_single_mode`, `mode_field_diameter`, `coupling_efficiency_gaussian`
- [x] **Additional dispersion models**: `HerzbergerCoefficients`, `SchottCoefficients` (BK7 preset), `ConradyCoefficients`

### P2 — Research & Future Features

- [ ] Gradient-index (GRIN) optics: curved ray paths through variable-n media
- [ ] Diffractive optical elements (DOE): phase gratings, holographic elements
- [ ] Vectorial diffraction (Richards-Wolf): high-NA focusing beyond scalar theory
- [ ] Hermite-Gaussian and Laguerre-Gaussian beam modes (HG_mn, LG_pl)
- [ ] Beam quality factor M² and propagation ratio
- [ ] GPU-friendly API: `f32` variants of hot-path functions for shader-side computation
- [ ] Higher-order (5th-order Buchdahl) monochromatic aberrations
- [ ] Wavefront aberration coefficients (W_040, W_131, etc.) from Seidel sums
- [ ] Aberrated MTF from generalized pupil function autocorrelation

### P3 — Advanced / Demand-Gated

- [ ] Fluorescence: Stokes shift, excitation/emission spectra
- [ ] Non-linear optics: second-harmonic generation, Kerr effect (if joshua needs it)
- [ ] Orbital angular momentum (OAM) of light — Laguerre-Gaussian modes carry OAM
- [ ] Metamaterials and negative refractive index
- [ ] Age-dependent CIE observer model (CIE 2006)

### P3 — Infrastructure

- [x] `cargo semver-checks` in CI (already in ci.yml)
- [x] Property-based testing (proptest) — 8 property tests for Snell, Fresnel, Planck, Beer-Lambert, Zernike, wavelength roundtrip
- [ ] Benchmark regression CI gate — fail on >10% regression
- [x] Doc tests for key functions — 8 doc tests across all modules (ray, wave, spectral, lens, pbr, atmosphere, zernike)
- [x] Coverage gate in CI — 85% threshold warning in ci.yml

## Consumers

| Consumer | What it uses |
|----------|-------------|
| **kiran** | PBR shading (Cook-Torrance, Fresnel-Schlick), lighting math |
| **aethersafta** | Color temperature for display calibration, compositing light math |
| **rasa** | Lens blur (DoF simulation), chromatic aberration filter, light leak effects |
| **joshua** | Simulation: laser/mirror puzzles, optics experiments, atmospheric effects |
| **tazama** | Color grading: color temperature adjustment, spectral analysis of footage |

## Boundary with Other Crates

| Feature | prakash | other |
|---------|---------|-------|
| Fresnel reflectance math | Yes (delegates to bijli) | bijli (EM foundation) |
| Pixel-level image filter | — | ranga |
| wgpu shader code | — | aethersafta |
| 3D scene graph | — | kiran |
| Vector/matrix math | — | hisab |
| Ray-geometry intersection | — | hisab (geo module) |
| Color space conversion (ICC) | — | ranga |
| Spectral → RGB conversion | Yes | — |
| Polarization formalism | Yes (bridge to bijli) | bijli (Jones/Mueller) |
| Gaussian beam / ABCD | Re-export from bijli | bijli (beam module) |
| Mie scattering (exact) | Re-export from bijli | bijli (scattering module) |
