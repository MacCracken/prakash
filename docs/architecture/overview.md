# Architecture

## Module Map

| Module | Feature | Files | Lines | Key Types | Purpose |
|--------|---------|-------|-------|-----------|---------|
| `error` | always | `error.rs` | 174 | `PrakashError`, `Result<T>` | Shared error types with `Cow<'static, str>` messages, bijli error bridge |
| `ray` | `ray` | `mod.rs`, `dispersion.rs`, `simulate.rs`, `system.rs`, `trace.rs` | 4144 | `Medium`, `ComplexMedium`, `CauchyCoefficients`, `SellmeierCoefficients`, `TraceRay`, `OpticalSurface`, `PolarizedTraceHit`, `ParaxialRay`, `Prescription` | Geometric optics: Snell, Fresnel (real + complex), dispersion, sequential/recursive tracing with polarization, ray fans, spot diagrams, OPD |
| `wave` | `wave` | `mod.rs`, `diffraction.rs`, `pattern.rs`, `polarization.rs`, `zernike.rs` | 4007 | `Polarization`, `StokesVector`, `MuellerMatrix`, `Pattern2D`, `ZernikeWavefront` | Wave optics: interference, diffraction, Fabry-Perot, AR coatings, Jones/Stokes/Mueller, Zernike polynomials, 2D patterns, PSF |
| `spectral` | `spectral` | `mod.rs`, `cie.rs` | 1496 | `Rgb`, `Xyz`, `Spd` | Color science: wavelength-to-RGB, Planck (numerically stable), Wien, CIE 1931 XYZ, SPD (Cow storage), illuminants, CRI |
| `lens` | `lens` | `lens.rs` | 1102 | `LensType`, `CardinalPoints`, `SeidelCoefficients` | Lens/mirror geometry: thin/thick lens, aberrations, MTF, DoF, Petzval, multi-element |
| `pbr` | `pbr` | `mod.rs`, `advanced.rs` | 1628 | (free functions) | PBR shading: Cook-Torrance, GGX, sheen, clearcoat, SSS, iridescence, volumetric, importance sampling, env maps |
| `atmosphere` | `atmosphere` | `atmosphere.rs` | 763 | (free functions + constants) | Atmospheric optics: Rayleigh/Mie scattering, sky color, air mass, optical depth, sunset model |
| `ai` | `ai` | `ai.rs` | 127 | `DaimonClient`, `DaimonConfig` | AI-assisted optics queries via reqwest/tokio |
| `logging` | `logging` | `logging.rs` | 34 | (init functions) | tracing-subscriber setup |

**Total**: ~13,540 lines across 20 source files.

## Design Principles

- Flat library crate -- no internal binaries
- Feature-gated modules -- consumers pull only what they need
- f64 precision throughout -- optics demands double precision
- No external FFT/linalg deps -- inline implementations for zero dependency cost
- `#[must_use]` on all pure functions, `#[non_exhaustive]` on all public enums
- `#![warn(missing_docs)]` -- compile-time doc coverage enforcement
- `Cow<'static, str>` in error variants -- zero alloc on static messages
- `Cow<'static, [f64]>` in SPD values -- zero alloc for standard illuminants
- Precomputed constants where possible (Rayleigh prefactor, 1/pi, CIE tables)
- `#[inline]` on hot-path functions -- every optics function that does math
- tracing instrumentation on functions that involve multi-step computation

## Data Flow

```
spectral ──> color science (wavelength <-> RGB, blackbody, CIE XYZ, SPD, CRI)
    |
    v
ray ──> geometric optics (Snell, Fresnel [real + complex], dispersion, sequential trace)
    |         |
    |         +──> ray::trace (polarization-aware sequential tracing, s/p tracking)
    |         +──> ray::system (paraxial trace, prescriptions, cardinal points)
    |         +──> ray::simulate (recursive trace, ray fans, spot diagrams, OPD)
    |
wave ──> wave optics (interference, diffraction, Fabry-Perot, AR coatings)
    |         |
    |         +──> wave::polarization (Stokes/Mueller formalism)
    |         +──> wave::zernike (Zernike polynomials, wavefront decomposition, Strehl ratio)
    |         +──> wave::pattern (2D diffraction patterns, PSF, spectrum strips)
    |
lens ──> lens geometry (thin/thick lens, aberrations, MTF, DoF, Petzval)
    |
pbr ──> PBR shading (Cook-Torrance, GGX, Fresnel-Schlick)
    |         |
    |         +──> pbr::advanced (sheen, clearcoat, SSS, iridescence, volumetric, importance sampling)
    |
atmosphere ──> sky models (Rayleigh/Mie, air mass, optical depth, sunset)
    |
    v
error ──> shared by all modules (PrakashError, Result<T>)
```

All physics modules import `crate::error`. No module imports another physics module. The flow is always: caller picks features, each module operates independently, error types unify them.

## Bijli Backend

When the `bijli-backend` feature is enabled (default), prakash delegates foundational EM physics to the bijli crate:

| Delegated | prakash function | bijli function |
|-----------|-----------------|----------------|
| Snell's law | `ray::snell` | `bijli::wave::snell_refraction_angle` |
| Critical angle | `ray::critical_angle` | `bijli::wave::critical_angle` |
| Brewster angle | `ray::brewster_angle` | `bijli::wave::brewster_angle` |
| Normal reflectance | `ray::fresnel_normal` | `bijli::wave::reflectance_normal` |
| Unpolarized reflectance | `ray::fresnel_unpolarized` | `bijli::wave::reflectance_unpolarized` |
| Speed of light | `spectral::SPEED_OF_LIGHT` | `bijli::field::SPEED_OF_LIGHT` |

New types re-exported from bijli: `JonesVector`, `JonesMatrix`, `Complex`, `GaussianBeam`, `AbcdMatrix`, `ResonatorStability`, `MieResult`.

Bidirectional `From` conversions: `StokesVector`, `MuellerMatrix`, `Polarization` → `JonesVector`.

All delegates have `#[cfg(not(feature = "bijli-backend"))]` fallback implementations, so prakash compiles and works without bijli.

## Feature Independence

Each module compiles independently. No cross-feature imports except `error.rs`.

Intentional duplications to avoid cross-feature dependencies:

- `wave::pattern` contains its own `wavelength_to_rgb` -- avoids requiring `spectral` feature for pattern visualization
- `atmosphere` defines its own `RGB_WAVELENGTHS` constant -- avoids requiring `spectral` feature for sky color computation

This means a consumer can use `--features wave` without pulling in `spectral`, or `--features atmosphere` without pulling in `spectral` or `ray`.

Default features: `ray`, `wave`, `spectral`, `lens`, `bijli-backend`. The `pbr`, `atmosphere`, `ai`, and `logging` features are opt-in.

## Dependencies

### Required

| Crate | Purpose |
|-------|---------|
| `hisab` | Math foundations (AGNOS crate) |
| `serde` | Serialization for all public types |
| `thiserror` | Derive `Error` for `PrakashError` |
| `tracing` | Structured instrumentation |

### Optional (feature-gated)

| Crate | Feature | Purpose |
|-------|---------|---------|
| `bijli` | `bijli-backend` | EM-correct Fresnel/Snell, polarization, Gaussian beams, Mie scattering (AGNOS crate) |
| `reqwest` | `ai` | HTTP client for AI queries |
| `tokio` | `ai` | Async runtime |
| `serde_json` | `ai` | JSON serialization for AI payloads |
| `tracing-subscriber` | `logging` | Log output formatting and filtering |

### Dev-only

| Crate | Purpose |
|-------|---------|
| `criterion` | Benchmarking |
| `serde_json` | Test serialization round-trips |
