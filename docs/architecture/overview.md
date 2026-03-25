# Architecture

## Module Map

| Module | Feature | Files | Lines | Key Types | Purpose |
|--------|---------|-------|-------|-----------|---------|
| `error` | always | `error.rs` | 126 | `PrakashError`, `Result<T>` | Shared error types with `Cow<'static, str>` messages |
| `ray` | `ray` | `mod.rs`, `dispersion.rs`, `simulate.rs`, `system.rs`, `trace.rs` | 3740 | `Medium`, `CauchyCoefficients`, `SellmeierCoefficients`, `TraceRay`, `OpticalSurface`, `ParaxialRay`, `Prescription` | Geometric optics: Snell, Fresnel, dispersion, sequential/recursive tracing, ray fans, spot diagrams, OPD |
| `wave` | `wave` | `mod.rs`, `diffraction.rs`, `pattern.rs`, `polarization.rs` | 3291 | `Polarization`, `StokesVector`, `MuellerMatrix`, `Pattern2D`, `PointSource`, `Spd` | Wave optics: interference, diffraction, Fabry-Perot, AR coatings, Jones/Stokes/Mueller, 2D patterns, PSF |
| `spectral` | `spectral` | `mod.rs`, `cie.rs` | 1338 | `Rgb`, `Xyz`, `Spd` | Color science: wavelength-to-RGB, Planck, Wien, CIE 1931 XYZ, SPD, illuminants, CRI |
| `lens` | `lens` | `lens.rs` | 1099 | `LensType`, `CardinalPoints`, `SeidelCoefficients` | Lens/mirror geometry: thin/thick lens, aberrations, MTF, DoF, Petzval, multi-element |
| `pbr` | `pbr` | `mod.rs`, `advanced.rs` | 1625 | (free functions) | PBR shading: Cook-Torrance, GGX, sheen, clearcoat, SSS, iridescence, volumetric, importance sampling, env maps |
| `atmosphere` | `atmosphere` | `atmosphere.rs` | 753 | (free functions + constants) | Atmospheric optics: Rayleigh/Mie scattering, sky color, air mass, optical depth, sunset model |
| `ai` | `ai` | `ai.rs` | 127 | (async functions) | AI-assisted optics queries via reqwest/tokio |
| `logging` | `logging` | `logging.rs` | 30 | (init function) | tracing-subscriber setup |

**Total**: ~11,840 lines across 19 source files.

## Design Principles

- Flat library crate -- no internal binaries
- Feature-gated modules -- consumers pull only what they need
- f64 precision throughout -- optics demands double precision
- No external FFT/linalg deps -- inline implementations for zero dependency cost
- `#[must_use]` on all pure functions, `#[non_exhaustive]` on all public enums
- `Cow<'static, str>` in error variants -- zero alloc on static messages
- Precomputed constants where possible (Rayleigh prefactor, 1/pi, CIE tables)
- `#[inline]` on hot-path functions -- every optics function that does math
- tracing instrumentation on functions that involve multi-step computation

## Data Flow

```
spectral ──> color science (wavelength <-> RGB, blackbody, CIE XYZ, SPD, CRI)
    |
    v
ray ──> geometric optics (Snell, Fresnel, dispersion, sequential trace)
    |         |
    |         +──> ray::system (paraxial trace, prescriptions, cardinal points)
    |         +──> ray::simulate (recursive trace, ray fans, spot diagrams, OPD)
    |
wave ──> wave optics (interference, diffraction, Fabry-Perot, AR coatings)
    |         |
    |         +──> wave::polarization (Stokes/Mueller formalism)
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

## Feature Independence

Each module compiles independently. No cross-feature imports except `error.rs`.

Intentional duplications to avoid cross-feature dependencies:

- `wave::pattern` contains its own `wavelength_to_rgb` -- avoids requiring `spectral` feature for pattern visualization
- `atmosphere` defines its own `RGB_WAVELENGTHS` constant -- avoids requiring `spectral` feature for sky color computation

This means a consumer can use `--features wave` without pulling in `spectral`, or `--features atmosphere` without pulling in `spectral` or `ray`.

Default features: `ray`, `wave`, `spectral`, `lens`. The `pbr`, `atmosphere`, `ai`, and `logging` features are opt-in.

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
| `reqwest` | `ai` | HTTP client for AI queries |
| `tokio` | `ai` | Async runtime |
| `serde_json` | `ai` | JSON serialization for AI payloads |
| `tracing-subscriber` | `logging` | Log output formatting and filtering |

### Dev-only

| Crate | Purpose |
|-------|---------|
| `criterion` | Benchmarking |
| `serde_json` | Test serialization round-trips |
