# Architecture

prakash is a flat Cyrius library: flat `src/*.cyr` modules with free functions,
no internal binaries (`src/main.cyr` is a build smoke only). Ported from the Rust
source retained under `rust-old/`.

## Module Map

| Module | Files | Tests | Key Types | Purpose |
|--------|-------|-------|-----------|---------|
| `error` | error.cyr | (smoke) | `PK_ERR_*` codes | Shared error codes + `prakash_set_log_level` / `_prk_trace` (sakshi logging) |
| `ray` | ray_core, ray_fresnel, ray_trace, ray_simulate, ray_system, ray_dispersion, ray_fiber | 593 | `Medium`, `ComplexMedium`, `*Coefficients`, `TraceRay`, `OpticalSurface`, `PolarizedTraceHit`, `ParaxialRay`, `Prescription` | Geometric optics: Snell, Fresnel (real + complex), dispersion (Sellmeier/Cauchy/Herzberger/Schott/Conrady), chromatic aberration, fiber optics, sequential/recursive tracing with polarization, ray fans, spot diagrams, OPD |
| `spectral` | spectral_core, spectral_cie, spectral_photometry | 1691 | `Rgb`, `Xyz`, `Spd`, `Observer` (tag constants) | Color science: wavelength↔RGB, Planck (numerically stable), Wien, CIE 1931/1964/2015 XYZ, SPD, illuminants, CRI, photometry (V(λ), luminous flux/efficacy) |
| `wave` | wave_core, wave_polarization, wave_coherence, wave_airy, wave_fabry_perot, wave_diffraction, wave_zernike, wave_pattern | 1568 | `Polarization`, `StokesVector`, `MuellerMatrix` (16-f64 buffer), `Pattern2D`, `ZernikeWavefront`, `ThinFilmResult` | Wave optics: interference, coherence, Airy/Bessel, Fabry-Pérot, Fraunhofer/Fresnel diffraction, TMM (oblique s/p), AR coatings, Jones/Stokes/Mueller, Zernike polynomials, 2D FFT patterns, PSF |
| `lens` | lens.cyr | 145 | `CardinalPoints`, `SeidelCoefficients` | Lens/mirror geometry: thin/thick lens, aberrations, MTF (mono + poly + through-focus), DoF, Petzval, multi-element |
| `pbr` | pbr_core, pbr_advanced | 865 | (free functions) | PBR shading: Cook-Torrance, GGX, sheen, clearcoat, SSS, iridescence, volumetric, importance sampling, split-sum IBL |
| `atmosphere` | atmosphere.cyr | 313 | (free functions + constants) | Rayleigh/Mie scattering, King correction, sky color, air mass, optical depth, sunset model |
| `bridge` | bridge.cyr | 28 | (free functions) | Primitive-value cross-crate hooks (bijli/tara/badal) — no dependency on sibling crates |
| `serialize` | serialize.cyr | 24 | (free functions) | JSON roundtrips (bayan) for the seven serde-tested types |
| `ai` | ai.cyr | 12 | `DaimonClient`, `DaimonConfig`, `HooshConfig` | AI-assisted optics queries via sandhi HTTP POST — **not in the core bundle** |

**Total**: 25 science modules + error, **5251 test assertions across 26 suites**,
27 benchmarks.

## Design Principles

- Flat library — free functions, heap-allocated structs (`alloc(sizeof(T))` +
  `#derive(accessors)`), enums as tag constants.
- `f64` precision throughout — optics demands double precision. Values are
  IEEE-754 bit patterns in `i64`; arithmetic goes through `f64_*`.
- **Bit-fidelity to the Rust original** — constants encoded as exact ratios or
  IEEE-754 hex; `powi` replicated as square-and-multiply; left-associative fold
  order preserved so results match to the ULP.
- `#must_use` on pure functions.
- Errors are `PK_ERR_*` codes returned via an `err_out` pointer — never
  `unwrap`/`panic` in library code.
- Precomputed constants where possible (Rayleigh prefactor, 1/π, CIE tables).
- Optional cost is opt-in: the AI client (and its TLS stack) lives in a separate
  bundle so math-only consumers pay nothing for it.

## Data Flow

```
spectral ──> color science (wavelength <-> RGB, blackbody, CIE XYZ, SPD, CRI)
    |
ray ──> geometric optics (Snell, Fresnel [real + complex], dispersion, trace)
    |         ├──> ray_trace (polarization-aware sequential tracing, s/p tracking)
    |         ├──> ray_system (paraxial trace, prescriptions, cardinal points)
    |         └──> ray_simulate (recursive trace, ray fans, spot diagrams, OPD)
    |
wave ──> wave optics (interference, diffraction, Fabry-Pérot, AR coatings)
    |         ├──> wave_polarization (Stokes/Mueller formalism)
    |         ├──> wave_zernike (Zernike polynomials, Strehl ratio)
    |         └──> wave_pattern (2D diffraction/PSF via hisab FFT)
    |
lens ──> lens geometry (thin/thick lens, aberrations, MTF, DoF, Petzval)
pbr  ──> PBR shading (Cook-Torrance, GGX; advanced: sheen/clearcoat/SSS/iridescence)
atmosphere ──> sky models (Rayleigh/Mie, air mass, optical depth, sunset)
bridge ──> primitive-value cross-crate hooks (bijli/tara/badal)
    |
error ──> PK_ERR_* codes + logging, shared by all modules
```

Every module includes `error.cyr` (for `PK_ERR_*` and `_prk_trace`). Within the
ray and wave groups, later files build on earlier ones (ray_trace/simulate/system
on ray_core+ray_fresnel; wave_polarization on wave_core; pbr_advanced on pbr_core).
Across module groups there are two cross-dependencies: `serialize`, which reads
accessors from the ray/spectral/wave/lens type modules for its JSON roundtrips
(hence it is bundled last), and `wave_pattern`, which depends on the hisab FFT dep.

## Module Independence

Each module is self-contained; the distlib bundler strips `include` lines and
resolves stdlib from the consumer's `[deps] stdlib`. Two intentional duplications
avoid cross-module coupling:

- `wave_pattern` carries its own `_pat_wl_to_rgb` — pattern visualization does not
  pull the full `spectral` module.
- `atmosphere` defines its own RGB wavelength constants — sky color does not pull
  `spectral`.

`[lib]` order (dependency order): `error` → the 7 ray modules → the 3 spectral →
the 8 wave → pbr_core, pbr_advanced → lens → atmosphere → bridge → serialize.
`[lib.ai]` is the same list plus `ai.cyr`.

### The `dist/*.deps` sidecars

Each bundle ships a `.deps` sidecar naming the stdlib folds a consumer must have
in scope; `cyrius deps` **validates** a consumer's `[deps] stdlib` against it and
hard-errors on anything missing, so it is a contract rather than documentation.

prakash's bundle layout is **inverted** relative to what `cyrius distlib` assumes:
the base bundle (`[lib]` → `dist/prakash.cyr`) is the *narrow* math-only one and
the profile (`[lib.ai]`) is the *wide* one, whereas the generator treats the base
as widest and prunes profiles. Left alone it therefore over-reports for the core
bundle (handing it the sandhi/TLS stack it never touches) and under-reports for
the ai bundle. prakash regenerates both from the manifest with
**`scripts/sync-deps-sidecar.sh`** — core = declared stdlib minus the AI-only
folds (`net http tls async random fdlopen dynlib chrono sandhi`), ai = the full
declared list. CI enforces both the sync and a core-bundle-is-TLS-free scan.

## Dependencies

| Dependency | Kind | Purpose |
|-----------|------|---------|
| `hisab` | git dep (tag 2.11.1) | Complex + FFT (`num_fft`/`num_ifft`) for `wave_pattern` |
| `ganita` | stdlib | Transcendentals (acos/asin/atan2/pow/sinh/…) + linear algebra |
| `math` | stdlib | Comparisons, clamp/lerp/min/max, `F64_PI` etc., aarch64 polyfills |
| `bayan` | stdlib | JSON (`serialize` module) |
| `sakshi` | stdlib | Logging (trace diagnostics) |
| `sandhi` + TLS stack | stdlib | HTTP POST for the `ai` bundle only (net/http/tls/async/random/fdlopen/dynlib/chrono) |

The Rust build (`Cargo.toml`, criterion, bijli-backend, reqwest/tokio/serde) is
gone — see the [2.0.0] CHANGELOG entry. The `bridge` module replaces the former
`bijli-backend` feature with dependency-free primitive-value hooks.
