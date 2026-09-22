# Architecture

prakash is a flat Cyrius library: flat `src/*.cyr` modules with free functions,
no internal binaries (`src/main.cyr` is a build smoke only). Ported from the Rust
source, which was retained under `rust-old/` through 2.2.2 and removed in 2.2.3
after the port-completeness review; recover it with
`git checkout 2.2.2 -- rust-old/`.

## Module Map

| Module | Files | Tests | Key Types | Purpose |
|--------|-------|-------|-----------|---------|
| `error` | error.cyr | (smoke) | `PK_ERR_*` codes | Shared error codes + `prakash_set_log_level` / `_prk_trace` (sakshi logging) |
| `ray` | ray_core, ray_fresnel, ray_trace, ray_simulate, ray_system, ray_dispersion, ray_fiber | 599 | `Medium`, `ComplexMedium`, `*Coefficients`, `TraceRay`, `OpticalSurface`, `PolarizedTraceHit`, `ParaxialRay`, `Prescription` | Geometric optics: Snell, Fresnel (real + complex), dispersion (Sellmeier/Cauchy/Herzberger/Schott/Conrady), chromatic aberration, fiber optics, sequential/recursive tracing with polarization, ray fans, spot diagrams, OPD |
| `spectral` | spectral_core, spectral_cie, spectral_photometry | 2025 | `Rgb`, `Xyz`, `Spd`, `Observer` (tag constants) | Color science: wavelength↔RGB, Planck (numerically stable), Wien, CIE 1931/1964/2015 XYZ, SPD, illuminants, CRI, photometry (V(λ), luminous flux/efficacy) |
| `wave` | wave_core, wave_polarization, wave_coherence, wave_airy, wave_fabry_perot, wave_diffraction, wave_zernike, wave_pattern | 1665 | `Polarization`, `StokesVector`, `MuellerMatrix` (16-f64 buffer), `Pattern2D`, `ZernikeWavefront`, `ThinFilmResult` | Wave optics: interference, coherence, Airy/Bessel, Fabry-Pérot, Fraunhofer/Fresnel diffraction, TMM (oblique s/p), AR coatings, Jones/Stokes/Mueller, Zernike polynomials, 2D FFT patterns, PSF |
| `lens` | lens.cyr | 239 | `CardinalPoints`, `SeidelCoefficients` | Lens/mirror geometry: thin/thick lens, aberrations, MTF (mono + poly + through-focus), DoF, Petzval, multi-element |
| `pbr` | pbr_core, pbr_advanced | 1404 | (free functions) | PBR shading: Cook-Torrance, GGX, sheen, clearcoat, SSS, iridescence, volumetric, importance sampling, split-sum IBL |
| `atmosphere` | atmosphere.cyr | 378 | (free functions + constants) | Rayleigh/Mie scattering, King correction, sky color, air mass, optical depth, sunset model |
| `bridge` | bridge.cyr | 28 | (free functions) | Primitive-value cross-crate hooks (bijli/tara/badal) — no dependency on sibling crates |
| `serialize` | serialize.cyr | 49 | (free functions) | JSON roundtrips for the seven serde-tested types. ⚠ **Encode goes straight into a `str_builder`** since 2.3.4/2.3.5 — bayan is used only for Grisu2 float rendering; **decode** still walks a bayan value tree. Floats are bit-exact; every `*_from_json` reports via `err_out`; every wire format is pinned to exact bytes |
| `ai` | ai.cyr | 30 | `DaimonClient`, `DaimonConfig`, `HooshConfig` | AI-assisted optics queries via sandhi HTTP POST — **not in the core bundle** |

**Total**: 25 science modules + error, **6757 test assertions across 31 suites**,
138 benchmarks. (`tests/hardening.tcyr` is the cross-module regression suite for the
2.0.2 audit repairs and the 2.1.0 error channels — see those CHANGELOG entries.)

## Design Principles

- Flat library — free functions, heap-allocated structs (`alloc(sizeof(T))` +
  `#derive(accessors)`), enums as tag constants.
- `f64` precision throughout — optics demands double precision. Values are
  IEEE-754 bit patterns in `i64`; arithmetic goes through `f64_*`.
- **Bit-fidelity to the Rust original** — constants encoded as exact ratios or
  IEEE-754 hex; `powi` replicated as square-and-multiply; left-associative fold
  order preserved so results match to the ULP.
  ⚠ **Deliberately NOT bit-faithful at six sites.** 2.2.6 found six wrong physics
  formulas, five character-for-character identical in the Rust archive, and fixed
  them — GGX at low roughness, both Fresnel-integral branches, the sphere vertex
  cap, the doublet prescription, third-order spherical, and `atm_air_mass` past
  the horizon. 2.2.8 added four more. Fidelity is the default, not the rule: the
  Rust original is a fidelity reference, never a correctness one.
- `#must_use` on pure functions.
- Errors are `PK_ERR_*` codes returned via an `err_out` pointer — never
  `unwrap`/`panic` in library code.
- Precomputed constants where possible (Rayleigh prefactor, 1/π).
  ⚠ The CIE tables are **not** static data — they are 21 lazily-built memoised
  allocations, which is why `prakash_reset_caches()` exists and why it MUST follow
  any `alloc_reset()`. See `docs/guides/allocation.md`.
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

⚠ **No module includes another — `src/*.cyr` contains no `include` statement at
all** (only `src/main.cyr`, the build entry, has any). The bundler concatenates
the `[lib]` list in dependency order and everything resolves at file scope, which
is why that order is load-bearing. Every module *uses* `error.cyr`'s `PK_ERR_*`
and `_prk_trace`, which is why it is bundled first. Within the
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

Each bundle ships a `.deps` sidecar naming the stdlib folds that bundle needs in
scope. It is published metadata describing the artifact — **not** an enforced
contract: measured against cyrius 6.5.33 (re-measured for 2.2.4; unchanged
from the 6.5.20 measurement), a consumer that declares fewer folds
than the sidecar lists still resolves cleanly (`cyrius deps` exits 0), a bogus
fold name in the sidecar raises no error, and in the git-dep flow the vendored
set is driven by prakash's own `cyrius.cyml`, not by the sidecar.

Even so it should state the truth, and by default it does not. prakash's bundle
layout is **inverted** relative to what `cyrius distlib` assumes: the base bundle
(`[lib]` → `dist/prakash.cyr`) is the *narrow* math-only one and the profile
(`[lib.ai]`) is the *wide* one, whereas since 6.5.10 the generator emits
`[deps] stdlib` ∪ include-scan for the base and a pruned inference for profiles —
i.e. it assumes the base is widest. Left alone, the core sidecar therefore
advertises the sandhi/TLS stack the math-only bundle never touches, and the ai
sidecar under-reports.

⚠ **The ai half of that was re-measured at 6.5.33, at 6.6.4 and at 6.6.6; the
figure keeps moving and the defect does not.** Through 6.5.20 the pruned inference
yielded literally `syscalls io`. At 6.5.33 it yielded ten folds — `string alloc
str vec math ganita tagged bayan sandhi sakshi`. At **6.6.4** (measured for 2.3.0)
it yielded **25** at 2.3.0 and yields **26** at 2.4.0 (the extra leaf is `fnptr`,
pulled in by `spd_from_function`'s `fncall1` — prakash's first callback site), and
the generator now reports its own patching — `sidecar:
re-added 15 leaf(s) the inference missed (compile-verified)` — but the inversion
survives all of it: the same run emits **29** leaves for the math-only base
against those 25 for ai, i.e. the narrow bundle is still advertised as the wider
one, complete with the sandhi/TLS stack it never touches. At **6.6.6** (2.4.1) the
same 29-for-base / 26-for-ai pair reproduces. The generator still
under-reports for the wide profile; only the size of the shortfall moved. The core
over-reporting reproduces unchanged. **`scripts/sync-deps-sidecar.sh`** regenerates
both from the manifest — core = declared stdlib minus the AI-only folds
(`net http tls async random fdlopen dynlib chrono sandhi`), ai = the full declared
list. CI enforces the sync plus a core-bundle-is-TLS-free symbol scan, which is
the check that actually has teeth.

## Dependencies

| Dependency | Kind | Purpose |
|-----------|------|---------|
| `hisab` | git dep (tag 3.2.1) | FFT (`num_fft`) for `wave_pattern` (the suite also exercises `num_ifft`). `RayVec3` is layout-identical to hisab's `HVec3`, so all 26 `hvec3_*` ops work on prakash vectors unconverted — contract pinned by `tests/hisab_interop.tcyr` |
| `ganita` | stdlib | Transcendentals (acos/asin/atan2/pow/sinh/…) + linear algebra |
| `math` | stdlib | Comparisons, clamp/lerp/min/max, `F64_PI` etc., aarch64 polyfills |
| `bayan` | stdlib | JSON (`serialize` module) |
| `sakshi` | stdlib | Logging (trace diagnostics) |
| `sandhi` + TLS stack | stdlib | HTTP POST for the `ai` bundle only (net/http/tls/async/random/fdlopen/dynlib/chrono) |

The Rust build (`Cargo.toml`, criterion, bijli-backend, reqwest/tokio/serde) is
gone — see the [2.0.0] CHANGELOG entry. The `bridge` module replaces the former
`bijli-backend` feature with dependency-free primitive-value hooks.
