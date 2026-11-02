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
| `ray` | ray_core, ray_fresnel, ray_trace, ray_prt, ray_simulate, ray_system, ray_dispersion, ray_fiber, ray_grin | 2128 | `Medium`, `ComplexMedium`, `*Coefficients`, `TraceRay`, `OpticalSurface`, `PolarizedTraceHit`, `FresnelAmplitudes`, `PrtStep`, `PrtTraceHit`, `ParaxialRay`, `Prescription`, `GrinProfile` | Geometric optics: Snell, Fresnel (real + complex), dispersion (Sellmeier/Cauchy/Herzberger/Schott/Conrady, analytic derivatives and group index), chromatic aberration, fiber optics and fiber chromatic dispersion (material, waveguide, total; LP01 exact within weak guidance), sequential/recursive tracing with polarization, ray fans, spot diagrams, OPD; gradient-index tracing (SELFOC, sech, polynomial, Luneburg, fish-eye) by RK-Nyström with optical path |
| `spectral` | spectral_core, spectral_cie, spectral_photometry, spectral_fluor | 2473 | `Rgb`, `Xyz`, `Spd`, `Observer` (tag constants), `Donaldson` | Color science: wavelength↔RGB, Planck (numerically stable), Wien, CIE 1931/1964/2015 XYZ, SPD, illuminants, CRI, photometry (V(λ), luminous flux/efficacy); fluorescence (peaks, Stokes shift, quantum yield, the Donaldson re-radiation matrix, Förster transfer) |
| `wave` | wave_core, wave_polarization, wave_coherence, wave_airy, wave_fabry_perot, wave_diffraction, wave_doe, wave_zernike, wave_pattern, wave_beam, wave_focus | 5428 | `Polarization`, `StokesVector`, `MuellerMatrix` (16-f64 buffer), `Pattern2D`, `ZernikeWavefront`, `ThinFilmResult`, `AbsorbingFilmResult`, `GaussianBeam`, `VolumeGrating`, `FocusLens` | Wave optics: interference, coherence, Airy/Bessel, Fabry-Pérot, Fraunhofer/Fresnel diffraction, TMM (oblique s/p; complex-index R/T/A for absorbing layers), AR coatings, Jones/Stokes/Mueller, Zernike polynomials, 2D FFT patterns, PSF; Gaussian beams (q/ABCD, resonator modes, HG/LG modes, M², OAM); diffractive optics (grating equation, blazed/multilevel/binary/Dammann efficiencies, diffractive lens, Kogelnik volume holograms); vectorial focusing (Richards-Wolf: linear/elliptical, radial, azimuthal; uniform, Gaussian, doughnut or `GaussianBeam` pupil) |
| `lens` | lens.cyr, lens_aberration.cyr | 1628 | `CardinalPoints`, `SeidelCoefficients`, `WavefrontCoefficients`, `AberrSums` | Lens/mirror geometry: thin/thick lens, aberrations, wavefront coefficients from the Seidel sums, MTF (mono + poly + through-focus), aberrated OTF/MTF by pupil autocorrelation (Zernike or Seidel, sagittal/tangential, polychromatic complex sum), DoF, Petzval, multi-element; the Seidel and Buchdahl fifth-order sums of a `Prescription`, surface by surface, each fifth-order term split into intrinsic and induced |
| `pbr` | pbr_core, pbr_advanced | 1461 | (free functions) | PBR shading: Cook-Torrance, GGX, sheen, clearcoat, SSS, iridescence, volumetric, importance sampling, split-sum IBL |
| `atmosphere` | atmosphere.cyr | 381 | (free functions + constants) | Rayleigh/Mie scattering, King correction, sky color, air mass, optical depth, sunset model |
| `bridge` | bridge.cyr | 35 | (free functions) | Primitive-value cross-crate hooks (bijli/tara/badal) — no dependency on sibling crates |
| `serialize` | serialize.cyr | 80 | (free functions) | JSON roundtrips for the seven serde-tested types. ⚠ **Encode goes straight into a `str_builder`** since 2.3.4/2.3.5 — bayan is used only for Grisu2 float rendering; **decode** still walks a bayan value tree. The encoder emits round-trip-correct floats; bayan's decoder mis-rounds ~1 in 10⁵ doubles by 1 ULP (see serialize.cyr); every `*_from_json` reports via `err_out`; every wire format is pinned to exact bytes |
| `ai` | ai.cyr | 74 | `DaimonClient`, `DaimonConfig`, `HooshConfig` | AI-assisted optics queries via sandhi HTTP POST — **not in the core bundle** |

**Total**: 32 science modules + error, **14187 test assertions across 39 suites**,
185 benchmarks. (`tests/hardening.tcyr` is the cross-module regression suite for the
2.0.2 audit repairs and the 2.1.0 error channels — see those CHANGELOG entries.)

## Design Principles

- Flat library — free functions, heap-allocated structs (`alloc(sizeof(T))` +
  `#derive(accessors)`), enums as tag constants.
- `f64` precision throughout — optics demands double precision. Values are
  IEEE-754 bit patterns in `i64`; arithmetic goes through `f64_*`.
- **Correctness from the literature; fidelity to the Rust only where it was right.**
  Constants are exact ratios or IEEE-754 hex, and the ported modules keep the Rust's
  arithmetic (`powi` as square-and-multiply, left-associative folds) where it was correct. It
  often was not: 2.2.6 found six wrong formulas, five character for character in the Rust
  archive, 2.2.8 four more, and the 2.15.1 audit replaced ported arithmetic outright where it
  was wrong — the thin-film and polarizer formulas, the Fresnel integrals (now Cephes),
  `zernike_radial` (now a Jacobi recurrence), the through-focus MTF (now the defocused OTF).
  The modules added since 2.6.0 were never ports. Every value is pinned against an
  independent reference, never against the code's own output.
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
    |         └──> spectral_fluor (Stokes shift, quantum yield, Donaldson matrix, FRET)
    |
ray ──> geometric optics (Snell, Fresnel [real + complex], dispersion, trace, 3×3 polarization)
    |         ├──> ray_trace (polarization-aware sequential tracing, s/p tracking)
    |         ├──> ray_system (paraxial trace, prescriptions, cardinal points)
    |         └──> ray_simulate (recursive trace, ray fans, spot diagrams, OPD)
    |
wave ──> wave optics (interference, diffraction, Fabry-Pérot, AR coatings)
    |         ├──> wave_polarization (Stokes/Mueller formalism)
    |         ├──> wave_zernike (Zernike polynomials, Strehl ratio)
    |         ├──> wave_pattern (2D diffraction/PSF via hisab FFT)
    |         ├──> wave_beam (Gaussian beams: q/ABCD, HG/LG modes, M², OAM)
    |         └──> wave_focus (Richards-Wolf vectorial focusing; fields as ray_prt's PrtField)
    |
lens ──> lens geometry (thin/thick lens, aberrations, MTF, DoF, Petzval)
    |         └──> lens_aberration (Seidel and Buchdahl sums of a ray_system Prescription)
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
on ray_core+ray_fresnel; ray_prt and ray_grin on ray_core+ray_trace; wave_polarization on wave_core; pbr_advanced on pbr_core).
Across module groups there are five cross-dependencies: `lens`, whose aberrated
OTF (2.7.0) evaluates Zernike wavefronts from `wave_zernike` (bundled before it);
`lens_aberration` (2.14.0), which expands `ray_system`'s `Prescription` with its
`ParaxialRay`s;
`wave_focus` (2.13.0), which writes its fields into `ray_prt`'s `PrtField` and reads
a `wave_beam` beam;
`serialize`, which reads accessors from the ray/spectral/wave/lens type modules for
its JSON roundtrips (hence it is bundled last); and `wave_pattern`, which depends on
the hisab FFT dep.

## Module Independence

Each module is self-contained; the distlib bundler strips `include` lines and
resolves stdlib from the consumer's `[deps] stdlib`. Two intentional duplications
avoid cross-module coupling:

- `wave_pattern` carries its own `_pat_wl_to_rgb` — pattern visualization does not
  pull the full `spectral` module.
- `atmosphere` defines its own RGB wavelength constants — sky color does not pull
  `spectral`.

`[lib]` order (dependency order): `error` → the 9 ray modules → the 4 spectral →
the 11 wave → pbr_core, pbr_advanced → lens, lens_aberration → atmosphere → bridge →
serialize.
`[lib.ai]` is the same list plus `ai.cyr`.

### The `dist/*.deps` sidecars

Each bundle ships a `.deps` sidecar naming the stdlib folds that bundle needs in
scope. It is published metadata describing the artifact — **not** an enforced
contract: measured against cyrius 6.5.33 (re-measured for 2.2.4; unchanged
from the 6.5.20 measurement), a consumer that declares fewer folds
than the sidecar lists still resolves cleanly (`cyrius deps` exits 0), a bogus
fold name in the sidecar raises no error, and in the git-dep flow the vendored
set is driven by prakash's own `cyrius.cyml`, not by the sidecar.

Even so it should state the truth, and what `cyrius distlib` generates does not.

Through cyrius 6.6.17 the generator got prakash's layout backwards. The base bundle
(`[lib]` → `dist/prakash.cyr`) is the *narrow* math-only one and the profile
(`[lib.ai]`) the *wide* one, and the generator assumed the opposite: the core sidecar
advertised the sandhi/TLS stack and the ai sidecar under-reported.

⚠ **6.6.18 fixed that and opened a different hole.** A sidecar is now what a
compile-verify fixpoint proves the bundle needs, and the bundle opens with a matching
`# Requires` include block. But the verify runs with the package's named deps in scope,
and hisab is optional for prakash's consumers: its own requires block supplies `math`
and `ganita`, and `[deps.sakshi]` supplies sakshi. So at 6.7.18 the generator reports
`bayan` for core and `sandhi bayan` for ai, "compile-verified", and that is false for
the hisab-free link the README promises: the bundle with those leaves alone is
refused, on `SK_TRACE` and the unresolved `ganita_*` and `sakshi_*` calls, and
`scripts/check-consumer-link.sh` fails on it. With `[deps.hisab]` removed
from a scratch manifest, the same generator derives ten leaves (`math sakshi ganita
string alloc vec fnptr result bayan str`) and warns on `num_fft`.

**`scripts/sync-deps-sidecar.sh`** therefore regenerates both from the manifest: core =
declared stdlib minus the AI-only folds (`net http tls async random fdlopen dynlib
chrono sandhi`), ai = the full declared list. Two consequences are expected:
`cyrius distlib --check` reads STALE afterwards, and the `# Requires` block still names
only the generated leaves, so `include "dist/prakash.cyr"` alone does not compile. CI
enforces the sync, a core-bundle-is-TLS-free symbol scan and the consumer link, which
are the checks with teeth. bayan's include closure is 35 files, none of them net or TLS.

## Dependencies

| Dependency | Kind | Purpose |
|-----------|------|---------|
| `hisab` | git dep (tag 3.4.1) | FFT (`num_fft`) for `wave_pattern` (the suite also exercises `num_ifft`). `RayVec3` is layout-identical to hisab's `HVec3`, so all 26 `hvec3_*` ops work on prakash vectors unconverted — contract pinned by `tests/hisab_interop.tcyr` |
| `ganita` | stdlib | Transcendentals (acos/asin/atan2/pow/sinh/…) + linear algebra |
| `math` | stdlib | Comparisons, clamp/lerp/min/max, `F64_PI` etc., aarch64 polyfills |
| `bayan` | stdlib | JSON (`serialize` module) |
| `sakshi` | stdlib | Logging (trace diagnostics) |
| `sandhi` + TLS stack | stdlib | HTTP POST for the `ai` bundle only (net/http/tls/async/random/fdlopen/dynlib/chrono) |

The Rust build (`Cargo.toml`, criterion, bijli-backend, reqwest/tokio/serde) is
gone — see the [2.0.0] CHANGELOG entry. The `bridge` module replaces the former
`bijli-backend` feature with dependency-free primitive-value hooks.

## Consumers

Surveyed across every local repo at 2.5.1 (and kiran on GitHub, which has no local
clone). Each consumer is in a different state, so there is no single "consumers
move to Cyrius" milestone.

| Consumer | Language | How it depends | What it uses |
|----------|----------|----------------|--------------|
| **ranga** | Cyrius | git dep on `dist/prakash.cyr`, tag **2.2.8**, optional `spectral` profile only; omits hisab on purpose | spectral / colour: the `Xyz` bridge (read through the accessors and copied; its tests pin the 0/8/16 layout with raw `load64`), `spd_to_xyz`, `cct_from_xy`, `cie_cmf_at`, `spd_blackbody`; its tests also call the white points and `illuminant_d65`. Re-exports ~15 more flat names (`color_rendering_index`, `wavelength_to_rgb`, …) without calling them |
| **tanmatra** | Rust | crate `prakash = "1.1"` (locked 1.2.0) behind its `optics` feature | `Spd` (it builds Gaussian line SPDs as prakash `Spd`s, mapped to `spd_new` in its port plan), `Rgb` / `PrakashError`, `wavelength_to_rgb`. Its Cyrius port plan (not started) targets an opt-in `[lib.optics]` profile against the bundle |
| **soorat** | Rust | crate `prakash = "1"` (locked 1.2.0) behind its `optics` feature | `Rgb`, `color_temperature_to_rgb`, `wavelength_to_rgb`, `ior_to_f0`, `integrate_brdf_lut` (bakes its IBL LUT). Its Cook-Torrance / GGX / Fresnel are hand-copied into WGSL, which no bundle can reach. No Cyrius plan |
| **kiran** | Rust | crate `prakash` 1.2 (`pbr`, `atmosphere`) behind its `rendering` feature | re-exports six prakash modules wholesale; calls none of them. No Cyrius plan |

- **ranga's upgrade path is verified, not assumed:** `docs/guides/upgrading.md`
  (what moves between 2.2.8 and now), and `scripts/check-consumer-link.sh`, which CI
  runs to prove the core bundle links without hisab over the surface ranga uses
  (pinned toolchain, prakash's `lib/`; `--root` reproduces ranga's own setup).
- ⚠ **The three Rust consumers are frozen on 1.x**, and a `"1"` / `"1.2"`
  requirement never resolves to 2.x. They carry every value 2.x corrected — the
  2.2.6 six, the 2.4.x coefficient and constant repairs, 2.5.0's Seidel/LSA/OPD —
  and soorat's WGSL `distribution_ggx` still adds an absolute epsilon of the kind
  2.2.6 removed (`+ 1e-7` in its f32 port; prakash's was `1e-15`), which at its
  0.04 roughness floor leaves the specular peak ~4,900x low (25.6 against
  1/(pi a^2) = 1.24e5). Nothing prakash releases reaches them.
