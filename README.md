# Prakash

> **Prakash** (Sanskrit: प्रकाश — light, illumination) — optics and light simulation for AGNOS

Physics of light: ray optics, wave optics, spectral math, lens geometry, atmospheric scattering, and physically-based rendering primitives. Written in [Cyrius](https://github.com/MacCracken/cyrius), ported from Rust (2.0.0). Math foundations (Complex + FFT) come from [hisab](https://github.com/MacCracken/hisab).

Consumed by [soorat](https://github.com/MacCracken/soorat) (PBR shading), [kiran](https://github.com/MacCracken/kiran) (lighting), and [ranga](https://github.com/MacCracken/ranga) (lens effects) — still Rust for now; the `dist/prakash*.cyr` bundle is the interop surface until they port.

## Modules

| Module | Files | Description |
|--------|-------|-------------|
| **ray** | ray_core, ray_fresnel, ray_trace, ray_simulate, ray_system, ray_dispersion, ray_fiber | Geometric optics: Snell, Fresnel (real + complex refractive index), reflection (2D/3D), critical/Brewster angle, Beer-Lambert, TIR, 12 materials. Dispersion (Sellmeier/Cauchy/Herzberger/Schott/Conrady), chromatic aberration, fiber optics. Sequential + recursive ray tracer with polarization, ray fans, spot diagrams, OPD, paraxial system builder, cardinal points |
| **spectral** | spectral_core, spectral_cie, spectral_photometry | Wavelength↔RGB, Planck (numerically stable), Wien, color-temperature→RGB, photon energy. CIE 1931/1964/2015 observers, XYZ, SPD integration, illuminants (A/D50/D65/F2/F11, all verified against published CIE chromaticities and CCTs), **CIE 13.3 CRI** (Ra plus per-sample R_i, validated against published Ra for F2/F11/D65/A), photometry (V(λ), CIE 1951 V'(λ), luminous flux/efficacy) |
| **wave** | wave_core, wave_polarization, wave_coherence, wave_airy, wave_fabry_perot, wave_diffraction, wave_zernike, wave_pattern | Interference, single/double-slit + grating diffraction, thin-film reflectance, Malus, Jones/Stokes/Mueller, coherence, Airy/Bessel, Fabry-Pérot, Fraunhofer/Fresnel/Huygens diffraction, AR coatings + TMM, Zernike wavefronts, 2D FFT patterns, PSF |
| **lens** | lens | Thin/thick lens, lensmaker, mirrors, f-number/NA, FOV, MTF (mono/poly/through-focus), Seidel aberrations, Petzval, depth of field, multi-element systems |
| **pbr** | pbr_core, pbr_advanced | Cook-Torrance, GGX/Beckmann NDF, Smith geometry, Fresnel-Schlick, Lambert; anisotropy, sheen, clearcoat, SSS, iridescence, volumetric scattering, importance sampling, split-sum IBL |
| **atmosphere** | atmosphere | Rayleigh/Mie scattering, King factor, sky color, air mass (Kasten-Young), optical depth, sunset gradient |
| **bridge** | bridge | Primitive-value cross-crate hooks: bijli (EM ↔ wavelength/index), tara (stellar temp → RGB), badal (density/humidity → scattering) |
| **serialize** | serialize | JSON roundtrips (via bayan) for rgb, medium, sellmeier, lens type, polarization, prescription, SPD |
| **ai** *(opt-in)* | ai | Daimon/Hoosh AI client — blocking HTTP POST via sandhi. **Not in the core bundle** (pulls the TLS stack); ships in `dist/prakash-ai.cyr` |
| **error** | error | `PK_ERR_*` codes + `prakash_set_log_level` / trace logging (sakshi) |

## Quick Start

```toml
# cyrius.cyml
[package]
name     = "your-project"
version  = "${file:VERSION}"
language = "cyrius"
cyrius   = "6.5.33"

[deps]
# ganita provides the transcendentals (acos/asin/atan2/pow/sinh/…) and subsumes
# matrix/linalg. math stays for comparisons/clamp/lerp/min/max + polyfills.
stdlib = ["string", "fmt", "alloc", "vec", "str", "math", "ganita", "tagged", "fnptr"]

[deps.prakash]
git     = "https://github.com/MacCracken/prakash.git"
tag     = "2.2.4"
modules = ["dist/prakash.cyr"]        # math-only core (no TLS)
# For the AI client instead, pull the ai bundle (adds the sandhi HTTP/TLS stack):
# modules = ["dist/prakash-ai.cyr"]
```

```cyrius
include "lib/prakash.cyr"    # the vendored bundle

alloc_init();
var err = alloc(8);          # Result<T> -> err_out pointer (PK_ERR_NONE == 0)

# Snell's law: air (n=1.0) into glass (n=1.52) at 30 degrees
var angle = f64_div(F64_PI, f64_from(6));                 # 30 deg in radians
var n_glass = f64_div(f64_from(152), f64_from(100));      # 1.52
var refracted = ray_snell(F64_ONE, n_glass, angle, err);

# Fresnel reflectance at normal incidence (~4% off glass)
var reflectance = fresnel_normal(F64_ONE, n_glass);

# 550 nm light -> sRGB (Rgb struct: Rgb_r/_g/_b accessors)
var green = wavelength_to_rgb(f64_from(550), err);

# Where does a 50 mm lens focus an object at 2 m?
var image_dist = lens_thin_image_distance(f64_from(50), f64_from(2000), err);

# Malus's law: light through a polarizer at 45 degrees (-> 50%)
var transmitted = malus_law(F64_ONE, F64_PI_4);
```

`f64` values are IEEE-754 bit patterns in `i64`; do arithmetic through `f64_add`/`f64_mul`/`f64_sqrt`/… and compare with `f64_lt`/`f64_gt` (which return `1`/`0`). `Result<T>` becomes an `err_out` pointer (check `prakash_is_ok(load64(err))`); `Option<T>` returns `1`/`0` plus an out-pointer; tuples and `[f64; N]` return through caller-supplied buffers.

## Two Bundles

`cyrius distlib` emits two self-contained bundles — a consumer links **one**:

| Bundle | Contents | Pulls TLS? |
|--------|----------|-----------|
| `dist/prakash.cyr` | all optics + serialize (math-only core) | no |
| `dist/prakash-ai.cyr` | the same core **plus** the daimon/hoosh AI client | yes (sandhi) |

## Architecture

```
prakash (Cyrius)
  ├── hisab   (git dep, tag 2.11.2) — FFT (num_fft) for wave_pattern
  ├── ganita  (stdlib) — transcendentals + linear algebra
  ├── bayan   (stdlib) — JSON (serialize module)
  ├── sakshi  (stdlib) — logging (trace diagnostics)
  └── sandhi  (stdlib) — HTTP/TLS, ai bundle only
```

## Building

```sh
cyrius deps                 # resolve the hisab git dep
for f in tests/*.tcyr; do cyrius test "$f"; done   # 5615 assertions, 29 suites
cyrius distlib              # regenerate dist/prakash.cyr
cyrius distlib ai           # regenerate dist/prakash-ai.cyr
cyrius bench tests/prakash.bcyr                     # 36 benchmarks
./scripts/bench-history.sh  # append to the CSV history + benchmarks.md
```

The original Rust implementation is no longer in-tree. It was retained under
`rust-old/` through 2.2.2 as the translation reference and removed in 2.2.3, once
the port-completeness review (2.2.1) had confirmed nothing was left behind. It
remains in git: `git checkout 2.2.2 -- rust-old/`, or tag `1.2.0` for the
pre-port project at the repo root. See
[docs/benchmarks-rust-v-cyrius.md](docs/benchmarks-rust-v-cyrius.md) for a
Rust-vs-Cyrius performance comparison.

## License

GPL-3.0-only — see [LICENSE](LICENSE).
