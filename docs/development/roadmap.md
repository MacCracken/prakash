# Prakash Roadmap

> **Prakash** is the optics/light-simulation library, written in [Cyrius](https://github.com/MacCracken/cyrius). Math foundations (Complex + FFT) come from [hisab](https://github.com/MacCracken/hisab). Consumed by soorat (PBR shading), kiran (lighting), and ranga (lens effects).

**This document is forward-facing.** It lists what is not done yet. Shipped work
lives in `CHANGELOG.md`, which is the release history and the place to look for
what a given version changed and why. Nothing here is checked off; an item that
is finished is deleted from this file and described there instead.

What a completed item *leaves behind* is a constraint on future work — a
measurement that says "this approach does not pay". Those are kept, attached to
the open item they govern, so the next attempt does not repeat a failed one.

## Scope

Prakash owns the **physics of light**: how light travels, bends, reflects, interferes, diffracts, polarizes, and interacts with materials. It provides the math; consumers decide what to do with it (render pixels, simulate experiments, process images).

Prakash does NOT own:
- **Rendering pipeline** → kiran/soorat (they consume prakash for lighting math)
- **Image processing** → ranga (pixel operations, color spaces, filters)
- **Math primitives** → hisab (vectors, geometry, calculus, Complex, FFT)
- **Color science beyond spectral** → ranga (ICC profiles, gamut mapping)

## Performance

Ranked by measured impact against the 2.0.1 Rust-vs-Cyrius baseline — see
`docs/benchmarks-rust-v-cyrius.md`. Both remaining rows carry a caveat that cost
a measurement to establish.

- [ ] **`serialize/rgb_to_json` ~2×** (floor-corrected ~1,520 ns → ~3,400 ns).
      ⚠ **Not a regression to undo.** 85% of the cost is `bayan_json_v_build` and
      ~720 ns of that is *per float*, because bayan 1.2.1 replaced a 6-decimal
      renderer with round-trip-correct Grisu2. That bought bit-exact f64
      round-trips and fixed a real data-loss bug (the old encoder flushed any
      |x| < 5e-7 to zero). **Do not "fix" the float path.** The remaining lever is
      the 14% spent on object construction.
      ⚠ **The obvious version of that lever does not work.** The `_a` allocator
      variants (`bayan_json_v_obj_set_a` / `_build_a`) do exist in the vendored
      bayan, so the swap looks mechanical — but prakash allocates through a **bump
      allocator that never frees**, so "pass an arena" means nothing unless the
      arena is reused across calls, and `bayan_json_v_obj_new` / `_float_new` /
      `str_from` on the same path have no `_a` form to route through it. This is a
      question about allocator lifetime, not a find-and-replace.
- [ ] **Inline expansion** of the tiny `f64_*` wrappers and `_pbr_*`/`_lens_*`
      helpers — per-call overhead is the dominant term on the cheap ops
      (`pbr/fresnel_schlick` is 16.3× Rust at 17 ns absolute).
      ⚠ **There is no `#inline` in Cyrius**, and `cycc` silently accepts unknown
      attributes, so writing one would compile, lint clean, and do nothing. This
      means hand-expansion at every site: real churn against rows of 16–22 ns,
      squarely inside this host's 40% sub-100 ns spread. Any attempt must be
      measured **same-binary with interleaved arms**, never as a whole-suite delta.

### Constraints established by measurement — read before optimizing

- **SIMD is exhausted.** `pattern2d_normalized` is the only loop in prakash the
  `simd` fold applies to. The typed `f64v2_*`/`f64v4_*` wrappers are **slower than
  scalar**; only the raw `f64v_*` builtins help, and only as ONE bulk call;
  `f64v_scale` **overruns by one element on an odd count**; `memcpy`/`memset` are
  byte-at-a-time and 7–8× slower than a scalar loop. The remaining hot spots are
  unreachable from this fold: `interference_pattern` and `spd_blackbody` are
  transcendental-bound and there are **no vector transcendentals**; `spot_diagram`
  is branchy; `max_intensity` needs an `f64v_max` that does not exist;
  `_spd_integrate` reads an **interleaved** CMF table (stride 24). Real speedups
  there need algorithmic work or vector transcendentals upstream.
- **Do not route prakash's internals through `hvec3_*`.** `RayVec3` and `HVec3`
  are layout-identical and the interop contract is pinned by
  `tests/hisab_interop.tcyr` — but sending the tracer's dot products through
  `hvec3_dot` cost **2–11%** on `ray/trace_surface`, `ray/trace_sequential` and
  `ray/fresnel_unpolarized` across three runs, because it is two nested calls
  where the inline form is straight-line arithmetic. Revisit only if hisab gains
  an inlinable form.
- **Look for loop-invariant transcendentals before reaching for allocators.**
  The 2.3.1 win on `atm_sky_color_rgb` (−14%) was hoisting `_prk_cos` and both
  phase functions out of a per-channel loop; the roadmap row that had stood there
  for three releases blamed allocation, and `src/atmosphere.cyr` performs no
  allocation at all. Check the premise before optimizing against it.

## Accuracy & completeness (demand-gated — build on request)

- [ ] Gradient-index (GRIN) optics: curved ray paths through variable-n media
- [ ] Diffractive optical elements (DOE): phase gratings, holographic elements
- [ ] Vectorial diffraction (Richards-Wolf): high-NA focusing beyond scalar theory
- [ ] Hermite-Gaussian / Laguerre-Gaussian beam modes; M² beam quality
- [ ] Higher-order (5th-order Buchdahl) aberrations; wavefront coefficients (W_040, …) from Seidel sums
- [ ] Aberrated MTF from generalized pupil-function autocorrelation
- [ ] `spd_from_function(f, start_nm, end_nm)` — build an SPD from a continuous
      spectral function via hisab `calc_integral_gauss5`. **The one place hisab's
      quadrature genuinely fits.** The other candidates were investigated and do
      not: `_spd_integrate` is the CIE-defined weighted sum over the tabulated
      81-entry CMFs (the standard's method, not an approximation to improve),
      `huygens_fresnel_1d` integrates a caller-supplied discrete buffer with no
      continuous integrand, and `spd_blackbody` samples rather than integrates.
      `src/` contains **zero `fncall` sites**, so this is a new capability — it
      would be prakash's first function taking a callable integrand.

## Advanced / demand-gated

- [ ] Fluorescence (Stokes shift, excitation/emission spectra)
- [ ] Non-linear optics (SHG, Kerr) — if joshua needs it
- [ ] Orbital angular momentum (Laguerre-Gaussian modes)
- [ ] Metamaterials / negative refractive index
- [ ] Age-dependent CIE observer (CIE 2006)

## Housekeeping

- [ ] **Benchmark parity with `rust-old/`** — 180 Rust benches against 36 here;
      131 subjects uncovered. Mostly trivial scalar micro-benchmarks, and the
      expensive composites are already covered. Bulk-porting them would add noise
      to `bench-history.csv` without changing a decision — **do it only if a
      specific regression needs the resolution.**
- [ ] **Tracked probe files in the repo root.** `_probe_ai_bundle.cyr`,
      `_probe_bundle_fft.cyr` and `_probe_pow.tcyr` are committed at the top level,
      which `CLAUDE.md` forbids. All three landed in `15d900c` (2026-09-11) during
      the 2.2.9 toolchain bump. Triaged in 2.3.1 — **none is a test**: they print or
      `SYS_EXIT(0)` and assert nothing, so they cannot fail and CI never runs them.
      - `_probe_pow.tcyr` prints `_prk_pow` results beside Rust's hex. Fully
        superseded — `tests/hardening.tcyr` asserts that entire table live (251
        assertions). Safe to delete.
      - `_probe_ai_bundle.cyr` / `_probe_bundle_fft.cyr` include the built bundles
        and touch one symbol each, i.e. a *consumer-side link check* on
        `dist/prakash-ai.cyr` and `dist/prakash.cyr`. That check has value the other
        gates lack, but as written it proves nothing automatically. Either promote to
        a real `tests/*.tcyr` with assertions and let CI run it, or delete.
      Deleting is the maintainer's call, not a side effect of another release.

## Consumer integration

Blocked on the consumers, not on prakash.

- [ ] soorat / kiran / ranga: consume `dist/prakash.cyr` directly once they move to Cyrius

## Consumers

| Consumer | What it uses |
|----------|-------------|
| **soorat** | PBR shading (Cook-Torrance, Fresnel-Schlick) |
| **kiran** | Physically-based lighting math |
| **ranga** | Lens effects (DoF, chromatic aberration) |

## Boundary with Other Crates

| Feature | prakash | other |
|---------|---------|-------|
| Fresnel/Snell reflectance math | Yes (self-contained) | — |
| EM ↔ optics primitive bridge | Yes (`bridge`) | bijli (EM foundation) |
| Pixel-level image filter | — | ranga |
| 3D scene graph | — | kiran |
| Vector/matrix math, Complex, FFT | — | hisab |
| Color space conversion (ICC) | — | ranga |
| Spectral → RGB conversion | Yes | — |
| Polarization formalism (Jones/Stokes/Mueller) | Yes | — |
