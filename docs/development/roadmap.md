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

## How this file is organised

Items are grouped by the **release class their change implies**, not by the order
anyone intends to do them:

| Bucket | Means | Rule |
|---|---|---|
| **2.3.x — patch** | no public API moves | internals, perf, tests, docs, tooling |
| **2.4.x — minor** | adds public API | new entry points, new capability |
| **2.x — demand-gated** | adds a subsystem | build when a consumer actually asks |
| **Blocked** | not prakash's move | waiting on something external |

⚠ **The bucket is a SemVer classification, not a queue.** Anything in 2.3.x can
ship in any order, in any patch release, in any combination — the bucket only
promises it will not force a minor bump. Same within 2.4.x. So reshuffling is
free by construction, and the only thing that moves an item between buckets is a
change in what it does to the public surface.

**Nothing here has a dependency on anything else here.** If that ever stops being
true, say so in the item itself rather than relying on list order.

## 2.3.x — patch: no public API change

### Performance

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

### Housekeeping

- [ ] **Benchmark parity with `rust-old/`** — 180 Rust benches against 36 here;
      131 subjects uncovered. Mostly trivial scalar micro-benchmarks, and the
      expensive composites are already covered. Bulk-porting them would add noise
      to `bench-history.csv` without changing a decision — **do it only if a
      specific regression needs the resolution.**

## 2.4.x — minor: adds public API

- [ ] `spd_from_function(f, start_nm, end_nm)` — build an SPD from a continuous
      spectral function via hisab `calc_integral_gauss5`. **The one place hisab's
      quadrature genuinely fits.** The other candidates were investigated and do
      not: `_spd_integrate` is the CIE-defined weighted sum over the tabulated
      81-entry CMFs (the standard's method, not an approximation to improve),
      `huygens_fresnel_1d` integrates a caller-supplied discrete buffer with no
      continuous integrand, and `spd_blackbody` samples rather than integrates.
      `src/` contains **zero `fncall` sites**, so this is a new capability — it
      would be prakash's first function taking a callable integrand.

## 2.x — demand-gated: build when a consumer asks

Each of these is a subsystem, not an afternoon. None is speculative work worth
doing before someone needs it; all are listed so the scope boundary stays visible.

### Optics capability

- [ ] Gradient-index (GRIN) optics: curved ray paths through variable-n media
- [ ] Diffractive optical elements (DOE): phase gratings, holographic elements
- [ ] Vectorial diffraction (Richards-Wolf): high-NA focusing beyond scalar theory
- [ ] Hermite-Gaussian / Laguerre-Gaussian beam modes; M² beam quality
- [ ] Higher-order (5th-order Buchdahl) aberrations; wavefront coefficients (W_040, …) from Seidel sums
- [ ] Aberrated MTF from generalized pupil-function autocorrelation

### Advanced

- [ ] Fluorescence (Stokes shift, excitation/emission spectra)
- [ ] Non-linear optics (SHG, Kerr) — if joshua needs it
- [ ] Orbital angular momentum (Laguerre-Gaussian modes)
- [ ] Metamaterials / negative refractive index
- [ ] Age-dependent CIE observer (CIE 2006)

## Blocked — not prakash's move

Blocked on the consumers, not on prakash.

- [ ] soorat / kiran / ranga: consume `dist/prakash.cyr` directly once they move to Cyrius

## Constraints established by measurement — read before optimizing

- **Hand-inlining pays at this scale, and only same-binary A/B can show it.**
  2.3.2 took `pbr/fresnel_schlick` 16.3 → 12 ns (−26%), `distribution_ggx`
  22 → 17 (−23%) and `cook_torrance` 93 → 79.7 (−14%) by writing out `f64_clamp`
  and `_pbr_pow5` at the hot sites and turning eighteen per-call constant
  divisions into hex literals. Cyrius has no `#inline`, so a helper call is a real
  call. ⚠ Two rules came out of it: **a constant helper must return a literal**
  (`_pbr_ln2`'s comment has always said so), and **every such change must be
  pinned by an assertion that recomputes the original expression bit-for-bit** —
  `tests/constants.tcyr` exists for that. Expanding a helper anywhere that is not
  measurably hot is churn; the delegating form stays the default.
- **Test the source, not the bundle, when pinning something in `src/`.**
  `tests/constants.tcyr` first included `dist/prakash.cyr` and would have happily
  checked a stale artifact's digits after a `src/` edit — passing on exactly the
  change it exists to catch. Bundle inclusion is right for one question only,
  *does the shipped artifact link*, which is `tests/ai_bundle.tcyr`'s job.

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
