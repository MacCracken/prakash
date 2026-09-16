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

- [ ] **Ten per-call scratch allocations, not four — in two tiers wanting opposite
      fixes.** Measured in 2.3.5 with `alloc_used()` deltas.
      ⭐ **Tier 1 — large, variable-size, all in `wave_pattern.cyr`** (`grid` at
      :269 and :549, `aperture` at :299, `col_buf` at :220). Scratch is the
      MAJORITY of these functions' allocation: `diffraction_pattern_2d(64×64)`
      hands out 99,352 B/call of which **66,560 B (67%)** is scratch;
      `diffraction_pattern_circular(32)` **75%**; `psf_from_wavefront(32×32)`
      **67%**. A grow-only module-scope cache cut those to 536 / 32,792 / 8,216
      B/call, bit-identical across every cell.
      ⚠ **And bought NO measurable time** — 3 interleaved same-binary rounds, all
      inside noise. **These are leak fixes, not speed fixes**; justify them as
      memory or not at all.
      ⭐ **Tier 2 — small, fixed-size** (`pbr_advanced.cyr:428/436`,
      `spectral_cie.cyr:3117/3126/3127` and `:3156/3157`,
      `wave_diffraction.cyr:497/501`, `wave_pattern.cyr:493` and `:517`). Only
      `multilayer_rt` moved measurably: **493 → 462 ns (−6.3%)**, reproduced 3/3.
      ⭐ **The stack-local budget is 122,880 BYTES**, bisected under cyrius 6.6.4:
      `var b[N]` compiles clean at N = 122,864 and trips *"oversized array local
      kept in shared global"* at N = 122,872 — so `var X[N]` is N **bytes**, not N
      slots. Every Tier 2 site is at 0.04% of budget, so stack locals are
      unconditionally safe there. `col_buf` is safe only for nh ≤ 7,679 while
      `_pat_dims_bad` admits nh up to 8,192, so it does **not** cover the domain;
      `grid` and `aperture` are never stack-eligible.
      ⚠ **Correction:** `trace_surface`'s 5 boxed structs are **not** scratch — all
      five escape through the returned `TraceHit`. That is the row above, and no
      stack local can fix it.

### Housekeeping

- [ ] **Benchmark parity with `rust-old/`** — 180 Rust benches against 36 here;
      131 subjects uncovered. Mostly trivial scalar micro-benchmarks, and the
      expensive composites are already covered. Bulk-porting them would add noise
      to `bench-history.csv` without changing a decision — **do it only if a
      specific regression needs the resolution.**

## 2.4.x — minor: adds public API

⚠ **Readiness checked in 2.3.3.** `calc_integral_gauss5` is **ready**: it is one of
the few integral forms hisab 3.x did **not** move onto `Result<T,E>`, and it links
from prakash's existing include set (compiled and run). So the plumbing is not the
gap. ⛔ **The real gap is a semantics question this row never asked:** whether an
SPD sample means the *average over its bin* or a *point sample at its wavelength*.
Integrating when the CIE convention wants point samples is wrong by a factor of
about the step width — roughly **5×** at 5 nm. Settle that before writing code.


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

- **Bayan's value tree costs far more in STRING bytes than in allocations.**
  `_jb_append_string` appends string content **one byte at a time** through
  `str_builder_add_cstr_a` with a 2-byte buffer — a strlen, grow-check and memcpy
  call per character, for every key and every string value. That, not the node
  allocations, is the dominant cost: allocation count explains only ~11–25% of the
  time 2.3.5 recovered. It is why the win tracked string content almost exactly —
  `spd_to_json` (3 short keys, 81 floats) −5.8%, `rgb_to_json` (3 one-char keys)
  −20.5%, `medium_to_json` (10 string bytes) −34.3%, `prescription_to_json`
  (255 string bytes at 6 surfaces) −34.5%. ⭐ If a future document is
  string-heavy, expect a large win; if it is float-heavy, expect a small one.
- **`str_builder_add_json_str` is byte-identical to bayan's escaping** — verified
  across all 255 reachable byte values, 0 mismatches. Use it rather than
  hand-rolling escapes.
- **CLAUDE.md's "write into a caller buffer" rule has a precondition the rule does
  not state: the caller must own the lifetime.** Applying it to `trace_surface`
  produces aliasing corruption, because `trace_sequential` retains hits in a
  returned vec. Check where the pointer ends up before reaching for a caller
  buffer.

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
