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
      ⚠ **2.3.2's note on this row was WRONG and is corrected here.** It claimed
      `bayan_json_v_obj_new` / `_float_new` "have no `_a` form to route through"
      an arena. They do. bayan ships **15** `_a` variants covering every
      constructor this path uses — `bayan_json_v_obj_new_a`,
      `bayan_json_v_float_new_a`, `bayan_json_v_obj_set_a`, `bayan_json_v_build_a`
      — and each plain form is literally a wrapper passing `default_alloc()`.
      **Nothing is missing upstream; there is no bayan issue to file.** This is an
      internal repair.
      ⛔ **MEASURED IN 2.3.3, AND THE ARENA IS NOT WORTH TAKING. Do not attempt it
      again without new evidence.** Two independent measurements agreed:
      construction is **10.9–11.1%** of the row (not 14%), an arena captures only
      ~28% of that, and the whole prize is **−1.8% to −4.4%** — at or below the
      per-arm jitter. Even if construction were FREE the ceiling is **−11%**.
      ⛔ **An intermediate claim of mine was also wrong and is retracted here.** I
      wrote that `alloc()` "takes a LOCK on every call". It does not:
      `_alloc_lock_acquire` (lib/alloc.cyr) early-returns while
      `_threads_active == 0`, and prakash never arms threads — `src/` contains no
      `thread_create`. There is no lock to remove, and the per-allocation delta
      between `alloc()` and `arena_alloc` is ~1–3 ns. Going through an allocator
      *handle* (`alloc_via`, an indirect `fncall2`) is measurably **worse** than
      calling `alloc()` directly.
      ⛔ **And the arena breaks the lifetime contract.** Reproduced, not theorised:
      `bayan_json_v_build_a` puts BOTH the `Str` header and its byte buffer in the
      passed allocator, so an `arena_reset` at the top of the next call hands back
      the identical pointers and a previously returned string silently becomes the
      new result. Silent corruption, no crash — and the existing benchmark discards
      its result, so it cannot see it. `tests/serialize.tcyr` already holds one
      `to_json` result across a later `to_json` call.
      ⭐ **The lever that IS large was never costed: skip the bayan value tree.**
      Emitting the same bytes straight through a `str_builder` measured **−19.6%**,
      and hand-assembling into a stack buffer **−28.4%** — the floor with Grisu2
      untouched. That is the shape any future attempt on this row should take.
      ⭐ **"Do not fix the float path" is now proven rather than asserted:**
      `bayan_f64_to_json` is 671 ns against 66 ns for the 6-decimal `fmt_float_buf`
      it replaced, same binary — 3 floats × ~605 ns IS the ~1.9× step, and that
      renderer is what makes f64 round-trips bit-exact.

- [ ] **`trace_surface` boxes 5 structs per call** (128 bytes/call, measured), and
      it runs **2× per `trace_sequential`** (635 ns) and **38× per `spot_diagram`**
      (19.0 µs). A controlled same-binary A/B on a faithful clone: `alloc()` 206 ns
      (control, matching the shipped 210 ns), a caller-supplied flat buffer
      **132 ns (−36%)** — projecting to roughly **−23% on `trace_sequential`** and
      **−15% on `spot_diagram`**. This is CLAUDE.md's "write into a caller buffer"
      rule applied to the hottest allocating path in the tracer.
      ⚠ It changes an internal calling convention, so it is a bite of its own, and
      the projections above are projections — re-measure on the real functions.
- [ ] **`_prk_trace` runs `strlen` on a compile-time literal**, 29 ns per traced
      entry (20 ns of it `strlen`) across **32 call sites** — about **14% of
      `trace_surface`**. The label is always a literal, so the length is known at
      the call site. ⚠ sakshi's `sakshi_trace(name, len)` already takes the length;
      the waste is entirely on prakash's side of the wrapper.
- [ ] **Four sites allocate a per-call scratch buffer**, which CLAUDE.md's DO-NOT
      list calls a leak under a bump allocator that never frees. Identified during
      the 2.3.3 investigation; each needs its own look, since the fix is either a
      caller buffer or a stack local depending on lifetime.
- [ ] **`spd_to_json` has no benchmark row and is the heaviest serializer.** It
      allocates *inside a loop* — one bayan float node plus one array push per
      sample, so ~81 for a standard SPD against `rgb_to_json`'s 3, on the same
      pattern the whole `rgb_to_json` row is about. Measured arena footprint 6,808
      bytes/call and unbounded in `Spd_len`. Bench it before optimising anything
      else in `serialize`; the row that is measured is not the row that costs.

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
