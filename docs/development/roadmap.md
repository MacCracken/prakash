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

## Open work at a glance

**12 items. ⭐ Nothing is scheduled. Every item is demand-gated or waiting on a consumer —
build one only when someone asks for it.**

| # | Bucket | Item | State |
|---|---|---|---|
| 1–6 | 2.x | GRIN, DOE, Richards-Wolf, HG/LG beams, Buchdahl, aberrated MTF | demand-gated |
| 7–11 | 2.x | Fluorescence, non-linear, OAM, metamaterials, CIE 2006 observer | demand-gated |
| 12 | Blocked | soorat / kiran / ranga consume `dist/prakash.cyr` | waiting on the consumers |

⚠ **Items 1–11 are not a backlog anyone is working through.** Each is a subsystem,
listed so the scope boundary stays visible. Do not start one without a consumer
asking for it.

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
promises it will not force a minor bump. Reshuffling is free by construction.

**Nothing here depends on anything else here**, with one exception stated in
item 2. ⚠ **Keep rows SHORT.** A row says what the work is and what would block
it. Measurements belong in `CHANGELOG.md`; a row that grows past ~8 lines has
started duplicating the release history and should be cut back.

## 2.x — demand-gated: build when a consumer asks

Each is a subsystem, not an afternoon. None is speculative work worth doing before
someone needs it; all are listed so the scope boundary stays visible.

### Optics capability

- [ ] Gradient-index (GRIN) optics: curved ray paths through variable-n media
- [ ] Diffractive optical elements (DOE): phase gratings, holographic elements
- [ ] Vectorial diffraction (Richards-Wolf): high-NA focusing beyond scalar theory
- [ ] Hermite-Gaussian / Laguerre-Gaussian beam modes; M² beam quality
- [ ] Higher-order (5th-order Buchdahl) aberrations; wavefront coefficients from Seidel sums
- [ ] Aberrated MTF from generalized pupil-function autocorrelation

### Advanced

- [ ] Fluorescence (Stokes shift, excitation/emission spectra)
- [ ] Non-linear optics (SHG, Kerr) — if joshua needs it
- [ ] Orbital angular momentum (Laguerre-Gaussian modes)
- [ ] Metamaterials / negative refractive index
- [ ] Age-dependent CIE observer (CIE 2006)

## Blocked — not prakash's move

- [ ] soorat / kiran / ranga: consume `dist/prakash.cyr` directly once they move to Cyrius

## Constraints established by measurement — read before optimizing

- **An `Spd` holds POINT SAMPLES, not bin averages or bin integrals.**
  `_spd_integrate` is `sum(power[i] * cmf[i]) * 5` — a rectangle rule where the
  `* 5` is the bin width, so the consumer applies the integral once and the stored
  values are the function's value AT each wavelength. `spd_blackbody` and
  `spd_from_function` both point-sample. ⛔ **Do not integrate when building an
  Spd.** A bin integral over 5 nm is ~5x a point sample, and because `spd_to_xyz`
  normalises by X+Y+Z the chromaticity would still look right while every absolute
  photometric value was 5x wrong. The roadmap specified this incorrectly for six
  releases; see the [2.4.0] CHANGELOG entry.
- **hisab's quadrature applies to no prakash path.** 2.1.2 established that for
  `_spd_integrate`, `huygens_fresnel_1d` and `spd_blackbody`; 2.4.0 closed the last
  candidate, `spd_from_function`. `calc_integral_gauss5` is not used and there is
  no known place for it.

- **The small fixed-size scratch allocations are not worth converting, measured.**
  `alloc(16)` costs 6.4–7.5 ns on this host, so removing one or two only registers
  on a row whose baseline is a few hundred ns. `multilayer_rt` was the only site
  that moved (−6.8%, taken in 2.3.7). The five `spectral_cie` CRI sites are the
  biggest *byte* win left — 1,528 → 1,224 B/call, −19.9% of CRI's allocation
  volume — and measured **+0.04% / +0.24%** on the clock, i.e. the stack-local arm
  was marginally *slower*. `spectrum_strip` −0.68% at 24 B/call. Revisit only if
  the bytes matter for a named consumer, never for speed.

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
