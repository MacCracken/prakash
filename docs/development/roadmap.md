# Prakash Roadmap

> **Prakash** is the optics/light-simulation library, written in [Cyrius](https://github.com/MacCracken/cyrius). Its one external math call is hisab's FFT (`num_fft`); complex arithmetic and everything else is its own. Consumed by ranga (spectral/colour, on Cyrius) and, through the frozen Rust 1.x crate, tanmatra, soorat and kiran — see the consumers section of `docs/architecture/overview.md`.

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
- **Rendering pipeline** → kiran/soorat (soorat bakes its IBL LUT with prakash and shades in its own WGSL; kiran re-exports prakash modules)
- **Image processing** → ranga (pixel operations, color spaces, filters)
- **Math primitives** → hisab (vectors, geometry, calculus, FFT). prakash calls only `num_fft`; its complex arithmetic is its own
- **Color science beyond spectral** → ranga (ICC profiles, gamut mapping)

## Open work at a glance

**12 items, all demand-gated or waiting on a consumer.** No known defect is open.

| # | Bucket | Item | State |
|---|---|---|---|
| 1–6 | 2.x | GRIN, DOE, Richards-Wolf, HG/LG beams, Buchdahl, aberrated MTF | demand-gated, no ask |
| 7–11 | 2.x | Fluorescence, non-linear, OAM, metamaterials, CIE 2006 observer | demand-gated, no ask |
| 12 | Consumers | ranga's bump from 2.2.8; tanmatra's planned `[lib.optics]` profile | their move; prakash's side done |

⭐ **Demand was measured at 2.5.1, not assumed.** A survey of every local repo (plus
kiran on GitHub) found **no consumer asking for any of items 1–11**. Item 7's only
trace is a stale registry line in agnosticos written about prakash 1.1; item 8's
named requester, joshua (GitHub-only, Rust 0.1.0), has no prakash dependency of its
own — it reaches the 1.x crate only through kiran, behind its optional `engine` feature — and
nothing in it mentions non-linear optics. Implicit signals, none of them a request:
tanmatra builds Gaussian line SPDs itself (see the `spd_to_xyz` caveat), ranga
inlines `xyz_to_xyy` to avoid an allocation, soorat carries a roughness-aware
Fresnel in WGSL. Re-survey before starting any item.

⚠ **Items 1–11 are not a backlog anyone is working through.** Each is a subsystem,
listed so the scope boundary stays visible. Do not start one without a consumer
asking for it.

## How this file is organised

Items are grouped by the **release class their change implies**, not by the order
anyone intends to do them:

| Bucket | Means | Rule |
|---|---|---|
| **2.5.x — patch** | no public API moves | internals, perf, tests, docs, tooling, data fixes |
| **2.6.x — minor** | adds public API | new entry points, new capability |
| **2.x — demand-gated** | adds a subsystem | build when a consumer actually asks |
| **Consumers** | not prakash's move | waiting on a consumer's own change |

⚠ **The bucket is a SemVer classification, not a queue.** Anything in 2.5.x can
ship in any order, in any patch release, in any combination — the bucket only
promises it will not force a minor bump. Reshuffling is free by construction.

**Nothing here depends on anything else here.** ⚠ **Keep rows SHORT.** A row says
what the work is and what would block it. Measurements belong in `CHANGELOG.md`;
a row that grows past ~8 lines has started duplicating the release history and
should be cut back.

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

## Consumers — their move

- [ ] **ranga** already consumes `dist/prakash.cyr` (Cyrius, tag 2.2.8, `spectral`
  profile). Bumping the tag is ranga's change; the path is verified in
  `docs/guides/upgrading.md`. CI holds the hisab-free link on the pinned toolchain
  via `scripts/check-consumer-link.sh`; ranga's exact 6.6.2 setup is `--root`.
- [ ] **tanmatra** plans an opt-in `[lib.optics]` profile in its Cyrius port (not
  started). If its M6 line-spectrum grid needs colour, a line-spectrum → XYZ entry
  point is a 2.6.x minor — only when it asks.
- soorat and kiran are Rust 1.x with no port plan, so there is no item for them.

## Constraints established by measurement — read before optimizing

- **`SellmeierCoefficients` holds exactly THREE terms, and that is now a known
  boundary rather than an assumption.** 2.4.3 needed Daimon & Masumura's four-term
  water fit and did not widen the struct: D&M's two adjacent UV resonances (135 and
  162 nm) are carried by one least-squares merged term at 138.3 nm, which tracks the
  published four-term curve to **3.5e-6 in n** across 380-780 nm. So a four-term
  published fit is not a reason to move public API — measure the merge first. The
  term that cannot be dropped is the **IR** one: absent, it costs ~0.0014 in n and
  17% of the dispersion (that was the pre-2.4.3 defect).
- ⚠ **Documentation counts drift silently and must be MEASURED, not read.** 2.4.7 found
  every per-module assertion count in `docs/architecture/overview.md` stale, pbr by 539
  and spectral by 334, and the totals in two files stale since 2.4.0. Nothing in the gate
  checks a number written in prose. Re-measure them whenever a release touches tests:
  `cyrius test tests/<suite>.tcyr | grep -oE '^[0-9]+ passed'` per suite, and the
  parenthesised-total filter from CLAUDE.md for the project sum.
- ⚠ **A recorded measurement is a claim like any other.** 2.4.7 found `_fresnel_fg`'s
  accuracy table reporting the PUBLISHED value in a column that reports what the code
  does — one cell, in an otherwise correct five-row table. A table of measurements earns
  trust that prose does not, so a wrong cell in one is worse than a vague sentence.
  Re-run the numbers a comment states when you touch the function.
- ⭐ **Sweep a defect CLASS mechanically, not by reading.** Three releases guarded null
  handles by inspection (2.2.7 encoders, 2.2.8 Spd, 2.4.8 Prescription) and each
  asserted the job was done. 2.4.9 called every public function (413) in its own
  process with handles, arrays and counts zeroed, then again with counts = 4 so that
  loop bodies were reached: **68 more SIGSEGVs**. The harness is a generator over
  `cyrius.cyml`'s module list plus one `cyrius build` per function; run builds in
  parallel only with a serial retry, because concurrent builds race on dependency
  resolution (cyrius refuses the racing write safely and the build just fails).
- ⚠ **A numerically-stable rewrite of a hot expression must be benchmarked same-binary,
  and the first stable form is not necessarily the right one.** 2.4.9's GGX repair:
  `ndh²·α² + (1−ndh)(1+ndh)` was exact and cost **+18%** (every `f64_*` op is a real
  call — Cyrius has no `#inline`); `ndh²·α² + (1−ndh²)` is equally exact and cost
  **nothing**, same op count as the original. Both were correct; only measurement
  told them apart.
- ⭐ **Pin a known UPSTREAM defect at its defective value, so the note retires itself.**
  bayan's decoder mis-rounds ~1 in 10⁵ doubles by 1 ULP; tests/hardening.tcyr
  asserts the WRONG result for 1.621274542797433e-9 on purpose. When bayan is fixed
  that assertion fails, which is the signal to delete it and the caveat in
  src/serialize.cyr. Same pattern the file already uses for the math-shim divergences.
  Filed upstream at 2.5.0 with a self-proving repro, 27 vectors and the root cause
  (bayan `docs/development/issues/2026-09-22-prakash-f64-parse-double-rounding-at-midpoint.md`).
- ⛔ **Two public functions that name the same state must be asserted against EACH
  OTHER, or they will contradict each other and both suites will stay green.** 2.4.8
  found `polarization_circular_right` (Jones, S3 = -1) and `stokes_circular_right`
  (Stokes, S3 = +1) disagreeing on the same `StokesVector` type. Each suite pinned its
  own side; 31 suites passed over a flat contradiction. When a second constructor for
  a named quantity appears, the cross-module assertion is part of the work.
- ⛔ **A documented null return is only safe if its CONSUMERS guard it — sweep by
  family, and do not take a previous sweep's word.** 2.2.8 guarded the `Spd` family
  and its own comment asserted null-guarding was "already the library-wide convention
  — ray_core, ray_trace, ray_system and the wave_* accessors all do it". It was false
  of `ray_system`: all seven Prescription consumers exited 139 on
  `prescription_from_json`'s documented 0. The remaining decoder families were then
  covered by 2.4.9's mechanical sweep (the ⭐ bullet above), not by reading.
- ⚠ **Sizing a numerical window from a measurement taken INSIDE that window proves
  nothing.** `_cri_cct_refine`'s ±3% bracket was justified by "McCamy's error, worst
  18.6 K at 2000 K, i.e. 0.93%" — measured only where the bracket already worked.
  McCamy is +50.1% out at 1200 K. A golden section then returns its endpoint silently,
  and Ra fell to 15.6 for a 1000 K blackbody that must score 100. Measure a tolerance
  across the DOMAIN it will be used on, including both edges.
- ⛔ **A shape check cannot see a multiplier, and a value pin derived from the code
  cannot see anything.** Four audits (2.2.6, 2.2.8, 2.4.6, 2.4.9) repaired the Seidel
  bracket against its best-form argmin and aplanatic zero — both invariant under any
  positive factor — while `spherical` sat 6x too large, `coma` 2x and LSA 16.7x too
  small, and the one value pin had been computed from the expression it tested.
  2.5.0 found all three by tracing a real singlet with `ray_trace` and comparing.
  **When a formula has an independent implementation in this library, pin against
  it**, and check that the residual shrinks where theory says it should (here: as
  the lens thins) — a wrong factor does not.
- ⚠ **A pin relating two functions catches them drifting apart, never both being
  wrong together.** 2.4.6's `LSA == h²·S₁/(2φ)` caught the second copy of the Seidel
  bracket, and then held both at the same wrong divisor until 2.5.0.
- ⛔ **When one quantity is computed in two places, a repair to one does not reach the
  other — and nothing will tell you.** 2.4.6 found `lens_longitudinal_spherical_aberration`
  still carrying BOTH halves of the defect 2.2.6 repaired in `lens_seidel_coefficients`
  thirty lines above it, and 2.2.8 audited that bracket again without looking down. The
  bracket was 6.15x too small at n = 1.5. **Before repairing a formula, grep for its other
  evaluation sites**, and leave behind an assertion that pins the copies against each
  other rather than trusting them to stay in step.
- ⭐ **An integral is the strongest check a distribution can be given.** NDFs
  (`integral D cos = 1` over the hemisphere), phase functions (`integral = 1` over the
  sphere), sampling PDFs and CMF tables all have one, and 2.4.5-2.4.6 found that none of
  them was pinned. A wrong 4pi, a dropped factor or a guard epsilon that has taken over
  the answer survives "is it positive", "is it symmetric" and a single spot value; it
  cannot survive the integral. The 2.2.6 GGX epsilon defect WAS this integral collapsing.
- ⭐ **A derived constant must be pinned to the RELATIONSHIP it comes from, not only to
  its output.** 2.4.5 found `_atm_prefactor` and `_atm_n_s` describing air at two
  different temperatures — each a plausible textbook number, together 12% wrong. The
  assertion that catches it states that `prefactor*3*N_S^2/(8pi^3)` reproduces
  `(n^2-1)^2` for the same air. ⚠ It is also the **only** pin that catches all three
  mutants: reverting `N_S` alone lands 2.9% from the published figure, inside the 3%
  band a physics assertion can justify. When two constants are derived from one
  another, pin the derivation.
- ⚠ **An order-of-magnitude range assertion is not a test of a physical quantity.**
  The atmosphere suite had 367 assertions and pinned Rayleigh scattering only as
  `beta > 1e-6 && beta < 1e-4` — two decades wide, so a 12% calibration error passed
  for the life of the module. Same shape as the too-wide tolerances the 2.4.2-2.4.4
  dispersion repairs found. A range check earns its place next to a value pin, not
  instead of one.
- ⭐ **Two independent models of one material is the sharpest instrument this library
  has, and it must be set to the model error, not to a round number.** prakash carries
  a Sellmeier fit AND a Schott series for N-BK7, and a Cauchy pair AND a Sellmeier fit
  for fused silica. 2.4.4 found real defects in both second members — and the Schott
  cross-check already existed, at `TOL_002`, **100x the 2.1e-5 defect it was built to
  catch**. Corrected, the pair agrees to 1.7e-6. When adding a preset that duplicates
  an existing model, pin the two against each other at their true agreement.
- **A fidelity repair below the instrument needs a pin on the LITERAL, or it reverts
  silently.** 2.4.4 corrected N-SF11's `b1` (2e-7 in n) and sapphire's three slipped
  digits (3.3e-5); the mutation check showed the `b1` revert passing all 79 assertions.
  Both are now pinned at the coefficient via the accessors at 1e-8, labelled in the
  test file as regression pins rather than physics checks.
- **Every Sellmeier preset carries an Abbe pin, and that is the assertion that bites.**
  Both 2.4.x coefficient defects (diamond, water) sat within ~0.4% on n_d while being
  8.9 and 9.4 out on V_d. A preset's index can look right while its dispersion is
  badly wrong; pin both.

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

Who uses what, how, and at which version is in the consumers section of
`docs/architecture/overview.md`. Until 2.5.1 this table said soorat used prakash's
Cook-Torrance (it hand-copies it into WGSL), kiran its lighting math (it re-exports
modules and calls none) and ranga its lens effects (ranga uses only spectral/colour),
and it omitted tanmatra.

## Boundary with Other Crates

| Feature | prakash | other |
|---------|---------|-------|
| Fresnel/Snell reflectance math | Yes (self-contained) | — |
| EM ↔ optics primitive bridge | Yes (`bridge`) | bijli (EM foundation) |
| Pixel-level image filter | — | ranga |
| 3D scene graph | — | kiran |
| Complex arithmetic (complex Fresnel for absorbing media) | Yes (own) | — |
| Vector/matrix math, FFT | — | hisab (prakash calls only `num_fft`) |
| Color space conversion (ICC) | — | ranga |
| Spectral → RGB conversion | Yes | — |
| Polarization formalism (Jones/Stokes/Mueller) | Yes | — |
