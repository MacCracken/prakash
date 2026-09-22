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

**13 items. ⭐ Item 1 is a filed defect backlog from the 2.4.8 audit and can ship in any
patch; the rest are demand-gated or waiting on a consumer.**

| # | Bucket | Item | State |
|---|---|---|---|
| 1 | 2.4.x patch | 20 findings from the 2.4.8 audit, filed with evidence, unverified | open — see below |
| 2–7 | 2.x | GRIN, DOE, Richards-Wolf, HG/LG beams, Buchdahl, aberrated MTF | demand-gated |
| 8–12 | 2.x | Fluorescence, non-linear, OAM, metamaterials, CIE 2006 observer | demand-gated |
| 13 | Blocked | soorat / kiran / ranga consume `dist/prakash.cyr` | waiting on the consumers |

⚠ **Items 2–12 are not a backlog anyone is working through.** Each is a subsystem,
listed so the scope boundary stays visible. Do not start one without a consumer
asking for it.

## How this file is organised

Items are grouped by the **release class their change implies**, not by the order
anyone intends to do them:

| Bucket | Means | Rule |
|---|---|---|
| **2.4.x — patch** | no public API moves | internals, perf, tests, docs, tooling, data fixes |
| **2.5.x — minor** | adds public API | new entry points, new capability |
| **2.x — demand-gated** | adds a subsystem | build when a consumer actually asks |
| **Blocked** | not prakash's move | waiting on something external |

⚠ **The bucket is a SemVer classification, not a queue.** Anything in 2.4.x can
ship in any order, in any patch release, in any combination — the bucket only
promises it will not force a minor bump. Reshuffling is free by construction.

**Nothing here depends on anything else here.** ⚠ **Keep rows SHORT.** A row says
what the work is and what would block it. Measurements belong in `CHANGELOG.md`;
a row that grows past ~8 lines has started duplicating the release history and
should be cut back.

## 2.4.x — patch: no public API moves

- [ ] **The 2.4.8 audit's 20 lower-severity findings.** A 23-agent adversarial sweep
  verified 8 findings and repaired 6 (see the [2.4.8] CHANGELOG entry). These 20 ranked
  below the verification cut and are recorded **unverified** — each carries the finder's
  evidence and none has been independently reproduced, so confirm before repairing.
  ⚠ They are NOT one bite. The clusters below are separable and CLAUDE.md's sizing rule
  says not to batch them.

  - **Rgb/Xyz accessor family dereferences the documented 0 from rgb_from_json** (`src/spectral_core.cyr`) —
    MEASURED, exit 139. `rgb_from_json` is documented to return 0 and tests/hardening.tcyr:270-272 pins three
    ways of getting it. Probe /tmp/pkprobe/chain2.cyr: `var c = rgb_from_json("{\"r\":\"x\"}", 9, &e);` -> c
    == 0 (printed), then `rgb_to_
  - **Medium / Polarization / Stokes / Sellmeier / Zernike consumers dereference documented 0 sentinels**
    (`src/ray_core.cyr`) — MEASURED, exit 139 on every one (probes /tmp/pkprobe/chain3-5.cyr and u1-u11.cyr).
    Three reachable chains, each starting from a 0 the library itself documents: (1) `medium_custom(0.5,
    "bad", &e)` returns 0 with PK_ERR_INVALID_INDEX per ray_
  - **ai.cyr's two JSON encoders were left out of the 2.2.7 encode-side null sweep** (`src/ai.cyr`) —
    MEASURED, exit 139. Probe /tmp/pkprobe/ai1.cyr: `var c = ai_daimon_config_from_json("not json", 8, &e);`
    -> c == 0 (printed; ai.cyr:130-131 documents that return), then `ai_daimon_config_to_json(c)` -> exit=139,
    because line 118 does `str_f
  - **The photopic V(λ) table has zero value pins — a 3.4% corruption passes the whole suite**
    (`tests/spectral_photometry.tcyr`) — The photopic table's entire coverage is two assertions:
    `_assert_near(_vl(vp, 35), F64_ONE, TOL6, "photopic V peaks 1.0 at 555nm")` (line 27) and an 81-iteration
    `_in_range(_vl(vp, pi), Z, F64_ONE)` loop (line 30). Not one of the 81 values
  - **planck_radiance has no absolute-value pin anywhere — c1 = 2hc² is unconstrained by the test suite**
    (`tests/spectral_core.tcyr`) — Every assertion in the `blackbody` group (tests/spectral_core.tcyr lines
    121-155) is either scale-invariant or sign-only: `planck positive`, `peak > shorter`, `peak > longer`,
    `hotter more radiance`, `T=0 gives 0`, `T<0 gives 0`, `extreme l
  - **The Planckian Ra = 100 invariant loop samples only 2000–4800 K — exactly the window where the CCT clamp
    above does not bite** (`tests/spectral_cie.tcyr`) — `for (var _bb = 2000; _bb <= 4800; _bb = _bb + 400)`
    asserts Ra(spd_blackbody(T)) = 100 to 1e-3. The invariant is exactly right and the tolerance is tight; the
    SAMPLE RANGE is the problem. Measured Ra at the untested points of the same inva
  - **The negative-zero wire-format pin tests +0.0, and the behaviour it documents is the opposite of what
    bayan actually does** (`tests/serialize.tcyr`) — Two measurements (/tmp/prkaudit/p2.cyr,
    /tmp/prkaudit/p3.cyr). (1) `f64_neg(f64_from(0))` returns bits 0x0, i.e. POSITIVE zero — so line 233
    stores +0.0 and the assertion at 235-237 pins nothing about a signed zero. (2) With the real bit pa
  - **spectral_band_similarity's doc table publishes pre-2.2.7 measurements for F2 and F11; both numbers are
    stale by ~1.1 and ~1.7 points, and the pin the test file claims exists was never written**
    (`src/spectral_cie.cyr`) — The table at lines 1938-1943 is headed "Measured against published Ra" and
    gives the CURRENT output of this metric as F2 98.289 and F11 92.001. MEASURED with the shipped CIE 15:2004
    tables: spectral_band_similarity(illuminant_f2()) = 97.196
  - **The paragraph above _cri_cct_refine still documents the superseded search (40 ternary iterations, ±8%
    bracket) that the code beneath it replaced** (`src/spectral_cie.cyr`) — Lines 3209-3211 read "40 ternary
    iterations over a +/-8% bracket close the seed error to well under 0.01 K". The code at 3218-3243 is a
    22-step golden-section over a ±3% bracket (`_cie_r(97, 100)` / `_cie_r(103, 100)`, `for (var it = 0; it
  - **The "floats survive a to_json/from_json cycle BIT-EXACTLY" claim is false; bayan's decimal->f64 parser
    loses 1 ULP on ~1 in 38,000 finite doubles** (`src/serialize.cyr`) — Measured (/tmp/prkaudit/p7.cyr):
    499,748 random finite doubles encoded with the module's own path and decoded with rgb_from_json -> 13
    values came back with different bits, all exactly +/-1 ULP (e.g. json=6.28282780197287e+197, in=0x6900cf5
  - **pbr_split_sum_scale_bias's accuracy table was measured against the pre-2.2.8 (energy-gaining) LUT and is
    wrong by up to 3.1x — including understating the bias error a consumer would size against**
    (`src/pbr_advanced.cyr`) — Comment claims (measured against pbr_integrate_brdf_lut over an 11x10 grid):
    `Karis: scale max err 0.5547, mean 0.1275 | bias max 0.0861, mean 0.0130`, and at line 368 `worst point r
    = 0.6, n.v = 0.1, where the integral is 1.158 and the fit
  - **src/error.cyr's math-shim tables describe a toolchain two minor versions old — 4 of the 5 divergences
    they document are closed, and the file ships verbatim in both dist/ bundles** (`src/error.cyr`) — Two
    present-tense claim blocks, both false on the pinned cyrius 6.6.6 / ganita 1.2.6 (measured with
    /tmp/pk_probe2.cyr, /tmp/pk_probe3.cyr, raw bits):

(a) Lines 155-166, `WHAT THIS DOES NOT FIX — the three remaining divergent rows`:
    cl
  - **docs/guides/allocation.md's per-call byte table understates medium_to_json by 12% and spd_to_json by up
    to 2.2x — the two numbers a consumer sizes JSON output from** (`docs/guides/allocation.md`) — Measured
    with alloc_used() deltas, caches warm, pre-built handles (/tmp/pk_probe8.cyr, /tmp/pk_probe7.cyr) — same
    method the table states:
  rgb_to_json               doc 136   measured 136  (control: the method reproduces)
  medium_to_json
  - **The independent textbook check on the Seidel bracket quotes 1.0162 for n = 1.7; 2(n^2-1)/(n+2) is 1.0216
    and the shipped code measures 1.022** (`src/lens.cyr`) — src/lens.cyr:349-352 presents the best-form
    singlet argmin as the INDEPENDENT check that replaced 2.2.6's self-referential one: 'the textbook
    best-form singlet at p = -1, argmin q = 2(n^2-1)/(n+2): 0.7143 at n = 1.5, 0.8667 at 1.6, 1.0162 a
  - **wavelength_to_rgb returns a null handle with err_out = PK_ERR_NONE when rgb_new's allocation fails**
    (`src/spectral_core.cyr`) — `wavelength_to_rgb` stores `PK_ERR_NONE` into `err_out` at entry (line 112)
    and its success path ends `return rgb_new(f64_mul(r, factor), ...)` (line 163). `rgb_new` returns 0 when
    `alloc(sizeof(Rgb))` fails (line 45: `if (c == 0) { return
  - **_cri_context builds a full D-series reference SPD then discards it for every source below 5000 K — 680
    of the 2240 bytes a CRI call burns on a never-freeing allocator** (`src/spectral_cie.cyr`) — `var ref_spd
    = _cri_d_illuminant(t_cct); if (t_cct < 5000) { ref_spd = spd_blackbody(t_cct); }` — the D-series
    reconstruction (81 S0/S1/S2 basis mixes plus an allocation) runs unconditionally and is overwritten for
    every warm source, which
  - **A NaN surface radius is silently laundered into a valid flat surface (+inf) by the prescription round
    trip** (`src/serialize.cyr`) — Measured (/tmp/prkaudit/p5.cyr): a PrescriptionSurface with radius =
    0x7FF8000000000000 (NaN) encodes to {"name":"nanr","surfaces":[{"radius":null,...}],...} and decodes back
    to radius bits 0x7ff0000000000000 — +infinity, i.e. a legitimate
  - **A zero-length Spd encodes successfully but its own decoder rejects the result** (`src/serialize.cyr`) —
    Measured (/tmp/prkaudit/p10.cyr): spd_to_json on an Spd with len = 0 emits
    {"start_nm":380.0,"step_nm":5.0,"values":[]}, and feeding those exact bytes to spd_from_json returns 0
    with err = -7 (PK_ERR_INVALID_PARAMETER), because of the `if (
  - **pbr_distribution_ggx_aniso keeps the denominator epsilon that was removed from the isotropic GGX, so the
    two disagree below roughness 1e-3** (`src/pbr_advanced.cyr`) — At ax = ay and h_dot_x = h_dot_y = 0 the
    anisotropic GGX reduces analytically to the isotropic one exactly (the algebra cancels to
    a2/(pi*(ndh^2(a2-1)+1)^2)), so any difference is implementation, not physics. Measured ratio aniso/iso at
    n_d
  - **Two call-site counts in shipped comments are stale: '_prk_pow ... other seven call sites' is 10, and
    '_pbr_eps15 ... called from eight sites' is 13** (`src/error.cyr`) — src/error.cyr:152 ('The other seven
    `_prk_pow` call sites pass a finite constant exponent and cannot reach here'): `grep -rn '_prk_pow('
    src/*.cyr` excluding the definition and comment lines gives 12 call sites — atmosphere.cyr:157, bridge.

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
  `prescription_from_json`'s documented 0. **The remaining decoder families (Rgb/Xyz,
  Medium, Polarization, Stokes, Sellmeier, Zernike, and ai.cyr's encoders) are filed
  above and have not been checked.**
- ⚠ **Sizing a numerical window from a measurement taken INSIDE that window proves
  nothing.** `_cri_cct_refine`'s ±3% bracket was justified by "McCamy's error, worst
  18.6 K at 2000 K, i.e. 0.93%" — measured only where the bracket already worked.
  McCamy is +50.1% out at 1200 K. A golden section then returns its endpoint silently,
  and Ra fell to 15.6 for a 1000 K blackbody that must score 100. Measure a tolerance
  across the DOMAIN it will be used on, including both edges.
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
