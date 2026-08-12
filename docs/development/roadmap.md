# Prakash Roadmap

> **Prakash** is the optics/light-simulation library, written in [Cyrius](https://github.com/MacCracken/cyrius). Math foundations (Complex + FFT) come from [hisab](https://github.com/MacCracken/hisab). Consumed by soorat (PBR shading), kiran (lighting), and ranga (lens effects).

## Scope

Prakash owns the **physics of light**: how light travels, bends, reflects, interferes, diffracts, polarizes, and interacts with materials. It provides the math; consumers decide what to do with it (render pixels, simulate experiments, process images).

Prakash does NOT own:
- **Rendering pipeline** → kiran/soorat (they consume prakash for lighting math)
- **Image processing** → ranga (pixel operations, color spaces, filters)
- **Math primitives** → hisab (vectors, geometry, calculus, Complex, FFT)
- **Color science beyond spectral** → ranga (ICC profiles, gamut mapping)

## Completed

| Phase | Release | Summary |
|-------|---------|---------|
| V0.1–V0.29 | Rust foundation → V1.0 | Ray/wave/spectral/lens/pbr/atmosphere physics, optical systems, CIE color science, advanced PBR, pattern computation (see git history / `rust-old/`) |
| V1.0 | Stable (Rust) | API review, `#[must_use]`/`#[non_exhaustive]` audit, doc coverage, architecture docs |
| V1.1 | P1 audit (Rust) | Polarization consolidation, zero-alloc SPD (`Cow`), missing-docs, wavelength-first ordering |
| V1.2 | Bridge (Rust) | Removed `bijli-backend`; primitive-value cross-crate bridges (bijli/tara/badal) |
| **V2.0** | **Cyrius port** | **Full Rust→Cyrius rewrite. All 25 science modules ported and adversarially verified; two distlib bundles (math-only + AI); sakshi logging; bayan serialization; benchmark suite + CSV history. 5251 test assertions, 27 benchmarks. See the [2.0.0] CHANGELOG entry.** |
| V2.0.1 | Toolchain refresh | Cyrius 6.4.10 → 6.5.20, hisab 2.6.8 → 2.11.1, sakshi 2.4.2 → 2.4.10. `bayan_json_v_parse_str` → `_parse_buf` rename absorbed. 5251 assertions unchanged. 6.5.20's harness subtracts its timer floor, so the Rust-vs-Cyrius comparison is measurable for the first time (~7× median, was reported as 60×–1,300×). See the [2.0.1] CHANGELOG entry. |
| V2.0.2 | Audit & hardening | Six-dimension P(-1) sweep with adversarial verification. Fixed 7 public entry points that SIGSEGVed on malformed input, unchecked/overflowing allocations, missing bounds checks on Pattern2D and Mueller accessors, an out-of-bounds read in `spd_at`, silent fabrication on bad JSON, flat-surface JSON corruption, and three math boundaries that diverged from Rust. 5251 → 5334 assertions, 26 → 27 suites. See the [2.0.2] CHANGELOG entry. |
| **V2.1.0** | **Error channels** | **Every `*_from_json` (9 functions) takes a trailing `err_out` and reports `PK_ERR_PARSE` (bad bytes) / `PK_ERR_INVALID_PARAMETER` (wrong schema) / `PK_ERR_ALLOCATION`, closing the gap 2.0.2 left open. Breaking signature change; migration table in the [2.1.0] CHANGELOG entry. 5361 assertions.** |
| **V2.1.1** | **Audit cleanup** | **Clears every remaining 2.0.2 audit item: trig out-of-domain (`f64_sin(1e300)` returned 1e300), `expm1` precision loss below 1e-16, the unusable `point_source_new`, and 46 unguarded constructor allocs. Closes all three coverage gaps and adds 8 composite benchmarks (27 → 35) — `interference_pattern` is 7× the slowest previously measured op. 5424 assertions.** |
| **V2.1.2** | **hisab interop** | **`tests/hisab_interop.tcyr` pins the `RayVec3` ↔ `HVec3` layout contract (identical `{x;y;z}`), so all 26 `hvec3_*` ops work on prakash vectors unconverted — incl. `ray_reflect_3d` agreeing with `hvec3_reflect` over 121 cases. Internal `hvec3_dot` adoption measured at 2–11% slower and reverted; hisab quadrature investigated and does not apply. No behaviour change. 5444 assertions, 28 suites.** |

### V2.0 exit criteria (met)

- [x] Every Rust module translated to Cyrius, bit-faithful (exact-ratio / hex constants, `powi` square-and-multiply, left-assoc fold order)
- [x] Each module adversarially verified against `rust-old/`
- [x] `PK_ERR_*` error namespace (avoids Cyrius last-one-wins global collisions)
- [x] Two bundles: `dist/prakash.cyr` (core) and `dist/prakash-ai.cyr` (+ AI client)
- [x] Logging ported onto sakshi (`prakash_set_log_level` + 31 trace markers)
- [x] Serialization reimplemented on bayan JSON
- [x] Benchmark suite (`tests/prakash.bcyr`) + `scripts/bench-history.sh` baseline
- [x] Cyrius CI/release workflows; VERSION = 2.0.0 (SSOT via `${file:VERSION}`)

## Cross-Crate Bridges (done in V2.0)

The `bridge` module ships primitive-value hooks — no dependency on sibling crates:

- [x] **bijli bridge**: EM frequency/wavelength ↔ nm; E-field ↔ intensity; Cauchy refractive index (BK7)
- [x] **tara bridge**: stellar temperature → blackbody RGB; B-V color index → temperature; spectral class → temperature; absolute magnitude → luminosity
- [x] **badal bridge**: density → Rayleigh scale; humidity → Mie scale; cloud cover → diffuse fraction; visibility → extinction

## Post-2.0 Backlog

### Consumer integration (genuinely waiting on the consumers)

- [ ] soorat / kiran / ranga: consume `dist/prakash.cyr` directly once they move to Cyrius

### Dependency adoption (settled in 2.1.2 — both investigated, one shipped)

- [x] **hisab geometry bridge — shipped as a tested interop contract, not an
      internal swap.** `RayVec3` and `HVec3` are both `{ x; y; z; }`, so a prakash
      vector is already an hisab vector and all 26 `hvec3_*` ops work on it
      unconverted. `tests/hisab_interop.tcyr` now pins the layout, the field
      order, and that `ray_reflect_3d` agrees with `hvec3_reflect` (121 cases).
      ⚠ **prakash does NOT call `hvec3_*` internally, and that is measured, not
      assumed:** routing the tracer's dot products through `hvec3_dot` cost
      **2–11%** on `ray/trace_surface`, `ray/trace_sequential` and
      `ray/fresnel_unpolarized` across three runs, because it is two nested calls
      (`hvec3_dot` → `f64v_dot`) where the inline form is straight-line
      arithmetic. Reverted. Revisit only if hisab gains an inlinable form.
- [x] **hisab numerical integration — investigated, DOES NOT APPLY.** hisab does
      ship quadrature (`calc_integral_gauss5(f, a, b)`, `calc_integral_simpson`,
      `calc_adaptive_simpson`), but none of the named targets is a quadrature:
      - `_spd_integrate` is the **CIE-defined weighted sum** over the tabulated
        81-entry CMFs at 5 nm. It is the standard's method, not an approximation
        to improve; replacing it would change published colour values and break
        1168 `spectral_cie` assertions.
      - `huygens_fresnel_1d` integrates a **caller-supplied discrete sample
        buffer** — there is no continuous integrand to hand a quadrature rule.
      - `spd_blackbody` **samples** 81 points; it does not integrate.

      `src/` contains **zero `fncall` sites**, i.e. no prakash function takes the
      callable integrand these routines require. If a continuous-integrand entry
      point is ever wanted (e.g. `spd_from_function(f, start, end)`), that is a
      new feature, not dependency adoption — filed under Accuracy below.

### Performance (from the 2.0.1 measurement)

Now that the bench harness reports true op cost, these are ranked by measured
impact rather than guessed at — see `docs/benchmarks-rust-v-cyrius.md`.

- [ ] **`serialize/rgb_to_json` ~2×** (floor-corrected ~1,520 ns → ~3,400 ns).
      ⚠ **Not a regression to undo.** Measured decomposition: 85% of the cost is
      `bayan_json_v_build` and ~720 ns of that is *per float*, because bayan 1.2.1
      replaced a 6-decimal renderer with round-trip-correct Grisu2. That bought
      bit-exact f64 round-trips and fixed a real data-loss bug (the old encoder
      flushed any |x| < 5e-7 to zero). The remaining lever is the 14% spent on
      object construction — the `_a` allocator variants
      (`bayan_json_v_obj_set_a` / `_build_a`) can take an arena instead of
      `default_alloc()` per append. Do not "fix" the float path.
- [ ] **Inline expansion** of the tiny `f64_*` wrappers and `_pbr_*`/`_lens_*`
      helpers — per-call overhead is now the dominant term on the cheap ops
      (`pbr/fresnel_schlick` is 16.3× Rust at 17 ns absolute).
- [ ] **Arena allocation** for struct-returning ops — `atmosphere/sky_color_rgb`
      is the worst scalar row (19.7×) and the only one that allocates.
- [ ] **SIMD** — available now, NOT blocked: `lib/simd.cyr` ships 127 functions
      (`f64v2_*` / `f64v4_*`, incl. `fmadd`/`dot`); it simply is not in prakash's
      `[deps] stdlib` yet. The compute-bound paths are the ones worth vectorising,
      now that 2.1.1 measures them: `wave/interference_pattern_16x16` (108 µs),
      `ray/spot_diagram` (28 µs), `wave/diffraction_pattern_2d_8x8` (15 µs).

### Correctness (from the 2.0.2 audit — ALL FIXED in 2.1.1)

Each was found and verified during the 2.0.2 sweep and deliberately left out of
that release. They are real and reproduced; none is a guess.

- [x] **`point_source_new` is unusable by its only consumer.** `interference_pattern`
      indexes a contiguous `PointSource` array (`sources + i*sizeof(PointSource)`),
      but the only public constructor returns standalone allocations that cannot be
      indexed that way — `tests/wave_pattern.tcyr` works around it with a local
      `_mksrc` helper. Needs either `point_source_array_new(n)` + a setter, or an
      init-into-caller-slot form.
- [x] **`f64_sin`/`f64_cos` return the argument unchanged for |x| >= 2^63**
      (x87 `fsin` out-of-range), silently corrupting `interference_intensity` and
      `single_slit_intensity` at extreme phase. Either range-reduce mod TAU in a
      prakash helper before the trig call, or document an explicit domain limit.
- [x] **`planck_radiance` substitutes `exp(x)-1` for Rust's `exp_m1`**, losing the
      small-x branch below x ~ 1e-16. Known deviation since the port; a proper
      `expm1` (series for |x| < 0.5) would close it.
- [x] **~78 `alloc(sizeof(T))` constructor sites do not check for 0.** Only
      reachable under heap exhaustion (the caller-sized ones were fixed in 2.0.2),
      but the constructors return an unguarded pointer, so `PK_ERR_ALLOCATION`
      cannot be reported from them. Decide whether to guard them all or document
      the limit.

### Coverage (from the 2.0.2 audit — ALL CLOSED in 2.1.1)

- [x] **Benchmark the expensive composites.** `tests/prakash.bcyr` covers 27 ops but
      misses every heavy one: `spd_to_xyz`, `luminous_flux`, `color_rendering_index`,
      `spd_blackbody`, `trace_surface`, `trace_sequential`, `spot_diagram`,
      `interference_pattern`. The audit measured `interference_pattern` at ~37× the
      slowest currently-benchmarked op, so the suite is not covering the real cost.
- [x] **`trace_surface` has 5 error branches and 1 is tested** — the four
      geometry-miss cases the whole tracing stack depends on are unexercised.
- [x] **`mueller_set`/`mueller_get` had zero direct tests** before 2.0.2 added the
      bounds cases; a full 16-element set/get roundtrip is still missing.

### Accuracy & completeness (demand-gated)

- [ ] Gradient-index (GRIN) optics: curved ray paths through variable-n media
- [ ] Diffractive optical elements (DOE): phase gratings, holographic elements
- [ ] Vectorial diffraction (Richards-Wolf): high-NA focusing beyond scalar theory
- [ ] Hermite-Gaussian / Laguerre-Gaussian beam modes; M² beam quality
- [ ] Higher-order (5th-order Buchdahl) aberrations; wavefront coefficients (W_040, …) from Seidel sums
- [ ] Aberrated MTF from generalized pupil-function autocorrelation
- [ ] `spd_from_function(f, start_nm, end_nm)` — build an SPD from a continuous
      spectral function via hisab `calc_integral_gauss5`. The one place hisab's
      quadrature would genuinely fit; a new capability rather than an adoption.

### Advanced / demand-gated

- [ ] Fluorescence (Stokes shift, excitation/emission spectra)
- [ ] Non-linear optics (SHG, Kerr) — if joshua needs it
- [ ] Orbital angular momentum (Laguerre-Gaussian modes)
- [ ] Metamaterials / negative refractive index
- [ ] Age-dependent CIE observer (CIE 2006)

### Housekeeping

- [ ] **Port-completeness review against `rust-old/`** — confirm nothing from the
      Rust original was left on the table (functions, tests, benchmarks, doc
      content) before anything is deleted.
- [ ] Remove `rust-old/` — **gated on the review above**, not on a version number
      (it is also preserved in the pre-2.0 git tags).

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
