# Benchmarks: Rust vs Cyrius

> prakash v2.0.0 benchmark comparison.
>
> - **Rust**: criterion v0.8, release mode. Latest run from `rust-old/bench-history.csv` (2026-03-26, commit `01684d4`). f64 throughout.
> - **Cyrius**: cyrius 6.4.10, `tests/prakash.bcyr`. Run 2026-07-06 (commit `6977621`). f64 (IEEE-754 bits in i64), heap-allocated structs.
> - **Platform**: x86_64 Linux.
> - The full Rust set lives in [`rust-old/benchmarks.md`](../rust-old/benchmarks.md); the Cyrius history in [`../benchmarks.md`](../benchmarks.md).

## ⚠ Read this first — the Cyrius numbers are measurement-floor-bound

Unlike a heap-heavy library, prakash's optics are **scalar f64 math** — a Snell
refraction or a Fresnel term is a handful of `f64_mul`/`f64_add` with no
allocation. Those operations are *faster than the Cyrius bench harness can
measure*. The harness has a fixed per-iteration floor (fnptr indirection + timer
read) of roughly **0.5–0.9 µs**, and every scalar optics bench bottoms out there.

The proof is in the `min` column below: it **quantizes to ~489 / 838 / 907 ns**
regardless of what the function does.

- `pbr/fresnel_schlick` — Rust **1.04 ns** (one `powi(5)` + lerp) → Cyrius min **488 ns**
- `pbr/iridescence_fresnel` — Rust **23.2 ns** (Airy thin-film, ~20× the work) → Cyrius min **907 ns**

A 23× difference in real work collapses to a <2× difference in the Cyrius
measurement. The harness simply cannot resolve below its floor. **So the "avg ÷
Rust" ratios below (60×–1300×) are dominated by harness overhead, not by a real
compute gap.** Treat them as an upper bound with most of the mass being fixed cost,
not as the cost of the optics math.

The one bench where real work exceeds the floor — the 2D FFT diffraction pattern —
tells the honest story: **~3×**, not 100–1000×.

## Head-to-Head (25 matched operations)

| Operation | Rust (ns) | Cyrius avg (ns) | Cyrius min (ns) | avg ÷ Rust |
|---|--:|--:|--:|--:|
| **Ray** | | | | |
| fresnel_unpolarized | 24.9 | 1,442 | 907 | 58× |
| sellmeier_n_at | 5.99 | 1,414 | 838 | 236× |
| fiber_na | 2.84 | 1,329 | 838 | 467× |
| v_number | 1.26 | 1,332 | 489 | 1,059× |
| **Spectral** | | | | |
| wavelength_to_rgb | 9.58 | 1,455 | 907 | 152× |
| planck_radiance | 11.2 | 1,434 | 838 | 128× |
| spd_at † | 992.5 | 1,414 | 489 | 1× |
| **Wave** | | | | |
| single_slit_intensity | 17.3 | 1,444 | 838 | 84× |
| malus_law | 7.83 | 1,378 | 489 | 176× |
| bessel_j1 | 10.8 | 1,426 | 907 | 132× |
| airy_pattern | 13.3 | 1,560 | 908 | 117× |
| fresnel_integral_c | 16.4 | 1,471 | 838 | 90× |
| coating_reflectance | 15.3 | 1,431 | 907 | 94× |
| zernike_poly | 18.1 | 1,441 | 907 | 80× |
| diffraction_pattern_2d ‡ | 5,639 | 15,462 | 14,108 | **3×** |
| **PBR** | | | | |
| fresnel_schlick | 1.04 | 1,340 | 488 | 1,285× |
| distribution_ggx | 1.71 | 1,343 | 838 | 787× |
| cook_torrance | 7.23 | 1,419 | 907 | 196× |
| henyey_greenstein | 3.70 | 1,359 | 489 | 368× |
| iridescence_fresnel | 23.2 | 1,457 | 907 | 63× |
| **Lens** | | | | |
| mtf_diffraction_limited | 7.68 | 1,417 | 838 | 185× |
| seidel_coefficients | 10.7 | 1,470 | 907 | 137× |
| **Atmosphere** | | | | |
| rayleigh_cross_section | 1.23 | 1,329 | 838 | 1,078× |
| mie_phase_cornette_shanks | 3.70 | 1,346 | 489 | 363× |
| sky_color_rgb | 51.5 | 2,329 | 1,397 | 45× |

† **spd_at** is not a like-for-like pair — the Cyrius bench is a single interpolated
sample lookup, while the nearest Rust bench (`spd_blackbody`) builds a full 81-sample
blackbody SPD. Ignore its ratio.

‡ **diffraction_pattern_2d** is the only compute-bound row (a real 2D FFT via hisab).
Sizes differ (Cyrius 8×8 = 64 points vs Rust 16-point), but both do genuine FFT work
that dwarfs the harness floor — so its **~3×** is the honest ceiling on the real gap.

## Cyrius-only benches (no Rust counterpart)

Added in the Cyrius era (`bridge` primitive hooks; bayan JSON `serialize`):

| Benchmark | Cyrius avg (ns) | Cyrius min (ns) | Iterations |
|---|--:|--:|--:|
| bridge/stellar_temperature_to_rgb | 1,453 | 489 | 100,000 |
| serialize/rgb_to_json | 2,840 | 1,396 | 50,000 |

## Where the real gap comes from

For prakash specifically, the genuine (measurable) Cyrius overhead is small and
narrow, because the optics are scalar:

| Factor | Cost | Applies to |
|--------|------|-----------|
| Bench harness floor | ~0.5–0.9 µs/iter | **every scalar bench** (dominates all rows but the FFT) |
| f64, no SIMD | ~1.5–3× | the actual math, once you subtract the floor |
| fn-call / no inlining | ~a few ns | wrappers and helpers |
| Heap structs | ~200–400 ns/alloc | only the struct-returning ops (Rgb, Pattern2D, CardinalPoints) |

Note prakash's Rust was **already f64** (optics demands double precision), so unlike
some ports there is no f32→f64 penalty to attribute — the honest compute gap is the
FFT row's **~3×**, plus per-op fn-call overhead that the harness can't isolate.

## Where Cyrius wins

| Metric | Rust | Cyrius |
|--------|------|--------|
| External crates | hisab + serde + thiserror + tracing (+ reqwest/tokio/serde_json/tracing-subscriber for `ai`/`logging`) | hisab (git) + stdlib — **0 external crates** |
| Library source | 17,904 lines / 27 files | **7,592 lines / 27 files** |
| Distribution | crate + transitive dep tree | one self-contained `dist/prakash.cyr` bundle |
| Build | cargo + dep compilation | single static binary, instant |
| Reproducibility | depends on dep versions | bit-identical (pinned toolchain + vendored `lib/`) |
| Precision | f64 | f64 (bit-faithful to the Rust original — constants match to the ULP) |

## Full Cyrius set (27 benchmarks, 2026-07-06)

| Benchmark | Avg (ns) | Min (ns) | Iterations |
|---|--:|--:|--:|
| atmosphere/mie_phase_cornette_shanks | 1,346 | 489 | 500,000 |
| atmosphere/rayleigh_cross_section | 1,329 | 838 | 500,000 |
| atmosphere/sky_color_rgb | 2,329 | 1,397 | 20,000 |
| bridge/stellar_temperature_to_rgb | 1,453 | 489 | 100,000 |
| lens/mtf_diffraction_limited | 1,417 | 838 | 200,000 |
| lens/seidel_coefficients | 1,470 | 907 | 200,000 |
| pbr/cook_torrance | 1,419 | 907 | 200,000 |
| pbr/distribution_ggx | 1,343 | 838 | 500,000 |
| pbr/fresnel_schlick | 1,340 | 488 | 500,000 |
| pbr/henyey_greenstein | 1,359 | 489 | 500,000 |
| pbr/iridescence_fresnel | 1,457 | 907 | 100,000 |
| ray/fiber_na | 1,329 | 838 | 500,000 |
| ray/fresnel_unpolarized | 1,442 | 907 | 500,000 |
| ray/sellmeier_n_at | 1,414 | 838 | 200,000 |
| ray/v_number | 1,332 | 489 | 500,000 |
| serialize/rgb_to_json | 2,840 | 1,396 | 50,000 |
| spectral/planck_radiance | 1,434 | 838 | 200,000 |
| spectral/spd_at | 1,414 | 489 | 200,000 |
| spectral/wavelength_to_rgb | 1,455 | 907 | 200,000 |
| wave/airy_pattern | 1,560 | 908 | 200,000 |
| wave/bessel_j1 | 1,426 | 907 | 500,000 |
| wave/coating_reflectance | 1,431 | 907 | 200,000 |
| wave/diffraction_pattern_2d_8x8 | 15,462 | 14,108 | 2,000 |
| wave/fresnel_integral_c | 1,471 | 838 | 500,000 |
| wave/malus_law | 1,378 | 489 | 500,000 |
| wave/single_slit_intensity | 1,444 | 838 | 200,000 |
| wave/zernike_poly | 1,441 | 907 | 200,000 |

## Optimization vectors

1. **A lower-overhead bench harness** — the single biggest win for *measurement*
   would be a loop-based microbench that amortizes the fnptr+timer floor, so scalar
   optics ops report their true (tens-of-ns) cost instead of the ~0.5 µs floor.
2. **Inline expansion** of the tiny `f64_*` wrappers and `_pbr_*`/`_lens_*` helpers —
   removes per-call overhead on the hot scalar paths.
3. **SIMD** (Cyrius roadmap) — would help the FFT and any batched pattern/PSF work,
   the only genuinely compute-bound paths.
4. **Arena allocation** for the struct-returning ops (Pattern2D grids, ray-fan
   buffers) — amortize the heap cost across a batch.
