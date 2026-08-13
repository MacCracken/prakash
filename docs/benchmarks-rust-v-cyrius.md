# Benchmarks: Rust vs Cyrius

> prakash v2.0.1 benchmark comparison.
>
> - **Rust**: criterion v0.8, release mode. Latest run from `rust-old/bench-history.csv` as of tag 2.2.2 (2026-03-26, commit `01684d4`). f64 throughout.
> - **Cyrius**: cyrius 6.5.20, `tests/prakash.bcyr`. Run 2026-08-12 (commit `7ec51c8`). f64 (IEEE-754 bits in i64), heap-allocated structs.
> - **Platform**: x86_64 Linux.
> - **Pre-SIMD baseline**: Cyrius 6.5.20 emits scalar f64 only — no SIMD (it's on the Cyrius roadmap). Every number here is scalar-vs-scalar (or scalar Cyrius vs Rust's autovectorized/glam paths). The gap is a current-toolchain artifact, not a language ceiling; SIMD closes the vector/FFT-heavy paths.
> - The Rust numbers below are the record: `rust-old/` was removed in 2.2.3, so
>   the full Rust set is in git rather than in-tree —
>   `git show 2.2.2:rust-old/benchmarks.md` (and `rust-old/bench-history.csv` for
>   the raw data). The Cyrius history is in [`../benchmarks.md`](../benchmarks.md).

## What changed in this revision

The v2.0.0 edition of this page carried a large caveat: the Cyrius numbers were
**measurement-floor-bound**. The bench harness had a fixed per-iteration cost
(fnptr indirection + timer read) of ~1.32 µs, every scalar optics bench bottomed
out there, and the resulting ratios read **60×–1,300×**. Those ratios were an
artifact of the instrument, not a property of the code — as the old page said in
its own warning.

**Cyrius 6.5.20 fixed the instrument.** The harness now measures the timer floor
and subtracts it from every sample:

```
[timer floor 1.347us per clock read, measured; subtracted from every sample]
```

That was optimization vector #1 on the old page, and it is now resolved upstream.
The evidence that the floor is genuinely gone is the `min` column: it used to
quantize to ~489 / 838 / 907 ns regardless of what the function did. It no longer
does — `pbr/fresnel_schlick` now reports a 17 ns min and
`atmosphere/rayleigh_cross_section` an 8 ns min, tracking their actual work.

**The honest gap is 3×–20×, median ~7×** — not the 60×–1,300× the instrument
previously suggested. No prakash code got faster between 2.0.0 and 2.0.1; only the
measurement did. Do not read the deltas in [`../benchmarks.md`](../benchmarks.md)
across the 6.5.20 boundary as improvements.

## Head-to-Head (24 matched operations)

| Operation | Rust (ns) | Cyrius avg (ns) | Cyrius min (ns) | avg ÷ Rust |
|---|--:|--:|--:|--:|
| **Ray** | | | | |
| fresnel_unpolarized | 24.91 | 108 | 103 | 4.3× |
| sellmeier_n_at | 5.99 | 51 | 49 | 8.5× |
| fiber_na | 2.84 | 9 | 9 | 3.2× |
| v_number | 1.26 | 9 | 7 | 7.2× |
| **Spectral** | | | | |
| wavelength_to_rgb | 9.58 | 89 | 84 | 9.3× |
| planck_radiance | 11.22 | 80 | 77 | 7.1× |
| **Wave** | | | | |
| single_slit_intensity | 17.28 | 121 | 116 | 7.0× |
| malus_law | 7.83 | 46 | 43 | 5.9× |
| bessel_j1 | 10.78 | 61 | 58 | 5.7× |
| airy_pattern | 13.31 | 217 | 206 | 16.3× |
| fresnel_integral_c | 16.42 | 157 | 149 | 9.6× |
| coating_reflectance | 15.26 | 103 | 98 | 6.7× |
| zernike_poly | 18.07 | 110 | 104 | 6.1× |
| diffraction_pattern_2d ‡ | 5,639.30 | 15,888 | 14,996 | **2.8×** |
| **PBR** | | | | |
| fresnel_schlick | 1.04 | 17 | 17 | 16.3× |
| distribution_ggx | 1.71 | 22 | 21 | 12.9× |
| cook_torrance | 7.23 | 101 | 97 | 14.0× |
| henyey_greenstein | 3.70 | 19 | 18 | 5.1× |
| iridescence_fresnel | 23.19 | 144 | 138 | 6.2× |
| **Lens** | | | | |
| mtf_diffraction_limited | 7.68 | 88 | 84 | 11.5× |
| seidel_coefficients | 10.71 | 117 | 102 | 10.9× |
| **Atmosphere** | | | | |
| rayleigh_cross_section | 1.23 | 9 | 8 | 7.3× |
| mie_phase_cornette_shanks | 3.70 | 23 | 22 | 6.2× |
| sky_color_rgb | 51.55 | 1,014 | 967 | 19.7× |

Spread across the 23 scalar rows: **min 3.2×** (`fiber_na`), **max 19.7×**
(`sky_color_rgb`), **median 7.2×**.

‡ **diffraction_pattern_2d** is the compute-bound row (a real 2D FFT via hisab)
and the only one that was never floor-bound. Sizes differ (Cyrius 8×8 = 64 points
vs Rust `diffraction_2d_16` = 16 points), so it is not strictly like-for-like —
but it was **~3×** under the old harness and is **2.8×** under the new one, which
is the useful cross-check: the row the floor never touched did not move.

`spectral/spd_at` is omitted from the table — the Cyrius bench is a single
interpolated sample lookup while the nearest Rust bench (`spd_blackbody`, 992.5 ns)
builds a full 81-sample blackbody SPD. It was never a like-for-like pair.

## Cyrius-only benches (no Rust counterpart)

Added in the Cyrius era (`bridge` primitive hooks; bayan JSON `serialize`):

| Benchmark | Cyrius avg (ns) | Cyrius min (ns) | Iterations |
|---|--:|--:|--:|
| bridge/stellar_temperature_to_rgb | 119 | 115 | 100,000 |
| serialize/rgb_to_json | 3,443 | 3,145 | 50,000 |

⚠ **`serialize/rgb_to_json` is the one row that genuinely regressed.** Correcting
its v2.0.0 reading for the old harness floor (2,840 − ~1,320 ≈ **1,520 ns**), it
roughly **doubled** to ~3,400 ns. The cause is upstream, not in prakash: the
function's body is unchanged, and bayan 1.4.1 rebuilt the serializer so every
append routes through an allocator-variant with per-append failure propagation
(`_jb_append_string` → `_jb_append_string_a(default_alloc(), …)`), adding a call
layer and an allocator lookup per append. That is a deliberate
robustness-for-speed trade upstream, landing on prakash's one allocation-dominated
bench. Tracked on the roadmap; it does not affect any math path.

## Where the real gap comes from

With the harness floor removed, the remaining gap is attributable:

| Factor | Cost | Applies to |
|--------|------|-----------|
| f64, no SIMD | ~1.5–3× | the scalar math on every row |
| fn-call / no inlining | ~a few ns/call | wrappers and `_pbr_*`/`_lens_*` helpers — significant on the 1–20 ns ops |
| Heap structs | ~200–400 ns/alloc | only the struct-returning ops (Rgb, Pattern2D, CardinalPoints) |

The pattern in the table follows from this. The **highest** ratios are the
*cheapest* Rust ops — `fresnel_schlick` (1.04 ns) and `distribution_ggx`
(1.71 ns) — where a fixed few-ns call overhead is most of the Cyrius number;
16× of ~1 ns is still only 17 ns. The **lowest** ratios are the ops with enough
real arithmetic to amortize that overhead. The outlier at the top,
`sky_color_rgb` (19.7×), is the one scalar row that allocates: it returns RGB
through composite scattering terms, so it pays the heap-struct cost the others
avoid.

Note prakash's Rust was **already f64** (optics demands double precision), so
unlike some ports there is no f32→f64 penalty to attribute.

## Where Cyrius wins

| Metric | Rust | Cyrius |
|--------|------|--------|
| External crates | hisab + serde + thiserror + tracing (+ reqwest/tokio/serde_json/tracing-subscriber for `ai`/`logging`) | hisab (git) + stdlib — **0 external crates** |
| Library source | 17,904 lines / 27 files | **7,592 lines / 27 files** |
| Distribution | crate + transitive dep tree | one self-contained `dist/prakash.cyr` bundle |
| Build | cargo + dep compilation | single static binary, instant |
| Reproducibility | depends on dep versions | bit-identical (pinned toolchain + vendored `lib/`) |
| Precision | f64 | f64 (bit-faithful to the Rust original — constants match to the ULP) |

## Full Cyrius set (the 27 benchmarks as of 2.0.1, cyrius 6.5.20)

⚠ This page is a **2.0.1-era snapshot**, kept because it is the like-for-like
Rust comparison and the Rust side can never be re-run now that `rust-old/` is
out of tree. The Cyrius side has moved since: the suite is 36 benchmarks as of
2.2.0, and 2.1.1 removed per-call allocations that shifted several rows (e.g.
`fresnel_integral_c` 157 → 136 ns). For current Cyrius numbers see
[`../benchmarks.md`](../benchmarks.md); the ratios here remain the best available
estimate of the language gap.

| Benchmark | Avg (ns) | Min (ns) | Iterations |
|---|--:|--:|--:|
| atmosphere/mie_phase_cornette_shanks | 23 | 22 | 500,000 |
| atmosphere/rayleigh_cross_section | 9 | 8 | 500,000 |
| atmosphere/sky_color_rgb | 1,014 | 967 | 20,000 |
| bridge/stellar_temperature_to_rgb | 119 | 115 | 100,000 |
| lens/mtf_diffraction_limited | 88 | 84 | 200,000 |
| lens/seidel_coefficients | 117 | 102 | 200,000 |
| pbr/cook_torrance | 101 | 97 | 200,000 |
| pbr/distribution_ggx | 22 | 21 | 500,000 |
| pbr/fresnel_schlick | 17 | 17 | 500,000 |
| pbr/henyey_greenstein | 19 | 18 | 500,000 |
| pbr/iridescence_fresnel | 144 | 138 | 100,000 |
| ray/fiber_na | 9 | 9 | 500,000 |
| ray/fresnel_unpolarized | 108 | 103 | 500,000 |
| ray/sellmeier_n_at | 51 | 49 | 200,000 |
| ray/v_number | 9 | 7 | 500,000 |
| serialize/rgb_to_json | 3,443 | 3,145 | 50,000 |
| spectral/planck_radiance | 80 | 77 | 200,000 |
| spectral/spd_at | 43 | 41 | 200,000 |
| spectral/wavelength_to_rgb | 89 | 84 | 200,000 |
| wave/airy_pattern | 217 | 206 | 200,000 |
| wave/bessel_j1 | 61 | 58 | 500,000 |
| wave/coating_reflectance | 103 | 98 | 200,000 |
| wave/diffraction_pattern_2d_8x8 | 15,888 | 14,996 | 2,000 |
| wave/fresnel_integral_c | 157 | 149 | 500,000 |
| wave/malus_law | 46 | 43 | 500,000 |
| wave/single_slit_intensity | 121 | 116 | 200,000 |
| wave/zernike_poly | 110 | 104 | 200,000 |

## Optimization vectors

1. ~~**A lower-overhead bench harness**~~ — **done upstream in Cyrius 6.5.20**,
   which measures and subtracts the timer floor. This page is the first edition
   that reports true op cost.
2. **Inline expansion** of the tiny `f64_*` wrappers and `_pbr_*`/`_lens_*`
   helpers — now the *largest* remaining lever, since per-call overhead is the
   dominant term on the cheap ops (see `fresnel_schlick` at 16.3×).
3. **SIMD** (Cyrius roadmap) — would help the FFT and any batched pattern/PSF
   work, the only genuinely compute-bound paths.
4. **Arena allocation** for the struct-returning ops (Pattern2D grids, ray-fan
   buffers) — amortize the heap cost across a batch. `sky_color_rgb`, the worst
   scalar row at 19.7×, is the case in point.
5. **`serialize/rgb_to_json`** — re-check against bayan upstream; the `_a`
   allocator variants (`bayan_json_v_obj_set_a`, `bayan_json_v_build_a`) allow
   passing an arena instead of `default_alloc()` per append.
