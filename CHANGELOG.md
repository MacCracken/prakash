# Changelog

## [2.2.4] - 2026-08-21 — Toolchain bump: fifteen cyrius releases, and a CI block that could not survive them

Cyrius **6.5.20 → 6.5.33**, hisab **2.11.1 → 2.11.2**, sakshi **2.4.10 → 2.4.11**,
and — vendored with the toolchain — sandhi **1.9.9 → 1.9.10** and ganita
**1.0.4 → 1.1.4**. No feature work.

Suite **5485 assertions across 29 suites, unchanged on both sides of the bump**,
and that is the point: **not one finding here was surfaced by a failing test.**
Every one is a claim that was measured once, written down, and then quietly
stopped being true — the same shape hisab recorded in its own 2.11.2, which
absorbed this identical toolchain range.

### Security

- **CI installed an unverified binary toolchain, and the step could not fail.**
  Both workflows fetched the release tarball with a bare `curl -sLO`, untarred
  it, and ended every copy in `2>/dev/null || true`. This repo ships a
  `SECURITY.md` and a security job while doing that — on the **release**
  workflow too, which builds the artifact that ships. Both now pipe the pin into
  the **upstream installer**, fetched from the pin's own immutable tag, which
  fails closed on a missing, unverifiable or mismatched `.sha256` (CVE-21) and
  verifies an **Ed25519 signature** over `SHA256SUMS` where one is published
  (CVE-13).
- **sandhi 1.9.10** (vendored) fixes a P1 in which `_sandhi_resp_frame_a`
  returned `SANDHI_OK` carrying a **NULL body pointer** when the allocator could
  not fit the body copy, so a caller that correctly checked
  `sandhi_http_err_kind(r) == SANDHI_OK` then read from address 0.
  ⚠ **prakash is on that code path** — traced, not assumed:
  `sandhi_http_post_opts` (used at `src/ai.cyr:209`) delegates to
  `sandhi_http_post_opts_a(default_alloc(), ...)`, so the arena-flavoured
  framing runs, backed by the **global** allocator. And `src/ai.cyr:211` does
  exactly what the bug defeated — checks `err_kind` and then trusts the body
  pointer. What prakash was **not** in is the regime that made it reachable in
  practice upstream: a *constrained* arena where `rbuf`'s eager allocation at
  the full `max_response_bytes` cap starves the body copy. prakash takes the
  global allocator and sandhi's default 262,144-byte cap, so reaching it needs a
  genuine global allocation failure. 2.2.3 shipped `lib/sandhi.cyr` at 1.9.9,
  **without** the guard; 2.2.4 ships it with.

### Fixed

- **CI could not have installed 6.5.33 at all — reproduced, `rc=1`.** The
  hand-rolled block laid the toolchain out **flat** as `~/.cyrius/{bin,lib}`,
  but cyrius resolves its stdlib snapshot from `~/.cyrius/versions/<pin>/lib`,
  so `cyrius deps` fails with *"pins version 6.5.33 but it is not installed
  at .../versions/6.5.33/lib"*. A new **`Verify toolchain matches the pin`**
  step now asserts this where it is diagnosable — the original failure surfaced
  two steps downstream as a path error.
- ⭐ **CI's own remediation advice would have truncated a contributor's source
  file to zero bytes.** On format drift the workflow printed
  `run 'cyrius fmt $f > tmp && mv tmp $f'`. `cyrius fmt` **rewrites in place as
  of cyrius 6.5.28** (it was stdout-only before), so that redirect captures
  **0 bytes** — measured — and moves an empty file over the source. The advice
  now reads `cyrius fmt $f`, and `CONTRIBUTING.md` carries the same warning.
  The instruction had been wrong since 6.5.28 and nothing executes it in CI, so
  only a human following it would have found out.
- **`fmt_float_buf` dropped the carry when a fraction rounded up to a full
  unit — and corrupted the fraction while doing it.** Measured either side, same
  source:

  | value | decimals | 6.5.20 | 6.5.33 |
  |---|---|---|---|
  | `wien_peak(5778)` = 501.5182 nm | 0 | **501.1** | 502.0 |
  | 2.96 | 1 | **2.10** | 3.0 |
  | 0.999 | 2 | **0.100** | 1.00 |
  | 9.99 | 1 | **9.10** | 10.0 |

  ⚠ **Blast radius is `examples/` only, and that is verified, not asserted:**
  `fmt_float_buf` appears in all four examples and in **no** file under `src/`
  or `tests/`. The library surface renders floats through bayan's Grisu2, whose
  6.5.20 → 6.5.33 diff has no hunk in the encoder — **no JSON byte moves**.
  `examples/basic_optics.cyr` now prints `Sun peaks at 502.0nm` where it printed
  `501.1nm`; the true value is 501.5182 nm, so the new output is the correct
  round and the old one was simply wrong.
- **The three `bayan.cyr` "assigning non-pointer to typed pointer" warnings are
  gone** — every prakash compile emitted them; the build is now warning-clean
  apart from the pre-existing `large static data` notes.

### Changed

- **`lib/` re-vendored from the 6.5.33 snapshot** — 30 files changed and
  **`async_macos.cyr` added** (new upstream; a partial refresh would have left
  `async.cyr` with a dangling include). Verified **byte-identical to
  `~/.cyrius/versions/6.5.33/lib`, file by file**, rather than assumed.
  ⚠ prakash's `lib/` was checked against **its own declared pin first** and had
  **zero drift** — hisab found its tree had arrived mid-sync and that comparing
  old-pin against new-pin would have shown a tidy diff and missed it. prakash
  was clean; the check is recorded so the next bump repeats it.
- **`ganita_f64_pow`'s domain moved, and `src/error.cyr`'s measured table was
  stamped "Measured on cyrius 6.5.20" with a row that had gone false.** ganita
  1.1.4 special-cases zero base, zero exponent and **negative base with an
  integral exponent** ahead of the `exp(n*ln(base))` path. Re-measured on this
  host, raw f64 bits:

  | expression | 6.5.20 | 6.5.33 |
  |---|---|---|
  | `pow(0, 2)` | NaN | `0x0` (0.0) |
  | `pow(0, 0)` | NaN | `0x3FF0000000000000` (1.0) |
  | `pow(0, -1)` | NaN | `0x7FF0000000000000` (+inf) |
  | `pow(-2, 3)` | NaN | `0xC01FFFFFFFFFFFFE` (−7.99999999999999982) |
  | `pow(2, 3)` | `0x401FFFFFFFFFFFFE` | `0x401FFFFFFFFFFFFE` *(control)* |
  | `pow(-2, 0.5)` | NaN | NaN *(control)* |
  | `pow(1, 0)` | 1.0 | 1.0 *(control)* |

  **Three controls hold on both sides**, so the flips are real and not an
  artefact of the probe. **Nothing prakash computes moved**: `_prk_pow` already
  produced those three zero-base values, which is precisely why 5485 assertions
  passed either side and the prose could rot unobserved.
  - **`_prk_pow` is KEPT, and its rationale rewritten** — the stdlib now
    produces the same bits, so the branch is a deliberate pin against a silent
    upstream revert, not a repair.
  - ⚠ **A negative base no longer fails loudly.** With an integral exponent it
    is now finite but inexact (`pow(-2,3)` vs Rust's exact `-8`), where it used
    to be a NaN somebody would notice. No prakash call site passes a negative
    base — verified across all nine `_prk_pow` / `ganita_f64_pow` sites — but a
    future one would now break bit-fidelity silently.
  - The same stale reasoning in `tests/hardening.tcyr` is corrected; its
    `(was NaN)` assertion messages are kept as history, now meaning *"before
    6.5.33"* rather than *"without the shim"*.
- **52 of 61 gated files reformatted.** ⚠ **This drift is PRE-EXISTING and is
  not caused by the bump** — measured with each toolchain's own wrapper against
  a pristine checkout: 6.5.20 flags **the same 52 files**, and `diff` of the two
  sorted lists is empty. CI's Format check was already failing. Proven
  **leading-whitespace-only**: 555 insertions against 555 deletions,
  `git diff -w` **empty**, and no file changed line count.
- **`bench-history.csv` gains derived `regime` and `floor_ns` columns, and the
  trend table now filters on them.** `scripts/bench-history.sh`'s header claimed
  the harness's per-call overhead "is included ... the min in the raw output is
  the better absolute proxy" — false since cyrius 6.5.19, which taught the
  harness to calibrate one clock read and subtract it. The script had been
  echoing the harness's own `[timer floor ...]` banner while its header denied
  it. `benchmarks.md` was consequently anchoring every Δ to a pre-floor
  baseline and printing artefacts like −98.3% as improvements.
  ⚠ **`regime` is DERIVED** — read from whether the harness printed that banner
  — so it cannot go stale the way the constant it replaces did. ⚠ **What it
  cannot do:** it separates "floor subtracted" from "floor not subtracted"; a
  future instrument change that kept printing the banner would not flip it.
  This run reports the **81 rows across 3 earlier runs** it excluded rather than
  quietly dropping them.
- **Stale version prose corrected** where it had become false rather than
  historical: `cyrius.cyml`'s `[deps.hisab]` block (which named 2.11.1, sakshi
  2.4.10 and "the 6.5.20 stdlib" in one sentence), `README.md`'s consumer
  snippet (which pinned cyrius 6.5.20 and **prakash 2.0.1, four releases
  behind**), `docs/architecture/overview.md`, and `tests/serialize.tcyr`'s
  header (which justified its ~1e-5 tolerance with "bayan renders floats to 6
  decimals" — untrue since bayan 1.2.1 replaced that renderer with round-trip-
  correct Grisu2; the tolerance is now loose, not forced).
  Historical narrative — prior CHANGELOG entries, the roadmap's per-release
  rows, `docs/benchmarks-rust-v-cyrius.md`'s 2.0.1-era numbers — was
  **deliberately left alone**: a figure stamped "as of X" is a record, not a
  claim. Where such a doc carried a *toolchain* claim that is now false
  ("Cyrius emits scalar f64 only — no SIMD"), a correction was added beside it
  rather than the record being rewritten.
- **The `.deps` sidecar behaviour was re-measured, and half of what prakash
  documented had changed.** The core bundle's **over**-reporting reproduces
  unchanged at 6.5.33. The ai bundle's **under**-reporting does not: through
  6.5.20 the pruned inference yielded literally `syscalls io`, and it now yields
  ten folds — still short of the stack sandhi needs, but a different shortfall.
  `scripts/sync-deps-sidecar.sh` and `overview.md` both corrected.

### Performance

- ⭐ **Binary 12,767,544 → 2,785,712 bytes (−78%)** and static data
  **10,792,936 → 798,472 bytes (−93%)**, from upstream restructuring of the
  large stdlib folds. Same source, same build command.
- ⚠ **No optics performance claim is made in 2.2.4.** 36 benchmarks compared
  against a pre-bump baseline captured on this host in the same regime: **median
  +0.9%**, and the single row moving ≥10% is `ray/fiber_na` at 9 ns → 8 ns —
  one nanosecond of quantisation at the harness resolution limit, not a signal.

### Notes

- **A known latent divergence is now documented rather than fixed.**
  `_prk_pow(0.0, NaN)` returns 1.0 (a NaN exponent is neither `> 0` nor `< 0`,
  so it falls through to the `F64_ONE` branch); the raw builtin returns 0.0;
  IEEE-754 / C99 `pow`, which Rust follows, says NaN. All three disagree. It is
  **pre-existing, not introduced here**, and unreachable from prakash today —
  every call site passes an exponent it computed itself. Fixing it is a
  behaviour change and belongs in its own release, not in a toolchain bump.
- `cyrius distlib --check` reports `dist/prakash-ai.cyr` STALE for prakash and
  always will: it regenerates the sidecar (10 folds) and compares against the
  28-fold override prakash deliberately owns. Bundle **content** is idempotent —
  verified by hashing across two consecutive generations. CI does not use
  `--check`; it regenerates, syncs, and diffs, which is correct.

## [2.2.3] - 2026-08-12 — Remove the `rust-old/` translation archive

The Rust original is out of tree. It was retained through 2.2.2 as the
translation reference; the 2.2.1 port-completeness review confirmed nothing was
left behind (all 282 `pub fn`, 37 structs, 5 enums accounted for; the 4 dropped
examples ported), so the gate it was waiting on is satisfied.

No code, test, benchmark or API change — the suite is **5485 assertions across 29
suites**, unchanged, and `rust-old/` was never a build input (no reference in
`cyrius.cyml`, CI or `scripts/`).

### Removed
- **`rust-old/`** — 42 files, 820 KB.

  **It is not lost.** Recover the archive as it stood with
  `git checkout 2.2.2 -- rust-old/`, or take the pre-port project from tag
  `1.2.0`, where the same 33 `.rs` files sit at the repo root. Individual files:
  `git show 2.2.2:rust-old/src/lens.rs`, `git show 1.2.0:src/lens.rs`.

### Changed
- **Four claims that would have become false are corrected**, rather than left to
  rot: `README.md` and `docs/architecture/overview.md` both said the Rust source
  was "retained in-tree"; `CONTRIBUTING.md` told contributors to "compare against
  `rust-old/`" when in doubt about bit-fidelity; and
  `docs/benchmarks-rust-v-cyrius.md` carried a markdown link to
  `../rust-old/benchmarks.md` that would have 404'd. Each now names the git
  incantation that retrieves what it refers to.
- **`docs/benchmarks-rust-v-cyrius.md` is marked as the 2.0.1-era snapshot it is.**
  It is kept because it is the only like-for-like Rust comparison and the Rust
  side can never be re-run now — but the Cyrius side has moved (36 benchmarks
  since 2.2.0; 2.1.1's allocation removals shifted rows such as
  `fresnel_integral_c`, 157 → 136 ns). Current numbers live in `benchmarks.md`.
- `.gitignore` drops the now-dead `rust-old/target/` and `rust-old/Cargo.lock`
  entries.

### Notes
- The per-file provenance comments across `src/` and `tests/` ("ported from
  `rust-old/src/ray/mod.rs`") are deliberately **kept**. They record where each
  module came from, which is still true and still useful; the paths resolve
  against tag 2.2.2.
- Bit-fidelity to the Rust original remains a stated property of the port. What
  changes is where the reference lives, not whether it exists — and the
  assertions that encode it (exact-ratio and hex constants, `powi`
  square-and-multiply, left-associative fold order) are all in the suite.

## [2.2.2] - 2026-08-12 — P(-1) audit: the crash every dimension found

A full audit/hardening/security sweep of the 2.2.1 tree across six dimensions,
every candidate handed to an independent skeptic instructed to refute it.
**14 findings confirmed, 4 refuted.**

Suite **5467 → 5485** assertions.

### Security
- **`pattern2d_max_intensity` SIGSEGVs on a null pattern — found independently by
  ALL SIX audit dimensions and reproduced as exit 139.** It was the only
  `Pattern2D` entry point without a null guard: `pattern2d_get`, `pattern2d_set`
  and `pattern2d_normalized` all got one in 2.0.2 and this one was missed. Null
  is not hypothetical here — it is the documented return of `pattern2d_new`,
  `diffraction_pattern_2d`, `diffraction_pattern_circular`, `psf_from_wavefront`,
  `psf_diffraction_limited` and `interference_pattern`, and `tests/hardening.tcyr`
  already asserted those return 0. So the suite pinned the null and then the next
  natural call killed the process. Now returns 0.0, which is the fold's own
  identity.
- **`_zern_single` stored through an unchecked `alloc`** — the one allocation site
  both the 2.0.2 and 2.1.1 sweeps missed. It backs every named aberration
  (`zernike_defocus`, `zernike_spherical`, `zernike_coma_x/y`, `zernike_astigmatism_*`).

### Fixed
- **Allocation failures that returned 0 without reporting it.** Four `err_out`-
  carrying functions returned the failure sentinel while leaving `err_out` at
  `PK_ERR_NONE`, so a caller following the documented `if (prakash_is_err(e))`
  pattern proceeded to dereference null: `trace_surface` (whose own in-library
  callers then deref it), `find_system_properties`, and the five
  `*_from_json` decoders on their terminal constructor. All now write
  `PK_ERR_ALLOCATION`.
- **`prescription_from_json` silently truncated a lens.** It discarded
  `prescription_add_surface`'s result, so a failed surface append produced a
  Prescription with **fewer surfaces than the document** and `PK_ERR_NONE` — a
  different lens, reported as success.
- **`rgb_to_u8` returned `INT64_MIN` for a NaN channel** where Rust's `as u8` is a
  saturating cast returning 0. `f64_clamp` passes NaN through (correctly matching
  `f64::clamp`), so the NaN had to be caught at the cast. New `_r_sat_u8` also
  saturates the high end.
- **`_ai_json_str_field` reported `PK_ERR_INVALID_PARAMETER` for a non-object
  document**, contradicting its own module contract and the sibling `serialize`
  family, which report `PK_ERR_PARSE` for bad shape. Bad shape is a parse
  failure; a missing field is a parameter failure.

### Documentation
- **`math.md`'s clearcoat formula double-counted the geometry term.** It read
  `f_coat = D F G / (4(n·v)(n·l))`, but `pbr_clearcoat_geometry` **is**
  `1/(4(n·v)(n·l))` — the implementation multiplies `D·F·G` and is right.
- **README's benchmark count** said 27; it is 36.
- **Corrected a 2.2.1 claim**: that entry said tests were added for "the 3 ported
  functions that had no test reference". The 2.2.1 sweep compared against
  `rust-old/` and so never examined Cyrius-only additions — 12 public functions
  still had no test, bench or example reference. Ten now do (see below); the
  wording is corrected in place.

### Tests
- Coverage for the previously unreferenced public functions: `spec_visible_min_nm`,
  `spec_visible_max_nm`, `photometry_k_m`, `observer_cmf` (distinct table per
  observer), `cie_2015_10deg_table`, `xyz_d50_white`, `stokes_intensity`,
  `stokes_diagonal_minus`, `spd_from_static` (asserts it *borrows* the caller's
  buffer rather than copying), and `prakash_set_log_level`.

### Refuted
Four claims did not survive verification and were dropped rather than acted on:
`ai_register_agent`'s unchecked sandhi handles (the failure modes are not
reachable as claimed), the README Quick Start manifest (built and run — it
works), and twice, `point_source_set`'s missing upper bound — the array is a
length-less raw buffer by the same contract as `spectrum_strip`, and
`interference_pattern` takes the count separately.

## [2.2.1] - 2026-08-12 — Port-completeness review against `rust-old/`

A systematic comparison of the Rust original against the Cyrius port, run
**before** anything gets deleted. The function-level port is complete; two real
gaps were found and closed, and one is reported rather than closed.

Suite **5457 → 5467** assertions; **4 examples ported** (a new `examples/` tree).

### Verified complete
- **All 282 Rust `pub fn` are accounted for.** 274 matched by name or the port's
  prefix convention; the remaining 8 are word-order renames that resolve:
  `classify_lens`→`lens_classify`, `radial_polynomial`→`zernike_radial`,
  `thin_lens_image_distance`→`lens_thin_image_distance`,
  `thick_lens_focal_length`→`lens_thick_focal_length`,
  `separated_lenses_focal_length`/`_bfd`→`lens_separated_*`,
  `cmf_table`→`observer_cmf`, and `init_with_level`→`prakash_set_log_level`
  (the `tracing_subscriber` setup has no Cyrius analogue; sakshi replaces it).
- **All 37 `pub struct` and 5 `pub enum` are present.**

### Added
- **`examples/` — the 4 Rust examples were silently dropped by the port.**
  Nothing referenced them: not the port plan, not CI, not the docs. Ported and
  now run by CI:
  - `basic_optics.cyr` — Snell, Fresnel, wavelength→RGB, Wien, thin lens
  - `rainbow.cyr` — dispersion through a raindrop; 9000-step minimum-deviation
    sweep per wavelength. Reproduces the physics: violet 40.29°, red 42.11°
  - `pbr_materials.cyr` — Fresnel-Schlick, Cook-Torrance dielectric vs metal,
    clearcoat, thin-film iridescence
  - `camera_lens.cyr` — a cemented doublet through the whole ray stack:
    paraxial properties (EFL 49.07 mm), marginal trace, sequential trace with
    per-surface reflectance, spot diagrams on/off-axis, OPD fan

  They `include "dist/prakash.cyr"` — the exact surface a consumer gets — so CI
  running them is an **end-to-end check of the published bundle**, not just of
  `src/`. Lint and fmt now cover `examples/` too.
- **Tests for the 3 RUST-ORIGIN functions that had no test reference**
  (the sweep compared against `rust-old/`; Cyrius-only additions were out of its
  scope and 2.2.2 closes those separately):
  `spd_end_nm`, `fresnel_integral_cs` (agreement with the separate `C`/`S`
  entry points across both branches and for negative x), and
  `pbr_clearcoat_distribution` (GGX lobe peaks at `n_dot_h = 1`, tighter when
  smoother).

### Notes
- ⚠ **An example silently redefined a library function.** A helper named `_r` in
  `rainbow.cyr` collided with `spectral_core`'s own `_r` inside the bundle —
  Cyrius globals are **last-one-wins with no diagnostic**, so the example
  overrode the library's helper. Harmless only because both happened to be
  identical. All example helpers are now prefixed (`_rb_`, `_pm_`, `_cl_`,
  `_ex_`) and each file says why.
- ⚠ **Benchmark coverage is the one gap left open**: the Rust suite had 180
  benchmarks against prakash's 36, and 131 Rust bench subjects have no Cyrius
  counterpart. Most are trivial scalar micro-benchmarks (`abbe_number`,
  `brewster_angle`, `critical_angle`, …) and 2.1.1 already added the expensive
  composites that actually dominate runtime. Filed on the roadmap rather than
  bulk-ported, because 131 more one-line scalar benches would add noise to
  `bench-history.csv` without changing a decision.

## [2.2.0] - 2026-08-12 — SIMD, and where it actually applies

Adds the `simd` stdlib fold and applies it to the one loop in prakash that
genuinely benefits. The larger result is negative and is recorded as such: most
of prakash's cost is transcendental and branchy, not elementwise-array, and the
stdlib has no vector transcendentals. Two of the three benchmark hot spots
cannot be vectorised at all with the available primitives.

Suite **5444 → 5457** assertions across **28 → 29** suites; benchmarks 35 → 36.

### Added
- **`simd` in `[deps] stdlib`**, and `pattern2d_normalized` now scales its grid
  with a single bulk `f64v_scale` instead of a per-element loop.
  Whole function on 64×64: **35.2 µs → 26.0 µs (−26%)**; the scale step itself
  **9.15 µs → 1.05 µs per 4096 elements (8.7×)**. The remainder is
  `pattern2d_max_intensity` and `pattern2d_new`'s zero-fill, neither of which has
  a SIMD form (see below).
- **`tests/simd.tcyr`** pins the contract that makes the swap safe: the result is
  **bit-identical** to the scalar loop — 512 mixed-magnitude values, plus ±0, ±inf,
  `DBL_MIN` and `DBL_MAX`, plus all 256 cells of a real pattern. Bit-fidelity to
  `rust-old/` is a core property of the port, so a toolchain change that made
  `f64v_scale` reassociate or fuse must fail a test rather than silently shift
  published optics results.
- **`wave/pattern2d_normalized_64x64`** benchmark (35 → 36).

### Fixed
- ⚠ **`f64v_scale` writes one f64 PAST the requested count when the count is
  odd** — it processes in pairs and rounds up. Measured for every odd n in 1..12.
  `pattern2d_normalized` passes `width * height`, which is odd for a 3×5 or 7×7
  pattern, so the naive call corrupted the 8 bytes after the grid. prakash now
  vectorises the even prefix and does the final element scalar; the suite pins
  both the overrun (so a toolchain fix is noticed) and that a 3×5 pattern
  normalizes correctly.

### Investigated and rejected
Each measured, not assumed — the same discipline as 2.1.2:

- **The typed `f64v2_*` / `f64v4_*` wrappers are SLOWER than scalar.** Per 4096
  elements: scalar 9.15 µs, `f64v2` **23.0 µs**, `f64v4` **20.3 µs**. Their ABI
  plus per-lane extraction costs more than the arithmetic saves. Only the raw
  `f64v_*` builtins are worth using.
- **Chunking loses most of the win**: `f64v_scale` called per 4 elements is
  2.83 µs/4096 against **1.05 µs** for one bulk call. It must span the range.
- **`memcpy` / `memset` are byte-at-a-time** in `lib/string.cyr` (`store8` in a
  `while`) — **7–8× SLOWER** than an 8-byte scalar loop (48.8 µs / 47.2 µs vs
  6.4 µs / 5.9 µs per 4096 f64). An earlier draft of this release used `memcpy`
  for the clone branch; it was a regression and was reverted.
- **The remaining hot spots cannot be vectorised** with what the stdlib provides:
  - `wave/interference_pattern_16x16` (108 µs) and `spectral/spd_blackbody`
    (9.8 µs) are dominated by per-element `hypot`/`sin`/`cos`/`exp`. There are
    **no vector transcendentals**.
  - `ray/spot_diagram` (28 µs) is branchy ray tracing, not an array kernel.
  - `pattern2d_max_intensity` needs a max reduction — **`f64v_max` does not
    exist** (nor `f64v_min`, `f64v_fill`, `f64v_sum`, `f64v_copy`).
  - `_spd_integrate` looks like a dot product, but the CIE CMF table is
    **interleaved** (x@0, y@8, z@16 per row, stride 24) so `f64v_dot` cannot read
    it, and a reassociated sum would not be bit-exact against `rust-old/` anyway.
  - The zero-fill loops have no `f64v_fill`, and `memset` is slower than the loop.

### Notes
- ⚠ **Toolchain constraint:** `f64v_scale` **segfaults when called from top-level
  statement scope**; it must be invoked from inside a function. Library code is
  unaffected (`pattern2d_normalized` is a function), but `.tcyr` suites run at
  top level, so `tests/simd.tcyr` routes every SIMD call through a helper.

## [2.1.2] - 2026-08-12 — hisab interop, measured

Settles both hisab adoption items. One shipped as a tested contract; the other
was investigated and does not apply. No behaviour changes, no signature changes —
`src/` is unchanged apart from documentation and one new test suite.

Suite **5424 → 5444** assertions across **27 → 28** suites.

### Added
- **`tests/hisab_interop.tcyr` — the `RayVec3` ↔ `HVec3` layout contract is now
  enforced.** Both are `{ x; y; z; }`, so a prakash vector **is** an hisab vector:
  all 26 of hisab's `hvec3_*` operations (dot, cross, normalize, reflect,
  project, angle, …) work on prakash vectors with **zero conversion**, and hisab
  results read back through the `RayVec3` accessors. The port has claimed this
  since 2.0.0 but nothing tested it — a field reorder or an added field in either
  project would have broken every consumer silently, at the point of a wrong
  number rather than a compile error. The suite pins:
  - `sizeof(RayVec3) == sizeof(HVec3)`, and the field **order** (x@0, y@8, z@16),
    which a size check alone would not catch
  - hisab ops on prakash vectors in both directions, with exact expected values
  - `hvec3_dot` bit-identical to prakash's inline dot
  - `ray_reflect_3d` agreeing with `hvec3_reflect` across **121 cases**, so a
    consumer may reach for either library's geometry and get the same answer

### Changed
- `src/ray_core.cyr` documents the interop guarantee at `struct RayVec3`, and why
  prakash does not use it internally.

### Investigated and rejected
- ⚠ **prakash does NOT call `hvec3_*` internally — measured, not assumed.**
  Routing the tracer's dot products through `hvec3_dot` (which is allocation-free
  and SIMD-backed, so it looked like a free win) cost **2–11%** across three runs:

  | Benchmark | inline | via `hvec3_dot` |
  |---|--:|--:|
  | `ray/trace_sequential` | 957–964 ns | 985–1068 ns |
  | `ray/trace_surface` | 329–334 ns | 340–360 ns |
  | `ray/fresnel_unpolarized` | 116–117 ns | 121–123 ns |

  `hvec3_dot` is two nested calls (`hvec3_dot` → `f64v_dot`) where the inline form
  is straight-line arithmetic; for a 3-element dot the call overhead exceeds any
  vectorisation gain. Reverted. `hvec3_reflect` is worse still for the hot path —
  it composes `hvec3_sub(hvec3_scale(…))`, **two allocations** where prakash's
  does one.
- ⚠ **hisab quadrature does not apply to the named targets.** hisab does ship
  `calc_integral_gauss5(f, a, b)`, `calc_integral_simpson` and
  `calc_adaptive_simpson`, but:
  - `_spd_integrate` is the **CIE-defined weighted sum** over the tabulated
    81-entry CMFs at 5 nm — the standard's method, not an approximation to
    improve. Replacing it would change published colour values and break 1168
    `spectral_cie` assertions.
  - `huygens_fresnel_1d` integrates a **caller-supplied discrete sample buffer**;
    there is no continuous integrand to hand a quadrature rule.
  - `spd_blackbody` **samples** 81 points, it does not integrate.

  `src/` contains **zero `fncall` sites** — no prakash function takes the callable
  integrand these routines require. A `spd_from_function(f, start, end)` entry
  point would be the one genuine fit; it is filed on the roadmap as a new
  capability, not as dependency adoption.

## [2.1.1] - 2026-08-12 — The audit's remaining defects and coverage gaps

Clears every item the 2.0.2 audit reproduced but deferred, plus the three coverage
gaps it identified. Suite **5361 → 5424** assertions; benchmarks **27 → 35**.
No signature changes to existing functions.

### Fixed
- **`f64_sin`/`f64_cos` returned the argument unchanged for |x| >= 2^63.** The x86
  `fsin`/`fcos` instructions are undefined past that point, so `f64_sin(1e300)`
  evaluated to `1e300` — an "intensity" far outside [-1, 1] that propagated
  silently through interference and slit calculations. Measured: correct through
  2^62, broken at 1e300 and at the infinities. New `_prk_sin` / `_prk_cos` in
  `error.cyr` return NaN past the limit so the loss is visible; every trig call
  site in `src/` routes through them and in-domain values are untouched.
  ⚠ Rust's libm would perform a Payne-Hanek reduction and return a real value —
  this is a deliberate divergence in favour of not fabricating one. At 2^63 a
  double's spacing exceeds 2048 against a period of 2π, so the argument no longer
  identifies a phase.
- **`planck_radiance` lost all precision for small arguments.** The port
  substituted a bare `exp(x) - 1` for Rust's `f64::exp_m1`, which returns exactly
  **0** once x drops below about 1e-16 — the entire result destroyed by
  cancellation. `_expm1` now uses Kahan's identity `(u-1)·x / ln(u)` with
  `u = exp(x)`, where the same `u` in numerator and denominator cancels the
  rounding error. Verified: `expm1(1e-18)` returns 1e-18 (was 0), `expm1(1e-14)`
  and `expm1(1e-6)` match the series, and `expm1(10)` is unchanged.
- **`point_source_new` could not build input for its only consumer.**
  `interference_pattern` indexes a contiguous `PointSource` array
  (`sources + i*sizeof(PointSource)`), but the constructor returns standalone
  allocations that are not contiguous — the test suite carried a private `_mksrc`
  helper to work around it. Added **`point_source_array_new(n)`**,
  **`point_source_set(arr, i, …)`** and **`point_source_at(arr, i)`**;
  `point_source_new` remains for a single standalone source and now says so.
- **46 constructor allocations stored through an unchecked result.** All are
  fixed-size (`alloc(sizeof(T))` and friends), so only heap exhaustion reaches
  them — but they returned an unguarded pointer, which meant a caller could not
  detect the failure at all. Each now returns 0.
  ⚠ The **10 private lazy-init table getters** in `spectral_cie` /
  `spectral_photometry` are deliberately left unguarded: their callers index the
  result immediately (`_cmf_x(t, i)`), so returning 0 would convert a clean
  segfault into a wild read. Documented rather than changed.

### Added
- **8 benchmarks for the expensive composites**, which the suite had never
  measured (27 → 35). The audit was right that this was hiding the real cost:

  | Benchmark | avg |
  |---|--:|
  | `wave/interference_pattern_16x16` | 107.9 µs |
  | `ray/spot_diagram` | 27.9 µs |
  | `spectral/color_rendering_index` | 11.1 µs |
  | `spectral/spd_blackbody` | 9.8 µs |
  | `spectral/luminous_flux` | 3.7 µs |
  | `spectral/spd_to_xyz` | 1.5 µs |
  | `ray/trace_sequential` | 932 ns |
  | `ray/trace_surface` | 327 ns |

  `interference_pattern` alone is **7× the slowest previously-benchmarked
  operation** (`diffraction_pattern_2d_8x8`, 15.4 µs). `tests/prakash.bcyr` now
  includes `spectral_photometry` and `ray_simulate`, which it had never compiled.

### Tests
- **`trace_surface`'s five error branches are now all exercised** — previously
  only the aperture case was, leaving the four geometry-miss paths the whole
  tracing stack depends on untested: plane-parallel ray, plane behind the ray,
  sphere missed (negative discriminant), and sphere entirely behind. Plus a
  control that a successful trace clears a stale `err_out`.
- **Mueller accessors: full 16-element set/get roundtrip**, including that every
  `(i, j)` addresses a distinct slot (no aliasing) — they had no direct tests
  before 2.0.2 added the bounds cases.
- New coverage for the `PointSource` array API, the trig domain guard (including
  a 2^62 control just below the limit), and `expm1` at 1e-18 / 1e-14 / 1e-6.

## [2.1.0] - 2026-08-12 — Error channels for the deserialization surface

2.0.2 stopped the `*_from_json` family from crashing or fabricating on malformed
input, but left it unable to say **why** a decode failed — it returned a bare `0`,
which broke the library's own convention in `error.cyr` ("functions return
`PK_ERR_NONE` on success or a negative `PK_ERR_*` code on failure"). This release
closes that gap.

Suite **5334 → 5361** assertions across 27 suites.

### Breaking
- **Every `*_from_json` takes a trailing `err_out` pointer.** Nine functions:

  | Function | Was | Now |
  |---|---|---|
  | `rgb_from_json` | `(json, len)` | `(json, len, err_out)` |
  | `medium_from_json` | `(json, len)` | `(json, len, err_out)` |
  | `sellmeier_from_json` | `(json, len)` | `(json, len, err_out)` |
  | `lens_type_from_json` | `(json, len)` | `(json, len, err_out)` |
  | `polarization_from_json` | `(json, len)` | `(json, len, err_out)` |
  | `prescription_from_json` | `(json, len)` | `(json, len, err_out)` |
  | `spd_from_json` | `(json, len)` | `(json, len, err_out)` |
  | `ai_daimon_config_from_json` | `(json, len)` | `(json, len, err_out)` |
  | `ai_hoosh_config_from_json` | `(json, len)` | `(json, len, err_out)` |

  **Migration** — pass the address of a local and test it, exactly as the rest of
  the library already works:

  ```cyrius
  # before (2.0.2)
  var c = rgb_from_json(data, len);
  if (c == 0) { ... }                  # failed, reason unknown

  # after (2.1.0)
  var e = 0;
  var c = rgb_from_json(data, len, &e);
  if (prakash_is_err(e) == 1) { ... }  # e says which of the three below
  ```

  The return value is unchanged in every case (`0` on failure;
  `lens_type_from_json` still returns `LENS_TYPE_INVALID`), so a caller that only
  checks the return value needs nothing but the extra argument.

### Added
- **`PK_ERR_PARSE`** (`-9`) — the input bytes are not well-formed JSON, or the
  document is not the type the decoder needs. Distinct from
  `PK_ERR_INVALID_PARAMETER`, which now means the document parsed but a required
  field is absent or of the wrong type. **Bad bytes vs. wrong schema** is the
  distinction a caller actually needs, and the suite asserts the two are
  reported separately.
- The three codes the family emits: `PK_ERR_PARSE`,
  `PK_ERR_INVALID_PARAMETER`, and `PK_ERR_ALLOCATION` (for a caller-sized buffer
  that could not be allocated — `spd_from_json`'s sample array,
  `prescription_new`, and the `ai` config handles).

### Changed
- `err_out` doubles as the failure accumulator inside a decoder: each field
  helper writes its own code and success never writes, so several fields are read
  in a row and the slot is tested once. Success always clears the slot to
  `PK_ERR_NONE` on entry, so a stale code from an earlier call cannot be mistaken
  for a fresh failure — asserted directly.
- **`ai_daimon_client_new` and `ai_register_agent` now report allocation failures
  as `PK_ERR_ALLOCATION`** rather than `PK_ERR_INVALID_PARAMETER`; the arguments
  were fine, the memory was not there. `ai_daimon_client_new`'s "never fails"
  comment was false as of 2.0.2 and is corrected.

## [2.0.2] - 2026-08-12 — Audit release: crashes, fabricated results, and a corrupted round-trip

A full P(-1) hardening sweep of the 2.0.1 tree across six dimensions — allocation
safety, error channels, untrusted input, numerical fidelity against `rust-old/`,
API/doc accuracy, and coverage — with every candidate finding handed to an
independent skeptic instructed to **refute** it. Two findings were refuted and
dropped; several had their severity or trigger corrected by the verifier.

No API signature changed, so this is a drop-in upgrade. Suite **5251 → 5334**
assertions across **26 → 27** suites (new `tests/hardening.tcyr` with 71, plus 12
added to `tests/ai.tcyr`).

**Every defect below was reproduced against 2.0.1 before it was fixed**, and the
crash cases were reproduced as an actual `SIGSEGV` (exit 139) from a public entry
point — not inferred from reading. `tests/hardening.tcyr` re-runs each one.

### Security
- **Seven public entry points aborted the process on malformed input.**
  `medium_from_json`, `lens_type_from_json`, `prescription_from_json` (all in
  `serialize`) and `ai_daimon_config_from_json` / `ai_hoosh_config_from_json`
  (in `ai`) dereferenced null and SIGSEGVed on JSON that merely lacked a required
  field or gave it the wrong type — `{}`, `{"endpoint":42}`, `not json` all
  sufficed. Every bayan accessor returns 0 for a missing or wrong-typed node and
  `str_data(0)` is `load64(0)`. These are the deserialization paths, i.e. exactly
  where untrusted bytes arrive. All now validate and fail safely.
- **Unchecked allocations turned a size request into a process abort.**
  `alloc` returns 0 past `ALLOC_MAX` (2 GiB) *and on heap exhaustion*, and none
  of the pattern/zernike allocation sites checked it before storing through the
  result. The cheapest trigger took two scalar arguments:
  `psf_diffraction_limited(8193, 100.0)` — the power-of-two padding rounds 8193
  to 16384, making the 16-byte-per-cell FFT grid 4 GiB. Guarded in
  `pattern2d_new`, `diffraction_pattern_2d`, `diffraction_pattern_circular`,
  `psf_from_wavefront`, `_pat_fft2d`, `spectrum_strip`, `spectrum_strip_range`,
  `zernike_wf_to_grid`, `point_source_new` and the `ai` constructors.
- **Size computations could overflow before reaching `alloc`.**
  `pattern2d_new(1, 2305843009213693953)` wrapped `width*height*8` to 8, allocated
  an 8-byte buffer, and then wrote 2^61 f64s through it. New `_pat_dims_bad`
  screens non-positive dimensions, multiplication overflow, and the `ALLOC_MAX`
  budget *before* the product is formed. `_pat_npot` no longer doubles into a
  negative.
- **`ai_register_agent` POSTed to a corrupted URL.** `str_cat` returns a `Str`
  whose buffer has no NUL terminator, but sandhi consumes the URL as a cstring
  and read off the end. Now terminated via `str_cstr`.

### Fixed
- **`pattern2d_get`/`pattern2d_set` and `mueller_get`/`mueller_set` did no bounds
  checking**, so an out-of-range index read or wrote outside the buffer — the
  Rust original used a checked `[[f64; 4]; 4]` and `Vec`. Reads now return 0.0 and
  writes are refused with **`PK_ERR_INVALID_PARAMETER`** — a named code, not a
  bare `-1`, which is already `PK_ERR_TIR` and would have reported "total
  internal reflection" for an out-of-range matrix index.
- **`spd_at` read `values[-1]` on a degenerate SPD.** With `len = 0` and
  `step = 0`, `idx_f` is `0/0 = NaN`, which reached the `idx >= len - 1` branch
  and loaded 8 bytes before the buffer. (The Rust original indexed a `Vec` and
  panicked rather than reading out of bounds.)
- **Malformed JSON silently fabricated a valid-looking result.** `rgb_from_json`
  on garbage returned `Rgb(0,0,0)` — indistinguishable from a legitimate black —
  because bayan returns 0 for a missing node and 0 *is* the bit pattern for 0.0.
  All `*_from_json` now return 0 on failure; `lens_type_from_json` returns the new
  `LENS_TYPE_INVALID`, since 0 was already `LENS_CONVERGING`.
- **A flat optical surface was corrupted by a JSON round-trip.** JSON has no
  infinity literal, so bayan encodes a non-finite f64 as `null`, and the decoder
  read it back as 0.0 — turning a plano surface (radius `+inf`, *zero* curvature)
  into a radius-0 surface (*infinite* curvature), the exact physical inverse.
  A plano-convex lens, one of the most common elements, hit this. `radius` now
  restores `+inf` from `null`.
- **Three math boundaries diverged from Rust**, because `ganita_f64_pow` is
  `exp(n*ln(base))` and the x87 `exp` returns NaN for an infinite argument:

  | expression | 2.0.1 | Rust | 2.0.2 |
  |---|---|---|---|
  | `pow(0.0, 2)` | NaN | 0.0 | 0.0 |
  | `exp(-inf)` | NaN | 0.0 | 0.0 |
  | `exp(+inf)` | NaN | +inf | +inf |

  New `_prk_pow` / `_prk_exp` shims in `error.cyr` carry Rust's semantics; all
  call sites route through them and finite arguments are delegated unchanged.
  The reachable consequences: `pbr_distribution_charlie` returned **NaN at the
  specular peak** (`n_dot_h = 1`, where `sin_theta` is 0) instead of 0.0, and
  `atm_atmospheric_transmittance` returned NaN instead of 0.0 at saturation.
- **`color_temperature_to_rgb`** tested `temp > 66` where Rust tests
  `temp <= 66` and takes the else branch. Identical for finite input, but on a
  NaN both comparisons are false and Rust computes while the port returned 1.0.
  Now negated to match (and to match `bridge.cyr`, which already did).

### Performance
- **Two hot paths allocated per call into an allocator that never frees.**
  `zernike_wf_evaluate` took a 16-byte scratch buffer on every call and is called
  once per grid cell by `zernike_wf_to_grid`, so a grid leaked more memory than it
  returned; `fresnel_integral_c`/`_s`/`_cs`/`fresnel_edge_intensity` did the same
  per call. Both now use stack locals via `&`-pointers — no allocation.
  - `wave/fresnel_integral_c`: 157 ns → **136 ns (−13%)**
  - `wave/zernike_poly`: 110 ns → 103 ns
- **The added guards cost nothing measurable.** Despite per-call bounds checks in
  `pattern2d_set` (invoked once per cell), `wave/diffraction_pattern_2d_8x8` is
  15.888 µs → 15.400 µs, and every other row is flat or slightly faster. Full
  27-benchmark run recorded in `bench-history.csv` / `benchmarks.md`.

### Added
- **`PK_ERR_ALLOCATION`** (`-8`) — an allocation was refused because the request
  exceeds `ALLOC_MAX`, the size computation overflowed, or the heap is exhausted.
  Distinct from `PK_ERR_INVALID_PARAMETER`: the arguments were sane, the memory
  was not there. Every status-returning site added in this release now returns a
  named `PK_ERR_*` code, per `error.cyr`'s stated convention.

### Changed
- **`#must_use`** added to ten pure functions that lacked it (`prakash_is_err`,
  `prakash_is_ok`, the five `spec_*` constants, the three `fraunhofer_*` lines)
  plus `mueller_get`.
- **Removed dead code**: `_atm_n_air` had no caller (the Rayleigh prefactor it
  documents is a compile-time constant).

### Documentation
- **`docs/architecture/math.md` had three wrong formulas**, each verified against
  both the implementation and `rust-old/`: the Petzval radius read `R_p = -1/(2Σ)`
  where the code computes `-1/Σ`; the McCamy CCT denominator was `(y - 0.1858)`
  where both the code and the published formula use `(0.1858 - y)`; and the
  double-slit intensity omitted the factor of 4.
- **README** advertised a CIE illuminant `E` that does not exist (only A, D50,
  D65, F2, F11 are implemented).
- **`num_ifft`** was described as part of prakash's hisab surface; the library
  calls only `num_fft` (the test suite exercises `num_ifft`).
- ⚠ **Corrected two claims that shipped in the tagged 2.0.1 entry**: that `cyrius
  deps` hard-errors on a sidecar mismatch, and that the over-reporting sidecar
  would force math-only consumers to compile TLS. Both came from the upstream
  issue report rather than measurement, and neither reproduces on cyrius 6.5.20
  (a consumer declaring fewer folds resolves cleanly; a bogus fold name in the
  sidecar raises no error). The 2.0.1 entry below carries an amendment note. The
  sidecar fix stands; its blast radius does not.

## [2.0.1] - 2026-08-12 — Toolchain and dependency refresh

Maintenance release: Cyrius `6.4.10` → **`6.5.20`**, hisab `2.6.8` → **`2.11.1`**,
sakshi `2.4.2` → **`2.4.10`**. No optics behaviour changes — all **5251 test
assertions across 26 suites** pass unchanged, and both bundles differ from 2.0.0
only by the version stamp and one dependency-driven rename.

### Fixed
- **serialize / ai** — `bayan_json_v_parse_str(buf, len)` was renamed upstream to
  **`bayan_json_v_parse_buf(buf, len)`** (bayan 1.3.0); the old name no longer
  exists, which broke the `serialize` and `ai` suites with `undefined function`.
  Updated all 10 call sites (7 in `src/serialize.cyr`, 3 in `src/ai.cyr`).
  Semantics are identical — the rename was forced by Cyrius's `X_str` overload
  dispatch, which silently rewrote `bayan_json_v_parse(someStr)` into a 1-arg
  call to the 2-arg `_str` function. prakash always called the explicit
  `(buf, len)` form, so it never hit that mis-dispatch; this is a pure rename.
- **`dist/*.deps` sidecars described the wrong bundles.** These files name the
  stdlib folds each published bundle needs. Cyrius 6.5.10 made the **base**
  bundle's sidecar `[deps] stdlib` ∪ include-scan while **profile** bundles keep
  a pruned inference — a rule that assumes the base bundle is the widest.
  prakash is inverted (`[lib]` is the narrow math-only core; `[lib.ai]` is the
  wide one), so `dist/prakash.deps` advertised the whole sandhi/net/TLS stack
  that `[deps] stdlib` must declare for `src/ai.cyr` — despite the core bundle
  referencing **none** of it (verified by symbol scan) — while
  `dist/prakash-ai.deps` listed only `syscalls io`, omitting sandhi entirely.

  prakash now owns both via **`scripts/sync-deps-sidecar.sh`** (core = declared
  stdlib − the AI-only folds = 18; ai = declared stdlib in full = 27), following
  the precedent setu set for the same upstream issue. CI runs it after `distlib`
  and fails on drift, plus a **core-bundle-is-TLS-free** symbol scan.
  `.gitignore` no longer lists `dist/*.deps` — they were ignored *and* tracked,
  an incoherent state; they are deliberately tracked.

  ⚠ **Scope, measured rather than assumed:** the sidecar is published metadata,
  **not** an enforced contract. Against cyrius 6.5.20, a consumer declaring fewer
  folds than the sidecar lists still resolves cleanly (exit 0); a bogus fold name
  in the sidecar raises no error; and the folds vendored into a git-dep consumer
  are driven by prakash's own `cyrius.cyml`, not by the sidecar. So no consumer
  was actually forced to compile TLS. The upstream issue report describes a
  hard-error behaviour that does not reproduce here.

  ⚠ **This paragraph is an amendment.** As tagged, 2.0.1 asserted that `cyrius
  deps` hard-errors on a sidecar mismatch and that the over-reporting sidecar
  would force math-only consumers to compile TLS. Both claims were taken from the
  upstream issue report rather than measured, and both are wrong for cyrius
  6.5.20; they were corrected while preparing 2.0.2. The sidecar fix itself
  stands — published metadata should be true — but its blast radius does not.

### Changed
- **toolchain** — Cyrius `6.5.20`. `lib/` re-vendored via `cyrius lib sync --full`
  (99 → 102 files; adds `async_win.cyr`, `thread_macos.cyr`, `unicode/`), and
  `cyrius.lock` regenerated (99 → 108 locked deps, all verifying).
- **hisab** `2.6.8` → `2.11.1`. prakash's surface is unchanged — `num_fft` /
  `num_ifft`, both still returning `HSB_ERR_NONE` (0). The 2.6.8
  collision-hardening properties this port depends on still hold. Scanned for
  global collisions against the larger bundle (509 prakash fns / 31 globals vs
  912 / 220 in hisab): **none**, in all four name cross-products — which matters
  because Cyrius globals are last-one-wins with no diagnostic.
- **sakshi** `2.4.2` → `2.4.10`, matching what the 6.5.20 stdlib bundles; this
  clears the `./lib/ shadows version-pinned …` warning.
- **`scripts/bench-history.sh`** — the generated `benchmarks.md` header now warns
  that the Δ column spans the 6.5.20 harness change and is not like-for-like.

### Performance
- **No prakash code changed, but every benchmark reports a much lower number.**
  Cyrius 6.5.20's harness measures its own timer floor (~1.32–1.35 µs) and
  subtracts it from every sample; previously that fixed cost was included in each
  reading. Scalar benches therefore drop ~90% (e.g. `pbr/fresnel_schlick`
  1.341 µs → 17 ns) **without getting faster**. Verified by reconstruction: the
  old overhead derived from the near-zero-cost rows is ~1,320 ns, and applying it
  to `atmosphere/sky_color_rgb` predicts 967 ns against 965 ns measured.
  To compare any pre-6.5.20 row, subtract ~1.32 µs from it first.
- **The real Rust-vs-Cyrius gap is now measurable: 3.2×–19.7×, median ~7.2×** —
  not the 60×–1,300× the old instrument implied.
  `docs/benchmarks-rust-v-cyrius.md` is rewritten against true op cost. The
  compute-bound FFT row, which the floor never touched, is unmoved (~3× → 2.8×),
  which cross-checks the correction.
- ⚠ **`serialize/rgb_to_json` regressed ~2×** (floor-corrected ~1,520 ns → ~3,400 ns),
  reproduced across two runs. Cause is upstream: bayan 1.4.1 rebuilt the
  serializer so each append routes through an allocator variant with per-append
  failure propagation (`_jb_append_string` → `_jb_append_string_a(default_alloc(), …)`).
  prakash's `rgb_to_json` body is unchanged. Math paths are unaffected; tracked on
  the roadmap.

## [2.0.0] - 2026-07-06 — Full port to the Cyrius systems language

prakash is rewritten from Rust to **Cyrius**, the AGNOS systems language. Every
optics module — ray, spectral, wave, lens, PBR, atmosphere, bridge, AI, plus
error/logging — was translated line-by-line from the original Rust source
(retained in-tree under `rust-old/` as the translation reference) into
hand-verified Cyrius. **25 science modules, 5251 test
assertions across 26 suites, and 27 benchmarks are green.** The math is
bit-faithful to the Rust original: f64 constants are encoded as exact ratios or
IEEE-754 hex bit patterns, `powi` replicates the compiler's square-and-multiply,
and left-associative fold order is preserved so results match to the ULP. Every
module was checked by an adversarial multi-agent verification pass.

This is a hard break: consumers no longer link a Rust crate — they include a
self-contained Cyrius bundle. soorat, kiran, and ranga remain Rust for now; the
`dist/prakash*.cyr` bundle is the interop surface until they port.

### Breaking
- **Language & ABI** — the entire library is reimplemented in Cyrius. There is no
  Rust crate, no `Cargo.toml`, no crates.io publish. The Rust API (methods,
  generics, traits, `Result<T>`, `Option<T>`, `Cow`, `Vec`) is replaced by flat
  Cyrius **free functions**. Migration: link a bundle and call the `*_` functions.
  - `Result<f64>` → an `err_out` pointer parameter (`store64(err_out, PK_ERR_*)`,
    returns 0 on error); `Option<T>` → returns `1`/`0` plus an out-pointer.
  - Tuples / `[f64; N]` returns → caller-supplied out-pointer buffers.
  - `f64` values are IEEE-754 bit patterns in `i64`; all arithmetic goes through
    `f64_add`/`f64_mul`/`f64_sqrt`/… (comparisons return `1`/`0`).
  - Structs are heap-allocated (`alloc(sizeof(T))` + `#derive(accessors)`); enums
    become tag constants.
- **error constants** — namespaced `PK_ERR_*` (`PK_ERR_NONE`, `PK_ERR_TIR`,
  `PK_ERR_DIVISION_BY_ZERO`, …) so they never collide with a dependency's globals
  under Cyrius's last-one-wins global rule. `PK_ERR_NONE` is 0.

### Added
- **Two distlib bundles** — `cyrius distlib` emits `dist/prakash.cyr` (math-only,
  no TLS) and `cyrius distlib ai` emits `dist/prakash-ai.cyr` (same core **plus**
  the daimon/hoosh AI client, which pulls the sandhi HTTP/TLS stack). Consumers
  pick one — mirroring the old optional `ai` cargo feature.
- **logging** — the Rust `tracing::trace!` diagnostics are ported onto **sakshi**:
  `prakash_set_log_level(level)` (SK_ERROR=1 … SK_TRACE=5, default SK_INFO leaves
  traces off) and 31 operation-entry `_prk_trace` markers. sakshi_trace is
  compile-elidable and runtime-guarded — zero measured overhead when disabled.
- **serialization** — `src/serialize.cyr` reimplements serde_json roundtrips on
  **bayan** (`*_to_json` / `*_from_json`) for the seven tested types (rgb, medium,
  sellmeier, lens_type, polarization, prescription, spd). Float rendering is
  6-decimal (roundtrips checked to ~1e-5).
- **benchmarking** — `tests/prakash.bcyr` (27 benches across every module) plus
  `scripts/bench-history.sh`, which records a CSV history and a 3-point trend
  `benchmarks.md`.

### Changed
- **Dependency stack** — first-party **hisab** (git dep, tag 2.6.8) supplies
  Complex + FFT (`num_fft`/`num_ifft`) for `wave/pattern`; **ganita** (stdlib)
  supplies transcendentals and linear algebra; **bayan** replaces serde_json;
  **sandhi** replaces reqwest/tokio for the AI client's blocking HTTP POST;
  **sakshi** replaces `tracing`. No external Rust crates remain.
- **toolchain** — Cyrius `6.4.10`; `VERSION` is the single source of truth
  (`cyrius.cyml` pulls it via `${file:VERSION}`).

### Removed
- **Rust build system** — `Cargo.toml`, `criterion` benches, `cargo`-based CI, and
  the crates.io publish flow. The Rust source itself is retained in-tree under
  `rust-old/` as a translation reference (to be removed a couple releases out) and
  is also preserved in pre-2.0 git tags.
- **bijli-backend feature** — already dropped in 1.2.0; the port keeps optics math
  self-contained (the `bridge` module exposes primitive-value cross-crate hooks).

## [1.2.0]

### Added
- **bridge** — cross-crate primitive-value bridges for bijli (EM frequency/wavelength conversion, E-field to intensity, Cauchy refractive index), tara (stellar temperature to RGB, B-V color index, spectral class, absolute magnitude to luminosity), badal (density to Rayleigh scale, humidity to Mie scale, cloud cover to diffuse fraction, visibility to extinction)

### Removed
- **bijli-backend feature** — removed direct bijli dependency and all `#[cfg(feature = "bijli-backend")]` code. Cross-crate coupling now uses primitive-value bridges instead. Removed: bijli type re-exports (JonesVector, GaussianBeam, AbcdMatrix, etc.), From impls between prakash/bijli types, Snell/Fresnel/Brewster delegation to bijli. All optics math is now self-contained in prakash.

### Updated
- iri-string 0.7.11 -> 0.7.12

## [1.1.0] - 2026-03-25

### Changed — P1 Post-V1.0 Audit
- **Polarization consolidation**: added `From<Polarization> for StokesVector` (Jones → Stokes conversion)
- **Duplicate function docs**: cross-referenced `diffraction_limit`/`rayleigh_criterion` and `beer_lambert`/`volume_transmittance`
- **Spd zero-alloc illuminants**: `Spd.values` changed from `Vec<f64>` to `Cow<'static, [f64]>`; illuminant functions now borrow static data via `Spd::from_static()`
- **`#![warn(missing_docs)]`** enabled — all public items now documented
- **Wavelength parameter ordering**: standardized wavelength-first convention across all wave/diffraction functions (`single_slit_intensity`, `double_slit_intensity`, `grating_maxima`, `airy_pattern`, `path_to_phase`, `fraunhofer_rect`, `fresnel_number`, `fresnel_parameter`, `coating_reflectance`, `multilayer_reflectance`)
- **Unit-suffix convention** documented in lib.rs: `_nm`/`_m`/`_um` for specific units, unsuffixed for generic

### Updated
- hisab 0.22 → 1.1
- criterion 0.5 → 0.8 (migrated from `criterion::black_box` to `std::hint::black_box`)
- reqwest 0.12 → 0.13 (switched from native-tls to rustls)

### Added — Bijli EM Backend (`bijli-backend` feature)
- **Dependency**: bijli 0.24 (path dep during development, optional `bijli-backend` feature, default on)
- **Scalar optics delegation**: `snell`, `critical_angle`, `brewster_angle`, `fresnel_normal`, `fresnel_unpolarized` delegate to bijli's EM-correct implementations with trig-free Fresnel variants
- **Constants**: `SPEED_OF_LIGHT` re-exported from `bijli::field` (single source of truth)
- **Polarization bridge**: bidirectional `From` impls for `StokesVector` and `MuellerMatrix`, `From<Polarization> for JonesVector`
- **New re-exports**: `JonesVector`, `JonesMatrix`, `Complex` (Jones formalism), `GaussianBeam`, `AbcdMatrix`, `ResonatorStability` (beam optics), `mie`, `MieResult`, `em_rayleigh_cross_section` (exact Mie scattering)
- **Material bridge**: `Medium::permittivity()` (n² = ε_r)
- **Error bridge**: `From<BijliError> for PrakashError`
- **Fallback**: all original implementations retained behind `#[cfg(not(feature = "bijli-backend"))]`

### Added — P1 Accuracy & Completeness (from physics research audit)
- **Planck numerical stability**: replaced `exp()-1.0` with `exp_m1()`, Wien approximation for x>500, T≤0 guard
- **Complex Fresnel equations**: `ComplexMedium` type, `fresnel_s_complex`/`fresnel_p_complex`/`fresnel_unpolarized_complex`/`fresnel_normal_complex` for metals and absorbing media. Presets: gold, silver, copper, aluminum at 550nm
- **Zernike polynomials** (`wave::zernike` module): radial polynomial R_n^m, full Z_n^m evaluation, Noll index ↔ (n,m) conversion, `ZernikeWavefront` (evaluate, to_grid, RMS error, peak-to-valley, Strehl ratio), named constructors (defocus, spherical, coma_x/y, astigmatism_0/45)
- **Polarization ray tracing**: `trace_sequential_polarized` tracks cumulative s/p Fresnel transmittance per surface, `PolarizedTraceHit` output with R_p/R_s ratio

### Added — Stack Integration
- **bijli 1.0**: upgraded from path dep to registry `bijli = "1"`, added `maxwell` feature
- **bijli wave re-exports**: `radiation_pressure_absorbed/reflected`, `poynting_vector`, `momentum_density`, `plane_wave_intensity`
- **bijli maxwell re-exports**: `wave_speed`, `impedance`, `free_space_impedance`, `refractive_index`
- **hisab FFT consolidation**: replaced inline Cooley-Tukey FFT in `wave::pattern` with `hisab::num::fft`

### Changed — Refactoring
- `ray/mod.rs` split: extracted `ray/fresnel.rs` (Fresnel + Brewster + ComplexMedium + Beer-Lambert)
- `wave/mod.rs` split: extracted `wave/coherence.rs`, `wave/airy.rs`, `wave/fabry_perot.rs`
- No file over 1,021 lines (down from 1,313 max)

### Tests
- 702 unit tests, 34 integration tests, 8 doc tests = 744 total (up from 598 + 10)
- Property-based tests (proptest): 8 properties covering Snell reversibility, Fresnel symmetry, Planck positivity, Beer-Lambert bounds, Zernike pupil, wavelength roundtrip, complex Fresnel dielectric equivalence
- Doc tests: 8 across all modules (ray, wave, spectral, lens, pbr, atmosphere, zernike)

## [1.0.0] - 2026-03-24

### Added — atmosphere V0.25: Atmospheric Optics
- Rayleigh scattering: cross-section (λ⁻⁴), sea-level coefficient, altitude-dependent coefficient, phase function
- Mie scattering: sea-level coefficient, altitude-dependent coefficient, Cornette-Shanks phase function (improved Henyey-Greenstein)
- Air mass: Kasten & Young (1989) formula (accurate through horizon)
- Optical depth: Rayleigh (wavelength-dependent), Mie, total atmospheric transmittance
- Sky color: single-scattering RGB model (Rayleigh + Mie), per-wavelength radiance
- Sunlight color: atmospheric transmittance filter (overhead white → sunset red)
- Sunset/sunrise gradient: view-angle-dependent sky color with sun/view path attenuation
- Scattering angle utility (zenith/azimuth → angular separation)
- New `atmosphere` feature gate

### Changed — V1.0: Stable Release
- API review: renamed `fresnel_c`/`fresnel_s`/`fresnel_cs` (diffraction integrals) to `fresnel_integral_c`/`fresnel_integral_s`/`fresnel_integral_cs` to resolve name collision with Fresnel reflectance functions in `ray` module
- API review: renamed `atmosphere::rayleigh_phase` to `atmosphere::phase_rayleigh` for consistency with `pbr::phase_rayleigh`
- API review: renamed `coating_reflectance` params from `n1/n2/n3` to `n_incident/n_coating/n_substrate`
- Documentation: added `///` doc comments to all ~24 previously undocumented public items (Medium constants, Rgb/Xyz constructors, StripColor fields, etc.)
- README: full rewrite with module table, feature flags, examples, architecture, consumer crates

### Added — examples V0.29: Examples
- `rainbow`: physically accurate rainbow via Sellmeier water dispersion, minimum deviation angles, ANSI color spectrum
- `camera_lens`: cemented doublet prescription, paraxial analysis, full trace, spot diagrams, OPD fan, recursive ghost reflections
- `pbr_materials`: Fresnel-Schlick vs angle, Cook-Torrance roughness sweep, clearcoat blend, iridescence, SSS profiles, sheen, phase functions, importance sampling

### Added — wave/pattern V0.28: Pattern Computation
- 2D diffraction pattern from arbitrary aperture via radix-2 Cooley-Tukey FFT (inline, no external dep)
- Circular aperture diffraction pattern (generates Airy pattern)
- N-source interference pattern on 2D grid (complex amplitude summation with 1/r falloff)
- Visible spectrum RGB strip generator (full 380–780nm range + custom wavelength range)
- PSF from wavefront error: pupil function × exp(i·k·OPD) → 2D FFT → intensity
- Diffraction-limited PSF convenience function
- `Pattern2D` grid type with get/set, max_intensity, normalized view

### Added — ray/system V0.27: Optical Bench
- Paraxial ray trace (y-nu method): `ParaxialRay`, `paraxial_refract`, `paraxial_transfer`, `paraxial_trace`
- Prescription builder: `Prescription` with fluent `.add_surface()` API, `to_trace_surfaces()`, `to_paraxial_surfaces()`
- System cardinal point finder: `find_system_properties()` — EFL, BFD, FFD, power, pupil positions from paraxial marginal ray trace
- Common prescriptions: `prescription_biconvex`, `prescription_planoconvex`, `prescription_doublet` (achromatic)
- Full serde support for `Prescription` and `PrescriptionSurface`

### Added — ray/simulate V0.26: Simulation Primitives
- Recursive ray tracer: follows both reflected and refracted paths at each interface, energy tracking, configurable max depth and min energy threshold
- Ray fan generators: meridional (y-z plane), sagittal (x-z plane), radial bundle (concentric rings), all with field angle support
- Spot diagram: traces ray bundle through optical system to image plane, computes RMS spot radius
- Optical path length/difference (OPD): per-ray OPL through system, OPD relative to chief ray, full fan OPD computation

### Changed — Code Quality & Performance Audit
- `#[must_use]` on all 190+ public pure functions across all modules
- `#[non_exhaustive]` on `SurfaceShape` and `LensType` enums (were missing)
- `tracing::trace!` instrumentation on key entry-point functions
- Error variants `DivisionByZero` and `InvalidParameter`: `String` → `Cow<'static, str>` (zero-alloc for static messages)
- `cook_torrance`: added explicit `h_dot_v` parameter (was using `n_dot_h` approximation)
- `MuellerMatrix`: added `Serialize`/`Deserialize` derives
- `Polarization::circular_right/left`: now `const fn` using `FRAC_1_SQRT_2`
- Performance: `snell_3d` inlined refraction (eliminated duplicate trig), `fresnel_unpolarized` cos(asin(x))→sqrt, `sin_cos()` in 5 functions, `fresnel_integral_cs` combined function, `lambert_diffuse` div→mul, precomputed Rayleigh prefactor, illuminant data as `pub const` arrays
- Canonical `RGB_WAVELENGTHS_NM`/`RGB_WAVELENGTHS_M` constants in spectral module (650/550/450nm)
- Test suite: 487 → 608 tests
- Benchmarks: 137 → 162 functions

### Quality
- 608 tests (598 unit + 10 integration), 162 criterion benchmarks
- Zero `unwrap()`/`panic!()` in library code
- `#[must_use]` on all 190+ pure functions, `#[non_exhaustive]` on all public enums
- `tracing` instrumentation on key entry points
- Full P(-1) scaffold hardening + development loop audit

## [0.23.3] - 2026-03-23

### Added — wave V0.3: Wave Optics Expansion
- Coherence: temporal coherence length/time, spatial coherence angle/area, coherence ratio
- Circular aperture diffraction: Bessel J₁ (Abramowitz & Stegun rational approximation), Airy pattern, first zero, Rayleigh criterion
- Fabry-Pérot interferometer: transmittance (Airy function), finesse, FSR (Hz and wavelength), resolving power
- Stokes vectors with 7 polarization state constructors and degree of polarization
- Mueller matrices: identity, horizontal/vertical polarizer, arbitrary-angle polarizer, quarter/half-wave plates, general retarder, rotation, matrix multiply, chain application
- Birefringent materials: calcite, quartz, rutile, mica presets; retardation, quarter/half-wave thickness, Mueller generation
- Fraunhofer diffraction: rectangular aperture, 1D arbitrary aperture (numerical DFT)
- Fresnel diffraction: Fresnel number, Fresnel C(x)/S(x) integrals, straight-edge intensity, Fresnel parameter, Huygens-Fresnel 1D numerical integral
- Anti-reflection coatings: ideal AR index, quarter-wave thickness, single-layer reflectance, V-coat, multi-layer transfer matrix method

### Added — pbr V0.4: Advanced PBR
- Anisotropic GGX NDF and Smith geometry (directional roughness)
- Sheen: Charlie distribution (Estevez & Kulla 2017) and Ashikhmin velvet model
- Clearcoat: GGX distribution, fixed-IOR Fresnel, Kelemen geometry, energy-conserving blend
- Subsurface scattering: Burley normalized diffusion profile, Gaussian profile, SSS diffuse term, thin-slab transmittance
- Iridescence: thin-film Fresnel (Airy formula) at single wavelength and RGB, angle/thickness/wavelength-dependent color
- Volumetric scattering: Henyey-Greenstein phase function, Rayleigh phase, isotropic phase, extinction, transmittance, single-scatter albedo, in-scattering
- GGX importance sampling: half-vector sampling, PDF, cosine-weighted hemisphere sampling
- Environment map: split-sum scale/bias (Lazarov analytical fit), mip LOD from roughness, numerical BRDF LUT integration (Hammersley sequence)

### Changed — Infrastructure
- Module refactoring: wave, ray, pbr, spectral split into submodules (largest file: 2174→956 lines)
- CI/CD: GitHub Actions workflows (ci.yml, release.yml) — fmt, clippy, test (Linux/macOS), MSRV 1.89, coverage, docs, security audit, deny, semver check
- Added: SECURITY.md, CONTRIBUTING.md, CODE_OF_CONDUCT.md, deny.toml, codecov.yml
- Added: scripts/bench-history.sh (CSV tracking + benchmarks.md generation), scripts/version-bump.sh
- Makefile: added `all`, `semver` targets
- Test suite: 330 → 487 tests
- Benchmarks: 95 → 137 functions
- Performance: SPD→XYZ 80% faster (aligned fast path), Fresnel edge 15%, CRI 41%, numerous #[inline] additions
- Bug fixes: Rgb::to_u8 rounding, depth_of_field hyperfocal, shape_factor div-by-zero, conjugate_factor div-by-zero

## [0.22.3] - 2026-03-23

### Added — ray V0.2: Optical Systems
- 3D vector refraction (`refract_3d`) and 3D Snell's law with Fresnel reflectance (`snell_3d`)
- Sequential ray tracing through multiple optical surfaces (`trace_surface`, `trace_sequential`)
- Optical surface types: `SurfaceShape` (Sphere/Plane), `OpticalSurface`, `TraceRay`, `TraceHit`
- Cauchy dispersion model (`CauchyCoefficients`) with BK7 and fused silica presets
- Sellmeier dispersion model (`SellmeierCoefficients`) with 6 presets (BK7, SF11, fused silica, sapphire, water, diamond)
- Abbe number calculation from Sellmeier coefficients
- Fraunhofer spectral line constants (D, F, C)
- Prism deviation, dispersion, and angular spread functions

### Added — lens V0.2: Advanced Optics
- Thick lens equation and cardinal points (FFD, BFD, principal planes)
- F-number, aperture diameter, numerical aperture, NA from f-number
- Diffraction limit (Rayleigh criterion) and Airy disk radius
- Field of view (horizontal and diagonal)
- Diffraction-limited MTF (cutoff frequency and modulation curve)
- Seidel aberration coefficients (spherical, coma, astigmatism, field curvature, distortion)
- Shape factor and conjugate factor for lens analysis
- Longitudinal spherical aberration, chromatic aberration
- Petzval sum and Petzval radius for field curvature analysis
- Separated two-lens system (focal length, BFD)
- System magnification for multi-element systems

### Added — spectral V0.2: Color Science
- CIE 1931 2° standard observer color matching functions (81 entries, 380–780nm @ 5nm)
- `Xyz` tristimulus type with XYZ↔xyY and XYZ↔sRGB conversions
- sRGB gamma correction (linear↔gamma) with proper IEC 61966-2-1 transfer function
- Linear sRGB↔XYZ matrix conversions (D65 white point)
- Correlated color temperature from xy chromaticity (McCamy's approximation)
- CIE CMF interpolation at arbitrary wavelengths
- `Spd` spectral power distribution type with interpolation, XYZ integration, and sRGB conversion
- Blackbody SPD generator
- Standard illuminants: D65, D50, A, F2, F11
- Color Rendering Index (CRI Ra) calculation

### Fixed
- ai: replaced `anyhow` dependency with native `PrakashError`
- ai: `DaimonClient::new()` returns `Result` instead of panicking
- spectral: `wavelength_to_rgb` rejects NaN input
- spectral: `Rgb::to_u8` now rounds instead of truncating (0.999→255)
- lens: `depth_of_field` returns `f64::INFINITY` at hyperfocal distance instead of negative
- lens: `shape_factor` guards against division by zero when r1 == r2
- lens: `optical_power` and `combined_focal_length` return `Result` for zero inputs
- ray: `critical_angle` uses static error string instead of allocating on hot path
- ray: Sellmeier `n_at` guards against resonance pole division-by-zero and negative n²
- pbr: fixed operator precedence in GGX, Beckmann, and geometry functions

### Changed
- Test suite: 92 → 330 tests
- Benchmarks: 24 → 95 benchmarked functions
- Performance improvements across all modules via `#[inline]`, precomputed constants, eliminated redundant computation
- `snell_3d` delegates to `refract_3d` (single source of truth)
- Sequential ray trace sphere intersection assumes normalized direction (eliminates redundant dot product)

## [0.1.0] - 2026-03-23

### Added
- ray: 12 material refractive indices, Snell's law with TIR, Fresnel equations (s/p/unpolarized/normal), reflection (2D/3D), critical angle, Brewster's angle, Beer-Lambert attenuation
- wave: interference (constructive/destructive detection), single/double-slit diffraction, diffraction grating, thin film interference, polarization (Jones vectors), Malus's law
- spectral: wavelength→RGB (CIE 1931), Planck blackbody radiation, Wien's displacement law, color temperature→RGB, wavelength↔frequency, photon energy (J and eV), physical constants
- lens: thin lens equation, magnification, lensmaker's equation, optical power, mirrors, combined focal length, lens classification, depth of field
- pbr: Fresnel-Schlick (scalar/RGB), GGX and Beckmann NDF, Schlick-GGX geometry, Smith geometry, Cook-Torrance specular BRDF, Lambert diffuse (scalar/RGB), IOR→F0
- error: PrakashError with #[non_exhaustive], 7 variants
