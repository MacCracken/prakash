# Upgrading a consumer of `dist/prakash.cyr`

For a project that vendors the bundle as a git dependency — today that is
[ranga](https://github.com/MacCracken/ranga), pinned at **2.2.8** for its optional
`spectral` profile. Everything below was measured at 2.5.1 against that pin; the
details of each change are in `CHANGELOG.md`.

## The short version (2.2.8 → 2.5.x)

- **Nothing ranga calls changed value or signature.** `xyz_new`, the `Xyz_*`
  accessors, `spd_to_xyz`, `cct_from_xy`, `cie_cmf_at`, `spd_blackbody`,
  `xyz_d65_white` / `xyz_d50_white`, `illuminant_d65`, `spec_visible_min_nm` /
  `_max_nm` and the `Spd_*` accessors have identical bodies apart from
  allocation-failure guards, and the CIE, illuminant and V(λ) tables have not moved
  a bit. ranga's pinned values (D65 CCT 6503.461953 K, D50 5002.115988 K) reproduce.
- **`struct Xyz { x; y; z; }` keeps its layout** (0 / 8 / 16). ranga's source reads
  it only through `Xyz_x` / `_y` / `_z` and copies it, but `tests/spectral.tcyr`
  pins the layout with raw `load64`, so a layout change would fail ranga's suite.
- **No name collisions.** 24 global names were added since 2.2.8 and none removed;
  none of them matches any of ranga's 784 source-level names or anything in its
  vendored `lib/`.
- **It builds without hisab, and on ranga's toolchain.** The bundle compiles and
  runs ranga's surface with one `undefined function 'num_fft'` warning.
  `scripts/check-consumer-link.sh` is the check. **CI runs it bare**, which covers
  the pinned cyrius 6.6.6 with prakash's own `lib/` on every commit. **ranga's exact
  configuration** — cyrius 6.6.2 and ranga's own vendored stdlib — was checked by
  hand at 2.5.1 and is not held by CI; to repeat it:
  `./scripts/check-consumer-link.sh --root ~/Repos/ranga ~/.cyrius/versions/6.6.2/bin/cycc`.

## What does move

### On ranga's re-exported surface

`src/spectral.cyr` in ranga promises that prakash's flat names are callable
through `ranga-spectral`, so these reach ranga's users even though ranga never
calls them:

| Function | Release | Change |
|---|---|---|
| `color_rendering_index` | 2.4.8 | Ra was wrong for sources whose CCT falls outside ~1630–11850 K: blackbody 1000 K 15.589 → 99.999, 1200 K 63.938 → 100.000, 1500 K 95.300 → 100.000. D65 / A / F2 / F11 unchanged. |
| `color_rendering_index` | 2.4.9 | Not a value change: it leaked 680 bytes per call for sources below 5000 K (#16). |
| `wavelength_to_rgb` | 2.4.9 | Reports `PK_ERR_ALLOCATION` through `err_out` on allocation failure instead of a null handle with `PK_ERR_NONE`. |

### Memory reclaim — new, and it fixes a live hazard

**At 2.2.8, calling `alloc_reset()` silently corrupted every colour value** that
followed: prakash's 21 lazily-built tables kept pointers into the reset heap, and
`cie_cmf_at(555)` returned ȳ = 0.0 instead of 1.0, with no error. 2.3.7 added
**`prakash_reset_caches()`**. Any code that reclaims the heap must call it after
`alloc_reset()` and before the next prakash call. ranga does not reset its
allocator itself; its users might.

⚠ **Use 2.5.1 or later if you reset.** From 2.3.8 through 2.5.0 the recipe was
incomplete: the 2D-diffraction scratch blocks (`diffraction_pattern_2d` /
`_circular`, `psf_*`) survived it and were then written into memory the allocator
had already handed back out — 12,416 of a fresh 16,384-word buffer overwritten,
and the pattern itself wrong. 2.5.1 covers them.

### Elsewhere in the same bundle

Linking the bundle links all of it. These moved published values in modules ranga
does not call — relevant only if a consumer starts to:

| Release | Module | Change |
|---|---|---|
| 2.3.0 | wave_pattern | hisab 3.x: `num_fft` returns a `Result`; the four FFT entry points need hisab ≥ 3.0 if linked |
| 2.4.2–2.4.4 | ray_dispersion | `sellmeier_diamond`, `sellmeier_water`, `schott_bk7`, `cauchy_fused_silica`, `sellmeier_sf11`, `sellmeier_sapphire` coefficients corrected |
| 2.4.5 | atmosphere | Rayleigh constants now describe one air (they were 12% inconsistent) |
| 2.4.6 | lens | `lens_longitudinal_spherical_aberration` bracket corrected: ×6.15 at n = 1.5 (4.46 / 3.45 / 2.79 at 1.6 / 1.7 / 1.8). The 2.5.0 factor applies on top |
| 2.4.8 | wave_core | **Breaking:** `polarization_circular_right` / `_left` swap S3 sign to the S3 > 0 = right convention |
| 2.4.8 | wave_diffraction | `coating_reflectance` gained a missing numerator term |
| 2.4.9 | pbr | isotropic GGX 19% low at its alpha floor; anisotropic GGX epsilon removed |
| 2.5.0 | lens, ray_simulate | **Breaking:** Seidel `spherical` ÷ n/(n−1)², `coma` × (n−1), LSA × f(n−1)²/n, OPD re-referenced (see the 2.5.0 entry) |

Across 2.3.5–2.4.9, 68+ public functions that dereferenced a null handle now return
their documented sentinel (0 / 0.0), reporting through `err_out` where the function
has one, instead of exiting 139.

## Toolchain and transitive dependencies

- prakash 2.5.x pins cyrius **6.6.6**, hisab **3.2.1** and sakshi **2.5.2**; 2.2.8
  pinned cyrius 6.5.33, hisab 2.11.2 and sakshi 2.4.11.
- The core bundle needs exactly the stdlib leaves in `dist/prakash.deps` (19 at
  2.5.1). ranga's `[deps] stdlib` already has all of them.
- ⚠ **Not verified:** whether `cyrius deps --features spectral` on a cyrius 6.6.2
  project enforces hisab 3.2.1's **cyrius ≥ 6.6.3** floor when it re-locks
  prakash's transitive deps. ranga never compiles hisab — it omits it on purpose —
  so the floor matters only if `deps` refuses the lock. Bumping ranga's own pin to
  6.6.6 sidesteps the question.

## Bundles

`dist/prakash.cyr` is math-only and TLS-free; `dist/prakash-ai.cyr` adds the sandhi
HTTP client for `ai.cyr`. Pick one. Both names, and the `.deps` sidecars beside
them, are unchanged since 2.2.8.
