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
  the pinned cyrius with prakash's own `lib/` on every commit. **ranga's exact
  configuration** — cyrius 6.6.2 and ranga's own vendored stdlib — was checked by
  hand at 2.5.1 and is not held by CI; to repeat it, name ranga's current pin:
  `./scripts/check-consumer-link.sh --root ~/Repos/ranga ~/.cyrius/versions/<pin>/bin/cycc`.
  ⚠ From 2.10.1 that passes in full only at cyrius ≥ 6.6.13; see below.

## 2.15.1: an audit, and what it moves

A repair release: results that were wrong change value, and inputs that used to come back as a
plausible number are refused. Nothing is renamed or removed, no signature changes, and the cyrius
pin stays at 6.7.18. `CHANGELOG.md` has the measurement behind each item.

**What ranga calls does not move on a valid handle.** `spd_to_xyz`, `spd_blackbody`,
`cct_from_xy`, `cie_cmf_at`, `xyz_new` and the illuminant tables give the same bits;
`spd_to_xyz` and `spd_blackbody` were refactored onto shared stack-buffer cores with the same
arithmetic. On a NULL handle, `spd_at` and `spd_end_nm` now answer NaN instead of 0.0 (see the
null-handle list below). That happens only after an `Spd` constructor or `spd_from_json` has
already reported its own failure. If you call `color_rendering_index` or `cri_special`, they now
follow CIE 13.3's whiteness limit: NaN for a source more than DC = 5.4e-3 from the Planckian locus
at its CCT, or below 1000 K. Every illuminant prakash ships passes.

**A null handle is NaN, not 0.0.** 41 f64 functions with no `err_out` answered 0.0 for a null
handle or a null buffer with a nonzero count. They are listed below. Each can also return 0.0 as a
computed answer (an index, a luminance, an R_i, a Mueller element), so the two could not be told
apart. They answer NaN now, which no valid input produces, and a caller that checked
`== 0` for failure should test `f64_eq(x, x) == 0`. The same goes for
three other non-answers:
- an index outside the object, in `mueller_get`, `pattern2d_get` and `cri_special`;
- an empty spot diagram, in `spot_rms_radius`, which said 0, a perfect focus;
- a polychromatic MTF with no weight.

The functions:
- `lens_mtf_polychromatic`, `lens_petzval_sum` and `lens_system_magnification`;
- `spot_rms_radius`;
- `fresnel_normal_complex`, `fresnel_s_complex`, `fresnel_p_complex` and
  `fresnel_unpolarized_complex`;
- `polarization_through_polarizer`, `polarization_intensity` and
  `stokes_degree_of_polarization`;
- `medium_permittivity` and `rgb_luminance`;
- the dispersion family: `cauchy_n_at`, `sellmeier_n_at`, `abbe_number`, `herzberger_n_at`,
  `schott_n_at`, `conrady_n_at` and `partial_dispersion`;
- `stokes_intensity`, `stokes_ellipticity_angle` and `stokes_orientation_angle`;
- `mueller_get`;
- the `biref_*` family: `biref_birefringence`, `biref_retardation`,
  `biref_quarter_wave_thickness` and `biref_half_wave_thickness`;
- `luminous_flux`, `luminous_efficacy` and `luminous_flux_scotopic` (also for an internal table
  that could not be allocated);
- `zernike_wf_evaluate`, `zernike_wf_rms_error` and `zernike_wf_strehl_ratio`;
- `spd_end_nm`, `spd_at`, `spectral_band_similarity`, `cri_special` and `color_rendering_index`;
- `pattern2d_get` and `pattern2d_max_intensity`.

Zeros that are physically correct stay: `spd_at` outside its sampled range (no power), a
wavefront outside its pupil, and the degree of polarization of unpolarized light.
Functions that return a pointer still return 0 for a null input, as a null pointer is
unambiguous; functions with an `err_out` still return 0 and write the code.

**Values that change, because the old one was wrong:**
- `thin_film_reflectance` is sin²(2πnd/λ). It was cos², so a vanishing film read 1.
- `polarization_through_polarizer` uses the Jones phase: circular light at 45° passes 1/2, not 1.
- `biref_mica` is n = 1.5936 / 1.5977 (Δn 0.0041). Every mica retardation was 8.2× too large.
- `airy_first_zero` is j₁,₁/π λ/D, 5e-7 relative above the old constant.
- `fresnel_integral_c`, `_s`, `_cs` and `fresnel_edge_intensity` are accurate to 3e-16. They
  were off by up to 2.2e-3.
- `multilayer_rt` with a layer at exactly its critical angle gives a number. s was NaN, and p
  dropped the layer.
- `trace_recursive` sends a ray crossing a surface toward −z into the medium in front of it.
  Trees with internal reflections gain the bounces they were missing.
- `ray_fan_sagittal` keeps the field in y and spreads its rays in x. It was the meridional fan,
  rotated.
- `find_system_properties`: `bfd` and `ffd` are true distances in the image and object media, a
  factor n′ and n from before. `entrance_pupil` and `exit_pupil` are NaN; a prescription has no
  stop.
- `optical_path_length` subtracts the leg back to an image plane behind the last surface. It
  used to drop that leg.
- `lens_mtf_through_focus` is the defocused pupil's OTF. It reads lower at every defocus than
  the sinc model did: at f/2.8, 200 c/mm and 10 µm, 0.452 rather than 0.606.
- `lens_mtf_diffraction_limited` is even in frequency. A negative ν read 1.
- `zernike_radial` is correct from n = 21 up, and a Zernike coefficient below 1e-15 counts.
- `pbr_iridescence_fresnel` / `_rgb` use signed amplitudes and average s and p. A thin soap film
  no longer reads at its quarter-wave maximum.
- `pbr_subsurface_diffuse` uses F_ss90 = roughness, Burley's value at l·h = 1. It was
  roughness², so values move wherever n·l or n·v < 1.
- A back-facing GGX half vector has D = 0.
- `beam_second_moment_radius` forms its moments in sample units, so values can differ in the
  last bit.

**Refused now, where a number came back before.**

Out-of-domain inputs:
- a wavelength ≤ 0 or NaN in the dispersion formulas, Planck, Strehl, interference and the PSF;
- a reflectance outside [0, 1) in Fabry-Pérot;
- a radius of 0 in the lens and system functions;
- |g| ≥ 1 in Henyey-Greenstein;
- an index that is ≤ 0, NaN or infinite;
- incidence on the back of a surface;
- |θ| ≥ π/2 in `multilayer_rt`;
- an unknown spectral class, or B − V ≤ −0.674, in `bridge`.

Results beyond the double range in `wave_beam`.

Malformed or unsafe input:
- JSON number literals past the double range;
- an SPD whose step is ≤ 0;
- a null name or endpoint;
- an endpoint or key containing a control byte, which closes a CR-LF injection.

Each refusal uses the function's existing channel: NaN from a scalar, 0 from a constructor, and
a `PK_ERR_*` code through `err_out`.

**Behaviour you may have relied on:**
- `sellmeier_n_at` no longer clamps n² < 1 to 1. Below 0 it is NaN, and between 0 and 1 you get
  the model's value.
- `ray_snell` and `fresnel_unpolarized` at the exact critical angle answer grazing (R = 1), not
  TIR.
- `lens_type_to_json` returns 0 for a value that is not one of the two lens types.
- `point_source_new` returns a one-element array with a length header, so it can be passed to
  `interference_pattern`.

**A hand-assembled build** that includes `src/lens.cyr` without `src/wave_zernike.cyr` now
fails to link `lens_mtf_through_focus`, as it already did for the `lens_otf_*` family. The
bundles carry both files.

## 2.15.0: fluorescence

New API only: nothing existing changes value or signature, and the cyrius pin stays at
6.7.18. `dist/prakash.cyr` gains `spectral_fluor`. A name of yours that now collides:
- one beginning `fluor_`, `fret_` or `FLUOR_`;
- `Donaldson`;
- `spd_energy_to_photons` or `spd_photons_to_energy`.

Its spectra are `Spd`s of spectral power, except where a function asks for photons.

## 2.14.0: fifth-order aberrations

New API only: nothing existing changes value or signature, and the cyrius pin stays at
6.7.18. `dist/prakash.cyr` gains `lens_aberration`, so a name of yours beginning `aberr_`,
`Aberr` or `ABERR` would now collide. Its third-order coefficients are a `Prescription`'s
surface sums, which `lens_seidel_coefficients` (one thin lens, stop at the lens) never was;
that function is unchanged. A coefficient here is the transverse aberration in lens units at
full aperture and field for the rays you pass, not a per-unit-height form.

## 2.13.0: vectorial focusing

New API only: nothing existing changes value or signature, and the cyrius pin stays at
6.7.18. `dist/prakash.cyr` gains `wave_focus`, so a name of yours beginning `focus_`,
`Focus` or `FOCUS_` would now collide. Its fields follow e^(-iωt), like `ray_prt`, and are
the complex conjugates of fields written in e^(+iωt), which `wave_beam` uses.
`bessel_j1` in `wave_airy` is unchanged and still a 1e-7 fit; the new module carries its
own J₀, J₁ and J₂.

## 2.12.0: polarization ray tracing

New API only: nothing existing changes value or signature, and the cyrius pin stays at
6.7.18. `dist/prakash.cyr` gains `ray_prt`, so a name of yours beginning `prt_` or `Prt`, or
one called `fresnel_amplitudes`, would now collide. A complex phase from the new API follows
e^(-iωt), the conjugate of a published matrix written in e^(+iωt), which Chipman's group's is.

## 2.11.1: cyrius 6.7.18 and stricter JSON

- **No optics result moves.** About 3 million samples are bit-identical between 6.7.9
  and 6.7.18, and JSON output is byte-identical.
- ⚠ **`*_from_json` refuses two malformed inputs.** A duplicate key used to take its first
  value and a leading zero (`01`) used to read as its digits. Both are `PK_ERR_PARSE`
  now. If a stored document starts failing, fix the producer. prakash's own encoder
  writes neither.
- **A hand-vendored `lib/` needs more files.** bayan now reaches `chrono`, `slice`,
  `hashseed` and `lib/unicode/`. `cyrius deps` vendors them for you.
- **The cyrius ≥ 6.6.13 floor** for the JSON surface is unchanged.

## 2.11.0: cyrius 6.7.9 and long optical paths

- **The bump itself moves nothing.** About 3 million samples are bit-identical between
  6.6.20 and 6.7.9, so a consumer can take either pin for the same numbers.
- **`interference_pattern` and `fabry_perot_transmittance` move, within rounding:** at
  most 3.0e-8 and 6.3e-8 relative, and their worst phase error is smaller. A pin on an
  old value of either may need its tolerance checked.
- **`fabry_perot_transmittance` is corrected at a high-finesse resonance:** it returns 1
  exactly at a whole order, where R = 1 − 2⁻³⁰ at 2²² orders used to give 0.451.
- **The 9× cost below is gone** from both functions, which now run ~9× faster on long
  paths. Other calls with a huge sin/cos argument still pay it.
- **hisab 3.4.0 retires 27 public names.** prakash calls none of them, but a consumer
  that calls hisab directly should read hisab's `docs/guides/migration-3.4.md`.

## 2.10.1: what the cyrius 6.6.20 bump moves

A consumer compiles the bundle against its **own** stdlib, so none of this reaches it
until its own pin moves. Measured by running one program under 6.6.6 and 6.6.20 and
comparing bits.

- **What ranga calls is bit-identical** over the sample: `cct_from_xy`, `cie_cmf_at`,
  `spd_to_xyz` on blackbodies. So is `wavelength_to_rgb`. On its re-exported surface,
  `color_temperature_to_rgb`'s red channel and `color_rendering_index` move by at most
  2 ulp.
- ⚠ **A zero can come back as −0.** `f64_neg(+0)` is −0 from cyrius 6.6.8, and
  `x == 0` is a bit compare, so `doe_lens_phase` on the axis or `fiber_dispersion` in
  vacuum now read false against it. Use `f64_eq(x, 0)` or `(x & 0x7FFFFFFFFFFFFFFF) == 0`.
- **`sin`/`cos` are right near multiples of π.** x87 `fsin` was up to 44% wrong there.
  Resonant results move: Fabry-Perot by up to 2.2e-11, slits and multilayers by up to
  6e-14. The new values are the correct ones.
- **FFT patterns** are bit-identical to n = 128, and move by at most 3.4e-14 of peak at
  n = 1024.
- **`num_modes`** saturates at `i64::MAX` past V ≈ 2³², where it used to answer 1.
- **JSON** round-trips bit-exactly, so a pin on an old decoded value goes red by design.
- **Cost:** `sin`/`cos` at |x| ≥ 2²⁰·π/2 ≈ 1.65e6 are ~10× slower (525 ns, full
  reduction), and `pow` is ~2× slower (198 ns). `interference_pattern` is the row that
  shows it, at about 9×; a phase k·r passes that threshold at r ≈ 0.14 m in visible light.
- ⚠ **New floor: bayan ≥ 1.5.10, which is cyrius ≥ 6.6.13, for the JSON surface.** Below
  it the bundle still links and runs everything else, with a second undefined-function
  warning, `bayan_json_v_obj_get_by_cstr`, beside `num_fft`. ranga's cyrius 6.6.9 is below
  it: its surface passes, and the FFT-refusal check reports the extra name.
- **`dist/prakash.cyr` opens with `include "lib/bayan.cyr"`.** It is include-once, so a
  consumer that already declares bayan is unaffected, but it does not make the bundle
  self-sufficient: declare `dist/prakash.deps`.

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
has one, instead of exiting 139. ⚠ From 2.15.1 the 41 f64 functions without an `err_out`
answer NaN rather than 0.0; see the 2.15.1 section.

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
