# Prakash → Cyrius Port Plan

Porting prakash from Rust (1.2.0) to the **Cyrius** systems language. The port
ships as **prakash 2.0.0** — mirroring how [`hisab`](https://github.com/MacCracken/hisab)
bumped 1.4.0 → 2.0.0 for its own Rust→Cyrius port. `hisab` is the canonical
reference template for this work.

Status: **planning complete, scaffolding underway.**

---

## 1. What the port is (and isn't)

Cyrius is a sovereign, C-level systems language. This is not a syntax reskin —
it is a re-expression of ~15,100 lines of Rust into a language with:

- **No native float.** `f64` is an IEEE-754 bit pattern held in an `i64`. Every
  arithmetic op is a function/builtin call: `a*b + c` → `f64_add(f64_mul(a,b), c)`.
- **No methods.** Free functions only: `medium_permittivity(m)`, not `m.permittivity()`.
- **Heap structs.** `alloc(sizeof(T))` + `#derive(accessors)` → `T_field(v)` / `T_set_field(v, x)`.
- **No `Result`/`Option` enums for our error model.** Integer error codes
  (`ERR_NONE = 0`, negative `ERR_*`); out-values via pointer args.
- **No generics, no traits, no closures-with-capture, no borrow checker.**
- **Manual memory** via a bump allocator (`alloc_init()` once at program start).

Because ~90% of prakash is dense f64 math with no `unsafe`, no SIMD, and almost
no trait/generic use, the dominant task is **mechanical arithmetic translation**
guarded by the ported test suite (600+ cases → `.tcyr` assertions).

---

## 2. Toolchain, layout, and file kinds

| Rust | Cyrius |
| --- | --- |
| `cargo build/test/bench` | `cyrius build` / `cyrius test` / `cyrius bench` |
| `cargo fmt` / `clippy` | `cyrius fmt` / `cyrius lint` / `cyrius vet` / `cyrius audit` |
| `Cargo.toml` | `cyrius.cyml` (`VERSION` file is SSOT via `${file:VERSION}`) |
| `Cargo.lock` | `cyrius.lock` |
| `src/lib.rs` + `pub mod` | flat `src/*.cyr` modules, bundled by `cyrius distlib` → `dist/prakash.cyr` |
| `src/foo.rs` inline `#[test]` | `tests/*.tcyr` (`test_group` + `assert_*`) |
| `benches/*.rs` (criterion) | `tests/*.bcyr` (`bench_new`/`bench_run`/`bench_report`) |
| — | `tests/*.fcyr` (fuzz targets) |

`cyrius distlib` strips `include` lines and concatenates the `[lib] modules` list
in dependency order into `dist/prakash.cyr`. Consumers depend on that bundle (or
symlink individual `lib/prakash_*.cyr` files). CI gets a **distlib-drift gate**
(rebuild must be a no-op), matching hisab.

---

## 3. Dependencies (`cyrius.cyml`)

```toml
[deps]
stdlib = [
    "syscalls", "string", "alloc", "str", "fmt", "vec", "io", "args",
    "assert", "math", "ganita", "tagged", "fnptr", "bench", "callback",
    "bayan",           # JSON — replaces serde
    "sandhi",          # HTTP client — replaces reqwest/tokio (ai module only)
]

[deps.hisab]
git = "https://github.com/MacCracken/hisab.git"
tag = "2.6.7"          # Complex + FFT for wave/pattern
modules = ["dist/hisab.cyr"]

[deps.sakshi]
git = "https://github.com/MacCracken/sakshi.git"
tag = "2.4.4"          # structured logging — replaces tracing
modules = ["dist/sakshi.cyr"]
```

**Note on `ai` + `sandhi`:** `sandhi` drags in TLS. To avoid forcing every
math-only consumer to compile it, `ai.cyr` is kept **out of the core
`dist/prakash.cyr` bundle** and shipped as a separate opt-in module — the Cyrius
equivalent of the old `ai` cargo feature. Consumers who want it include
`ai.cyr` explicitly and add `sandhi` to their own stdlib deps.

---

## 4. The capability map (verified — nothing is dropped)

| prakash need | Cyrius source | Symbols |
| --- | --- | --- |
| `+ − × ÷ neg`, `sqrt`, `sin`/`cos`, `exp`/`ln`, `floor`/`ceil`, `abs`, `f64_lt/gt/eq` | compiler f64 builtins (aarch64 polyfills in stdlib `math`) | `f64_add`, `f64_sqrt`, `f64_sin`, … |
| `f64_le/ge`, `min`/`max`/`clamp`/`lerp`/`hypot`/`round`/`trunc`/`fract`/`sign`, `F64_PI` & friends, `f64_from`/`f64_to`/`f64_parse` | stdlib `math` | `f64_min`, `f64_clamp`, `F64_PI`, … |
| `asin`/`acos`/`atan2`/`pow`/`sinh`/`cosh`/`tanh`/`asinh`/`acosh`/`atanh`/`hypot` | **ganita** `math_advanced` | `ganita_f64_asin`, `ganita_f64_atan2`, `ganita_f64_pow`, … |
| matrix, **SVD**, **least-squares**, eigen(sym), LU, QR, inv, rank | **ganita** `linalg` | `ganita_mat_svd`, `ganita_mat_least_squares`, `ganita_mat_eigen_sym`, … |
| Complex numbers, FFT/IFFT | **hisab** | `cx_new`/`cx_mul`/`cx_add`/`cx_norm_sq`/`cx_exp`, `num_fft`/`num_ifft` |
| JSON encode/decode (was serde) | **bayan** | `bayan_json_v_obj_new`/`_obj_set`/`_str_new`/`_int_new`/`_float_new`/`_arr_push`, `_v_build`(`_pretty`), `_parse`, `_v_obj_get`/`_v_str`/`_v_float` |
| HTTP POST + timeout + JSON (was reqwest) | **sandhi** | `sandhi_http_post`, `sandhi_http_options_total_ms`, `sandhi_http_status`/`_body`, `sandhi_json_obj_new`/`_add_string`/`_add_raw`/`_build`, `sandhi_json_get_string` |
| structured logging (was tracing) | **sakshi** | `sakshi_info`/`_warn`/`_error`/`_debug`/`_trace`, `sakshi_span_enter`/`_exit`, `sakshi_log_kv` |

`atan` (unused by hisab) and `powf` are covered by `ganita_f64_atan2`/`ganita_f64_pow`;
any missing thin wrapper we add to a local `constants.cyr`/helper as needed.

---

## 5. Idiom translation cheat-sheet

| Rust | Cyrius |
| --- | --- |
| `struct Rgb { r: f64, g: f64, b: f64 }` | `#derive(accessors) struct Rgb { r; g; b; }` + `rgb_new(r,g,b)` allocating `sizeof(Rgb)` |
| `impl Rgb { fn luminance(&self) -> f64 }` | `fn rgb_luminance(c): i64` (returns f64 bits) |
| `Result<f64, PrakashError>` | return f64 bits on success; signal failure via error code out-param or a sentinel + `err` code (see §6) |
| `Option<T>` | `tagged.cyr` (`Tag_Some`/`Tag_None`) or sentinel |
| `enum SurfaceShape { Sphere{radius}, Plane }` | tag constant + payload fields in a heap struct, or `enum` tag + branch |
| `Vec<TraceSegment>` | stdlib `vec` (holds `i64`/pointers) |
| `[f64; 3]` direction/normal | **hisab `HVec3`** (reuse tested `hvec3_dot`/`_cross`/`_normalize`/`reflect`) or a raw 3×f64 alloc |
| `#[derive(Serialize)]` | hand-written `*_to_json` / `*_from_json` via bayan |
| `tracing::info!("x", k=v)` | `sakshi_info("x", 1)` / `sakshi_log_kv(SK_INFO, msg, mlen, k, klen, v, vlen)` |
| doc `///` / `//!` | leading `#` comment block per file/function |
| `const EPS: f64 = 1e-6` | `var EPS = 0x3EB0C6F7A0B5ED8D;  # 1e-6` (bit pattern, with decimal comment) |
| negative literal `-1` | `(0 - 1)` |
| private `fn helper` | `fn _helper` (underscore-prefixed) |

**Vectors:** prefer hisab `HVec2`/`HVec3` for `[f64;2]`/`[f64;3]` so we reuse
tested dot/cross/normalize/reflect and stay consistent with the ecosystem.

---

## 6. Error model

Port `PrakashError`'s 8 variants to negative integer codes in `error.cyr`,
following hisab's `error.cyr`. Detailed message formatting (angles, indices) is
dropped from the return path — callers that need a message log it via sakshi at
the failure site. Proposed codes:

```
ERR_NONE                    =  0
PRK_ERR_TIR                 = (0 - 1)   # total internal reflection
PRK_ERR_INVALID_INDEX       = (0 - 2)   # refractive index < 1.0
PRK_ERR_WAVELENGTH_RANGE    = (0 - 3)   # outside 380–780 nm
PRK_ERR_INVALID_ANGLE       = (0 - 4)   # outside 0–90°
PRK_ERR_INVALID_FOCAL       = (0 - 5)   # negative focal length
PRK_ERR_DIVISION_BY_ZERO    = (0 - 6)
PRK_ERR_INVALID_PARAMETER   = (0 - 7)
```

Functions that returned `Result<f64>` return the f64 bits and take an
`err_out` pointer (write `ERR_NONE` or a code), or — where a NaN sentinel is
unambiguous — return NaN and expose a separate validating variant. Decide
per-function; default to the `err_out` pointer for clarity. `prakash_is_err(code)`
mirrors `hisab_is_err`.

---

## 7. Module map & sequencing

Cyrius `src/` is flat (no subdirs), so Rust's `ray/mod.rs` etc. become prefixed
files. Ship order is dependency order; each module is ported → `fmt/lint/vet`
clean → tests ported → `cyrius test` green → benched → next.

| # | Cyrius module | From (Rust) | LOC | fn prefix | Notes / deps |
| --- | --- | --- | --- | --- | --- |
| 0 | `error.cyr` | `error.rs` | 157 | `PRK_ERR_*` | integer codes; first, needed everywhere |
| 0 | `constants.cyr` | scattered consts | — | `PRK_*` | physical consts (c, h, k_B), visible range, shared EPS |
| 1 | `ray_core.cyr` | `ray/mod.rs` | 988 | `ray_`, `medium_` | Snell, Fresnel-basic, reflect/refract 2d/3d, Medium DB |
| 1 | `ray_fresnel.cyr` | `ray/fresnel.rs` | 271 | `fresnel_` | ComplexMedium, s/p/unpolarized, Brewster, Beer-Lambert |
| 1 | `ray_trace.cyr` | `ray/trace.rs` | 531 | `trace_` | surfaces, sequential trace |
| 1 | `ray_simulate.cyr` | `ray/simulate.rs` | 999 | `sim_` | recursive trace, ray fans, spot diagram, OPD (`vec`) |
| 1 | `ray_system.cyr` | `ray/system.rs` | 723 | `sys_` | multi-element, Abbe, dispersion trend |
| 1 | `ray_dispersion.cyr` | `ray/dispersion.rs` | 832 | `disp_` | Cauchy, Sellmeier, Helmholtz |
| 1 | `ray_fiber.cyr` | `ray/fiber.rs` | 174 | `fiber_` | NA, V-number, modes, coupling |
| 2 | `spectral_core.cyr` | `spectral/mod.rs` | 562 | `spec_`, `rgb_` | wavelength↔RGB, Planck, color temp |
| 2 | `spectral_cie.cyr` | `spectral/cie.rs` | 1350 | `cie_` | CIE tables, XYZ/Lab, ΔE, CRI, CCT (largest module) |
| 2 | `spectral_photometry.cyr` | `spectral/photometry.rs` | 198 | `photo_` | V(λ), lumens |
| 3 | `wave_core.cyr` | `wave/mod.rs` | 892 | `wave_` | interference, diffraction, Jones/Stokes, Malus |
| 3 | `wave_polarization.cyr` | `wave/polarization.rs` | 695 | `pol_` | Mueller/Jones matrices (ganita mat) |
| 3 | `wave_coherence.cyr` | `wave/coherence.rs` | 61 | `coh_` | includes a Bessel J1 approx |
| 3 | `wave_airy.cyr` | `wave/airy.rs` | 82 | `airy_` | Airy disk, Rayleigh criterion |
| 3 | `wave_fabry_perot.cyr` | `wave/fabry_perot.rs` | 89 | `fp_` | finesse, FSR, transmittance |
| 3 | `wave_diffraction.cyr` | `wave/diffraction.rs` | 841 | `diff_` | gratings, Bragg, Raman-Nath |
| 3 | `wave_zernike.cyr` | `wave/zernike.rs` | 678 | `zern_` | Zernike; fit via `ganita_mat_least_squares` |
| 3 | `wave_pattern.cyr` | `wave/pattern.rs` | 748 | `patt_` | **uses hisab `cx_*` + `num_fft`/`num_ifft`** |
| 4 | `lens.cyr` | `lens.rs` | 1185 | `lens_` | thin/thick, Seidel, MTF, DoF |
| 4 | `pbr_core.cyr` | `pbr/mod.rs` | 457 | `pbr_` | Cook-Torrance, GGX, Smith, Schlick |
| 4 | `pbr_advanced.cyr` | `pbr/advanced.rs` | 1179 | `pbr_` | aniso, sheen, clearcoat, SSS, iridescence |
| 4 | `atmosphere.cyr` | `atmosphere.rs` | 827 | `atm_` | Rayleigh/Mie, sky color, air mass |
| 5 | `bridge.cyr` | `bridge.rs` | 378 | `br_` | bijli/tara/badal converters (primitive f64 in/out) |
| 6 | `ai.cyr` | `ai.rs` | 135 | `ai_` | **out of core bundle**; sandhi POST + bayan JSON |

`logging.rs` has no dedicated module — internal `tracing` calls become `sakshi_*`
calls inline; consumers call `sakshi_set_level`. A thin `ai_*`/global init helper
may wrap `sakshi` setup if convenient.

---

## 8. Serialization (bayan) — replacing serde

Types that were `#[derive(Serialize, Deserialize)]` — `Medium`, `TraceRay`,
`TraceConfig`, `Polarization`, `OpticalSurface`, `DaimonConfig`, `HooshConfig`,
etc. — get hand-written pairs:

```
fn medium_to_json(m): i64        # returns bayan json value (obj)
fn medium_from_json(v): i64      # returns Medium ptr (or err via out-param)
```

Build with `bayan_json_v_obj_new` + `_obj_set`/`_str_new`/`_float_new`; read with
`bayan_json_v_obj_get` + `_v_str`/`_v_float`. Only add these where a Rust test
exercised serde or a consumer needs wire format — not blanket. The serde
round-trip unit tests become bayan build→parse→compare `.tcyr` assertions.

---

## 9. The `ai` module (sandhi) — replacing reqwest/tokio

`register_agent` (async in Rust) becomes a **blocking** sandhi call:

1. `ai_daimon_config_new(endpoint, api_key)` → heap struct.
2. Build body: `sandhi_json_obj_new()` → `sandhi_json_add_string(o,"name","prakash")`
   → `sandhi_json_add_raw(o,"capabilities","[\"optics\",\"ray_tracing\",…]")`
   → `body = sandhi_json_build(o)`.
3. `opts = sandhi_http_options_new(); sandhi_http_options_total_ms(opts, 30000);`
4. `resp = sandhi_http_post(url, headers, body, body_len /*, opts */)`.
5. Guard `sandhi_http_status(resp)` in 200–299; else error + `sakshi_error`.
6. `agent_id = sandhi_json_get_string(sandhi_http_body(resp), "agent_id")`.

Config round-trip tests port to bayan (§8). Kept out of the core distlib bundle.

---

## 10. Testing, benchmarks, fuzz

- **Tests** (`tests/*.tcyr`): split like hisab — e.g. `prakash.tcyr` (smoke),
  `ray.tcyr`, `wave.tcyr`, `spectral.tcyr`, `optics_edge.tcyr`. Each `include`s the
  modules it exercises, calls `alloc_init()`, uses `test_group` + `assert_*`.
  Float asserts compare `f64_to(...)` or an `f64_approx_eq(a, b, eps)` helper
  (prakash uses `EPS = 1e-6`). Target ≥80% coverage, port all 600+ Rust cases.
- **Benchmarks** (`tests/prakash.bcyr`): mirror the criterion set (ray-trace
  throughput, FFT diffraction, Zernike fit, Cook-Torrance). Use amplification
  (x64) for sub-timer-floor ops. Record to `bench-history.csv`.
- **Fuzz** (`tests/prakash.fcyr`): invariants — Snell round-trips, Fresnel ∈ [0,1],
  energy conservation ≤ 1, RGB clamps, no crash on extreme inputs.

---

## 11. Docs & release

- Update `README.md` (Cyrius quick-start, build/test/bench, consumers),
  `CHANGELOG.md` (`[2.0.0]` Breaking: language change — Rust→Cyrius, feature
  gates → module includes, serde→bayan, tracing→sakshi, reqwest→sandhi),
  `docs/architecture/overview.md`, `docs/architecture/math.md`, roadmap.
- `CLAUDE.md` process stays; add Cyrius toolchain commands.
- Bump `VERSION` → `2.0.0`; keep `${file:VERSION}` as SSOT in `cyrius.cyml`.
- Archive the Rust source in a pre-2.0 git tag (user handles git); repo goes
  Cyrius-native like hisab.
- CI: version-verify gate, `fmt`/`lint`/`vet`, `deps --verify` (lock), distlib
  drift gate, test/fuzz/bench gates.

---

## 12. Risks & open items

1. **Volume.** ~15k LOC of arithmetic to translate by hand; the test suite is the
   safety net — port tests alongside each module, never after.
2. **`[f64;N]` ↔ HVec.** Decide HVec vs raw alloc per function (default HVec).
3. **Error ergonomics.** `err_out` pointer vs NaN sentinel — pick per function; be
   consistent within a module.
4. **CIE / V(λ) / atmosphere constant tables.** Large `[(f64,f64,...)]` literals →
   flat f64-bit arrays; generate with a helper rather than hand-typing bit patterns.
5. **bayan/sandhi API drift.** Pin exact tags; the versions here are today's.
6. **Numerical parity.** f32→f64 already done in spirit (Rust is f64); watch
   transcendental impl differences (ganita vs libm) at tight tolerances.
7. **Consumers still Rust** (soorat, ranga, kiran) — no Cyrius consumer to satisfy
   yet, so we follow hisab idioms freely; their ports come later.

---

## 13. Milestone checklist

- [x] **M0 Scaffold** — `cyrius.cyml`, `src/main.cyr` entry, `src/error.cyr`,
      `tests/prakash.tcyr` smoke. `cyrius build` OK, `cyrius test` 12/0, `lint`
      clean, `cyrfmt --check` clean, `cyrius distlib` → `dist/prakash.cyr` OK.
      (`VERSION`→2.0.0 deferred to M8 per "version 2.0.0 when port completed".
      `constants.cyr` lands with the first module that needs it — spectral.)
- [x] **M1 Ray** — 7 modules DONE & verified: ray_core (96), ray_fresnel (175),
      ray_trace (28), ray_simulate (180), ray_system (46), ray_dispersion (54),
      ray_fiber (14) = **593 test assertions**, all lint/fmt clean, each with a
      passing adversarial-verify workflow. (benches deferred to a later pass.)
- [ ] **M2 Spectral** — 3 modules (incl. CIE tables) + tests.
- [ ] **M3 Wave** — 8 modules (pattern last, on hisab FFT) + tests.
- [ ] **M4 Lens / PBR / Atmosphere** — + tests.
- [ ] **M5 Bridge** — + tests.
- [ ] **M6 AI** — sandhi port (out-of-bundle) + tests.
- [ ] **M7 Serialization** — bayan `*_to_json`/`*_from_json` where needed.
- [ ] **M8 Release** — distlib, docs, CHANGELOG, benchmarks, CI gates, 2.0.0.
