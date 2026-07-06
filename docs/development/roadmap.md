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

### Consumer integration (as the Rust consumers port)

- [ ] soorat / kiran / ranga: consume `dist/prakash.cyr` directly once they move to Cyrius
- [ ] hisab geometry bridge: adopt `hisab` `HVec3`/ray types in the tracer (currently local `RayVec3`)
- [ ] hisab numerical integration: use `hisab` Gauss-Legendre for Huygens-Fresnel / SPD generation

### Accuracy & completeness (demand-gated)

- [ ] Gradient-index (GRIN) optics: curved ray paths through variable-n media
- [ ] Diffractive optical elements (DOE): phase gratings, holographic elements
- [ ] Vectorial diffraction (Richards-Wolf): high-NA focusing beyond scalar theory
- [ ] Hermite-Gaussian / Laguerre-Gaussian beam modes; M² beam quality
- [ ] Higher-order (5th-order Buchdahl) aberrations; wavefront coefficients (W_040, …) from Seidel sums
- [ ] Aberrated MTF from generalized pupil-function autocorrelation

### Advanced / demand-gated

- [ ] Fluorescence (Stokes shift, excitation/emission spectra)
- [ ] Non-linear optics (SHG, Kerr) — if joshua needs it
- [ ] Orbital angular momentum (Laguerre-Gaussian modes)
- [ ] Metamaterials / negative refractive index
- [ ] Age-dependent CIE observer (CIE 2006)

### Housekeeping

- [ ] Remove `rust-old/` a couple releases out (retained now as the translation reference; also in pre-2.0 git tags)

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
