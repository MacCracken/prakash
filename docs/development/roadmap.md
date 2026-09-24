# Prakash Roadmap

> **Prakash** is the optics/light-simulation library, written in [Cyrius](https://github.com/MacCracken/cyrius). Its one external math call is hisab's FFT (`num_fft`). Consumers are listed in `docs/architecture/overview.md`.

**Forward-facing only.** This lists what is not done yet. Shipped work lives in
`CHANGELOG.md`; a finished item is deleted from here, not checked off.

## Scope

Prakash owns the **physics of light**: how light travels, bends, reflects, interferes, diffracts, polarizes, and interacts with materials. It provides the math; consumers decide what to do with it.

| Not prakash's | Owner |
|---|---|
| Rendering pipeline | kiran / soorat |
| Image processing, ICC profiles, gamut mapping | ranga |
| Vectors, geometry, calculus, FFT | hisab (prakash calls only `num_fft`) |
| 3D scene graph | kiran |

## Release plan

Each open item is pinned to the minor release that ships it. Order is build order,
but nothing here depends on anything else here, so any release can be swapped.
Fixes found along the way ship as patches of the current minor.

| Version | Release | Pinned against |
|---|---|---|
| **2.11.0** | 3×3 polarization ray tracing | Chipman's PRT matrix |
| 2.12.0 | Vectorial focusing | Richards-Wolf, scalar Airy limit |
| 2.13.0 | Fifth-order aberrations | Buchdahl, real-ray residual |
| 2.14.0 | Fluorescence | Stokes shift, Donaldson matrix |
| 2.15.0 | Nonlinear optics | SHG, Kerr |
| 2.16.0 | Negative-index media | Veselago-Pendry lens |
| 2.17.0 | CIE 2006 observer | CIE 170-1:2006 tables |

Keep each section below to a few lines: what ships and what it is pinned against.

### 2.11.0 — 3×3 polarization ray tracing
Chipman's polarization ray-tracing matrix per surface, accumulated through the
sequential tracer with the complex Fresnel coefficients. Today
`trace_sequential_polarized` tracks scalar s/p transmittance from the real
coefficients only.

### 2.12.0 — Vectorial focusing (Richards-Wolf)
The I₀, I₁, I₂ integrals for high-NA focusing of linear, radial and azimuthal
polarization, with a `GaussianBeam` pupil fill. It must reproduce the scalar Airy
pattern at low NA.

### 2.13.0 — Fifth-order aberrations (Buchdahl)
Buchdahl's fifth-order sums for a surface sequence. They are pinned against the
real-ray residual left after the third-order wavefront coefficients are subtracted.

### 2.14.0 — Fluorescence
Excitation and emission spectra as `Spd`s, Stokes shift, quantum yield, the
re-radiation (Donaldson) matrix for spectral rendering, and FRET efficiency.

### 2.15.0 — Nonlinear optics
- **SHG:** phase mismatch sinc²(ΔkL/2), coherence length, undepleted-pump
  efficiency.
- **Kerr:** n = n₀ + n₂I, self-phase modulation, B-integral, critical power for
  self-focusing.

### 2.16.0 — Negative-index media
Snell and Fresnel for ε, μ < 0, Drude-Lorentz ε(ω) and μ(ω), and the
Veselago-Pendry flat lens.

### 2.17.0 — CIE 2006 physiological observer
Age- and field-size-dependent cone fundamentals: lens and macular pigment optical
densities, and LMS → XYZ. Pinned to the published 2° and 10° tables at age 32.

## Constraints — read before optimizing or testing

Each of these rules is a measurement. The CHANGELOG has the story behind each one.

**Performance**
- **Cost model (cyrius 6.6.6):**
  - `f64_mul` (~2 ns) and `f64_div` (~4.5 ns) are inline; a user-function call adds
    1–3 ns.
  - The expensive calls are transcendentals: ln ~43 ns, atan2 ~53 ns, exp ~20 ns.
  - So remove transcendental and function calls, not arithmetic, and hoist
    loop-invariant transcendentals first.
- **Benchmarking:** run-to-run spread reaches 40%. Prove any change same-binary,
  with both copies in one bench binary plus an A/A pair.
- **Numerically stable rewrites:** the first stable form is not always the cheapest
  one. Benchmark the alternatives.
- **Hand-inlining:** only at measurably hot sites. A constant helper returns a
  literal, pinned bit for bit in `tests/constants.tcyr`.
- **Small allocations:** fixed-size scratch `alloc(16)` costs ~7 ns. It is worth
  converting only on rows of a few hundred ns.
- **Caller buffers:** a caller buffer needs a caller that owns the lifetime.
  `trace_sequential` retains its hits.
- **Dead ends:**
  - SIMD is exhausted: only raw `f64v_*` bulk calls help, and `f64v_scale` overruns
    on odd counts.
  - `hvec3_*` in internals costs 2–11%.
  - hisab's quadrature fits no prakash path.
- **JSON cost:** bayan's cost is string bytes, appended one at a time. Use
  `str_builder_add_json_str`, which escapes byte-identically.

**Data and API**
- **`SellmeierCoefficients`** holds three terms. Measure a least-squares merge
  before widening it; water's four-term fit merges to 3.5e-6 in n.
- **An `Spd`** holds point samples, not bin integrals. Integrating while building
  one makes photometry 5× wrong while chromaticity still looks right.

**Testing**
- **Pin against an independent implementation or reference**, never a value derived
  from the code. A shape or ratio check cannot see a constant factor.
- **Pin relationships:**
  - a derived constant against its derivation;
  - two functions naming one quantity against each other;
  - two models of one material at their true agreement.
- **Every Sellmeier preset** carries an Abbe pin. A fidelity repair below the
  physics tolerance needs a pin on the literal.
- **Before repairing a formula,** grep for its other evaluation sites.
- **Distributions:** the integral is the strongest check. Range assertions
  supplement value pins; they never replace them.
- **Numeric windows and tolerances:** measure them across the whole domain,
  including both edges.
- **Defect classes** (null handles, zero counts) are swept mechanically, one build
  per public function, not by reading.
- **Upstream defects** are pinned at their defective value, so the assertion fails
  when upstream fixes them. Example: bayan's 1-ULP decoder mis-round.
- **Counts and measurements in prose drift:** re-measure them whenever a release
  touches the code or the tests.
- **Test the source, not the bundle.** Bundle inclusion only answers "does it link"
  (`tests/ai_bundle.tcyr`).
