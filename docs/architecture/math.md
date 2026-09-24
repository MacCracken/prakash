# Mathematical Reference

Quick reference for the core ray, wave, spectral, lens, PBR and atmosphere
formulas in prakash.

⚠ **This file is a subset, not a complete index, and 2.2.8 narrowed the claim
rather than pretending otherwise.** It previously said "all physics formulas".
Not covered here: CIE 13.3 colour rendering (R_i = 100 - 4.6*dE_i in CIE 1964
W*U*V* after the von Kries c/d adaptation, with the reference illuminant chosen
Planckian below 5000 K and D-series above), the thin-film transfer matrix, the
Zernike/Strehl family, Malus/Jones/Mueller polarisation, the Conrady and
Herzberger dispersion models, fibre NA / V-number / mode-field diameter, and
Huygens-Fresnel propagation. For those, the source comments are the reference —
they carry the derivations and the measured validation tables.

## Ray Optics

### Snell's Law
n₁ sin θ₁ = n₂ sin θ₂

### Fresnel Equations
Rs = ((n₁ cos θᵢ - n₂ cos θₜ) / (n₁ cos θᵢ + n₂ cos θₜ))²
Rp = ((n₂ cos θᵢ - n₁ cos θₜ) / (n₂ cos θᵢ + n₁ cos θₜ))²

### Fresnel at Normal Incidence
R = ((n₁ - n₂) / (n₁ + n₂))²

### Brewster's Angle
θ_B = arctan(n₂ / n₁)

### Critical Angle
θ_c = arcsin(n₂ / n₁)

### Beer-Lambert
I = I₀ exp(-αd)

### Cauchy Dispersion
n(λ) = B + C/λ²

### Sellmeier Dispersion
n²(λ) = 1 + Σ Bᵢλ² / (λ² - Cᵢ)

### Abbe Number
V = (n_d - 1) / (n_F - n_C)

### Prism Deviation
δ = 2 arcsin(n sin(α/2)) - α

## GRIN Ray Tracing (`ray_grin`, 2.8.0)

### The ray equation (Sharma, Kumar & Ghatak 1982)
d/ds(n dr/ds) = ∇n, rewritten with dt = ds/n and the optical direction T = n dr/ds:

d²R/dt² = D(R),  D = n∇n = ½∇(n²),  T = dR/dt,  |T| = n.

Integrated with their Runge-Kutta-Nyström scheme, fourth order, three D
evaluations per step (the end point's D is the next step's first):

A = Δt D(R),  B = Δt D(R + (Δt/2)T + (Δt/8)A),  C = Δt D(R + Δt T + (Δt/2)B)
R′ = R + Δt [T + (A + 2B)/6],  T′ = T + (A + 4B + C)/6

Optical path: OPL = ∫n ds = ∫n² dt, by the two-point Hermite rule
Δt/2 (f₀ + f₁) + Δt²/12 (f₀′ − f₁′) with f = n², f′ = 2 D·T. That is also fourth order and
costs nothing, because n² and D at both ends are already known.

**The step is an arc length, every step:** Δt = min(step/n(start), step/n(here)). The cap
stops Δt growing where n falls towards 0. The global error is O(step⁴), with measured
ratios of 16.3–16.9 per halving.

A step is re-taken at half the length, up to 40 times, when it does not behave like a
step of a ray:
- a stage point has no index;
- T turns by more than 0.5 rad;
- n² changes by more than 2×;
- the optical path does not grow.

A trace stops at once if |T| drifts 10% from n. It is refused if |T|²/n² is more than
2·10⁻³ from 1 at its end. That end check is **relative to the index at the end**:
drift made where n is large reads (n_max/n_end)² larger where the ray lands. So a trace
landing in much lower index than it crossed needs a proportionally smaller step, and
may be refused even when its result would have been accurate. Rescaling the check by
n_max instead let 0.29 rad direction errors through. A fish-eye trace is capped at an
optical path of 2πn₀R: every fish-eye ray is a closed circle of optical length πn₀R.
None of this is error control: it keeps a coarse step from returning garbage as success.

A trace ends exactly on its boundary (a plane z = z_end, or a sphere left from inside):
- the step before is clamped to 9/8 of the crossing that the step's own quadratic
  g + g′t + g″t²/2 predicts (curvature included, so a ray skimming a sphere's rim
  from inside is not slowed to a crawl);
- a crossing inside a step whose ends are both inside, such as a thin cap grazed on
  the way out, is found as a peak of the cubic Hermite interpolant of g, from g and
  dg/dt at the step ends;
- the crossing step is then re-taken with Δt = τ, and τ is found by safeguarded
  Newton on g(step(τ)) = 0.

The sphere's tolerances scale with R(R + |c|), the rounding of |X − c|² − R². A rod's
wall is a cylinder r = a, checked the same way: the largest r² within a step is the
peak of the Hermite cubic of r², with d(r²)/dt = 2(xT_x + yT_y).

### Profiles
| Profile | n | D = ½∇n² |
|---|---|---|
| SELFOC (NSG form) | n² = n₀²[1 − (gr)² + h₄(gr)⁴ + h₆(gr)⁶] | n₀²g²(−1 + 2h₄u + 3h₆u²)(x, y, 0), u = (gr)² |
| sech | n₀ sech(gr) | −n² g² (tanh(gr)/(gr)) (x, y, 0) |
| polynomial ("gradient 3") | n₀ + n_r2 r² + n_r4 r⁴ + n_r6 r⁶ + n_z1 z + n_z2 z² + n_z3 z³ | n∇n |
| Luneburg | n_s √(2 − ρ²/R²) | −n_s² (R − c)/R² |
| Maxwell fish-eye | n₀ / (1 + ρ²/R²) | −2n₀² (R − c) / (R² (1 + ρ²/R²)³) |

### Closed forms the tests pin
- **Parabolic-n² SELFOC** (h₄ = h₆ = 0): D is linear, so x and y are harmonic in t
  with ω = n₀g and T_z is constant. Every ray, meridional or skew, is exact.
  Pitch 2πβ/(n₀g) with β = T_z, and 2π/g paraxially.
- **sech:** meridional rays satisfy sinh(gy) = C sin(gz + φ). Every one is periodic
  with period 2π/g, so a collimated fan focuses exactly on axis at a quarter pitch.
  SELFOC with h₄ = 2/3, h₆ = −17/45 matches it through (gr)⁶.
- **Linear axial n = n₀ + az:** transverse T is conserved (Snell) and
  y(z) = (β/a)[acosh(n/β)]. OPL = [n√(n² − β²) + β² acosh(n/β)]/(2a) between the
  end indices.
- **Luneburg:** a harmonic oscillator with ω = n_s/R. A collimated beam, at any height
  and in any direction d, focuses at c + Rd and leaves along −P/R (P the entry point).
  OPL from the tangent plane is n_s R(1 + π/2) on every ray.
- **Maxwell fish-eye:** every ray from P meets again at −R²P/|P|² (relative to the
  centre), with OPL n₀πR/2. A point on ρ = R images to its antipode.

### SELFOC paraxial optics
Pitch 2π/g. A rod of length L between flat faces, index n_ext outside, on
(height, real slope):

[ cos gL, sin gL · n_ext/(n₀g) ; −sin gL · n₀g/n_ext, cos gL ]

So a quarter-pitch rod has A = D = 0 and focal length 1/(n₀g) in air.

## Wave Optics

### Interference
I = A₁² + A₂² + 2A₁A₂ cos(δ)

### Phase from Path Difference
δ = 2π Δ/λ

### Thin Film Reflectance
R = sin²(δ/2) where δ = 2π(2nt)/λ + π

### Single-Slit Diffraction
I(θ) = I₀ (sin β / β)² where β = πa sin θ / λ

### Double-Slit Diffraction
I(θ) = 4 I₀ (sin β / β)² cos²(γ) where γ = πd sin θ / λ

(the two-slit term peaks at 4I₀, matching the implementation and the Rust original)

### Diffraction Grating
d sin θ = mλ

### Airy Pattern
I(θ) = I₀ [2J₁(x)/x]² where x = πD sin θ / λ

### Rayleigh Criterion
θ_min = 1.22 λ/D

### Fabry-Perot Transmittance
T = 1 / (1 + F sin²(δ/2))

### Fabry-Perot Finesse
F = π√R / (1 - R)

### Fabry-Perot Free Spectral Range
FSR = c / (2nt)

### Coherence Length
L_c = λ² / Δλ

### Coherence Time
τ_c = λ² / (c Δλ)

### Fresnel Number
N_F = a² / (λz)

### Fresnel Integrals
C(x) = ∫₀ˣ cos(πt²/2) dt
S(x) = ∫₀ˣ sin(πt²/2) dt

### Anti-Reflection Ideal Index
n_coat = √(n₁ n₂)

### Quarter-Wave Thickness
t = λ / (4n)

### Mueller/Stokes Formalism
S' = M · S (4x4 matrix times 4-vector)

## Diffractive Optics (`wave_doe`, 2.9.0)

### Grating equation
n_out sin θ_m = n_in sin θ_in + mλ/Λ. An order with |sin θ_m| > 1 is evanescent,
reported as `PK_ERR_TIR`: for order 0 that is total internal reflection, which is how
`ray_snell` reports it. A grazing order (|sin| = 1) counts as propagating.

### Scalar efficiencies (thin elements, the Fourier coefficients of the profile)
| Element | η_m |
|---|---|
| Blazed (sawtooth, phase depth α waves) | sinc²(α − m) |
| N-level staircase (levels at 2πα j/N) | sinc²(m/N) · [sin π(α − m) / (N sin(π(α − m)/N))]² |
| Binary phase (depth δ waves over duty D) | η₀ = 1 − 4D(1 − D) sin²πδ where 4D(1 − D) ≤ ½, else cos²πδ + (1 − 2D)² sin²πδ;  η_m = [2 sin πδ · sin πmD / (πm)]² |
| Binary amplitude (open fraction D) | η₀ = D²,  η_m = [sin πmD / (πm)]² |
| Dammann (0/π, even, transitions x_k ∈ (0, ½)) | c₀ = (−1)^J + 4Σ(−1)^{k+1}x_k,  c_m = (2/πm) Σ(−1)^{k+1} sin 2πm x_k,  η = c² |

sinc x = sin πx/(πx). The off-design depth of a surface-relief blaze in air is
α = M(λ₀/λ)(n(λ) − 1)/(n(λ₀) − 1). At design, a staircase puts sinc²(1/N) in order 1:
40.5%, 81.1%, 95.0% and 98.7% for N = 2, 4, 8 and 16. The published Dammann designs
reach 66.4% (1×3) and 77.4% (1×5). sin πx is evaluated with x reduced exactly, so
integer x gives exactly 0. The numerator sin π(α − m) is formed as ±sin πα, and the
Dirichlet denominator carries α − m with its exact rounding error (TwoSum). That way a
weak grating (α ≈ 10⁻⁶) and a staircase near α = kN keep their relative accuracy.
m is reduced modulo N in integers first. α and m past 2⁵² are refused.

### Diffractive lens
φ(r) = −sign(f)(2π/λ₀)(√(f² + r²) − |f|) exactly (paraxially −πr²/(λ₀f)), formed
without cancellation. Zone j ends at r_j = √(2jλ₀|f| + (jλ₀)²).
f(λ) = f₀λ₀/(mλ). Abbe number V = λ_d/(λ_F − λ_C) = −3.4534 for the d, F and C lines.
Hybrid achromat: φ_r = P V_r/(V_r − V_d), φ_d = −P V_d/(V_r − V_d), about 5% diffractive
with a crown.

### Volume holograms (Kogelnik 1969)
β = 2πn/λ, K = 2π/Λ, c_R = cos θ, c_S = cos θ − (K/β) cos φ,
ϑ = K cos(φ − θ) − K²/(2β), κ = πn₁/λ. TM multiplies κ by |r̂·ŝ|, which is
|β − K cos(φ − θ)| / √(β² − 2βK cos(φ − θ) + K²) with σ = ρ − K. That is Kogelnik's
eqs. (78)–(90), equal to cos 2(θ_B − φ) at Bragg. Taking cos 2(θ − φ) at the incidence
angle instead is 0.4% of η off at 4 mrad of detuning.

K and −K are the same fringes. The order nearer Bragg is taken: the sign of K that
makes K cos(φ − θ) ≥ 0. So slant φ + π and the mirror incidence −θ_B behave as the
physics says.
- Transmission (c_S > 0): ν = κd/√(c_R c_S), ξ = ϑd/(2c_S), η = [(ν/w) sin w]² with
  w = √(ν² + ξ²). Since ν/w ≤ 1 in rounding too, η ≤ 1.
- Reflection (c_S < 0): ν = κd/√(c_R|c_S|), ξ = ϑd/(2|c_S|), η = ν²/(ν² + F), with
  F = (s/sinh s)² for s² = ν² − ξ² > 0, and (q/sin q)² for q² = ξ² − ν².

Bragg: |cos(φ − θ_B)| = λ/(2nΛ). The entering root nearest the normal is returned;
a tie goes to the side of sin φ, so θ_B(−φ) = −θ_B(φ).

Refused:
- grazing geometry (c_R ≤ 2⁻⁵⁰, c_S = 0, σ = 0);
- a sine argument the doubles do not resolve to 2⁻¹¹ rad against max(w, 1), where the
  uncertainty of w² is ν²ε_ν + (|ξ| + δξ)δξ;
- inside a reflection band, an η whose propagated relative error passes 2⁻¹¹, or a band
  edge within the uncertainty of s². Validity: Klein Q = 2πλd/(nΛ²) ≫ 1, n₁ ≪ n, lossless,
two waves.

## Gaussian Beams (`wave_beam`, 2.6.0)

⚠ **Conventions**, stated because sources disagree on every one of them:
time dependence e^{+iωt} (Siegman; Saleh & Teich), so the omitted carrier is
e^{−ikz}; q = z + i·z_R with z from the waist, positive downstream; w is the 1/e²
**intensity** radius; `wavelength` is the wavelength in the medium, λ₀/n; ABCD
matrices act on (height, **real** slope), so det = n_in/n_out.

### TEM00
z_R = π w₀² / λ,  w(z) = w₀ √(1 + (z/z_R)²),  1/R(z) = z / (z² + z_R²),
ψ(z) = atan(z/z_R),  θ = λ / (π w₀),  I(r, z) = (2P / π w²) e^{−2r²/w²}.
Power inside radius a: 1 − e^{−2a²/w²} (computed with Kahan's expm1 identity).

### Real beams (ISO 11146)
W(z) = W₀ √(1 + (z M² λ / π W₀²)²),  M² = π W₀ θ / λ,  BPP = W₀θ = M² λ / π.
Second-moment radius W = 2σ, σ² = Σ I (x − x̄)² / Σ I.

### q Parameter and ABCD Propagation
1/q = 1/R − i λ / (π w²),  q' = (A q + B) / (C q + D),  λ' = λ · (AD − BC).
Gouy phase through an element: Δψ = −arg(A + B/q) — 0 for B = 0 with A > 0 (lens, mirror,
interface), π for B = 0 with A < 0 (an inverting relay: 4f telescope, 2f–2f imaging).
Elements: free space [[1, d], [0, 1]]; thin lens [[1, 0], [−1/f, 1]]; mirror
[[1, 0], [−2/R, 1]] (R > 0 concave); curved interface [[1, 0], [(n₁ − n₂)/(n₂R), n₁/n₂]]
(R > 0 with the centre downstream).

### Resonator Eigenmode (Kogelnik & Li 1966)
For a round trip [[A, B], [C, D]] with AD − BC = 1 and |A + D| < 2:
1/q = (D − A)/(2B) − i √(4 − (A + D)²) / (2|B|).
Edges of stability (|A + D| → 2: confocal, concentric, planar) select no unique mode — a
symmetric confocal round trip is −I — and are refused within 1e-12·(|A| + |D|).

### Mode Overlap
η = 4 / ((w₁/w₂ + w₂/w₁)² + (π w₁ w₂ / λ)² (1/R₁ − 1/R₂)²) = 4 Im q₁ Im q₂ / |q₁ − q₂*|².

### Hermite-Gaussian and Laguerre-Gaussian Modes (unit power)
HG_mn: u = (√2/w) φ_m(√2x/w) φ_n(√2y/w) e^{−ik(x²+y²)/2R} e^{i(m+n+1)Ψ},
φ_k(t) = H_k(t) e^{−t²/2} / √(2^k k! √π).
LG_pl: u = √(2 p! / π (p+|l|)!) (1/w) s^{|l|/2} L_p^{|l|}(s) e^{−s/2}
e^{−ikr²/2R} e^{−ilφ} e^{i(2p+|l|+1)Ψ},  s = 2r²/w².
Ψ is the accumulated Gouy phase. LG₀,±₁ = (HG₁₀ ∓ i HG₀₁)/√2;
LG₁₀ = −(HG₂₀ + HG₀₂)/√2. M²: HG_mn → (2m+1, 2n+1); LG_pl → 2p + |l| + 1.
Recurrences: H_{k+1} = 2x H_k − 2k H_{k−1} (A&S 22.7.13);
(k+1) L_{k+1}^α = (2k+1+α−x) L_k^α − (k+α) L_{k−1}^α (A&S 22.7.12).

### Orbital Angular Momentum
e^{−ilφ} in this time convention carries +lħ per photon. Angular-momentum flux
l P / ω = l P λ₀ / (2π c). Spiral phase plate step h = l λ / (n_rel − 1).

## Spectral

### Planck's Law
L(λ,T) = 2hc² / (λ⁵(exp(hc/λkT) - 1))

### Wien's Law
λ_max = b/T where b = 2.898 x 10⁻³ m K

### Photon Energy
E = hc/λ

### CIE XYZ
X = ∫ S(λ) x̄(λ) dλ (similarly Y, Z)

### XYZ to xyY
x = X/(X+Y+Z), y = Y/(X+Y+Z)

### CCT from Chromaticity (McCamy)
CCT = 449n³ + 3525n² + 6823.3n + 5520.33 where n = (x - 0.3320)/(0.1858 - y)

### Luminance (Rec. 709)
Y = 0.2126R + 0.7152G + 0.0722B

## Lens

### Thin Lens
1/f = 1/dₒ + 1/dᵢ

### Magnification
M = -dᵢ/dₒ

### Lensmaker's Equation
1/f = (n-1)(1/R₁ - 1/R₂)

### Thick Lens
1/f = (n-1)[1/R₁ - 1/R₂ + (n-1)d/(nR₁R₂)]

### Optical Power
P = 1/f (diopters when f in meters)

### f-Number
N = f/D

### Numerical Aperture
NA = n sin θ

### Diffraction Limit
θ = 1.22 λ/D

### Airy Disk Radius
r = 1.22 λ N

### Field of View
FOV = 2 arctan(s / 2f)

### MTF Cutoff
f_c = 1 / (λN)

### Seidel Aberrations
S₁ (spherical), S₂ (coma), S₃ (astigmatism), S₄ (field curvature), S₅ (distortion)

For a thin lens in air with the stop at the lens — shape factor q and conjugate
factor p below, h the marginal ray height, H the Lagrange invariant — these are
Welford's sums (*Aberrations of Optical Systems*, thin-lens chapter). They are the
most-repaired formulas in the library (2.2.6, 2.2.8, 2.4.6, 2.5.0):

S₁ = (h⁴φ³/4) · [ (n+2)/(n(n-1)²) q² + 4(n+1)/(n(n-1)) qp + (3n+2)/n p² + n²/(n-1)² ]

S₂ = −(h²φ²H/2) · [ (n+1)/(n(n-1)) q + (2n+1)/n p ]

S₃ = H²φ,  S₄ = H²φ/n,  S₅ = 0

`lens_seidel_coefficients` returns them per unit aperture and field:
`spherical` = S₁/h⁴, `coma` = −S₂/(h²H) = (φ²/2)·[…], `astigmatism` = S₃/H² = φ,
`field_curvature` = S₄/H² = φ/n, `distortion` = 0. The sign on `coma` is the one
the function has always had; physically, **`coma` < 0 is a flare pointing away
from the axis**, the ordinary coma of a biconvex singlet.

What they predict, each checked against `trace_sequential` (below):

| Quantity | Formula | Where |
|----------|---------|-------|
| Longitudinal SA, z_paraxial − z_marginal | S₁/(2u′²), u′ = hφ → h²𝔅/(8f) | `lens_longitudinal_spherical_aberration` |
| Wavefront aberration at paraxial focus | −S₁/8 → −h⁴𝔅/(32f³) | `optical_path_difference`, `opd_fan` |
| Tangential coma, field angle θ | −(3/2) · `coma` · h²θf | — |

with 𝔅 = (3n+2)/n + n²/(n-1)², the bracket at q = 0, p = −1 (equiconvex, object at
infinity): 13.333 at n = 1.5. **`lens_longitudinal_spherical_aberration` is the
same bracket, not a separate formula** — it disagreed with S₁ until 2.4.6.

**Two kinds of check, and each is blind to what the other sees.**

- *The bracket's shape.* With p = −1 the spherical-minimum "best-form" shape is
  q = 2(n²−1)/(n+2) = **0.7143 at n = 1.5** (0.8667 / 1.0216 / 1.1789 at 1.6–1.8),
  and the coma-free "aplanatic" shape is q = **0.80 at n = 1.5**. A bracket that
  does not reproduce both is wrong, whatever it was copied from. But an argmin and
  a zero crossing depend only on the RATIOS of the terms inside the bracket, so
  both are blind to the multiplier in front of it.
- *The magnitude.* tests/lens.tcyr traces an equiconvex singlet (radii ±R, 1
  thick) and compares code/traced. The residual is lens thickness, which thin-lens
  theory leaves out, so it shrinks as R grows; a wrong multiplier would not:

  | n = 1.5 | R = 100 | R = 200 | R = 400 |
  |---------|---------|---------|---------|
  | LSA | 1.0047 | 1.0024 | 1.0012 |
  | S₁/h⁴ | 1.0047 | 1.0024 | 1.0012 |
  | tangential coma | 1.0159 | 1.0082 | 1.0041 |

⛔ **From 2.4.7 to 2.4.9 this page published φ³·n/(4(n−1)²) and φ²·1/(2(n−1)) as
the multipliers**, and called the two shape checks "the reason to trust these
rather than the algebra". They were wrong: `spherical` came out n/(n−1)² too large
(6× at n = 1.5), `coma` 1/(n−1) too large (2×), and LSA — which also divided by φ
where it needed φ² — f(n−1)²/n too small (16.7× at f = 100, and scaling as 1/f²).
Four audits of the bracket passed over them because every check was a shape check
and every value pin was derived from the code. Fixed in 2.5.0, found by the 2.4.8
audit's critic, confirmed against the tracer before any change was made.

### Wavefront Aberration (OPD)
W(ray) = [OPL to the last surface + n′·|Q − P|] − [OPL of the chief ray to P]

P is the chief ray's intercept with the image plane and Q is where the ray leaves
the last surface. The leg Q → P counts negative when the ray travels away from P
(a virtual image). **Every ray must end at the SAME point.** Through 2.4.9 each
ended at its own image-plane intercept, which subtracts h·dW/dh from W: −3W for
spherical aberration, sign flipped. The straight leg Q → P is not the ray's own
path; the error is second order in the ray's angular deviation from QP.

Independent check (tests/ray_simulate.tcyr): the transverse ray error at the image
plane equals f·dW/dh (Hamilton), measured to 0.01–0.32% for h = 1..4 on the
singlet above. The 2.4.9 form gives −0.333 there.

### Wavefront Aberration Coefficients (2.7.0)
From the Seidel sums at marginal height h and Lagrange invariant H (H = hθ > 0 for an
object at infinity with the field on +y, stop at the lens), in
`optical_path_difference`'s sign (ray OPL − chief OPL):

W₀₄₀ = −S_I/8,  W₁₃₁ = +S_II/2,  W₂₂₂ = −S_III/2,  W₂₂₀ = −(S_III + S_IV)/4,  W₃₁₁ = +S_V/2

Welford's W has the opposite overall sign, and his H = n′u′η′ is NEGATIVE for an
image on +y, so the terms even in H flip against his relations and the terms odd in H
(W₁₃₁, W₃₁₁) keep his sign. The first 2.7.0 cut flipped all five; W₁₃₁ then disagreed
in sign with the traced OPD, and an OTF phase conjugated the other way hid it end to end.

W(x, y) = W₀₄₀r⁴ + W₁₃₁yr² + W₂₂₂y² + W₂₂₀r² + W₃₁₁y at the field edge, (x, y)
normalised to the pupil radius. W₂₂₀ carries the sagittal field curvature
(S_III + S_IV); W₂₂₂ + W₂₂₀ the tangential (3S_III + S_IV)/4. tests/lens_otf.tcyr pins
W₀₄₀ to the traced OPD at the pupil edge (0.2% at R = 400) and W₁₃₁, signed, to the odd
part [W(+h) − W(−h)]/2 of a traced tilted fan (0.982 of the Seidel value at θ = 0.01).

### Aberrated OTF (pupil autocorrelation, 2.7.0)
OTF(v) = (1/π) ∬_overlap exp(ik[W(a − v) − W(a + v)]) dA, v = ν/ν_c, ν_c = 1/(λN)

The standard OTF (Goodman eq. 6-25, ISO 9334): the Fourier transform with kernel
e^(−2πiνx) of the PSF |FT{P}|² that `psf_from_wavefront` forms, P = exp(+ikW). A tilt
W₃₁₁y gives MTF_dl·exp(−2πi·W₃₁₁·2v) along y. (Goodman's eq. 6-31 form
W(a + v) − W(a − v) returns the conjugate; the first 2.7.0 cut used it.) The two unit
pupils are shifted by ∓v along the frequency direction (x: sagittal, y: tangential),
and the overlap is integrated exactly — y = √(1−v²)·sin t, then x across the chord —
with n-point Gauss-Legendre in each. That converges exponentially once n resolves the
phase: ~1e-12 at 8 nodes per wave of max|W(a − v) − W(a + v)| for v ≥ 0.02, and v < 0.02
wants 128. The caller's n is a minimum, raised from a bound on that difference:
Seidel min(Σ|Wᵢ|·span, 2v·Σ|Wᵢ|·max|∇|) with spans 1, 2, 1, 1, 2 and gradients 4, 3, 2,
2, 1; Zernike Σ 2|cⱼ|·√(n+1) (m = 0) or √(2(n+1)). Past 256 nodes (32 waves) it is
refused. With W = 0 it is (2/π)(acos v − v√(1−v²)).

Polychromatic: OTF(ν) = Σ wᵢ·OTFᵢ(νλᵢN) / Σ wᵢ, summing the COMPLEX OTFs, so bands
where some wavelengths' contrast reverses cancel. `lens_mtf_polychromatic`'s real
mean of diffraction-limited MTFs equals it only while every OTF is real and ≥ 0.

### Shape Factor
q = (R₂ + R₁) / (R₂ - R₁)

### Conjugate Factor
p = (dᵢ - dₒ) / (dᵢ + dₒ)

### Petzval Sum
Σ 1/(nᵢfᵢ)

### Petzval Radius
R_p = -1/Σ

### Separated Thin Lenses
1/f = 1/f₁ + 1/f₂ - d/(f₁f₂)

### Depth of Field
DoF_near = Hd/(H + d - f), DoF_far = Hd/(H - d + f)

### Paraxial y-nu Trace
Refraction: nu' = nu - y phi
Transfer: y' = y + (t/n) nu

### Chromatic Aberration
Δf = f/V

## PBR

### Fresnel-Schlick
F(θ) = F₀ + (1 - F₀)(1 - cos θ)⁵

### IOR to F0
F₀ = ((n - 1) / (n + 1))²

### GGX/Trowbridge-Reitz NDF
D(h) = α² / (π((n h)²(α² - 1) + 1)²)

### Beckmann NDF
D(h) = exp(-tan²θ / α²) / (πα² cos⁴θ)

### Schlick-GGX Geometry
G₁(v) = (n v) / ((n v)(1 - k) + k)

### Smith Geometry
G(v, l) = G₁(v) G₁(l)

### Cook-Torrance
f = DFG / (4(n v)(n l))

### Lambert Diffuse
f_d = albedo / π

### Anisotropic GGX
D(h) = 1 / (π αₓ αᵧ ((h x/αₓ)² + (h y/αᵧ)² + (n h)²)²)

### Charlie Sheen NDF
D(h) = (2 + 1/α) sin^(1/α) θ / (2π)

### Clearcoat (GGX with fixed F0 = 0.04)
f_coat = D_coat F_coat G_coat,  where G_coat = 1/(4(n v)(n l))

(the 1/(4(n·v)(n·l)) factor lives inside G_coat — writing it again outside would
double-count it, which the implementation does not do)

### Subsurface Scattering (Burley Profile)
R(r) = (e^(-r/d) + e^(-r/3d)) / (8πdr)

### SSS Gaussian Profile
R(r) = (1 / (2πσ²)) exp(-r² / (2σ²))

### SSS Transmittance
T = exp(-σ_t d)

### Iridescence (Thin-Film Interference)
OPD = 2 n_film d cos θ_film, phase shift per channel

### Henyey-Greenstein Phase Function
p(cos θ) = (1 - g²) / (4π(1 + g² - 2g cos θ)^(3/2))

### Rayleigh Phase Function
p(cos θ) = 3(1 + cos²θ) / (16π)

### Volume Transmittance
T = exp(-σ_t d)

### Single Scatter Albedo
ω = σ_s / σ_t

### Importance Sampling GGX
θ = arctan(α √(ξ₁ / (1 - ξ₁)))
φ = 2π ξ₂

### Cosine-Weighted Hemisphere Sampling
pdf = (n l) / π

### Split-Sum Approximation
∫ f(v,l) Lᵢ dω = (∫ Lᵢ dω)(∫ f dω) -- scale + bias from BRDF LUT

## Atmosphere

### Rayleigh Cross-Section
σ(λ) = (8π³/3)(n² - 1)² / (N²λ⁴)

⚠ **n and N must describe the SAME air.** `(n - 1)` is proportional to number
density, so the two do not cancel — pairing an `n` quoted at one temperature with
an `N` at another compounds the error. prakash uses standard air throughout,
**15 °C and 101.325 kPa**: N = 2.5469e25 m⁻³ and n = 1.000277824, the index at
550 nm. ⛔ The textbook `n_air = 1.000293` is air at **0 °C**; prakash paired it
with a ~20 °C N until 2.4.5 and was 12% high in β and τ as a result.

### Rayleigh Scattering Coefficient
β(λ) = N σ(λ)

### Rayleigh at Altitude
β(λ,h) = β₀(λ) exp(-h / H_R)

### Rayleigh Phase Function
p(cos θ) = 3(1 + cos²θ) / (16π)

### Mie Scattering at Altitude
β_M(h) = β_M₀ exp(-h / H_M)

### Mie Phase (Cornette-Shanks)
p(cos θ) = (3(1 - g²)(1 + cos²θ)) / (8π(2 + g²)(1 + g² - 2g cos θ)^(3/2))

### Air Mass (Kasten & Young)
m = 1 / (cos θ + 0.50572(96.07995 - θ)^(-1.6364))

### Optical Depth (Rayleigh)
τ_R(λ) = β_R(λ) H_R

### Atmospheric Transmittance
T(λ,θ) = exp(-(τ_R + τ_M) m(θ))

### Scattering Angle
cos Θ = cos θ_s cos θ_v + sin θ_s sin θ_v cos(Δφ)
