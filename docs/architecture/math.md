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
