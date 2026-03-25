# Mathematical Reference

Quick reference for all physics formulas implemented in prakash.

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
I(θ) = I₀ (sin β / β)² cos²(γ) where γ = πd sin θ / λ

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
CCT = 449n³ + 3525n² + 6823.3n + 5520.33 where n = (x - 0.3320)/(y - 0.1858)

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

### Shape Factor
q = (R₂ + R₁) / (R₂ - R₁)

### Conjugate Factor
p = (dᵢ - dₒ) / (dᵢ + dₒ)

### Petzval Sum
Σ 1/(nᵢfᵢ)

### Petzval Radius
R_p = -1/(2 Σ)

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
f_coat = D_coat F_coat G_coat / (4(n v)(n l))

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
