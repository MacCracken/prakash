# Prakash Roadmap

> **Prakash** is the optics/light simulation crate. Math foundations come from [hisab](https://github.com/MacCracken/hisab). Rendering integration is via [aethersafta](https://github.com/MacCracken/aethersafta) and [kiran](https://github.com/MacCracken/kiran).

## Scope

Prakash owns the **physics of light**: how light travels, bends, reflects, interferes, diffracts, polarizes, and interacts with materials. It provides the math; consumers decide what to do with it (render pixels, simulate experiments, process images).

Prakash does NOT own:
- **Rendering pipeline** → aethersafta/kiran (they consume prakash for lighting math)
- **Image processing** → ranga (pixel operations, color spaces, filters)
- **Math primitives** → hisab (vectors, geometry, calculus)
- **Color science beyond spectral** → ranga (ICC profiles, gamut mapping)

## V0.1 — Foundation (done)

### ray
- [x] 12 material refractive indices (vacuum, air, water, glass, diamond, etc.)
- [x] Snell's law with total internal reflection detection
- [x] Critical angle calculation
- [x] Reflection (2D and 3D vector reflection)
- [x] Fresnel equations (s-polarized, p-polarized, unpolarized, normal incidence)
- [x] Brewster's angle
- [x] Beer-Lambert attenuation

### wave
- [x] Two-wave interference intensity
- [x] Constructive/destructive interference detection
- [x] Path difference ↔ phase difference
- [x] Thin film interference
- [x] Single-slit diffraction (sinc² envelope)
- [x] Double-slit diffraction (envelope × interference)
- [x] Diffraction grating maxima (multiple orders)
- [x] Polarization state (Jones vector)
- [x] Malus's law

### spectral
- [x] Wavelength → RGB (CIE 1931 piecewise approximation)
- [x] Planck's blackbody radiation
- [x] Wien's displacement law (peak wavelength)
- [x] Color temperature → RGB (Tanner Helland algorithm)
- [x] Wavelength ↔ frequency conversion
- [x] Photon energy (Joules and eV)
- [x] Physical constants (c, h, k)

### lens
- [x] Thin lens equation (image distance from focal + object distance)
- [x] Magnification
- [x] Lensmaker's equation (focal from radii + refractive index)
- [x] Optical power (diopters)
- [x] Mirror focal length and image distance
- [x] Combined focal length (two lenses in contact)
- [x] Lens classification (converging/diverging)
- [x] Depth of field

### pbr
- [x] Fresnel-Schlick approximation (scalar and RGB)
- [x] GGX/Trowbridge-Reitz normal distribution function
- [x] Beckmann normal distribution function
- [x] Schlick-GGX geometry function
- [x] Smith's geometry function
- [x] Cook-Torrance specular BRDF
- [x] Lambert diffuse BRDF (scalar and RGB)
- [x] IOR → F0 conversion

## V0.2 — Optical Systems (done)

### ray
- [x] Snell's law for 3D vectors (not just angles)
- [x] Ray tracing through multiple surfaces (sequential ray trace)
- [x] Dispersion: wavelength-dependent refractive index (Cauchy/Sellmeier equations)
- [x] Prism dispersion (angular spread by wavelength)
- [x] Abbe number calculation
- [x] Sellmeier presets: BK7, SF11, fused silica, sapphire, water, diamond

### lens
- [x] Thick lens equation
- [x] Cardinal points (FFD, BFD, principal planes)
- [x] Lens aberrations: spherical, chromatic, coma, astigmatism, distortion, field curvature (Seidel coefficients)
- [x] Multi-element lens system (separated lenses, system magnification)
- [x] Aperture, f-number, numerical aperture
- [x] MTF (Modulation Transfer Function) calculation
- [x] Diffraction limit, Airy disk radius
- [x] Field of view (horizontal and diagonal)
- [x] Petzval sum and field curvature

### spectral
- [x] CIE 1931 2° standard observer (81-entry table, interpolated)
- [x] XYZ ↔ sRGB conversion matrix (D65 white point)
- [x] sRGB gamma correction (IEC 61966-2-1)
- [x] XYZ ↔ xyY chromaticity
- [x] Spectral power distribution (SPD) type with integration
- [x] Standard illuminants (D50, D65, A, F2, F11)
- [x] Color rendering index (CRI Ra) calculation
- [x] Correlated color temperature from chromaticity (McCamy)

## V0.3 — Wave Optics Expansion (done)

### wave
- [x] Coherence length and time (temporal + spatial)
- [x] Circular aperture diffraction (Bessel J₁, Airy pattern, Rayleigh criterion)
- [x] Fabry-Pérot interferometer (transmittance, finesse, FSR, resolving power)
- [x] Stokes parameters (7 constructors, degree of polarization, ellipticity, orientation)
- [x] Mueller matrices (polarizers, wave plates, retarder, rotation, chain application)
- [x] Birefringence (calcite, quartz, rutile, mica; retardation, wave plate thickness)
- [x] Fraunhofer diffraction (rectangular aperture, 1D arbitrary aperture)
- [x] Fresnel diffraction (Fresnel integrals C/S, straight-edge, Fresnel number)
- [x] Huygens-Fresnel diffraction integral (1D numerical)
- [x] Anti-reflection coatings (quarter-wave, V-coat, multi-layer transfer matrix)

## V0.4 — Advanced PBR (done)

### pbr
- [x] Anisotropic GGX NDF and Smith geometry
- [x] Sheen: Charlie distribution + Ashikhmin velvet
- [x] Clearcoat: distribution, Fresnel, geometry, energy-conserving blend
- [x] Subsurface scattering: Burley + Gaussian profiles, SSS diffuse, transmittance
- [x] Iridescence: thin-film Fresnel (Airy formula), RGB color shift
- [x] Volumetric scattering: Henyey-Greenstein, Rayleigh, isotropic phase; extinction, transmittance, albedo, in-scattering
- [x] Importance sampling: GGX half-vector + PDF, cosine hemisphere + PDF
- [x] Environment map: split-sum scale/bias, mip LOD, BRDF LUT integration (Hammersley)

## V0.25 — Atmospheric Optics (done)

- [x] Rayleigh scattering coefficient (λ⁻⁴ wavelength dependence)
- [x] Rayleigh sky color at arbitrary sun angle
- [x] Mie scattering for aerosols (extinction + phase function)
- [x] Optical depth and air mass calculation
- [x] Sunset/sunrise color gradient model

## V0.26 — Simulation Primitives

- [ ] Recursive ray tracer (reflect/refract through hisab geometry)
- [ ] Ray fan generator (marginal, chief, meridional rays)
- [ ] Spot diagram computation (ray bundle through optical system)
- [ ] Optical path difference (OPD) calculation

## V0.27 — Optical Bench

- [ ] Optical system builder (chain of surfaces with spacings)
- [ ] System cardinal point finder (trace marginal + chief rays)
- [ ] Paraxial ray trace (y-nu method)
- [ ] System prescription format (serialize/deserialize)

## V0.28 — Pattern Computation

- [ ] 2D diffraction pattern from arbitrary aperture (FFT-based)
- [ ] Interference pattern generator (N-source, 2D grid)
- [ ] SPD → RGB strip visualization data
- [ ] PSF (point spread function) from wavefront

## V0.29 — Examples & Documentation

- [ ] Example: physically accurate rainbow simulation
- [ ] Example: camera lens simulator (multi-element trace + spot diagram)
- [ ] Example: PBR material preview (Cook-Torrance + clearcoat + iridescence)
- [ ] Physics explanations in module-level documentation

## V1.0 — Stable Release

- [ ] API review: naming consistency, parameter ordering
- [ ] `#[must_use]` on pure functions and accessors
- [ ] Feature gate audit (minimize compile time for partial usage)
- [ ] Documentation coverage check (all public items documented)
- [ ] README update with full feature matrix

## Consumers

| Consumer | What it uses |
|----------|-------------|
| **kiran** | PBR shading (Cook-Torrance, Fresnel-Schlick), lighting math |
| **aethersafta** | Color temperature for display calibration, compositing light math |
| **rasa** | Lens blur (DoF simulation), chromatic aberration filter, light leak effects |
| **joshua** | Simulation: laser/mirror puzzles, optics experiments, atmospheric effects |
| **tazama** | Color grading: color temperature adjustment, spectral analysis of footage |

## Boundary with Other Crates

| Feature | prakash | other |
|---------|---------|-------|
| Fresnel reflectance math | Yes | — |
| Pixel-level image filter | — | ranga |
| wgpu shader code | — | aethersafta |
| 3D scene graph | — | kiran |
| Vector/matrix math | — | hisab |
| Ray-geometry intersection | — | hisab (geo module) |
| Color space conversion (ICC) | — | ranga |
| Spectral → RGB conversion | Yes | — |
