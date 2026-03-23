# Prakash Roadmap

> **Prakash** is the optics/light simulation crate. Math foundations come from [hisab](https://github.com/MacCracken/hisab). Rendering integration is via [aethersafta](https://github.com/MacCracken/aethersafta) and [kiran](https://github.com/MacCracken/kiran).

## Scope

Prakash owns the **physics of light**: how light travels, bends, reflects, interferes, diffracts, polarizes, and interacts with materials. It provides the math; consumers decide what to do with it (render pixels, simulate experiments, process images).

Prakash does NOT own:
- **Rendering pipeline** → aethersafta/kiran (they consume prakash for lighting math)
- **Image processing** → ranga (pixel operations, color spaces, filters)
- **Math primitives** → hisab (vectors, geometry, calculus)
- **Color science beyond spectral** → ranga (ICC profiles, gamut mapping)

## V0.1 — Foundation (done, 2026-03-22)

### ray
- 12 material prefractive indices (vacuum, air, water, glass, diamond, etc.)
- Snell's law with total internal reflection detection
- Critical angle calculation
- Reflection (2D and 3D vector reflection)
- Fresnel equations (s-polarized, p-polarized, unpolarized, normal incidence)
- Brewster's angle
- Beer-Lambert attenuation

### wave
- Two-wave interference intensity
- Constructive/destructive interference detection
- Path difference ↔ phase difference
- Thin film interference
- Single-slit diffraction (sinc² envelope)
- Double-slit diffraction (envelope × interference)
- Diffraction grating maxima (multiple orders)
- Polarization state (Jones vector)
- Malus's law

### spectral
- Wavelength → RGB (CIE 1931 piecewise approximation)
- Planck's blackbody radiation
- Wien's displacement law (peak wavelength)
- Color temperature → RGB (Tanner Helland algorithm)
- Wavelength ↔ frequency conversion
- Photon energy (Joules and eV)
- Physical constants (c, h, k)

### lens
- Thin lens equation (image distance from focal + object distance)
- Magnification
- Lensmaker's equation (focal from radii + refractive index)
- Optical power (diopters)
- Mirror focal length and image distance
- Combined focal length (two lenses in contact)
- Lens classification (converging/diverging)
- Depth of field

### pbr
- Fresnel-Schlick approximation (scalar and RGB)
- GGX/Trowbridge-Reitz normal distribution function
- Beckmann normal distribution function
- Schlick-GGX geometry function
- Smith's geometry function
- Cook-Torrance specular BRDF
- Lambert diffuse BRDF (scalar and RGB)
- IOR → F0 conversion

## V0.2 — Optical Systems

### ray
- [ ] Snell's law for 3D vectors (not just angles)
- [ ] Ray tracing through multiple surfaces (sequential ray trace)
- [ ] Dispersion: wavelength-dependent refractive index (Cauchy/Sellmeier equations)
- [ ] Prism dispersion (angular spread by wavelength)

### lens
- [ ] Thick lens equation
- [ ] Lens aberrations: spherical, chromatic, coma, astigmatism, distortion, field curvature
- [ ] Multi-element lens system (sequential surfaces)
- [ ] Aperture and f-number effects
- [ ] MTF (Modulation Transfer Function) calculation

### spectral
- [ ] CIE XYZ color matching functions (tabulated, interpolated)
- [ ] XYZ ↔ sRGB conversion matrix
- [ ] Spectral power distribution (SPD) type with integration
- [ ] Standard illuminants (D50, D65, A, F-series)
- [ ] Color rendering index (CRI) calculation

## V0.3 — Wave Optics Expansion

### wave
- [ ] Huygens-Fresnel diffraction integral (numerical)
- [ ] Fraunhofer diffraction (far-field patterns)
- [ ] Fresnel diffraction (near-field patterns)
- [ ] Circular aperture diffraction (Airy disk, Rayleigh criterion)
- [ ] Coherence length and time
- [ ] Fabry-Pérot interferometer (cavity resonance)
- [ ] Anti-reflection coating design (quarter-wave, multi-layer)
- [ ] Stokes parameters and Mueller matrices (full polarization calculus)
- [ ] Birefringence (ordinary/extraordinary rays)

## V0.4 — Advanced PBR

### pbr
- [ ] Subsurface scattering (BSSRDF approximation)
- [ ] Anisotropic GGX (for brushed metals, hair)
- [ ] Sheen (for fabric, velvet)
- [ ] Clearcoat (for car paint, lacquered surfaces)
- [ ] Iridescence (thin-film on surface — pearls, soap bubbles, beetle shells)
- [ ] Volumetric scattering (Henyey-Greenstein phase function, participating media)
- [ ] Importance sampling for GGX distribution
- [ ] Environment map sampling (pre-filtered, split-sum)

## V0.5 — Simulation & Visualization

- [ ] Ray tracer: recursive ray tracing through hisab geometry (spheres, planes, triangles)
- [ ] Optical bench simulator: define surfaces, trace rays, compute image
- [ ] Spectrum visualization: SPD → RGB image strip
- [ ] Interference pattern visualization: 2D intensity maps
- [ ] Diffraction pattern computation: far-field intensity from arbitrary aperture shape
- [ ] Atmospheric scattering (Rayleigh for sky color, Mie for clouds/fog)

## V1.0 — Stable API

- [ ] API stabilization and review
- [ ] Comprehensive documentation with physics explanations
- [ ] Criterion benchmarks with CSV history
- [ ] Feature gates optimized (minimize compile time for partial usage)
- [ ] Published to crates.io
- [ ] Example: physically accurate rainbow simulation
- [ ] Example: camera lens simulator
- [ ] Example: PBR material preview

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
