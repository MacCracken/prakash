# Changelog

## [0.23.3] - 2026-03-23

### Added — ray V0.2: Optical Systems
- 3D vector refraction (`refract_3d`) and 3D Snell's law with Fresnel reflectance (`snell_3d`)
- Sequential ray tracing through multiple optical surfaces (`trace_surface`, `trace_sequential`)
- Optical surface types: `SurfaceShape` (Sphere/Plane), `OpticalSurface`, `TraceRay`, `TraceHit`
- Cauchy dispersion model (`CauchyCoefficients`) with BK7 and fused silica presets
- Sellmeier dispersion model (`SellmeierCoefficients`) with 6 presets (BK7, SF11, fused silica, sapphire, water, diamond)
- Abbe number calculation from Sellmeier coefficients
- Fraunhofer spectral line constants (D, F, C)
- Prism deviation, dispersion, and angular spread functions

### Added — lens V0.2: Advanced Optics
- Thick lens equation and cardinal points (FFD, BFD, principal planes)
- F-number, aperture diameter, numerical aperture, NA from f-number
- Diffraction limit (Rayleigh criterion) and Airy disk radius
- Field of view (horizontal and diagonal)
- Diffraction-limited MTF (cutoff frequency and modulation curve)
- Seidel aberration coefficients (spherical, coma, astigmatism, field curvature, distortion)
- Shape factor and conjugate factor for lens analysis
- Longitudinal spherical aberration, chromatic aberration
- Petzval sum and Petzval radius for field curvature analysis
- Separated two-lens system (focal length, BFD)
- System magnification for multi-element systems

### Added — spectral V0.2: Color Science
- CIE 1931 2° standard observer color matching functions (81 entries, 380–780nm @ 5nm)
- `Xyz` tristimulus type with XYZ↔xyY and XYZ↔sRGB conversions
- sRGB gamma correction (linear↔gamma) with proper IEC 61966-2-1 transfer function
- Linear sRGB↔XYZ matrix conversions (D65 white point)
- Correlated color temperature from xy chromaticity (McCamy's approximation)
- CIE CMF interpolation at arbitrary wavelengths
- `Spd` spectral power distribution type with interpolation, XYZ integration, and sRGB conversion
- Blackbody SPD generator
- Standard illuminants: D65, D50, A, F2, F11
- Color Rendering Index (CRI Ra) calculation

### Fixed
- ai: replaced `anyhow` dependency with native `PrakashError`
- ai: `DaimonClient::new()` returns `Result` instead of panicking
- spectral: `wavelength_to_rgb` rejects NaN input
- spectral: `Rgb::to_u8` now rounds instead of truncating (0.999→255)
- lens: `depth_of_field` returns `f64::INFINITY` at hyperfocal distance instead of negative
- lens: `shape_factor` guards against division by zero when r1 == r2
- lens: `optical_power` and `combined_focal_length` return `Result` for zero inputs
- ray: `critical_angle` uses static error string instead of allocating on hot path
- ray: Sellmeier `n_at` guards against resonance pole division-by-zero and negative n²
- pbr: fixed operator precedence in GGX, Beckmann, and geometry functions

### Changed
- Test suite: 92 → 330 tests
- Benchmarks: 24 → 95 benchmarked functions
- Performance improvements across all modules via `#[inline]`, precomputed constants, eliminated redundant computation
- `snell_3d` delegates to `refract_3d` (single source of truth)
- Sequential ray trace sphere intersection assumes normalized direction (eliminates redundant dot product)

## [0.1.0] - 2026-03-23

### Added
- ray: 12 material refractive indices, Snell's law with TIR, Fresnel equations (s/p/unpolarized/normal), reflection (2D/3D), critical angle, Brewster's angle, Beer-Lambert attenuation
- wave: interference (constructive/destructive detection), single/double-slit diffraction, diffraction grating, thin film interference, polarization (Jones vectors), Malus's law
- spectral: wavelength→RGB (CIE 1931), Planck blackbody radiation, Wien's displacement law, color temperature→RGB, wavelength↔frequency, photon energy (J and eV), physical constants
- lens: thin lens equation, magnification, lensmaker's equation, optical power, mirrors, combined focal length, lens classification, depth of field
- pbr: Fresnel-Schlick (scalar/RGB), GGX and Beckmann NDF, Schlick-GGX geometry, Smith geometry, Cook-Torrance specular BRDF, Lambert diffuse (scalar/RGB), IOR→F0
- error: PrakashError with #[non_exhaustive], 7 variants
