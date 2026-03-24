# Changelog

## [0.24.3] - 2026-03-23

### Added — wave V0.3: Wave Optics Expansion
- Coherence: temporal coherence length/time, spatial coherence angle/area, coherence ratio
- Circular aperture diffraction: Bessel J₁ (Abramowitz & Stegun rational approximation), Airy pattern, first zero, Rayleigh criterion
- Fabry-Pérot interferometer: transmittance (Airy function), finesse, FSR (Hz and wavelength), resolving power
- Stokes vectors with 7 polarization state constructors and degree of polarization
- Mueller matrices: identity, horizontal/vertical polarizer, arbitrary-angle polarizer, quarter/half-wave plates, general retarder, rotation, matrix multiply, chain application
- Birefringent materials: calcite, quartz, rutile, mica presets; retardation, quarter/half-wave thickness, Mueller generation
- Fraunhofer diffraction: rectangular aperture, 1D arbitrary aperture (numerical DFT)
- Fresnel diffraction: Fresnel number, Fresnel C(x)/S(x) integrals, straight-edge intensity, Fresnel parameter, Huygens-Fresnel 1D numerical integral
- Anti-reflection coatings: ideal AR index, quarter-wave thickness, single-layer reflectance, V-coat, multi-layer transfer matrix method

### Added — pbr V0.4: Advanced PBR
- Anisotropic GGX NDF and Smith geometry (directional roughness)
- Sheen: Charlie distribution (Estevez & Kulla 2017) and Ashikhmin velvet model
- Clearcoat: GGX distribution, fixed-IOR Fresnel, Kelemen geometry, energy-conserving blend
- Subsurface scattering: Burley normalized diffusion profile, Gaussian profile, SSS diffuse term, thin-slab transmittance
- Iridescence: thin-film Fresnel (Airy formula) at single wavelength and RGB, angle/thickness/wavelength-dependent color
- Volumetric scattering: Henyey-Greenstein phase function, Rayleigh phase, isotropic phase, extinction, transmittance, single-scatter albedo, in-scattering
- GGX importance sampling: half-vector sampling, PDF, cosine-weighted hemisphere sampling
- Environment map: split-sum scale/bias (Lazarov analytical fit), mip LOD from roughness, numerical BRDF LUT integration (Hammersley sequence)

### Changed — Infrastructure
- Module refactoring: wave, ray, pbr, spectral split into submodules (largest file: 2174→956 lines)
- CI/CD: GitHub Actions workflows (ci.yml, release.yml) — fmt, clippy, test (Linux/macOS), MSRV 1.89, coverage, docs, security audit, deny, semver check
- Added: SECURITY.md, CONTRIBUTING.md, CODE_OF_CONDUCT.md, deny.toml, codecov.yml
- Added: scripts/bench-history.sh (CSV tracking + benchmarks.md generation), scripts/version-bump.sh
- Makefile: added `all`, `semver` targets
- Test suite: 330 → 487 tests
- Benchmarks: 95 → 137 functions
- Performance: SPD→XYZ 80% faster (aligned fast path), Fresnel edge 15%, CRI 41%, numerous #[inline] additions
- Bug fixes: Rgb::to_u8 rounding, depth_of_field hyperfocal, shape_factor div-by-zero, conjugate_factor div-by-zero

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
