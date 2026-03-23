# Changelog

## [0.22.3] - 2026-03-23

### Fixed
- ai: replaced `anyhow` dependency with native `PrakashError` — `DaimonClient::new()` now returns `Result` instead of panicking
- ai: `register_agent()` now returns `Result` with proper error context instead of `anyhow::Result`

### Changed
- Expanded test suite: 92 → 212 tests with edge cases, boundary conditions, property checks, cross-module validation, and serde roundtrips
- Expanded benchmarks: 24 → 47 benchmarked functions covering all public APIs
- spectral: precomputed Planck radiation constants (PLANCK_C1, PLANCK_C2) — eliminates per-call multiply chain
- ray: `deg_to_rad`/`rad_to_deg` now use stdlib `to_radians()`/`to_degrees()` for precision
- wave: `grating_maxima` preallocates Vec capacity
- pbr: fixed operator precedence in `distribution_ggx`, `distribution_beckmann`, and `geometry_schlick_ggx` — `.max()` was binding to wrong subexpression, replaced with `+ epsilon` guards
- lens: `optical_power` and `combined_focal_length` now return `Result` to prevent division-by-zero panics
- spectral: `wavelength_to_rgb` now correctly rejects NaN input (was passing through range check due to NaN comparison semantics)
- Integration tests now include cross-module consistency checks (ray↔pbr F0 agreement, Snell transitivity, spectral energy ordering)

## [0.1.0] - 2026-03-23

### Added
- ray: 12 material refractive indices, Snell's law with TIR, Fresnel equations (s/p/unpolarized/normal), reflection (2D/3D), critical angle, Brewster's angle, Beer-Lambert attenuation
- wave: interference (constructive/destructive detection), single/double-slit diffraction, diffraction grating, thin film interference, polarization (Jones vectors), Malus's law
- spectral: wavelength→RGB (CIE 1931), Planck blackbody radiation, Wien's displacement law, color temperature→RGB, wavelength↔frequency, photon energy (J and eV), physical constants
- lens: thin lens equation, magnification, lensmaker's equation, optical power, mirrors, combined focal length, lens classification, depth of field
- pbr: Fresnel-Schlick (scalar/RGB), GGX and Beckmann NDF, Schlick-GGX geometry, Smith geometry, Cook-Torrance specular BRDF, Lambert diffuse (scalar/RGB), IOR→F0
- error: PrakashError with #[non_exhaustive], 7 variants
