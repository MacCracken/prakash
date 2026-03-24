# Prakash

> **Prakash** (Sanskrit: light, illumination) -- optics and light simulation for AGNOS

Physics of light: ray optics, wave optics, spectral math, lens geometry, atmospheric scattering, and physically-based rendering primitives. Built on [hisab](https://crates.io/crates/hisab) for math foundations.

## Quick Start

```rust
use prakash::ray::{Medium, snell, fresnel_normal, brewster_angle, beer_lambert};
use prakash::spectral::{wavelength_to_rgb, color_temperature_to_rgb, wien_peak};
use prakash::lens::{thin_lens_image_distance, magnification};
use prakash::wave::{malus_law, single_slit_intensity};

// Snell's law: light entering glass from air at 30 degrees
let refracted = snell(Medium::AIR.n, Medium::GLASS.n, 0.5236).unwrap();

// How much light reflects off glass at normal incidence?
let reflectance = fresnel_normal(1.0, 1.52); // ~4%

// What color is 550nm light?
let green = wavelength_to_rgb(550.0).unwrap();

// What color is a 2700K warm bulb?
let warm = color_temperature_to_rgb(2700.0);

// Where does a 50mm lens focus an object at 2 meters?
let image_dist = thin_lens_image_distance(50.0, 2000.0).unwrap();
let mag = magnification(2000.0, image_dist);

// Malus's law: light through a polarizer at 45 degrees
let transmitted = malus_law(1.0, std::f64::consts::FRAC_PI_4); // 50%
```

## Modules

| Module | Description |
|--------|-------------|
| `ray` | Geometric optics: Snell's law, Fresnel equations (s/p/normal/unpolarized), reflection (2D/3D), critical angle, Brewster angle, Beer-Lambert, total internal reflection. 12 built-in materials with refractive indices. Sequential ray tracer, optical system builder, paraxial analysis, cardinal point finder. |
| `wave` | Wave optics: interference intensity, single/double slit diffraction, diffraction grating, thin film reflectance, Malus's law, polarization (Jones vectors/calculus). 2D Fraunhofer diffraction (FFT), multi-source interference patterns, PSF from wavefront error. |
| `spectral` | Wavelength to RGB conversion, Planck blackbody radiance, Wien displacement law, color temperature to RGB, photon energy (J and eV), wavelength/frequency conversion. CIE 1931 2-degree observer, XYZ tristimulus, SPD integration, standard illuminants (A/D50/D65/E/F2/F11), color rendering index (CRI). Physical constants (c, h, k_B). |
| `lens` | Thin lens equation, magnification, lensmaker's equation, combined focal length, depth of field, f-number/NA, mirror focal length, lens classification. |
| `pbr` | Physically-based rendering: Fresnel-Schlick, GGX/Beckmann NDF, geometry Smith, Cook-Torrance specular BRDF, Lambert diffuse, IOR-to-F0 conversion. Full metallic-roughness pipeline. |
| `atmosphere` | Atmospheric optics: Rayleigh and Mie scattering, sky color computation, air mass, optical depth. |
| `error` | Unified error types for all modules. |

## Feature Flags

| Feature | Default | Description |
|---------|---------|-------------|
| `ray` | yes | Geometric optics, sequential ray tracing, optical system builder |
| `wave` | yes | Wave optics, interference, diffraction, polarization, FFT patterns |
| `spectral` | yes | Spectral math, CIE color science, blackbody, color temperature |
| `lens` | yes | Thin lens, lensmaker, mirrors, depth of field |
| `pbr` | no | Cook-Torrance BRDF, Fresnel-Schlick, GGX/Beckmann NDF |
| `atmosphere` | no | Rayleigh/Mie scattering, sky color, air mass |
| `ai` | no | Daimon/hoosh integration (network deps) |
| `logging` | no | Structured logging via `PRAKASH_LOG` env var |
| `full` | -- | Enables all features |

## Examples

| Example | Features | Description |
|---------|----------|-------------|
| [`basic_optics`](examples/basic_optics.rs) | `ray`, `spectral`, `lens` | Snell's law, Fresnel, wavelength-to-RGB, thin lens |
| [`rainbow`](examples/rainbow.rs) | `ray`, `spectral` | Rainbow simulation via raindrop refraction |
| [`camera_lens`](examples/camera_lens.rs) | `ray`, `spectral`, `lens` | Camera lens system analysis with prescriptions |
| [`pbr_materials`](examples/pbr_materials.rs) | `pbr` | PBR material shading with Cook-Torrance BRDF |

Run an example:

```sh
cargo run --example basic_optics --features ray,spectral,lens
```

## Architecture

```
prakash (this crate)
  +-- hisab (math: vectors, geometry, calculus, numerical methods)
```

## Consumer Crates

| Crate | Usage |
|-------|-------|
| soorat | PBR shading |
| kiran | Physically-based lighting |
| ranga | Lens effects |

## Building

```sh
cargo build --all-features
cargo test --all-features
```

## License

GPL-3.0 -- see [LICENSE](LICENSE).
