# Prakash

> **Prakash** (Sanskrit: प्रकाश — light, illumination, to make visible) — optics and light simulation for AGNOS

Physics of light: ray optics, wave optics, spectral math, lens geometry, and physically-based rendering primitives. Built on [hisab](https://crates.io/crates/hisab) for math foundations.

## Features

| Feature | Default | Description |
|---------|---------|-------------|
| `ray` | yes | Geometric optics: Snell's law, Fresnel equations, reflection, refraction, TIR, Beer-Lambert |
| `wave` | yes | Wave optics: interference, diffraction (single/double slit, grating), polarization, thin film |
| `spectral` | yes | Wavelength ↔ RGB, blackbody radiation (Planck), Wien's law, color temperature, photon energy |
| `lens` | yes | Thin lens equation, magnification, lensmaker's equation, mirrors, depth of field |
| `pbr` | no | Physically-based rendering: Fresnel-Schlick, GGX/Beckmann NDF, Cook-Torrance BRDF, Lambert |
| `ai` | no | Daimon/hoosh integration (network deps) |
| `logging` | no | Structured logging via `PRAKASH_LOG` env var |
| `full` | — | Enables all features |

## Architecture

```
prakash (this crate)
  └── hisab (math: vectors, geometry, calculus, numerical methods)

Consumers:
  kiran       ──→ prakash (physically-based lighting, lens effects)
  aethersafta ──→ prakash (compositor light math, color temperature)
  rasa        ──→ prakash (image filter optics: lens blur, chromatic aberration)
  joshua      ──→ prakash (simulation: laser puzzles, optics experiments)
```

## Quick Start

```rust
use prakash::ray::{Medium, snell, fresnel_normal, brewster_angle, beer_lambert};
use prakash::spectral::{wavelength_to_rgb, color_temperature_to_rgb, wien_peak};
use prakash::lens::{thin_lens_image_distance, magnification};
use prakash::wave::{malus_law, single_slit_intensity};

// Snell's law: light entering glass from air at 30°
let refracted = snell(Medium::AIR.n, Medium::GLASS.n, 0.5236).unwrap(); // 30° in radians

// How much light reflects off glass at normal incidence?
let reflectance = fresnel_normal(1.0, 1.52); // ~4%

// What color is 550nm light?
let green = wavelength_to_rgb(550.0).unwrap();

// What color is a 2700K warm bulb?
let warm = color_temperature_to_rgb(2700.0);

// Where does a 50mm lens focus an object at 2 meters?
let image_dist = thin_lens_image_distance(50.0, 2000.0).unwrap();
let mag = magnification(2000.0, image_dist);

// Malus's law: light through a polarizer at 45°
let transmitted = malus_law(1.0, std::f64::consts::FRAC_PI_4); // 50%
```

## Modules

| Module | Functions | Description |
|--------|-----------|-------------|
| `ray` | `snell`, `fresnel_*`, `reflect_2d/3d`, `critical_angle`, `brewster_angle`, `beer_lambert` | 12 common materials with refractive indices |
| `wave` | `interference_intensity`, `single_slit_intensity`, `double_slit_intensity`, `grating_maxima`, `thin_film_reflectance`, `malus_law` | Polarization state (Jones vectors) |
| `spectral` | `wavelength_to_rgb`, `planck_radiance`, `wien_peak`, `color_temperature_to_rgb`, `photon_energy` | Physical constants (c, h, k) |
| `lens` | `thin_lens_image_distance`, `magnification`, `lensmaker_focal_length`, `depth_of_field`, `combined_focal_length` | Mirror equations, lens classification |
| `pbr` | `fresnel_schlick`, `distribution_ggx`, `geometry_smith`, `cook_torrance`, `lambert_diffuse`, `ior_to_f0` | Full metallic-roughness PBR pipeline |

## Building

```sh
cargo build --all-features
cargo test --all-features
```

## Roadmap

See [docs/development/roadmap.md](docs/development/roadmap.md).

## License

GPL-3.0 — see [LICENSE](LICENSE).
