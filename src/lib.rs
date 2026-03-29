#![warn(missing_docs)]

//! Prakash — Optics and light simulation for AGNOS
//!
//! Sanskrit: प्रकाश (prakash) — light, illumination, to make visible
//!
//! Provides ray optics, wave optics, spectral math, lens geometry,
//! physically-based rendering, and atmospheric optics. Built on
//! [hisab](https://crates.io/crates/hisab) for math foundations.
//!
//! # Modules
//!
//! - [`ray`] — Geometric optics: reflection, refraction, Snell's law, Fresnel, dispersion, sequential/recursive ray tracing, ray fans, spot diagrams, OPD
//! - [`wave`] — Wave optics: interference, diffraction, polarization (Jones/Stokes/Mueller), coherence, Fabry-Pérot, anti-reflection coatings, 2D patterns, PSF
//! - [`spectral`] — Wavelength ↔ RGB, blackbody radiation, CIE 1931, SPD, illuminants, CRI, color temperature
//! - [`lens`] — Lens/mirror geometry: thin/thick lens, aberrations (Seidel), MTF, DoF, Petzval, multi-element systems
//! - [`pbr`] — PBR shading: Cook-Torrance, GGX, sheen, clearcoat, SSS, iridescence, volumetric, importance sampling, environment maps
//! - [`atmosphere`] — Atmospheric optics: Rayleigh/Mie scattering, sky color, air mass, optical depth, sunset model
//! - [`error`] — Error types
//!
//! # Wavelength parameter convention
//!
//! - **`_nm`** suffix — wavelength in nanometers (spectral, CIE, visible-range functions)
//! - **`_m`** suffix — wavelength in meters (Planck, atmospheric scattering, coherence)
//! - **`_um`** suffix — wavelength in micrometers (dispersion: Sellmeier, Cauchy)
//! - **unsuffixed** — unit-agnostic; caller provides consistent units (wave optics, lens)
//!
//! Within wave-optics functions, **wavelength is always the first parameter**,
//! followed by geometry (aperture, slit width), then observation angle/position.

pub mod error;

#[cfg(feature = "ray")]
pub mod ray;

#[cfg(feature = "wave")]
pub mod wave;

#[cfg(feature = "spectral")]
pub mod spectral;

#[cfg(feature = "lens")]
pub mod lens;

#[cfg(feature = "pbr")]
pub mod pbr;

#[cfg(feature = "atmosphere")]
pub mod atmosphere;

#[cfg(feature = "logging")]
pub mod logging;

#[cfg(feature = "ai")]
pub mod ai;

pub mod bridge;

pub use error::PrakashError;
