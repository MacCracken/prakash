//! Prakash — Optics and light simulation for AGNOS
//!
//! Sanskrit: प्रकाश (prakash) — light, illumination, to make visible
//!
//! Provides ray optics, wave optics, spectral math, lens geometry, and
//! physically-based rendering primitives. Built on [hisab](https://crates.io/crates/hisab)
//! for math foundations.
//!
//! # Modules
//!
//! - [`ray`] — Geometric optics: reflection, refraction, Snell's law, Fresnel, total internal reflection
//! - [`wave`] — Wave optics: interference, diffraction, polarization
//! - [`spectral`] — Wavelength ↔ RGB, blackbody radiation, color temperature, spectral power
//! - [`lens`] — Lens/mirror geometry: focal length, magnification, thin lens, aberrations
//! - [`pbr`] — Physically-based rendering: BRDF, Cook-Torrance, Fresnel-Schlick
//! - [`error`] — Error types

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

#[cfg(feature = "logging")]
pub mod logging;

#[cfg(feature = "ai")]
pub mod ai;

pub use error::PrakashError;
