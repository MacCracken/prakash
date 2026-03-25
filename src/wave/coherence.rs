//! Temporal and spatial coherence functions.

// ── Coherence ─────────────────────────────────────────────────────────────

/// Temporal coherence length.
///
/// l_c = λ² / Δλ
///
/// `center_wavelength` and `bandwidth` in same units (e.g., nm or m).
/// Returns coherence length in the same units.
#[must_use]
#[inline]
pub fn coherence_length(center_wavelength: f64, bandwidth: f64) -> f64 {
    center_wavelength * center_wavelength / bandwidth
}

/// Temporal coherence time.
///
/// τ_c = 1 / Δν = λ² / (c · Δλ)
///
/// `center_wavelength_m` and `bandwidth_m` in meters.
/// Returns coherence time in seconds.
#[must_use]
#[inline]
pub fn coherence_time(center_wavelength_m: f64, bandwidth_m: f64) -> f64 {
    coherence_length(center_wavelength_m, bandwidth_m) / 299_792_458.0
}

/// Spatial coherence angle (van Cittert-Zernike theorem).
///
/// θ_c = λ / d
///
/// `wavelength` and `source_diameter` in same units.
/// Returns the coherence angle in radians.
#[must_use]
#[inline]
pub fn spatial_coherence_angle(wavelength: f64, source_diameter: f64) -> f64 {
    wavelength / source_diameter
}

/// Spatial coherence area at distance R from a source.
///
/// A_c = (λ · R / d)²
///
/// `wavelength`, `distance`, and `source_diameter` in same units.
/// Returns the coherence area in the same units squared.
#[must_use]
#[inline]
pub fn coherence_area(wavelength: f64, distance: f64, source_diameter: f64) -> f64 {
    let l = wavelength * distance / source_diameter;
    l * l
}

/// Number of coherence lengths that fit in a given path difference.
///
/// A visibility metric: values >> 1 mean the source is incoherent at this path difference.
#[must_use]
#[inline]
pub fn coherence_ratio(path_difference: f64, center_wavelength: f64, bandwidth: f64) -> f64 {
    path_difference / coherence_length(center_wavelength, bandwidth)
}
