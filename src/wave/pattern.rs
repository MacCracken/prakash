//! Pattern computation — 2D diffraction, interference patterns, SPD visualization, PSF.
//!
//! Tools for generating 2D optical patterns: diffraction from arbitrary apertures
//! (via FFT), multi-source interference, spectral color strips, and point spread
//! functions from wavefront data.

use tracing::trace;

// ── Complex Arithmetic (minimal, inline) ────────────────────────────────────

/// Minimal complex number for FFT computation.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Complex {
    re: f64,
    im: f64,
}

impl Complex {
    #[inline]
    const fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    #[inline]
    const fn zero() -> Self {
        Self { re: 0.0, im: 0.0 }
    }

    #[inline]
    fn norm_sq(self) -> f64 {
        self.re * self.re + self.im * self.im
    }

    #[inline]
    fn mul(self, other: Self) -> Self {
        Self {
            re: self.re * other.re - self.im * other.im,
            im: self.re * other.im + self.im * other.re,
        }
    }

    #[inline]
    fn add(self, other: Self) -> Self {
        Self {
            re: self.re + other.re,
            im: self.im + other.im,
        }
    }

    #[inline]
    fn sub(self, other: Self) -> Self {
        Self {
            re: self.re - other.re,
            im: self.im - other.im,
        }
    }

    /// e^(iθ)
    #[inline]
    fn from_phase(theta: f64) -> Self {
        let (s, c) = theta.sin_cos();
        Self { re: c, im: s }
    }
}

// ── FFT (Cooley-Tukey radix-2) ──────────────────────────────────────────────

/// In-place radix-2 Cooley-Tukey FFT.
///
/// `data` length must be a power of 2. `inverse` = true for IFFT.
fn fft_radix2(data: &mut [Complex], inverse: bool) {
    let n = data.len();
    debug_assert!(n.is_power_of_two(), "FFT length must be power of 2");

    // Bit-reversal permutation
    let mut j = 0;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            data.swap(i, j);
        }
    }

    // Butterfly passes
    let sign = if inverse { 1.0 } else { -1.0 };
    let mut len = 2;
    while len <= n {
        let half = len / 2;
        let angle = sign * std::f64::consts::TAU / len as f64;
        let wn = Complex::from_phase(angle);

        let mut start = 0;
        while start < n {
            let mut w = Complex::new(1.0, 0.0);
            for k in 0..half {
                let u = data[start + k];
                let v = data[start + k + half].mul(w);
                data[start + k] = u.add(v);
                data[start + k + half] = u.sub(v);
                w = w.mul(wn);
            }
            start += len;
        }
        len <<= 1;
    }

    // Normalize for inverse
    if inverse {
        let inv_n = 1.0 / n as f64;
        for c in data.iter_mut() {
            c.re *= inv_n;
            c.im *= inv_n;
        }
    }
}

/// Round up to the next power of 2.
#[inline]
fn next_power_of_2(n: usize) -> usize {
    n.next_power_of_two()
}

// ── 2D Diffraction Pattern ──────────────────────────────────────────────────

/// A 2D grid of intensity values.
#[derive(Debug, Clone, PartialEq)]
pub struct Pattern2D {
    /// Intensity values in row-major order.
    pub data: Vec<f64>,
    /// Number of columns (width).
    pub width: usize,
    /// Number of rows (height).
    pub height: usize,
}

impl Pattern2D {
    /// Create a new pattern filled with zeros.
    #[must_use]
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            data: vec![0.0; width * height],
            width,
            height,
        }
    }

    /// Get the intensity at (col, row).
    #[must_use]
    #[inline]
    pub fn get(&self, col: usize, row: usize) -> f64 {
        self.data[row * self.width + col]
    }

    /// Set the intensity at (col, row).
    #[inline]
    pub fn set(&mut self, col: usize, row: usize, value: f64) {
        self.data[row * self.width + col] = value;
    }

    /// Maximum intensity in the pattern.
    #[must_use]
    pub fn max_intensity(&self) -> f64 {
        self.data.iter().copied().fold(0.0_f64, f64::max)
    }

    /// Normalize all values to [0, 1] relative to the maximum.
    #[must_use]
    pub fn normalized(&self) -> Self {
        let max = self.max_intensity();
        if max < 1e-30 {
            return self.clone();
        }
        let inv = 1.0 / max;
        Self {
            data: self.data.iter().map(|&v| v * inv).collect(),
            width: self.width,
            height: self.height,
        }
    }
}

/// Compute the 2D Fraunhofer diffraction pattern from an arbitrary aperture.
///
/// The aperture is defined as a 2D grid of transmission values (0.0 = opaque,
/// 1.0 = fully transparent). The pattern is computed via 2D FFT.
///
/// `aperture` is a row-major grid with dimensions `width × height`.
/// These dimensions will be padded to the next power of 2.
///
/// Returns the far-field intensity pattern (|E|²).
#[must_use]
pub fn diffraction_pattern_2d(aperture: &[f64], width: usize, height: usize) -> Pattern2D {
    trace!(width, height, "diffraction_pattern_2d");
    let nw = next_power_of_2(width);
    let nh = next_power_of_2(height);

    // Create zero-padded complex array
    let mut grid = vec![Complex::zero(); nw * nh];
    for row in 0..height {
        for col in 0..width {
            grid[row * nw + col] = Complex::new(aperture[row * width + col], 0.0);
        }
    }

    // 2D FFT: rows then columns
    // Row FFTs
    for row in 0..nh {
        let start = row * nw;
        fft_radix2(&mut grid[start..start + nw], false);
    }
    // Column FFTs
    let mut col_buf = vec![Complex::zero(); nh];
    for col in 0..nw {
        for row in 0..nh {
            col_buf[row] = grid[row * nw + col];
        }
        fft_radix2(&mut col_buf, false);
        for row in 0..nh {
            grid[row * nw + col] = col_buf[row];
        }
    }

    // Extract |E|² and fftshift (move DC to center)
    let mut pattern = Pattern2D::new(nw, nh);
    for row in 0..nh {
        for col in 0..nw {
            // fftshift: swap quadrants
            let sr = (row + nh / 2) % nh;
            let sc = (col + nw / 2) % nw;
            pattern.set(col, row, grid[sr * nw + sc].norm_sq());
        }
    }
    pattern
}

/// Compute the 2D diffraction pattern of a circular aperture.
///
/// Creates an aperture mask for a circle of given radius centered in a grid,
/// then computes the far-field pattern via FFT. The result is the Airy pattern.
///
/// `grid_size` = size of the computation grid (will be rounded to power of 2).
/// `radius` = aperture radius in grid pixels.
#[must_use]
pub fn diffraction_pattern_circular(grid_size: usize, radius: f64) -> Pattern2D {
    trace!(grid_size, radius, "diffraction_pattern_circular");
    let n = next_power_of_2(grid_size);
    let cx = n as f64 / 2.0;
    let cy = n as f64 / 2.0;
    let r2 = radius * radius;

    let mut aperture = vec![0.0; n * n];
    for row in 0..n {
        for col in 0..n {
            let dx = col as f64 - cx;
            let dy = row as f64 - cy;
            if dx * dx + dy * dy <= r2 {
                aperture[row * n + col] = 1.0;
            }
        }
    }
    diffraction_pattern_2d(&aperture, n, n)
}

// ── Interference Pattern ────────────────────────────────────────────────────

/// A point source for interference pattern computation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointSource {
    /// Position x-coordinate in the source plane.
    pub x: f64,
    /// Position y-coordinate in the source plane.
    pub y: f64,
    /// Amplitude.
    pub amplitude: f64,
    /// Initial phase (radians).
    pub phase: f64,
}

impl PointSource {
    /// Create a new point source.
    #[must_use]
    #[inline]
    pub const fn new(x: f64, y: f64, amplitude: f64, phase: f64) -> Self {
        Self {
            x,
            y,
            amplitude,
            phase,
        }
    }
}

/// Compute the interference pattern from N point sources on a 2D grid.
///
/// For each grid point, sums the complex amplitudes from all sources
/// (accounting for path differences) and returns |Σ E|².
///
/// `sources` = list of point sources.
/// `wavelength` = wavelength (same units as source positions and grid extents).
/// `distance` = propagation distance from source plane to observation plane.
/// `grid_x` = (x_min, x_max) range of the observation grid.
/// `grid_y` = (y_min, y_max) range.
/// `nx`, `ny` = grid resolution.
#[must_use]
pub fn interference_pattern(
    sources: &[PointSource],
    wavelength: f64,
    distance: f64,
    grid_x: (f64, f64),
    grid_y: (f64, f64),
    nx: usize,
    ny: usize,
) -> Pattern2D {
    trace!(num_sources = sources.len(), nx, ny, "interference_pattern");
    let k = std::f64::consts::TAU / wavelength;
    let dx = if nx > 1 {
        (grid_x.1 - grid_x.0) / (nx as f64 - 1.0)
    } else {
        0.0
    };
    let dy = if ny > 1 {
        (grid_y.1 - grid_y.0) / (ny as f64 - 1.0)
    } else {
        0.0
    };

    let mut pattern = Pattern2D::new(nx, ny);

    for row in 0..ny {
        let obs_y = grid_y.0 + row as f64 * dy;
        for col in 0..nx {
            let obs_x = grid_x.0 + col as f64 * dx;

            let mut e_re = 0.0;
            let mut e_im = 0.0;
            for src in sources {
                let sx = obs_x - src.x;
                let sy = obs_y - src.y;
                let r = (sx * sx + sy * sy + distance * distance).sqrt();
                let phase = k * r + src.phase;
                let (sin_p, cos_p) = phase.sin_cos();
                let amp = src.amplitude / r; // 1/r falloff
                e_re += amp * cos_p;
                e_im += amp * sin_p;
            }
            pattern.set(col, row, e_re * e_re + e_im * e_im);
        }
    }
    pattern
}

// ── SPD → RGB Strip ─────────────────────────────────────────────────────────

/// An RGB color triplet for spectrum visualization.
///
/// Separate from [`crate::spectral::Rgb`] to avoid requiring the `spectral` feature.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StripColor {
    /// Red component (0.0--1.0).
    pub r: f64,
    /// Green component (0.0--1.0).
    pub g: f64,
    /// Blue component (0.0--1.0).
    pub b: f64,
}

/// Generate a visible spectrum strip as a sequence of RGB colors.
///
/// Produces `num_samples` colors spanning the visible range (380–780 nm)
/// using the piecewise CIE 1931 approximation from the spectral module.
///
/// Returns colors ready for display (gamma not applied — linear RGB).
#[must_use]
pub fn spectrum_strip(num_samples: usize) -> Vec<StripColor> {
    trace!(num_samples, "spectrum_strip");
    let mut colors = Vec::with_capacity(num_samples);
    for i in 0..num_samples {
        let t = i as f64 / (num_samples as f64 - 1.0).max(1.0);
        let wl = 380.0 + t * 400.0; // 380–780 nm
        let (r, g, b) = wavelength_to_rgb_approx(wl);
        colors.push(StripColor { r, g, b });
    }
    colors
}

/// Generate a spectrum strip for a custom wavelength range.
///
/// `start_nm` and `end_nm` define the range. Wavelengths outside 380–780nm
/// will produce black.
#[must_use]
pub fn spectrum_strip_range(start_nm: f64, end_nm: f64, num_samples: usize) -> Vec<StripColor> {
    trace!(start_nm, end_nm, num_samples, "spectrum_strip_range");
    let mut colors = Vec::with_capacity(num_samples);
    for i in 0..num_samples {
        let t = i as f64 / (num_samples as f64 - 1.0).max(1.0);
        let wl = start_nm + t * (end_nm - start_nm);
        let (r, g, b) = wavelength_to_rgb_approx(wl);
        colors.push(StripColor { r, g, b });
    }
    colors
}

/// Piecewise wavelength→RGB (inline, no spectral feature dependency).
///
/// Same algorithm as `spectral::wavelength_to_rgb` but returns (r,g,b) tuple
/// and doesn't require the spectral feature. Returns (0,0,0) outside visible range.
fn wavelength_to_rgb_approx(nm: f64) -> (f64, f64, f64) {
    if !(380.0..=780.0).contains(&nm) {
        return (0.0, 0.0, 0.0);
    }
    let (r, g, b) = if nm < 440.0 {
        (-(nm - 440.0) / (440.0 - 380.0), 0.0, 1.0)
    } else if nm < 490.0 {
        (0.0, (nm - 440.0) / (490.0 - 440.0), 1.0)
    } else if nm < 510.0 {
        (0.0, 1.0, -(nm - 510.0) / (510.0 - 490.0))
    } else if nm < 580.0 {
        ((nm - 510.0) / (580.0 - 510.0), 1.0, 0.0)
    } else if nm < 645.0 {
        (1.0, -(nm - 645.0) / (645.0 - 580.0), 0.0)
    } else {
        (1.0, 0.0, 0.0)
    };
    let factor = if nm < 420.0 {
        0.3 + 0.7 * (nm - 380.0) / 40.0
    } else if nm > 700.0 {
        0.3 + 0.7 * (780.0 - nm) / 80.0
    } else {
        1.0
    };
    (r * factor, g * factor, b * factor)
}

// ── PSF (Point Spread Function) ─────────────────────────────────────────────

/// Compute the Point Spread Function from wavefront error data.
///
/// The PSF is the squared magnitude of the Fourier transform of the pupil function:
///   PSF = |FT{P(x,y) · exp(i·2π·W(x,y)/λ)}|²
///
/// where P is the pupil mask (1 inside aperture, 0 outside) and W is the
/// wavefront error (OPD) in the same units as `wavelength`.
///
/// `pupil` = aperture mask (row-major, `width × height`): 1.0 inside, 0.0 outside.
/// `wavefront_error` = OPD at each pupil point (same grid, same units as wavelength).
/// `wavelength` = wavelength for phase computation (same units as wavefront_error).
///
/// Returns the PSF as a 2D pattern (power of 2 dimensions, centered).
#[must_use]
pub fn psf_from_wavefront(
    pupil: &[f64],
    wavefront_error: &[f64],
    width: usize,
    height: usize,
    wavelength: f64,
) -> Pattern2D {
    trace!(width, height, wavelength, "psf_from_wavefront");
    let nw = next_power_of_2(width);
    let nh = next_power_of_2(height);
    let k = std::f64::consts::TAU / wavelength;

    // Build complex pupil function: P(x,y) · exp(i·k·W(x,y))
    let mut grid = vec![Complex::zero(); nw * nh];
    for row in 0..height {
        for col in 0..width {
            let idx = row * width + col;
            let p = pupil[idx];
            if p > 0.0 {
                let phase = k * wavefront_error[idx];
                let e = Complex::from_phase(phase);
                grid[row * nw + col] = Complex::new(p * e.re, p * e.im);
            }
        }
    }

    // 2D FFT
    for row in 0..nh {
        let start = row * nw;
        fft_radix2(&mut grid[start..start + nw], false);
    }
    let mut col_buf = vec![Complex::zero(); nh];
    for col in 0..nw {
        for row in 0..nh {
            col_buf[row] = grid[row * nw + col];
        }
        fft_radix2(&mut col_buf, false);
        for row in 0..nh {
            grid[row * nw + col] = col_buf[row];
        }
    }

    // |E|² with fftshift
    let mut pattern = Pattern2D::new(nw, nh);
    for row in 0..nh {
        for col in 0..nw {
            let sr = (row + nh / 2) % nh;
            let sc = (col + nw / 2) % nw;
            pattern.set(col, row, grid[sr * nw + sc].norm_sq());
        }
    }
    pattern
}

/// Compute the diffraction-limited PSF (perfect wavefront, zero OPD).
///
/// Equivalent to `psf_from_wavefront` with all-zero wavefront error.
/// The result is the Airy pattern for a circular aperture.
///
/// `grid_size` = computation grid size, `radius` = aperture radius in pixels.
#[must_use]
pub fn psf_diffraction_limited(grid_size: usize, radius: f64) -> Pattern2D {
    trace!(grid_size, radius, "psf_diffraction_limited");
    diffraction_pattern_circular(grid_size, radius)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-6;

    // ── FFT ──────────────────────────────────────────────────────────────

    #[test]
    fn test_fft_single_element() {
        let mut data = [Complex::new(5.0, 0.0)];
        fft_radix2(&mut data, false);
        assert!((data[0].re - 5.0).abs() < EPS);
    }

    #[test]
    fn test_fft_two_elements() {
        let mut data = [Complex::new(1.0, 0.0), Complex::new(1.0, 0.0)];
        fft_radix2(&mut data, false);
        assert!((data[0].re - 2.0).abs() < EPS); // DC = sum
        assert!(data[1].norm_sq() < EPS); // Nyquist = 0 for constant
    }

    #[test]
    fn test_fft_roundtrip() {
        let original = [
            Complex::new(1.0, 0.0),
            Complex::new(2.0, 0.0),
            Complex::new(3.0, 0.0),
            Complex::new(4.0, 0.0),
        ];
        let mut data = original;
        fft_radix2(&mut data, false);
        fft_radix2(&mut data, true);
        for (a, b) in data.iter().zip(original.iter()) {
            assert!(
                (a.re - b.re).abs() < EPS && (a.im - b.im).abs() < EPS,
                "FFT roundtrip failed"
            );
        }
    }

    #[test]
    fn test_fft_parseval() {
        // Parseval's theorem: sum|x|² = (1/N)·sum|X|²
        let mut data = [
            Complex::new(1.0, 0.0),
            Complex::new(0.0, 1.0),
            Complex::new(-1.0, 0.0),
            Complex::new(0.0, -1.0),
        ];
        let energy_time: f64 = data.iter().map(|c| c.norm_sq()).sum();
        fft_radix2(&mut data, false);
        let energy_freq: f64 = data.iter().map(|c| c.norm_sq()).sum();
        assert!(
            (energy_time * 4.0 - energy_freq).abs() < EPS,
            "Parseval: {energy_time}*4 ≠ {energy_freq}"
        );
    }

    // ── 2D Diffraction ───────────────────────────────────────────────────

    #[test]
    fn test_diffraction_2d_dimensions() {
        let aperture = vec![1.0; 8 * 8];
        let pattern = diffraction_pattern_2d(&aperture, 8, 8);
        assert_eq!(pattern.width, 8);
        assert_eq!(pattern.height, 8);
    }

    #[test]
    fn test_diffraction_2d_non_power_of_2() {
        let aperture = vec![1.0; 6 * 6];
        let pattern = diffraction_pattern_2d(&aperture, 6, 6);
        assert_eq!(pattern.width, 8); // padded to next power of 2
        assert_eq!(pattern.height, 8);
    }

    #[test]
    fn test_diffraction_2d_center_peak() {
        let aperture = vec![1.0; 16 * 16];
        let pattern = diffraction_pattern_2d(&aperture, 16, 16);
        let center = pattern.get(8, 8);
        let corner = pattern.get(0, 0);
        assert!(center > corner, "Center should be brightest");
    }

    #[test]
    fn test_diffraction_2d_all_non_negative() {
        let aperture = vec![1.0; 16 * 16];
        let pattern = diffraction_pattern_2d(&aperture, 16, 16);
        for &v in &pattern.data {
            assert!(v >= 0.0, "Intensity should be non-negative");
        }
    }

    #[test]
    fn test_diffraction_circular_symmetric() {
        let pattern = diffraction_pattern_circular(32, 8.0);
        // Check approximate symmetry: top vs bottom of center
        let above = pattern.get(16, 12);
        let below = pattern.get(16, 20);
        assert!(
            (above - below).abs() / (above + below + 1e-30) < 0.1,
            "Circular pattern should be symmetric"
        );
    }

    #[test]
    fn test_diffraction_circular_center_brightest() {
        let pattern = diffraction_pattern_circular(32, 8.0);
        let max = pattern.max_intensity();
        let center = pattern.get(16, 16);
        assert!(
            (center - max).abs() / max < 0.01,
            "Center should be the maximum"
        );
    }

    // ── Interference pattern ─────────────────────────────────────────────

    #[test]
    fn test_interference_single_source() {
        let sources = [PointSource::new(0.0, 0.0, 1.0, 0.0)];
        let pattern =
            interference_pattern(&sources, 0.5e-6, 1.0, (-1e-3, 1e-3), (-1e-3, 1e-3), 16, 16);
        assert_eq!(pattern.width, 16);
        assert_eq!(pattern.height, 16);
        // Single source → smooth falloff, all positive
        for &v in &pattern.data {
            assert!(v >= 0.0);
        }
    }

    #[test]
    fn test_interference_two_sources_fringes() {
        let sources = [
            PointSource::new(-0.5e-3, 0.0, 1.0, 0.0),
            PointSource::new(0.5e-3, 0.0, 1.0, 0.0),
        ];
        let pattern = interference_pattern(&sources, 0.5e-6, 1.0, (-5e-3, 5e-3), (0.0, 0.0), 64, 1);
        // Should show intensity variation (fringes) along x
        let min = pattern.data.iter().copied().fold(f64::INFINITY, f64::min);
        let max = pattern.data.iter().copied().fold(0.0_f64, f64::max);
        assert!(
            max > min * 1.5,
            "Two sources should produce fringes: min={min}, max={max}"
        );
    }

    #[test]
    fn test_interference_constructive() {
        // Two in-phase sources at same location → 4× intensity of one
        let sources = [
            PointSource::new(0.0, 0.0, 1.0, 0.0),
            PointSource::new(0.0, 0.0, 1.0, 0.0),
        ];
        let one = [PointSource::new(0.0, 0.0, 1.0, 0.0)];
        let p2 = interference_pattern(&sources, 0.5e-6, 1.0, (0.0, 0.0), (0.0, 0.0), 1, 1);
        let p1 = interference_pattern(&one, 0.5e-6, 1.0, (0.0, 0.0), (0.0, 0.0), 1, 1);
        assert!(
            (p2.data[0] / p1.data[0] - 4.0).abs() < 0.01,
            "Constructive: 2 sources → 4× intensity"
        );
    }

    #[test]
    fn test_interference_all_non_negative() {
        let sources = [
            PointSource::new(-1e-3, 0.0, 1.0, 0.0),
            PointSource::new(1e-3, 0.0, 1.0, 0.0),
            PointSource::new(0.0, 1e-3, 0.5, 0.5),
        ];
        let pattern =
            interference_pattern(&sources, 0.5e-6, 1.0, (-5e-3, 5e-3), (-5e-3, 5e-3), 16, 16);
        for &v in &pattern.data {
            assert!(v >= 0.0, "Intensity must be non-negative");
        }
    }

    // ── Spectrum strip ───────────────────────────────────────────────────

    #[test]
    fn test_spectrum_strip_count() {
        let strip = spectrum_strip(100);
        assert_eq!(strip.len(), 100);
    }

    #[test]
    fn test_spectrum_strip_range() {
        let strip = spectrum_strip(50);
        for c in &strip {
            assert!((0.0..=1.0).contains(&c.r), "R out of range");
            assert!((0.0..=1.0).contains(&c.g), "G out of range");
            assert!((0.0..=1.0).contains(&c.b), "B out of range");
        }
    }

    #[test]
    fn test_spectrum_strip_red_end() {
        let strip = spectrum_strip(100);
        let last = &strip[strip.len() - 1];
        assert!(last.r > last.b, "Red end should be red-dominant");
    }

    #[test]
    fn test_spectrum_strip_blue_start() {
        let strip = spectrum_strip(100);
        // Sample a point clearly in the blue range (not the violet edge at 380nm)
        let blue_region = &strip[10]; // ~420nm
        assert!(
            blue_region.b > blue_region.r,
            "Blue region should be blue-dominant"
        );
    }

    #[test]
    fn test_spectrum_strip_custom_range() {
        let strip = spectrum_strip_range(500.0, 600.0, 20);
        assert_eq!(strip.len(), 20);
        // 500-600nm is green-yellow
        for c in &strip {
            assert!(c.g > 0.3, "500-600nm should have green: g={}", c.g);
        }
    }

    // ── PSF ──────────────────────────────────────────────────────────────

    #[test]
    fn test_psf_perfect_wavefront() {
        let n = 32;
        let radius = 8.0;
        let cx = n as f64 / 2.0;
        let cy = n as f64 / 2.0;
        let r2 = radius * radius;

        let mut pupil = vec![0.0; n * n];
        let wavefront = vec![0.0; n * n]; // perfect

        for row in 0..n {
            for col in 0..n {
                let dx = col as f64 - cx;
                let dy = row as f64 - cy;
                if dx * dx + dy * dy <= r2 {
                    pupil[row * n + col] = 1.0;
                }
            }
        }

        let psf = psf_from_wavefront(&pupil, &wavefront, n, n, 0.55e-6);
        assert!(psf.max_intensity() > 0.0);
        // Peak should be at center
        let center = psf.get(n / 2, n / 2);
        assert!(center > psf.get(0, 0), "PSF peak should be at center");
    }

    #[test]
    fn test_psf_diffraction_limited() {
        let psf = psf_diffraction_limited(32, 8.0);
        assert_eq!(psf.width, 32);
        assert!(psf.max_intensity() > 0.0);
    }

    #[test]
    fn test_psf_aberrated_broader() {
        let n = 32;
        let radius = 8.0;
        let cx = n as f64 / 2.0;
        let cy = n as f64 / 2.0;
        let r2 = radius * radius;

        let mut pupil = vec![0.0; n * n];
        let wfe_zero = vec![0.0; n * n];
        let mut wfe_aberr = vec![0.0; n * n];

        for row in 0..n {
            for col in 0..n {
                let dx = col as f64 - cx;
                let dy = row as f64 - cy;
                let rho2 = dx * dx + dy * dy;
                if rho2 <= r2 {
                    pupil[row * n + col] = 1.0;
                    // Spherical aberration: W = rho^4
                    wfe_aberr[row * n + col] = 0.5e-6 * (rho2 / r2) * (rho2 / r2);
                }
            }
        }

        let psf_perfect = psf_from_wavefront(&pupil, &wfe_zero, n, n, 0.55e-6);
        let psf_aberr = psf_from_wavefront(&pupil, &wfe_aberr, n, n, 0.55e-6);

        // Aberrated PSF should have lower peak (energy spread out)
        assert!(
            psf_aberr.max_intensity() < psf_perfect.max_intensity(),
            "Aberrated PSF should have lower peak"
        );
    }

    // ── Pattern2D ────────────────────────────────────────────────────────

    #[test]
    fn test_pattern_normalized() {
        let mut p = Pattern2D::new(4, 4);
        p.set(0, 0, 10.0);
        p.set(1, 1, 5.0);
        let norm = p.normalized();
        assert!((norm.get(0, 0) - 1.0).abs() < EPS);
        assert!((norm.get(1, 1) - 0.5).abs() < EPS);
    }

    #[test]
    fn test_pattern_normalized_zero() {
        let p = Pattern2D::new(4, 4);
        let norm = p.normalized();
        assert!(norm.max_intensity().abs() < EPS);
    }

    #[test]
    fn test_pattern_max_intensity() {
        let mut p = Pattern2D::new(4, 4);
        p.set(2, 3, 42.0);
        assert!((p.max_intensity() - 42.0).abs() < EPS);
    }
}
