//! Integration tests — cross-module optics scenarios.

#[test]
fn test_glass_window_reflection_and_transmission() {
    use prakash::ray::{Medium, fresnel_normal, snell, deg_to_rad};

    // Light hitting a glass window at 30°
    let angle_i = deg_to_rad(30.0);
    let angle_t = snell(Medium::AIR.n, Medium::GLASS.n, angle_i).unwrap();

    // Some light reflects, some transmits
    let reflectance = fresnel_normal(Medium::AIR.n, Medium::GLASS.n);
    assert!(reflectance > 0.0 && reflectance < 1.0);
    assert!(angle_t < angle_i); // bends toward normal entering denser medium
}

#[test]
fn test_sun_peak_wavelength_is_green() {
    use prakash::spectral::{wien_peak, wavelength_to_rgb};

    let peak_m = wien_peak(5778.0); // Sun's surface temperature
    let peak_nm = peak_m * 1e9;

    // Sun peaks in green (~502nm)
    assert!(peak_nm > 480.0 && peak_nm < 530.0);

    // That wavelength should be greenish
    let rgb = wavelength_to_rgb(peak_nm).unwrap();
    assert!(rgb.g > rgb.r);
    assert!(rgb.g > rgb.b);
}

#[test]
fn test_pbr_dielectric_vs_metal() {
    use prakash::pbr::{fresnel_schlick, ior_to_f0};

    // Dielectric (glass): low F0, high Fresnel at grazing
    let f0_glass = ior_to_f0(1.5);
    assert!(f0_glass < 0.1);
    assert!(fresnel_schlick(f0_glass, 0.0) > 0.9); // grazing → high

    // Metal (gold-like): high F0, stays high everywhere
    let f0_metal = 0.95;
    assert!(fresnel_schlick(f0_metal, 1.0) > 0.9); // even at normal
}

#[test]
fn test_polarizer_chain() {
    use prakash::wave::malus_law;
    use std::f64::consts::FRAC_PI_4;

    // Two crossed polarizers: 0% transmission
    let i = malus_law(1.0, std::f64::consts::FRAC_PI_2);
    assert!(i < 1e-10);

    // Insert 45° polarizer between them: some light gets through
    let after_first = malus_law(1.0, FRAC_PI_4); // 50%
    let after_second = malus_law(after_first, FRAC_PI_4); // 25% of original
    assert!((after_second - 0.25).abs() < 0.01);
}

#[test]
fn test_camera_lens_magnification() {
    use prakash::lens::{thin_lens_image_distance, magnification};

    // 50mm portrait lens, subject at 2 meters
    let di = thin_lens_image_distance(50.0, 2000.0).unwrap();
    let m = magnification(2000.0, di);

    // Should produce small inverted real image
    assert!(di > 0.0); // real image
    assert!(m < 0.0); // inverted
    assert!(m.abs() < 1.0); // diminished
}
