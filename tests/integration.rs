//! Integration tests — cross-module optics scenarios.

#[test]
fn test_glass_window_reflection_and_transmission() {
    use prakash::ray::{Medium, deg_to_rad, fresnel_normal, fresnel_unpolarized, snell};

    let angle_i = deg_to_rad(30.0);
    let angle_t = snell(Medium::AIR.n, Medium::GLASS.n, angle_i).unwrap();

    // Some light reflects, some transmits
    let reflectance = fresnel_normal(Medium::AIR.n, Medium::GLASS.n);
    assert!(reflectance > 0.0 && reflectance < 1.0);
    assert!(angle_t < angle_i); // bends toward normal entering denser medium

    // Unpolarized reflectance at this angle should be > normal incidence
    let r_30 = fresnel_unpolarized(Medium::AIR.n, Medium::GLASS.n, angle_i).unwrap();
    assert!(r_30 >= reflectance - 0.01);
}

#[test]
fn test_sun_peak_wavelength_is_green() {
    use prakash::spectral::{wavelength_to_rgb, wien_peak};

    let peak_m = wien_peak(5778.0);
    let peak_nm = peak_m * 1e9;

    assert!(peak_nm > 480.0 && peak_nm < 530.0);

    let rgb = wavelength_to_rgb(peak_nm).unwrap();
    assert!(rgb.g > rgb.r);
    assert!(rgb.g > rgb.b);
}

#[test]
fn test_pbr_dielectric_vs_metal() {
    use prakash::pbr::{fresnel_schlick, ior_to_f0};

    let f0_glass = ior_to_f0(1.5);
    assert!(f0_glass < 0.1);
    assert!(fresnel_schlick(f0_glass, 0.0) > 0.9);

    let f0_metal = 0.95;
    assert!(fresnel_schlick(f0_metal, 1.0) > 0.9);
}

#[test]
fn test_polarizer_chain() {
    use prakash::wave::malus_law;
    use std::f64::consts::FRAC_PI_4;

    let i = malus_law(1.0, std::f64::consts::FRAC_PI_2);
    assert!(i < 1e-10);

    let after_first = malus_law(1.0, FRAC_PI_4);
    let after_second = malus_law(after_first, FRAC_PI_4);
    assert!((after_second - 0.25).abs() < 0.01);
}

#[test]
fn test_camera_lens_magnification() {
    use prakash::lens::{magnification, thin_lens_image_distance};

    let di = thin_lens_image_distance(50.0, 2000.0).unwrap();
    let m = magnification(2000.0, di);

    assert!(di > 0.0);
    assert!(m < 0.0);
    assert!(m.abs() < 1.0);
}

#[test]
fn test_ray_fresnel_matches_pbr_ior_to_f0() {
    use prakash::pbr::ior_to_f0;
    use prakash::ray::{Medium, fresnel_normal};

    // fresnel_normal and ior_to_f0 should agree for air→material
    let materials = [Medium::WATER, Medium::GLASS, Medium::DIAMOND];
    for m in &materials {
        let ray_r = fresnel_normal(1.0, m.n);
        let pbr_f0 = ior_to_f0(m.n);
        assert!(
            (ray_r - pbr_f0).abs() < 0.001,
            "Mismatch for {}: ray={ray_r}, pbr={pbr_f0}",
            m.name
        );
    }
}

#[test]
fn test_spectral_energy_conservation_intuition() {
    use prakash::spectral::{SPEED_OF_LIGHT, photon_energy_ev, wavelength_to_frequency};

    // Higher frequency → higher energy
    let freq_blue = wavelength_to_frequency(450.0);
    let freq_red = wavelength_to_frequency(700.0);
    assert!(freq_blue > freq_red);

    let e_blue = photon_energy_ev(450.0);
    let e_red = photon_energy_ev(700.0);
    assert!(e_blue > e_red);

    // Verify c = λ·f
    let lambda_m = 550e-9;
    let freq = wavelength_to_frequency(550.0);
    assert!((lambda_m * freq - SPEED_OF_LIGHT).abs() < 1.0);
}

#[test]
fn test_lens_system_two_element() {
    use prakash::lens::{
        LensType, classify_lens, combined_focal_length, magnification, thin_lens_image_distance,
    };

    // Two-element system: 50mm + 100mm
    let f_combined = combined_focal_length(50.0, 100.0).unwrap();
    assert!(f_combined > 0.0);
    assert!(f_combined < 50.0); // shorter than either element
    assert_eq!(classify_lens(f_combined), LensType::Converging);

    // Image through the combined system
    let di = thin_lens_image_distance(f_combined, 500.0).unwrap();
    let m = magnification(500.0, di);
    assert!(di > 0.0);
    assert!(m.abs() < 1.0); // diminished at this distance
}

#[test]
fn test_snell_through_multiple_media() {
    use prakash::ray::{Medium, deg_to_rad, snell};

    // Air → Glass → Water: trace a ray through two boundaries
    let angle_air = deg_to_rad(30.0);
    let angle_glass = snell(Medium::AIR.n, Medium::GLASS.n, angle_air).unwrap();
    let angle_water = snell(Medium::GLASS.n, Medium::WATER.n, angle_glass).unwrap();

    // Direct air → water should give same result
    let angle_water_direct = snell(Medium::AIR.n, Medium::WATER.n, angle_air).unwrap();
    assert!(
        (angle_water - angle_water_direct).abs() < 1e-6,
        "Snell's law should be transitive through media"
    );
}

#[test]
fn test_wave_and_spectral_consistency() {
    use prakash::spectral::wavelength_to_frequency;
    use prakash::wave::path_to_phase;

    // Phase from one full wavelength should be 2π
    let wl_nm = 550.0;
    let phase = path_to_phase(wl_nm, wl_nm);
    assert!((phase - 2.0 * std::f64::consts::PI).abs() < 1e-6);

    // Half wavelength → π
    let phase_half = path_to_phase(wl_nm / 2.0, wl_nm);
    assert!((phase_half - std::f64::consts::PI).abs() < 1e-6);

    // Frequency should be positive and finite
    let freq = wavelength_to_frequency(wl_nm);
    assert!(freq > 0.0);
    assert!(freq.is_finite());
}
