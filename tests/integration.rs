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
fn test_jones_to_stokes_and_mueller_chain() {
    use prakash::wave::{MuellerMatrix, Polarization, StokesVector};

    // Convert Jones horizontal to Stokes, pass through a horizontal polarizer
    let jones = Polarization::HORIZONTAL;
    let stokes = StokesVector::from(jones);
    let result = MuellerMatrix::POLARIZER_HORIZONTAL.apply(&stokes);

    // Should pass through undiminished (horizontal through horizontal polarizer)
    assert!(
        (result.s0 - 1.0).abs() < 0.01,
        "Horizontal through H polarizer should pass"
    );

    // Vertical through horizontal polarizer should be blocked
    let jones_v = Polarization::VERTICAL;
    let stokes_v = StokesVector::from(jones_v);
    let result_v = MuellerMatrix::POLARIZER_HORIZONTAL.apply(&stokes_v);
    assert!(
        result_v.s0.abs() < 0.01,
        "Vertical through H polarizer should be blocked"
    );
}

#[test]
fn test_spd_cow_illuminant_works_in_cri() {
    use prakash::spectral::{color_rendering_index, illuminant_d65};

    // illuminant_d65 now uses Cow::Borrowed — verify it still works correctly
    let d65 = illuminant_d65();
    let cri = color_rendering_index(&d65);
    assert!(
        cri > 90.0,
        "D65 CRI should still be high with Cow storage, got {cri}"
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
    let phase_half = path_to_phase(wl_nm, wl_nm / 2.0);
    assert!((phase_half - std::f64::consts::PI).abs() < 1e-6);

    // Frequency should be positive and finite
    let freq = wavelength_to_frequency(wl_nm);
    assert!(freq > 0.0);
    assert!(freq.is_finite());
}

// ── Bijli integration tests ─────────────────────────────────────────────────

#[cfg(feature = "bijli-backend")]
mod bijli_integration {
    const EPS: f64 = 1e-6;

    #[test]
    fn snell_delegates_correctly() {
        use prakash::ray::{deg_to_rad, snell};
        // Air → glass at 30° — verify the bijli-backed path gives correct physics
        let angle_t = snell(1.0, 1.52, deg_to_rad(30.0)).unwrap();
        assert!((angle_t.to_degrees() - 19.2).abs() < 0.5);
    }

    #[test]
    fn snell_tir_still_errors() {
        use prakash::ray::{deg_to_rad, snell};
        // Glass → air at 45° should still be TIR
        assert!(snell(1.52, 1.0, deg_to_rad(45.0)).is_err());
    }

    #[test]
    fn critical_angle_delegates_correctly() {
        use prakash::ray::{Medium, critical_angle};
        let ca = critical_angle(Medium::GLASS.n, Medium::AIR.n).unwrap();
        assert!((ca.to_degrees() - 41.1).abs() < 0.5);
    }

    #[test]
    fn critical_angle_rejects_n1_le_n2() {
        use prakash::ray::critical_angle;
        assert!(critical_angle(1.0, 1.5).is_err());
    }

    #[test]
    fn brewster_angle_delegates_correctly() {
        use prakash::ray::brewster_angle;
        let ba = brewster_angle(1.0, 1.5);
        assert!((ba.to_degrees() - 56.3).abs() < 0.5);
    }

    #[test]
    fn fresnel_normal_delegates_correctly() {
        use prakash::ray::{Medium, fresnel_normal};
        let r = fresnel_normal(Medium::AIR.n, Medium::GLASS.n);
        assert!((r - 0.04).abs() < 0.01);
    }

    #[test]
    fn fresnel_unpolarized_delegates_correctly() {
        use prakash::ray::{deg_to_rad, fresnel_unpolarized};
        let r = fresnel_unpolarized(1.0, 1.5, deg_to_rad(30.0)).unwrap();
        assert!(r > 0.0 && r < 1.0);
    }

    #[test]
    fn speed_of_light_matches_bijli() {
        assert!(
            (prakash::spectral::SPEED_OF_LIGHT - bijli::field::SPEED_OF_LIGHT).abs() < EPS,
            "prakash and bijli should share the same c"
        );
    }

    #[test]
    fn stokes_roundtrip() {
        use prakash::wave::StokesVector;
        let p = StokesVector::horizontal(1.0);
        let b: bijli::polarization::StokesVector = p.into();
        let back: StokesVector = b.into();
        assert!((back.s0 - p.s0).abs() < EPS);
        assert!((back.s1 - p.s1).abs() < EPS);
        assert!((back.s2 - p.s2).abs() < EPS);
        assert!((back.s3 - p.s3).abs() < EPS);
    }

    #[test]
    fn mueller_roundtrip() {
        use prakash::wave::MuellerMatrix;
        let p = MuellerMatrix::POLARIZER_HORIZONTAL;
        let b: bijli::polarization::MuellerMatrix = p.into();
        let back: MuellerMatrix = b.into();
        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    (back.m[i][j] - p.m[i][j]).abs() < EPS,
                    "Mueller roundtrip mismatch at [{i}][{j}]"
                );
            }
        }
    }

    #[test]
    fn polarization_to_jones() {
        use prakash::wave::{JonesVector, Polarization};
        let p = Polarization::HORIZONTAL;
        let j: JonesVector = p.into();
        // Horizontal should have x=1, y=0
        assert!((j.intensity() - 1.0).abs() < EPS);
    }

    #[test]
    fn gaussian_beam_available() {
        use prakash::wave::GaussianBeam;
        let beam = GaussianBeam::new(632.8e-9, 1e-3).unwrap();
        assert!(beam.divergence() > 0.0);
        assert!(beam.rayleigh_range > 0.0);
    }

    #[test]
    fn abcd_thin_lens() {
        use prakash::wave::AbcdMatrix;
        let lens = AbcdMatrix::thin_lens(0.1).unwrap();
        // Thin lens should have a=1, d=1, c=-1/f
        assert!((lens.a - 1.0).abs() < EPS);
        assert!((lens.d - 1.0).abs() < EPS);
        assert!((lens.c + 10.0).abs() < EPS); // -1/0.1 = -10
    }

    #[test]
    fn medium_permittivity() {
        use prakash::ray::Medium;
        // n=1.5 → ε_r = 2.25
        assert!((Medium::GLASS.permittivity() - 2.3104).abs() < 0.01); // 1.52² = 2.3104
    }
}

// ── Property-based tests ────────────────────────────────────────────────────

mod property_tests {
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn snell_reversible(angle_deg in 1.0f64..40.0) {
            use prakash::ray::snell;
            let angle = angle_deg.to_radians();
            let angle_t = snell(1.0, 1.52, angle).unwrap();
            let angle_back = snell(1.52, 1.0, angle_t).unwrap();
            prop_assert!((angle_back - angle).abs() < 1e-10,
                "Snell should be reversible: {angle} → {angle_t} → {angle_back}");
        }

        #[test]
        fn fresnel_normal_symmetric(n1 in 1.0f64..3.0, n2 in 1.0f64..3.0) {
            use prakash::ray::fresnel_normal;
            let r1 = fresnel_normal(n1, n2);
            let r2 = fresnel_normal(n2, n1);
            prop_assert!((r1 - r2).abs() < 1e-10,
                "Fresnel normal should be symmetric: R({n1},{n2})={r1} vs R({n2},{n1})={r2}");
        }

        #[test]
        fn fresnel_normal_range(n1 in 1.0f64..4.0, n2 in 1.0f64..4.0) {
            use prakash::ray::fresnel_normal;
            let r = fresnel_normal(n1, n2);
            prop_assert!((0.0..=1.0).contains(&r), "Fresnel in [0,1]: R={r}");
        }

        #[test]
        fn wavelength_frequency_roundtrip(nm in 380.0f64..780.0) {
            use prakash::spectral::{wavelength_to_frequency, frequency_to_wavelength};
            let freq = wavelength_to_frequency(nm);
            let back = frequency_to_wavelength(freq);
            prop_assert!((back - nm).abs() < 1e-6,
                "Roundtrip failed: {nm} → {freq} → {back}");
        }

        #[test]
        fn planck_always_positive(wl_nm in 200.0f64..2000.0, temp in 100.0f64..50000.0) {
            use prakash::spectral::planck_radiance;
            let r = planck_radiance(wl_nm * 1e-9, temp);
            prop_assert!(r >= 0.0 && r.is_finite(),
                "Planck should be non-negative and finite: {r} at {wl_nm}nm, {temp}K");
        }

        #[test]
        fn beer_lambert_monotonic_decay(alpha in 0.01f64..10.0, d in 0.0f64..100.0) {
            use prakash::ray::beer_lambert;
            let i = beer_lambert(1.0, alpha, d);
            prop_assert!((0.0..=1.0).contains(&i),
                "Beer-Lambert in [0,1]: I={i} for α={alpha}, d={d}");
        }

        #[test]
        fn zernike_outside_pupil_is_zero(rho in 1.01f64..5.0, theta in 0.0f64..std::f64::consts::TAU) {
            use prakash::wave::zernike::zernike;
            let z = zernike(2, 0, rho, theta);
            prop_assert!(z.abs() < 1e-10,
                "Zernike outside pupil should be 0: Z={z} at ρ={rho}");
        }

        #[test]
        fn complex_fresnel_dielectric_matches_real(
            n in 1.0f64..3.0,
            angle_deg in 0.0f64..85.0
        ) {
            use prakash::ray::{ComplexMedium, fresnel_unpolarized, fresnel_unpolarized_complex};
            let m = ComplexMedium::dielectric(n, "test");
            let angle = angle_deg.to_radians();
            let r_complex = fresnel_unpolarized_complex(1.0, &m, angle);
            if let Ok(r_real) = fresnel_unpolarized(1.0, n, angle) {
                prop_assert!((r_complex - r_real).abs() < 0.01,
                    "Complex should match real for dielectric: {r_complex} vs {r_real}");
            }
        }
    }
}
