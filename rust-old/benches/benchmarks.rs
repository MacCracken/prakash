use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use std::f64::consts::FRAC_PI_6;

fn bench_ray(c: &mut Criterion) {
    use prakash::ray::*;
    let mut group = c.benchmark_group("ray");

    group.bench_function("snell_air_glass", |b| {
        b.iter(|| snell(black_box(1.0), black_box(1.52), black_box(FRAC_PI_6)))
    });
    group.bench_function("snell_tir_check", |b| {
        b.iter(|| snell(black_box(1.52), black_box(1.0), black_box(0.8)))
    });
    group.bench_function("fresnel_normal", |b| {
        b.iter(|| fresnel_normal(black_box(1.0), black_box(1.52)))
    });
    group.bench_function("fresnel_s", |b| {
        b.iter(|| {
            fresnel_s(
                black_box(1.0),
                black_box(1.52),
                black_box(0.866),
                black_box(0.94),
            )
        })
    });
    group.bench_function("fresnel_p", |b| {
        b.iter(|| {
            fresnel_p(
                black_box(1.0),
                black_box(1.52),
                black_box(0.866),
                black_box(0.94),
            )
        })
    });
    group.bench_function("fresnel_unpolarized", |b| {
        b.iter(|| fresnel_unpolarized(black_box(1.0), black_box(1.52), black_box(0.5)))
    });
    group.bench_function("critical_angle", |b| {
        b.iter(|| critical_angle(black_box(1.52), black_box(1.0)))
    });
    group.bench_function("reflect_2d", |b| {
        b.iter(|| reflect_2d(black_box([0.707, -0.707]), black_box([0.0, 1.0])))
    });
    group.bench_function("reflect_3d", |b| {
        b.iter(|| reflect_3d(black_box([0.707, -0.707, 0.0]), black_box([0.0, 1.0, 0.0])))
    });
    group.bench_function("beer_lambert", |b| {
        b.iter(|| beer_lambert(black_box(1.0), black_box(0.1), black_box(5.0)))
    });
    group.bench_function("brewster_angle", |b| {
        b.iter(|| brewster_angle(black_box(1.0), black_box(1.52)))
    });
    group.bench_function("refract_3d", |b| {
        let dir = [0.5f64.sqrt(), 0.0, -(0.5f64.sqrt())];
        b.iter(|| {
            refract_3d(
                black_box(dir),
                black_box([0.0, 0.0, 1.0]),
                black_box(1.0),
                black_box(1.52),
            )
        })
    });
    group.bench_function("snell_3d", |b| {
        let dir = [0.5f64.sqrt(), 0.0, -(0.5f64.sqrt())];
        b.iter(|| {
            snell_3d(
                black_box(dir),
                black_box([0.0, 0.0, 1.0]),
                black_box(1.0),
                black_box(1.52),
            )
        })
    });
    group.bench_function("sellmeier_n_at", |b| {
        b.iter(|| SellmeierCoefficients::BK7.n_at(black_box(0.5876)))
    });
    group.bench_function("cauchy_n_at", |b| {
        b.iter(|| CauchyCoefficients::BK7.n_at(black_box(0.5876)))
    });
    group.bench_function("abbe_number", |b| {
        b.iter(|| abbe_number(black_box(&SellmeierCoefficients::BK7)))
    });
    group.bench_function("prism_deviation", |b| {
        b.iter(|| prism_deviation(black_box(std::f64::consts::FRAC_PI_3), black_box(1.5168)))
    });
    group.bench_function("trace_plane_surface", |b| {
        let ray = TraceRay {
            position: [0.0, 0.0, 0.0],
            direction: [0.0, 0.0, 1.0],
            n: 1.0,
        };
        let surface = OpticalSurface {
            shape: SurfaceShape::Plane,
            z_position: 10.0,
            n_after: 1.5,
            aperture_radius: 25.0,
        };
        b.iter(|| trace_surface(black_box(&ray), black_box(&surface)))
    });
    group.bench_function("trace_sphere_surface", |b| {
        let ray = TraceRay {
            position: [2.0, 0.0, -50.0],
            direction: [0.0, 0.0, 1.0],
            n: 1.0,
        };
        let surface = OpticalSurface {
            shape: SurfaceShape::Sphere { radius: 50.0 },
            z_position: 0.0,
            n_after: 1.5,
            aperture_radius: 25.0,
        };
        b.iter(|| trace_surface(black_box(&ray), black_box(&surface)))
    });
    group.bench_function("trace_sequential_2", |b| {
        let ray = TraceRay {
            position: [2.0, 0.0, -50.0],
            direction: [0.0, 0.0, 1.0],
            n: 1.0,
        };
        let surfaces = [
            OpticalSurface {
                shape: SurfaceShape::Sphere { radius: 50.0 },
                z_position: 0.0,
                n_after: 1.5,
                aperture_radius: 25.0,
            },
            OpticalSurface {
                shape: SurfaceShape::Sphere { radius: -50.0 },
                z_position: 5.0,
                n_after: 1.0,
                aperture_radius: 25.0,
            },
        ];
        b.iter(|| trace_sequential(black_box(&ray), black_box(&surfaces)))
    });
    group.bench_function("prism_dispersion", |b| {
        let s = SellmeierCoefficients::BK7;
        b.iter(|| {
            prism_dispersion(
                black_box(std::f64::consts::FRAC_PI_3),
                black_box(&s),
                black_box(0.55),
            )
        })
    });
    group.bench_function("prism_angular_spread", |b| {
        let s = SellmeierCoefficients::BK7;
        b.iter(|| {
            prism_angular_spread(
                black_box(std::f64::consts::FRAC_PI_3),
                black_box(&s),
                black_box(0.486),
                black_box(0.656),
            )
        })
    });

    // ── Simulation primitives ────────────────────────────────────────────
    let biconvex = [
        OpticalSurface {
            shape: SurfaceShape::Sphere { radius: 100.0 },
            z_position: 0.0,
            n_after: 1.5,
            aperture_radius: 25.0,
        },
        OpticalSurface {
            shape: SurfaceShape::Sphere { radius: -100.0 },
            z_position: 5.0,
            n_after: 1.0,
            aperture_radius: 25.0,
        },
    ];
    let on_axis = TraceRay {
        position: [0.0, 0.0, -50.0],
        direction: [0.0, 0.0, 1.0],
        n: 1.0,
    };
    group.bench_function("trace_recursive_2", |b| {
        b.iter(|| {
            trace_recursive(
                black_box(&on_axis),
                black_box(&biconvex),
                black_box(&TraceConfig::default()),
            )
        })
    });
    group.bench_function("ray_fan_meridional_21", |b| {
        b.iter(|| {
            ray_fan_meridional(
                black_box(10.0),
                black_box(21),
                black_box(0.0),
                black_box(-50.0),
            )
        })
    });
    group.bench_function("ray_bundle_3x8", |b| {
        b.iter(|| {
            ray_bundle(
                black_box(10.0),
                black_box(3),
                black_box(8),
                black_box(0.0),
                black_box(-50.0),
            )
        })
    });
    group.bench_function("spot_diagram_2x6", |b| {
        b.iter(|| {
            spot_diagram(
                black_box(&biconvex),
                black_box(5.0),
                black_box(2),
                black_box(6),
                black_box(0.0),
                black_box(-50.0),
                black_box(200.0),
            )
        })
    });
    group.bench_function("optical_path_length", |b| {
        b.iter(|| optical_path_length(black_box(&on_axis), black_box(&biconvex), black_box(200.0)))
    });
    let prescription = Prescription::new("bench")
        .add_surface(100.0, 5.0, 1.5, 25.0)
        .add_surface(-100.0, 0.0, 1.0, 25.0);
    group.bench_function("paraxial_trace_2", |b| {
        let surfaces = prescription.to_paraxial_surfaces();
        let ray = ParaxialRay::marginal(10.0);
        b.iter(|| paraxial_trace(black_box(&ray), black_box(&surfaces)))
    });
    group.bench_function("find_system_properties", |b| {
        b.iter(|| find_system_properties(black_box(&prescription)))
    });
    group.bench_function("fresnel_normal_complex_gold", |b| {
        b.iter(|| fresnel_normal_complex(black_box(1.0), black_box(&ComplexMedium::GOLD_550NM)))
    });
    group.bench_function("fresnel_s_complex_gold", |b| {
        b.iter(|| {
            fresnel_s_complex(
                black_box(1.0),
                black_box(&ComplexMedium::GOLD_550NM),
                black_box(0.5),
            )
        })
    });
    group.bench_function("fiber_na", |b| {
        b.iter(|| fiber::fiber_na(black_box(1.4682), black_box(1.4629)))
    });
    group.bench_function("fiber_v_number", |b| {
        b.iter(|| fiber::v_number(black_box(4.1e-6), black_box(0.12), black_box(1.55e-6)))
    });
    group.bench_function("fiber_mfd", |b| {
        b.iter(|| fiber::mode_field_diameter(black_box(4.1e-6), black_box(2.0)))
    });
    group.bench_function("fiber_coupling", |b| {
        b.iter(|| fiber::coupling_efficiency_gaussian(black_box(5e-6), black_box(5.2e-6)))
    });
    group.bench_function("partial_dispersion_bk7", |b| {
        b.iter(|| partial_dispersion(black_box(&SellmeierCoefficients::BK7)))
    });
    group.bench_function("schott_n_at", |b| {
        b.iter(|| SchottCoefficients::BK7.n_at(black_box(0.5876)))
    });

    group.finish();
}

fn bench_spectral(c: &mut Criterion) {
    use prakash::spectral::*;
    let mut group = c.benchmark_group("spectral");

    group.bench_function("wavelength_to_rgb", |b| {
        b.iter(|| wavelength_to_rgb(black_box(550.0)))
    });
    group.bench_function("wavelength_to_rgb_edge", |b| {
        b.iter(|| wavelength_to_rgb(black_box(400.0)))
    });
    group.bench_function("planck_radiance", |b| {
        b.iter(|| planck_radiance(black_box(500e-9), black_box(5778.0)))
    });
    group.bench_function("wien_peak", |b| b.iter(|| wien_peak(black_box(5778.0))));
    group.bench_function("color_temp_to_rgb", |b| {
        b.iter(|| color_temperature_to_rgb(black_box(6500.0)))
    });
    group.bench_function("color_temp_warm", |b| {
        b.iter(|| color_temperature_to_rgb(black_box(2700.0)))
    });
    group.bench_function("photon_energy_ev", |b| {
        b.iter(|| photon_energy_ev(black_box(550.0)))
    });
    group.bench_function("wavelength_to_frequency", |b| {
        b.iter(|| wavelength_to_frequency(black_box(550.0)))
    });
    group.bench_function("frequency_to_wavelength", |b| {
        b.iter(|| frequency_to_wavelength(black_box(5.45e14)))
    });
    group.bench_function("photon_energy", |b| {
        b.iter(|| photon_energy(black_box(550.0)))
    });
    group.bench_function("rgb_luminance", |b| {
        let rgb = Rgb::new(0.5, 0.7, 0.3);
        b.iter(|| black_box(rgb).luminance())
    });
    group.bench_function("rgb_to_u8", |b| {
        let rgb = Rgb::new(0.5, 0.7, 0.3);
        b.iter(|| black_box(rgb).to_u8())
    });
    group.bench_function("cie_cmf_at", |b| b.iter(|| cie_cmf_at(black_box(555.0))));
    group.bench_function("cie_cmf_at_interp", |b| {
        b.iter(|| cie_cmf_at(black_box(557.5)))
    });
    group.bench_function("xyz_to_xyy", |b| {
        let xyz = Xyz::D65_WHITE;
        b.iter(|| black_box(xyz).to_xyy())
    });
    group.bench_function("xyz_to_srgb", |b| {
        let xyz = Xyz::D65_WHITE;
        b.iter(|| black_box(xyz).to_srgb())
    });
    group.bench_function("xyz_to_linear_srgb", |b| {
        let xyz = Xyz::D65_WHITE;
        b.iter(|| black_box(xyz).to_linear_srgb())
    });
    group.bench_function("linear_srgb_to_xyz", |b| {
        let rgb = Rgb::new(0.5, 0.7, 0.3);
        b.iter(|| linear_srgb_to_xyz(black_box(&rgb)))
    });
    group.bench_function("srgb_gamma", |b| {
        b.iter(|| linear_to_srgb_gamma(black_box(0.5)))
    });
    group.bench_function("srgb_gamma_inv", |b| {
        b.iter(|| srgb_gamma_to_linear(black_box(0.5)))
    });
    group.bench_function("cct_from_xy", |b| {
        b.iter(|| cct_from_xy(black_box(0.3127), black_box(0.3290)))
    });
    group.bench_function("illuminant_d65_cow", |b| b.iter(illuminant_d65));
    group.bench_function("spd_blackbody", |b| {
        b.iter(|| Spd::blackbody(black_box(5778.0)))
    });
    group.bench_function("spd_to_xyz", |b| {
        let spd = Spd::blackbody(5778.0);
        b.iter(|| black_box(&spd).to_xyz())
    });
    group.bench_function("spd_to_srgb", |b| {
        let spd = Spd::blackbody(5778.0);
        b.iter(|| black_box(&spd).to_srgb())
    });
    group.bench_function("cri", |b| {
        let spd = illuminant_d65();
        b.iter(|| color_rendering_index(black_box(&spd)))
    });
    group.bench_function("luminous_flux_d65", |b| {
        let spd = illuminant_d65();
        b.iter(|| luminous_flux(black_box(&spd)))
    });
    group.bench_function("luminous_efficacy_d65", |b| {
        let spd = illuminant_d65();
        b.iter(|| luminous_efficacy(black_box(&spd)))
    });
    group.bench_function("spd_to_xyz_2015_2deg", |b| {
        let spd = illuminant_d65();
        b.iter(|| black_box(&spd).to_xyz_observer(black_box(Observer::Cie2015_2deg)))
    });

    group.finish();
}

fn bench_wave(c: &mut Criterion) {
    use prakash::wave::*;
    let mut group = c.benchmark_group("wave");

    group.bench_function("interference", |b| {
        b.iter(|| interference_intensity(black_box(1.0), black_box(1.0), black_box(0.5)))
    });
    group.bench_function("single_slit", |b| {
        b.iter(|| {
            single_slit_intensity(
                black_box(500e-9),
                black_box(1e-3),
                black_box(0.01),
                black_box(1.0),
            )
        })
    });
    group.bench_function("double_slit", |b| {
        b.iter(|| {
            double_slit_intensity(
                black_box(500e-9),
                black_box(0.1e-3),
                black_box(0.5e-3),
                black_box(0.01),
                black_box(1.0),
            )
        })
    });
    group.bench_function("malus_law", |b| {
        b.iter(|| malus_law(black_box(1.0), black_box(0.785)))
    });
    group.bench_function("thin_film", |b| {
        b.iter(|| thin_film_reflectance(black_box(550.0), black_box(100.0), black_box(1.5)))
    });
    group.bench_function("path_to_phase", |b| {
        b.iter(|| path_to_phase(black_box(500.0), black_box(500.0)))
    });
    group.bench_function("grating_maxima_3", |b| {
        b.iter(|| grating_maxima(black_box(500e-9), black_box(1e-6), black_box(3)))
    });
    group.bench_function("polarization_intensity", |b| {
        let p = Polarization::circular_right();
        b.iter(|| black_box(p).intensity())
    });
    group.bench_function("polarization_through", |b| {
        let p = Polarization::HORIZONTAL;
        b.iter(|| black_box(p).through_polarizer(black_box(0.5)))
    });
    group.bench_function("is_constructive", |b| {
        b.iter(|| is_constructive(black_box(500.0), black_box(500.0)))
    });
    group.bench_function("is_destructive", |b| {
        b.iter(|| is_destructive(black_box(250.0), black_box(500.0)))
    });
    group.bench_function("coherence_length", |b| {
        b.iter(|| coherence_length(black_box(550.0), black_box(1.0)))
    });
    group.bench_function("coherence_time", |b| {
        b.iter(|| coherence_time(black_box(550e-9), black_box(1e-9)))
    });
    group.bench_function("bessel_j1", |b| b.iter(|| bessel_j1(black_box(3.0))));
    group.bench_function("bessel_j1_large", |b| b.iter(|| bessel_j1(black_box(15.0))));
    group.bench_function("airy_pattern", |b| {
        b.iter(|| {
            airy_pattern(
                black_box(550e-9),
                black_box(10e-3),
                black_box(1e-5),
                black_box(1.0),
            )
        })
    });
    group.bench_function("fp_transmittance", |b| {
        b.iter(|| {
            fabry_perot_transmittance(
                black_box(550e-9),
                black_box(1e-3),
                black_box(1.0),
                black_box(0.0),
                black_box(0.9),
            )
        })
    });
    group.bench_function("fp_finesse", |b| {
        b.iter(|| fabry_perot_finesse(black_box(0.9)))
    });
    group.bench_function("fp_resolving_power", |b| {
        b.iter(|| {
            fabry_perot_resolving_power(
                black_box(550e-9),
                black_box(1e-3),
                black_box(1.0),
                black_box(0.9),
            )
        })
    });
    group.bench_function("stokes_dop", |b| {
        let s = StokesVector::new(1.0, 0.5, 0.3, 0.1);
        b.iter(|| black_box(s).degree_of_polarization())
    });
    group.bench_function("mueller_apply", |b| {
        let m = MuellerMatrix::POLARIZER_HORIZONTAL;
        let s = StokesVector::unpolarized(1.0);
        b.iter(|| black_box(m).apply(black_box(&s)))
    });
    group.bench_function("mueller_multiply", |b| {
        let a = MuellerMatrix::POLARIZER_HORIZONTAL;
        let b_mat = MuellerMatrix::rotation(0.5);
        b.iter(|| black_box(a).multiply(black_box(&b_mat)))
    });
    group.bench_function("jones_to_stokes", |b| {
        let p = Polarization::circular_right();
        b.iter(|| StokesVector::from(black_box(p)))
    });
    group.bench_function("mueller_polarizer", |b| {
        b.iter(|| MuellerMatrix::polarizer(black_box(0.785)))
    });
    group.bench_function("mueller_retarder", |b| {
        b.iter(|| MuellerMatrix::retarder(black_box(std::f64::consts::FRAC_PI_2)))
    });
    group.bench_function("mueller_chain_3", |b| {
        let s = StokesVector::unpolarized(1.0);
        let elements = [
            MuellerMatrix::POLARIZER_HORIZONTAL,
            MuellerMatrix::QUARTER_WAVE_HORIZONTAL,
            MuellerMatrix::POLARIZER_VERTICAL,
        ];
        b.iter(|| mueller_chain(black_box(&s), black_box(&elements)))
    });
    group.bench_function("birefringent_retardation", |b| {
        let mat = BirefringentMaterial::QUARTZ;
        b.iter(|| mat.retardation(black_box(100.0), black_box(550.0)))
    });
    group.bench_function("fraunhofer_rect", |b| {
        b.iter(|| {
            fraunhofer_rect(
                black_box(550e-9),
                black_box(1e-3),
                black_box(1e-3),
                black_box(0.001),
                black_box(0.001),
                black_box(1.0),
            )
        })
    });
    group.bench_function("fresnel_integral_c", |b| {
        b.iter(|| fresnel_integral_c(black_box(2.0)))
    });
    group.bench_function("fresnel_integral_s", |b| {
        b.iter(|| fresnel_integral_s(black_box(2.0)))
    });
    group.bench_function("fresnel_edge", |b| {
        b.iter(|| fresnel_edge_intensity(black_box(1.0)))
    });
    group.bench_function("huygens_fresnel_50", |b| {
        let aperture: Vec<(f64, f64)> = (0..50).map(|i| ((i as f64 - 25.0) * 1e-5, 1.0)).collect();
        b.iter(|| {
            huygens_fresnel_1d(
                black_box(&aperture),
                black_box(550e-9),
                black_box(0.1),
                black_box(0.0),
            )
        })
    });
    group.bench_function("coating_reflectance", |b| {
        b.iter(|| {
            coating_reflectance(
                black_box(550.0),
                black_box(1.0),
                black_box(1.38),
                black_box(1.52),
                black_box(99.6),
            )
        })
    });
    group.bench_function("multilayer_2", |b| {
        let layers = [(1.38, 99.6), (2.1, 65.5)];
        b.iter(|| {
            multilayer_reflectance(
                black_box(550.0),
                black_box(1.0),
                black_box(1.52),
                black_box(&layers),
            )
        })
    });

    // ── Pattern computation ─────────────────────────────────────────────
    group.bench_function("diffraction_2d_16", |b| {
        let aperture = vec![1.0; 16 * 16];
        b.iter(|| diffraction_pattern_2d(black_box(&aperture), black_box(16), black_box(16)))
    });
    group.bench_function("diffraction_circular_32", |b| {
        b.iter(|| diffraction_pattern_circular(black_box(32), black_box(8.0)))
    });
    group.bench_function("interference_2src_16x16", |b| {
        let sources = [
            PointSource::new(-0.5e-3, 0.0, 1.0, 0.0),
            PointSource::new(0.5e-3, 0.0, 1.0, 0.0),
        ];
        b.iter(|| {
            interference_pattern(
                black_box(&sources),
                black_box(0.5e-6),
                black_box(1.0),
                black_box((-5e-3, 5e-3)),
                black_box((-5e-3, 5e-3)),
                black_box(16),
                black_box(16),
            )
        })
    });
    group.bench_function("spectrum_strip_100", |b| {
        b.iter(|| spectrum_strip(black_box(100)))
    });
    group.bench_function("psf_32", |b| {
        b.iter(|| psf_diffraction_limited(black_box(32), black_box(8.0)))
    });
    group.bench_function("zernike_eval", |b| {
        b.iter(|| zernike::zernike(black_box(4), black_box(0), black_box(0.7), black_box(0.5)))
    });
    group.bench_function("zernike_wavefront_eval", |b| {
        let wf = zernike::defocus(0.1);
        b.iter(|| wf.evaluate(black_box(0.7), black_box(0.5)))
    });
    group.bench_function("zernike_rms", |b| {
        let wf = zernike::ZernikeWavefront::new(vec![
            0.0, 0.01, 0.02, 0.1, 0.05, 0.03, 0.02, 0.01, 0.005, 0.002, 0.05,
        ]);
        b.iter(|| black_box(&wf).rms_error())
    });
    group.bench_function("multilayer_rt_oblique", |b| {
        b.iter(|| {
            multilayer_rt(
                black_box(550.0),
                black_box(0.5),
                black_box(1.0),
                black_box(1.52),
                black_box(&[(1.38, 99.6), (2.1, 65.5)]),
            )
        })
    });
    group.bench_function("king_factor", |b| {
        use prakash::atmosphere::king_factor;
        b.iter(|| king_factor(black_box(550e-9)))
    });

    group.finish();
}

fn bench_lens(c: &mut Criterion) {
    use prakash::lens::*;
    let mut group = c.benchmark_group("lens");

    group.bench_function("thin_lens", |b| {
        b.iter(|| thin_lens_image_distance(black_box(50.0), black_box(2000.0)))
    });
    group.bench_function("magnification", |b| {
        b.iter(|| magnification(black_box(2000.0), black_box(51.28)))
    });
    group.bench_function("lensmaker", |b| {
        b.iter(|| lensmaker_focal_length(black_box(1.5), black_box(100.0), black_box(-100.0)))
    });
    group.bench_function("depth_of_field", |b| {
        b.iter(|| {
            depth_of_field(
                black_box(50.0),
                black_box(2.8),
                black_box(0.03),
                black_box(2000.0),
            )
        })
    });
    group.bench_function("optical_power", |b| {
        b.iter(|| optical_power(black_box(0.05)).unwrap())
    });
    group.bench_function("combined_focal", |b| {
        b.iter(|| combined_focal_length(black_box(50.0), black_box(100.0)).unwrap())
    });
    group.bench_function("mirror_focal", |b| {
        b.iter(|| mirror_focal_length(black_box(200.0)))
    });
    group.bench_function("thick_lens", |b| {
        b.iter(|| {
            thick_lens_focal_length(
                black_box(1.5),
                black_box(100.0),
                black_box(-100.0),
                black_box(10.0),
            )
        })
    });
    group.bench_function("cardinal_points", |b| {
        b.iter(|| {
            cardinal_points(
                black_box(1.5),
                black_box(100.0),
                black_box(-100.0),
                black_box(10.0),
            )
        })
    });
    group.bench_function("f_number", |b| {
        b.iter(|| f_number(black_box(50.0), black_box(25.0)))
    });
    group.bench_function("numerical_aperture", |b| {
        b.iter(|| numerical_aperture(black_box(1.0), black_box(0.5)))
    });
    group.bench_function("airy_disk_radius", |b| {
        b.iter(|| airy_disk_radius(black_box(0.00055), black_box(2.8)))
    });
    group.bench_function("field_of_view", |b| {
        b.iter(|| field_of_view(black_box(36.0), black_box(50.0)))
    });
    group.bench_function("mtf_diffraction", |b| {
        b.iter(|| mtf_diffraction_limited(black_box(500.0), black_box(1000.0)))
    });
    group.bench_function("seidel_coefficients", |b| {
        b.iter(|| {
            seidel_coefficients(
                black_box(1.5),
                black_box(100.0),
                black_box(0.0),
                black_box(-1.0),
            )
        })
    });
    group.bench_function("separated_lenses", |b| {
        b.iter(|| separated_lenses_focal_length(black_box(50.0), black_box(-80.0), black_box(30.0)))
    });
    group.bench_function("aperture_from_f_number", |b| {
        b.iter(|| aperture_from_f_number(black_box(50.0), black_box(2.8)))
    });
    group.bench_function("na_from_f_number", |b| {
        b.iter(|| na_from_f_number(black_box(2.8)))
    });
    group.bench_function("diffraction_limit", |b| {
        b.iter(|| diffraction_limit(black_box(550e-6), black_box(10.0)))
    });
    group.bench_function("field_of_view_diagonal", |b| {
        b.iter(|| field_of_view_diagonal(black_box(36.0), black_box(24.0), black_box(50.0)))
    });
    group.bench_function("mtf_cutoff_frequency", |b| {
        b.iter(|| mtf_cutoff_frequency(black_box(0.00055), black_box(2.8)))
    });
    group.bench_function("shape_factor", |b| {
        b.iter(|| shape_factor(black_box(100.0), black_box(-100.0)))
    });
    group.bench_function("conjugate_factor", |b| {
        b.iter(|| conjugate_factor(black_box(200.0), black_box(100.0)))
    });
    group.bench_function("lsa", |b| {
        b.iter(|| {
            longitudinal_spherical_aberration(black_box(10.0), black_box(100.0), black_box(1.5))
        })
    });
    group.bench_function("chromatic_aberration", |b| {
        b.iter(|| chromatic_aberration(black_box(100.0), black_box(64.0)))
    });
    group.bench_function("petzval_sum", |b| {
        let elements = [(1.5, 100.0), (1.7, -150.0)];
        b.iter(|| petzval_sum(black_box(&elements)))
    });
    group.bench_function("separated_lenses_bfd", |b| {
        b.iter(|| separated_lenses_bfd(black_box(50.0), black_box(-80.0), black_box(30.0)))
    });
    group.bench_function("mtf_polychromatic_3wl", |b| {
        let wls = [486.1e-6, 587.6e-6, 656.3e-6];
        let wts = [1.0, 2.0, 1.0];
        b.iter(|| {
            mtf_polychromatic(
                black_box(&wls),
                black_box(&wts),
                black_box(5.6),
                black_box(50.0),
            )
        })
    });

    group.finish();
}

fn bench_pbr(c: &mut Criterion) {
    use prakash::pbr::*;
    let mut group = c.benchmark_group("pbr");

    group.bench_function("fresnel_schlick", |b| {
        b.iter(|| fresnel_schlick(black_box(0.04), black_box(0.8)))
    });
    group.bench_function("fresnel_schlick_rgb", |b| {
        b.iter(|| fresnel_schlick_rgb(black_box([0.56, 0.57, 0.58]), black_box(0.8)))
    });
    group.bench_function("distribution_ggx", |b| {
        b.iter(|| distribution_ggx(black_box(0.9), black_box(0.3)))
    });
    group.bench_function("distribution_beckmann", |b| {
        b.iter(|| distribution_beckmann(black_box(0.9), black_box(0.3)))
    });
    group.bench_function("geometry_schlick_ggx", |b| {
        b.iter(|| geometry_schlick_ggx(black_box(0.8), black_box(0.3)))
    });
    group.bench_function("geometry_smith", |b| {
        b.iter(|| geometry_smith(black_box(0.8), black_box(0.7), black_box(0.3)))
    });
    group.bench_function("cook_torrance", |b| {
        b.iter(|| {
            cook_torrance(
                black_box(0.9),
                black_box(0.85),
                black_box(0.8),
                black_box(0.7),
                black_box(0.3),
                black_box(0.04),
            )
        })
    });
    group.bench_function("lambert_diffuse", |b| {
        b.iter(|| lambert_diffuse(black_box(0.8)))
    });
    group.bench_function("lambert_diffuse_rgb", |b| {
        b.iter(|| lambert_diffuse_rgb(black_box([0.8, 0.6, 0.4])))
    });
    group.bench_function("ior_to_f0", |b| b.iter(|| ior_to_f0(black_box(1.5))));
    group.bench_function("ggx_aniso", |b| {
        b.iter(|| {
            distribution_ggx_aniso(
                black_box(0.9),
                black_box(0.3),
                black_box(0.2),
                black_box(0.3),
                black_box(0.7),
            )
        })
    });
    group.bench_function("geometry_smith_aniso", |b| {
        b.iter(|| {
            geometry_smith_aniso(
                black_box(0.8),
                black_box(0.7),
                black_box(0.3),
                black_box(0.2),
                black_box(0.3),
                black_box(0.2),
                black_box(0.3),
                black_box(0.7),
            )
        })
    });
    group.bench_function("sheen_charlie", |b| {
        b.iter(|| {
            sheen_charlie(
                black_box(0.9),
                black_box(0.8),
                black_box(0.7),
                black_box(0.5),
            )
        })
    });
    group.bench_function("sheen_ashikhmin", |b| {
        b.iter(|| sheen_ashikhmin(black_box(0.8), black_box(0.5)))
    });
    group.bench_function("clearcoat_brdf", |b| {
        b.iter(|| {
            clearcoat_brdf(
                black_box(0.9),
                black_box(0.8),
                black_box(0.7),
                black_box(0.85),
                black_box(0.05),
            )
        })
    });
    group.bench_function("clearcoat_blend", |b| {
        b.iter(|| {
            clearcoat_blend(
                black_box(0.5),
                black_box(0.3),
                black_box(0.8),
                black_box(0.85),
            )
        })
    });
    group.bench_function("sss_burley", |b| {
        b.iter(|| sss_profile_burley(black_box(0.5), black_box(1.0)))
    });
    group.bench_function("sss_gaussian", |b| {
        b.iter(|| sss_profile_gaussian(black_box(0.5), black_box(1.0)))
    });
    group.bench_function("subsurface_diffuse", |b| {
        b.iter(|| subsurface_diffuse(black_box(0.8), black_box(0.7), black_box(0.5)))
    });
    group.bench_function("iridescence_fresnel", |b| {
        b.iter(|| {
            iridescence_fresnel(
                black_box(1.0),
                black_box(1.3),
                black_box(1.5),
                black_box(300.0),
                black_box(550.0),
                black_box(0.8),
            )
        })
    });
    group.bench_function("iridescence_rgb", |b| {
        b.iter(|| {
            iridescence_rgb(
                black_box(1.0),
                black_box(1.3),
                black_box(1.5),
                black_box(300.0),
                black_box(0.8),
            )
        })
    });
    group.bench_function("henyey_greenstein", |b| {
        b.iter(|| henyey_greenstein(black_box(0.5), black_box(0.7)))
    });
    group.bench_function("phase_rayleigh", |b| {
        b.iter(|| phase_rayleigh(black_box(0.5)))
    });
    group.bench_function("volume_transmittance", |b| {
        b.iter(|| volume_transmittance(black_box(0.5), black_box(2.0)))
    });
    group.bench_function("inscattering", |b| {
        b.iter(|| {
            single_scatter_inscattering(
                black_box(0.5),
                black_box(0.8),
                black_box(0.5),
                black_box(1.0),
                black_box(0.5),
                black_box(1.0),
            )
        })
    });
    group.bench_function("sample_ggx", |b| {
        b.iter(|| sample_ggx(black_box(0.3), black_box(0.5), black_box(0.5)))
    });
    group.bench_function("sample_ggx_pdf", |b| {
        b.iter(|| sample_ggx_pdf(black_box(0.9), black_box(0.8), black_box(0.3)))
    });
    group.bench_function("sample_cosine", |b| {
        b.iter(|| sample_cosine_hemisphere(black_box(0.5), black_box(0.5)))
    });
    group.bench_function("split_sum", |b| {
        b.iter(|| split_sum_scale_bias(black_box(0.8), black_box(0.3)))
    });
    group.bench_function("integrate_brdf_64", |b| {
        b.iter(|| integrate_brdf_lut(black_box(0.5), black_box(0.3), black_box(64)))
    });

    group.finish();
}

fn bench_atmosphere(c: &mut Criterion) {
    use prakash::atmosphere::*;
    use std::f64::consts::PI;
    let mut group = c.benchmark_group("atmosphere");

    group.bench_function("rayleigh_cross_section", |b| {
        b.iter(|| rayleigh_cross_section(black_box(550e-9)))
    });
    group.bench_function("rayleigh_scattering_coeff", |b| {
        b.iter(|| rayleigh_scattering_coefficient(black_box(550e-9)))
    });
    group.bench_function("rayleigh_at_altitude", |b| {
        b.iter(|| rayleigh_scattering_at_altitude(black_box(550e-9), black_box(5000.0)))
    });
    group.bench_function("phase_rayleigh", |b| {
        b.iter(|| phase_rayleigh(black_box(0.5)))
    });
    group.bench_function("mie_phase_cornette_shanks", |b| {
        b.iter(|| mie_phase_cornette_shanks(black_box(0.5), black_box(0.76)))
    });
    group.bench_function("air_mass", |b| b.iter(|| air_mass(black_box(1.2))));
    group.bench_function("air_mass_horizon", |b| {
        b.iter(|| air_mass(black_box(PI / 2.0)))
    });
    group.bench_function("optical_depth_rayleigh", |b| {
        b.iter(|| optical_depth_rayleigh(black_box(550e-9)))
    });
    group.bench_function("atmospheric_transmittance", |b| {
        b.iter(|| atmospheric_transmittance(black_box(550e-9), black_box(1.0)))
    });
    group.bench_function("sky_color_rgb", |b| {
        b.iter(|| sky_color_rgb(black_box(0.5), black_box(PI / 2.0)))
    });
    group.bench_function("sunlight_color", |b| {
        b.iter(|| sunlight_color(black_box(85.0_f64.to_radians())))
    });
    group.bench_function("sunset_gradient", |b| {
        b.iter(|| {
            sunset_gradient(
                black_box(85.0_f64.to_radians()),
                black_box(10.0_f64.to_radians()),
                black_box(0.3),
            )
        })
    });
    group.bench_function("scattering_angle", |b| {
        b.iter(|| scattering_angle(black_box(0.5), black_box(1.0), black_box(0.8)))
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_ray,
    bench_spectral,
    bench_wave,
    bench_lens,
    bench_pbr,
    bench_atmosphere
);
criterion_main!(benches);
