use criterion::{Criterion, black_box, criterion_group, criterion_main};
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
                black_box(1e-3),
                black_box(500e-9),
                black_box(0.01),
                black_box(1.0),
            )
        })
    });
    group.bench_function("double_slit", |b| {
        b.iter(|| {
            double_slit_intensity(
                black_box(0.1e-3),
                black_box(0.5e-3),
                black_box(500e-9),
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
        b.iter(|| grating_maxima(black_box(1e-6), black_box(500e-9), black_box(3)))
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

    group.finish();
}

criterion_group!(
    benches,
    bench_ray,
    bench_spectral,
    bench_wave,
    bench_lens,
    bench_pbr
);
criterion_main!(benches);
