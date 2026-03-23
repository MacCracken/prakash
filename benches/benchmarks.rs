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
    group.bench_function("rgb_luminance", |b| {
        let rgb = Rgb::new(0.5, 0.7, 0.3);
        b.iter(|| black_box(rgb).luminance())
    });
    group.bench_function("rgb_to_u8", |b| {
        let rgb = Rgb::new(0.5, 0.7, 0.3);
        b.iter(|| black_box(rgb).to_u8())
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
