use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn bench_ray(c: &mut Criterion) {
    use prakash::ray::*;
    let mut group = c.benchmark_group("ray");

    group.bench_function("snell_air_glass", |b| {
        b.iter(|| snell(black_box(1.0), black_box(1.52), black_box(0.5236)))
    });
    group.bench_function("fresnel_normal", |b| {
        b.iter(|| fresnel_normal(black_box(1.0), black_box(1.52)))
    });
    group.bench_function("fresnel_unpolarized", |b| {
        b.iter(|| fresnel_unpolarized(black_box(1.0), black_box(1.52), black_box(0.5)))
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

    group.finish();
}

fn bench_spectral(c: &mut Criterion) {
    use prakash::spectral::*;
    let mut group = c.benchmark_group("spectral");

    group.bench_function("wavelength_to_rgb", |b| {
        b.iter(|| wavelength_to_rgb(black_box(550.0)))
    });
    group.bench_function("planck_radiance", |b| {
        b.iter(|| planck_radiance(black_box(500e-9), black_box(5778.0)))
    });
    group.bench_function("wien_peak", |b| {
        b.iter(|| wien_peak(black_box(5778.0)))
    });
    group.bench_function("color_temp_to_rgb", |b| {
        b.iter(|| color_temperature_to_rgb(black_box(6500.0)))
    });
    group.bench_function("photon_energy_ev", |b| {
        b.iter(|| photon_energy_ev(black_box(550.0)))
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
        b.iter(|| single_slit_intensity(black_box(1e-3), black_box(500e-9), black_box(0.01), black_box(1.0)))
    });
    group.bench_function("malus_law", |b| {
        b.iter(|| malus_law(black_box(1.0), black_box(0.785)))
    });
    group.bench_function("thin_film", |b| {
        b.iter(|| thin_film_reflectance(black_box(550.0), black_box(100.0), black_box(1.5)))
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
        b.iter(|| depth_of_field(black_box(50.0), black_box(2.8), black_box(0.03), black_box(2000.0)))
    });

    group.finish();
}

fn bench_pbr(c: &mut Criterion) {
    use prakash::pbr::*;
    let mut group = c.benchmark_group("pbr");

    group.bench_function("fresnel_schlick", |b| {
        b.iter(|| fresnel_schlick(black_box(0.04), black_box(0.8)))
    });
    group.bench_function("distribution_ggx", |b| {
        b.iter(|| distribution_ggx(black_box(0.9), black_box(0.3)))
    });
    group.bench_function("geometry_smith", |b| {
        b.iter(|| geometry_smith(black_box(0.8), black_box(0.7), black_box(0.3)))
    });
    group.bench_function("cook_torrance", |b| {
        b.iter(|| cook_torrance(black_box(0.9), black_box(0.8), black_box(0.7), black_box(0.3), black_box(0.04)))
    });
    group.bench_function("lambert_diffuse", |b| {
        b.iter(|| lambert_diffuse(black_box(0.8)))
    });

    group.finish();
}

criterion_group!(benches, bench_ray, bench_spectral, bench_wave, bench_lens, bench_pbr);
criterion_main!(benches);
