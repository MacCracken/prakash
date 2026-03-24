# Benchmarks

Three-point tracking: **baseline** (first run) / **previous** / **latest**

| Point | Date | Commit |
|-------|------|--------|
| Baseline | 2026-03-24T18:16:28Z | `e6cfd50` |
| Latest | 2026-03-24T18:48:48Z | `e6cfd50` |

## ray

| Benchmark | Baseline | Latest |
|-----------|----------|--------|
| `snell_air_glass` | 14.50 ns | 14.90 ns |
| `snell_tir_check` | 16.85 ns | 19.17 ns |
| `fresnel_normal` | 1.23 ns | 1.32 ns |
| `fresnel_s` | 2.06 ns | 2.42 ns |
| `fresnel_p` | 2.04 ns | 2.19 ns |
| `fresnel_unpolarized` | 34.56 ns | 35.55 ns |
| `critical_angle` | 7.23 ns | 7.39 ns |
| `reflect_2d` | 1.80 ns | 1.92 ns |
| `reflect_3d` | 2.54 ns | 2.61 ns |
| `beer_lambert` | 4.78 ns | 5.09 ns |
| `brewster_angle` | 8.50 ns | 8.80 ns |
| `refract_3d` | 9.02 ns | 9.52 ns |
| `snell_3d` | 17.99 ns | 18.23 ns |
| `sellmeier_n_at` | 6.01 ns | 6.15 ns |
| `cauchy_n_at` | 1.16 ns | 1.20 ns |
| `abbe_number` | 13.63 ns | 13.86 ns |
| `prism_deviation` | 17.02 ns | 17.52 ns |
| `trace_plane_surface` | 13.92 ns | 14.29 ns |
| `trace_sphere_surface` | 32.13 ns | 33.60 ns |
| `trace_sequential_2` | 141.54 ns | 144.82 ns |
| `prism_dispersion` | 22.18 ns | 23.00 ns |
| `prism_angular_spread` | 43.13 ns | 43.68 ns |

## spectral

| Benchmark | Baseline | Latest |
|-----------|----------|--------|
| `wavelength_to_rgb` | 3.22 ns | 3.35 ns |
| `wavelength_to_rgb_edge` | 2.81 ns | 2.95 ns |
| `planck_radiance` | 6.87 ns | 7.08 ns |
| `wien_peak` | 1.16 ns | 1.23 ns |
| `color_temp_to_rgb` | 12.94 ns | 13.96 ns |
| `color_temp_warm` | 13.01 ns | 16.91 ns |
| `photon_energy_ev` | 2.33 ns | 2.63 ns |
| `wavelength_to_frequency` | 1.17 ns | 1.29 ns |
| `frequency_to_wavelength` | 1.18 ns | 1.26 ns |
| `photon_energy` | 1.18 ns | 1.28 ns |
| `rgb_luminance` | 1.28 ns | 1.55 ns |
| `rgb_to_u8` | 3.52 ns | 4.16 ns |
| `cie_cmf_at` | 3.93 ns | 4.64 ns |
| `cie_cmf_at_interp` | 3.90 ns | 4.24 ns |
| `xyz_to_xyy` | 2.30 ns | 2.85 ns |
| `xyz_to_srgb` | 54.97 ns | 67.11 ns |
| `xyz_to_linear_srgb` | 2.21 ns | 2.75 ns |
| `linear_srgb_to_xyz` | 1.97 ns | 1.93 ns |
| `srgb_gamma` | 12.69 ns | 14.64 ns |
| `srgb_gamma_inv` | 13.29 ns | 15.71 ns |
| `cct_from_xy` | 1.97 ns | 2.36 ns |
| `spd_blackbody` | 617.18 ns | 675.67 ns |
| `spd_to_xyz` | 65.55 ns | 75.68 ns |
| `spd_to_srgb` | 116.85 ns | 131.53 ns |
| `cri` | 894.16 ns | 1.02 µs |

## wave

| Benchmark | Baseline | Latest |
|-----------|----------|--------|
| `interference` | 8.01 ns | 9.00 ns |
| `single_slit` | 18.25 ns | 20.54 ns |
| `double_slit` | 28.13 ns | 31.13 ns |
| `malus_law` | 7.92 ns | 8.49 ns |
| `thin_film` | 11.94 ns | 13.19 ns |
| `path_to_phase` | 1.18 ns | 1.24 ns |
| `grating_maxima_3` | 32.42 ns | 32.29 ns |
| `polarization_intensity` | 1.26 ns | 1.31 ns |
| `polarization_through` | 11.79 ns | 12.20 ns |
| `is_constructive` | 2.85 ns | 3.20 ns |
| `is_destructive` | 2.85 ns | 3.18 ns |
| `coherence_length` | 1.20 ns | 1.30 ns |
| `coherence_time` | 2.33 ns | 2.58 ns |
| `bessel_j1` | 11.81 ns | 11.59 ns |
| `bessel_j1_large` | 14.82 ns | 19.78 ns |
| `airy_pattern` | 13.84 ns | 16.18 ns |
| `fp_transmittance` | 19.11 ns | 22.01 ns |
| `fp_finesse` | 3.44 ns | 3.81 ns |
| `fp_resolving_power` | 4.05 ns | 5.04 ns |
| `stokes_dop` | 3.65 ns | 4.42 ns |
| `mueller_apply` | 6.66 ns | 8.08 ns |
| `mueller_multiply` | 13.19 ns | 15.14 ns |
| `mueller_polarizer` | 14.04 ns | 17.86 ns |
| `mueller_retarder` | 12.08 ns | 15.61 ns |
| `mueller_chain_3` | 9.90 ns | 13.80 ns |
| `birefringent_retardation` | 1.21 ns | 1.73 ns |
| `fraunhofer_rect` | 39.66 ns | 49.00 ns |
| `fresnel_c` | 16.61 ns | 22.81 ns |
| `fresnel_s` | 16.86 ns | 20.24 ns |
| `fresnel_edge` | 19.18 ns | 24.18 ns |
| `huygens_fresnel_50` | 1.24 µs | 1.92 µs |
| `coating_reflectance` | 16.41 ns | 19.94 ns |
| `multilayer_2` | 35.64 ns | 45.38 ns |

## lens

| Benchmark | Baseline | Latest |
|-----------|----------|--------|
| `thin_lens` | 1.67 ns | 2.15 ns |
| `magnification` | 1.26 ns | 1.47 ns |
| `lensmaker` | 2.91 ns | 3.35 ns |
| `depth_of_field` | 3.72 ns | 4.43 ns |
| `optical_power` | 2.02 ns | 2.05 ns |
| `combined_focal` | 1.83 ns | 1.97 ns |
| `mirror_focal` | 521.60 ps | 804.10 ps |
| `thick_lens` | 3.80 ns | 5.03 ns |
| `cardinal_points` | 6.03 ns | 8.48 ns |
| `f_number` | 1.18 ns | 1.37 ns |
| `numerical_aperture` | 8.10 ns | 10.07 ns |
| `airy_disk_radius` | 780.90 ps | 988.10 ps |
| `field_of_view` | 7.84 ns | 9.93 ns |
| `mtf_diffraction` | 8.13 ns | 10.26 ns |
| `seidel_coefficients` | 14.33 ns | 15.39 ns |
| `separated_lenses` | 4.17 ns | 4.20 ns |
| `aperture_from_f_number` | 1.33 ns | 1.39 ns |
| `na_from_f_number` | 1.26 ns | 1.28 ns |
| `diffraction_limit` | 1.20 ns | 1.43 ns |
| `field_of_view_diagonal` | 15.97 ns | 18.14 ns |
| `mtf_cutoff_frequency` | 1.36 ns | 1.47 ns |
| `shape_factor` | 1.61 ns | 1.24 ns |
| `conjugate_factor` | 1.71 ns | 1.31 ns |
| `lsa` | 5.23 ns | 4.68 ns |
| `chromatic_aberration` | 1.27 ns | 1.37 ns |
| `petzval_sum` | 2.44 ns | 2.70 ns |
| `separated_lenses_bfd` | 4.87 ns | 6.35 ns |

## pbr

| Benchmark | Baseline | Latest |
|-----------|----------|--------|
| `fresnel_schlick` | 1.28 ns | 1.36 ns |
| `fresnel_schlick_rgb` | 2.32 ns | 2.33 ns |
| `distribution_ggx` | 3.79 ns | 2.49 ns |
| `distribution_beckmann` | 15.12 ns | 14.14 ns |
| `geometry_schlick_ggx` | 1.71 ns | 3.11 ns |
| `geometry_smith` | 2.91 ns | 3.83 ns |
| `cook_torrance` | 7.63 ns | 8.42 ns |
| `lambert_diffuse` | 1.20 ns | 1.45 ns |
| `lambert_diffuse_rgb` | 2.43 ns | 3.12 ns |
| `ior_to_f0` | 1.24 ns | 1.68 ns |
| `ggx_aniso` | 4.05 ns | 8.43 ns |
| `geometry_smith_aniso` | 6.99 ns | 13.80 ns |
| `sheen_charlie` | 19.27 ns | 30.86 ns |
| `sheen_ashikhmin` | 2.33 ns | 1.77 ns |
| `clearcoat_brdf` | 8.69 ns | 5.80 ns |
| `clearcoat_blend` | 2.71 ns | 2.93 ns |
| `sss_burley` | 15.40 ns | 15.12 ns |
| `sss_gaussian` | 7.21 ns | 7.44 ns |
| `subsurface_diffuse` | 5.68 ns | 6.38 ns |
| `iridescence_fresnel` | 38.98 ns | 31.05 ns |
| `iridescence_rgb` | 62.04 ns | 64.11 ns |
| `henyey_greenstein` | 4.38 ns | 5.07 ns |
| `phase_rayleigh` | 1.20 ns | 1.47 ns |
| `volume_transmittance` | 4.63 ns | 5.92 ns |
| `inscattering` | 8.04 ns | 9.20 ns |
| `sample_ggx` | 22.97 ns | 27.88 ns |
| `sample_ggx_pdf` | 3.71 ns | 4.32 ns |
| `sample_cosine` | 14.45 ns | 18.44 ns |
| `split_sum` | 3.40 ns | 4.09 ns |
| `integrate_brdf_64` | 1.63 µs | 2.19 µs |

---

Generated by `./scripts/bench-history.sh`. Full history in `bench-history.csv`.
