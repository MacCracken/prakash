//! Physically accurate rainbow simulation.
//!
//! Traces sunlight through a spherical raindrop to show how dispersion
//! separates white light into the rainbow spectrum. Uses Sellmeier water
//! coefficients for wavelength-dependent refraction and computes the
//! deviation angle for each color.

fn main() {
    use prakash::ray::SellmeierCoefficients;

    println!("╔══════════════════════════════════════════╗");
    println!("║     Rainbow Simulation (Prakash)         ║");
    println!("╚══════════════════════════════════════════╝\n");

    // Sunlight enters a raindrop, refracts, reflects off the back, refracts out.
    // The minimum deviation angle determines the rainbow angle.
    //
    // For a sphere: δ_min = 2(θ_i - θ_r) + (π - 2θ_r)
    // where θ_i is incidence, θ_r is refraction angle.
    // The rainbow angle is 180° - δ_min.

    // Sample visible wavelengths
    let wavelengths_nm = [
        (380.0, "Violet"),
        (420.0, "Indigo"),
        (450.0, "Blue"),
        (490.0, "Cyan"),
        (530.0, "Green"),
        (570.0, "Yellow"),
        (600.0, "Orange"),
        (650.0, "Red"),
        (700.0, "Deep Red"),
    ];

    println!("Primary Rainbow (single internal reflection):");
    println!("┌──────────┬──────────┬──────────┬────────────┐");
    println!("│ Color    │ λ (nm)   │ n(water) │ Angle (°)  │");
    println!("├──────────┼──────────┼──────────┼────────────┤");

    for (wl_nm, name) in &wavelengths_nm {
        let wl_um = wl_nm / 1000.0;
        let n = SellmeierCoefficients::WATER.n_at(wl_um);

        // Find minimum deviation by scanning incidence angles
        let mut min_dev = f64::MAX;
        for i in 0..9000 {
            let theta_i = (i as f64) * 0.01_f64.to_radians();
            let sin_r = theta_i.sin() / n;
            if sin_r.abs() > 1.0 {
                continue;
            }
            let theta_r = sin_r.asin();
            // Deviation for primary bow: δ = 2(θi - θr) + (π - 2θr)
            let dev = 2.0 * (theta_i - theta_r) + (std::f64::consts::PI - 2.0 * theta_r);
            if dev < min_dev {
                min_dev = dev;
            }
        }
        let rainbow_angle = 180.0 - min_dev.to_degrees();

        println!(
            "│ {:<8} │ {:>6.0}   │ {:.5}  │ {:>8.2}°  │",
            name, wl_nm, n, rainbow_angle
        );
    }

    println!("└──────────┴──────────┴──────────┴────────────┘");
    println!("\nRed light bends least → appears on the outside of the bow (~42.4°)");
    println!("Violet bends most → appears on the inside (~40.5°)");

    // Show the RGB colors along the rainbow
    println!("\nRainbow spectrum:");
    for wl in (400..=700).step_by(10) {
        if let Ok(rgb) = prakash::spectral::wavelength_to_rgb(wl as f64) {
            let u8 = rgb.to_u8();
            print!("\x1b[48;2;{};{};{}m  \x1b[0m", u8[0], u8[1], u8[2]);
        }
    }
    println!("  (380nm → 700nm)");

    // Dispersion statistics
    println!("\nDispersion analysis:");
    let n_red = SellmeierCoefficients::WATER.n_at(0.656);
    let n_blue = SellmeierCoefficients::WATER.n_at(0.486);
    println!("  n(red,  656nm) = {:.6}", n_red);
    println!("  n(blue, 486nm) = {:.6}", n_blue);
    println!("  Δn = {:.6} (dispersion)", n_blue - n_red);
}
