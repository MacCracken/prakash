//! PBR material preview — Cook-Torrance + clearcoat + iridescence.
//!
//! Evaluates physically-based rendering BRDFs at various angles to show
//! how different material properties affect light reflection.

fn main() {
    use prakash::pbr::*;

    println!("╔══════════════════════════════════════════╗");
    println!("║     PBR Material Preview (Prakash)       ║");
    println!("╚══════════════════════════════════════════╝\n");

    // ── Fresnel-Schlick at various angles ────────────────────────────────
    println!("Fresnel-Schlick reflectance (F0=0.04, dielectric):");
    println!("┌────────────┬────────────┐");
    println!("│ Angle (°)  │ Reflectance│");
    println!("├────────────┼────────────┤");
    for deg in [0, 15, 30, 45, 60, 75, 85, 89] {
        let cos_theta = (deg as f64).to_radians().cos();
        let f = fresnel_schlick(0.04, cos_theta);
        println!("│ {:>8}°  │ {:>9.4}  │", deg, f);
    }
    println!("└────────────┴────────────┘\n");

    // ── Cook-Torrance BRDF for different roughness ──────────────────────
    println!("Cook-Torrance specular (n·h=0.95, h·v=0.9, n·v=0.8, n·l=0.7):");
    println!("┌───────────┬──────────────┬──────────────┐");
    println!("│ Roughness │  Dielectric  │    Metal     │");
    println!("├───────────┼──────────────┼──────────────┤");
    for r in [0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0] {
        let d = cook_torrance(0.95, 0.9, 0.8, 0.7, r, 0.04);
        let m = cook_torrance(0.95, 0.9, 0.8, 0.7, r, 0.95);
        println!("│ {:>7.2}   │ {:>10.4}   │ {:>10.4}   │", r, d, m);
    }
    println!("└───────────┴──────────────┴──────────────┘\n");

    // ── Clearcoat effect ─────────────────────────────────────────────────
    println!("Clearcoat blend (base=0.5, coat_roughness=0.05):");
    let base = 0.5;
    let coat = clearcoat_brdf(0.95, 0.8, 0.7, 0.9, 0.05);
    println!("  Coat BRDF:     {:.4}", coat);
    for intensity in [0.0, 0.25, 0.5, 0.75, 1.0] {
        let blended = clearcoat_blend(base, coat, intensity, 0.9);
        println!("  Clearcoat={:.2}: blended={:.4}", intensity, blended);
    }
    println!();

    // ── Iridescence ──────────────────────────────────────────────────────
    println!("Iridescence (thin film, n_film=1.3, n_base=1.5):");
    println!("┌────────────┬────────────┬────────────┬────────────┐");
    println!("│ Thickness  │     R      │     G      │     B      │");
    println!("├────────────┼────────────┼────────────┼────────────┤");
    for t in [100.0, 200.0, 300.0, 400.0, 500.0, 600.0] {
        let rgb = iridescence_rgb(1.0, 1.3, 1.5, t, 0.8);
        println!(
            "│ {:>6.0} nm  │ {:>8.4}   │ {:>8.4}   │ {:>8.4}   │",
            t, rgb[0], rgb[1], rgb[2]
        );
    }
    println!("└────────────┴────────────┴────────────┴────────────┘\n");

    // ── Subsurface scattering profiles ───────────────────────────────────
    println!("SSS profiles (d=1.0):");
    println!("  Distance │ Burley     │ Gaussian(σ=1)");
    println!("  ─────────┼────────────┼─────────────");
    for dist in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0] {
        let b = sss_profile_burley(dist, 1.0);
        let g = sss_profile_gaussian(dist, 1.0);
        println!("  {:>7.2}  │ {:>8.5}   │ {:>8.5}", dist, b, g);
    }
    println!();

    // ── Sheen comparison ─────────────────────────────────────────────────
    println!("Sheen (fabric/velvet):");
    println!("  cos(θ) │ Charlie(r=0.5) │ Ashikhmin(i=0.8)");
    println!("  ───────┼────────────────┼─────────────────");
    for ct_i in [0, 2, 4, 6, 8, 10] {
        let ct = ct_i as f64 / 10.0;
        let charlie = sheen_charlie(ct.max(0.01), 0.7, 0.8, 0.5);
        let ashikhmin = sheen_ashikhmin(ct, 0.8);
        println!("  {:>5.1}  │ {:>12.6}   │ {:>12.6}", ct, charlie, ashikhmin);
    }
    println!();

    // ── Volumetric scattering ────────────────────────────────────────────
    println!("Phase functions at various scattering angles:");
    println!("  cos(θ) │ HG(g=0.76) │ Rayleigh   │ Isotropic");
    println!("  ───────┼────────────┼────────────┼──────────");
    for ct_i in [-10, -5, 0, 5, 8, 10] {
        let ct = ct_i as f64 / 10.0;
        let hg = henyey_greenstein(ct, 0.76);
        let ray = phase_rayleigh(ct);
        let iso = phase_isotropic();
        println!(
            "  {:>+5.1}  │ {:>8.5}   │ {:>8.5}   │ {:>8.5}",
            ct, hg, ray, iso
        );
    }

    // ── Importance sampling validation ───────────────────────────────────
    println!();
    println!("GGX importance sampling (roughness=0.3):");
    let n_samples = 8;
    for i in 0..n_samples {
        let xi1 = (i as f64 + 0.5) / n_samples as f64;
        let h = sample_ggx(0.3, xi1, 0.0);
        let pdf = sample_ggx_pdf(h[2], h[2], 0.3);
        println!(
            "  ξ₁={:.3}: h=({:+.4}, {:+.4}, {:.4}), pdf={:.4}",
            xi1, h[0], h[1], h[2], pdf
        );
    }
}
