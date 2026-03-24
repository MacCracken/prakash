//! Camera lens simulator — multi-element trace + spot diagram.
//!
//! Defines a simple camera lens as an optical prescription, traces rays
//! through it, computes spot diagrams, and analyzes image quality.

fn main() {
    use prakash::ray::*;

    println!("╔══════════════════════════════════════════╗");
    println!("║   Camera Lens Simulator (Prakash)        ║");
    println!("╚══════════════════════════════════════════╝\n");

    // ── Define a cemented doublet (50mm f/4) ────────────────────────────
    let lens = Prescription::new("50mm f/4 doublet")
        .add_surface(28.0, 4.0, 1.517, 6.25) // BK7 crown, R1
        .add_surface(-22.0, 2.0, 1.648, 6.25) // SF2 flint, R2 (cemented)
        .add_surface(-80.0, 0.0, 1.0, 6.25); // R3 → air

    println!("Lens: {}", lens.name);
    println!("Surfaces: {}", lens.len());
    println!();

    // ── Paraxial analysis ────────────────────────────────────────────────
    let props = find_system_properties(&lens).unwrap();
    println!("System Properties (paraxial):");
    println!("  Effective focal length: {:.2} mm", props.focal_length);
    println!("  Back focal distance:    {:.2} mm", props.bfd);
    println!("  Front focal distance:   {:.2} mm", props.ffd);
    println!("  Optical power:          {:.4} mm⁻¹", props.power);
    println!();

    // ── Paraxial ray trace ───────────────────────────────────────────────
    let paraxial_surfaces = lens.to_paraxial_surfaces();
    let marginal = ParaxialRay::marginal(6.25);
    let states = paraxial_trace(&marginal, &paraxial_surfaces);
    println!("Paraxial marginal ray (y=6.25mm):");
    for (i, s) in states.iter().enumerate() {
        println!("  Surface {}: y={:.4} mm, nu={:.6}", i, s.y, s.nu);
    }
    println!();

    // ── Full ray trace ───────────────────────────────────────────────────
    let trace_surfaces = lens.to_trace_surfaces();
    let image_z = props.bfd + 6.0; // approximate image plane

    println!("Full ray trace (on-axis, height=3mm):");
    let ray = TraceRay {
        position: [0.0, 3.0, -100.0],
        direction: [0.0, 0.0, 1.0],
        n: 1.0,
    };
    match trace_sequential(&ray, &trace_surfaces) {
        Ok(hits) => {
            for (i, hit) in hits.iter().enumerate() {
                println!(
                    "  Surface {}: hit=({:.3}, {:.3}, {:.3}), R={:.4}",
                    i, hit.hit_point[0], hit.hit_point[1], hit.hit_point[2], hit.reflectance
                );
            }
        }
        Err(e) => println!("  Ray vignetted: {e}"),
    }
    println!();

    // ── Spot diagram ─────────────────────────────────────────────────────
    println!("Spot diagrams:");
    for (label, field_deg) in [
        ("On-axis", 0.0_f64),
        ("5° off-axis", 5.0),
        ("10° off-axis", 10.0),
    ] {
        let field_rad = field_deg.to_radians();
        let spots = spot_diagram(
            &trace_surfaces,
            5.0, // aperture radius
            3,   // rings
            8,   // arms per ring
            field_rad,
            -100.0, // start z
            image_z,
        );
        let rms = spot_rms_radius(&spots);
        println!(
            "  {:<15}: {} rays landed, RMS radius = {:.3} mm",
            label,
            spots.len(),
            rms
        );
    }
    println!();

    // ── OPD analysis ─────────────────────────────────────────────────────
    println!("OPD fan (meridional, on-axis):");
    let fan = ray_fan_meridional(5.0, 11, 0.0, -100.0);
    let opd_data = opd_fan(&fan, &trace_surfaces, image_z);
    for (h, opd) in &opd_data {
        let bar_len = (opd.abs() * 1000.0).min(40.0) as usize;
        let bar = if *opd >= 0.0 {
            format!("{:>20}{}", "", "█".repeat(bar_len))
        } else {
            let offset = 20 - bar_len.min(20);
            format!("{:>offset$}{}", "", "█".repeat(bar_len))
        };
        println!("  h={:+.2}: OPD={:+.6} mm {}", h, opd, bar);
    }

    // ── Recursive trace ──────────────────────────────────────────────────
    println!();
    println!("Recursive trace (showing ghost reflections):");
    let ray = TraceRay {
        position: [0.0, 0.0, -100.0],
        direction: [0.0, 0.0, 1.0],
        n: 1.0,
    };
    let config = TraceConfig {
        max_depth: 4,
        min_energy: 0.001,
    };
    let tree = trace_recursive(&ray, &trace_surfaces, &config);
    println!(
        "  {} segments, {} interactions",
        tree.segments.len(),
        tree.interactions
    );
    for seg in &tree.segments {
        let event_str = match seg.event {
            TraceEvent::Refraction {
                reflectance,
                surface_idx,
            } => {
                format!("refract(S{}, R={:.4})", surface_idx, reflectance)
            }
            TraceEvent::Reflection { surface_idx } => format!("reflect(S{})", surface_idx),
            TraceEvent::Escaped => "escaped".to_string(),
            TraceEvent::Terminated => "terminated".to_string(),
            _ => "unknown".to_string(),
        };
        println!(
            "  depth={}, energy={:.4}, OPL={:.2}, {}",
            seg.depth, seg.energy, seg.optical_path, event_str
        );
    }
}
