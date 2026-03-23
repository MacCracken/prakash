//! Example: basic optics demonstrations.

fn main() {
    // Snell's law
    let angle_i = prakash::ray::deg_to_rad(30.0);
    let angle_t = prakash::ray::snell(
        prakash::ray::Medium::AIR.n,
        prakash::ray::Medium::GLASS.n,
        angle_i,
    ).unwrap();
    println!("Light entering glass at 30°: refracted to {:.1}°", angle_t.to_degrees());

    // Fresnel reflectance
    let r = prakash::ray::fresnel_normal(1.0, 1.52);
    println!("Glass reflects {:.1}% at normal incidence", r * 100.0);

    // Color of light
    let rgb = prakash::spectral::wavelength_to_rgb(550.0).unwrap();
    println!("550nm light: RGB({:.2}, {:.2}, {:.2})", rgb.r, rgb.g, rgb.b);

    // Sun's peak wavelength
    let peak = prakash::spectral::wien_peak(5778.0) * 1e9;
    println!("Sun peaks at {:.0}nm", peak);

    // Camera lens
    let di = prakash::lens::thin_lens_image_distance(50.0, 2000.0).unwrap();
    let mag = prakash::lens::magnification(2000.0, di);
    println!("50mm lens, subject at 2m: image at {:.1}mm, mag={:.3}", di, mag);
}
