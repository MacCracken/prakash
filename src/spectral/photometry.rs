//! Photometry — radiometric-to-photometric conversion via V(λ).
//!
//! Bridges spectral radiance/power to human-perceived quantities (lumens, lux, candela).

use super::Spd;

/// Maximum luminous efficacy at 555nm: 683 lm/W.
pub const K_M: f64 = 683.0;

/// CIE 1924 photopic luminous efficiency function V(λ).
///
/// Tabulated at 5nm intervals, 380–780nm (81 values).
/// Represents the relative sensitivity of the human eye under daylight conditions.
pub const V_LAMBDA_PHOTOPIC: [f64; 81] = [
    // 380-400
    0.000039, 0.000064, 0.000120, 0.000217, 0.000396, // 405-425
    0.000640, 0.001210, 0.002180, 0.004000, 0.007300, // 430-450
    0.011600, 0.016840, 0.023000, 0.029800, 0.038000, // 455-475
    0.048000, 0.060000, 0.073900, 0.090980, 0.112600, // 480-500
    0.139020, 0.169300, 0.208020, 0.258600, 0.323000, // 505-525
    0.407300, 0.503000, 0.608200, 0.710000, 0.793200, // 530-550
    0.862000, 0.914850, 0.954000, 0.980300, 0.994950, // 555-575
    1.000000, 0.995000, 0.978600, 0.952000, 0.915400, // 580-600
    0.870000, 0.816300, 0.757000, 0.694900, 0.631000, // 605-625
    0.566800, 0.503000, 0.441200, 0.381000, 0.321000, // 630-650
    0.265000, 0.217000, 0.175000, 0.138200, 0.107000, // 655-675
    0.081600, 0.061000, 0.044580, 0.032000, 0.023200, // 680-700
    0.017000, 0.011920, 0.008210, 0.005723, 0.004102, // 705-725
    0.002929, 0.002091, 0.001484, 0.001047, 0.000740, // 730-750
    0.000520, 0.000361, 0.000249, 0.000172, 0.000120, // 755-780
    0.000085, 0.000060, 0.000042, 0.000030, 0.000021, 0.000015,
];

/// CIE 1951 scotopic luminous efficiency function V'(λ).
///
/// Tabulated at 5nm intervals, 380–780nm (81 values).
/// Represents rod-mediated vision (night/dark-adapted).
/// Peaks at 507nm (vs 555nm for photopic).
pub const V_LAMBDA_SCOTOPIC: [f64; 81] = [
    // 380-400
    0.000589, 0.001108, 0.002209, 0.004530, 0.009290, // 405-425
    0.018480, 0.034840, 0.060400, 0.096600, 0.143900, // 430-450
    0.199800, 0.262500, 0.328100, 0.405300, 0.484900, // 455-475
    0.567600, 0.649400, 0.726000, 0.793200, 0.851400, // 480-500
    0.899200, 0.937600, 0.963100, 0.978300, 0.982000, // 505-525
    0.980600, 0.970600, 0.948600, 0.913600, 0.869500, // 530-550
    0.811000, 0.733000, 0.650000, 0.564000, 0.481000, // 555-575
    0.402000, 0.328800, 0.264500, 0.207600, 0.160200, // 580-600
    0.121200, 0.089900, 0.065500, 0.046900, 0.033100, // 605-625
    0.022700, 0.015700, 0.010400, 0.006900, 0.004530, // 630-650
    0.002900, 0.001840, 0.001140, 0.000690, 0.000420, // 655-675
    0.000250, 0.000149, 0.000088, 0.000052, 0.000030, // 680-700
    0.000018, 0.000010, 0.000006, 0.000003, 0.000002, // 705-725
    0.000001, 0.000001, 0.000000, 0.000000, 0.000000, // 730-750
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, // 755-780
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
];

/// Luminous flux from a spectral power distribution (photopic).
///
/// Φ_v = K_m · ∫ Φ(λ) · V(λ) dλ
///
/// where K_m = 683 lm/W. Uses 5nm integration matching the V(λ) table.
/// Returns luminous flux in lumens (relative to the SPD's power units).
#[must_use]
pub fn luminous_flux(spd: &Spd) -> f64 {
    let mut sum = 0.0;
    for (i, &v) in V_LAMBDA_PHOTOPIC.iter().enumerate() {
        let wl = 380.0 + i as f64 * 5.0;
        sum += spd.at(wl) * v;
    }
    K_M * sum * 5.0 // 5nm step
}

/// Luminous efficacy of a source.
///
/// η = Φ_v / Φ_e (lm/W)
///
/// Computes the ratio of luminous flux to total radiant flux.
/// A source emitting only at 555nm has η = 683 lm/W (maximum).
///
/// `spd` is the spectral power distribution. Returns efficacy in lm/W.
#[must_use]
pub fn luminous_efficacy(spd: &Spd) -> f64 {
    let mut radiant_sum = 0.0;
    let mut luminous_sum = 0.0;
    for (i, &v) in V_LAMBDA_PHOTOPIC.iter().enumerate() {
        let wl = 380.0 + i as f64 * 5.0;
        let power = spd.at(wl);
        radiant_sum += power;
        luminous_sum += power * v;
    }
    if radiant_sum < 1e-15 {
        return 0.0;
    }
    K_M * luminous_sum / radiant_sum
}

/// Luminous flux from a spectral power distribution (scotopic / night vision).
///
/// Uses V'(λ) with K'_m = 1700 lm/W.
#[must_use]
pub fn luminous_flux_scotopic(spd: &Spd) -> f64 {
    const K_M_SCOTOPIC: f64 = 1700.0;
    let mut sum = 0.0;
    for (i, &v) in V_LAMBDA_SCOTOPIC.iter().enumerate() {
        let wl = 380.0 + i as f64 * 5.0;
        sum += spd.at(wl) * v;
    }
    K_M_SCOTOPIC * sum * 5.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spectral::cie::{Spd, illuminant_d65};

    #[test]
    fn test_v_lambda_photopic_peak() {
        // V(λ) peaks at 555nm (index 35)
        assert!((V_LAMBDA_PHOTOPIC[35] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_v_lambda_photopic_range() {
        for &v in &V_LAMBDA_PHOTOPIC {
            assert!((0.0..=1.0).contains(&v));
        }
    }

    #[test]
    fn test_v_lambda_scotopic_peak_near_507() {
        // Scotopic peak near 507nm (index ~25)
        let peak_idx = V_LAMBDA_SCOTOPIC
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .unwrap()
            .0;
        let peak_wl = 380.0 + peak_idx as f64 * 5.0;
        assert!(
            (peak_wl - 505.0).abs() < 10.0,
            "Scotopic peak at {peak_wl}nm, expected ~507nm"
        );
    }

    #[test]
    fn test_luminous_flux_d65_positive() {
        let d65 = illuminant_d65();
        let flux = luminous_flux(&d65);
        assert!(flux > 0.0, "D65 luminous flux should be positive");
    }

    #[test]
    fn test_luminous_efficacy_d65() {
        let d65 = illuminant_d65();
        let eff = luminous_efficacy(&d65);
        // D65 efficacy should be around 200-300 lm/W
        assert!(
            eff > 100.0 && eff < 400.0,
            "D65 efficacy should be ~200-300, got {eff}"
        );
    }

    #[test]
    fn test_luminous_efficacy_max_at_555() {
        // Monochromatic 555nm source → maximum efficacy = 683 lm/W
        let spd = Spd::new(555.0, 5.0, vec![1.0]);
        let eff = luminous_efficacy(&spd);
        assert!(
            (eff - 683.0).abs() < 10.0,
            "555nm efficacy should be ~683, got {eff}"
        );
    }

    #[test]
    fn test_luminous_flux_scotopic_positive() {
        let d65 = illuminant_d65();
        let flux = luminous_flux_scotopic(&d65);
        assert!(flux > 0.0);
    }

    #[test]
    fn test_scotopic_higher_blue_sensitivity() {
        // Scotopic should weight blue more than photopic
        let spd_blue = Spd::new(380.0, 5.0, {
            let mut v = vec![0.0; 81];
            v[10] = 1.0; // 430nm
            v
        });
        let photopic = luminous_flux(&spd_blue);
        let scotopic = luminous_flux_scotopic(&spd_blue);
        assert!(
            scotopic > photopic,
            "Scotopic should be more sensitive to blue: scot={scotopic}, phot={photopic}"
        );
    }
}
