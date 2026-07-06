//! CIE color science: XYZ tristimulus, color matching functions, SPD, illuminants, CRI.

use std::borrow::Cow;

use serde::{Deserialize, Serialize};
use tracing::trace;

use super::{Rgb, planck_radiance};

/// CIE 1931 2° standard observer color matching functions.
///
/// Tabulated at 5nm intervals from 380nm to 780nm (81 entries).
/// Each entry is (x̄, ȳ, z̄).
pub const CIE_1931_2DEG: [(f64, f64, f64); 81] = [
    // 380-400
    (0.001368, 0.000039, 0.006450),
    (0.002236, 0.000064, 0.010550),
    (0.004243, 0.000120, 0.020050),
    (0.007650, 0.000217, 0.036210),
    (0.014310, 0.000396, 0.067850),
    // 405-425
    (0.023190, 0.000640, 0.110200),
    (0.043510, 0.001210, 0.207400),
    (0.077630, 0.002180, 0.371300),
    (0.134380, 0.004000, 0.645600),
    (0.214770, 0.007300, 1.039050),
    // 430-450
    (0.283900, 0.011600, 1.385600),
    (0.328500, 0.016840, 1.622960),
    (0.348280, 0.023000, 1.747060),
    (0.348060, 0.029800, 1.782600),
    (0.336200, 0.038000, 1.772110),
    // 455-475
    (0.318700, 0.048000, 1.744100),
    (0.290800, 0.060000, 1.669200),
    (0.251100, 0.073900, 1.528100),
    (0.195360, 0.090980, 1.287640),
    (0.142100, 0.112600, 1.041900),
    // 480-500
    (0.095640, 0.139020, 0.812950),
    (0.058010, 0.169300, 0.616200),
    (0.032010, 0.208020, 0.465180),
    (0.014700, 0.258600, 0.353300),
    (0.004900, 0.323000, 0.272000),
    // 505-525
    (0.002400, 0.407300, 0.212300),
    (0.009300, 0.503000, 0.158200),
    (0.029100, 0.608200, 0.111700),
    (0.063270, 0.710000, 0.078250),
    (0.109600, 0.793200, 0.057250),
    // 530-550
    (0.165500, 0.862000, 0.042160),
    (0.225750, 0.914850, 0.029840),
    (0.290400, 0.954000, 0.020300),
    (0.359700, 0.980300, 0.013400),
    (0.433450, 0.994950, 0.008750),
    // 555-575
    (0.512050, 1.000000, 0.005750),
    (0.594500, 0.995000, 0.003900),
    (0.678400, 0.978600, 0.002750),
    (0.762100, 0.952000, 0.002100),
    (0.842500, 0.915400, 0.001800),
    // 580-600
    (0.916300, 0.870000, 0.001650),
    (0.978600, 0.816300, 0.001400),
    (1.026300, 0.757000, 0.001100),
    (1.056700, 0.694900, 0.001000),
    (1.062200, 0.631000, 0.000800),
    // 605-625
    (1.045600, 0.566800, 0.000600),
    (1.002600, 0.503000, 0.000340),
    (0.938400, 0.441200, 0.000240),
    (0.854450, 0.381000, 0.000190),
    (0.751400, 0.321000, 0.000100),
    // 630-650
    (0.642400, 0.265000, 0.000050),
    (0.541900, 0.217000, 0.000030),
    (0.447900, 0.175000, 0.000020),
    (0.360800, 0.138200, 0.000010),
    (0.283500, 0.107000, 0.000000),
    // 655-675
    (0.218700, 0.081600, 0.000000),
    (0.164900, 0.061000, 0.000000),
    (0.121200, 0.044580, 0.000000),
    (0.087400, 0.032000, 0.000000),
    (0.063600, 0.023200, 0.000000),
    // 680-700
    (0.046770, 0.017000, 0.000000),
    (0.032900, 0.011920, 0.000000),
    (0.022700, 0.008210, 0.000000),
    (0.015840, 0.005723, 0.000000),
    (0.011359, 0.004102, 0.000000),
    // 705-725
    (0.008111, 0.002929, 0.000000),
    (0.005790, 0.002091, 0.000000),
    (0.004109, 0.001484, 0.000000),
    (0.002899, 0.001047, 0.000000),
    (0.002049, 0.000740, 0.000000),
    // 730-750
    (0.001440, 0.000520, 0.000000),
    (0.001000, 0.000361, 0.000000),
    (0.000690, 0.000249, 0.000000),
    (0.000476, 0.000172, 0.000000),
    (0.000332, 0.000120, 0.000000),
    // 755-780
    (0.000235, 0.000085, 0.000000),
    (0.000166, 0.000060, 0.000000),
    (0.000117, 0.000042, 0.000000),
    (0.000083, 0.000030, 0.000000),
    (0.000059, 0.000021, 0.000000),
    (0.000042, 0.000015, 0.000000),
];

/// CIE 1964 10-degree standard observer color matching functions.
///
/// Tabulated at 5nm intervals from 380nm to 780nm (81 entries).
/// Each entry is (x̄₁₀, ȳ₁₀, z̄₁₀). Source: CIE S 014-1 / cvrl.org.
pub const CIE_1964_10DEG: [(f64, f64, f64); 81] = [
    (0.000160, 0.000017, 0.000705),
    (0.000662, 0.000072, 0.002928),
    (0.002362, 0.000253, 0.010482),
    (0.007242, 0.000769, 0.032344),
    (0.019110, 0.002004, 0.086011),
    (0.043400, 0.004509, 0.197120),
    (0.084736, 0.008756, 0.389366),
    (0.140638, 0.014456, 0.656760),
    (0.204492, 0.021391, 0.972542),
    (0.264737, 0.029497, 1.282500),
    (0.314679, 0.038676, 1.553480),
    (0.357719, 0.049602, 1.798500),
    (0.383734, 0.062077, 1.967280),
    (0.386726, 0.074704, 2.027300),
    (0.370702, 0.089456, 1.994800),
    (0.342957, 0.106256, 1.900700),
    (0.302273, 0.128201, 1.745370),
    (0.254085, 0.152761, 1.554900),
    (0.195618, 0.185190, 1.317560),
    (0.132349, 0.219940, 1.030200),
    (0.080507, 0.253589, 0.772125),
    (0.041072, 0.297665, 0.570060),
    (0.016172, 0.339133, 0.415254),
    (0.005132, 0.395379, 0.302356),
    (0.003816, 0.460777, 0.218502),
    (0.015444, 0.531360, 0.159249),
    (0.037465, 0.606741, 0.112044),
    (0.071358, 0.685660, 0.082248),
    (0.117749, 0.761757, 0.060709),
    (0.172953, 0.823330, 0.043050),
    (0.236491, 0.875211, 0.030451),
    (0.304213, 0.923810, 0.020584),
    (0.376772, 0.961988, 0.013676),
    (0.451584, 0.982200, 0.007918),
    (0.529826, 0.991761, 0.003988),
    (0.616053, 0.999110, 0.001091),
    (0.705224, 0.997340, 0.000000),
    (0.793832, 0.982380, 0.000000),
    (0.878655, 0.955552, 0.000000),
    (0.951162, 0.915175, 0.000000),
    (1.014160, 0.868934, 0.000000),
    (1.074300, 0.825623, 0.000000),
    (1.118520, 0.777405, 0.000000),
    (1.134300, 0.720353, 0.000000),
    (1.123990, 0.658341, 0.000000),
    (1.089100, 0.593878, 0.000000),
    (1.030480, 0.527963, 0.000000),
    (0.950740, 0.461834, 0.000000),
    (0.856297, 0.398057, 0.000000),
    (0.754930, 0.339554, 0.000000),
    (0.647467, 0.283493, 0.000000),
    (0.535110, 0.228254, 0.000000),
    (0.431567, 0.179828, 0.000000),
    (0.343690, 0.140211, 0.000000),
    (0.268329, 0.107633, 0.000000),
    (0.204300, 0.081187, 0.000000),
    (0.152568, 0.060281, 0.000000),
    (0.112210, 0.044096, 0.000000),
    (0.081261, 0.031800, 0.000000),
    (0.057930, 0.022602, 0.000000),
    (0.040851, 0.015905, 0.000000),
    (0.028623, 0.011130, 0.000000),
    (0.019941, 0.007749, 0.000000),
    (0.013842, 0.005375, 0.000000),
    (0.009577, 0.003718, 0.000000),
    (0.006605, 0.002565, 0.000000),
    (0.004553, 0.001768, 0.000000),
    (0.003145, 0.001222, 0.000000),
    (0.002175, 0.000846, 0.000000),
    (0.001506, 0.000586, 0.000000),
    (0.001045, 0.000407, 0.000000),
    (0.000727, 0.000284, 0.000000),
    (0.000508, 0.000199, 0.000000),
    (0.000356, 0.000140, 0.000000),
    (0.000251, 0.000098, 0.000000),
    (0.000178, 0.000070, 0.000000),
    (0.000126, 0.000050, 0.000000),
    (0.000090, 0.000036, 0.000000),
    (0.000065, 0.000025, 0.000000),
    (0.000046, 0.000018, 0.000000),
    (0.000033, 0.000013, 0.000000),
];

/// CIE 2015 2-degree cone-fundamental-based observer.
///
/// Tabulated at 5nm intervals from 380nm to 780nm (81 entries).
/// Source: CIE 170-2:2015 (DOI: 10.25039/CIE.DS.548rw69q).
/// More physiologically accurate than CIE 1931, especially below 460nm.
/// 380nm and 385nm are zero (official data starts at 390nm).
pub const CIE_2015_2DEG: [(f64, f64, f64); 81] = [
    (0.000000, 0.000000, 0.000000),
    (0.000000, 0.000000, 0.000000),
    (0.003770, 0.000415, 0.018473),
    (0.009383, 0.001060, 0.046098),
    (0.022143, 0.002452, 0.109609),
    (0.047430, 0.004972, 0.236925),
    (0.089538, 0.009080, 0.450837),
    (0.144621, 0.014294, 0.737882),
    (0.203573, 0.020274, 1.051821),
    (0.248852, 0.026121, 1.305008),
    (0.291825, 0.033190, 1.552826),
    (0.322709, 0.041579, 1.748280),
    (0.348255, 0.050337, 1.917479),
    (0.341848, 0.057434, 1.918437),
    (0.322464, 0.064724, 1.848545),
    (0.282665, 0.072383, 1.664439),
    (0.248525, 0.085148, 1.522157),
    (0.221978, 0.106015, 1.428440),
    (0.180691, 0.129896, 1.250610),
    (0.129192, 0.153507, 0.999179),
    (0.081829, 0.178805, 0.755238),
    (0.046009, 0.206483, 0.561731),
    (0.020840, 0.237916, 0.409931),
    (0.007098, 0.285068, 0.310594),
    (0.002462, 0.348354, 0.237675),
    (0.003649, 0.427760, 0.172002),
    (0.015570, 0.520497, 0.117680),
    (0.043152, 0.620626, 0.082835),
    (0.079629, 0.718089, 0.056504),
    (0.126847, 0.794645, 0.037519),
    (0.181803, 0.857580, 0.024382),
    (0.240502, 0.907135, 0.015662),
    (0.309812, 0.954468, 0.009846),
    (0.380424, 0.981411, 0.006131),
    (0.449421, 0.989023, 0.003790),
    (0.528023, 0.999461, 0.002327),
    (0.613378, 0.996774, 0.001432),
    (0.701677, 0.990255, 0.000882),
    (0.796775, 0.973261, 0.000545),
    (0.885338, 0.942457, 0.000339),
    (0.963839, 0.896361, 0.000212),
    (1.051011, 0.858720, 0.000134),
    (1.109767, 0.811587, 0.000085),
    (1.143620, 0.754479, 0.000055),
    (1.151033, 0.691855, 0.000035),
    (1.134757, 0.627007, 0.000023),
    (1.083928, 0.558375, 0.000016),
    (1.007344, 0.489595, 0.000010),
    (0.914288, 0.422990, 0.000000),
    (0.813557, 0.360925, 0.000000),
    (0.692472, 0.298087, 0.000000),
    (0.575541, 0.241690, 0.000000),
    (0.473122, 0.194312, 0.000000),
    (0.384499, 0.154740, 0.000000),
    (0.299737, 0.119312, 0.000000),
    (0.227779, 0.089796, 0.000000),
    (0.170791, 0.066710, 0.000000),
    (0.126381, 0.048997, 0.000000),
    (0.092246, 0.035600, 0.000000),
    (0.066400, 0.025542, 0.000000),
    (0.047106, 0.018079, 0.000000),
    (0.032921, 0.012616, 0.000000),
    (0.022623, 0.008661, 0.000000),
    (0.015754, 0.006028, 0.000000),
    (0.010968, 0.004196, 0.000000),
    (0.007609, 0.002911, 0.000000),
    (0.005215, 0.001996, 0.000000),
    (0.003569, 0.001367, 0.000000),
    (0.002465, 0.000945, 0.000000),
    (0.001704, 0.000654, 0.000000),
    (0.001186, 0.000456, 0.000000),
    (0.000827, 0.000318, 0.000000),
    (0.000576, 0.000222, 0.000000),
    (0.000406, 0.000157, 0.000000),
    (0.000286, 0.000110, 0.000000),
    (0.000202, 0.000078, 0.000000),
    (0.000144, 0.000056, 0.000000),
    (0.000102, 0.000040, 0.000000),
    (0.000073, 0.000029, 0.000000),
    (0.000053, 0.000021, 0.000000),
    (0.000038, 0.000015, 0.000000),
];

/// CIE 2015 10-degree cone-fundamental-based observer.
///
/// Tabulated at 5nm intervals from 380nm to 780nm (81 entries).
/// Source: CIE 170-2:2015 (DOI: 10.25039/CIE.DS.dm6qiig7).
/// 380nm and 385nm are zero (official data starts at 390nm).
pub const CIE_2015_10DEG: [(f64, f64, f64); 81] = [
    (0.000000, 0.000000, 0.000000),
    (0.000000, 0.000000, 0.000000),
    (0.002952, 0.000408, 0.013188),
    (0.007641, 0.001078, 0.034246),
    (0.018793, 0.002590, 0.085083),
    (0.042050, 0.005474, 0.192707),
    (0.082773, 0.010413, 0.383282),
    (0.139513, 0.017130, 0.656819),
    (0.207765, 0.025761, 0.993344),
    (0.268899, 0.035296, 1.308674),
    (0.328180, 0.046982, 1.624940),
    (0.369308, 0.060474, 1.867751),
    (0.402619, 0.074683, 2.075946),
    (0.404253, 0.088205, 2.132574),
    (0.393214, 0.103903, 2.128264),
    (0.348221, 0.119539, 1.946651),
    (0.301311, 0.141459, 1.768440),
    (0.253422, 0.170137, 1.582342),
    (0.191418, 0.199986, 1.310576),
    (0.128317, 0.231243, 1.010952),
    (0.075931, 0.268227, 0.751639),
    (0.038368, 0.310944, 0.554962),
    (0.014007, 0.355402, 0.397811),
    (0.003447, 0.414823, 0.290582),
    (0.005652, 0.478048, 0.207816),
    (0.015620, 0.549134, 0.139464),
    (0.037782, 0.624830, 0.088524),
    (0.075389, 0.701229, 0.058245),
    (0.120151, 0.778820, 0.037849),
    (0.175683, 0.837636, 0.024314),
    (0.238025, 0.882955, 0.015395),
    (0.304699, 0.923386, 0.009753),
    (0.384186, 0.966533, 0.006083),
    (0.463311, 0.988689, 0.003769),
    (0.537417, 0.990750, 0.002324),
    (0.623089, 0.999778, 0.001427),
    (0.712385, 0.994430, 0.000878),
    (0.801628, 0.984813, 0.000541),
    (0.893341, 0.964055, 0.000334),
    (0.972130, 0.928650, 0.000208),
    (1.034327, 0.877536, 0.000130),
    (1.106886, 0.837084, 0.000082),
    (1.147304, 0.786995, 0.000052),
    (1.160477, 0.727231, 0.000033),
    (1.148163, 0.662904, 0.000022),
    (1.113846, 0.597038, 0.000014),
    (1.048485, 0.528230, 0.000010),
    (0.961711, 0.460131, 0.000006),
    (0.862958, 0.395076, 0.000000),
    (0.760350, 0.335179, 0.000000),
    (0.641398, 0.275181, 0.000000),
    (0.529098, 0.221956, 0.000000),
    (0.432313, 0.177688, 0.000000),
    (0.349636, 0.141020, 0.000000),
    (0.271490, 0.108400, 0.000000),
    (0.205651, 0.081377, 0.000000),
    (0.153816, 0.060340, 0.000000),
    (0.113607, 0.044254, 0.000000),
    (0.082810, 0.032119, 0.000000),
    (0.059548, 0.023026, 0.000000),
    (0.042215, 0.016288, 0.000000),
    (0.029488, 0.011361, 0.000000),
    (0.020256, 0.007797, 0.000000),
    (0.014102, 0.005425, 0.000000),
    (0.009816, 0.003776, 0.000000),
    (0.006809, 0.002619, 0.000000),
    (0.004666, 0.001796, 0.000000),
    (0.003194, 0.001230, 0.000000),
    (0.002206, 0.000850, 0.000000),
    (0.001525, 0.000588, 0.000000),
    (0.001061, 0.000410, 0.000000),
    (0.000740, 0.000286, 0.000000),
    (0.000515, 0.000199, 0.000000),
    (0.000363, 0.000141, 0.000000),
    (0.000256, 0.000099, 0.000000),
    (0.000181, 0.000070, 0.000000),
    (0.000129, 0.000050, 0.000000),
    (0.000092, 0.000036, 0.000000),
    (0.000066, 0.000026, 0.000000),
    (0.000047, 0.000018, 0.000000),
    (0.000034, 0.000013, 0.000000),
];

/// Standard observer type for generic SPD→XYZ integration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Observer {
    /// CIE 1931 2-degree standard observer (default).
    Cie1931,
    /// CIE 1964 10-degree standard observer.
    Cie1964,
    /// CIE 2015 2-degree cone-fundamental-based observer.
    Cie2015_2deg,
    /// CIE 2015 10-degree cone-fundamental-based observer.
    Cie2015_10deg,
}

impl Observer {
    /// Get the color matching function table for this observer.
    #[must_use]
    pub fn cmf_table(&self) -> &'static [(f64, f64, f64); 81] {
        match self {
            Observer::Cie1931 => &CIE_1931_2DEG,
            Observer::Cie1964 => &CIE_1964_10DEG,
            Observer::Cie2015_2deg => &CIE_2015_2DEG,
            Observer::Cie2015_10deg => &CIE_2015_10DEG,
        }
    }
}

/// CIE XYZ tristimulus values.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Xyz {
    /// X tristimulus value.
    pub x: f64,
    /// Y tristimulus value (luminance).
    pub y: f64,
    /// Z tristimulus value.
    pub z: f64,
}

impl Xyz {
    /// Create new XYZ tristimulus values.
    #[must_use]
    #[inline]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Convert to CIE xyY chromaticity coordinates.
    ///
    /// Returns (x, y, Y) where x,y are chromaticity and Y is luminance.
    #[must_use]
    #[inline]
    pub fn to_xyy(self) -> (f64, f64, f64) {
        let sum = self.x + self.y + self.z;
        if sum < 1e-15 {
            return (0.0, 0.0, 0.0);
        }
        (self.x / sum, self.y / sum, self.y)
    }

    /// Create from xyY chromaticity coordinates.
    #[must_use]
    #[inline]
    pub fn from_xyy(cx: f64, cy: f64, big_y: f64) -> Self {
        if cy < 1e-15 {
            return Self::new(0.0, 0.0, 0.0);
        }
        Self {
            x: (cx / cy) * big_y,
            y: big_y,
            z: ((1.0 - cx - cy) / cy) * big_y,
        }
    }

    /// Convert to linear sRGB (not gamma-corrected).
    ///
    /// Uses the sRGB/Rec.709 conversion matrix (D65 white point).
    #[must_use]
    #[inline]
    pub fn to_linear_srgb(self) -> Rgb {
        Rgb::new(
            3.2404542 * self.x - 1.5371385 * self.y - 0.4985314 * self.z,
            -0.9692660 * self.x + 1.8760108 * self.y + 0.0415560 * self.z,
            0.0556434 * self.x - 0.2040259 * self.y + 1.0572252 * self.z,
        )
    }

    /// Convert to sRGB with gamma correction and clamping.
    #[must_use]
    pub fn to_srgb(self) -> Rgb {
        let linear = self.to_linear_srgb();
        Rgb::new(
            linear_to_srgb_gamma(linear.r),
            linear_to_srgb_gamma(linear.g),
            linear_to_srgb_gamma(linear.b),
        )
        .clamp()
    }

    /// D65 standard illuminant white point.
    pub const D65_WHITE: Xyz = Xyz::new(0.95047, 1.0, 1.08883);
    /// D50 standard illuminant white point.
    pub const D50_WHITE: Xyz = Xyz::new(0.96422, 1.0, 0.82521);
}

/// Convert linear sRGB to XYZ (inverse of the sRGB matrix).
#[must_use]
#[inline]
pub fn linear_srgb_to_xyz(rgb: &Rgb) -> Xyz {
    Xyz::new(
        0.4124564 * rgb.r + 0.3575761 * rgb.g + 0.1804375 * rgb.b,
        0.2126729 * rgb.r + 0.7151522 * rgb.g + 0.0721750 * rgb.b,
        0.0193339 * rgb.r + 0.1191920 * rgb.g + 0.9503041 * rgb.b,
    )
}

/// sRGB gamma correction: linear → sRGB.
#[must_use]
#[inline]
pub fn linear_to_srgb_gamma(c: f64) -> f64 {
    if c <= 0.003_130_8 {
        12.92 * c
    } else {
        1.055 * c.powf(1.0 / 2.4) - 0.055
    }
}

/// sRGB inverse gamma: sRGB → linear.
#[must_use]
#[inline]
pub fn srgb_gamma_to_linear(c: f64) -> f64 {
    if c <= 0.040_45 {
        c / 12.92
    } else {
        ((c + 0.055) / 1.055).powf(2.4)
    }
}

/// Correlated color temperature from CIE xy chromaticity (McCamy's approximation).
///
/// Valid for ~3000K–50000K. Input is CIE 1931 (x, y) chromaticity.
/// Returns temperature in Kelvin.
#[must_use]
#[inline]
pub fn cct_from_xy(cx: f64, cy: f64) -> f64 {
    let n = (cx - 0.3320) / (0.1858 - cy);
    449.0 * n * n * n + 3525.0 * n * n + 6823.3 * n + 5520.33
}

/// Interpolate CIE 1931 CMF at an arbitrary wavelength (nm).
///
/// Uses linear interpolation between the 5nm tabulated values.
/// Returns (x̄, ȳ, z̄). Returns (0,0,0) outside 380–780nm.
#[must_use]
#[inline]
pub fn cie_cmf_at(wavelength_nm: f64) -> (f64, f64, f64) {
    if !(380.0..=780.0).contains(&wavelength_nm) {
        return (0.0, 0.0, 0.0);
    }
    let idx_f = (wavelength_nm - 380.0) / 5.0;
    let idx = idx_f as usize;
    if idx >= 80 {
        return CIE_1931_2DEG[80];
    }
    let frac = idx_f - idx as f64;
    let (x0, y0, z0) = CIE_1931_2DEG[idx];
    let (x1, y1, z1) = CIE_1931_2DEG[idx + 1];
    (
        x0 + frac * (x1 - x0),
        y0 + frac * (y1 - y0),
        z0 + frac * (z1 - z0),
    )
}

// ── Spectral Power Distribution ───────────────────────────────────────────

/// A spectral power distribution (SPD).
///
/// Stores relative power at uniform wavelength intervals. Uses
/// [`Cow`] storage so that standard illuminant data can be borrowed
/// from static arrays without allocating.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Spd {
    /// Start wavelength in nm.
    pub start_nm: f64,
    /// Wavelength step in nm.
    pub step_nm: f64,
    /// Relative power values (borrowed from static data or owned).
    pub values: Cow<'static, [f64]>,
}

impl Spd {
    /// Create a new SPD from owned samples.
    #[must_use]
    pub fn new(start_nm: f64, step_nm: f64, values: Vec<f64>) -> Self {
        Self {
            start_nm,
            step_nm,
            values: Cow::Owned(values),
        }
    }

    /// Create a new SPD that borrows from a static slice (zero allocation).
    #[must_use]
    pub const fn from_static(start_nm: f64, step_nm: f64, values: &'static [f64]) -> Self {
        Self {
            start_nm,
            step_nm,
            values: Cow::Borrowed(values),
        }
    }

    /// End wavelength (inclusive).
    #[must_use]
    #[inline]
    pub fn end_nm(&self) -> f64 {
        self.start_nm + self.step_nm * (self.values.len() as f64 - 1.0)
    }

    /// Interpolate the SPD at an arbitrary wavelength.
    #[must_use]
    #[inline]
    pub fn at(&self, wavelength_nm: f64) -> f64 {
        if wavelength_nm < self.start_nm || wavelength_nm > self.end_nm() {
            return 0.0;
        }
        let idx_f = (wavelength_nm - self.start_nm) / self.step_nm;
        let idx = idx_f as usize;
        if idx >= self.values.len() - 1 {
            return *self.values.last().unwrap_or(&0.0);
        }
        let frac = idx_f - idx as f64;
        self.values[idx] + frac * (self.values[idx + 1] - self.values[idx])
    }

    /// Integrate this SPD against the CIE 1931 CMFs to get XYZ tristimulus values.
    ///
    /// Uses 5nm integration steps over the visible range.
    #[must_use]
    pub fn to_xyz(&self) -> Xyz {
        trace!(
            start_nm = self.start_nm,
            step_nm = self.step_nm,
            samples = self.values.len(),
            "spd_to_xyz"
        );
        let mut x = 0.0;
        let mut y = 0.0;
        let mut z = 0.0;

        // Fast path: SPD starts at 380nm with 5nm step (aligned with CMF table)
        let aligned = (self.start_nm - 380.0).abs() < 1e-10 && (self.step_nm - 5.0).abs() < 1e-10;

        for (i, &(cx, cy, cz)) in CIE_1931_2DEG.iter().enumerate() {
            let power = if aligned && i < self.values.len() {
                self.values[i]
            } else {
                self.at(380.0 + i as f64 * 5.0)
            };
            x += power * cx;
            y += power * cy;
            z += power * cz;
        }
        // Integration weight = step size
        Xyz::new(x * 5.0, y * 5.0, z * 5.0)
    }

    /// Integrate this SPD against a specific observer's CMFs to get XYZ tristimulus values.
    ///
    /// Supports CIE 1931, 1964, and 2015 observers via the [`Observer`] enum.
    #[must_use]
    pub fn to_xyz_observer(&self, observer: Observer) -> Xyz {
        let cmf = observer.cmf_table();
        let mut x = 0.0;
        let mut y = 0.0;
        let mut z = 0.0;
        let aligned = (self.start_nm - 380.0).abs() < 1e-10 && (self.step_nm - 5.0).abs() < 1e-10;
        for (i, &(cx, cy, cz)) in cmf.iter().enumerate() {
            let power = if aligned && i < self.values.len() {
                self.values[i]
            } else {
                self.at(380.0 + i as f64 * 5.0)
            };
            x += power * cx;
            y += power * cy;
            z += power * cz;
        }
        Xyz::new(x * 5.0, y * 5.0, z * 5.0)
    }

    /// Convert this SPD to sRGB color.
    #[must_use]
    pub fn to_srgb(&self) -> Rgb {
        let xyz = self.to_xyz();
        // Normalize to Y=1 for white point
        let norm = if xyz.y > 1e-15 { 1.0 / xyz.y } else { 1.0 };
        let normalized = Xyz::new(xyz.x * norm, xyz.y * norm, xyz.z * norm);
        normalized.to_srgb()
    }

    /// Blackbody SPD at a given temperature, sampled at 5nm from 380–780nm.
    #[must_use]
    pub fn blackbody(temperature_k: f64) -> Self {
        let mut values = Vec::with_capacity(81);
        let mut wl = 380.0;
        while wl <= 780.0 {
            values.push(planck_radiance(wl * 1e-9, temperature_k));
            wl += 5.0;
        }
        Self::new(380.0, 5.0, values)
    }
}

// ── Standard Illuminants ──────────────────────────────────────────────────

/// D65 illuminant SPD data (5nm intervals, 380–780nm, 81 values).
pub const ILLUMINANT_D65_DATA: [f64; 81] = [
    49.98, 52.31, 54.65, 68.70, 82.75, 87.12, 91.49, 92.46, 93.43, 90.06, 86.68, 95.77, 104.86,
    110.94, 117.01, 117.41, 117.81, 116.34, 114.86, 115.39, 115.92, 112.37, 108.81, 109.08, 109.35,
    108.58, 107.80, 106.30, 104.79, 106.24, 107.69, 106.05, 104.41, 104.23, 104.05, 102.02, 100.00,
    98.17, 96.33, 96.06, 95.79, 92.24, 88.69, 89.35, 90.01, 89.80, 89.60, 88.65, 87.70, 85.49,
    83.29, 83.49, 83.69, 81.86, 80.03, 80.12, 80.21, 81.25, 82.28, 80.28, 78.28, 74.00, 69.72,
    70.67, 71.61, 72.98, 74.35, 67.98, 61.60, 65.74, 69.89, 72.49, 75.09, 69.34, 63.59, 55.01,
    46.42, 56.61, 66.81, 65.09, 63.38,
];

/// D50 illuminant SPD data (5nm intervals, 380–780nm, 81 values).
pub const ILLUMINANT_D50_DATA: [f64; 81] = [
    24.49, 27.18, 29.87, 39.59, 49.31, 52.91, 56.51, 58.27, 60.03, 58.93, 57.82, 66.32, 74.82,
    81.04, 87.25, 88.93, 90.61, 90.99, 91.37, 93.24, 95.11, 93.54, 91.96, 93.84, 95.72, 96.17,
    96.61, 96.87, 97.13, 99.61, 102.10, 101.43, 100.75, 101.54, 102.32, 101.16, 100.00, 98.87,
    97.74, 98.33, 98.92, 96.21, 93.50, 95.59, 97.69, 98.48, 99.27, 99.16, 99.04, 97.38, 95.72,
    97.29, 98.86, 97.26, 95.67, 96.93, 98.19, 100.60, 103.00, 101.07, 99.13, 93.26, 87.38, 89.49,
    91.60, 92.25, 92.89, 84.87, 76.85, 81.68, 86.51, 89.55, 92.58, 85.40, 78.23, 67.96, 57.69,
    70.31, 82.92, 80.60, 78.27,
];

/// F2 illuminant SPD data (5nm intervals, 380–780nm, 81 values).
pub const ILLUMINANT_F2_DATA: [f64; 81] = [
    1.18, 1.48, 1.84, 2.15, 3.44, 15.69, 3.85, 3.74, 4.19, 4.62, 5.06, 34.98, 11.81, 6.27, 6.63,
    6.93, 7.19, 7.40, 7.54, 7.62, 7.65, 7.62, 7.62, 7.45, 7.28, 7.15, 7.05, 7.04, 7.16, 7.47, 8.04,
    8.88, 10.01, 24.24, 16.55, 14.71, 16.73, 19.88, 28.75, 16.52, 14.36, 13.45, 12.72, 12.27,
    12.26, 31.85, 12.75, 12.68, 12.00, 11.17, 10.41, 9.87, 9.26, 8.59, 7.96, 7.41, 6.87, 6.30,
    5.82, 5.27, 4.81, 4.43, 4.07, 3.71, 3.38, 3.09, 2.84, 2.64, 2.43, 2.28, 2.08, 1.91, 1.77, 1.63,
    1.53, 1.42, 1.32, 1.24, 1.15, 1.08, 1.01,
];

/// F11 illuminant SPD data (5nm intervals, 380–780nm, 81 values).
pub const ILLUMINANT_F11_DATA: [f64; 81] = [
    0.91, 0.63, 0.46, 0.37, 1.29, 12.68, 0.96, 0.64, 0.45, 0.37, 1.00, 7.77, 16.46, 1.86, 0.76,
    0.49, 0.38, 0.36, 0.39, 0.45, 0.55, 0.66, 0.77, 0.90, 1.05, 1.18, 1.32, 1.47, 1.64, 1.86, 2.14,
    2.52, 3.02, 3.77, 5.04, 15.39, 40.98, 72.84, 80.01, 45.35, 19.44, 8.99, 4.70, 3.22, 2.82, 3.73,
    7.39, 74.89, 21.44, 11.69, 4.72, 2.82, 2.00, 1.52, 1.22, 1.10, 1.09, 1.07, 0.99, 0.88, 0.76,
    0.68, 0.63, 0.59, 0.55, 0.52, 0.49, 0.46, 0.44, 0.42, 0.39, 0.37, 0.35, 0.33, 0.31, 0.30, 0.29,
    0.28, 0.27, 0.26, 0.25,
];

/// CIE Standard Illuminant D65 (daylight, ~6504K).
///
/// Relative SPD at 5nm intervals, 380–780nm. Borrows static data (zero allocation).
#[must_use]
pub const fn illuminant_d65() -> Spd {
    Spd::from_static(380.0, 5.0, &ILLUMINANT_D65_DATA)
}

/// CIE Standard Illuminant D50 (horizon daylight, ~5003K).
///
/// Borrows static data (zero allocation).
#[must_use]
pub const fn illuminant_d50() -> Spd {
    Spd::from_static(380.0, 5.0, &ILLUMINANT_D50_DATA)
}

/// CIE Standard Illuminant A (incandescent, ~2856K).
pub fn illuminant_a() -> Spd {
    // Planck formula normalized to 100 at 560nm
    let mut values = Vec::with_capacity(81);
    let norm = planck_radiance(560e-9, 2856.0);
    let mut wl = 380.0;
    while wl <= 780.0 {
        values.push(100.0 * planck_radiance(wl * 1e-9, 2856.0) / norm);
        wl += 5.0;
    }
    Spd::new(380.0, 5.0, values)
}

/// CIE Standard Illuminant F2 (cool white fluorescent).
///
/// Borrows static data (zero allocation).
#[must_use]
pub const fn illuminant_f2() -> Spd {
    Spd::from_static(380.0, 5.0, &ILLUMINANT_F2_DATA)
}

/// CIE Standard Illuminant F11 (narrow-band fluorescent, ~4000K).
///
/// Borrows static data (zero allocation).
#[must_use]
pub const fn illuminant_f11() -> Spd {
    Spd::from_static(380.0, 5.0, &ILLUMINANT_F11_DATA)
}

// ── CRI (Color Rendering Index) ───────────────────────────────────────────

/// Compute the Color Rendering Index (CRI Ra) of a test illuminant.
///
/// Compares the test illuminant against a reference (blackbody for CCT < 5000K,
/// daylight for CCT >= 5000K). Returns Ra (0–100, 100 = perfect rendering).
///
/// This is a simplified CRI that uses the chromaticity shift method.
#[must_use]
pub fn color_rendering_index(test_spd: &Spd) -> f64 {
    trace!("color_rendering_index");
    let test_xyz = test_spd.to_xyz();
    let (tx, ty, _) = test_xyz.to_xyy();
    let test_cct = cct_from_xy(tx, ty);

    // Reference illuminant: blackbody for warm, D-illuminant for cool
    let ref_spd = if test_cct < 5000.0 {
        Spd::blackbody(test_cct)
    } else {
        // Use D65 as approximation for daylight reference
        illuminant_d65()
    };

    let ref_xyz = ref_spd.to_xyz();

    // Simplified CRI: compare XYZ ratios for 8 spectral bands
    // (This is a simplified version — full CRI uses 14 TCS reflectance samples)
    let mut delta_e_sum = 0.0;
    let bands: [(f64, f64); 8] = [
        (400.0, 450.0),
        (450.0, 490.0),
        (490.0, 520.0),
        (520.0, 560.0),
        (560.0, 600.0),
        (600.0, 640.0),
        (640.0, 690.0),
        (690.0, 750.0),
    ];

    for (start, end) in &bands {
        let mut test_y = 0.0;
        let mut ref_y = 0.0;
        let mut wl = *start;
        while wl <= *end {
            // Use table directly for 5nm-aligned wavelengths in visible range
            let idx = ((wl - 380.0) / 5.0) as usize;
            let cy = if idx < CIE_1931_2DEG.len() {
                CIE_1931_2DEG[idx].1
            } else {
                0.0
            };
            test_y += test_spd.at(wl) * cy;
            ref_y += ref_spd.at(wl) * cy;
            wl += 5.0;
        }
        // Normalize
        let test_norm = if test_xyz.y > 1e-15 {
            test_y / test_xyz.y
        } else {
            0.0
        };
        let ref_norm = if ref_xyz.y > 1e-15 {
            ref_y / ref_xyz.y
        } else {
            0.0
        };
        delta_e_sum += (test_norm - ref_norm).abs();
    }

    let delta_e_avg = delta_e_sum / 8.0;
    // Scale: CRI Ra = 100 - 4.6 * ΔE (approximate)
    (100.0 - 4.6 * delta_e_avg * 100.0).clamp(0.0, 100.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-6;

    #[test]
    fn test_cie_cmf_table_length() {
        assert_eq!(CIE_1931_2DEG.len(), 81); // 380..=780 at 5nm
    }

    #[test]
    fn test_cie_cmf_peak_y_at_555nm() {
        // ȳ peaks near 555nm (index 35)
        let (_, y_peak, _) = CIE_1931_2DEG[35]; // 555nm
        assert!(
            (y_peak - 1.0).abs() < 0.01,
            "ȳ should peak at ~1.0 near 555nm, got {y_peak}"
        );
    }

    #[test]
    fn test_cie_cmf_all_non_negative() {
        for (i, (x, y, z)) in CIE_1931_2DEG.iter().enumerate() {
            let wl = 380.0 + i as f64 * 5.0;
            assert!(*x >= 0.0, "x̄ negative at {wl}nm");
            assert!(*y >= 0.0, "ȳ negative at {wl}nm");
            assert!(*z >= 0.0, "z̄ negative at {wl}nm");
        }
    }

    #[test]
    fn test_cie_cmf_interpolation() {
        // At exact 5nm boundaries, should match table
        let (x, y, z) = cie_cmf_at(555.0);
        let (tx, ty, tz) = CIE_1931_2DEG[35];
        assert!((x - tx).abs() < EPS);
        assert!((y - ty).abs() < EPS);
        assert!((z - tz).abs() < EPS);
    }

    #[test]
    fn test_cie_cmf_interpolation_midpoint() {
        // Between table points, should interpolate
        let (x, _, _) = cie_cmf_at(557.5);
        let (x0, _, _) = CIE_1931_2DEG[35]; // 555
        let (x1, _, _) = CIE_1931_2DEG[36]; // 560
        let expected = (x0 + x1) / 2.0;
        assert!((x - expected).abs() < 0.001);
    }

    #[test]
    fn test_cie_cmf_out_of_range() {
        assert_eq!(cie_cmf_at(300.0), (0.0, 0.0, 0.0));
        assert_eq!(cie_cmf_at(800.0), (0.0, 0.0, 0.0));
    }

    // ── XYZ conversion tests ──────────────────────────────────────────────

    #[test]
    fn test_xyz_to_xyy_roundtrip() {
        let xyz = Xyz::new(0.5, 0.4, 0.3);
        let (cx, cy, big_y) = xyz.to_xyy();
        let back = Xyz::from_xyy(cx, cy, big_y);
        assert!((back.x - xyz.x).abs() < 0.001);
        assert!((back.y - xyz.y).abs() < 0.001);
        assert!((back.z - xyz.z).abs() < 0.001);
    }

    #[test]
    fn test_xyz_to_xyy_white() {
        let (cx, cy, _) = Xyz::D65_WHITE.to_xyy();
        // D65 chromaticity: x ≈ 0.3127, y ≈ 0.3290
        assert!((cx - 0.3127).abs() < 0.001);
        assert!((cy - 0.3290).abs() < 0.001);
    }

    #[test]
    fn test_xyz_to_xyy_black() {
        let (cx, cy, big_y) = Xyz::new(0.0, 0.0, 0.0).to_xyy();
        assert_eq!(cx, 0.0);
        assert_eq!(cy, 0.0);
        assert_eq!(big_y, 0.0);
    }

    #[test]
    fn test_xyz_to_srgb_white() {
        let rgb = Xyz::D65_WHITE.to_srgb();
        // D65 should map to approximately white
        assert!((rgb.r - 1.0).abs() < 0.05, "r={}", rgb.r);
        assert!((rgb.g - 1.0).abs() < 0.05, "g={}", rgb.g);
        assert!((rgb.b - 1.0).abs() < 0.05, "b={}", rgb.b);
    }

    #[test]
    fn test_xyz_to_srgb_black() {
        let rgb = Xyz::new(0.0, 0.0, 0.0).to_srgb();
        assert!(rgb.r.abs() < EPS);
        assert!(rgb.g.abs() < EPS);
        assert!(rgb.b.abs() < EPS);
    }

    #[test]
    fn test_linear_srgb_roundtrip() {
        let rgb = Rgb::new(0.5, 0.3, 0.8);
        let xyz = linear_srgb_to_xyz(&rgb);
        let back = xyz.to_linear_srgb();
        assert!((back.r - rgb.r).abs() < 0.001);
        assert!((back.g - rgb.g).abs() < 0.001);
        assert!((back.b - rgb.b).abs() < 0.001);
    }

    #[test]
    fn test_srgb_gamma_roundtrip() {
        for v in [0.0, 0.01, 0.1, 0.5, 0.9, 1.0] {
            let gamma = linear_to_srgb_gamma(v);
            let back = srgb_gamma_to_linear(gamma);
            assert!((back - v).abs() < 0.001, "Gamma roundtrip failed at {v}");
        }
    }

    #[test]
    fn test_srgb_gamma_monotonic() {
        let mut prev = 0.0;
        for i in 0..=100 {
            let v = i as f64 / 100.0;
            let g = linear_to_srgb_gamma(v);
            assert!(g >= prev - EPS, "Gamma not monotonic at {v}");
            prev = g;
        }
    }

    // ── CCT tests ─────────────────────────────────────────────────────────

    #[test]
    fn test_cct_d65() {
        let (cx, cy, _) = Xyz::D65_WHITE.to_xyy();
        let cct = cct_from_xy(cx, cy);
        assert!(
            (cct - 6504.0).abs() < 100.0,
            "D65 CCT should be ≈6504K, got {cct}"
        );
    }

    #[test]
    fn test_cct_warm() {
        // Approximate warm white chromaticity (x≈0.45, y≈0.41)
        let cct = cct_from_xy(0.45, 0.41);
        assert!(
            cct < 4000.0 && cct > 2000.0,
            "Warm chromaticity CCT ≈ 2800K, got {cct}"
        );
    }

    // ── SPD tests ─────────────────────────────────────────────────────────

    #[test]
    fn test_spd_blackbody_length() {
        let spd = Spd::blackbody(5778.0);
        assert_eq!(spd.values.len(), 81);
        assert!((spd.start_nm - 380.0).abs() < EPS);
        assert!((spd.step_nm - 5.0).abs() < EPS);
    }

    #[test]
    fn test_spd_blackbody_all_positive() {
        let spd = Spd::blackbody(5778.0);
        for (i, v) in spd.values.iter().enumerate() {
            assert!(*v > 0.0, "Blackbody power negative at index {i}");
        }
    }

    #[test]
    fn test_spd_interpolation() {
        let spd = Spd::new(380.0, 5.0, vec![1.0, 2.0, 3.0]);
        assert!((spd.at(380.0) - 1.0).abs() < EPS);
        assert!((spd.at(382.5) - 1.5).abs() < EPS);
        assert!((spd.at(385.0) - 2.0).abs() < EPS);
        assert!((spd.at(390.0) - 3.0).abs() < EPS);
    }

    #[test]
    fn test_spd_out_of_range() {
        let spd = Spd::new(380.0, 5.0, vec![1.0, 2.0, 3.0]);
        assert!((spd.at(370.0)).abs() < EPS);
        assert!((spd.at(400.0)).abs() < EPS);
    }

    #[test]
    fn test_spd_to_xyz_positive() {
        let spd = Spd::blackbody(5778.0);
        let xyz = spd.to_xyz();
        assert!(xyz.x > 0.0);
        assert!(xyz.y > 0.0);
        assert!(xyz.z > 0.0);
    }

    #[test]
    fn test_spd_to_srgb_blackbody_reasonable() {
        let rgb = Spd::blackbody(5778.0).to_srgb();
        // Sun should be roughly warm white
        assert!((0.0..=1.0).contains(&rgb.r), "r={}", rgb.r);
        assert!((0.0..=1.0).contains(&rgb.g), "g={}", rgb.g);
        assert!((0.0..=1.0).contains(&rgb.b), "b={}", rgb.b);
    }

    // ── Illuminant tests ──────────────────────────────────────────────────

    #[test]
    fn test_illuminant_d65_length() {
        let d65 = illuminant_d65();
        assert_eq!(d65.values.len(), 81);
    }

    #[test]
    fn test_illuminant_d65_normalized_at_560() {
        let d65 = illuminant_d65();
        let val = d65.at(560.0);
        assert!(
            (val - 100.0).abs() < 2.0,
            "D65 should be ~100 at 560nm, got {val}"
        );
    }

    #[test]
    fn test_illuminant_d50_length() {
        assert_eq!(illuminant_d50().values.len(), 81);
    }

    #[test]
    fn test_illuminant_a_warm() {
        let a = illuminant_a();
        // Illuminant A should be red-heavy (more power at long wavelengths)
        let blue = a.at(450.0);
        let red = a.at(650.0);
        assert!(red > blue, "Illuminant A should have more red than blue");
    }

    #[test]
    fn test_illuminant_f2_length() {
        assert_eq!(illuminant_f2().values.len(), 81);
    }

    #[test]
    fn test_illuminant_f11_length() {
        assert_eq!(illuminant_f11().values.len(), 81);
    }

    #[test]
    fn test_all_illuminants_positive() {
        let illuminants: Vec<(&str, Spd)> = vec![
            ("D65", illuminant_d65()),
            ("D50", illuminant_d50()),
            ("A", illuminant_a()),
            ("F2", illuminant_f2()),
            ("F11", illuminant_f11()),
        ];
        for (name, spd) in &illuminants {
            for (i, v) in spd.values.iter().enumerate() {
                assert!(*v >= 0.0, "{name} has negative value at index {i}: {v}");
            }
        }
    }

    #[test]
    fn test_all_illuminants_to_xyz_positive() {
        let illuminants: Vec<(&str, Spd)> = vec![
            ("D65", illuminant_d65()),
            ("D50", illuminant_d50()),
            ("A", illuminant_a()),
            ("F2", illuminant_f2()),
            ("F11", illuminant_f11()),
        ];
        for (name, spd) in &illuminants {
            let xyz = spd.to_xyz();
            assert!(xyz.x > 0.0, "{name}: X should be positive");
            assert!(xyz.y > 0.0, "{name}: Y should be positive");
            assert!(xyz.z > 0.0, "{name}: Z should be positive");
        }
    }

    // ── CRI tests ─────────────────────────────────────────────────────────

    #[test]
    fn test_cri_d65_perfect() {
        // D65 compared to itself (or daylight reference) should score high
        let cri = color_rendering_index(&illuminant_d65());
        assert!(cri > 90.0, "D65 CRI should be very high, got {cri}");
    }

    #[test]
    fn test_cri_blackbody_high() {
        // Blackbody at 3000K should score ~100 (it IS the reference)
        let cri = color_rendering_index(&Spd::blackbody(3000.0));
        assert!(cri > 85.0, "Blackbody CRI should be high, got {cri}");
    }

    #[test]
    fn test_cri_range() {
        let illuminants = [
            illuminant_d65(),
            illuminant_d50(),
            illuminant_a(),
            illuminant_f2(),
            illuminant_f11(),
        ];
        for spd in &illuminants {
            let cri = color_rendering_index(spd);
            assert!((0.0..=100.0).contains(&cri), "CRI out of range: {cri}");
        }
    }

    #[test]
    fn test_cri_fluorescent_lower_than_daylight() {
        let cri_d65 = color_rendering_index(&illuminant_d65());
        let cri_f11 = color_rendering_index(&illuminant_f11());
        assert!(
            cri_d65 > cri_f11,
            "Daylight CRI ({cri_d65}) should exceed narrow-band fluorescent ({cri_f11})"
        );
    }

    // ── Spd Cow storage tests ─────────────────────────────────────────────

    #[test]
    fn test_spd_from_static_borrows() {
        let spd = Spd::from_static(380.0, 5.0, &ILLUMINANT_D65_DATA);
        assert!(matches!(spd.values, std::borrow::Cow::Borrowed(_)));
        assert_eq!(spd.values.len(), 81);
    }

    #[test]
    fn test_spd_new_owns() {
        let spd = Spd::new(380.0, 5.0, vec![1.0, 2.0, 3.0]);
        assert!(matches!(spd.values, std::borrow::Cow::Owned(_)));
    }

    #[test]
    fn test_illuminant_d65_is_borrowed() {
        let d65 = illuminant_d65();
        assert!(
            matches!(d65.values, std::borrow::Cow::Borrowed(_)),
            "illuminant_d65 should borrow static data"
        );
    }

    #[test]
    fn test_illuminant_d50_is_borrowed() {
        let d50 = illuminant_d50();
        assert!(matches!(d50.values, std::borrow::Cow::Borrowed(_)));
    }

    #[test]
    fn test_illuminant_f2_is_borrowed() {
        let f2 = illuminant_f2();
        assert!(matches!(f2.values, std::borrow::Cow::Borrowed(_)));
    }

    #[test]
    fn test_illuminant_f11_is_borrowed() {
        let f11 = illuminant_f11();
        assert!(matches!(f11.values, std::borrow::Cow::Borrowed(_)));
    }

    #[test]
    fn test_spd_blackbody_is_owned() {
        let bb = Spd::blackbody(5778.0);
        assert!(matches!(bb.values, std::borrow::Cow::Owned(_)));
    }

    #[test]
    fn test_spd_cow_clone_borrowed_stays_borrowed() {
        let d65 = illuminant_d65();
        let cloned = d65.clone();
        assert!(matches!(cloned.values, std::borrow::Cow::Borrowed(_)));
    }

    #[test]
    fn test_spd_cow_serde_roundtrip() {
        let d65 = illuminant_d65();
        let json = serde_json::to_string(&d65).unwrap();
        let back: Spd = serde_json::from_str(&json).unwrap();
        assert_eq!(back.values.len(), 81);
        assert!((back.values[0] - d65.values[0]).abs() < EPS);
    }

    #[test]
    fn test_spd_cow_to_xyz_identical() {
        // Borrowed and owned should produce identical XYZ
        let borrowed = illuminant_d65();
        let owned = Spd::new(380.0, 5.0, ILLUMINANT_D65_DATA.to_vec());
        let xyz_b = borrowed.to_xyz();
        let xyz_o = owned.to_xyz();
        assert!((xyz_b.x - xyz_o.x).abs() < EPS);
        assert!((xyz_b.y - xyz_o.y).abs() < EPS);
        assert!((xyz_b.z - xyz_o.z).abs() < EPS);
    }

    // ── Observer table tests ──────────────────────────────────────────────

    #[test]
    fn test_cie_1964_table_length() {
        assert_eq!(CIE_1964_10DEG.len(), 81);
    }

    #[test]
    fn test_cie_1964_all_non_negative() {
        for (i, (x, y, z)) in CIE_1964_10DEG.iter().enumerate() {
            let wl = 380.0 + i as f64 * 5.0;
            assert!(*x >= 0.0, "x̄₁₀ negative at {wl}nm");
            assert!(*y >= 0.0, "ȳ₁₀ negative at {wl}nm");
            assert!(*z >= 0.0, "z̄₁₀ negative at {wl}nm");
        }
    }

    #[test]
    fn test_cie_1964_y_peaks_near_555() {
        let peak_idx = CIE_1964_10DEG
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.1.partial_cmp(&b.1.1).unwrap())
            .unwrap()
            .0;
        let peak_wl = 380.0 + peak_idx as f64 * 5.0;
        assert!(
            (peak_wl - 555.0).abs() < 10.0,
            "10-deg ȳ peak at {peak_wl}nm, expected ~555nm"
        );
    }

    #[test]
    fn test_cie_2015_2deg_table_length() {
        assert_eq!(CIE_2015_2DEG.len(), 81);
    }

    #[test]
    fn test_cie_2015_2deg_starts_zero() {
        // 380nm and 385nm should be zero (below CIE 170-2 range)
        let (x, y, z) = CIE_2015_2DEG[0];
        assert!(x.abs() < EPS && y.abs() < EPS && z.abs() < EPS);
    }

    #[test]
    fn test_cie_2015_10deg_table_length() {
        assert_eq!(CIE_2015_10DEG.len(), 81);
    }

    #[test]
    fn test_observer_enum_returns_correct_table() {
        assert_eq!(Observer::Cie1931.cmf_table().len(), 81);
        assert_eq!(Observer::Cie1964.cmf_table().len(), 81);
        assert_eq!(Observer::Cie2015_2deg.cmf_table().len(), 81);
        assert_eq!(Observer::Cie2015_10deg.cmf_table().len(), 81);
    }

    #[test]
    fn test_spd_to_xyz_observer_default_matches() {
        // Observer::Cie1931 should produce same result as to_xyz()
        let d65 = illuminant_d65();
        let xyz_default = d65.to_xyz();
        let xyz_observer = d65.to_xyz_observer(Observer::Cie1931);
        assert!((xyz_default.x - xyz_observer.x).abs() < EPS);
        assert!((xyz_default.y - xyz_observer.y).abs() < EPS);
        assert!((xyz_default.z - xyz_observer.z).abs() < EPS);
    }

    #[test]
    fn test_spd_to_xyz_1964_differs_from_1931() {
        let d65 = illuminant_d65();
        let xyz_1931 = d65.to_xyz_observer(Observer::Cie1931);
        let xyz_1964 = d65.to_xyz_observer(Observer::Cie1964);
        // Should be similar but not identical
        assert!(
            (xyz_1931.x - xyz_1964.x).abs() > 0.1,
            "1931 and 1964 should differ"
        );
    }

    #[test]
    fn test_all_observers_produce_positive_xyz_for_d65() {
        let d65 = illuminant_d65();
        let observers = [
            Observer::Cie1931,
            Observer::Cie1964,
            Observer::Cie2015_2deg,
            Observer::Cie2015_10deg,
        ];
        for obs in &observers {
            let xyz = d65.to_xyz_observer(*obs);
            assert!(xyz.x > 0.0, "{obs:?}: X should be positive");
            assert!(xyz.y > 0.0, "{obs:?}: Y should be positive");
            assert!(xyz.z > 0.0, "{obs:?}: Z should be positive");
        }
    }
}
