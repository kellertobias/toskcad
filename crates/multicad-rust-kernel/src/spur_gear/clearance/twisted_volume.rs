//! Common material volume of two twisted gears on parallel axes.
//!
//! At height `z` each helical or herringbone gear is its exact spur section
//! turned about its own axis by `phase_at(z)`, so the section of the common
//! material is the exact common region of the two turned sections
//! (`intersect_sections`) and its area is the closed-form Green integral of
//! the exact trimmed pieces. The volume is the integral of that area over
//! the aligned face width by adaptive Gauss-Kronrod quadrature, split where
//! a chevron reverses its hand; its estimate is the summed Kronrod-Gauss
//! difference.
//!
//! Every evaluated section also bounds the true volume. Turning a section
//! by `dφ` about its centre sweeps each boundary element over the area
//! `|dφ| |(p - c) · dp|`, and the symmetric difference of the turned and
//! original section lies in that sweep. Only boundary edges that can enter
//! the other gear's tip disk matter, so the common area is Lipschitz in `z`
//! with `L = Σ |twist| ∫ |d(|p - c|² / 2)|` over those edges. Each node
//! value therefore encloses the integral over the cell around it within
//! `L (w₋² + w₊²) / 2`, and the sum of cells is a certified enclosure of
//! the exact volume that needs no derivative estimate. Intervals are halved
//! until the quadrature estimate meets the volume tolerance and the
//! enclosure of a nonempty common material is within one percent of its
//! upper bound.
use super::chevron::placed_at;
use super::overlap::loops;
use super::{Error, Placement, Section, checkpoint, intersect_sections};
use crate::normal_sweep;
use crate::revolution::length;
use crate::spur_gear::{Edge, Point, herringbone::Herringbone, sub};

const MAX_INTERVALS: usize = 4096;
/// Width of the certified enclosure relative to its upper bound.
const ENCLOSURE_FRACTION: f64 = 0.01;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TwistedVolume {
    /// Adaptive Gauss-Kronrod value of the exact section areas.
    pub volume_mm3: f64,
    /// Summed Kronrod-Gauss difference of that quadrature.
    pub error_mm3: f64,
    /// Certified lower and upper bounds of the exact common volume.
    pub enclosure_mm3: [f64; 2],
}

/// The exact common section area at `z_mm`, `None` when a contact or a
/// vertex crossing leaves that section unresolved.
fn area(
    gears: [&Herringbone; 2],
    placements: [Placement; 2],
    z_mm: f64,
    tolerance_mm: f64,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Option<f64>, Error> {
    let placed = placed_at(gears, placements, z_mm, 0.);
    let sections = gears.map(|gear| Section::Spur(gear.profile()));
    Ok(
        intersect_sections(sections, placed, tolerance_mm, interrupt)?.map(|regions| {
            regions
                .iter()
                .map(|region| region.area_mm2(placed))
                .sum::<f64>()
                .max(0.)
        }),
    )
}

/// `∫ |d(|p|² / 2)|` along a local edge, or an upper bound of it.
fn radial_variation(edge: Edge) -> f64 {
    let squared = |p: Point| (p[0] * p[0] + p[1] * p[1]) / 2.;
    match edge {
        // The involute radius is monotone along the flank.
        Edge::Involute(flank) => (squared(flank.point(1.)) - squared(flank.point(0.))).abs(),
        Edge::Planar(
            normal_sweep::Edge::Arc { center, .. } | normal_sweep::Edge::Circle { center, .. },
        ) if center == [0.; 2] => 0.,
        // |(p - c) · dp| ≤ |p - c| |dp|, and |p - c| ≤ |start| + length.
        _ => {
            let reach = edge.length_mm();
            reach * (length(edge.point(0.)) + reach)
        }
    }
}

/// The Lipschitz constant of the common section area in `z`.
fn lipschitz(gears: [&Herringbone; 2], placements: [Placement; 2], width_mm: f64) -> f64 {
    let middle = placed_at(gears, placements, width_mm / 2., 0.);
    (0..2)
        .map(|i| {
            let other = 1 - i;
            let rate = gears[i].twist_rad_per_mm().abs();
            let reach = gears[other].tip_radius_mm();
            let Ok(edges) = loops(Section::Spur(gears[i].profile())) else {
                return f64::INFINITY;
            };
            let variation: f64 = edges
                .iter()
                .flatten()
                .filter(|&&edge| {
                    // Any edge point stays within its length of the start,
                    // and every phase in the width lies within the half
                    // twist of the middle phase.
                    let start = edge.point(0.);
                    let size = edge.length_mm();
                    let swing = (length(start) + size) * rate * width_mm / 2.;
                    let at = middle[i].point(start);
                    length(sub(at, middle[other].center_mm)) - size - swing
                        <= reach + placements[i].uncertainty_mm + placements[other].uncertainty_mm
                })
                .map(|&edge| radial_variation(edge))
                .sum();
            rate * variation
        })
        .sum()
}

struct Interval {
    a: f64,
    b: f64,
    value: f64,
    error: f64,
    bounds: [f64; 2],
}

const X: [f64; 8] = [
    0.991_455_371_120_812_6,
    0.949_107_912_342_758_5,
    0.864_864_423_359_769_1,
    0.741_531_185_599_394_4,
    0.586_087_235_467_691_1,
    0.405_845_151_377_397_2,
    0.207_784_955_007_898_5,
    0.,
];
const K: [f64; 8] = [
    0.022_935_322_010_529_22,
    0.063_092_092_629_978_55,
    0.104_790_010_322_250_18,
    0.140_653_259_715_525_92,
    0.169_004_726_639_267_9,
    0.190_350_578_064_785_4,
    0.204_432_940_075_298_9,
    0.209_482_141_084_727_83,
];
const G: [f64; 4] = [
    0.129_484_966_168_869_7,
    0.279_705_391_489_276_7,
    0.381_830_050_505_118_9,
    0.417_959_183_673_469_4,
];

/// Seven-point Gauss and fifteen-point Kronrod rules on `[a, b]` and the
/// Lipschitz enclosure of the cells around their fifteen nodes. `None`
/// when a node section is unresolved.
fn kronrod(
    f: &mut dyn FnMut(f64) -> Result<Option<f64>, Error>,
    a: f64,
    b: f64,
    lipschitz: f64,
    ceiling: f64,
) -> Result<Option<Interval>, Error> {
    let (mid, half) = ((a + b) / 2., (b - a) / 2.);
    // Nodes in increasing order with their Kronrod and Gauss weights.
    let mut nodes = [(0., 0., 0., 0.); 15];
    for i in 0..8 {
        // Odd Kronrod nodes, and the centre, are the Gauss nodes.
        let gauss = if i % 2 == 1 { G[i / 2] } else { 0. };
        nodes[i] = (mid - half * X[i], K[i], gauss, 0.);
        nodes[14 - i] = (mid + half * X[i], K[i], gauss, 0.);
    }
    for node in &mut nodes {
        let Some(value) = f(node.0)? else {
            return Ok(None);
        };
        node.3 = value;
    }
    let kronrod: f64 = nodes.iter().map(|n| n.1 * n.3).sum();
    let gauss: f64 = nodes.iter().map(|n| n.2 * n.3).sum();
    let mut bounds = [0., 0.];
    for (k, &(x, _, _, value)) in nodes.iter().enumerate() {
        let left = if k == 0 { a } else { (nodes[k - 1].0 + x) / 2. };
        let right = if k == 14 {
            b
        } else {
            (x + nodes[k + 1].0) / 2.
        };
        let slack = lipschitz * ((x - left).powi(2) + (right - x).powi(2)) / 2.;
        let width = right - left;
        bounds[0] += (value * width - slack).max(0.);
        bounds[1] += (value * width + slack).min(ceiling * width);
    }
    Ok(Some(Interval {
        a,
        b,
        value: kronrod * half,
        error: ((kronrod - gauss) * half).abs(),
        bounds,
    }))
}

/// The common volume of two twisted gears whose aligned sections start on
/// one transverse plane, in a common orthonormal chart. Gaps and hubs are
/// `Invalid`; `Ok(None)` leaves an unresolved section to other providers.
pub fn twisted_interference_volume(
    gears: [&Herringbone; 2],
    placements: [Placement; 2],
    tolerance_mm: f64,
    volume_tolerance_mm3: f64,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Option<TwistedVolume>, Error> {
    let profiles = gears.map(Herringbone::profile);
    let width = profiles[0].height_mm();
    if width != profiles[1].height_mm()
        || gears
            .iter()
            .any(|g| g.is_gapped() || g.profile().hub().is_some())
        || volume_tolerance_mm3.is_nan()
        || volume_tolerance_mm3 <= 0.
        || placements.into_iter().any(|p| !p.valid(tolerance_mm))
        || !(1e-7..=1e-3).contains(&tolerance_mm)
    {
        return Err(Error::Invalid);
    }
    let lipschitz = lipschitz(gears, placements, width);
    let ceiling = profiles[0].area_mm2().min(profiles[1].area_mm2());
    // A chevron reverses its hand at mid-width, where the area has a kink.
    let mut splits = vec![0., width];
    if gears.iter().any(|g| !g.is_helical()) {
        splits.insert(1, width / 2.);
    }
    let mut f = |z: f64| {
        checkpoint(interrupt)?;
        area(gears, placements, z, tolerance_mm, interrupt)
    };
    let mut intervals = Vec::new();
    let mut pending: Vec<[f64; 2]> = splits.windows(2).map(|w| [w[0], w[1]]).collect();
    loop {
        while let Some([a, b]) = pending.pop() {
            if intervals.len() + pending.len() >= MAX_INTERVALS {
                return Err(Error::Limit);
            }
            match kronrod(&mut f, a, b, lipschitz, ceiling)? {
                Some(interval) => intervals.push(interval),
                // A node on a contact event: other nodes avoid it.
                None if b - a > tolerance_mm => {
                    let mid = (a + b) / 2.;
                    pending.extend([[a, mid], [mid, b]]);
                }
                None => return Ok(None),
            }
        }
        let error: f64 = intervals.iter().map(|i| i.error).sum();
        let low: f64 = intervals.iter().map(|i| i.bounds[0]).sum();
        let high: f64 = intervals.iter().map(|i| i.bounds[1]).sum();
        let volume: f64 = intervals.iter().map(|i| i.value).sum();
        // An empty common material has no relative enclosure to reach.
        let loose = volume > volume_tolerance_mm3 && high - low > ENCLOSURE_FRACTION * high;
        if error <= volume_tolerance_mm3 && !loose {
            return Ok(Some(TwistedVolume {
                volume_mm3: volume.max(0.),
                error_mm3: error,
                enclosure_mm3: [low, high],
            }));
        }
        // Refine the quadrature estimate first, then the enclosure.
        let key = |i: &Interval| {
            if error > volume_tolerance_mm3 {
                i.error
            } else {
                i.bounds[1] - i.bounds[0]
            }
        };
        let (worst, _) = intervals
            .iter()
            .enumerate()
            .max_by(|a, b| key(a.1).total_cmp(&key(b.1)))
            .expect("nonempty");
        let Interval { a, b, .. } = intervals.swap_remove(worst);
        let mid = (a + b) / 2.;
        pending.extend([[a, mid], [mid, b]]);
    }
}

#[cfg(test)]
mod tests;
