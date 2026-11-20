//! Certified minimum distance between two bevels sharing one pitch apex.
//!
//! Each solid is the set of points `r u` whose unit direction `u` lies in
//! its spherical domain and whose distance `r` lies between the toe and
//! heel back cones, `[T(γ), H(γ)] = [d_toe, d_heel] / cos(γ - δ)`, less the
//! axial bore. Hence the distance between the solids is the minimum, over
//! direction pairs, of the distance between two radial segments. For an
//! angle `θ ≥ θ₀` between the directions and radii in boxes `R₁ × R₂`, the
//! squared distance `r₁² + r₂² - 2 r₁ r₂ cos θ` is convex in the radii and
//! is smallest on a lower edge of the box, which bounds the pair below.
//!
//! The domain of one sector is covered by cells. A tooth cell is a range of
//! the right flank chain parameter `s` (fillet, meridian, involute) times a
//! range of the azimuth fraction `t ∈ [-1, 1]`: the direction
//! `polar(γ(s), t φ(s))`. Moving `t` at fixed `s` travels at `sin γ |φ|`,
//! and moving `s` at fixed `t` travels no faster than the flank itself, at
//! most `|u'(s_c)| + |s - s_c| · max|u''|` (`curvature_bound`). A core cell
//! is a polar box below the root cone. These speeds bound the angular radius
//! of every cell about its center direction, and the polar angle changes no
//! more than that radius.
//!
//! Every cell center is a material direction. Its optimal radii give a pair
//! of actual material points, outside the bore, whose distance bounds the
//! minimum from above. A depth-first search prunes a cell pair once its
//! lower bound is within the requested gap of the best actual pair.
use super::clearance::Placement;
use super::*;

/// Cell-pair evaluations before the search reports a limit.
const MAX_EVALUATIONS: usize = 40_000_000;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Error {
    Invalid,
    Limit,
    Cancelled,
    /// A cell pair could not be decided at the smallest representable size.
    Unresolved,
}

/// Actual points on both solids and a proven lower bound of their distance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Found {
    /// Distance between `points_mm`, an upper bound of the minimum.
    pub distance_mm: f64,
    pub lower_mm: f64,
    pub points_mm: [Point; 2],
}

#[derive(Clone, Copy, Debug)]
enum Kind {
    Tooth(Piece),
    Core,
}

/// Parameter box of one sector: `(s, t)` for teeth, `(γ, φ)` for the core.
#[derive(Clone, Copy, Debug)]
struct Cell {
    sector: u32,
    kind: Kind,
    ranges: [[f64; 2]; 2],
}

#[derive(Clone, Copy, Debug)]
struct Bound {
    /// Sector-frame center direction and its polar angle.
    local: Point,
    gamma: f64,
    /// Placed unit center direction.
    world: Point,
    /// Angular radius of the cell about `world`, and the share of each
    /// parameter range in it.
    radius: f64,
    share: [f64; 2],
}

struct Side<'a> {
    gear: &'a Bevel,
    placement: Placement,
    /// Bound of the placement's departure from an orthogonal map.
    defect: f64,
}

fn middle([a, b]: [f64; 2]) -> f64 {
    (a + b) / 2.
}

fn angle(a: Point, b: Point) -> f64 {
    norm(cross(a, b)).atan2(dot(a, b))
}

/// Smallest `r₁² + r₂² - 2 r₁ r₂ c` over a positive radius box, and its
/// arguments. The minimizer lies on `r₁ = a₁` or on `r₂ = a₂`.
fn radial(r1: [f64; 2], r2: [f64; 2], c: f64) -> (f64, [f64; 2]) {
    let value = |p: f64, q: f64| (p * p + q * q - 2. * p * q * c).max(0.);
    let first = [r1[0], (r1[0] * c).clamp(r2[0], r2[1])];
    let second = [(r2[0] * c).clamp(r1[0], r1[1]), r2[0]];
    let [a, b] = [first, second].map(|[p, q]| value(p, q));
    if a <= b { (a, first) } else { (b, second) }
}

impl<'a> Side<'a> {
    fn new(gear: &'a Bevel, placement: Placement) -> Self {
        let mut defect = 0_f64;
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1. } else { 0. };
                defect = defect.max((dot(placement.axes[i], placement.axes[j]) - expected).abs());
            }
        }
        Self {
            gear,
            placement,
            // A 3×3 Gram defect `e` bounds the spectral defect by `3e`; the
            // placed angles and lengths then move by at most twice that.
            defect: 6. * defect,
        }
    }

    fn cells(&self) -> Vec<Cell> {
        let gear = self.gear;
        let spacing = PI / f64::from(gear.spec.teeth);
        let mut cells = vec![];
        for sector in 0..gear.spec.teeth {
            cells.push(Cell {
                sector,
                kind: Kind::Core,
                ranges: [[0., gear.root_gamma()], [-spacing, spacing]],
            });
            for &piece in gear.pieces() {
                if matches!(
                    piece,
                    Piece::FilletRight | Piece::MeridianRight | Piece::InvoluteRight
                ) {
                    cells.push(Cell {
                        sector,
                        kind: Kind::Tooth(piece),
                        ranges: [[0., 1.], [-1., 1.]],
                    });
                }
            }
        }
        cells
    }

    fn bound(&self, cell: &Cell) -> Bound {
        let gear = self.gear;
        let half = cell.ranges.map(|[a, b]| (b - a) / 2.);
        let (local, gamma, share) = match cell.kind {
            Kind::Tooth(piece) => {
                let sample = gear.sample(piece, middle(cell.ranges[0]));
                let t = middle(cell.ranges[1]);
                let speed = norm(sample.derivative) + half[0] * gear.curvature_bound(piece);
                (
                    polar(sample.gamma, t * sample.phi),
                    sample.gamma,
                    [
                        half[0] * speed,
                        half[1] * sample.gamma.sin() * sample.phi.abs(),
                    ],
                )
            }
            Kind::Core => {
                let gamma = middle(cell.ranges[0]);
                (
                    polar(gamma, middle(cell.ranges[1])),
                    gamma,
                    [half[0], half[1] * gamma.sin()],
                )
            }
        };
        let world = self.placement.world(gear.rotate(cell.sector, local));
        let length = norm(world);
        // Rounding of the trigonometric evaluation and the placed map.
        let slack = self.defect + 64. * f64::EPSILON;
        Bound {
            local,
            gamma,
            world: scale(world, 1. / length),
            radius: (share[0] + share[1]) * (1. + 1e-12) + slack,
            share,
        }
    }

    /// Back-cone radii over a polar range, widened for the placed map.
    fn radii(&self, gamma: [f64; 2]) -> [f64; 2] {
        let gear = self.gear;
        let pitch = gear.spec.pitch_cone_angle_rad;
        let scale = |g: f64| 1. / (g - pitch).cos();
        let low = scale(pitch.clamp(gamma[0], gamma[1]));
        let high = scale(gamma[0]).max(scale(gamma[1]));
        [
            gear.end_distance_mm(End::Toe) * low * (1. - self.defect - 1e-12),
            gear.end_distance_mm(End::Heel) * high * (1. + self.defect + 1e-12),
        ]
    }

    /// Material radii along one direction of polar angle `gamma`.
    fn material(&self, gamma: f64) -> Option<[f64; 2]> {
        let gear = self.gear;
        let mut low = gear.end_distance_along(End::Toe, gamma);
        let high = gear.end_distance_along(End::Heel, gamma);
        if let Some(bore) = gear.bore_radius_mm() {
            // Strictly outside the bore, so the point is retained material.
            low = low.max(bore / gamma.sin() * (1. + 1e-12));
        }
        (low.is_finite() && low <= high).then_some([low, high])
    }

    fn point(&self, cell: &Cell, bound: &Bound, r: f64) -> Point {
        let p = self
            .placement
            .world(self.gear.rotate(cell.sector, scale(bound.local, r)));
        std::array::from_fn(|i| self.placement.apex_mm[i] + p[i])
    }
}

/// Minimum distance between two bevels whose placed pitch apexes coincide
/// within `tolerance_mm` and whose spherical domains are disjoint, so the
/// solids are. The result's lower bound is within `gap_mm` of its distance.
pub fn minimum(
    gears: [&Bevel; 2],
    placements: [Placement; 2],
    tolerance_mm: f64,
    gap_mm: f64,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Found, Error> {
    if !placements.iter().all(Placement::valid)
        || tolerance_mm.is_nan()
        || tolerance_mm <= 0.
        || gap_mm.is_nan()
        || gap_mm <= 0.
        || norm(sub(placements[0].apex_mm, placements[1].apex_mm)) > tolerance_mm
    {
        return Err(Error::Invalid);
    }
    let sides = [0, 1].map(|i| Side::new(gears[i], placements[i]));
    // The apexes differ by `apart`; measuring both from their own apex moves
    // every pair distance by at most that much.
    let apart = norm(sub(placements[0].apex_mm, placements[1].apex_mm));
    let mut best: Option<(f64, [Point; 2])> = None;
    let mut lower = f64::INFINITY;
    let mut evaluations = 0_usize;
    let mut evaluate = |pair: [Cell; 2],
                        best: &mut Option<(f64, [Point; 2])>|
     -> Result<(f64, [Bound; 2]), Error> {
        evaluations += 1;
        if evaluations > MAX_EVALUATIONS {
            return Err(Error::Limit);
        }
        if evaluations.is_multiple_of(4096) && interrupt.is_some_and(|stop| stop()) {
            return Err(Error::Cancelled);
        }
        let bounds = [0, 1].map(|i| sides[i].bound(&pair[i]));
        let theta = angle(bounds[0].world, bounds[1].world);
        // Lower bound over the cells.
        let spread = bounds[0].radius + bounds[1].radius;
        let floor = (theta - spread).max(0.);
        let ranges = [0, 1].map(|i| {
            let b = &bounds[i];
            sides[i].radii([(b.gamma - b.radius).max(0.), (b.gamma + b.radius).min(PI)])
        });
        let (square, _) = radial(ranges[0], ranges[1], floor.cos());
        let bound = (square.sqrt() * (1. - 1e-12) - apart - 1e-12).max(0.);
        // Actual material pair at the centers.
        if let (Some(a), Some(b)) = (
            sides[0].material(bounds[0].gamma),
            sides[1].material(bounds[1].gamma),
        ) {
            let (_, [r1, r2]) = radial(a, b, theta.cos());
            let points = [
                sides[0].point(&pair[0], &bounds[0], r1),
                sides[1].point(&pair[1], &bounds[1], r2),
            ];
            let distance = norm(sub(points[0], points[1]));
            if distance.is_finite() && best.is_none_or(|(d, _)| distance < d) {
                *best = Some((distance, points));
            }
        }
        Ok((bound, bounds))
    };
    let initial: Vec<[Cell; 2]> = {
        let cells = sides.each_ref().map(Side::cells);
        cells[0]
            .iter()
            .flat_map(|a| cells[1].iter().map(move |b| [*a, *b]))
            .collect()
    };
    let mut stack = Vec::with_capacity(initial.len());
    for pair in initial {
        let (bound, bounds) = evaluate(pair, &mut best)?;
        stack.push((bound, pair, bounds));
    }
    // Visit the most promising pairs last pushed, so first popped.
    stack.sort_by(|a, b| b.0.total_cmp(&a.0));
    while let Some((bound, pair, bounds)) = stack.pop() {
        let (distance, _) = best.ok_or(Error::Unresolved)?;
        if bound >= distance - gap_mm {
            lower = lower.min(bound);
            continue;
        }
        // Split the cell with the larger contribution along its larger share.
        let weight = |i: usize| bounds[i].radius * sides[i].radii([bounds[i].gamma; 2])[1];
        let which = usize::from(weight(1) > weight(0));
        let cell = pair[which];
        let axis = usize::from(bounds[which].share[1] > bounds[which].share[0]);
        let [lo, hi] = cell.ranges[axis];
        let mid = (lo + hi) / 2.;
        if mid <= lo || mid >= hi {
            return Err(Error::Unresolved);
        }
        let mut children = vec![];
        for range in [[lo, mid], [mid, hi]] {
            let mut child = cell;
            child.ranges[axis] = range;
            let mut next = pair;
            next[which] = child;
            let (b, bs) = evaluate(next, &mut best)?;
            children.push((b, next, bs));
        }
        children.sort_by(|a, b| b.0.total_cmp(&a.0));
        stack.extend(children);
    }
    let (distance, points) = best.ok_or(Error::Unresolved)?;
    Ok(Found {
        distance_mm: distance,
        lower_mm: lower.min(distance),
        points_mm: points,
    })
}

#[cfg(test)]
mod tests {
    use super::super::clearance::{ResultKind, certify};
    use super::super::tests::spec;
    use super::*;

    fn pair(counts: [u32; 2], drive: f64) -> ([Bevel; 2], [Placement; 2]) {
        let gears = [0, 1].map(|i| {
            Bevel::new(spec(1., counts[i], counts[1 - i], 2., 1., 0.2), 1e-3, 1e7).unwrap()
        });
        let [za, zb] = counts.map(f64::from);
        let phase = PI - PI / zb - drive * za / zb;
        let (s0, c0) = drive.sin_cos();
        let (s1, c1) = phase.sin_cos();
        let placements = [
            Placement {
                apex_mm: [0.; 3],
                axes: [[c0, s0, 0.], [-s0, c0, 0.], [0., 0., 1.]],
            },
            Placement {
                apex_mm: [0.; 3],
                axes: [[0., s1, -c1], [0., c1, s1], [1., 0., 0.]],
            },
        ];
        (gears, placements)
    }

    #[test]
    fn meshing_pairs_measure_a_certified_gap_below_the_flank_witness() {
        for (counts, drive) in [([24, 32], 0.), ([24, 32], 0.05), ([32, 24], -0.11)] {
            let (gears, placements) = pair(counts, drive);
            let refs = [&gears[0], &gears[1]];
            let clear = certify(refs, placements, 1e-6, 0.2, None).unwrap();
            assert_eq!(clear.result, ResultKind::Clear);
            let found = minimum(refs, placements, 1e-6, 1e-5, None).unwrap();
            assert!(found.lower_mm <= found.distance_mm);
            assert!(found.distance_mm - found.lower_mm <= 1e-5);
            assert!(found.lower_mm > 0.);
            let measured = norm(sub(found.points_mm[0], found.points_mm[1]));
            assert_eq!(measured, found.distance_mm);
            // Both witnesses are material points of their own gear.
            for (i, point) in found.points_mm.into_iter().enumerate() {
                let local = placements[i].axes.map(|axis| dot(axis, point));
                assert_ne!(gears[i].classify(local, 1e-9), Containment::Outside);
            }
            // The flank witness of the clearance check is an actual pair, so
            // the certified minimum lies at or below it.
            let witness = clear.nearest.unwrap();
            assert!(found.lower_mm <= witness.distance_mm + 1e-12);
            assert!(found.distance_mm < 0.1, "{}", found.distance_mm);
        }
        let (gears, placements) = pair([24, 32], 0.);
        let refs = [&gears[0], &gears[1]];
        let mut moved = placements;
        moved[1].apex_mm = [1., 0., 0.];
        assert_eq!(minimum(refs, moved, 1e-6, 1e-5, None), Err(Error::Invalid));
        assert_eq!(
            minimum(refs, placements, 1e-6, 0., None),
            Err(Error::Invalid)
        );
        let stop = || true;
        assert_eq!(
            minimum(refs, placements, 1e-6, 1e-5, Some(&stop)),
            Err(Error::Cancelled)
        );
    }
}
