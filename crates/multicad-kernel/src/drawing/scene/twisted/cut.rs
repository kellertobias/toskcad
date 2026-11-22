//! Exact intersection of a twisted gear with a plane facing the viewer.
use super::super::paths::{Lift, ProfilePath, Segment, add, scaled};
use super::super::*;
use super::flank::Shape;
use super::flank::rotate;
use super::solve::{Search, isolate};
use super::trace::trace;
use super::{EndFace, Gear, Level, Patch, Region, Traced};
use multicad_precision::PortableMath as _;
use multicad_rust_kernel::spur_gear::Edge;
use std::f64::consts::TAU;

/// Whether a path runs with the material on its left as the viewer sees it,
/// given the solid's outward normal at the path's middle.
fn leftward(path: &dyn Path, normal: P3) -> bool {
    let [p, q] = [path.point(0.5 - 1e-4), path.point(0.5 + 1e-4)];
    (q[1] - p[1]) * normal[0] - (q[0] - p[0]) * normal[1] > 0.
}

fn closed(edge: Edge) -> bool {
    matches!(edge, Edge::Planar(e) if e.start() == e.end())
}

/// A transverse section: the profile itself, turned by the local phase.
fn transverse(gear: &Gear, offset: f64) -> Vec<Loop> {
    let rise = gear.toward[2];
    let level = offset / rise;
    // Material immediately behind the plane decides at a shared level.
    let behind = level - rise.signum() * 2. * gear.tolerance;
    let mut rings: Vec<Vec<Edge>> = Vec::new();
    let mut phase = 0.;
    match gear.material(behind, -rise.signum()) {
        Some(Region::Profile(_)) => {
            let ends = [gear.bands[0].z[0], gear.bands[gear.bands.len() - 1].z[1]];
            phase = gear.phase_at(level.clamp(ends[0], ends[1]));
            rings.push(gear.flanks.iter().map(|flank| flank.edge).collect());
        }
        Some(Region::Disk(radius)) => rings.push(vec![Gear::circle(radius, false)]),
        None => return Vec::new(),
    }
    if gear.bore > 0. {
        rings.push(vec![Gear::circle(gear.bore, true)]);
    }
    let lift = gear.frame.level(level, phase);
    let mut output: Loop = Vec::new();
    for edges in rings {
        let mut ring: Loop = Vec::new();
        for edge in edges {
            let stops: &[f64] = if closed(edge) {
                &[0., 0.5, 1.]
            } else {
                &[0., 1.]
            };
            for pair in stops.windows(2) {
                ring.push(Arc::new(ProfilePath {
                    edge,
                    range: if rise < 0. {
                        [pair[1], pair[0]]
                    } else {
                        [pair[0], pair[1]]
                    },
                    lift,
                    straight: false,
                }));
            }
        }
        if rise < 0. {
            ring.reverse();
        }
        output.extend(ring);
    }
    vec![output]
}

/// The plane's curve on a coaxial flank (a tip or root arc swept through a
/// band), in closed form: generators of the cylinder for a plane along the
/// axis, otherwise arcs of the ellipse where the plane cuts the cylinder.
/// Each is restricted to where the twisted arc covers the cylinder: the
/// profile angle `theta - phase(z)` within the arc and `z` within the band.
/// A plane tangent to the cylinder touches it along a generator, which
/// bounds no area and contributes nothing.
fn coaxial(
    gear: &Gear,
    patch: &Patch,
    offset: f64,
    radius: f64,
    output: &mut Loop,
) -> Result<(), GeometryError> {
    let Shape::Arc { start, sweep, .. } = patch.flank.shape else {
        return Ok(());
    };
    let band = patch.band;
    let [rise, slope, azimuth] = [gear.toward[2], gear.slope, gear.azimuth];
    // The arc's profile angles as an ascending span.
    let low = start.min(start + sweep);
    let width = sweep.abs();
    let covered = |profile: f64| (profile - low).rem_euclid(TAU) <= width;
    // Outward normal: radial, toward the side the profile normal points.
    let probe = patch.flank.point(0.5);
    let tangent = patch.flank.tangent(0.5);
    let sign = (tangent[1] * probe[0] - tangent[0] * probe[1]).signum();
    let normal_at = |theta: f64| {
        let (sin, cos) = theta.portable_sin_cos();
        gear.frame.vector([sign * cos, sign * sin, 0.])
    };
    if gear.edge_on() {
        let ratio = offset / (slope * radius);
        // Within rounding of tangent the plane touches the cylinder along a
        // generator, which bounds no area.
        if ratio.abs() >= 1. - 1e-12 {
            return Ok(());
        }
        for side in [1., -1.] {
            let theta = azimuth + side * ratio.portable_acos();
            // phase(z) = theta - profile is linear in z; the covered
            // profile angles give z intervals, at most two per turn.
            let profile = |z: f64| theta - band.phase_at(z);
            let ends = band.z.map(profile);
            let [lo, hi] = [ends[0].min(ends[1]), ends[0].max(ends[1])];
            let mut k = ((lo - low - width) / TAU).ceil();
            while low + k * TAU <= hi {
                let [a, b] = [(low + k * TAU).max(lo), (low + k * TAU + width).min(hi)];
                if b > a {
                    let z =
                        [a, b].map(|angle| band.z[0] + (theta - angle - band.phase) / band.twist);
                    let at = |z: f64| {
                        let (sin, cos) = theta.portable_sin_cos();
                        gear.frame.place([radius * cos, radius * sin, z])
                    };
                    let mut path = Segment(at(z[0]), at(z[1]));
                    if !leftward(&path, normal_at(theta)) {
                        path = Segment(at(z[1]), at(z[0]));
                    }
                    output.push(Arc::new(path));
                }
                k += 1.;
            }
        }
        return Ok(());
    }
    // On the ellipse, z(theta) = (offset - slope radius cos(theta - azimuth)) / rise.
    let height = |theta: f64| (offset - slope * radius * (theta - azimuth).portable_cos()) / rise;
    // Angles at which the ellipse crosses the band's two levels.
    let mut cuts: Vec<f64> = Vec::new();
    for z in band.z {
        let cosine = (offset - rise * z) / (slope * radius);
        if cosine.abs() < 1. {
            let angle = cosine.portable_acos();
            cuts.extend([azimuth + angle, azimuth - angle]);
        }
    }
    // Angles at which the profile angle reaches an end of the arc:
    // theta - phase(z(theta)) - end is a multiple of a full turn.
    let lead = band.twist * slope * radius / rise;
    let mut budget = 1_000_000_usize;
    for end in [low, low + width] {
        let mut roots = Vec::new();
        // f(theta) = theta - phase(z(theta)) - end, wrapped to half a turn.
        isolate(
            &mut |theta| {
                let value = theta - band.phase_at(height(theta)) - end;
                (value + std::f64::consts::PI).rem_euclid(TAU) - std::f64::consts::PI
            },
            &Search {
                slope: 1. + lead.abs(),
                leaf: 1e-12,
                noise: 0.,
                second: lead.abs(),
                lift: 0.,
                wrapped: true,
            },
            [
                azimuth - std::f64::consts::PI,
                azimuth + std::f64::consts::PI,
            ],
            &mut budget,
            &mut roots,
        )?;
        cuts.extend(roots);
    }
    let base = azimuth - std::f64::consts::PI;
    let mut cuts: Vec<f64> = cuts
        .into_iter()
        .map(|angle| base + (angle - base).rem_euclid(TAU))
        .collect();
    cuts.extend([base, base + TAU]);
    cuts.sort_by(f64::total_cmp);
    cuts.dedup_by(|a, b| (*a - *b).abs() <= 1e-13);
    let lift = Lift {
        origin: gear.frame.place([0., 0., offset / rise]),
        x: add(
            gear.frame.axes[0],
            scaled(gear.frame.axes[2], -gear.toward[0] / rise),
        ),
        y: add(
            gear.frame.axes[1],
            scaled(gear.frame.axes[2], -gear.toward[1] / rise),
        ),
    };
    // Join consecutive covered pieces into maximal arcs.
    let inside = |theta: f64| {
        let z = height(theta);
        z > band.z[0] && z < band.z[1] && covered(theta - band.phase_at(z))
    };
    let mut runs: Vec<[f64; 2]> = Vec::new();
    for pair in cuts.windows(2) {
        if pair[1] - pair[0] <= 1e-13 || !inside(0.5 * (pair[0] + pair[1])) {
            continue;
        }
        match runs.last_mut() {
            Some(last) if last[1] == pair[0] => last[1] = pair[1],
            _ => runs.push([pair[0], pair[1]]),
        }
    }
    // A run through the seam of the turn continues the first run.
    if runs.len() > 1 && runs[0][0] == base && runs[runs.len() - 1][1] == base + TAU {
        let first = runs.remove(0);
        let last = runs.len() - 1;
        runs[last][1] = first[1] + TAU;
    }
    for [a, b] in runs {
        let mut path = ProfilePath {
            edge: Gear::circle(radius, false),
            range: [a / TAU, b / TAU],
            lift,
            straight: false,
        };
        if !leftward(&path, normal_at(0.5 * (a + b))) {
            path.range = [path.range[1], path.range[0]];
        }
        output.push(Arc::new(path));
    }
    Ok(())
}

pub(super) fn section(gear: &Gear, depth: f64) -> Result<Vec<Loop>, GeometryError> {
    let offset = depth - gear.frame.origin[2];
    if gear.axial() {
        return Ok(transverse(gear, offset));
    }
    let rise = gear.toward[2];
    let slope = gear.slope;
    let along = [gear.toward[0] / slope, gear.toward[1] / slope];
    let across = [-along[1], along[0]];
    let mut output: Loop = Vec::new();

    // Flanks: the curve of each helicoid on the plane.
    let level = Level::Cut(offset);
    let mut budget = 50_000_000_usize;
    for band in 0..gear.bands.len() {
        for flank in 0..gear.flanks.len() {
            let patch = gear.patch(band, flank);
            // Skip flanks whose support along the sight misses the plane.
            let reach = slope * patch.flank.radius[1];
            let levels = patch.band.z.map(|z| rise * z);
            if offset > levels[0].max(levels[1]) + reach
                || offset < levels[0].min(levels[1]) - reach
            {
                continue;
            }
            if let Some(radius) = patch.flank.coaxial() {
                coaxial(gear, &patch, offset, radius, &mut output)?;
                continue;
            }
            let chains = trace(
                &|s, zeta| patch.field(level, s, zeta),
                &patch.bounds(level),
                &mut budget,
            )?;
            for chain in chains {
                if chain.bridge {
                    let [a, b] = chain.ends.map(|[s, zeta]| patch.surface(s, zeta));
                    let [s, zeta] = [0, 1].map(|i| 0.5 * (chain.ends[0][i] + chain.ends[1][i]));
                    let mut path = Segment(a, b);
                    if !leftward(&path, patch.normal(s, zeta)) {
                        path = Segment(b, a);
                    }
                    output.push(Arc::new(path));
                    continue;
                }
                let mut path = Traced {
                    patch,
                    level,
                    steps: Arc::new(chain.steps),
                    range: [0., 1.],
                };
                let [s, zeta] = path.parameters(0.5);
                if !leftward(&path, patch.normal(s, zeta)) {
                    path.range = [1., 0.];
                }
                output.push(Arc::new(path));
            }
        }
    }

    // Bore and hub walls: generators or arcs of an ellipse.
    let mut cylinders: Vec<(f64, [f64; 2], f64)> = Vec::new();
    if gear.bore > 0. {
        cylinders.push((gear.bore, gear.axis, -1.));
    }
    for cylinder in &gear.cylinders {
        cylinders.push((cylinder.radius, cylinder.z, 1.));
    }
    for (radius, z, outward) in cylinders {
        let normal_at = |angle: f64| {
            let (sin, cos) = angle.portable_sin_cos();
            gear.frame.vector([outward * cos, outward * sin, 0.])
        };
        if gear.edge_on() {
            let line = offset / slope;
            if radius * (1. - 1e-12) <= line.abs() {
                continue;
            }
            for side in [-1., 1.] {
                let reach = side * (radius * radius - line * line).sqrt();
                let [x, y] = [0, 1].map(|i| line * along[i] + reach * across[i]);
                let [a, b] = z.map(|z| gear.frame.place([x, y, z]));
                let mut path = Segment(a, b);
                if !leftward(&path, normal_at(y.portable_atan2(x))) {
                    path = Segment(b, a);
                }
                output.push(Arc::new(path));
            }
            continue;
        }
        // cos(angle - azimuth) at the two ends of the wall.
        let ends = z.map(|z| (offset - rise * z) / (slope * radius));
        let [lo, hi] = [ends[0].min(ends[1]).max(-1.), ends[0].max(ends[1]).min(1.)];
        if hi - lo <= 1e-12 {
            continue;
        }
        let lift = Lift {
            origin: gear.frame.place([0., 0., offset / rise]),
            x: add(
                gear.frame.axes[0],
                scaled(gear.frame.axes[2], -gear.toward[0] / rise),
            ),
            y: add(
                gear.frame.axes[1],
                scaled(gear.frame.axes[2], -gear.toward[1] / rise),
            ),
        };
        for side in [1., -1.] {
            let angles = [hi.portable_acos(), lo.portable_acos()].map(|a| gear.azimuth + side * a);
            let mut path = ProfilePath {
                edge: Gear::circle(radius, false),
                range: angles.map(|a| a / TAU),
                lift,
                straight: false,
            };
            if !leftward(&path, normal_at(0.5 * (angles[0] + angles[1]))) {
                path.range = [path.range[1], path.range[0]];
            }
            output.push(Arc::new(path));
        }
    }

    // End faces: chords of the plane's line in each transverse face.
    for EndFace {
        z,
        outward,
        phase,
        inner,
        outer,
    } in gear.end_faces()
    {
        let line = (offset - rise * z) / slope;
        let base = [line * along[0], line * along[1], z];
        let direction = [across[0], across[1], 0.];
        let inside = |p: P3| {
            let radius = p[0].portable_hypot(p[1]);
            radius > inner
                && match phase {
                    Some(phase) => gear.contains_profile(rotate([p[0], p[1]], -phase)),
                    None => radius < outer,
                }
        };
        let normal = gear.frame.vector([0., 0., outward]);
        for span in gear.spans(base, direction, &inside)? {
            let [a, b] = span.map(|t| {
                gear.frame
                    .place(std::array::from_fn(|i| base[i] + t * direction[i]))
            });
            let mut path = Segment(a, b);
            if !leftward(&path, normal) {
                path = Segment(b, a);
            }
            output.push(Arc::new(path));
        }
    }
    Ok(if output.is_empty() {
        Vec::new()
    } else {
        vec![output]
    })
}
