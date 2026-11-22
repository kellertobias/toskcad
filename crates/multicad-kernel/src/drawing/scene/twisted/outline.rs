//! Edges, smooth seams and outlines of a twisted gear in one view.
use super::super::paths::{ProfilePath, Segment, dot2};
use super::super::profile::folds;
use super::super::*;
use super::flank::rotate;
use super::solve::{Search, isolate};
use super::trace::trace;
use super::{Gear, Helix, Level, Patch, Region, Traced};
use multicad_precision::PortableMath as _;
use multicad_rust_kernel::{normal_sweep, spur_gear::Edge};
use std::f64::consts::{FRAC_PI_2, TAU};

/// Where outlines end on the rims and rails of the flanks, so that every
/// edge is split where the adjacent flank turns away from the viewer.
struct Marks {
    /// By band, band end and profile edge: edge fractions.
    rims: Vec<[Vec<Vec<f64>>; 2]>,
    /// By band and profile vertex: band fractions.
    rails: Vec<Vec<Vec<f64>>>,
}

impl Marks {
    fn add(&mut self, band: usize, flank: usize, [s, zeta]: [f64; 2]) {
        let count = self.rims[0][0].len();
        const EDGE: f64 = 1e-7;
        if zeta <= EDGE {
            self.rims[band][0][flank].push(s);
        } else if zeta >= 1. - EDGE {
            self.rims[band][1][flank].push(s);
        }
        if s <= EDGE {
            self.rails[band][flank].push(zeta);
        } else if s >= 1. - EDGE {
            self.rails[band][(flank + 1) % count].push(zeta);
        }
    }
}

fn stops(mut marks: Vec<f64>) -> Vec<f64> {
    marks.retain(|m| *m > 1e-9 && *m < 1. - 1e-9);
    marks.sort_by(f64::total_cmp);
    marks.dedup_by(|a, b| (*a - *b).abs() <= 1e-9);
    let mut output = vec![0.];
    output.extend(marks);
    output.push(1.);
    output
}

/// Straight outlines of a flank whose normal turns uniformly along the edge
/// while `c . tau` stays constant: involute helicoids and coaxial cylinders.
/// With the normal's azimuth `theta(s) + phase(z)` relative to the sight
/// azimuth called `omega`, the outline is `slope cos(omega) + lead m = 0`.
fn rulings(
    gear: &Gear,
    patch: &Patch,
    [start, turn]: [f64; 2],
    mut emit: impl FnMut([[f64; 2]; 2]),
) {
    let moment = dot2(patch.flank.point(0.5), patch.flank.tangent(0.5));
    let ratio = -patch.band.twist * gear.toward[2] * moment / gear.slope;
    if ratio.abs() >= 1. - 1e-12 {
        // The flank faces one way throughout, or only grazes the sight.
        return;
    }
    let sweep = patch.band.twist * (patch.band.z[1] - patch.band.z[0]);
    let origin = start + patch.band.phase - gear.azimuth;
    let [lo, hi] = [
        origin + turn.min(0.) + sweep.min(0.),
        origin + turn.max(0.) + sweep.max(0.),
    ];
    for omega in [ratio.portable_acos(), -ratio.portable_acos()] {
        let mut k = ((lo - omega) / TAU).ceil();
        while omega + k * TAU <= hi {
            // origin + turn s + sweep zeta = omega + k TAU
            let target = omega + k * TAU - origin;
            let ends = [target / turn, (target - sweep) / turn];
            let [a, b] = [ends[0].min(ends[1]).max(0.), ends[0].max(ends[1]).min(1.)];
            if b - a > 1e-12 {
                emit([a, b].map(|s| [s, ((target - turn * s) / sweep).clamp(0., 1.)]));
            }
            k += 1.;
        }
    }
}

/// Parameters of a traced outline at which its paper image has a cusp: its
/// tangent passes through the sight direction, where its visibility can
/// change. With the curve's surface tangent `T = f_s S_z - f_z S_s` of the
/// fold function `f` and the normal `n`, the paper tangent vanishes where
/// `h = T . (n x d)` does. Every sign change of `h` along each chart is
/// isolated with a bound of its slope from the chart's slope bound and the
/// derivative bounds of `f`, the surface and its normal.
fn cusps(path: &Traced, budget: &mut usize) -> Result<Vec<f64>, GeometryError> {
    let patch = path.patch;
    let level = patch.bounds(Level::Fold);
    let [g1, g2] = [
        level.first[0].max(level.first[1]),
        level.second.iter().copied().fold(0., f64::max),
    ];
    let [s1, s2] = patch.surface_bounds();
    let [c0, c1, c2, _] = patch.flank.size;
    let twist = patch.band.twist.abs();
    let height = patch.band.z[1] - patch.band.z[0];
    // Traced flanks are lines and arcs: |c'| is the constant c1.
    let speed = patch.flank.at(0.5)[1];
    let speed = speed[0].portable_hypot(speed[1]);
    let normal = 1. + twist * c0;
    let turning = (c2 / speed + twist * (c1 + c0 * c2 / speed)).max(height * twist);
    let tangent = 2. * g1 * s1;
    let change = 2. * (g2 * s1 + g1 * s2);
    let gradient = change * normal + tangent * turning;
    let value = |s: f64, zeta: f64| {
        let (_, g) = patch.field(Level::Fold, s, zeta);
        let [ss, sz] = patch.tangents(s, zeta);
        let t: P3 = std::array::from_fn(|i| g[1] * ss[i] - g[0] * sz[i]);
        let n = patch.normal(s, zeta);
        t[0] * n[1] - t[1] * n[0]
    };
    let count = path.steps.len() as f64;
    let mut output = Vec::new();
    for (index, step) in path.steps.iter().enumerate() {
        let eval = |s: f64, zeta: f64| patch.field(Level::Fold, s, zeta);
        let mut roots = Vec::new();
        isolate(
            &mut |fraction| {
                let [s, zeta] = step.point(&eval, fraction);
                value(s, zeta)
            },
            &Search {
                slope: gradient * (1. + step.rate) * (step.u[1] - step.u[0]).abs(),
                leaf: 1e-12,
                noise: 1e-12 * gradient,
                second: f64::INFINITY,
                lift: 0.,
                wrapped: false,
            },
            [0., 1.],
            budget,
            &mut roots,
        )?;
        output.extend(roots.into_iter().map(|f| (index as f64 + f) / count));
    }
    output.retain(|t| *t > 1e-12 && *t < 1. - 1e-12);
    output.sort_by(f64::total_cmp);
    Ok(output)
}

/// The arcs of the circle of a radius about the axis, in profile
/// coordinates, that no coaxial profile arc of the same radius covers.
fn uncovered(gear: &Gear, radius: f64) -> Vec<Edge> {
    let mut covered: Vec<[f64; 2]> = gear
        .flanks
        .iter()
        .filter(|flank| {
            flank
                .coaxial()
                .is_some_and(|own| (own - radius).abs() <= 1e-9 * radius)
        })
        .filter_map(|flank| {
            let super::flank::Shape::Arc { start, sweep, .. } = flank.shape else {
                return None;
            };
            let [a, b] = [start, start + sweep];
            let lo = a.min(b).rem_euclid(TAU);
            Some([lo, lo + sweep.abs()])
        })
        .collect();
    if covered.is_empty() {
        return vec![Gear::circle(radius, false)];
    }
    covered.sort_by(|a, b| a[0].total_cmp(&b[0]));
    let point = |angle: f64| {
        let (sin, cos) = angle.portable_sin_cos();
        [radius * cos, radius * sin]
    };
    let mut output = Vec::new();
    for (index, span) in covered.iter().enumerate() {
        let next = covered
            .get(index + 1)
            .map_or(covered[0][0] + TAU, |next| next[0]);
        if next - span[1] > 1e-12 {
            output.push(Edge::Planar(normal_sweep::Edge::Arc {
                center: [0.; 2],
                a: point(span[1]),
                b: point(next),
                clockwise: false,
            }));
        }
    }
    output
}

pub(super) fn candidates(gear: &Gear) -> Result<Vec<Candidate>, GeometryError> {
    let count = gear.flanks.len();
    let axial = gear.axial();
    let mut output: Vec<Candidate> = Vec::new();
    let mut marks = Marks {
        rims: vec![[vec![Vec::new(); count], vec![Vec::new(); count]]; gear.bands.len()],
        rails: vec![vec![Vec::new(); count]; gear.bands.len()],
    };
    let mut budget = 50_000_000_usize;
    if !axial {
        for band in 0..gear.bands.len() {
            for flank in 0..count {
                let patch = gear.patch(band, flank);
                let ruled = patch.flank.coaxial().is_some()
                    || matches!(patch.flank.shape, super::flank::Shape::Involute { .. });
                if ruled && let Some(turn) = patch.flank.normal_turn() {
                    rulings(gear, &patch, turn, |ends| {
                        for end in ends {
                            marks.add(band, flank, end);
                        }
                        let [a, b] = ends.map(|[s, zeta]| patch.surface(s, zeta));
                        output.push(Candidate {
                            kind: Kind::Outline,
                            path: Arc::new(Segment(a, b)),
                        });
                    });
                    continue;
                }
                // Radial lines and fillets: trace the fold of the flank.
                let chains = trace(
                    &|s, zeta| patch.field(Level::Fold, s, zeta),
                    &patch.bounds(Level::Fold),
                    &mut budget,
                )?;
                for chain in chains {
                    for end in chain.ends {
                        marks.add(band, flank, end);
                    }
                    if chain.bridge {
                        let [a, b] = chain.ends.map(|[s, zeta]| patch.surface(s, zeta));
                        output.push(Candidate {
                            kind: Kind::Outline,
                            path: Arc::new(Segment(a, b)),
                        });
                        continue;
                    }
                    let steps = Arc::new(chain.steps);
                    let whole = Traced {
                        patch,
                        level: Level::Fold,
                        steps: steps.clone(),
                        range: [0., 1.],
                    };
                    let mut cuts = vec![0.];
                    cuts.extend(cusps(&whole, &mut budget)?);
                    cuts.push(1.);
                    for pair in cuts.windows(2) {
                        if pair[1] - pair[0] <= 1e-9 {
                            continue;
                        }
                        let path = Traced {
                            patch,
                            level: Level::Fold,
                            steps: steps.clone(),
                            range: [pair[0], pair[1]],
                        };
                        output.push(Candidate {
                            kind: Kind::Outline,
                            path: Arc::new(path),
                        });
                    }
                }
            }
        }
    }

    // Rims of the two tooth caps, and the apex ring of a herringbone.
    let across = [-gear.toward[1] / gear.slope, gear.toward[0] / gear.slope];
    let ring = |edge: Edge,
                z: f64,
                phase: f64,
                kind: Kind,
                extra: Vec<f64>,
                output: &mut Vec<Candidate>| {
        let mut marks = extra;
        if !axial {
            marks.extend(folds(edge, rotate(across, -phase)));
        }
        if matches!(edge, Edge::Planar(e) if e.start() == e.end()) {
            marks.push(0.5);
        }
        let lift = gear.frame.level(z, phase);
        for pair in stops(marks).windows(2) {
            output.push(Candidate {
                kind,
                path: Arc::new(ProfilePath {
                    edge,
                    range: [pair[0], pair[1]],
                    lift,
                    straight: gear.edge_on(),
                }),
            });
        }
    };
    // Band ends: a ring shared by two bands (the apex) turns sharply except
    // on the coaxial arcs; at a coaxial cylinder of the same radius (the
    // root core of a gap) the root arcs continue smoothly into it.
    let same = |a: f64, b: f64| (a - b).abs() <= 1e-9 * a.max(b);
    for (index, band) in gear.bands.iter().enumerate() {
        for end in 0..2 {
            let z = band.z[end];
            if end == 0 && index > 0 && gear.bands[index - 1].z[1] == z {
                continue;
            }
            let next = (end == 1)
                .then(|| gear.bands.get(index + 1))
                .flatten()
                .filter(|next| next.z[0] == z);
            let beyond = gear.material(z, if end == 0 { -1. } else { 1. });
            for (position, flank) in gear.flanks.iter().enumerate() {
                let smooth = flank.coaxial().is_some_and(|radius| {
                    next.is_some()
                        || beyond == Some(Region::Disk(radius))
                        || matches!(beyond, Some(Region::Disk(other)) if same(other, radius))
                });
                let mut extra = std::mem::take(&mut marks.rims[index][end][position]);
                if next.is_some() {
                    extra.extend(std::mem::take(&mut marks.rims[index + 1][0][position]));
                }
                ring(
                    flank.edge,
                    z,
                    band.phase_at(z),
                    if smooth { Kind::Smooth } else { Kind::Edge },
                    extra,
                    &mut output,
                );
            }
        }
    }

    // Rails: the helix of every profile vertex through every band.
    for (band_index, band) in gear.bands.iter().enumerate() {
        for (index, flank) in gear.flanks.iter().enumerate() {
            let before = &gear.flanks[(index + count - 1) % count];
            let [a, b] = [before.tangent(1.), flank.tangent(0.)];
            let sharp = (a[0] * b[1] - a[1] * b[0]).abs() > 1e-9 || dot2(a, b) < 0.;
            let kind = if sharp {
                Kind::Edge
            } else if axial && (before.coaxial().is_some() != flank.coaxial().is_some()) {
                // Seen along the axis the cylinder beside this seam is
                // edge-on, so the seam bounds the image of the other flank.
                Kind::Outline
            } else {
                Kind::Smooth
            };
            let height = band.z[1] - band.z[0];
            for pair in stops(std::mem::take(&mut marks.rails[band_index][index])).windows(2) {
                output.push(Candidate {
                    kind,
                    path: Arc::new(Helix {
                        frame: gear.frame,
                        band: *band,
                        point: flank.point(0.),
                        z: [band.z[0] + pair[0] * height, band.z[0] + pair[1] * height],
                        axial,
                    }),
                });
            }
        }
    }

    // Bore and hub: rims, and the two generators perpendicular to the sight.
    let mut cylinders: Vec<(f64, [f64; 2])> = Vec::new();
    if gear.bore > 0. {
        cylinders.push((gear.bore, gear.axis));
        for z in gear.axis {
            ring(
                Gear::circle(gear.bore, true),
                z,
                0.,
                Kind::Edge,
                Vec::new(),
                &mut output,
            );
        }
    }
    for cylinder in &gear.cylinders {
        cylinders.push((cylinder.radius, cylinder.z));
        for (end, z) in cylinder.z.into_iter().enumerate() {
            let beyond = gear.material(z, if end == 0 { -1. } else { 1. });
            let (arcs, phase) = match beyond {
                Some(Region::Disk(other)) if same(other, cylinder.radius) => continue,
                // Where the cylinder meets a tooth band, the band's own
                // coaxial arcs of the same radius continue it; the rest of
                // its rim is the edge under each tooth.
                Some(Region::Profile(phase)) => (uncovered(gear, cylinder.radius), phase),
                _ => (vec![Gear::circle(cylinder.radius, false)], 0.),
            };
            for edge in arcs {
                ring(edge, z, phase, Kind::Edge, Vec::new(), &mut output);
            }
        }
    }
    if !axial {
        for (radius, z) in cylinders {
            for side in [FRAC_PI_2, -FRAC_PI_2] {
                let (sin, cos) = (gear.azimuth + side).portable_sin_cos();
                let at = |z: f64| gear.frame.place([radius * cos, radius * sin, z]);
                output.push(Candidate {
                    kind: Kind::Outline,
                    path: Arc::new(Segment(at(z[0]), at(z[1]))),
                });
            }
        }
    }
    Ok(output)
}
