//! Edges, smooth seams and outlines of a twisted gear in one view.
use super::super::paths::{ProfilePath, Segment, dot2};
use super::super::profile::folds;
use super::super::*;
use super::flank::rotate;
use super::trace::trace;
use super::{Gear, Helix, Level, Patch, Traced};
use std::f64::consts::{FRAC_PI_2, TAU};

/// Where outlines end on the rims and rails of the flanks, so that every
/// edge is split where the adjacent flank turns away from the viewer.
struct Marks {
    /// By axial level and profile edge: edge fractions.
    rims: Vec<Vec<Vec<f64>>>,
    /// By band and profile vertex: band fractions.
    rails: Vec<Vec<Vec<f64>>>,
}

impl Marks {
    fn add(&mut self, band: usize, flank: usize, [s, zeta]: [f64; 2]) {
        let count = self.rims[0].len();
        const EDGE: f64 = 1e-7;
        if zeta <= EDGE {
            self.rims[band][flank].push(s);
        } else if zeta >= 1. - EDGE {
            self.rims[band + 1][flank].push(s);
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
    for omega in [ratio.acos(), -ratio.acos()] {
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

/// Parameters at which the paper image of a traced outline has a cusp: its
/// tangent passes through the sight direction, where its visibility changes.
fn cusps(path: &Traced) -> Vec<f64> {
    let measure = |t: f64| {
        let h = 1e-6;
        let [a, b] = [(t - h).max(0.), (t + h).min(1.)];
        let [p, q] = [path.point(a), path.point(b)];
        let [s, zeta] = path.parameters(t);
        let n = path.patch.normal(s, zeta);
        ((q[0] - p[0]) * n[1] - (q[1] - p[1]) * n[0]) / (b - a)
    };
    let count = 8 * path.steps.len();
    let mut output = Vec::new();
    let mut previous = (0., measure(0.));
    for i in 1..=count {
        let t = i as f64 / count as f64;
        let value = measure(t);
        if (value < 0.) != (previous.1 < 0.) {
            let [mut a, mut b] = [previous.0, t];
            for _ in 0..60 {
                let mid = 0.5 * (a + b);
                if (measure(mid) < 0.) == (previous.1 < 0.) {
                    a = mid;
                } else {
                    b = mid;
                }
            }
            output.push(0.5 * (a + b));
        }
        previous = (t, value);
    }
    output
}

pub(super) fn candidates(gear: &Gear) -> Result<Vec<Candidate>, GeometryError> {
    let count = gear.flanks.len();
    let axial = gear.axial();
    let mut output: Vec<Candidate> = Vec::new();
    let mut marks = Marks {
        rims: vec![vec![Vec::new(); count]; gear.bands.len() + 1],
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
                    let steps = Arc::new(chain.steps);
                    let whole = Traced {
                        patch,
                        level: Level::Fold,
                        steps: steps.clone(),
                        range: [0., 1.],
                    };
                    let mut cuts = vec![0.];
                    cuts.extend(cusps(&whole));
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
    let ring = |edge: multicad_rust_kernel::spur_gear::Edge,
                z: f64,
                phase: f64,
                kind: Kind,
                extra: Vec<f64>,
                output: &mut Vec<Candidate>| {
        let mut marks = extra;
        if !axial {
            marks.extend(folds(edge, rotate(across, -phase)));
        }
        if matches!(edge, multicad_rust_kernel::spur_gear::Edge::Planar(e) if e.start() == e.end())
        {
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
    for level in 0..=gear.bands.len() {
        let band = &gear.bands[level.min(gear.bands.len() - 1)];
        let z = if level == gear.bands.len() {
            band.z[1]
        } else {
            band.z[0]
        };
        let interior = level != 0 && level != gear.bands.len();
        for (index, flank) in gear.flanks.iter().enumerate() {
            let kind = if interior && flank.coaxial().is_some() {
                Kind::Smooth
            } else {
                Kind::Edge
            };
            ring(
                flank.edge,
                z,
                band.phase_at(z),
                kind,
                std::mem::take(&mut marks.rims[level][index]),
                &mut output,
            );
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
    for extension in &gear.extensions {
        cylinders.push((extension.radius, extension.z));
        for z in extension.z {
            ring(
                Gear::circle(extension.radius, false),
                z,
                0.,
                Kind::Edge,
                Vec::new(),
                &mut output,
            );
        }
    }
    if !axial {
        for (radius, z) in cylinders {
            for side in [FRAC_PI_2, -FRAC_PI_2] {
                let (sin, cos) = (gear.azimuth + side).sin_cos();
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
