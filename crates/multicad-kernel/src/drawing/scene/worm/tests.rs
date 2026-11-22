//! The worm members proven against the provider's own membership predicate
//! and surface laws, independently of the scene engine.
use super::*;
use crate::{GeneratedZaWormWheel, RigidTransform, ZaWorm};
use multicad_core::document::EntityId;
use multicad_rust_kernel::worm_gear::Containment;

fn worm(m: f64, length: f64, bore: f64) -> Solid {
    Solid::rust_za_worm_named(
        ZaWorm {
            axial_module_mm: m,
            threaded_length_mm: length * m,
            bore_radius_mm: bore,
            axial_thinning_mm: 0.02 * m,
        },
        EntityId::new(),
        GeometryPolicy::default(),
    )
    .unwrap()
}

fn wheel(m: f64, teeth: u32, bore: f64) -> Solid {
    Solid::rust_generated_za_worm_wheel_named(
        GeneratedZaWormWheel {
            teeth,
            axial_module_mm: m,
            face_width_mm: 1.6 * m,
            bore_radius_mm: bore,
            pitch_thinning_mm: 0.03 * m,
        },
        EntityId::new(),
        GeometryPolicy::default(),
    )
    .unwrap()
}

fn posed(solid: Solid) -> Solid {
    let (s, c) = 0.3_f64.sin_cos();
    solid
        .transformed(RigidTransform {
            rotation: [[c, -s, 0.], [s, c, 0.], [0., 0., 1.]],
            translation_mm: [1.5, -2., 0.75],
        })
        .unwrap()
}

fn solids() -> Vec<(&'static str, Solid)> {
    vec![
        ("worm", worm(1., 20., 0.)),
        ("bored worm", posed(worm(1.5, 16., 2.5))),
        ("wheel", wheel(1., 32, 0.)),
        ("bored wheel", posed(wheel(1.25, 40, 6.))),
    ]
}

fn views() -> Vec<(&'static str, PlaneFrame)> {
    let origin = [0.3, -0.2, 0.1];
    let s = 0.5_f64.sqrt();
    vec![
        (
            "axial",
            PlaneFrame {
                origin,
                u: [1., 0., 0.],
                v: [0., 1., 0.],
            },
        ),
        (
            "reversed",
            PlaneFrame {
                origin,
                u: [-1., 0., 0.],
                v: [0., 1., 0.],
            },
        ),
        (
            "side",
            PlaneFrame {
                origin,
                u: [0., 1., 0.],
                v: [0., 0., 1.],
            },
        ),
        (
            "top",
            PlaneFrame {
                origin,
                u: [1., 0., 0.],
                v: [0., 0., 1.],
            },
        ),
        (
            "oblique",
            PlaneFrame {
                origin,
                u: [s, -s, 0.],
                v: [s / 3_f64.sqrt(), s / 3_f64.sqrt(), -2. * s / 3_f64.sqrt()],
            },
        ),
        (
            "tilted",
            PlaneFrame {
                origin,
                u: [1., 0., 0.],
                v: [0., 0.96, 0.28],
            },
        ),
    ]
}

fn built(solid: &Solid, frame: PlaneFrame) -> Body {
    body(solid, &View::new(frame)).unwrap().unwrap()
}

fn classify(solid: &Solid, local: V3, tolerance: f64) -> Containment {
    match solid.rust_worm() {
        Some((worm, _)) => worm.classify(local, tolerance),
        None => solid
            .rust_worm_wheel()
            .unwrap()
            .0
            .classify(local, tolerance),
    }
}

struct Random(u64);

impl Random {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1_u64 << 53) as f64
    }
}

#[test]
fn sight_intervals_are_the_membership_predicate_with_roots_on_the_surfaces() {
    let mut random = Random(7);
    for (name, solid) in solids() {
        for (label, frame) in views() {
            let body = built(&solid, frame);
            let inner = &body.0;
            let [lo, hi] = body.bounds();
            let mut inside = 0;
            let mut ends = 0;
            for _ in 0..400 {
                let paper = [0, 1].map(|k| lo[k] + (hi[k] - lo[k]) * random.next());
                let spans = body.sight(paper).unwrap();
                for pair in spans.windows(2) {
                    assert!(pair[0][1] < pair[1][0], "{name} {label} {spans:?}");
                }
                for span in &spans {
                    assert!(span[1] > span[0]);
                    for depth in span {
                        let local = inner.local([paper[0], paper[1], *depth]);
                        assert_eq!(
                            classify(&solid, local, 1e-9),
                            Containment::Boundary,
                            "{name} {label} {paper:?} {depth}"
                        );
                        ends += 1;
                    }
                }
                for _ in 0..24 {
                    let depth = lo[2] - 0.5 + (hi[2] - lo[2] + 1.) * random.next();
                    let local = inner.local([paper[0], paper[1], depth]);
                    let found = spans.iter().any(|[a, b]| *a < depth && depth < *b);
                    let expected = match classify(&solid, local, 1e-6) {
                        Containment::Inside => true,
                        Containment::Outside => false,
                        Containment::Boundary => continue,
                    };
                    assert_eq!(
                        found, expected,
                        "{name} {label} {paper:?} {depth} {spans:?}"
                    );
                    inside += usize::from(found);
                }
            }
            assert!(inside > 200 && ends > 100, "{name} {label} {inside} {ends}");
        }
    }
}

/// Distance from a point to a polyline.
fn distance(point: P3, line: &[P3]) -> f64 {
    line.windows(2)
        .map(|pair| {
            let d = sub(pair[1], pair[0]);
            let s = (dot3(sub(point, pair[0]), d) / dot3(d, d).max(1e-300)).clamp(0., 1.);
            norm(sub(point, add(pair[0], scaled(d, s))))
        })
        .fold(f64::INFINITY, f64::min)
}

#[test]
fn candidates_lie_on_the_boundary_and_outlines_are_complete() {
    for (name, solid) in solids() {
        for (label, frame) in views() {
            let body = built(&solid, frame);
            let inner = &body.0;
            let candidates = body.candidates(0.002).unwrap();
            let mut outlines: Vec<Vec<P3>> = Vec::new();
            for candidate in &candidates {
                let stops = candidate.path.flatten(0.002);
                assert!(stops[0] == 0. && *stops.last().unwrap() == 1.);
                let points: Vec<P3> = stops.iter().map(|&t| candidate.path.point(t)).collect();
                for (i, pair) in stops.windows(2).enumerate() {
                    assert!(pair[1] > pair[0]);
                    for s in [0., 0.37, 1.] {
                        let p = candidate.path.point(pair[0] + (pair[1] - pair[0]) * s);
                        assert_eq!(
                            classify(&solid, inner.local(p), 1e-8),
                            Containment::Boundary,
                            "{name} {label} {:?} {i}",
                            candidate.kind
                        );
                    }
                    // Chords stay within the tolerance of the path.
                    let mid = candidate.path.point(0.5 * (pair[0] + pair[1]));
                    let chord = scaled(add(points[i], points[i + 1]), 0.5);
                    assert!(norm(sub(mid, chord)) <= 0.002, "{name} {label}");
                }
                // Every part is injective in the paper along one paper axis.
                let monotone = |k: usize| {
                    points.windows(2).all(|p| p[1][k] >= p[0][k] - 1e-9)
                        || points.windows(2).all(|p| p[1][k] <= p[0][k] + 1e-9)
                };
                assert!(
                    monotone(0) || monotone(1),
                    "{name} {label} {:?}",
                    candidate.kind
                );
                if candidate.kind == Kind::Outline {
                    outlines.push(points);
                }
            }
            // The exact normal is perpendicular to the sight line along
            // every traced outline.
            let shell = &inner.shell;
            let half = shell.half();
            let mut traced = 0;
            let mut crossings = 0;
            for patch in 0..shell.patches() {
                for path in body.traced(patch, Level::Outline).unwrap() {
                    for i in 0..=20 {
                        let jet = path.jet(f64::from(i) / 20.);
                        let facing = dot3(unit(jet.normal()), inner.toward);
                        assert!(facing.abs() < 1e-9, "{name} {label} {patch} {facing}");
                    }
                    traced += 1;
                }
                // Dense sign changes of `n · d` on the patch all lie on a
                // candidate outline.
                let [nu, nz] = shell.grid(patch).map(|n| 5 * n + 3);
                let facing =
                    |u: f64, z: f64| dot3(unit(shell.eval(patch, u, z).normal()), inner.toward);
                let node = |i: usize, j: usize| {
                    (
                        i as f64 / nu as f64,
                        -half + 2. * half * j as f64 / nz as f64,
                    )
                };
                for i in 0..nu {
                    for j in 0..=nz {
                        let [(ua, z), (ub, _)] = [node(i, j), node(i + 1, j)];
                        let [fa, fb] = [facing(ua, z), facing(ub, z)];
                        if fa.abs() < 1e-7 || fb.abs() < 1e-7 || (fa < 0.) == (fb < 0.) {
                            continue;
                        }
                        let u = root(|u| facing(u, z), ua, ub);
                        let point = inner.view(shell.eval(patch, u, z).p);
                        let gap = outlines
                            .iter()
                            .map(|line| distance(point, line))
                            .fold(f64::INFINITY, f64::min);
                        assert!(gap < 0.004, "{name} {label} {patch} {u} {z} {gap}");
                        crossings += 1;
                    }
                }
            }
            let axial = inner.toward[0].hypot(inner.toward[1]) < 1e-9;
            assert!(axial || (traced > 0 && crossings > 0), "{name} {label}");
        }
    }
}

/// Signed area of oriented pieces by Green's theorem on fine chords.
fn area(pieces: &[Arc<dyn Path>]) -> f64 {
    let mut total = 0.;
    for piece in pieces {
        let points: Vec<P3> = piece
            .flatten(2e-5)
            .into_iter()
            .map(|t| piece.point(t))
            .collect();
        for pair in points.windows(2) {
            total += 0.5 * (pair[0][0] * pair[1][1] - pair[1][0] * pair[0][1]);
        }
    }
    total
}

/// Area of the membership predicate over the plane at a depth: exact
/// transitions along rows, midpoint rule across them. The rows are turned
/// in the paper so that no straight boundary runs along a row.
fn quadrature(solid: &Solid, body: &Body, depth: f64) -> f64 {
    let inner = &body.0;
    let [lo, hi] = body.bounds();
    // The wheel predicate bisects its section at every call.
    let (rows, columns) = if solid.rust_worm().is_some() {
        (1500, 3000)
    } else {
        (400, 800)
    };
    let centre = [0, 1].map(|k| 0.5 * (lo[k] + hi[k]));
    let reach = 0.5 * (hi[0] - lo[0]).hypot(hi[1] - lo[1]) + 1e-3;
    let (sin, cos) = 0.37_f64.sin_cos();
    let inside = |x: f64, y: f64| {
        let paper = [
            centre[0] + cos * x - sin * y,
            centre[1] + sin * x + cos * y,
            depth,
        ];
        classify(solid, inner.local(paper), 0.) == Containment::Inside
    };
    let mut total = 0.;
    let height = 2. * reach / rows as f64;
    let width = 2. * reach / columns as f64;
    for row in 0..rows {
        let y = -reach + (row as f64 + 0.5) * height;
        let mut previous = (-reach - width, false);
        let mut entered = 0.;
        for column in 0..=columns + 1 {
            let x = -reach + column as f64 * width;
            let state = inside(x, y);
            if state != previous.1 {
                let (mut a, mut b) = (previous.0, x);
                for _ in 0..50 {
                    let mid = 0.5 * (a + b);
                    if inside(mid, y) == previous.1 {
                        a = mid;
                    } else {
                        b = mid;
                    }
                }
                if state {
                    entered = a;
                } else {
                    total += (a - entered) * height;
                }
            }
            previous = (x, state);
        }
    }
    total
}

#[test]
fn sections_are_closed_oriented_boundaries_of_the_membership_predicate() {
    for (name, solid) in solids() {
        for (label, frame) in views() {
            let body = built(&solid, frame);
            let inner = &body.0;
            let [lo, hi] = body.bounds();
            for fraction in [0.31, 0.5, 0.83] {
                let depth = lo[2] + (hi[2] - lo[2]) * fraction;
                let pieces = body.cut(depth).unwrap();
                assert!(!pieces.is_empty(), "{name} {label} {fraction}");
                let ends: Vec<[P3; 2]> =
                    pieces.iter().map(|p| [p.point(0.), p.point(1.)]).collect();
                for (i, piece) in pieces.iter().enumerate() {
                    for step in 0..=64 {
                        let t = f64::from(step) / 64.;
                        let p = piece.point(t);
                        assert!((p[2] - depth).abs() < 1e-9, "{name} {label}");
                        assert_eq!(
                            classify(&solid, inner.local(p), 1e-8),
                            Containment::Boundary,
                            "{name} {label} {fraction} {i} {t}"
                        );
                    }
                    // Every end continues in exactly one other start.
                    let next = ends
                        .iter()
                        .filter(|other| norm(sub(other[0], ends[i][1])) < 1e-7)
                        .count();
                    assert_eq!(next, 1, "{name} {label} {fraction} {i} {:?}", ends[i]);
                    // Material lies on the left of the travel direction.
                    let [a, b] = [piece.point(0.5 - 1e-4), piece.point(0.5 + 1e-4)];
                    let travel = unit(sub(b, a));
                    let mid = piece.point(0.5);
                    let side = |sign: f64| {
                        let p = [
                            mid[0] - sign * 2e-4 * travel[1],
                            mid[1] + sign * 2e-4 * travel[0],
                            depth,
                        ];
                        classify(&solid, inner.local(p), 0.)
                    };
                    assert_eq!(
                        side(1.),
                        Containment::Inside,
                        "{name} {label} {fraction} {i}"
                    );
                    assert_eq!(
                        side(-1.),
                        Containment::Outside,
                        "{name} {label} {fraction} {i}"
                    );
                }
                let exact = area(&pieces);
                let sampled = quadrature(&solid, &body, depth);
                assert!(
                    (exact - sampled).abs() < 1e-3 * sampled,
                    "{name} {label} {fraction} {exact} {sampled}"
                );
            }
        }
    }
}

#[test]
fn the_axial_plane_section_of_a_worm_is_its_straight_sided_profile() {
    let solid = worm(1., 20., 2.);
    let (source, _) = solid.rust_worm().unwrap();
    let thread = source.thread();
    let frame = PlaneFrame {
        origin: [0.; 3],
        u: [0., 0., 1.],
        v: [1., 0., 0.],
    };
    let body = built(&solid, frame);
    let pieces = body.cut(0.).unwrap();
    let slopes = [1, 3].map(|band| thread.slope_at_phase(0., band));
    let mut flanks = 0;
    for piece in &pieces {
        // Paper x is the axis and paper y the signed radius.
        let points: Vec<P3> = (0..=8).map(|i| piece.point(f64::from(i) / 8.)).collect();
        let [a, b] = [points[0], points[8]];
        for (i, p) in points.iter().enumerate() {
            let s = i as f64 / 8.;
            let chord = add(scaled(a, 1. - s), scaled(b, s));
            let off = (p[0] - chord[0]) * (b[1] - a[1]) - (p[1] - chord[1]) * (b[0] - a[0]);
            assert!(off.abs() < 1e-12 * 36., "{off}");
            assert!(p[2].abs() < 5e-12);
        }
        let radial = (b[1].abs() - a[1].abs()).abs();
        if radial > 1. && (b[0] - a[0]).abs() > 0.1 {
            // A flank: the radius changes by the profile slope per axial mm.
            let slope = (b[1].abs() - a[1].abs()) / (b[0] - a[0]);
            assert!(
                slopes.iter().any(|s| (s - slope).abs() < 1e-11 * s.abs()),
                "{slope} {slopes:?}"
            );
            flanks += 1;
        }
    }
    // Two flanks per turn on both sides of the axis, less the cut turns.
    assert!(flanks >= 4 * 5, "{flanks}");
}

/// The generated wheel and oblique frame of `tests/rust_scene_worm.rs`.
fn compared_wheel() -> (Solid, PlaneFrame) {
    let (s, c) = 0.1_f64.sin_cos();
    let solid = Solid::rust_generated_za_worm_wheel_named(
        GeneratedZaWormWheel {
            teeth: 32,
            axial_module_mm: 1.,
            face_width_mm: 1.6,
            bore_radius_mm: 5.,
            pitch_thinning_mm: 0.03,
        },
        EntityId::new(),
        GeometryPolicy::default(),
    )
    .unwrap()
    .transformed(RigidTransform {
        rotation: [[c, -s, 0.], [s, c, 0.], [0., 0., 1.]],
        translation_mm: [0.; 3],
    })
    .unwrap();
    let r = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    (
        solid,
        PlaneFrame {
            origin: [0.3, -0.2, 0.1],
            u: [r, -r, 0.],
            v: [t, t, -2. * t],
        },
    )
}

/// Boundary crossings of the exact membership predicate along the sight line
/// through `paper`, between depths `from` and `to` in `steps` samples.
fn crossings(inside: &dyn Fn(P3) -> bool, paper: P2, [from, to]: [f64; 2], steps: usize) -> usize {
    let mut last = inside([paper[0], paper[1], from]);
    let mut count = 0;
    for i in 1..=steps {
        let z = from + (to - from) * i as f64 / steps as f64;
        let now = inside([paper[0], paper[1], z]);
        count += usize::from(now != last);
        last = now;
    }
    count
}

/// Roles of the published vectors that run along a curve through a paper
/// point: within the publishing chord tolerance of it and parallel to the
/// curve's paper tangent there, so that crossing curves do not count.
fn roles_at(vectors: &[ProjectedVector], paper: P2, tangent: P2) -> Vec<VectorRole> {
    // Published chords stay within half of the default 0.01 mm chord
    // tolerance of the exact curve.
    let reach = 0.0055;
    let mut roles = Vec::new();
    for vector in vectors {
        let near = vector.preview_mm.windows(2).any(|pair| {
            let [a, b] = [pair[0], pair[1]];
            let d = [b[0] - a[0], b[1] - a[1]];
            let length = d[0].hypot(d[1]);
            let s = (((paper[0] - a[0]) * d[0] + (paper[1] - a[1]) * d[1])
                / (length * length).max(1e-300))
            .clamp(0., 1.);
            let gap = (a[0] + s * d[0] - paper[0]).hypot(a[1] + s * d[1] - paper[1]);
            let sine = (d[0] * tangent[1] - d[1] * tangent[0]).abs() / length.max(1e-300);
            gap <= reach && sine < 0.2
        });
        if near && !roles.contains(&vector.role) {
            roles.push(vector.role);
        }
    }
    roles
}

/// Material angle in degrees around a curve point, from the exact predicate
/// on a small circle in the plane normal to the curve: 180 on a smooth
/// surface, otherwise the dihedral angle of the crease.
fn material_angle(own: &dyn Fn(P3) -> bool, path: &dyn Path, t: f64, radius: f64) -> f64 {
    let p = path.point(t);
    let tangent = unit(sub(
        path.point((t + 1e-5).min(1.)),
        path.point((t - 1e-5).max(0.)),
    ));
    let seed = if tangent[0].abs() < 0.9 {
        [1., 0., 0.]
    } else {
        [0., 1., 0.]
    };
    let e1 = unit(cross(tangent, seed));
    let e2 = cross(tangent, e1);
    let samples = 1440;
    let inside = (0..samples)
        .filter(|&i| {
            let (s, c) = (TAU * f64::from(i) / f64::from(samples)).sin_cos();
            own(add(p, add(scaled(e1, radius * c), scaled(e2, radius * s))))
        })
        .count();
    360. * inside as f64 / f64::from(samples)
}

/// Material counts for the proofs only where the provider's margin exceeds
/// this many millimetres: a march along a surface it grazes (a bore wall or
/// the root cylinder seen along the axis) meets no material, as the member
/// contract requires, while real material lies at least 0.01 mm deep.
const STRICT: f64 = 1e-9;

/// One member of a scene with its candidates and its own exact predicate.
type ScenePart<'a> = (&'a str, Vec<Candidate>, &'a dyn Fn(P3) -> bool);

/// What one sample of a candidate proves.
#[derive(Default)]
struct Proof {
    folds: usize,
    creases: usize,
    seams: usize,
    roles: usize,
    unpublished: usize,
}

/// A sample point of a candidate is proved against the exact predicates:
///
/// - an outline is a fold: sight lines just either side of it cross the
///   boundary of its solid a number of times differing by two near its
///   depth (offsets of 1e-5 mm keep the fold sliver inside the window even
///   on the 21 mm tip cylinder), so `n · d` changes sign across it;
/// - an edge is a crease: the material angle around it differs from a flat
///   180 degrees by more than five, and a smooth seam is flat;
/// - the published role is hidden exactly when marching from just in front
///   of the point toward the viewer meets material of any solid.
#[allow(clippy::too_many_arguments)]
fn prove(
    kind: Kind,
    path: &dyn Path,
    t: f64,
    own: &dyn Fn(P3) -> bool,
    scene: &dyn Fn(P3) -> bool,
    top: f64,
    vectors: &[ProjectedVector],
    proof: &mut Proof,
    label: &str,
) {
    let p = path.point(t);
    match kind {
        Kind::Outline => {
            let [a, b] = [path.point(t - 1e-4), path.point(t + 1e-4)];
            let tangent = [b[0] - a[0], b[1] - a[1]];
            let length = tangent[0].hypot(tangent[1]);
            assert!(
                length > 1e-9,
                "{label} t={t}: outline without a paper tangent"
            );
            let side = [-tangent[1] / length, tangent[0] / length];
            let [left, right] = [1., -1.].map(|sign: f64| {
                crossings(
                    own,
                    [p[0] + sign * 1e-5 * side[0], p[1] + sign * 1e-5 * side[1]],
                    [p[2] - 0.25, p[2] + 0.25],
                    10_000,
                )
            });
            assert_eq!(
                left.abs_diff(right),
                2,
                "{label} t={t}: no fold across the outline at {p:?} ({left} vs {right})"
            );
            proof.folds += 1;
        }
        Kind::Edge => {
            let angle = material_angle(own, path, t, 1e-3);
            assert!(
                (angle - 180.).abs() > 5.,
                "{label} t={t}: edge at {p:?} is flat ({angle} degrees)"
            );
            proof.creases += 1;
        }
        Kind::Smooth => {
            let angle = material_angle(own, path, t, 1e-3);
            assert!(
                (angle - 180.).abs() <= 1.,
                "{label} t={t}: smooth seam at {p:?} is a crease ({angle} degrees)"
            );
            proof.seams += 1;
            return;
        }
    }
    let [a, b] = [
        path.point((t - 1e-4).max(0.)),
        path.point((t + 1e-4).min(1.)),
    ];
    let tangent = [b[0] - a[0], b[1] - a[1]];
    let length = tangent[0].hypot(tangent[1]);
    if length < 1e-12 {
        // The curve runs along the sight line here.
        proof.unpublished += 1;
        return;
    }
    let roles = roles_at(
        vectors,
        [p[0], p[1]],
        [tangent[0] / length, tangent[1] / length],
    );
    if roles.is_empty() {
        // Pieces shorter than the engine's merge distance near a visibility
        // transition or a crossing publish without this exact point.
        proof.unpublished += 1;
        return;
    }
    let start = p[2] + 1e-3;
    let steps = (((top - start) / 2e-3).ceil() as usize).max(1);
    let covered = (0..=steps).any(|i| {
        let z = start + (top - start) * i as f64 / steps as f64;
        scene([p[0], p[1], z])
    });
    let expected = if covered {
        VectorRole::Hidden
    } else {
        VectorRole::Visible
    };
    assert!(
        roles.contains(&expected),
        "{label} t={t}: published {roles:?} at {p:?}, marching says {expected:?}"
    );
    proof.roles += 1;
}

/// Whether a candidate's paper image passes within `reach` of a point.
fn passes(c: &Candidate, target: P2, reach: f64) -> bool {
    c.path.flatten(0.0025).windows(2).any(|pair| {
        let [a, b] = [c.path.point(pair[0]), c.path.point(pair[1])];
        let d = [b[0] - a[0], b[1] - a[1]];
        let s = (((target[0] - a[0]) * d[0] + (target[1] - a[1]) * d[1])
            / (d[0] * d[0] + d[1] * d[1]).max(1e-300))
        .clamp(0., 1.);
        (a[0] + s * d[0] - target[0]).hypot(a[1] + s * d[1] - target[1]) < reach
    })
}

/// The curves behind the measured OpenCascade deviations of the oblique wheel
/// view (a hidden edge meeting a corner 0.60 mm from any native curve, a
/// visibility change at a vertex 0.035 mm from the native one) and a sample
/// of every other edge and outline of the view are proved against the
/// provider's predicate.
#[test]
fn oblique_wheel_curves_are_creases_and_folds_with_marched_visibility() {
    let (solid, frame) = compared_wheel();
    let body = built(&solid, frame);
    let inner = &body.0;
    let (wheel, _) = solid.rust_worm_wheel().unwrap();
    let own = |p: P3| wheel.classify(inner.local(p), STRICT) == Containment::Inside;
    let top = body.bounds()[1][2];
    let vectors = crate::project_drawing(
        &[crate::DrawingSolid {
            solid: &solid,
            hatch: None,
        }],
        crate::DrawingProjection {
            frame,
            hidden_lines: true,
            ..Default::default()
        },
    )
    .unwrap();
    let candidates = body.candidates(0.0025).unwrap();
    let mut proof = Proof::default();
    for (label, target) in [
        ("hidden corner", [6.568499986795132, 9.658492277809218]),
        ("visible end", [-4.217254641741771, 9.149364469823032]),
    ] {
        let found: Vec<&Candidate> = candidates
            .iter()
            .filter(|c| passes(c, target, 1e-6))
            .collect();
        assert!(
            found.len() >= 3,
            "{label}: {} curves through {target:?}",
            found.len()
        );
        for c in found {
            for i in 1..20 {
                let t = f64::from(i) / 20.;
                prove(
                    c.kind,
                    c.path.as_ref(),
                    t,
                    &own,
                    &own,
                    top,
                    &vectors,
                    &mut proof,
                    label,
                );
            }
        }
    }
    let outlines: Vec<&Candidate> = candidates
        .iter()
        .filter(|c| c.kind == Kind::Outline)
        .collect();
    let mut random = Random(11);
    for _ in 0..40 {
        let c = outlines[(random.next() * outlines.len() as f64) as usize];
        let t = 0.05 + 0.9 * random.next();
        prove(
            c.kind,
            c.path.as_ref(),
            t,
            &own,
            &own,
            top,
            &vectors,
            &mut proof,
            "outline",
        );
    }
    for _ in 0..80 {
        let c = &candidates[(random.next() * candidates.len() as f64) as usize];
        let t = 0.05 + 0.9 * random.next();
        prove(
            c.kind,
            c.path.as_ref(),
            t,
            &own,
            &own,
            top,
            &vectors,
            &mut proof,
            "sample",
        );
    }
    println!(
        "folds {} creases {} seams {} roles {} unpublished {}",
        proof.folds, proof.creases, proof.seams, proof.roles, proof.unpublished
    );
    assert!(
        proof.folds >= 30 && proof.creases >= 50 && proof.roles >= 120,
        "too few proved samples"
    );
}

/// The box edges behind the measured deviations of the oblique box-and-wheel
/// scene (0.28 mm visible, 1.27 mm hidden) are proved the same way, with the
/// interpenetrating wheel and the exact box both covering them.
#[test]
fn oblique_box_and_wheel_edges_have_marched_visibility() {
    let (wheel_solid, _) = compared_wheel();
    let wheel_solid = wheel_solid
        .transformed(RigidTransform {
            rotation: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            translation_mm: [3., 2., 1.],
        })
        .unwrap();
    let r = |v: i64| multicad_rust_kernel::Rational::new(v, 1).unwrap();
    let block = Solid::rust_box(
        EntityId::new().uuid_bits(),
        [r(8), r(6), r(5)],
        PlaneFrame::default(),
        GeometryPolicy::default(),
    )
    .unwrap();
    let s = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    let frame = PlaneFrame {
        origin: [1., 1., 1.],
        u: [s, -s, 0.],
        v: [t, t, -2. * t],
    };
    let view = View::new(frame);
    let gear = built(&wheel_solid, frame);
    let (wheel, _) = wheel_solid.rust_worm_wheel().unwrap();
    let bounds = block.bounds().unwrap();
    let world = |p: P3| -> [f64; 3] {
        std::array::from_fn(|i| {
            view.origin[i]
                + p[0] * view.axes[0][i]
                + p[1] * view.axes[1][i]
                + p[2] * view.axes[2][i]
        })
    };
    let in_box = |p: P3| {
        let w = world(p);
        (0..3).all(|i| bounds.min[i] + STRICT < w[i] && w[i] < bounds.max[i] - STRICT)
    };
    let in_wheel = |p: P3| wheel.classify(gear.0.local(p), STRICT) == Containment::Inside;
    let scene = |p: P3| in_box(p) || in_wheel(p);
    let member = super::super::member(&block, &view).unwrap().unwrap();
    let top = member.bounds()[1][2].max(gear.bounds()[1][2]);
    let mut proof = Proof::default();
    for (label, target, hidden) in [
        ("visible box edge", [0., 3.429285639896448], false),
        (
            "hidden box edge",
            [std::f64::consts::SQRT_2, 4.24578222082414],
            true,
        ),
    ] {
        let vectors = crate::project_drawing(
            &[
                crate::DrawingSolid {
                    solid: &block,
                    hatch: None,
                },
                crate::DrawingSolid {
                    solid: &wheel_solid,
                    hatch: None,
                },
            ],
            crate::DrawingProjection {
                frame,
                hidden_lines: hidden,
                ..Default::default()
            },
        )
        .unwrap();
        let mut found = 0;
        let members: [ScenePart<'_>; 2] = [
            ("box", member.candidates(0.0025).unwrap(), &in_box),
            ("wheel", gear.candidates(0.0025).unwrap(), &in_wheel),
        ];
        for (name, candidates, own) in &members {
            for c in candidates.iter().filter(|c| passes(c, target, 3e-3)) {
                println!("{label}: {name} {:?}", c.kind);
                found += 1;
                for i in 1..40 {
                    let t = f64::from(i) / 40.;
                    // Without hidden lines a covered point publishes nothing
                    // and counts as unpublished.
                    prove(
                        c.kind,
                        c.path.as_ref(),
                        t,
                        *own,
                        &scene,
                        top,
                        &vectors,
                        &mut proof,
                        label,
                    );
                }
            }
        }
        if found == 0 {
            for (name, candidates, _) in &members {
                for c in candidates.iter().filter(|c| passes(c, target, 5e-2)) {
                    println!(
                        "{label}: near {name} {:?} {:?} {:?}",
                        c.kind,
                        c.path.point(0.),
                        c.path.point(1.)
                    );
                }
            }
            for v in &vectors {
                if v.preview_mm
                    .iter()
                    .any(|p| (p[0] - target[0]).hypot(p[1] - target[1]) < 5e-2)
                {
                    println!(
                        "{label}: vector {} {:?} {:?}",
                        v.source_index, v.role, v.preview_mm
                    );
                }
            }
        }
        assert!(found > 0, "{label}: no curve through {target:?}");
    }
    println!(
        "box scene folds {} creases {} roles {} unpublished {}",
        proof.folds, proof.creases, proof.roles, proof.unpublished
    );
    assert!(
        proof.creases + proof.folds >= 40 && proof.roles >= 20,
        "too few proved samples"
    );
}

/// The axial worm detail of `tests/rust_scene_worm.rs`, where OpenCascade
/// publishes the hidden root circle beyond the detail circle: every published
/// piece lies inside the 3 mm circle, and the role of every curve inside it
/// (the hidden root land and the visible bore and crest among them) matches
/// marching through the provider's predicate.
#[test]
fn axial_worm_detail_is_cropped_and_its_roles_marched() {
    let (s, c) = 0.1_f64.sin_cos();
    let solid = worm(1., 20., 2.)
        .transformed(RigidTransform {
            rotation: [[c, -s, 0.], [s, c, 0.], [0., 0., 1.]],
            translation_mm: [0.; 3],
        })
        .unwrap();
    let frame = PlaneFrame {
        origin: [0.3, -0.2, 0.4],
        u: [1., 0., 0.],
        v: [0., 1., 0.],
    };
    let (centre, radius) = ([0.5, 0.25], 3.);
    let vectors = crate::project_drawing(
        &[crate::DrawingSolid {
            solid: &solid,
            hatch: None,
        }],
        crate::DrawingProjection {
            frame,
            detail: Some((centre, radius)),
            hidden_lines: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(vectors.iter().any(|v| v.role == VectorRole::Hidden));
    assert!(vectors.iter().any(|v| v.role == VectorRole::Visible));
    for vector in &vectors {
        for p in &vector.preview_mm {
            let gap = (p[0] - centre[0]).hypot(p[1] - centre[1]);
            assert!(
                gap <= radius + 1e-6,
                "published {p:?} is {gap} mm from the centre"
            );
        }
    }
    let body = built(&solid, frame);
    let inner = &body.0;
    let (screw, _) = solid.rust_worm().unwrap();
    let own = |p: P3| screw.classify(inner.local(p), STRICT) == Containment::Inside;
    let top = body.bounds()[1][2];
    let mut proof = Proof::default();
    let mut root_hidden = 0;
    for candidate in body.candidates(0.0025).unwrap() {
        for i in 1..64 {
            let t = f64::from(i) / 64.;
            let p = candidate.path.point(t);
            if (p[0] - centre[0]).hypot(p[1] - centre[1]) > radius - 0.05 {
                continue;
            }
            let before = proof.roles;
            prove(
                candidate.kind,
                candidate.path.as_ref(),
                t,
                &own,
                &own,
                top,
                &vectors,
                &mut proof,
                "axial detail",
            );
            let local = inner.local(p);
            let on_root = (local[0].hypot(local[1]) - screw.root_radius_mm()).abs() < 1e-9;
            if on_root && proof.roles > before {
                root_hidden += 1;
            }
        }
    }
    println!(
        "axial detail creases {} roles {} unpublished {} root samples {root_hidden}",
        proof.creases, proof.roles, proof.unpublished
    );
    assert!(
        proof.roles >= 40 && root_hidden >= 5,
        "too few proved samples"
    );
}

/// `line` on arbitrary lines agrees with the provider's predicate inside and
/// between its intervals, and every interval end lies on the boundary.
#[test]
fn line_intervals_are_the_membership_predicate_on_arbitrary_lines() {
    let mut random = Random(23);
    for (name, solid) in solids() {
        for (label, frame) in views().into_iter().take(3) {
            let body = built(&solid, frame);
            let inner = &body.0;
            let [lo, hi] = body.bounds();
            let mut inside = 0;
            for _ in 0..150 {
                let start: P3 = std::array::from_fn(|k| lo[k] + (hi[k] - lo[k]) * random.next());
                let direction: P3 = std::array::from_fn(|_| 4. * random.next() - 2.);
                let spans = body.line(start, direction).unwrap().unwrap();
                for pair in spans.windows(2) {
                    assert!(pair[0][1] < pair[1][0], "{name} {label} {spans:?}");
                }
                let at = |s: f64| add(start, scaled(direction, s));
                for span in &spans {
                    for s in span {
                        assert_eq!(
                            classify(&solid, inner.local(at(*s)), 1e-9),
                            Containment::Boundary,
                            "{name} {label} {start:?} {direction:?} {s}"
                        );
                    }
                }
                let reach = 2. * norm(sub(hi, lo)) / norm(direction);
                for _ in 0..24 {
                    let s = -reach + 2. * reach * random.next();
                    let expected = match classify(&solid, inner.local(at(s)), 1e-6) {
                        Containment::Inside => true,
                        Containment::Outside => false,
                        Containment::Boundary => continue,
                    };
                    let found = spans.iter().any(|[a, b]| *a < s && s < *b);
                    assert_eq!(
                        found, expected,
                        "{name} {label} {start:?} {direction:?} {s}"
                    );
                    inside += usize::from(found);
                }
            }
            assert!(inside > 20, "{name} {label} {inside}");
        }
    }
    let body = built(&solids()[0].1, views()[0].1);
    for (start, direction) in [([0.; 3], [0.; 3]), ([f64::NAN, 0., 0.], [1., 0., 0.])] {
        let error = body.line(start, direction).unwrap().unwrap_err();
        assert_eq!(error.kind, GeometryErrorKind::Input);
    }
}

/// Largest distance from the preview points of `edges` in `role` to the
/// nearest preview chord of `reference` in `roles`.
fn reach(
    edges: &[ProjectedVector],
    role: VectorRole,
    reference: &[ProjectedVector],
    roles: &[VectorRole],
) -> f64 {
    let chords: Vec<[P2; 2]> = reference
        .iter()
        .filter(|edge| roles.contains(&edge.role))
        .flat_map(|edge| edge.preview_mm.windows(2).map(|w| [w[0], w[1]]))
        .collect();
    edges
        .iter()
        .filter(|edge| edge.role == role)
        .flat_map(|edge| edge.preview_mm.iter())
        .map(|&p| {
            chords
                .iter()
                .map(|&[a, b]| {
                    let d = [b[0] - a[0], b[1] - a[1]];
                    let w = [p[0] - a[0], p[1] - a[1]];
                    let t = ((w[0] * d[0] + w[1] * d[1]) / (d[0] * d[0] + d[1] * d[1]).max(1e-300))
                        .clamp(0., 1.);
                    (w[0] - t * d[0]).hypot(w[1] - t * d[1])
                })
                .fold(f64::INFINITY, f64::min)
        })
        .fold(0., f64::max)
}

fn drawn_length(edges: &[ProjectedVector], role: VectorRole) -> f64 {
    edges
        .iter()
        .filter(|edge| edge.role == role)
        .flat_map(|edge| edge.preview_mm.windows(2))
        .map(|w| (w[1][0] - w[0][0]).hypot(w[1][1] - w[0][1]))
        .sum()
}

/// Seen along the axis every helical rail of the thread is an arc of its
/// circle, turn after turn. The scene engine then compares the coincident
/// turns as arcs of one conic, and draws the worm as the closed-form axial
/// path does: the same curves in both roles, each published once.
#[test]
fn axial_worm_rails_are_arcs_and_the_scene_matches_the_axial_law() {
    for (name, solid) in solids().into_iter().take(2) {
        for (label, frame) in views().into_iter().take(2) {
            let body = built(&solid, frame);
            let mut rails = 0;
            for candidate in body.candidates(0.0005).unwrap() {
                let Some(Curve::Ellipse(arc)) = candidate.path.paper([0., 1.]) else {
                    continue;
                };
                // Every chord stop of the exact path lies on the paper
                // circle, and the arc ends are the path's ends.
                for t in candidate.path.flatten(0.0005) {
                    let p = candidate.path.point(t);
                    let r = (p[0] - arc.center[0]).hypot(p[1] - arc.center[1]);
                    assert!((r - arc.major).abs() < 1e-9, "{name} {label} {r}");
                }
                let ends = [0., 1.].map(|t| {
                    let p = candidate.path.point(t);
                    [p[0], p[1]]
                });
                let trim = arc.trim.map(|angle| arc.point(angle));
                let apart = |a: P2, b: P2| (a[0] - b[0]).hypot(a[1] - b[1]);
                assert!(
                    (apart(ends[0], trim[0]) < 1e-9 && apart(ends[1], trim[1]) < 1e-9)
                        || (apart(ends[0], trim[1]) < 1e-9 && apart(ends[1], trim[0]) < 1e-9),
                    "{name} {label} {ends:?} {trim:?}"
                );
                rails += 1;
            }
            // Two rails at each of the root and crest radii in half turns,
            // with the arcs of both end planes.
            assert!(rails > 50, "{name} {label} {rails}");
            let options = crate::DrawingProjection {
                frame,
                hidden_lines: true,
                chord_tolerance_mm: 0.002,
                ..Default::default()
            };
            let input = [crate::DrawingSolid {
                solid: &solid,
                hatch: None,
            }];
            let axial = crate::project_drawing(&input, options).unwrap();
            let start = std::time::Instant::now();
            let scene = super::super::project_solids(&input, options)
                .unwrap()
                .unwrap();
            let elapsed = start.elapsed();
            println!(
                "{{\"case\":\"axial {name} {label} through the scene\",\"vectors\":{},\"axial_vectors\":{},\"ms\":{}}}",
                scene.len(),
                axial.len(),
                elapsed.as_millis()
            );
            assert!(scene.len() < 2 * axial.len(), "{name} {label}");
            for role in [VectorRole::Visible, VectorRole::Hidden] {
                let forward = reach(&scene, role, &axial, &[role]);
                let backward = reach(&axial, role, &scene, &[role]);
                assert!(
                    forward < 0.004 && backward < 0.004,
                    "{name} {label} {role:?} {forward} {backward}"
                );
                // Coincident turns publish once.
                let [a, b] = [&scene, &axial].map(|edges| drawn_length(edges, role));
                assert!((a - b).abs() < 1e-3 * b, "{name} {label} {role:?} {a} {b}");
            }
        }
    }
}

/// An axial line on the waist of a throated root touches the boundary at
/// the mid plane from inside the material: the line is material on both
/// sides of that contact, wherever the contact falls in its interval.
#[test]
fn axial_lines_on_the_root_waist_are_material_through_the_contact() {
    for (name, solid) in solids().into_iter().skip(2) {
        let frame = views()[0].1;
        let body = built(&solid, frame);
        let inner = &body.0;
        let [waist, _] = inner.shell.annulus();
        let half = inner.shell.half();
        let mut inside = 0;
        for k in 0..64 {
            let angle = TAU * f64::from(k) / 64.;
            let at = |z: f64| inner.view([waist * angle.cos(), waist * angle.sin(), z]);
            let spans = body.sight([at(0.)[0], at(0.)[1]]).unwrap();
            for z in [-0.5 * half, 0.5 * half] {
                let depth = at(z)[2];
                let expected = match classify(&solid, inner.local(at(z)), 1e-9) {
                    Containment::Inside => true,
                    Containment::Outside => false,
                    Containment::Boundary => continue,
                };
                let found = spans.iter().any(|[a, b]| *a < depth && depth < *b);
                assert_eq!(found, expected, "{name} {angle} {z} {spans:?}");
                inside += usize::from(found);
            }
        }
        assert!(inside > 32, "{name} {inside}");
    }
}
