//! The bevel member against the provider's exact membership predicate and
//! independently evaluated surface normals.
use super::*;

mod drawn;
use crate::{GeometryPolicy, RigidTransform, SphericalBevelGear};
use multicad_core::document::EntityId;

fn specification(module: f64, teeth: u32, mate: u32, width: f64, bore: f64) -> SphericalBevelGear {
    let (z, mate) = (f64::from(teeth), f64::from(mate));
    let alpha = 20_f64.to_radians();
    let delta = 1_f64.atan2(mate / z);
    let distance = module * z / (2. * delta.sin());
    let beta = (delta.sin() * alpha.cos()).asin();
    let inv = |g: f64| {
        let sigma = (g.cos() / beta.cos()).acos();
        sigma / beta.sin() - sigma.sin().atan2(beta.sin() * sigma.cos())
    };
    SphericalBevelGear {
        teeth,
        pitch_cone_angle_rad: delta,
        base_cone_angle_rad: beta,
        tip_cone_angle_rad: delta + (module / distance).atan(),
        root_cone_angle_rad: delta - (1.25 * module / distance).atan(),
        outer_cone_distance_mm: distance,
        face_width_mm: width,
        bore_radius_mm: bore,
        half_base_tooth_angle_rad: (PI * module / 2. - 0.05) / (module * z) + inv(delta),
        root_fillet_radius_mm: module * 0.2,
    }
}

const PLACED: RigidTransform = RigidTransform {
    rotation: [[0.8, 0., 0.6], [0., 1., 0.], [-0.6, 0., 0.8]],
    translation_mm: [10., -20., 30.],
};

fn gears() -> Vec<(&'static str, Solid)> {
    let make = |spec, placed: bool| {
        let solid = Solid::rust_spherical_bevel_gear_named(
            spec,
            EntityId::new(),
            GeometryPolicy::default(),
        )
        .unwrap();
        if placed {
            solid.transformed(PLACED).unwrap()
        } else {
            solid
        }
    };
    vec![
        (
            "20/30 bored",
            make(specification(2., 20, 30, 8., 3.), false),
        ),
        (
            "16/16 solid placed",
            make(specification(2., 16, 16, 6., 0.), true),
        ),
        (
            "12/36 bored placed",
            make(specification(1.5, 12, 36, 5., 2.), true),
        ),
        (
            "40/20 solid",
            make(specification(1., 40, 20, 6., 0.), false),
        ),
    ]
}

/// Views of a gear: along its axis from both sides, across it, and oblique.
fn views(solid: &Solid) -> Vec<(&'static str, Gear)> {
    let (body, frame) = solid.bevel_body().unwrap();
    let [x, y, z] = [frame.u, frame.v, frame.normal()];
    let mix = |a: f64, b: f64, c: f64| -> [f64; 3] {
        std::array::from_fn(|i| a * x[i] + b * y[i] + c * z[i])
    };
    let s = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    [
        ("axial", [x, y]),
        ("reversed", [y, x]),
        ("side", [x, z]),
        ("oblique", [mix(s, -s, 0.), mix(t, t, -2. * t)]),
        ("skew", [mix(0.36, 0.48, 0.8), mix(0.8, -0.6, 0.)]),
    ]
    .into_iter()
    .map(|(name, [u, v])| {
        let view = View::new(PlaneFrame {
            origin: [1., 2., 3.],
            u,
            v,
        });
        (
            name,
            Gear {
                core: Arc::new(Core::new(body.geometry.clone(), frame, &view)),
            },
        )
    })
    .collect()
}

struct Random(u64);

impl Random {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }

    fn between(&mut self, a: f64, b: f64) -> f64 {
        a + (b - a) * self.next()
    }
}

fn on_boundary(core: &Core, view_point: P3, band: f64) -> bool {
    core.bevel.classify(core.to_local(view_point), band) == Containment::Boundary
}

#[test]
fn derivative_bounds_certify_the_root_isolation() {
    for (name, solid) in gears() {
        let gear = &views(&solid)[4].1;
        let core = &gear.core;
        for piece in 0..core.pieces.len() {
            let shape = core.shapes[piece];
            let h = 1e-4;
            for i in 1..40 {
                let s = f64::from(i) / 40.;
                let [a, m, b] = [s - h, s, s + h].map(|s| core.direction(2, piece, s));
                let second = norm(std::array::from_fn(|k| a.0[k] - 2. * m.0[k] + b.0[k])) / (h * h);
                assert!(
                    norm(m.1) <= shape.speed * (1. + 1e-9),
                    "{name} speed {piece}"
                );
                assert!(
                    second <= shape.curvature * (1. + 1e-5) + 1e-6,
                    "{name} curvature {piece}"
                );
                let rate = norm(std::array::from_fn(|k| {
                    (b.0[k] - a.0[k]) / (2. * h) - m.1[k]
                }));
                assert!(rate <= 1e-6, "{name} direction rate {piece}: {rate}");
                for end in End::BOTH {
                    let [a, m, b] = [s - h, s, s + h].map(|s| core.rim(end, 2, piece, s));
                    let second =
                        norm(std::array::from_fn(|k| a.0[k] - 2. * m.0[k] + b.0[k])) / (h * h);
                    assert!(
                        second <= core.rim_bound(end, piece) * (1. + 1e-5) + 1e-5,
                        "{name} rim bound {piece}"
                    );
                    let rate = norm(std::array::from_fn(|k| {
                        (b.0[k] - a.0[k]) / (2. * h) - m.1[k]
                    }));
                    assert!(rate <= 1e-5, "{name} rim rate {piece}: {rate}");
                }
                if flank(core.pieces[piece]) && !meridian(core.pieces[piece]) {
                    let [a, m, b] = [s - h, s, s + h].map(|s| core.limb(2, piece, s));
                    assert!(
                        ((b.0 - a.0) / (2. * h) - m.1).abs() <= 1e-6,
                        "{name} limb rate {piece}: {} {}",
                        (b.0 - a.0) / (2. * h),
                        m.1
                    );
                    assert!(
                        ((a.0 - 2. * m.0 + b.0) / (h * h)).abs()
                            <= 2. * shape.turn * shape.turn * (1. + 1e-5) + 1e-6,
                        "{name} limb bound {piece}"
                    );
                }
            }
        }
    }
}

#[test]
fn sight_intervals_are_the_material_along_the_line() {
    let mut random = Random(7);
    for (name, solid) in gears() {
        for (label, gear) in views(&solid) {
            let core = &gear.core;
            let [lo, hi] = gear.bounds();
            let mut spans = 0;
            let mut several = 0;
            for _ in 0..1500 {
                let paper = [random.between(lo[0], hi[0]), random.between(lo[1], hi[1])];
                let intervals = gear.sight(paper).unwrap();
                spans += intervals.len();
                several += usize::from(intervals.len() > 1);
                for pair in intervals.windows(2) {
                    assert!(pair[0][1] < pair[1][0], "{name} {label}: order");
                }
                for &[a, b] in &intervals {
                    assert!(b > a);
                    for z in [a, b] {
                        assert!(
                            on_boundary(core, [paper[0], paper[1], z], 1e-9),
                            "{name} {label}: end {z} at {paper:?} is off the boundary"
                        );
                    }
                }
                for _ in 0..12 {
                    let z = random.between(lo[2] - 1., hi[2] + 1.);
                    if intervals.iter().flatten().any(|end| (end - z).abs() < 1e-7) {
                        continue;
                    }
                    let expected = core.inside(core.to_local([paper[0], paper[1], z]));
                    let reported = intervals.iter().any(|[a, b]| *a < z && z < *b);
                    assert_eq!(expected, reported, "{name} {label}: {paper:?} depth {z}");
                }
            }
            assert!(spans > 200, "{name} {label}: {spans} intervals");
            if label != "axial" && label != "reversed" {
                assert!(several > 0, "{name} {label}: no line through several teeth");
            }
        }
    }
}

#[test]
fn candidates_lie_on_the_boundary_and_outlines_graze() {
    for (name, solid) in gears() {
        for (label, gear) in views(&solid) {
            let core = &gear.core;
            let candidates = gear.candidates(0.01).unwrap();
            let mut outlines = 0;
            for candidate in &candidates {
                let stops = candidate.path.flatten(0.01);
                assert_eq!(stops[0], 0.);
                assert_eq!(*stops.last().unwrap(), 1.);
                for pair in stops.windows(2) {
                    assert!(pair[1] > pair[0]);
                    for t in [pair[0], 0.5 * (pair[0] + pair[1]), pair[1]] {
                        assert!(
                            on_boundary(core, candidate.path.point(t), 1e-9),
                            "{name} {label}: {:?} point off the boundary",
                            candidate.kind
                        );
                    }
                    // Chords stay within the requested tolerance.
                    let [a, m, b] = [pair[0], 0.5 * (pair[0] + pair[1]), pair[1]]
                        .map(|t| candidate.path.point(t));
                    let sag = norm(std::array::from_fn(|k| m[k] - 0.5 * (a[k] + b[k])));
                    assert!(sag <= 0.01, "{name} {label}: sagitta {sag}");
                }
                if let Some(Curve::Line([a, b])) = candidate.path.paper([0., 1.]) {
                    let [p, q] = [0., 1.].map(|t| candidate.path.point(t));
                    assert!((a[0] - p[0]).hypot(a[1] - p[1]) < 1e-12);
                    assert!((b[0] - q[0]).hypot(b[1] - q[1]) < 1e-12);
                }
                if candidate.kind == Kind::Outline {
                    outlines += 1;
                    // The line of sight grazes the surface there: a transversal
                    // crossing would end exactly one material interval at
                    // the point.
                    // Along the axis the meridian flank planes contain the
                    // sight itself, which then has no crossing to count.
                    let along_axis = label == "axial" || label == "reversed";
                    for t in [0.3, 0.5, 0.7].into_iter().filter(|_| !along_axis) {
                        let p = candidate.path.point(t);
                        let ends = gear
                            .sight([p[0], p[1]])
                            .unwrap()
                            .into_iter()
                            .flatten()
                            .filter(|z| (z - p[2]).abs() < 2e-5)
                            .count();
                        assert!(
                            ends % 2 == 0,
                            "{name} {label}: the sight crosses the surface at an outline"
                        );
                    }
                }
            }
            if label == "side" || label == "oblique" || label == "skew" {
                assert!(outlines >= 4, "{name} {label}: {outlines} outlines");
            }
        }
    }
}

/// Distance from a gear-frame point to the nearest outline candidate.
fn outline_distance(gear: &Gear, candidates: &[Candidate], local: P3) -> f64 {
    let target = gear.core.to_view(local);
    candidates
        .iter()
        .filter(|candidate| candidate.kind == Kind::Outline)
        .map(|candidate| {
            let [a, b] = [0., 1.].map(|t| candidate.path.point(t));
            let chord = sub(b, a);
            let along = (dot3(sub(target, a), chord) / dot3(chord, chord)).clamp(0., 1.);
            norm(sub(sub(target, a), scaled(chord, along)))
        })
        .fold(f64::INFINITY, f64::min)
}

#[test]
fn every_turn_of_a_face_normal_has_an_outline() {
    for (name, solid) in gears() {
        for (label, gear) in views(&solid) {
            let core = &gear.core;
            let bevel = &core.bevel;
            let d = core.sight;
            let candidates = gear.candidates(0.01).unwrap();
            let mut turns = 0;
            let mut check = |before: f64, after: f64, point: P3, slack: f64, face: &str| {
                if (before > 0.) != (after > 0.) && before.abs() > 1e-9 && after.abs() > 1e-9 {
                    turns += 1;
                    let distance = outline_distance(&gear, &candidates, point);
                    assert!(
                        distance <= slack,
                        "{name} {label}: {face} normal turns without an outline ({distance} mm)"
                    );
                }
            };
            let grid = 400;
            // Side faces: the normal of the cone over a boundary piece,
            // from the provider's samples by central differences.
            for sector in 0..core.teeth {
                let turn = |p: P3| core.turn(sector, p);
                for &piece in bevel.pieces() {
                    let facing = |i: usize| {
                        let s = (i as f64 + 0.5) / grid as f64;
                        let [a, b] = [s - 1e-5, s + 1e-5].map(|s| bevel.sample(piece, s).direction);
                        let u = bevel.sample(piece, s).direction;
                        let normal = cross(u, sub(b, a));
                        (dot3(d, turn(scaled(normal, 1. / norm(normal)))), turn(u))
                    };
                    for i in 0..grid - 1 {
                        let [(before, u), (after, v)] = [facing(i), facing(i + 1)];
                        let middle = 0.5 * (core.reach(End::Toe, u) + core.reach(End::Heel, u));
                        check(
                            before,
                            after,
                            scaled(u, middle),
                            1.5 * middle * norm(sub(v, u)),
                            "side",
                        );
                    }
                }
            }
            // Caps and bore by azimuth.
            let steps = 4000;
            let [sin, cos] = core.pitch;
            for end in End::BOTH {
                let facing = |i: usize| {
                    let phi = TAU * (i as f64 + 0.5) / steps as f64;
                    (dot3(d, [sin * phi.cos(), sin * phi.sin(), cos]), phi)
                };
                for i in 0..steps {
                    let [(before, phi), (after, _)] = [facing(i), facing(i + 1)];
                    // A cap point just above the bottom of the generator.
                    let bottom = core.bottom(end, phi);
                    let point = add(
                        bottom,
                        [cos * phi.cos() * 1e-3, cos * phi.sin() * 1e-3, -sin * 1e-3],
                    );
                    let radius = bottom[0].hypot(bottom[1]) + 1e-3;
                    check(
                        before,
                        after,
                        point,
                        1.5 * radius * TAU / steps as f64 + 1e-6,
                        "cap",
                    );
                }
            }
            if let Some(bore) = core.bore {
                let facing = |i: usize| {
                    let phi = TAU * (i as f64 + 0.5) / steps as f64;
                    (d[0] * phi.cos() + d[1] * phi.sin(), phi)
                };
                for i in 0..steps {
                    let [(before, phi), (after, _)] = [facing(i), facing(i + 1)];
                    let [a, b] = End::BOTH.map(|end| core.bottom(end, phi));
                    check(
                        before,
                        after,
                        scaled(add(a, b), 0.5),
                        1.5 * bore * TAU / steps as f64,
                        "bore",
                    );
                }
            }
            if label == "side" || label == "oblique" || label == "skew" {
                assert!(turns >= 4, "{name} {label}: {turns} turns");
            }
        }
    }
}

/// Section paths published as exact ellipse arcs.
static CONICS: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

/// Signed area of chained section loops from fine chords, checking closure,
/// the plane, the boundary and the material side on the way.
fn loop_area(name: &str, gear: &Gear, depth: f64, loops: &[Loop]) -> f64 {
    let core = &gear.core;
    let mut area = 0.;
    for paths in loops {
        let mut previous = paths.last().unwrap().point(1.);
        for path in paths {
            let first = path.point(0.);
            let gap = norm(sub(first, previous));
            assert!(gap <= 1e-9 * core.scale, "{name}: loop gap {gap}");
            let stops = path.flatten(1e-5);
            let points: Vec<P3> = stops.iter().map(|&t| path.point(t)).collect();
            for (index, point) in points.iter().enumerate() {
                assert!((point[2] - depth).abs() <= 1e-9, "{name}: off the plane");
                assert!(
                    on_boundary(core, *point, 1e-8),
                    "{name}: section point {index} of {} off the boundary",
                    points.len()
                );
            }
            for pair in points.windows(2) {
                area += 0.5 * (pair[0][0] * pair[1][1] - pair[1][0] * pair[0][1]);
            }
            // Certified chords: the curve stays within the tolerance of
            // every chord, for a coarse and a fine request.
            for tolerance in [2e-2, 1e-4] {
                let stops = path.flatten(tolerance);
                assert!(stops.len() < 20_000, "{name}: {} stops", stops.len());
                assert_eq!((stops[0], *stops.last().unwrap()), (0., 1.));
                for pair in stops.windows(2) {
                    assert!(pair[1] > pair[0]);
                    let [a, b] = [pair[0], pair[1]].map(|t| path.point(t));
                    let chord = sub(b, a);
                    for k in 1..8 {
                        let offset = sub(
                            path.point(pair[0] + (pair[1] - pair[0]) * f64::from(k) / 8.),
                            a,
                        );
                        let along = if dot3(chord, chord) > 0. {
                            (dot3(offset, chord) / dot3(chord, chord)).clamp(0., 1.)
                        } else {
                            0.
                        };
                        let off = norm(sub(offset, scaled(chord, along)));
                        assert!(
                            off <= tolerance * (1. + 1e-9) + 1e-12,
                            "{name}: chord sag {off} beyond {tolerance}"
                        );
                    }
                }
            }
            // An exact paper curve carries the whole path.
            match path.paper([0., 1.]) {
                Some(Curve::Ellipse(ellipse)) => {
                    let [first, last] = ellipse.trim.map(|angle| ellipse.point(angle));
                    let [p, q] = [0., 1.].map(|t| path.point(t));
                    let gap = |a: P2, b: P3| (a[0] - b[0]).hypot(a[1] - b[1]);
                    assert!(
                        (gap(first, p) + gap(last, q)).min(gap(first, q) + gap(last, p)) <= 1e-9,
                        "{name}: ellipse ends"
                    );
                    for k in 0..=32 {
                        let point = path.point(f64::from(k) / 32.);
                        let offset = [point[0] - ellipse.center[0], point[1] - ellipse.center[1]];
                        let x =
                            (offset[0] * ellipse.u[0] + offset[1] * ellipse.u[1]) / ellipse.major;
                        let y =
                            (offset[0] * ellipse.v[0] + offset[1] * ellipse.v[1]) / ellipse.minor;
                        assert!(
                            (x.hypot(y) - 1.).abs() * ellipse.minor <= 1e-8,
                            "{name}: point off its ellipse"
                        );
                        let angle = (y.atan2(x) - ellipse.trim[0]).rem_euclid(TAU);
                        let span = ellipse.trim[1] - ellipse.trim[0];
                        assert!(
                            angle <= span + 1e-7 || angle >= TAU - 1e-7,
                            "{name}: point outside the ellipse trim"
                        );
                    }
                    CONICS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                }
                Some(Curve::Line([a, b])) => {
                    for k in 0..=8 {
                        let point = path.point(f64::from(k) / 8.);
                        let chord = [b[0] - a[0], b[1] - a[1]];
                        let off = ((point[0] - a[0]) * chord[1] - (point[1] - a[1]) * chord[0])
                            .abs()
                            / chord[0].hypot(chord[1]).max(1e-300);
                        assert!(off <= 1e-9, "{name}: point off its line");
                    }
                }
                None => {}
            }
            let length: f64 = points.windows(2).map(|w| norm(sub(w[1], w[0]))).sum();
            if length > 0.05 {
                // Material on the left, as the viewer sees the loop.
                let [a, m, b] = [0.37 - 1e-4, 0.37, 0.37 + 1e-4].map(|t| path.point(t));
                let tangent = sub(b, a);
                let size = tangent[0].hypot(tangent[1]);
                let left = [-tangent[1] / size, tangent[0] / size, 0.];
                let probe =
                    |sign: f64| core.inside(core.to_local(add(m, scaled(left, sign * 1e-5))));
                assert!(
                    probe(1.) && !probe(-1.),
                    "{name}: material is not on the left"
                );
            }
            previous = path.point(1.);
        }
    }
    area
}

/// Area of the membership predicate over the plane `z = depth` by scan
/// lines with bisected transitions.
fn scanned_area(gear: &Gear, depth: f64) -> f64 {
    let core = &gear.core;
    let [lo, hi] = gear.bounds();
    let (rows, columns) = (700, 900);
    let dy = (hi[1] - lo[1]) / rows as f64;
    let dx = (hi[0] - lo[0] + 2e-3) / columns as f64;
    let inside = |x: f64, y: f64| core.inside(core.to_local([x, y, depth]));
    let mut area = 0.;
    for row in 0..rows {
        let y = lo[1] + dy * (row as f64 + 0.5);
        let mut previous = (lo[0] - 1e-3, false);
        for column in 1..=columns {
            let x = lo[0] - 1e-3 + dx * column as f64;
            let state = inside(x, y);
            if state != previous.1 {
                let (mut a, mut b) = (previous.0, x);
                for _ in 0..40 {
                    let mid = 0.5 * (a + b);
                    if inside(mid, y) == previous.1 {
                        a = mid;
                    } else {
                        b = mid;
                    }
                }
                area += if state { -a } else { a } * dy;
            }
            previous = (x, state);
        }
    }
    area
}

#[test]
fn sections_bound_the_material_in_the_plane() {
    for (name, solid) in gears() {
        for (label, gear) in views(&solid) {
            let core = &gear.core;
            let [lo, hi] = gear.bounds();
            let apex = core.origin[2];
            let mut cut = 0;
            for depth in [
                apex,
                apex + 1e-13,
                apex - 3e-7,
                apex + 0.37,
                lo[2] + 0.2 * (hi[2] - lo[2]),
                lo[2] + 0.55 * (hi[2] - lo[2]),
                lo[2] + 0.93 * (hi[2] - lo[2]),
            ] {
                let case = format!("{name} {label} depth {depth}");
                let loops = gear.section(depth, 1e-3).unwrap();
                let area = loop_area(&case, &gear, depth, &loops);
                let scanned = scanned_area(&gear, depth);
                assert!(
                    (area - scanned).abs() <= 2e-3 * scanned + 2e-2,
                    "{case}: loop area {area} against scanned {scanned}"
                );
                cut += usize::from(!loops.is_empty());
            }
            assert!(cut >= 3, "{name} {label}: only {cut} planes cut the gear");
            assert!(gear.section(hi[2] + 1., 1e-3).unwrap().is_empty());
        }
    }
    // Bore cuts and closed cuts of the cones are ellipse arcs.
    let conics = CONICS.load(std::sync::atomic::Ordering::Relaxed);
    assert!(conics >= 40, "{conics} elliptic section paths");
}

/// The engine hides a point when a sight interval extends in front of it.
/// OpenCascade draws some rear tooth edges of the oblique view as visible
/// (see `rust_scene_bevel.rs`); marching the membership predicate along the
/// sight from the midpoint of every edge decides each such case here.
#[test]
fn material_in_front_of_an_edge_is_found_by_marching_the_membership() {
    for (name, solid) in gears() {
        // Along the axis the sight through a meridian edge lies in its
        // flank plane, where membership is decided by rounding.
        for (label, gear) in views(&solid).into_iter().skip(2) {
            let core = &gear.core;
            let top = gear.bounds()[1][2] + 0.1;
            let (mut covered, mut open) = (0, 0);
            for candidate in gear.candidates(0.01).unwrap() {
                let p = candidate.path.point(0.5);
                let ahead = gear
                    .sight([p[0], p[1]])
                    .unwrap()
                    .into_iter()
                    .map(|[a, b]| b - a.max(p[2]))
                    .fold(0., f64::max);
                let mut marched = 0.;
                let mut z = p[2] + 5e-3;
                while z < top {
                    if core.inside(core.to_local([p[0], p[1], z])) {
                        marched += 2e-3;
                    }
                    z += 2e-3;
                }
                // Compare only clear cases: grazing lines sit on the surface.
                if ahead > 0.05 {
                    assert!(
                        marched > 0.02,
                        "{name} {label}: {ahead} mm reported, {marched} marched"
                    );
                    covered += 1;
                } else if ahead < 1e-6 {
                    assert!(
                        marched < 0.02,
                        "{name} {label}: {marched} mm marched, none reported"
                    );
                    open += 1;
                }
            }
            assert!(
                covered > 0 && open > 0,
                "{name} {label}: {covered} covered, {open} open"
            );
        }
    }
}

/// A nine-tooth gear with meridian flanks; no meridian plane carries a
/// second meridian face on the far side of the axis.
fn meridian_gear(placed: bool) -> Solid {
    let solid = Solid::rust_spherical_bevel_gear_named(
        specification(2., 9, 13, 8., 3.),
        EntityId::new(),
        GeometryPolicy::default(),
    )
    .unwrap();
    assert!(solid.bevel_body().unwrap().0.geometry.has_meridian());
    if placed {
        solid.transformed(PLACED).unwrap()
    } else {
        solid
    }
}

fn gear_in(solid: &Solid, frame: PlaneFrame) -> Gear {
    let (body, placement) = solid.bevel_body().unwrap();
    Gear {
        core: Arc::new(Core::new(
            body.geometry.clone(),
            placement,
            &View::new(frame),
        )),
    }
}

/// Strictly inside: the point and its neighbors turned about the gear axis
/// by a small angle either way are material. This open-set test does not
/// depend on the rounding of a point lying on a meridian face.
fn interior(core: &Core, local: P3) -> bool {
    [-1e-9, 0., 1e-9].into_iter().all(|turn: f64| {
        let (sin, cos) = turn.sin_cos();
        let p = [
            cos * local[0] - sin * local[1],
            sin * local[0] + cos * local[1],
            local[2],
        ];
        core.inside(p)
    })
}

#[test]
fn sight_lines_in_meridian_planes_are_decided_exactly() {
    for placed in [false, true] {
        let solid = meridian_gear(placed);
        let (_, frame) = solid.bevel_body().unwrap();
        for (label, [u, v]) in [
            ("axial", [frame.u, frame.v]),
            ("reversed", [frame.v, frame.u]),
        ] {
            let gear = gear_in(
                &solid,
                PlaneFrame {
                    origin: [1., 2., 3.],
                    u,
                    v,
                },
            );
            let core = &gear.core;
            let planes = core.meridians();
            assert_eq!(planes.len(), 18);
            let mut checked = 0;
            for &(sector, piece, ..) in planes.iter().step_by(5) {
                // Lines through the face, through its rims and through the
                // rays at both ends of the face.
                let mut points = Vec::new();
                for s in [0., 0.3, 0.7, 1.] {
                    let u = core.direction(sector, piece, s).0;
                    let [toe, heel] = [End::Toe, End::Heel].map(|end| core.reach(end, u));
                    for f in [0., 0.4, 1.] {
                        points.push(scaled(u, toe + (heel - toe) * f));
                    }
                }
                for local in points {
                    let target = core.to_view(local);
                    let reference = gear.sight([target[0], target[1]]).unwrap();
                    // Perturbing the paper point by rounding-sized amounts
                    // leaves the answer unchanged.
                    for (dx, dy) in [(1e-14, 0.), (0., -1e-14), (-3e-15, 2e-15)] {
                        let other = gear.sight([target[0] + dx, target[1] + dy]).unwrap();
                        assert_eq!(
                            other.len(),
                            reference.len(),
                            "{label}: {other:?} {reference:?}"
                        );
                        for (a, b) in other.iter().zip(&reference) {
                            assert!((a[0] - b[0]).abs() + (a[1] - b[1]).abs() <= 1e-9);
                        }
                    }
                    let [lo, hi] = gear.bounds();
                    let steps = 4000;
                    for k in 0..=steps {
                        let z =
                            lo[2] - 0.5 + (hi[2] - lo[2] + 1.) * f64::from(k) / f64::from(steps);
                        if reference.iter().flatten().any(|end| (end - z).abs() < 0.01) {
                            continue;
                        }
                        let expected = interior(core, core.to_local([target[0], target[1], z]));
                        let reported = reference.iter().any(|[a, b]| *a < z && z < *b);
                        assert_eq!(
                            expected, reported,
                            "{label} placed {placed}: depth {z} on {target:?}"
                        );
                        checked += 1;
                    }
                }
            }
            assert!(checked > 100_000, "{checked}");
        }
    }
}

/// Winding number of chained loops about a paper point.
fn winding(loops: &[Loop], point: P2) -> i32 {
    let mut angle = 0.;
    for paths in loops {
        for path in paths {
            let stops = path.flatten(1e-4);
            for pair in stops.windows(2) {
                let [a, b] = [pair[0], pair[1]].map(|t| {
                    let p = path.point(t);
                    (p[1] - point[1]).atan2(p[0] - point[0])
                });
                angle += (b - a + PI).rem_euclid(TAU) - PI;
            }
        }
    }
    (angle / TAU).round() as i32
}

#[test]
fn section_planes_containing_a_meridian_face_keep_the_material_behind() {
    let solid = meridian_gear(false);
    let gear = gear_in(&solid, PlaneFrame::default());
    let (sector, piece, radial, _) = gear.core.meridians()[0];
    let normal = [-radial[1], radial[0], 0.];
    let face = gear.core.bevel.side_area_mm2(gear.core.pieces[piece]);
    let mut behind_sides = 0;
    for flip in [false, true] {
        // The viewer looks along ∓ the meridian normal: u × v = ±normal.
        let [u, v] = if flip {
            [radial, [0., 0., 1.]]
        } else {
            [[0., 0., 1.], radial]
        };
        let gear = gear_in(
            &solid,
            PlaneFrame {
                origin: [1., 2., 3.],
                u,
                v,
            },
        );
        let core = &gear.core;
        let toward = core.sight;
        assert!((dot3(toward, normal).abs() - 1.).abs() <= 1e-12);
        let depth = core.origin[2];
        let label = format!("flip {flip}");
        let loops = gear.section(depth, 1e-3).unwrap();
        let area = loop_area(&label, &gear, depth, &loops);
        let [behind, front] = [-1e-7, 1e-7].map(|h| scanned_area(&gear, depth + h));
        println!(
            "{{\"case\":\"{label}\",\"loops\":{area},\"behind\":{behind},\"front\":{front},\"face\":{face}}}"
        );
        assert!(
            (area - behind).abs() <= 2e-3 * behind + 2e-2,
            "{label}: {area} against {behind} behind"
        );
        // The face is the whole difference between the two sides.
        assert!(
            ((behind - front).abs() - face).abs() <= 0.02 * face,
            "{label}"
        );
        // A point in the middle of the face lies in the section exactly
        // when the tooth is behind the plane.
        let u = core.direction(sector, piece, 0.5).0;
        let r = 0.5 * (core.reach(End::Toe, u) + core.reach(End::Heel, u));
        let middle = core.to_view(scaled(u, r));
        let tooth_behind = core.inside(core.to_local([middle[0], middle[1], depth - 1e-6]));
        behind_sides += usize::from(tooth_behind);
        assert_eq!(
            winding(&loops, [middle[0], middle[1]]),
            i32::from(tooth_behind),
            "{label}"
        );
        assert_eq!(behind > front, tooth_behind, "{label}");
    }
    assert_eq!(behind_sides, 1);
}

#[test]
fn line_intervals_are_the_material_along_any_line() {
    let mut random = Random(11);
    for (name, solid) in gears()
        .into_iter()
        .chain([("meridian", meridian_gear(true))])
    {
        for (label, gear) in views(&solid).into_iter().take(4).skip(2) {
            let core = &gear.core;
            let [lo, hi] = gear.bounds();
            let mut spans = 0;
            for _ in 0..600 {
                let start: P3 = std::array::from_fn(|k| random.between(lo[k] - 2., hi[k] + 2.));
                let mut direction: P3 = std::array::from_fn(|_| random.between(-1., 1.));
                let length = norm(direction);
                direction = scaled(direction, random.between(0.3, 3.) / length);
                let intervals = gear.line(start, direction).unwrap().unwrap();
                spans += intervals.len();
                for pair in intervals.windows(2) {
                    assert!(pair[0][1] < pair[1][0], "{name} {label}: order");
                }
                let at = |s: f64| add(start, scaled(direction, s));
                for &[a, b] in &intervals {
                    assert!(b > a);
                    for s in [a, b] {
                        assert!(
                            on_boundary(core, at(s), 1e-9),
                            "{name} {label}: end {s} off the boundary"
                        );
                    }
                }
                // The parameter range crossing the bounds box.
                let reach = 2. * norm(sub(hi, lo)) / norm(direction) + 10.;
                for _ in 0..16 {
                    let s = random.between(-reach, reach);
                    if intervals.iter().flatten().any(|end| (end - s).abs() < 1e-7) {
                        continue;
                    }
                    let expected = core.inside(core.to_local(at(s)));
                    let reported = intervals.iter().any(|[a, b]| *a < s && s < *b);
                    assert_eq!(
                        expected, reported,
                        "{name} {label}: {start:?} {direction:?} at {s}"
                    );
                }
                // A start moved far along the line keeps the intervals.
                let shift = 1e3;
                let moved = gear.line(at(-shift), direction).unwrap().unwrap();
                assert_eq!(moved.len(), intervals.len(), "{name} {label}: moved start");
                for (a, b) in moved.iter().zip(&intervals) {
                    assert!(
                        (a[0] - shift - b[0]).abs() + (a[1] - shift - b[1]).abs() <= 1e-9 * shift
                    );
                }
            }
            assert!(spans > 100, "{name} {label}: {spans} intervals");
            // Along the view normal the line is the line of sight.
            for _ in 0..200 {
                let paper = [random.between(lo[0], hi[0]), random.between(lo[1], hi[1])];
                let sight = gear.sight(paper).unwrap();
                let line = gear
                    .line([paper[0], paper[1], 0.], [0., 0., 1.])
                    .unwrap()
                    .unwrap();
                assert_eq!(sight.len(), line.len());
                for (a, b) in sight.iter().zip(&line) {
                    assert!((a[0] - b[0]).abs() + (a[1] - b[1]).abs() <= 1e-9);
                }
            }
            assert!(gear.line([0.; 3], [0.; 3]).unwrap().is_err());
        }
    }
}
