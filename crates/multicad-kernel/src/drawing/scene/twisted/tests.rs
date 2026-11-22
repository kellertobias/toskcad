//! The twisted member against an independent exact membership predicate.
//!
//! The oracle never touches the member's code: it rebuilds the transverse
//! profile from the gear specification with the pure provider and decides
//! membership by the winding number of the exact profile curves, refined
//! until every chord provably subtends the same angle as its curve.
use super::super::*;
use super::Gear;
use crate::{GearHub, GearSweep, InvoluteGear, RigidTransform};
use multicad_core::document::EntityId;
use multicad_rust_kernel::spur_gear::{Edge, Profile, Spec};
use std::f64::consts::PI;

pub(super) fn specification(sweep: GearSweep, twist: f64, bore: f64, hub: bool) -> InvoluteGear {
    let a = 20_f64.to_radians();
    let pitch = 12.;
    InvoluteGear {
        hub: hub.then_some(GearHub {
            radius_mm: 6.,
            axial_start_mm: -2.,
            axial_end_mm: 5.5,
        }),
        teeth: 24,
        pitch_radius_mm: pitch,
        base_radius_mm: pitch * a.cos(),
        tip_radius_mm: pitch + 1.,
        root_radius_mm: pitch - 1.25,
        root_fillet_radius_mm: 0.2,
        half_base_angle_rad: (PI / 2. - 0.05) / (2. * pitch) + a.tan() - a,
        face_width_mm: 4.,
        bore_radius_mm: bore,
        twist_rad_per_mm: twist,
        sweep,
    }
}

pub(super) fn pose() -> RigidTransform {
    RigidTransform {
        rotation: [[0.8, -0.36, 0.48], [0.6, 0.48, -0.64], [0., 0.8, 0.6]],
        translation_mm: [3., -2., 1.5],
    }
}

pub(super) fn solid(spec: InvoluteGear, placed: bool) -> Solid {
    let build = match spec.sweep {
        GearSweep::Helical => Solid::rust_helical_gear_named,
        GearSweep::GappedChevron { .. } => Solid::rust_double_helical_named,
        _ => Solid::rust_herringbone_named,
    };
    let solid = build(spec, EntityId::new(), Default::default()).unwrap();
    if placed {
        solid.transformed(pose()).unwrap()
    } else {
        solid
    }
}

pub(super) fn frames() -> Vec<(&'static str, PlaneFrame)> {
    let unit = |v: [f64; 3]| {
        let l = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        v.map(|c| c / l)
    };
    let u = unit([2., 1., -0.5]);
    let w = unit([0.3, 1., 2.]);
    let dot = u[0] * w[0] + u[1] * w[1] + u[2] * w[2];
    let v = unit(std::array::from_fn(|i| w[i] - dot * u[i]));
    vec![
        ("axial", PlaneFrame::default()),
        (
            "reversed",
            PlaneFrame {
                origin: [0.; 3],
                u: [1., 0., 0.],
                v: [0., -1., 0.],
            },
        ),
        (
            "side",
            PlaneFrame {
                origin: [0.; 3],
                u: [1., 0., 0.],
                v: [0., 0., 1.],
            },
        ),
        (
            "oblique",
            PlaneFrame {
                origin: [1., 2., 3.],
                u,
                v,
            },
        ),
    ]
}

pub(super) struct Oracle {
    edges: Vec<Edge>,
    spec: InvoluteGear,
    placement: PlaneFrame,
    view: View,
}

fn swept(edge: Edge, w: P2, [a, b]: [f64; 2], [pa, pb]: [P2; 2], depth: u32) -> f64 {
    let chord = (pb[0] - pa[0]).hypot(pb[1] - pa[1]);
    let mid = [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1])];
    if (w[0] - mid[0]).hypot(w[1] - mid[1]) > chord || depth > 60 {
        let [p, q] = [[pa[0] - w[0], pa[1] - w[1]], [pb[0] - w[0], pb[1] - w[1]]];
        return (p[0] * q[1] - p[1] * q[0]).atan2(p[0] * q[0] + p[1] * q[1]);
    }
    let m = 0.5 * (a + b);
    let pm = edge.point(m);
    swept(edge, w, [a, m], [pa, pm], depth + 1) + swept(edge, w, [m, b], [pm, pb], depth + 1)
}

impl Oracle {
    pub(super) fn new(spec: InvoluteGear, placement: PlaneFrame, view: View) -> Self {
        let profile = Profile::new(
            Spec {
                teeth: spec.teeth,
                pitch_radius_mm: spec.pitch_radius_mm,
                base_radius_mm: spec.base_radius_mm,
                tip_radius_mm: spec.tip_radius_mm,
                root_radius_mm: spec.root_radius_mm,
                root_fillet_radius_mm: spec.root_fillet_radius_mm,
                half_base_angle_rad: spec.half_base_angle_rad,
                face_width_mm: spec.face_width_mm,
                bore_radius_mm: spec.bore_radius_mm,
            },
            1e-6,
            1e6,
        )
        .unwrap();
        Self {
            edges: profile.outline().iter().map(|edge| edge.edge).collect(),
            spec,
            placement,
            view,
        }
    }

    fn gap(&self) -> f64 {
        match self.spec.sweep {
            GearSweep::GappedChevron { gap_mm } => gap_mm,
            _ => 0.,
        }
    }

    /// Length of each tooth band of a chevron.
    fn band(&self) -> f64 {
        0.5 * (self.spec.face_width_mm - self.gap())
    }

    /// Inside the central gap of a double-helical gear.
    fn in_gap(&self, z: f64) -> bool {
        self.gap() > 0. && z > self.band() && z < self.spec.face_width_mm - self.band()
    }

    pub(super) fn phase(&self, z: f64) -> f64 {
        let h = self.spec.face_width_mm;
        self.spec.twist_rad_per_mm
            * if self.spec.sweep == GearSweep::Helical {
                z
            } else {
                z.min(h - z).min(self.band())
            }
    }

    pub(super) fn in_profile(&self, w: P2) -> bool {
        let r = w[0].hypot(w[1]);
        if r < self.spec.root_radius_mm - 1e-3 {
            return true;
        }
        if r > self.spec.tip_radius_mm + 1e-3 {
            return false;
        }
        let mut total = 0.;
        for &edge in &self.edges {
            let mut previous = edge.point(0.);
            for i in 1..=8 {
                let [a, b] = [f64::from(i - 1) / 8., f64::from(i) / 8.];
                let next = edge.point(b);
                total += swept(edge, w, [a, b], [previous, next], 0);
                previous = next;
            }
        }
        total.abs() > PI
    }

    pub(super) fn local(&self, p: P3) -> P3 {
        let world: P3 = std::array::from_fn(|i| {
            self.view.origin[i]
                + p[0] * self.view.axes[0][i]
                + p[1] * self.view.axes[1][i]
                + p[2] * self.view.axes[2][i]
        });
        let offset: P3 = std::array::from_fn(|i| world[i] - self.placement.origin[i]);
        [self.placement.u, self.placement.v, self.placement.normal()]
            .map(|axis| axis[0] * offset[0] + axis[1] * offset[1] + axis[2] * offset[2])
    }

    pub(super) fn inside_local(&self, [x, y, z]: P3) -> bool {
        let r = x.hypot(y);
        if r <= self.spec.bore_radius_mm {
            return false;
        }
        if z > 0. && z < self.spec.face_width_mm {
            // The central gap of a double-helical gear keeps the root core.
            if self.in_gap(z) {
                return r < self.spec.root_radius_mm;
            }
            let (sin, cos) = (-self.phase(z)).sin_cos();
            return self.in_profile([x * cos - y * sin, x * sin + y * cos]);
        }
        self.spec
            .hub
            .is_some_and(|hub| r < hub.radius_mm && z > hub.axial_start_mm && z < hub.axial_end_mm)
    }

    pub(super) fn inside(&self, p: P3) -> bool {
        self.inside_local(self.local(p))
    }
}

pub(super) struct Random(u64);

impl Random {
    pub(super) fn unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 11) as f64 / (1_u64 << 53) as f64
    }
}

pub(super) struct Case {
    pub label: String,
    pub gear: Gear,
    pub oracle: Oracle,
}

pub(super) fn cases(views: &[&str]) -> Vec<Case> {
    let mut output = Vec::new();
    for (sweep, name) in [
        (GearSweep::Helical, "helical"),
        (GearSweep::Chevron, "herringbone"),
        (GearSweep::GappedChevron { gap_mm: 1. }, "double-helical"),
    ] {
        for (twist, bore, hub, placed) in [
            (0.05, 3., false, false),
            (-0.05, 0., false, true),
            (0.05, 3., true, true),
            (-0.03, 0., true, false),
        ] {
            let spec = specification(sweep, twist, bore, hub);
            let solid = solid(spec, placed);
            let (_, placement) = solid.herringbone_body().unwrap();
            for (view_name, frame) in frames() {
                if !views.contains(&view_name) {
                    continue;
                }
                let view = View::new(frame);
                output.push(Case {
                    label: format!(
                        "{name} twist {twist} bore {bore} hub {hub} placed {placed} {view_name}"
                    ),
                    gear: Gear::new(&solid, &view).unwrap().unwrap(),
                    oracle: Oracle::new(spec, placement, view),
                });
            }
        }
    }
    output
}

const ALL: [&str; 4] = ["axial", "reversed", "side", "oblique"];

#[test]
fn sight_intervals_match_exact_membership_and_end_on_the_surface() {
    let mut random = Random(7);
    let mut worst: f64 = 0.;
    let mut ends = 0;
    for case in cases(&ALL) {
        let [lo, hi] = case.gear.bounds();
        let scale = lo.iter().chain(&hi).fold(1_f64, |s, v| s.max(v.abs()));
        for _ in 0..60 {
            let paper = [0, 1].map(|k| lo[k] + (hi[k] - lo[k]) * random.unit());
            let spans = case.gear.sight(paper).unwrap();
            assert!(
                spans.windows(2).all(|pair| pair[0][1] < pair[1][0])
                    && spans.iter().all(|[a, b]| a < b),
                "{}: {spans:?}",
                case.label
            );
            for _ in 0..24 {
                let depth = lo[2] - 1. + (hi[2] - lo[2] + 2.) * random.unit();
                if spans.iter().flatten().any(|end| (end - depth).abs() < 1e-6) {
                    continue;
                }
                assert_eq!(
                    spans.iter().any(|[a, b]| *a < depth && depth < *b),
                    case.oracle.inside([paper[0], paper[1], depth]),
                    "{} at {paper:?} depth {depth}: {spans:?}",
                    case.label
                );
            }
            // Each end is the exact boundary: bisect the oracle around it.
            let stops: Vec<f64> = spans.iter().flatten().copied().collect();
            for (index, &end) in stops.iter().enumerate() {
                let h = 1e-5;
                if stops
                    .iter()
                    .enumerate()
                    .any(|(other, stop)| other != index && (stop - end).abs() < 4. * h)
                {
                    continue;
                }
                let at = |depth: f64| case.oracle.inside([paper[0], paper[1], depth]);
                let [mut a, mut b] = [end - h, end + h];
                let first = at(a);
                if first == at(b) {
                    // Another boundary passes within the bracket.
                    continue;
                }
                for _ in 0..60 {
                    let mid = 0.5 * (a + b);
                    if at(mid) == first {
                        a = mid;
                    } else {
                        b = mid;
                    }
                }
                let residual = (0.5 * (a + b) - end).abs() / scale;
                worst = worst.max(residual);
                ends += 1;
                assert!(
                    residual <= 1e-9,
                    "{} at {paper:?}: end {end} is {residual} from the surface",
                    case.label
                );
            }
        }
    }
    assert!(ends > 2000, "only {ends} interval ends were certified");
    println!("{{\"sight_ends\":{ends},\"worst_relative_residual\":{worst:e}}}");
}

fn directions() -> Vec<P3> {
    let mut output = Vec::new();
    for x in [-1., 0., 1.] {
        for y in [-1., 0., 1.] {
            for z in [-1., 0., 1.] {
                if x != 0. || y != 0. || z != 0. {
                    output.push([x, y, z]);
                }
            }
        }
    }
    output
}

/// Material and void both occur within `reach` of a boundary point.
fn on_boundary(oracle: &Oracle, p: P3, reach: f64) -> bool {
    let sides: Vec<bool> = directions()
        .into_iter()
        .map(|d| oracle.inside(std::array::from_fn(|i| p[i] + reach * d[i])))
        .collect();
    sides.contains(&true) && sides.contains(&false)
}

impl Oracle {
    fn twist(&self, z: f64) -> f64 {
        let k = self.spec.twist_rad_per_mm;
        if self.spec.sweep == GearSweep::Helical || z < 0.5 * self.spec.face_width_mm {
            k
        } else {
            -k
        }
    }

    /// Outward flank normal in gear coordinates at an edge fraction and level.
    fn flank_normal(&self, edge: Edge, s: f64, z: f64) -> P3 {
        let [tx, ty] = edge.tangent(s);
        let c = edge.point(s);
        let (sin, cos) = self.phase(z).sin_cos();
        let [nx, ny] = [ty, -tx];
        let n = [
            nx * cos - ny * sin,
            nx * sin + ny * cos,
            self.twist(z) * (c[0] * tx + c[1] * ty),
        ];
        let l = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        n.map(|v| v / l)
    }

    fn flank_point(&self, edge: Edge, s: f64, z: f64) -> P3 {
        let [x, y] = edge.point(s);
        let (sin, cos) = self.phase(z).sin_cos();
        [x * cos - y * sin, x * sin + y * cos, z]
    }

    /// The surface normal at a boundary point away from the end faces,
    /// found on the nearest exact profile curve.
    fn normal(&self, [x, y, z]: P3) -> P3 {
        let r = x.hypot(y);
        let radial = [x / r, y / r, 0.];
        if !(z > 0. && z < self.spec.face_width_mm)
            || (r - self.spec.bore_radius_mm).abs() < 1e-6
            || self.in_gap(z)
        {
            return radial;
        }
        let (sin, cos) = (-self.phase(z)).sin_cos();
        let w = [x * cos - y * sin, x * sin + y * cos];
        let distance = |edge: Edge, s: f64| {
            let p = edge.point(s);
            (p[0] - w[0]).hypot(p[1] - w[1])
        };
        let mut nearest: Vec<(f64, Edge, f64)> = self
            .edges
            .iter()
            .map(|&edge| {
                (0..=32)
                    .map(|i| f64::from(i) / 32.)
                    .map(|s| (distance(edge, s), edge, s))
                    .min_by(|a, b| a.0.total_cmp(&b.0))
                    .unwrap()
            })
            .collect();
        nearest.sort_by(|a, b| a.0.total_cmp(&b.0));
        let (found, edge, s) = nearest[..4]
            .iter()
            .map(|&(_, edge, s)| {
                let [mut a, mut b] = [(s - 1. / 32.).max(0.), (s + 1. / 32.).min(1.)];
                for _ in 0..200 {
                    let [m1, m2] = [a + (b - a) / 3., b - (b - a) / 3.];
                    if distance(edge, m1) < distance(edge, m2) {
                        b = m2;
                    } else {
                        a = m1;
                    }
                }
                let s = 0.5 * (a + b);
                (distance(edge, s), edge, s)
            })
            .min_by(|a, b| a.0.total_cmp(&b.0))
            .unwrap();
        assert!(
            found < 1e-6,
            "point [{x}, {y}, {z}] is {found} off every flank; nearest {edge:?} at {s}"
        );
        self.flank_normal(edge, s, z)
    }

    fn toward(&self) -> P3 {
        let [a, b] = [self.local([0.; 3]), self.local([0., 0., 1.])];
        std::array::from_fn(|i| b[i] - a[i])
    }
}

fn segment_distance(p: P3, [a, b]: [P3; 2]) -> f64 {
    let d: P3 = std::array::from_fn(|i| b[i] - a[i]);
    let o: P3 = std::array::from_fn(|i| p[i] - a[i]);
    let dd = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    let t = if dd > 0. {
        ((o[0] * d[0] + o[1] * d[1] + o[2] * d[2]) / dd).clamp(0., 1.)
    } else {
        0.
    };
    (0..3)
        .map(|i| (o[i] - t * d[i]).powi(2))
        .sum::<f64>()
        .sqrt()
}

fn chords_of(path: &dyn Path, tolerance: f64) -> Vec<[P3; 2]> {
    let ts = path.flatten(tolerance);
    assert!(ts[0] == 0. && *ts.last().unwrap() == 1. && ts.windows(2).all(|p| p[0] < p[1]));
    ts.windows(2)
        .map(|pair| [path.point(pair[0]), path.point(pair[1])])
        .collect()
}

#[test]
fn candidates_lie_on_the_boundary_and_outlines_are_folds() {
    let mut worst: f64 = 0.;
    let mut outlines = 0;
    // One hand and hub arrangement per view keeps the oracle affordable.
    for (index, case) in cases(&ALL).into_iter().enumerate() {
        if index % 5 != 0 {
            continue;
        }
        let toward = case.oracle.toward();
        for candidate in case.gear.candidates(1e-3).unwrap() {
            // Chords stay within the tolerance of the exact path.
            let ts = candidate.path.flatten(1e-3);
            for pair in ts.windows(2) {
                let [a, b] = [candidate.path.point(pair[0]), candidate.path.point(pair[1])];
                for f in [0.25, 0.5, 0.75] {
                    let p = candidate.path.point(pair[0] + (pair[1] - pair[0]) * f);
                    assert!(segment_distance(p, [a, b]) <= 1e-3, "{}", case.label);
                }
            }
            for t in [0.07, 0.5, 0.93] {
                let p = candidate.path.point(t);
                assert!(
                    on_boundary(&case.oracle, p, 1e-7),
                    "{}: {:?} candidate at {p:?} is off the boundary",
                    case.label,
                    candidate.kind
                );
                if candidate.kind == Kind::Outline {
                    let n = case.oracle.normal(case.oracle.local(p));
                    let residual = (n[0] * toward[0] + n[1] * toward[1] + n[2] * toward[2]).abs();
                    worst = worst.max(residual);
                    outlines += 1;
                    assert!(
                        residual <= 1e-6,
                        "{}: outline normal is {residual} off perpendicular at {p:?}",
                        case.label
                    );
                }
            }
        }
    }
    assert!(outlines > 300);
    println!("{{\"outline_points\":{outlines},\"worst_normal_residual\":{worst:e}}}");
}

#[test]
fn outlines_are_complete_on_every_flank() {
    let mut found = 0;
    let mut worst: f64 = 0.;
    for case in cases(&["side", "oblique"]) {
        let toward = case.oracle.toward();
        let outlines: Vec<[P3; 2]> = case
            .gear
            .candidates(1e-4)
            .unwrap()
            .into_iter()
            .filter(|candidate| candidate.kind == Kind::Outline)
            .flat_map(|candidate| chords_of(candidate.path.as_ref(), 1e-4))
            .collect();
        let h = case.oracle.spec.face_width_mm;
        let bands: Vec<[f64; 2]> = if case.oracle.spec.sweep == GearSweep::Helical {
            vec![[0., h]]
        } else {
            let band = case.oracle.band();
            vec![[0., band], [h - band, h]]
        };
        const N: usize = 20;
        for &edge in &case.oracle.edges {
            for band in &bands {
                // Stay strictly inside the band so that its twist applies.
                let z_at = |j: usize| {
                    band[0] + (band[1] - band[0]) * (0.001 + 0.998 * j as f64 / N as f64)
                };
                let s_at = |i: usize| 0.001 + 0.998 * i as f64 / N as f64;
                let value = |s: f64, z: f64| {
                    let n = case.oracle.flank_normal(edge, s, z);
                    n[0] * toward[0] + n[1] * toward[1] + n[2] * toward[2]
                };
                let mut check = |from: [f64; 2], to: [f64; 2]| {
                    let [fa, fb] = [value(from[0], from[1]), value(to[0], to[1])];
                    if fa.abs() < 1e-4 || fb.abs() < 1e-4 || (fa < 0.) == (fb < 0.) {
                        return;
                    }
                    let [mut a, mut b] = [0_f64, 1.];
                    let at = |t: f64| {
                        [
                            from[0] + (to[0] - from[0]) * t,
                            from[1] + (to[1] - from[1]) * t,
                        ]
                    };
                    for _ in 0..60 {
                        let mid = 0.5 * (a + b);
                        let [s, z] = at(mid);
                        if (value(s, z) < 0.) == (fa < 0.) {
                            a = mid;
                        } else {
                            b = mid;
                        }
                    }
                    let [s, z] = at(0.5 * (a + b));
                    let local = case.oracle.flank_point(edge, s, z);
                    let p = case.gear.frame.place(local);
                    let distance = outlines
                        .iter()
                        .map(|chord| segment_distance(p, *chord))
                        .fold(f64::INFINITY, f64::min);
                    worst = worst.max(distance);
                    found += 1;
                    assert!(
                        distance <= 3e-4,
                        "{}: a fold at {p:?} is {distance} mm from every outline",
                        case.label
                    );
                };
                for i in 0..=N {
                    for j in 0..=N {
                        if i < N {
                            check([s_at(i), z_at(j)], [s_at(i + 1), z_at(j)]);
                        }
                        if j < N {
                            check([s_at(i), z_at(j)], [s_at(i), z_at(j + 1)]);
                        }
                    }
                }
            }
        }
    }
    assert!(found > 500, "only {found} folds were sampled");
    println!("{{\"folds\":{found},\"worst_distance_to_outline_mm\":{worst:e}}}");
}

/// Area of the oracle's material in the plane `z = depth`, by exact chords
/// along scan lines.
fn oracle_area(oracle: &Oracle, [lo, hi]: [P3; 2], depth: f64) -> f64 {
    const LINES: usize = 240;
    const SAMPLES: usize = 600;
    let dy = (hi[1] - lo[1]) / LINES as f64;
    let dx = (hi[0] - lo[0] + 0.2) / SAMPLES as f64;
    let mut area = 0.;
    for line in 0..LINES {
        let y = lo[1] + (line as f64 + 0.5) * dy;
        let at = |x: f64| oracle.inside([x, y, depth]);
        let mut previous = (lo[0] - 0.1, false);
        let mut entered = 0.;
        for i in 1..=SAMPLES {
            let x = lo[0] - 0.1 + i as f64 * dx;
            let state = at(x);
            if state != previous.1 {
                let [mut a, mut b] = [previous.0, x];
                for _ in 0..50 {
                    let mid = 0.5 * (a + b);
                    if at(mid) == previous.1 {
                        a = mid;
                    } else {
                        b = mid;
                    }
                }
                if state {
                    entered = 0.5 * (a + b);
                } else {
                    area += (0.5 * (a + b) - entered) * dy;
                }
            }
            previous = (x, state);
        }
    }
    area
}

fn section_cases() -> Vec<Case> {
    let mut output = Vec::new();
    for (sweep, name, twist, bore, hub, placed, view_name) in [
        (
            GearSweep::Helical,
            "helical",
            0.05,
            3.,
            false,
            false,
            "side",
        ),
        (
            GearSweep::Helical,
            "helical",
            -0.05,
            3.,
            true,
            true,
            "oblique",
        ),
        (
            GearSweep::Chevron,
            "herringbone",
            0.05,
            3.,
            true,
            false,
            "side",
        ),
        (
            GearSweep::Chevron,
            "herringbone",
            -0.05,
            0.,
            false,
            true,
            "oblique",
        ),
        (
            GearSweep::Chevron,
            "herringbone",
            0.05,
            3.,
            false,
            false,
            "oblique",
        ),
        (
            GearSweep::GappedChevron { gap_mm: 1. },
            "double-helical",
            0.05,
            3.,
            true,
            false,
            "side",
        ),
        (
            GearSweep::GappedChevron { gap_mm: 1. },
            "double-helical",
            -0.05,
            0.,
            false,
            true,
            "oblique",
        ),
    ] {
        let spec = specification(sweep, twist, bore, hub);
        let solid = solid(spec, placed);
        let (_, placement) = solid.herringbone_body().unwrap();
        let frame = frames()
            .into_iter()
            .find(|(n, _)| *n == view_name)
            .unwrap()
            .1;
        let view = View::new(frame);
        output.push(Case {
            label: format!(
                "{name} twist {twist} bore {bore} hub {hub} placed {placed} {view_name}"
            ),
            gear: Gear::new(&solid, &view).unwrap().unwrap(),
            oracle: Oracle::new(spec, placement, view),
        });
    }
    output
}

/// Closed, oriented section boundary whose enclosed area matches a
/// scan-line quadrature of the oracle; returns the relative area error and
/// the largest gap between consecutive pieces.
fn check_section(case: &Case, depth: f64, loops: Vec<Loop>) -> (f64, f64) {
    let [lo, hi] = case.gear.bounds();
    let mut worst_gap: f64 = 0.;
    let pieces: Vec<Arc<dyn Path>> = loops.into_iter().flatten().collect();
    let starts: Vec<P3> = pieces.iter().map(|p| p.point(0.)).collect();
    let ends: Vec<P3> = pieces.iter().map(|p| p.point(1.)).collect();
    // Closed loops: every end continues in a start at the same point, and
    // no start continues two ends.
    let mut taken = vec![false; starts.len()];
    for p in &ends {
        assert!((p[2] - depth).abs() <= 1e-9, "{}", case.label);
        let (index, gap) = starts
            .iter()
            .enumerate()
            .filter(|(index, _)| !taken[*index])
            .map(|(index, q)| (index, (p[0] - q[0]).hypot(p[1] - q[1])))
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .unwrap_or_else(|| panic!("{} at {depth}: an end without a start", case.label));
        worst_gap = worst_gap.max(gap);
        assert!(
            gap <= 1e-7,
            "{} at {depth}: boundary end {p:?} continues {gap} mm away",
            case.label
        );
        taken[index] = true;
    }
    let mut area = 0.;
    for piece in &pieces {
        for [a, b] in chords_of(piece.as_ref(), 1e-6) {
            assert!((a[2] - depth).abs() <= 1e-9 && (b[2] - depth).abs() <= 1e-9);
            area += 0.5 * (a[0] + b[0]) * (b[1] - a[1]);
        }
        // Material on the left, void on the right.
        let [p, q] = [piece.point(0.5 - 1e-4), piece.point(0.5 + 1e-4)];
        let m = piece.point(0.5);
        let l = (q[0] - p[0]).hypot(q[1] - p[1]);
        let left = [-(q[1] - p[1]) / l, (q[0] - p[0]) / l];
        let at = |side: f64| {
            case.oracle.inside([
                m[0] + side * 1e-6 * left[0],
                m[1] + side * 1e-6 * left[1],
                depth,
            ])
        };
        assert!(
            at(1.) && !at(-1.),
            "{} at {depth}: orientation at {m:?}",
            case.label
        );
    }
    let expected = oracle_area(&case.oracle, [lo, hi], depth);
    let error = (area - expected).abs() / expected.max(1.);
    println!(
        "{{\"case\":\"{}\",\"depth\":{depth},\"pieces\":{},\"area\":{area},\"quadrature\":{expected}}}",
        case.label,
        pieces.len()
    );
    assert!(
        error <= 2e-3,
        "{} at {depth}: {area} against {expected}",
        case.label
    );
    (error, worst_gap)
}

#[test]
fn sections_are_closed_oriented_and_enclose_the_exact_area() {
    let mut worst_area: f64 = 0.;
    let mut worst_gap: f64 = 0.;
    for case in section_cases() {
        let [lo, hi] = case.gear.bounds();
        for fraction in [0.37, 0.5, 0.81] {
            let depth = lo[2] + (hi[2] - lo[2]) * fraction;
            let pieces = case.gear.section(depth, 1e-3).unwrap();
            assert!(!pieces.is_empty(), "{} at {depth}", case.label);
            let (error, gap) = check_section(&case, depth, pieces);
            worst_area = worst_area.max(error);
            worst_gap = worst_gap.max(gap);
        }
    }
    println!(
        "{{\"worst_relative_area_error\":{worst_area:e},\"worst_loop_gap_mm\":{worst_gap:e}}}"
    );
}

/// Sight lines tangent to the bore, the hub, the root and the tip cylinders
/// and through the axis keep the exact material on both sides of the
/// contact, which contributes no interval end of its own.
#[test]
fn sight_lines_tangent_to_coaxial_cylinders_keep_the_material() {
    let mut random = Random(11);
    let mut checked = 0;
    for case in cases(&["side"]) {
        if case.label.contains("placed true") {
            continue;
        }
        let spec = case.oracle.spec;
        let mut stations = vec![0., spec.root_radius_mm, spec.tip_radius_mm];
        stations.push(spec.bore_radius_mm);
        stations.extend(spec.hub.map(|hub| hub.radius_mm));
        let [lo, hi] = case.gear.bounds();
        for station in stations {
            for sign in [-1., 1.] {
                for _ in 0..6 {
                    let paper = [sign * station, lo[1] + (hi[1] - lo[1]) * random.unit()];
                    let spans = case.gear.sight(paper).unwrap();
                    for _ in 0..40 {
                        let depth = lo[2] + (hi[2] - lo[2]) * random.unit();
                        let p = [paper[0], paper[1], depth];
                        // Skip depths at the contact itself or at an end.
                        let local = case.oracle.local(p);
                        if (local[0].hypot(local[1]) - station).abs() < 1e-6
                            || spans.iter().flatten().any(|end| (end - depth).abs() < 1e-6)
                        {
                            continue;
                        }
                        checked += 1;
                        assert_eq!(
                            spans.iter().any(|[a, b]| *a < depth && depth < *b),
                            case.oracle.inside(p),
                            "{} at {paper:?} depth {depth}: {spans:?}",
                            case.label
                        );
                    }
                }
            }
        }
    }
    assert!(checked > 2000);
}

fn case_in(label: &str, spec: InvoluteGear, frame: PlaneFrame) -> Case {
    let solid = solid(spec, false);
    let (_, placement) = solid.herringbone_body().unwrap();
    let view = View::new(frame);
    Case {
        label: label.to_string(),
        gear: Gear::new(&solid, &view).unwrap().unwrap(),
        oracle: Oracle::new(spec, placement, view),
    }
}

/// A view frame whose direction toward the viewer is `toward`.
fn facing(toward: P3, origin: P3) -> PlaneFrame {
    let length = dot3(toward, toward).sqrt();
    let n = toward.map(|v| v / length);
    let seed = if n[0].abs() < 0.9 {
        [1., 0., 0.]
    } else {
        [0., 1., 0.]
    };
    let along = dot3(seed, n);
    let u: P3 = std::array::from_fn(|i| seed[i] - along * n[i]);
    let ul = dot3(u, u).sqrt();
    let u = u.map(|v| v / ul);
    let v: P3 =
        std::array::from_fn(|i| n[(i + 1) % 3] * u[(i + 2) % 3] - n[(i + 2) % 3] * u[(i + 1) % 3]);
    PlaneFrame { origin, u, v }
}

fn dot3(a: P3, b: P3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// A point inside an involute flank's patch, the ruling of the involute
/// helicoid through it (the helicoid is developable: each ruling is a
/// straight line on the surface) and the surface normal there, all in gear
/// coordinates of an unplaced gear.
fn ruling(gear: &Gear) -> (P3, P3, P3, [f64; 2]) {
    let (index, flank) = gear
        .flanks
        .iter()
        .enumerate()
        .find(|(_, flank)| matches!(flank.shape, super::flank::Shape::Involute { .. }))
        .unwrap();
    let patch = gear.patch(0, index);
    let super::flank::Shape::Involute { start, span, .. } = flank.shape else {
        unreachable!()
    };
    let height = patch.band.z[1] - patch.band.z[0];
    let [s0, zeta0] = [0.5, 0.5];
    let t0 = start + s0 * span;
    let twist = patch.band.twist;
    // Along a ruling the involute parameter changes by -twist per unit
    // height for one hand and +twist for the other; keep the hand on which
    // three points are collinear.
    let at = |hand: f64, lambda: f64| {
        let s = (t0 - hand * twist * lambda - start) / span;
        patch.local(s, zeta0 + lambda / height)
    };
    let hand = [1., -1.]
        .into_iter()
        .min_by(|a, b| {
            let bend = |hand: f64| {
                let [p, q, r] = [-0.2, 0., 0.2].map(|l| at(hand, l));
                (0..3)
                    .map(|i| (p[i] + r[i] - 2. * q[i]).powi(2))
                    .sum::<f64>()
                    .sqrt()
            };
            bend(*a).total_cmp(&bend(*b))
        })
        .unwrap();
    let [p, q, r] = [-0.2, 0., 0.2].map(|l| at(hand, l));
    let bend = (0..3)
        .map(|i| (p[i] + r[i] - 2. * q[i]).powi(2))
        .sum::<f64>()
        .sqrt();
    assert!(
        bend < 1e-12,
        "the involute helicoid is not ruled here: {bend}"
    );
    let direction: P3 = std::array::from_fn(|i| r[i] - q[i]);
    let normal = {
        let tangent = flank.tangent(s0);
        let c = flank.point(s0);
        let [x, y] = super::flank::rotate([tangent[1], -tangent[0]], patch.band.phase_at(q[2]));
        [x, y, twist * (c[0] * tangent[0] + c[1] * tangent[1])]
    };
    // The range of lambda (per unit of `direction`) on the patch.
    let unit = 0.2;
    let mut range = [f64::INFINITY, f64::NEG_INFINITY];
    for k in -10_000..=10_000 {
        let lambda = f64::from(k) * 0.01;
        let s = (t0 - hand * twist * lambda * unit - start) / span;
        let zeta = zeta0 + lambda * unit / height;
        if (0. ..=1.).contains(&s) && (0. ..=1.).contains(&zeta) {
            range = [range[0].min(lambda), range[1].max(lambda)];
        }
    }
    (q, direction, normal, range)
}

/// A sight line along a ruling of an involute flank runs in the surface: it
/// only touches the boundary there, and the material elsewhere on the line
/// is exact. The drawing in that view draws.
#[test]
fn a_sight_line_along_an_involute_ruling_touches_without_material() {
    for sweep in [
        GearSweep::Helical,
        GearSweep::Chevron,
        GearSweep::GappedChevron { gap_mm: 1. },
    ] {
        let spec = specification(sweep, 0.05, 3., false);
        let probe = case_in("probe", spec, PlaneFrame::default());
        let (point, direction, _, [lo, hi]) = ruling(&probe.gear);
        let frame = facing(direction, [0.; 3]);
        let case = case_in(&format!("{sweep:?} along a ruling"), spec, frame);
        let view = View::new(frame);
        let p = view.point(point);
        let spans = case.gear.sight([p[0], p[1]]).unwrap();
        // The stretch on the flank lies in no interval.
        let step = dot3(direction, direction).sqrt();
        for k in 1..20 {
            let lambda = lo + (hi - lo) * f64::from(k) / 20.;
            let depth = p[2] + lambda * step;
            assert!(
                !spans.iter().any(|[a, b]| *a < depth && depth < *b),
                "{}: the ruling at depth {depth} is inside {spans:?}",
                case.label
            );
        }
        // Elsewhere the line agrees with the oracle.
        let [blo, bhi] = case.gear.bounds();
        let mut checked = 0;
        for k in 0..400 {
            let depth = blo[2] - 1. + (bhi[2] - blo[2] + 2.) * f64::from(k) / 400.;
            let lambda = (depth - p[2]) / step;
            if (lo - 0.02..=hi + 0.02).contains(&lambda)
                || spans.iter().flatten().any(|end| (end - depth).abs() < 1e-6)
            {
                continue;
            }
            checked += 1;
            assert_eq!(
                spans.iter().any(|[a, b]| *a < depth && depth < *b),
                case.oracle.inside([p[0], p[1], depth]),
                "{} at depth {depth}: {spans:?}",
                case.label
            );
        }
        assert!(checked > 100);
        let solid = solid(spec, false);
        for options in [
            DrawingProjection {
                frame: facing(direction, point),
                hidden_lines: true,
                ..Default::default()
            },
            DrawingProjection {
                frame: facing(direction, point),
                section_offset_mm: Some(0.),
                ..Default::default()
            },
        ] {
            let vectors = crate::project_drawing(
                &[DrawingSolid {
                    solid: &solid,
                    hatch: None,
                }],
                options,
            )
            .unwrap();
            assert!(!vectors.is_empty(), "{}", case.label);
        }
    }
}

/// Planes tangent to a flank: along the axis at the tip and root radii
/// (touching the tip and root cylinders along generators) and an oblique
/// plane containing a ruling of an involute flank (the tangent plane of a
/// developable surface is constant along its ruling). The sections are
/// closed, oriented and enclose the exact area.
#[test]
fn sections_by_planes_tangent_to_a_flank_are_exact() {
    for sweep in [
        GearSweep::Helical,
        GearSweep::Chevron,
        GearSweep::GappedChevron { gap_mm: 1. },
    ] {
        let spec = specification(sweep, 0.05, 3., false);
        let side = frames().into_iter().find(|(n, _)| *n == "side").unwrap().1;
        let case = case_in(&format!("{sweep:?} side"), spec, side);
        // The side view looks along +y: depth is -y.
        for radius in [spec.tip_radius_mm, spec.root_radius_mm] {
            for depth in [radius, -radius] {
                let loops = case
                    .gear
                    .section(depth, 1e-3)
                    .unwrap_or_else(|error| panic!("{} at {depth}: {error:?}", case.label));
                check_section(&case, depth, loops);
            }
        }
        let probe = case_in("probe", spec, PlaneFrame::default());
        let (point, _, normal, _) = ruling(&probe.gear);
        for sign in [1., -1.] {
            let frame = facing(normal.map(|v| sign * v), [0.; 3]);
            let case = case_in(&format!("{sweep:?} tangent to an involute"), spec, frame);
            let depth = View::new(frame).point(point)[2];
            let loops = case.gear.section(depth, 1e-3).unwrap();
            check_section(&case, depth, loops);
        }
    }
}

/// `line` gives the material along arbitrary lines of view space, not only
/// sight lines: random lines through the solid against the oracle, with
/// interval ends on the actual surface, for every family and placement.
#[test]
fn material_along_arbitrary_lines_matches_the_oracle() {
    let mut random = Random(23);
    let mut ends = 0;
    let mut worst: f64 = 0.;
    for case in cases(&["oblique", "side"]) {
        let [lo, hi] = case.gear.bounds();
        let scale = lo.iter().chain(&hi).fold(1_f64, |s, v| s.max(v.abs()));
        for _ in 0..40 {
            let start: P3 = std::array::from_fn(|k| lo[k] + (hi[k] - lo[k]) * random.unit());
            let direction: P3 = std::array::from_fn(|_| 2. * random.unit() - 1.);
            let spans = case.gear.line(start, direction).unwrap().unwrap();
            assert!(spans.windows(2).all(|pair| pair[0][1] < pair[1][0]));
            let at = |s: f64| std::array::from_fn(|k| start[k] + s * direction[k]);
            for k in 0..60 {
                let s = -3. + 6. * f64::from(k) / 60. + 0.0123;
                if spans.iter().flatten().any(|end| (end - s).abs() < 1e-6) {
                    continue;
                }
                assert_eq!(
                    spans.iter().any(|[a, b]| *a < s && s < *b),
                    case.oracle.inside(at(s)),
                    "{} along {start:?} + s {direction:?} at {s}: {spans:?}",
                    case.label
                );
            }
            for (index, &end) in spans.iter().flatten().enumerate() {
                let h = 1e-5;
                if spans
                    .iter()
                    .flatten()
                    .enumerate()
                    .any(|(other, stop)| other != index && (stop - end).abs() < 4. * h)
                {
                    continue;
                }
                let inside = |s: f64| case.oracle.inside(at(s));
                let [mut a, mut b] = [end - h, end + h];
                let first = inside(a);
                if first == inside(b) {
                    continue;
                }
                for _ in 0..60 {
                    let mid = 0.5 * (a + b);
                    if inside(mid) == first {
                        a = mid;
                    } else {
                        b = mid;
                    }
                }
                let length = dot3(direction, direction).sqrt();
                let residual = (0.5 * (a + b) - end).abs() * length / scale;
                worst = worst.max(residual);
                ends += 1;
                assert!(
                    residual <= 1e-9,
                    "{}: end {end} off by {residual}",
                    case.label
                );
            }
        }
    }
    assert!(ends > 500, "only {ends} ends");
    println!("{{\"line_ends\":{ends},\"worst_relative_residual\":{worst:e}}}");
}
