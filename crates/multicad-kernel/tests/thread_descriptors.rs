//! Read-only thread descriptors: reading one preserves the solid and the
//! descriptor alone reproduces the exact modeled material.
use multicad_core::document::EntityId;
use multicad_kernel::{
    Bounds, GeometryPolicy, HelicalThread, PlaneFrame, RigidTransform, Solid, rust_kernel::Rational,
};
use std::f64::consts::PI;

fn iso(length_mm: f64, left_hand: bool, starts: u32) -> HelicalThread {
    HelicalThread {
        major_radius_mm: 3.,
        pitch_mm: 1.,
        length_mm,
        depth_mm: 5. * 3_f64.sqrt() / 16.,
        crest_width_mm: 0.125,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        radial_taper: 0.,
        starts,
        left_hand,
    }
}

fn tilted() -> PlaneFrame {
    PlaneFrame {
        origin: [5., -3., 2.],
        u: [0.8, 0.6, 0.],
        v: [-0.36, 0.48, 0.8],
    }
}

fn normal(frame: PlaneFrame) -> [f64; 3] {
    let [u, v] = [frame.u, frame.v];
    [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ]
}

fn local(frame: PlaneFrame, p: [f64; 3]) -> [f64; 3] {
    let d: [f64; 3] = std::array::from_fn(|i| p[i] - frame.origin[i]);
    [frame.u, frame.v, normal(frame)].map(|a| a[0] * d[0] + a[1] * d[1] + a[2] * d[2])
}

/// Everything observable about a solid that a read must leave unchanged.
#[derive(Debug, PartialEq)]
struct Fingerprint {
    volume: f64,
    bounds: Bounds,
    topology: multicad_core::references::TopologySnapshot,
    step: Vec<u8>,
}

fn fingerprint(solid: &Solid, owner: EntityId, name: &str) -> Fingerprint {
    let path = std::env::temp_dir().join(format!(
        "multicad-thread-descriptor-{}-{name}.step",
        std::process::id()
    ));
    solid.write_step(&path).unwrap();
    let step = std::fs::read(&path).unwrap();
    std::fs::remove_file(&path).unwrap();
    Fingerprint {
        volume: solid.volume_mm3().unwrap(),
        bounds: solid.bounds().unwrap(),
        topology: solid.topology_snapshot(owner, 3).unwrap(),
        step,
    }
}

/// World sample points around a frame's axis, inside and just outside a
/// cylinder of `radius` over `[-0.5, height + 0.5]`.
fn samples(frame: PlaneFrame, radius: f64, height: f64) -> Vec<[f64; 3]> {
    let n = normal(frame);
    let mut out = Vec::new();
    for i in 0..23 {
        let z = -0.5 + (height + 1.) * f64::from(i) / 22.;
        for j in 0..29 {
            let r = 1.2 * radius * f64::from(j) / 28.;
            for k in 0..11 {
                let a = 2. * PI * f64::from(k) / 11. + 0.013 * f64::from(i);
                let [x, y] = [r * a.cos(), r * a.sin()];
                out.push(std::array::from_fn(|c| {
                    frame.origin[c] + x * frame.u[c] + y * frame.v[c] + z * n[c]
                }));
            }
        }
    }
    out
}

fn pose() -> RigidTransform {
    RigidTransform {
        rotation: [[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]],
        translation_mm: [100., -30., 50.],
    }
}

#[test]
fn helical_thread_descriptor_preserves_and_reproduces_the_plain_thread() {
    let policy = GeometryPolicy::default();
    let tolerance = policy.linear_tolerance_mm;
    for (starts, left, root) in [(1, false, Some(1. / (4. * 3_f64.sqrt()))), (2, true, None)] {
        let spec = match root {
            Some(root) => HelicalThread {
                depth_mm: (1. - 0.125) / (2. * 30_f64.to_radians().tan()) - root,
                ..iso(6.6, left, starts)
            },
            None => iso(6.6, left, starts),
        };
        let owner = EntityId::new();
        let solid = Solid::rust_helical_thread(spec, root, owner, tilted(), policy).unwrap();
        let before = fingerprint(&solid, owner, "helical-before");
        let (thread, frame) = solid.rust_helical_thread_data().unwrap();
        assert_eq!(fingerprint(&solid, owner, "helical-after"), before);
        assert_eq!(solid.rust_helical_thread_data(), Some((thread, frame)));
        assert_eq!(frame, tilted());
        let s = thread.spec();
        assert_eq!(
            [s.major_radius_mm, s.pitch_mm, s.length_mm, s.crest_width_mm],
            [
                spec.major_radius_mm,
                spec.pitch_mm,
                spec.length_mm,
                spec.crest_width_mm
            ]
        );
        assert_eq!(
            (s.starts, s.left_hand, s.root_radius_mm),
            (starts, left, root)
        );
        assert_eq!(before.volume, thread.volume_mm3());
        for p in samples(frame, spec.major_radius_mm, spec.length_mm) {
            assert_eq!(
                solid.rust_thread_contains(p).unwrap(),
                thread.contains(local(frame, p), tolerance),
                "{p:?}"
            );
        }
        let moved = solid.transformed(pose()).unwrap();
        let (moved_thread, placed) = moved.rust_helical_thread_data().unwrap();
        assert_eq!(moved_thread, thread);
        for p in samples(placed, spec.major_radius_mm, spec.length_mm)
            .into_iter()
            .step_by(7)
        {
            assert_eq!(
                moved.rust_thread_contains(p).unwrap(),
                thread.contains(local(placed, p), tolerance)
            );
        }
        assert!(solid.rust_internal_thread_data().is_none());
        assert!(solid.rust_hex_nut_geometry().is_none());
    }
}

/// The stock less the void, with the full stock disc beyond the void as
/// a partial-depth cut's floor, written from the descriptor alone.
fn internal_contains(data: &multicad_kernel::InternalThreadData, p: [f64; 3], t: f64) -> bool {
    let [x, y, z] = local(data.frame, p);
    let length = data.thread.spec().length_mm;
    let r = x.hypot(y);
    z >= -t
        && z <= data.stock_length_mm + t
        && r <= data.stock_radius_mm + t
        && ((data.stock_length_mm > length && z >= length - t)
            || r >= data.thread.radius(y.atan2(x), z.clamp(0., length)) - t)
}

#[test]
fn internal_thread_descriptor_preserves_and_reproduces_stock_floor_and_overrun() {
    let policy = GeometryPolicy::default();
    let tolerance = policy.linear_tolerance_mm;
    let stock_frame = tilted();
    let n = normal(stock_frame);
    // The cap at the far end, entered against the stock normal.
    let top = PlaneFrame {
        origin: std::array::from_fn(|i| stock_frame.origin[i] + 4. * n[i]),
        u: stock_frame.u,
        v: stock_frame.v.map(|x| -x),
    };
    for (length, entry, left) in [
        (5.2, stock_frame, false),
        (4., stock_frame, true),
        (2.5, stock_frame, false),
        (2.5, top, true),
        (5.2, top, false),
    ] {
        let stock = Solid::cylinder(stock_frame, 5., 4., policy).unwrap();
        let owner = EntityId::new();
        let spec = iso(length, left, 1);
        let solid = stock
            .rust_cut_internal_thread(spec, 0.2, owner, entry, policy)
            .unwrap();
        let before = fingerprint(&solid, owner, "internal-before");
        let data = solid.rust_internal_thread_data().unwrap();
        assert_eq!(fingerprint(&solid, owner, "internal-after"), before);
        assert_eq!(solid.rust_internal_thread_data(), Some(data));
        assert_eq!(data.stock_radius_mm, 5.);
        assert_eq!(data.stock_length_mm, 4.);
        assert_eq!(data.tool.spec().length_mm, length);
        assert_eq!(data.thread.spec().length_mm, length.min(4.));
        if length <= 4. {
            assert_eq!(data.thread, data.tool);
        } else {
            assert_ne!(data.thread, data.tool);
            assert_eq!(data.thread.helical_rate(), data.tool.helical_rate());
        }
        let intervals = solid.rust_receiving_thread_intervals().unwrap();
        assert_eq!(intervals.frame, data.frame);
        assert_eq!(intervals.requested_tool_length_mm, length);
        assert_eq!(intervals.physical_thread_length_mm, length.min(4.));
        let close = |a: [f64; 3], b: [f64; 3]| (0..3).all(|i| (a[i] - b[i]).abs() < 1e-12);
        assert!(close(data.frame.origin, entry.origin));
        assert!(close(normal(data.frame), normal(entry)));
        assert!(close(data.stock_frame.origin, stock_frame.origin));
        assert!(close(normal(data.stock_frame), n));
        let expected = PI * 25. * 4. - data.thread.volume_mm3();
        assert!((before.volume - expected).abs() <= 1e-9 * expected);
        for p in samples(data.frame, 5., 4.) {
            assert_eq!(
                solid.rust_thread_contains(p).unwrap(),
                internal_contains(&data, p, tolerance),
                "{length} {p:?}"
            );
        }
        let moved = solid.transformed(pose()).unwrap();
        let placed = moved.rust_internal_thread_data().unwrap();
        assert_eq!(
            (
                placed.thread,
                placed.tool,
                placed.stock_radius_mm,
                placed.stock_length_mm
            ),
            (data.thread, data.tool, 5., 4.)
        );
        for p in samples(placed.frame, 5., 4.).into_iter().step_by(7) {
            assert_eq!(
                moved.rust_thread_contains(p).unwrap(),
                internal_contains(&placed, p, tolerance)
            );
        }
        assert!(solid.rust_helical_thread_data().is_none());
        assert!(solid.rust_hex_nut_geometry().is_none());
    }
}

#[test]
fn hex_nut_descriptor_preserves_and_reproduces_the_plain_hex_nut() {
    let policy = GeometryPolicy::default();
    let tolerance = policy.linear_tolerance_mm;
    let spec = HelicalThread {
        major_radius_mm: 1.55,
        pitch_mm: 0.5,
        length_mm: 2.4,
        depth_mm: 0.28,
        crest_width_mm: 0.08,
        ..iso(2.4, false, 1)
    };
    let owner = EntityId::new();
    let solid = Solid::rust_hex_nut(spec, 5.5, 0.2, owner, tilted(), policy).unwrap();
    let before = fingerprint(&solid, owner, "nut-before");
    let (nut, frame) = solid.rust_hex_nut_geometry().map(|(n, f)| (*n, f)).unwrap();
    assert_eq!(fingerprint(&solid, owner, "nut-after"), before);
    assert_eq!(frame, tilted());
    assert_eq!(nut.across_flats_mm(), 5.5);
    assert_eq!(before.volume, nut.volume_mm3());
    for p in samples(frame, 5.5 / 3_f64.sqrt(), 2.4) {
        assert_eq!(
            solid.rust_thread_contains(p).unwrap(),
            nut.contains(local(frame, p), tolerance)
        );
    }
    assert!(solid.rust_helical_thread_data().is_none());
    assert!(solid.rust_internal_thread_data().is_none());
}

#[test]
fn descriptors_refuse_other_bodies() {
    let policy = GeometryPolicy::default();
    let boxed = Solid::rust_box(
        EntityId::new().uuid_bits(),
        [Rational::new(20, 1).unwrap(); 3],
        PlaneFrame::default(),
        policy,
    )
    .unwrap();
    let relief = Solid::rust_helical_thread_end_relief(
        iso(6.6, false, 1),
        None,
        0.6,
        0.1,
        EntityId::new(),
        PlaneFrame::default(),
        policy,
    )
    .unwrap();
    let cylinder = Solid::cylinder(PlaneFrame::default(), 5., 4., policy).unwrap();
    for solid in [&boxed, &relief, &cylinder] {
        assert!(solid.rust_helical_thread_data().is_none());
        assert!(solid.rust_internal_thread_data().is_none());
    }
}
