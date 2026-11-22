use super::*;
use std::cell::Cell;
fn source(scale: f64, left: bool, posed: bool) -> Rc<OwnedDerivedMaterial> {
    super::super::tests::source(scale, left, posed)
}
fn feature() -> EntityId {
    EntityId::from_uuid_bits(0x00ca_d352_0000_0000_0000_0000_0000_0002)
}
fn solid(source: &Rc<OwnedDerivedMaterial>, frame: PlaneFrame) -> Solid {
    Solid::rust_repaired_phillips_material(
        Rc::clone(source),
        feature(),
        frame,
        Default::default(),
        None,
    )
    .unwrap()
}
// Independent source-datum map. It uses the original four chart columns and
// the actual occurrence normal, never Placement::point or its inverse.
fn world(source: &OwnedDerivedMaterial, frame: PlaneFrame, p: Point) -> [[f64; 2]; 3] {
    let chart = source.chart();
    let z = I::point(p[2]) - I::point(chart.outer_top_z);
    let constructed: [I; 3] = std::array::from_fn(|i| {
        I::point(chart.columns[0][i])
            + I::point(chart.columns[1][i]) * I::point(p[0])
            + I::point(chart.columns[2][i]) * I::point(p[1])
            + I::point(chart.columns[3][i]) * z
    });
    let n = frame.normal();
    std::array::from_fn(|i| {
        let b = I::point(frame.origin[i])
            + I::point(frame.u[i]) * constructed[0]
            + I::point(frame.v[i]) * constructed[1]
            + I::point(n[i]) * constructed[2];
        [b.lo, b.hi]
    })
}
fn point(bounds: [[f64; 2]; 3]) -> Point {
    bounds.map(|b| b[0].midpoint(b[1]))
}
#[test]
fn actual_owned_world_queries_keep_station_pose_and_real_material() {
    for (scale, left, posed, frame) in [
        (1., false, false, PlaneFrame::default()),
        (
            1.137,
            true,
            true,
            PlaneFrame {
                origin: [10., 20., 30.],
                u: [0., 1., 0.],
                v: [-1., 0., 0.],
            },
        ),
    ] {
        let source = source(scale, left, posed);
        let part = solid(&source, frame);
        let original = source.original_input().to_owned();
        let top = source.head().view().station(0).unwrap();
        let cap = source.thread().cap_axes();
        let z = (cap[0].origin[2] + cap[1].origin[2]) * 0.5;
        let r = source.local_enclosure_bounds()[0][1];
        for (local, expected) in [
            ([0., 0., top[1] + 0.1 * scale], Classification::Outside),
            (
                [0.99 * top[0], 0., top[1] + 0.001 * scale],
                Classification::Inside,
            ),
            ([0., 0., z], Classification::Inside),
            ([r + scale, 0., z], Classification::Outside),
        ] {
            let bounds = world(&source, frame, local);
            let receipt = part
                .rust_repaired_phillips_classify_world_box(bounds, Default::default(), None)
                .unwrap();
            assert_eq!(receipt.location(), expected, "local {local:?}");
            assert!(std::ptr::eq(receipt.source(), source.as_ref()));
            assert_eq!(receipt.world_bounds_mm(), bounds);
            assert_eq!(receipt.occurrence(), frame);
            for (b, x) in receipt.local_bounds_mm().into_iter().zip(local) {
                assert!(b[0] <= x && x <= b[1]);
            }
            assert_eq!(receipt.work_used() + receipt.remaining_work(), MAXIMUM_WORK);
            assert_eq!(
                part.rust_thread_contains(point(bounds)).unwrap(),
                expected == Classification::Inside
            );
            let exact = RepairedPhillipsQueryOptions {
                maximum_work: receipt.work_used(),
                ..Default::default()
            };
            assert_eq!(
                part.rust_repaired_phillips_classify_world_box(bounds, exact, None)
                    .unwrap()
                    .remaining_work(),
                0
            );
            assert!(
                part.rust_repaired_phillips_classify_world_box(
                    bounds,
                    RepairedPhillipsQueryOptions {
                        maximum_work: exact.maximum_work - 1,
                        ..exact
                    },
                    None
                )
                .is_err()
            );
        }
        // Raw mouth origin is exactly the construction datum. Its material
        // station is topZ. The actual finite thread-tip axis is independent.
        for local in [[0., 0., top[1]], cap[1].origin] {
            let receipt = part
                .rust_repaired_phillips_classify_world_box(
                    world(&source, frame, local),
                    Default::default(),
                    None,
                )
                .unwrap();
            for (b, x) in receipt.local_bounds_mm().into_iter().zip(local) {
                assert!(b[0] <= x && x <= b[1]);
            }
        }
        // This shared closure disk is not an exposed-boundary oracle. The
        // whole tolerance box crosses BOTH source component caps.
        let join = world(&source, frame, cap[0].origin);
        assert_eq!(
            part.rust_thread_contains(point(join)).unwrap_err().kind,
            GeometryErrorKind::Unsupported
        );
        assert_eq!(source.original_input(), original);
    }
}
#[test]
fn actual_query_cancellation_and_provider_limits_keep_source() {
    let source = source(1., false, false);
    let part = solid(&source, Default::default());
    let caps = source.thread().cap_axes();
    let bounds = world(
        &source,
        Default::default(),
        [0., 0., (caps[0].origin[2] + caps[1].origin[2]) * 0.5],
    );
    let polls = Cell::new(0usize);
    let count = || {
        polls.set(polls.get() + 1);
        false
    };
    let baseline = part
        .rust_repaired_phillips_classify_world_box(bounds, Default::default(), Some(&count))
        .unwrap();
    let original = source.original_input().to_owned();
    for stop in 0..polls.get() {
        let n = Cell::new(0usize);
        let cancel = || {
            let k = n.get();
            n.set(k + 1);
            k == stop
        };
        assert_eq!(
            part.rust_repaired_phillips_classify_world_box(
                bounds,
                Default::default(),
                Some(&cancel)
            )
            .err()
            .unwrap()
            .kind,
            GeometryErrorKind::Cancelled
        );
    }
    for options in [
        RepairedPhillipsQueryOptions {
            maximum_work: MAXIMUM_WORK + 1,
            ..Default::default()
        },
        RepairedPhillipsQueryOptions {
            maximum_head_cells: 0,
            ..Default::default()
        },
        RepairedPhillipsQueryOptions {
            maximum_head_cells: 1_000_001,
            ..Default::default()
        },
    ] {
        assert_eq!(
            part.rust_repaired_phillips_classify_world_box(bounds, options, None)
                .err()
                .unwrap()
                .kind,
            GeometryErrorKind::Input
        );
    }
    for (maximum_work, origin) in [
        (0, "work exhausted"),
        (INVERSE_WORK + 4160, "Head(Limit)"),
        (INVERSE_WORK + 4162, "Thread(Limit)"),
    ] {
        let error = part
            .rust_repaired_phillips_classify_world_box(
                bounds,
                RepairedPhillipsQueryOptions {
                    maximum_work,
                    ..Default::default()
                },
                None,
            )
            .err()
            .unwrap();
        assert_eq!(error.kind, GeometryErrorKind::Limit);
        assert!(error.message.contains(origin));
    }
    for bad in [[1., -1.], [f64::NAN, 0.], [0., f64::INFINITY]] {
        let mut invalid = bounds;
        invalid[0] = bad;
        assert_eq!(
            part.rust_repaired_phillips_classify_world_box(invalid, Default::default(), None)
                .err()
                .unwrap()
                .kind,
            GeometryErrorKind::Input
        );
    }
    assert_eq!(
        part.rust_repaired_phillips_classify_world_box(
            [[f64::MAX; 2]; 3],
            RepairedPhillipsQueryOptions {
                maximum_work: INVERSE_WORK,
                ..Default::default()
            },
            None
        )
        .err()
        .unwrap()
        .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        part.rust_repaired_phillips_classify_world_box(bounds, Default::default(), None)
            .unwrap()
            .location(),
        baseline.location()
    );
    assert_eq!(source.original_input(), original);
    let cube = Solid::box_solid([1.; 3], Default::default()).unwrap();
    assert_eq!(
        cube.rust_thread_contains([0.5; 3]).unwrap_err().kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        cube.rust_repaired_phillips_classify_world_box(bounds, Default::default(), None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
}
#[test]
fn actual_fourth_column_inverse_and_source_error_origins_are_preserved() {
    let source = source(1., false, false);
    let mut chart = source.chart();
    chart.columns[3][0] = 2_f64.powi(-30);
    let frame = PlaneFrame {
        origin: [10., 20., 30.],
        u: [0., 1., 0.],
        v: [-1., 0., 0.],
    };
    let placement = Placement::new(chart, frame).unwrap();
    let p = [0.25, 0.5, 0.4];
    let dz = I::point(p[2]) - I::point(chart.outer_top_z);
    let constructed = [
        I::point(p[0]) + I::point(chart.columns[3][0]) * dz,
        I::point(p[1]),
        I::point(p[2]),
    ];
    let image = [
        I::point(10.) - constructed[1],
        I::point(20.) + constructed[0],
        I::point(30.) + constructed[2],
    ];
    let inverse = placement
        .inverse_bounds(image.map(|b| [b.lo, b.hi]))
        .unwrap();
    for (b, x) in inverse.into_iter().zip(p) {
        assert!(b[0] <= x && x <= b[1]);
    }
    // Assuming N=cross(U,V) would miss a real nonzero horizontal shear.
    assert!((constructed[0] - I::point(p[0])).lo > 0.);
    chart.columns[3] = [0.; 3];
    assert!(Placement::new(chart, frame).is_err());
    for (cause, expected) in [
        (E::Invalid, GeometryErrorKind::Input),
        (E::Cancelled, GeometryErrorKind::Cancelled),
        (E::Limit, GeometryErrorKind::Limit),
        (E::Head(H::Unresolved), GeometryErrorKind::Unsupported),
        (E::Thread(T::Unresolved), GeometryErrorKind::Unsupported),
        (E::Cavity(C::Unresolved), GeometryErrorKind::Unsupported),
    ] {
        let error = query_failure(cause);
        assert_eq!(error.kind, expected);
        assert!(error.message.contains(&format!("{cause:?}")));
    }
}
