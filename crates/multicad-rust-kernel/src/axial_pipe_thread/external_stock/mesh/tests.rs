use super::*;
use crate::axial_pipe_thread::authored::{AuthoredSection, AuthoredSpec, AuthoredThread};
use multicad_precision::PortableMath as _;

fn engineering_cut(flat: f64) -> ExternalStockCut {
    let taper = 1. / 32.;
    let slope = 3_f64.sqrt();
    let slopes = 1. / (slope - taper) + 1. / (slope + taper);
    let source = AuthoredThread::new(
        AuthoredSpec {
            major_radius_mm: 5.,
            pitch_mm: 1.,
            length_mm: 6.,
            radial_taper: taper,
            section: AuthoredSection::NptFlat {
                depth_mm: (1. - 2. * flat) / slopes,
                flat_width_mm: flat,
            },
            crest_radius_limit_mm: Some(5.1),
        },
        1e-7,
        1e-5,
    )
    .unwrap();
    let thread = source.stock_clip_plan(Default::default(), None).unwrap();
    let zero = Rational::new(0, 1).unwrap();
    let stock = crate::cylinder::AxialCylinder::new(
        73,
        crate::Circle::new([zero; 2], Rational::new(51, 10).unwrap()).unwrap(),
        zero,
        Rational::new(7, 1).unwrap(),
    )
    .unwrap();
    ExternalStockCut::new(thread, stock, Default::default(), None).unwrap()
}
fn quality() -> Tessellation {
    Tessellation {
        deflection: 0.15,
        angle: 0.5,
        triangle_limit: 500_000,
    }
}
fn same_physical_point(a: [[f64; 2]; 3], b: [[f64; 2]; 3]) {
    for (a, b) in a.into_iter().zip(b) {
        assert!(
            a[0] <= b[1] + 1e-12 && b[0] <= a[1] + 1e-12,
            "{a:?} and {b:?}"
        );
        assert!((a[0].midpoint(a[1]) - b[0].midpoint(b[1])).abs() < 1e-10);
    }
}

#[test]
fn actual_end_rail_pairs_select_their_physical_ring_curve_in_both_directions() {
    let cut = engineering_cut(0.001);
    let plan = cut.thread();
    let spec = plan.source().spec();
    let base = plan
        .vertex_bounds(
            StockClipVertex {
                station: StockClipStation::Start,
                rail: StockClipRail::RootLower,
            },
            Default::default(),
            None,
        )
        .unwrap();
    for band in BANDS {
        let curve = StockClipCurve::Ring {
            station: StockClipStation::End,
            band,
        };
        let [a, b] = curve.endpoints();
        let fields = plan
            .curve_law(curve, Default::default(), None)
            .unwrap()
            .affine_fields_mm()
            .map(|p| p.map(interval));
        for (vertex, t) in [(a, 0.), (b, 1.)] {
            let v = plan
                .vertex_bounds(vertex, Default::default(), None)
                .unwrap();
            let c = plan
                .curve_bounds(curve, [t, t], Default::default(), None)
                .unwrap();
            same_physical_point(v.position_mm(), c.position_mm());
            let value = parameter(
                StockClipPointOwner::Vertex(vertex),
                band,
                fields,
                interval(base.phase_mm()),
                p(spec.pitch_mm),
                p(spec.length_mm),
                128,
            )
            .unwrap();
            assert_eq!([value.lo, value.hi], [t, t]);
        }
        assert_eq!(
            owner_band(
                StockClipPointOwner::Vertex(a),
                StockClipPointOwner::Vertex(b)
            )
            .unwrap(),
            band
        );
        assert_eq!(
            owner_band(
                StockClipPointOwner::Vertex(b),
                StockClipPointOwner::Vertex(a)
            )
            .unwrap(),
            band
        );
        // A third actual end vertex is not an endpoint of this finite curve.
        let other = RAILS
            .into_iter()
            .map(|rail| StockClipVertex {
                station: StockClipStation::End,
                rail,
            })
            .find(|v| *v != a && *v != b)
            .unwrap();
        let other_point = plan
            .vertex_bounds(other, Default::default(), None)
            .unwrap()
            .position_mm();
        for endpoint in [a, b] {
            let endpoint = plan
                .vertex_bounds(endpoint, Default::default(), None)
                .unwrap()
                .position_mm();
            assert!(
                other_point
                    .into_iter()
                    .zip(endpoint)
                    .any(|(a, b)| (a[0].midpoint(a[1]) - b[0].midpoint(b[1])).abs() > 1e-4)
            );
        }
        assert!(matches!(
            parameter(
                StockClipPointOwner::Vertex(other),
                band,
                fields,
                interval(base.phase_mm()),
                p(spec.pitch_mm),
                p(spec.length_mm),
                128
            ),
            Err(Error::Unresolved)
        ));
        let start = StockClipVertex {
            station: StockClipStation::Start,
            rail: a.rail,
        };
        assert!(
            plan.vertex_bounds(start, Default::default(), None)
                .unwrap()
                .axial_mm()[1]
                < plan
                    .vertex_bounds(a, Default::default(), None)
                    .unwrap()
                    .axial_mm()[0]
        );
        assert!(matches!(
            owner_band(
                StockClipPointOwner::Vertex(start),
                StockClipPointOwner::Vertex(b)
            ),
            Err(Error::Unresolved)
        ));
        assert!(matches!(
            owner_band(
                StockClipPointOwner::Vertex(a),
                StockClipPointOwner::Vertex(start)
            ),
            Err(Error::Unresolved)
        ));
        assert!(matches!(
            parameter(
                StockClipPointOwner::Vertex(start),
                band,
                fields,
                interval(base.phase_mm()),
                p(spec.pitch_mm),
                p(spec.length_mm),
                128
            ),
            Err(Error::Unresolved)
        ));
    }
    let a = StockClipVertex {
        station: StockClipStation::End,
        rail: StockClipRail::RootLower,
    };
    let opposite = StockClipVertex {
        station: StockClipStation::End,
        rail: StockClipRail::CrestUpper,
    };
    assert!(matches!(
        owner_band(
            StockClipPointOwner::Vertex(a),
            StockClipPointOwner::Vertex(opposite)
        ),
        Err(Error::Unresolved)
    ));
}

#[test]
fn genuine_thin_stock_mesh_keeps_source_and_reports_actual_cap_owner_routes() {
    let cut = engineering_cut(0.01);
    let before = cut.thread().source().spec();
    let q = quality();
    let cells = cut
        .thread()
        .cell_plan(128, 32, 512_000_000, 64 * 1024 * 1024, None)
        .unwrap();
    let old = cells
        .mesh_candidate(q, 512_000_000, 64 * 1024 * 1024, None)
        .unwrap();
    let mut vertex_pairs = 0;
    let mut curve_pairs = 0;
    for (ids, face) in old
        .mesh()
        .indices
        .as_chunks::<3>()
        .0
        .iter()
        .zip(&old.mesh().face_indices)
    {
        if *face != 6 {
            continue;
        }
        let a = cells.vertex_owner(ids[1] as usize).unwrap();
        let b = cells.vertex_owner(ids[2] as usize).unwrap();
        let actual = owner_band(a, b).unwrap();
        if matches!(a, StockClipPointOwner::Vertex(_))
            && matches!(b, StockClipPointOwner::Vertex(_))
        {
            vertex_pairs += 1;
        } else {
            curve_pairs += 1;
        }
        let curve = StockClipCurve::Ring {
            station: StockClipStation::End,
            band: actual,
        };
        let whole = cut
            .thread()
            .curve_bounds(curve, [0., 1.], Default::default(), None)
            .unwrap()
            .position_mm();
        for id in [ids[1], ids[2]] {
            let literal = &old.mesh().vertices[id as usize * 3..id as usize * 3 + 3];
            for k in 0..3 {
                assert!(literal[k] >= whole[k][0] - 1e-9 && literal[k] <= whole[k][1] + 1e-9);
            }
        }
    }
    eprintln!(
        "THIN_ACTUAL_CAP_ROUTES vertex_pairs={vertex_pairs} curve_pairs={curve_pairs} grid={:?}",
        cells.subdivisions()
    );
    assert!(curve_pairs > 0);
    let candidate = cut
        .mesh_candidate(q, 512_000_000, 64 * 1024 * 1024, None)
        .unwrap();
    assert!(candidate.connectivity().signed_volume_mm3 > 0.);
    assert_eq!(
        manufacturing::certify(candidate.mesh()).unwrap(),
        candidate.connectivity()
    );
    assert!(candidate.source_distance_mm() <= q.deflection);
    assert!(candidate.normal_chord() <= 2. * (q.angle / 2.).portable_sin());
    assert_eq!(candidate.source().thread().source().spec(), before);
    eprintln!(
        "THIN_ACTUAL_MESH triangles={} work={} bytes={}",
        candidate.mesh().indices.len() / 3,
        candidate.work_used(),
        candidate.retained_budget_bytes()
    );
    let tail_reserve = old.mesh().indices.len() / 3 + 5 * cells.cap_edge_count();
    assert!(matches!(
        cut.mesh_candidate(
            Tessellation {
                triangle_limit: tail_reserve - 1,
                ..q
            },
            512_000_000,
            64 * 1024 * 1024,
            None
        ),
        Err(Error::Limit)
    ));
    for (work, bytes) in [
        (candidate.work_used() - 1, 64 * 1024 * 1024),
        (512_000_000, candidate.retained_budget_bytes() - 1),
    ] {
        assert!(matches!(
            cut.mesh_candidate(q, work, bytes, None),
            Err(Error::Limit)
        ));
    }
    let callbacks = std::cell::Cell::new(0);
    let observe = || {
        callbacks.set(callbacks.get() + 1);
        false
    };
    let expected = cut
        .mesh_candidate(q, 512_000_000, 64 * 1024 * 1024, Some(&observe))
        .unwrap();
    let last = callbacks.get();
    callbacks.set(0);
    let cancel = || {
        callbacks.set(callbacks.get() + 1);
        callbacks.get() >= last
    };
    assert!(matches!(
        cut.mesh_candidate(q, 512_000_000, 64 * 1024 * 1024, Some(&cancel)),
        Err(Error::Cancelled)
    ));
    let retry = cut
        .mesh_candidate(q, 512_000_000, 64 * 1024 * 1024, None)
        .unwrap();
    assert_eq!(retry.mesh(), expected.mesh());
    assert_eq!(retry.connectivity(), expected.connectivity());
    assert_eq!(cut.thread().source().spec(), before);
}

#[test]
fn genuine_stock_mesh_rejects_invalid_resources_unresolvable_quality_and_too_thin_bands() {
    let cut = engineering_cut(0.01);
    let q = quality();
    for (work, bytes) in [
        (0, 64 * 1024 * 1024),
        (512_000_001, 64 * 1024 * 1024),
        (512_000_000, 0),
        (512_000_000, 64 * 1024 * 1024 + 1),
    ] {
        assert!(matches!(
            cut.mesh_candidate(q, work, bytes, None),
            Err(Error::Invalid)
        ));
    }
    assert!(matches!(
        cut.mesh_candidate(q, cut.work + 1024, 64 * 1024 * 1024, None),
        Err(Error::Limit)
    ));
    assert!(matches!(
        cut.mesh_candidate(q, 512_000_000, 4095, None),
        Err(Error::Limit)
    ));
    assert!(matches!(
        cut.mesh_candidate(
            Tessellation {
                deflection: 1e-12,
                ..q
            },
            512_000_000,
            64 * 1024 * 1024,
            None
        ),
        Err(Error::Limit)
    ));
    for flat in [0.001, 0.0001] {
        let thin = engineering_cut(flat);
        let cells = thin
            .thread()
            .cell_plan(128, 32, 512_000_000, 64 * 1024 * 1024, None)
            .unwrap();
        match cells.mesh_candidate(q, 512_000_000, 64 * 1024 * 1024, None) {
            Ok(old) => {
                let pairs = old
                    .mesh()
                    .indices
                    .as_chunks::<3>()
                    .0
                    .iter()
                    .zip(&old.mesh().face_indices)
                    .filter(|(_, face)| **face == 6)
                    .filter(|(ids, _)| {
                        matches!(
                            cells.vertex_owner(ids[1] as usize),
                            Some(StockClipPointOwner::Vertex(_))
                        ) && matches!(
                            cells.vertex_owner(ids[2] as usize),
                            Some(StockClipPointOwner::Vertex(_))
                        )
                    })
                    .count();
                eprintln!(
                    "THIN_REFUSAL_STAGE flat={flat} original_candidate=PASS end_vertex_pairs={pairs}"
                );
            }
            Err(error) => eprintln!("THIN_REFUSAL_STAGE flat={flat} original_candidate={error:?}"),
        }
        assert!(matches!(
            thin.mesh_candidate(q, 512_000_000, 64 * 1024 * 1024, None),
            Err(Error::Unresolved)
        ));
        assert_eq!(
            thin.thread().source().spec().crest_radius_limit_mm,
            Some(5.1)
        );
    }
}

#[test]
fn genuine_exact_stock_precision_and_large_radial_transfer_refuse_without_changing_source() {
    let base = engineering_cut(0.01);
    let zero = Rational::new(0, 1).unwrap();
    let far = Rational::new(1_i64 << 55, 1).unwrap();
    for (label, center, bottom, radius) in [
        (
            "radial-transfer",
            [zero; 2],
            zero,
            Rational::new(51, 5).unwrap(),
        ),
        (
            "unresolvable-join",
            [zero; 2],
            far,
            Rational::new(51, 10).unwrap(),
        ),
        (
            "unresolvable-center",
            [far, zero],
            zero,
            Rational::new(51, 10).unwrap(),
        ),
    ] {
        let stock = crate::cylinder::AxialCylinder::new(
            73,
            crate::Circle::new(center, radius).unwrap(),
            bottom,
            Rational::new(7, 1).unwrap(),
        )
        .unwrap();
        let cut = ExternalStockCut::new(*base.thread(), stock, Default::default(), None).unwrap();
        assert_eq!(cut.stock().circle().center(), center);
        assert_eq!(cut.stock().bottom(), bottom);
        assert_eq!(cut.stock().circle().radius(), radius);
        let before = cut.thread().source().spec();
        assert!(
            matches!(
                cut.mesh_candidate(quality(), 512_000_000, 64 * 1024 * 1024, None),
                Err(Error::Unresolved)
            ),
            "{label}"
        );
        assert_eq!(cut.thread().source().spec(), before);
        eprintln!("EXACT_STOCK_PUBLIC_REFUSAL {label}");
    }
    assert!(
        base.mesh_candidate(quality(), 512_000_000, 64 * 1024 * 1024, None)
            .is_ok()
    );
}
