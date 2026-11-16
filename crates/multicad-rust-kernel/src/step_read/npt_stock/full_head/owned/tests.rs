use super::*;
use num_bigint::BigInt;
use num_rational::BigRational as Q;
use num_traits::{One, Signed, Zero};
fn q(x: f64) -> Q {
    Q::from_float(x).unwrap()
}
fn decimal(s: &str) -> Q {
    let (m, e) = s.split_once('E').unwrap_or((s, "0"));
    let exponent = e.parse::<i32>().unwrap();
    let fraction = m.split_once('.').map_or(0, |(_, v)| v.len()) as i32;
    let n = m.replace('.', "").parse::<BigInt>().unwrap();
    let power = exponent - fraction;
    if power >= 0 {
        Q::from_integer(n * BigInt::from(10).pow(power as u32))
    } else {
        Q::new(n, BigInt::from(10).pow((-power) as u32))
    }
}
fn contains(bounds: [f64; 2], value: &Q) {
    assert!(
        q(bounds[0]) <= *value && *value <= q(bounds[1]),
        "{bounds:?} excludes {value}"
    );
}
fn candidate(quarter: bool) -> String {
    let input = super::super::super::tests::fixture(quarter);
    let doc = CompoundDocument::parse(&input, DocumentOptions::default(), None).unwrap();
    let parsed = Parsed::parse(
        &doc,
        doc.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    let cap = parsed
        .faces()
        .iter()
        .filter_map(|f| match f.support() {
            Support::Plane { height } => {
                Some((parsed.points()[height].coordinates()[2].value(), f.entity()))
            }
            _ => None,
        })
        .max_by(|a, b| a.0.total_cmp(&b.0))
        .unwrap()
        .1;
    let head = HeadInput::regular_hex(&parsed, cap, 18., 3.).unwrap();
    let proof = head::certify(&parsed, head, 1e-7, 1e-6, Default::default(), None).unwrap();
    crate::step::npt_stock::full_head::propose(&proof, Default::default(), None)
        .unwrap()
        .as_str()
        .to_owned()
}
#[test]
fn owned_actual_head_survives_all_original_arenas_and_has_actual_roles_and_sewn_measure() {
    for quarter in [false, true] {
        let owned = {
            let input: Arc<str> = Arc::from(candidate(quarter));
            OwnedActualHead::read(Arc::clone(&input), Default::default(), None).unwrap()
        };
        // Every original input/CompoundDocument/Selected/Parsed/head proof has
        // already been destroyed. Only actual Arc records and earned owners live.
        let doc =
            CompoundDocument::parse(owned.actual_file(), DocumentOptions::default(), None).unwrap();
        assert_eq!(owned.original_record_count(), doc.records().len());
        for index in 0..owned.original_record_count() {
            let id = owned.original_entity_at(index).unwrap();
            for (component, e) in doc.records()[&id].iter().enumerate() {
                assert_eq!(owned.original_component(id, component), e.source);
            }
        }
        let roles = owned.physical_roles();
        assert_eq!(
            doc.records()[&roles.solid_entity][0].name,
            "MANIFOLD_SOLID_BREP"
        );
        let annulus = &doc.records()[&roles.annulus_face_entity][0];
        assert!(matches!(&annulus.args[1],Value::List(v) if v.len()==2));
        assert!(
            !doc.records()
                .contains_key(&owned.head().auxiliary_cap_entity)
        );
        for (i, &id) in roles.head_edge_entities.iter().enumerate() {
            let e = &doc.records()[&id][0];
            assert_eq!(e.name, "EDGE_CURVE");
            let [a, z] = owned.head().edges[i];
            assert!(matches!(e.args[1],Value::Reference(v) if v==roles.head_vertex_entities[a]));
            assert!(matches!(e.args[2],Value::Reference(v) if v==roles.head_vertex_entities[z]));
        }
        drop(doc);
        let thread = owned.auxiliary_thread();
        let mut caps = Vec::new();
        for f in thread.faces() {
            if let Support::Plane { height } = f.support() {
                caps.push(decimal(
                    thread.point(height).unwrap().coordinates()[2].spelling(),
                ));
            }
        }
        caps.sort();
        let lower = caps[0].clone();
        let upper = caps[1].clone();
        let transfer = owned.axial_sewing_transfer();
        let scale = (q(owned.head_stations_mm()[0]) - &lower) / (&upper - &lower);
        assert!(scale > Q::zero());
        contains(transfer.positive_axial_scale_bounds(), &scale);
        // Independently enclose affine displacement at BOTH actual hull
        // endpoints. Convexity then bounds every interior raw image point,
        // including controls slightly outside the nominal cap interval.
        for z in thread.local_enclosure_mm()[2] {
            let delta = ((&scale - Q::one()) * (q(z) - &lower)).abs();
            assert!(delta <= q(transfer.domain_displacement_mm()));
        }
        let mut area = Q::zero();
        let poly = owned.head_polygon_mm();
        for i in 0..6 {
            let a = poly[i];
            let b = poly[(i + 1) % 6];
            area += q(a[0]) * q(b[1]) - q(a[1]) * q(b[0]);
        }
        area /= Q::from_integer(BigInt::from(2));
        let prism = area * (q(owned.head_stations_mm()[1]) - q(owned.head_stations_mm()[0]));
        let raw = thread.local_volume_mm3();
        contains(
            owned.local_sewn_volume_mm3(),
            &(q(raw[0]) * &scale + &prism),
        );
        contains(
            owned.local_sewn_volume_mm3(),
            &(q(raw[1]) * &scale + &prism),
        );
        // Exact decimal rational controls, not cached values, independently
        // verify the published affine transfer and whole sewn support bounds.
        let bounds = owned.local_sewn_enclosure_mm();
        for tensor in thread.tensors() {
            for &index in tensor.point_indices() {
                let coordinates = thread.point(index).unwrap().coordinates();
                let point = coordinates.map(|s| decimal(s.spelling()));
                contains(bounds[0], &point[0]);
                contains(bounds[1], &point[1]);
                let mapped_z = &lower + (&point[2] - &lower) * &scale;
                contains(bounds[2], &mapped_z);
                assert!((mapped_z - &point[2]).abs() <= q(transfer.domain_displacement_mm()));
            }
        }
        let v = [0.25, -0.75, 1.25];
        let boxes = v.map(|x| [x, x]);
        let mapped = transfer.map_raw_vector_bounds(boxes).unwrap();
        let normal = transfer.map_raw_normal_bounds(boxes).unwrap();
        for i in 0..3 {
            contains(
                mapped[i],
                &(q(v[i]) * if i == 2 { scale.clone() } else { Q::one() }),
            );
            contains(
                normal[i],
                &(q(v[i]) * if i < 2 { scale.clone() } else { Q::one() }),
            );
        }
        assert!(owned.minimum_wall_mm() > 0. && owned.join_sewing_mm() <= 1e-7 / 8.);
        eprintln!(
            "OWNED_ACTUAL_HEAD quarter={quarter} {:?} volume={:?}",
            owned.receipt(),
            owned.local_sewn_volume_mm3()
        );
    }
}

#[test]
fn owning_actual_head_work_storage_and_cancel_frontiers_preserve_immutable_input() {
    use std::cell::Cell;
    let input: Arc<str> = Arc::from(candidate(false));
    let good = OwnedActualHead::read(Arc::clone(&input), Default::default(), None).unwrap();
    let r = good.receipt();
    for maximum_work in [r.aggregate_work - 1, r.aggregate_work] {
        let options = OwnedHeadOptions {
            resources: Options {
                maximum_work,
                ..Default::default()
            },
            ..Default::default()
        };
        let value = OwnedActualHead::read(Arc::clone(&input), options, None);
        if maximum_work < r.aggregate_work {
            assert!(matches!(value, Err(Error::Limit)));
        } else {
            assert_eq!(
                value.unwrap().local_sewn_volume_mm3(),
                good.local_sewn_volume_mm3()
            );
        }
    }
    for maximum_bytes in [r.accounted_peak_bytes - 1, r.accounted_peak_bytes] {
        let options = OwnedHeadOptions {
            resources: Options {
                maximum_bytes,
                ..Default::default()
            },
            ..Default::default()
        };
        let value = OwnedActualHead::read(Arc::clone(&input), options, None);
        if maximum_bytes < r.accounted_peak_bytes {
            assert!(matches!(value, Err(Error::Limit)));
        } else {
            assert_eq!(value.unwrap().receipt().retained_bytes, r.retained_bytes);
        }
    }
    let calls = Cell::new(0usize);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    OwnedActualHead::read(Arc::clone(&input), Default::default(), Some(&count)).unwrap();
    let total = calls.get();
    for at in [1, total / 2, total] {
        let calls = Cell::new(0usize);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() >= at
        };
        assert!(matches!(
            OwnedActualHead::read(Arc::clone(&input), Default::default(), Some(&cancel)),
            Err(Error::Cancelled)
        ));
        assert_eq!(
            OwnedActualHead::read(Arc::clone(&input), Default::default(), None)
                .unwrap()
                .receipt()
                .aggregate_work,
            r.aggregate_work
        );
    }
    let other = OwnedActualHead::read(
        Arc::clone(&input),
        OwnedHeadOptions {
            resources: Options {
                additional_resident_bytes: 4096,
                ..Default::default()
            },
            ..Default::default()
        },
        None,
    )
    .unwrap();
    assert_eq!(other.receipt().aggregate_work, r.aggregate_work);
    assert_eq!(
        other.receipt().accounted_peak_bytes,
        r.accounted_peak_bytes + 4096
    );
    assert_eq!(good.actual_file(), input.as_ref());
}
