use super::*;
use std::hash::{Hash, Hasher};
fn source(left: bool, nonbinary: bool) -> Insert {
    let original = super::super::super::tests::fixture(left, 1, nonbinary);
    let mut outer = original.outer().spec();
    outer.crest_width_mm = 0.4;
    outer.lower_flank_deg = 30.;
    outer.upper_flank_deg = 30.;
    let mut inner = original.inner().spec();
    inner.left_hand = left;
    let options = crate::authored_flat_thread::SourceOptions {
        position_limit_mm: 1e-7,
        normal_angle_limit_rad: 1e-7,
    };
    Insert::new(
        AuthoredThread::new(outer, options).unwrap(),
        AuthoredThread::new(inner, options).unwrap(),
        0.5,
        None,
    )
    .unwrap()
}
fn options() -> Tessellation {
    Tessellation {
        deflection: 0.1,
        angle: 0.3,
        triangle_limit: 1_000_000,
    }
}
#[test]
fn actual_default_insert_candidate_has_closed_annular_material_and_paired_volume() {
    let body = source(false, false);
    let before = body;
    let time = std::time::Instant::now();
    let c = body.tessellate_candidate(options(), None).unwrap();
    eprintln!(
        "actual wood insert {}facets {:?} {}work source{:?} literal{:?}",
        c.mesh().face_indices.len(),
        time.elapsed(),
        c.work_used(),
        c.source_volume_mm3(),
        c.literal_volume_mm3()
    );
    let mut hash = std::collections::hash_map::DefaultHasher::new();
    for v in &c.mesh().vertices {
        v.to_bits().hash(&mut hash);
    }
    c.mesh().indices.hash(&mut hash);
    c.mesh().face_indices.hash(&mut hash);
    eprintln!("default actualmesh hash {:016x}", hash.finish());
    assert_eq!(*c.source(), before);
    assert!(c.source_chord_bound_mm() <= options().deflection);
    assert!(c.source_normal_angle_bound_rad() <= options().angle);
    let proof = crate::manufacturing::certify(c.mesh()).unwrap();
    assert!(proof.edge_pairs > 0 && proof.vertex_links > 0 && proof.signed_volume_mm3 > 0.);
    let v = c.literal_volume_mm3();
    assert!(v[0] <= proof.signed_volume_mm3 && proof.signed_volume_mm3 <= v[1]);
    let source = c.source_volume_mm3();
    assert!((proof.signed_volume_mm3 - source[0]).abs() < source[0] * 0.01);
    assert_eq!(body, before);
}

#[test]
fn actual_left_hand_and_nonbinary_insert_materials_keep_source_and_closed_manifold() {
    for (left, nonbinary) in [(true, false), (false, true)] {
        let body = source(left, nonbinary);
        let timer = std::time::Instant::now();
        let c = body.tessellate_candidate(options(), None).unwrap();
        let proof = crate::manufacturing::certify(c.mesh()).unwrap();
        let literal = c.literal_volume_mm3();
        assert!(literal[0] <= proof.signed_volume_mm3 && proof.signed_volume_mm3 <= literal[1]);
        assert_eq!(*c.source(), body);
        assert!(c.source_chord_bound_mm() <= 0.1 && c.source_normal_angle_bound_rad() <= 0.3);
        eprintln!(
            "insert LH{left}/nonbinary{nonbinary}: {}facets {:?} {}work",
            c.mesh().face_indices.len(),
            timer.elapsed(),
            c.work_used()
        );
    }
}
