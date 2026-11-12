use super::*;
use crate::hexalobular_drive::{Cone, TangentSpec};

fn source() -> Drive {
    Drive::from_authored_tangent(
        TangentSpec {
            major_diameter_mm: 2.8,
            minor_diameter_mm: 2.05,
            outer_lobe_radius_mm: 0.28,
        },
        Cone {
            tip_depth_mm: 1.1 + 1.4 / 3_f64.sqrt(),
            depth_per_radius: 1. / 3_f64.sqrt(),
        },
    )
    .unwrap()
}
#[test]
fn torx_encoder_retains_real_surfaces_contacts_and_apex_without_opening_public_admission() {
    let source = source();
    let p = Placement::new([10., -20., 30.], [0.8, 0.6, 0.], [0., 0., 1.]).unwrap();
    let mut file = StepFile::new();
    let (_mouth, faces) = encode(&mut file, &p, &source, -2.4, 1e-7).unwrap();
    assert_eq!(faces.len(), 13);
    assert_eq!(file.data.matches("CYLINDRICAL_SURFACE(").count(), 12);
    assert_eq!(file.data.matches("CONICAL_SURFACE(").count(), 1);
    assert_eq!(file.data.matches("CIRCLE(").count(), 12);
    assert_eq!(file.data.matches("INTERSECTION_CURVE(").count(), 12);
    assert_eq!(file.data.matches("B_SPLINE_CURVE_WITH_KNOTS(").count(), 12);
    assert_eq!(file.data.matches("VERTEX_POINT(").count(), 25);
    assert_eq!(file.data.matches("EDGE_CURVE(").count(), 37);
    assert_eq!(file.data.matches("EVALUATED_DEGENERATE_PCURVE(").count(), 1);
    assert_eq!(file.data.matches("SEAM_CURVE(").count(), 1);
    assert_eq!(file.data.matches("SEAM_EDGE(").count(), 2);
    for arc in source.arcs() {
        let xy = arc.start_mm;
        for z in [
            -2.4,
            -2.4 + source.cone().tip_depth_mm - source.cone().depth_per_radius * xy[0].hypot(xy[1]),
        ] {
            let point = p.point_f64(xy, z).map(|v| if v == 0. { 0. } else { v });
            assert!(file.points.contains_key(&point.map(f64::to_bits)));
        }
    }
}

#[test]
fn retained_pan_head_bundle_uses_curved_cavity_while_public_screw_export_stays_gated() {
    let thread = Thread::new(crate::helical_thread::Spec {
        major_radius_mm: 1.5,
        pitch_mm: 0.5,
        length_mm: 3.,
        depth_mm: 0.2,
        crest_width_mm: 0.0625,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: false,
        root_radius_mm: None,
    })
    .unwrap();
    let recess = crate::screw_recess::Recess::conical_hexalobular(source(), 1e-7, 100_000).unwrap();
    let screw = Screw::new_with_recess(
        thread,
        1.5,
        crate::head_profile::Head::Crowned {
            radius_mm: 2.8,
            height_mm: 2.4,
            crown_radius_mm: 5.,
            top_land_radius_mm: 2.05,
            neck_radius_mm: 0.1,
        },
        1.,
        0.1,
        recess,
    )
    .unwrap();
    let p = Placement::new([0.; 3], [1., 0., 0.], [0., 1., 0.]).unwrap();
    let mut file = StepFile::new();
    let rim = CircularRim::new(&mut file, [0.; 3], p.normal, p.u, 1.5).unwrap();
    let join = Junction {
        crests: vec![],
        complements: vec![],
        vertex: rim.vertex,
        reference: p.u,
        positive: vec![(rim.edge, true)],
        outer: Some(rim),
    };
    let faces = head(&mut file, &p, &screw, join, 1e-7).unwrap();
    assert_eq!(faces.len(), screw.profile().segments().count() + 14);
    assert_eq!(file.data.matches("INTERSECTION_CURVE(").count(), 12);
    assert!(!crate::step::supports_profiled_screw(&screw, [0.; 3], 1e-7));
    assert_eq!(
        crate::step::write_profiled_screw(&screw, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7),
        Err(StepError::Unsupported)
    );
}

#[test]
fn torx_encoder_refuses_unresolvable_precision_and_collapsed_world_geometry() {
    let source = source();
    let snapshot = source.clone();
    let p = Placement::new([0.; 3], [1., 0., 0.], [0., 1., 0.]).unwrap();
    for tolerance in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(matches!(
            encode(&mut StepFile::new(), &p, &source, -2.4, tolerance),
            Err(StepError::InvalidGeometry)
        ));
    }
    assert!(matches!(
        encode(&mut StepFile::new(), &p, &source, f64::NAN, 1e-7),
        Err(StepError::InvalidGeometry)
    ));
    assert!(matches!(
        encode(&mut StepFile::new(), &p, &source, -2.4, 1e-20),
        Err(StepError::Unsupported)
    ));
    let shallow = Drive::from_authored_tangent(
        TangentSpec {
            major_diameter_mm: 2.8,
            minor_diameter_mm: 2.05,
            outer_lobe_radius_mm: 0.28,
        },
        Cone {
            tip_depth_mm: 1.,
            depth_per_radius: f64::MIN_POSITIVE,
        },
    )
    .unwrap();
    assert!(matches!(
        encode(&mut StepFile::new(), &p, &shallow, -2.4, 1e-7),
        Err(StepError::Unsupported)
    ));
    let collapsed = Placement::new([1e300; 3], [1., 0., 0.], [0., 1., 0.]).unwrap();
    assert!(matches!(
        encode(&mut StepFile::new(), &collapsed, &source, -2.4, 1e-7),
        Err(StepError::InvalidGeometry)
    ));
    let skewed = Placement {
        origin: [0.; 3],
        u: [0.5, 0., 0.],
        v: [0., 1., 0.],
        normal: [0., 0., 1.],
    };
    assert!(matches!(
        encode(&mut StepFile::new(), &skewed, &source, -2.4, 1e-7),
        Err(StepError::InvalidGeometry)
    ));
    for error in [
        Error::Invalid,
        Error::Unresolved,
        Error::Limit,
        Error::Cancelled,
    ] {
        let actual = source_error(error);
        assert_eq!(
            actual,
            match error {
                Error::Invalid => StepError::InvalidGeometry,
                Error::Unresolved => StepError::Unsupported,
                Error::Limit | Error::Cancelled => StepError::Limit,
            }
        );
    }
    assert_eq!(source, snapshot);
    assert!(encode(&mut StepFile::new(), &p, &source, -2.4, 1e-7).is_ok());
}
