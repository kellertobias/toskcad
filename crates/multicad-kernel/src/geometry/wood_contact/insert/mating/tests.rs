use super::*;
use crate::HelicalThread;
use multicad_rust_kernel::{
    authored_flat_thread::{AuthoredThread, SourceOptions},
    headed_screw::{ScrewHead, ScrewSpec},
    helical_thread::Spec,
    wood_insert::Insert,
};
pub(in crate::geometry::wood_contact) fn fixture(
    left: bool,
    length: f64,
    plain: f64,
    frame: PlaneFrame,
) -> (Solid, Solid) {
    let p = 0.5;
    let major = (2.874 + 2.980) / 4.;
    let pitchmid = (2.580 + 2.655) / 2.;
    let crest = p / 2. - (2. * major - pitchmid) / 3_f64.sqrt();
    let r = p / (4. * 3_f64.sqrt());
    let source = |major, pitch, depth, crest| {
        AuthoredThread::new(
            Spec {
                major_radius_mm: major,
                pitch_mm: pitch,
                length_mm: 4.,
                depth_mm: depth,
                crest_width_mm: crest,
                lower_flank_deg: 30.,
                upper_flank_deg: 30.,
                starts: 1,
                left_hand: left,
                root_radius_mm: None,
            },
            SourceOptions {
                position_limit_mm: 1e-8,
                normal_angle_limit_rad: 1e-6,
            },
        )
        .unwrap()
    };
    let insert = Insert::new(
        source(5., 2., 1., 0.4),
        source(
            1.5,
            p,
            (3. - (2.459 + 2.599) / 2.) / 2.,
            p / 2. - (3. - (2.675 + 2.775) / 2.) / 3_f64.sqrt(),
        ),
        0.5,
        None,
    )
    .unwrap();
    let policy = GeometryPolicy::default();
    let insert = Solid::rust_wood_insert(insert, EntityId::new(), frame, policy, None).unwrap();
    let screw = Solid::rust_machine_screw_at_thread_frame(
        HelicalThread {
            major_radius_mm: major,
            pitch_mm: p,
            length_mm: length,
            depth_mm: (p - crest) / (2. / 3_f64.sqrt()) - r,
            crest_width_mm: crest,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            radial_taper: 0.,
            starts: 1,
            left_hand: left,
        },
        Some(r),
        ScrewSpec {
            head: ScrewHead::SocketCap {
                radius_mm: 2.75,
                height_mm: 3.,
                socket_across_flats_mm: 2.,
                socket_depth_mm: 1.,
            },
            shank_length_mm: plain,
            shank_radius_mm: 1.5,
        },
        EntityId::new(),
        frame,
        policy,
    )
    .unwrap();
    (screw, insert)
}
fn identity() -> PlaneFrame {
    PlaneFrame {
        origin: [0.; 3],
        u: [1., 0., 0.],
        v: [0., 1., 0.],
    }
}
#[test]
fn current_native_m3_screw_sources_mate_actual_6h_insert_without_recipe_aliases() {
    for left in [false, true] {
        let (screw, insert) = fixture(left, 3., 0., identity());
        let proof = screw.rust_wood_insert_mating(&insert, 2., None).unwrap();
        assert!(proof.engagement_mm[0] <= 3. && proof.engagement_mm[1] >= 3.);
        assert!(proof.thread_exit_clearance_mm[0] > 0.99);
        assert!(proof.minimum_radial_clearance_lower_mm > 0.005);
    }
}
#[test]
fn finite_continuation_and_source_frame_failure_controls_keep_original_geometry() {
    let (screw, insert) = fixture(false, 4.5, 0., identity());
    let proof = screw.rust_wood_insert_mating(&insert, 2., None).unwrap();
    assert!(!proof.entire_thread_within_insert);
    assert!(proof.thread_exit_clearance_mm[1] < -0.49);
    assert!(proof.engagement_mm[0] <= 4. && proof.engagement_mm[1] >= 4.);
    assert_eq!(
        screw
            .rust_wood_insert_mating(&insert, 2., Some(&|| true))
            .unwrap_err()
            .kind,
        GeometryErrorKind::Cancelled
    );
    assert_eq!(
        screw
            .rust_wood_insert_mating(&insert, 0., None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Input
    );
    let plain = Solid::box_solid([20.; 3], GeometryPolicy::default()).unwrap();
    assert_eq!(
        plain
            .rust_wood_insert_mating(&insert, 2., None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        screw
            .rust_wood_insert_mating(&plain, 2., None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
    let mut policy = GeometryPolicy::default();
    policy.linear_tolerance_mm *= 2.;
    let unequal = Solid::rust_wood_insert(
        insert.rust_wood_insert_data().unwrap().source,
        EntityId::new(),
        identity(),
        policy,
        None,
    )
    .unwrap();
    assert_eq!(
        screw
            .rust_wood_insert_mating(&unequal, 2., None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Input
    );
    let quarter = PlaneFrame {
        origin: [0.; 3],
        u: [0., 1., 0.],
        v: [-1., 0., 0.],
    };
    let (wrong, _) = fixture(false, 3., 0., quarter);
    assert_eq!(
        wrong
            .rust_wood_insert_mating(&insert, 2., None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
    let mut penetrating = identity();
    penetrating.origin[2] = 0.1;
    let (wrong, _) = fixture(false, 3., 1., penetrating);
    assert_eq!(
        wrong
            .rust_wood_insert_mating(&insert, 2., None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
    let world = PlaneFrame {
        origin: [105., -20., 59.],
        u: [0.8, 0.6, 0.],
        v: [0.36, -0.48, -0.8],
    };
    let (posed, receiver) = fixture(true, 3., 0., world);
    assert!(posed.rust_wood_insert_mating(&receiver, 2., None).is_ok());
}
#[test]
fn actual_profiled_and_phillips_native_sources_preserve_their_complete_prefixes() {
    use multicad_rust_kernel::{
        head_profile::{Head, Profile},
        phillips_drive::{self, Cavity, Cutter, PhillipsScrew},
        profiled_screw,
    };
    let (screw, insert) = fixture(false, 3., 0., identity());
    let thread = *screw.rust_machine_screw_geometry().unwrap().0.thread();
    let head = Head::Countersunk {
        radius_mm: 2.75,
        height_mm: 1.65,
        included_angle_deg: 90.,
        edge_chamfer_mm: 0.05,
        neck_radius_mm: 0.1,
    };
    let spec = profiled_screw::Spec {
        head,
        shank_length_mm: 1.,
        socket_across_flats_mm: 2.,
        socket_depth_mm: 0.6,
        minimum_wall_mm: 1e-5,
    };
    let actual = profiled_screw::Screw::new_with_shank_radius(thread, 1.5, spec).unwrap();
    let mut above = identity();
    above.origin[2] = -0.5;
    let profiled = Solid::rust_profiled_screw_from_source(
        actual,
        EntityId::new(),
        above,
        GeometryPolicy::default(),
    )
    .unwrap();
    assert!(profiled.rust_wood_insert_mating(&insert, 2., None).is_ok());
    let s = phillips_drive::Spec {
        mouth_diameter_mm: 3.2,
        mouth_radius_mm: 0.5,
        b_width_mm: 0.955,
        g_diameter_mm: 1.27,
        tilt_angle_deg: 7.,
        outer_flank_angle_deg: 26.5,
        floor_angle_deg: 28.,
        outer_corner_angle_deg: 92.,
        inner: phillips_drive::Inner::ShallowV {
            angle_deg: 138.,
            e_width_mm: 0.435,
        },
    };
    let tip = (s.mouth_diameter_mm - s.g_diameter_mm) / (2. * 26.5_f64.to_radians().tan())
        + s.g_diameter_mm / 2. * 28_f64.to_radians().tan();
    let domain = phillips_drive::Domain {
        corner_start_depth_mm: -6.4,
        corner_end_depth_mm: tip + 3.2,
        corner_reach_mm: 4. * (2.75 + 3.2 + (tip - 1.65).abs() + 1.),
    };
    let cavity = Cavity::new(Cutter::new(s, domain).unwrap(), 1e-5, 100_000).unwrap();
    let actual = PhillipsScrew::new(
        thread,
        Profile::new(head, 1.5, 1.).unwrap(),
        profiled_screw::SourceSpec {
            head: Some(head),
            shank_length_mm: 1.,
            minimum_wall_mm: 1e-5,
        },
        cavity,
        100_000,
        &|| false,
    )
    .unwrap();
    let phillips = Solid::rust_phillips_screw_from_source(
        actual,
        EntityId::new(),
        above,
        GeometryPolicy::default(),
    )
    .unwrap();
    assert!(phillips.rust_wood_insert_mating(&insert, 2., None).is_ok());
}
