use super::*;
use crate::HelicalThread;
use multicad_rust_kernel::{
    authored_flat_thread::{AuthoredThread, SourceOptions},
    headed_screw::{ScrewHead, ScrewSpec},
    helical_thread::Spec,
    wood_insert::Insert,
};
fn fixture(left: bool, length: f64, plain: f64, frame: PlaneFrame) -> (Solid, Solid) {
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
