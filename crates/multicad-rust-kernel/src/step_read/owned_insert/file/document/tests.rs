use super::*;
use crate::step::{
    Writer,
    insert::{CurveOptions, DomainOptions, DomainPlan, SurfacePlan, TopologyPlan},
};
fn mixed(left: bool) -> String {
    let body = crate::step_read::owned_insert::tests::graph(left, 4.);
    let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
    let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
    let t = TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    let gear = gear();
    let mut w = Writer::new();
    w.append_herringbone(&gear, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
        .unwrap();
    w.append_owned_insert_mapped_candidate(
        &t,
        [31., -17., 9.],
        [0.8, 0.6, 0.],
        [-0.36, 0.48, 0.8],
        16_000_000,
        64 * 1024 * 1024,
        None,
    )
    .unwrap();
    w.finish().unwrap()
}
#[test]
fn actual_twenty_four_tooth_compound_has_a_separate_global_and_selected_budget() {
    for left in [false, true] {
        let text = mixed(left);
        assert!(matches!(
            File::parse(&text, Options::default(), None),
            Err(Error::Limit)
        ));
        let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
        eprintln!(
            "DOCUMENT input={} work={} bytes={}",
            text.len(),
            document.work_used(),
            document.retained_budget_bytes()
        );
        assert!(document.retained_budget_bytes() > 64 * 1024 * 1024);
        assert!(document.work_used() < 16_000_000); // Proves old rejection was retained bytes, not work.
        assert_eq!(document.literal_input(), text);
        let member = document.represented_members()[1];
        assert!(member.mapping_link().is_some());
        let packet = document
            .member_packet(member.item_entity(), Options::default(), None)
            .unwrap();
        assert!(std::ptr::eq(packet.document(), &document));
        assert_eq!(packet.parsed().mapping_link(), member.mapping_link());
        assert_eq!(packet.parsed().literal_input(), text);
        assert!(packet.query_budget_bytes() < 64 * 1024 * 1024);
        let parsed = packet.parsed();
        let bindings = parsed.support_bindings(Options::default(), None).unwrap();
        let material = bindings.material_bounds(Options::default(), None).unwrap();
        let volume = material
            .occupied_volume(FluxOptions::default(), None)
            .unwrap();
        assert!(volume.volume_mm3()[0] > 0.);
        let link = member.mapping_link().unwrap();
        let mapping = parsed
            .physical_mapping(link.operator_entity(), Options::default(), None)
            .unwrap();
        assert_eq!(mapping.handedness().unwrap(), 1);
        let view = document
            .local_member_document(
                document.represented_members()[0].item_entity(),
                Options::default(),
                None,
            )
            .unwrap();
        assert!(view.query_budget_bytes() < 64 * 1024 * 1024);
        assert!(view.work_used() < 16_000_000);
        let ordinary = crate::step_read::verify_herringbone(
            view.as_str(),
            &gear(),
            [0.; 3],
            [1., 0., 0.],
            [0., 1., 0.],
            1e-7,
        );
        assert!(ordinary.is_ok(), "{ordinary:?}");
    }
}
fn gear() -> crate::spur_gear::herringbone::Herringbone {
    let beta = 15_f64.to_radians();
    let pressure = (20_f64.to_radians().tan() / beta.cos()).atan();
    let module = 2. / beta.cos();
    let pitch = module * 24. / 2.;
    crate::spur_gear::herringbone::Herringbone::new_helical(
        crate::spur_gear::Spec {
            teeth: 24,
            pitch_radius_mm: pitch,
            base_radius_mm: pitch * pressure.cos(),
            tip_radius_mm: pitch + 2.,
            root_radius_mm: pitch - 2.5,
            root_fillet_radius_mm: 0.4,
            half_base_angle_rad: (std::f64::consts::PI * module / 2. - 0.05) / (2. * pitch)
                + pressure.tan()
                - pressure,
            face_width_mm: 4.,
            bore_radius_mm: 4.,
        },
        beta.tan() / pitch,
        1e-5,
        1e7,
    )
    .unwrap()
}
