use super::*;
fn owned_source(left: bool, extra: f64) -> OwnedPeriodicMaterial {
    let original = degree11::tests::fixture::actual_source();
    let thread = crate::helical_thread::Thread::new(crate::helical_thread::Spec {
        left_hand: left,
        length_mm: original.thread().spec().length_mm + extra,
        ..original.thread().spec()
    })
    .unwrap();
    let source = crate::phillips_drive::PhillipsScrew::new(
        thread,
        *original.profile(),
        original.spec(),
        original.cavity().clone(),
        1_000_000,
        &|| false,
    )
    .unwrap();
    OwnedPeriodicMaterial::from_authored(&source, Default::default(), None).unwrap()
}
fn file(source: &OwnedPeriodicMaterial, posed: bool) -> (String, [[f64; 3]; 4]) {
    let proposal = source
        .authored_mesh_approximation(16_000_000, None)
        .unwrap();
    let (o, u, v) = if posed {
        ([105., -20., 60.], [0.8, 0.6, 0.], [0.36, -0.48, -0.8])
    } else {
        ([0.; 3], [1., 0., 0.], [0., 1., 0.])
    };
    let normal = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    (
        crate::step::write_phillips_screw(proposal.proposal(), o, u, v, 1e-7).unwrap(),
        [o, u, v, normal],
    )
}
#[test]
fn actual_metric_files_match_complete_boundary_caps_world_maps_and_independent_mass() {
    for (left, posed) in [(false, false), (true, true)] {
        let source = owned_source(left, 0.);
        let (text, occurrence) = file(&source, posed);
        let decoded =
            crate::step_read::read_periodic_phillips_material(&text, Default::default(), None)
                .unwrap();
        let now = std::time::Instant::now();
        let proof = source
            .match_material(decoded.material(), occurrence, Default::default(), None)
            .unwrap();
        assert!(
            std::ptr::eq(proof.source(), &source)
                && std::ptr::eq(proof.target(), decoded.material())
        );
        assert!(
            proof.maximum_world_distance_mm() <= 1e-7
                && proof.maximum_world_chart_derivative_mm() <= 1e-7
        );
        assert!(proof.independent_volume_discrepancy_mm3() <= 1e-5);
        assert!(proof.fixed_work_used() < 16_000_000 && proof.cavity_cells_used() > 0);
        eprintln!(
            "complete source match left{left} pose{posed}: {:?} C0{} chartD{} massgap{} work{} cells{}",
            now.elapsed(),
            proof.maximum_world_distance_mm(),
            proof.maximum_world_chart_derivative_mm(),
            proof.independent_volume_discrepancy_mm3(),
            proof.fixed_work_used(),
            proof.cavity_cells_used()
        );
    }
}
#[test]
fn actual_geometry_and_occurrence_mismatches_refuse_without_scalar_or_header_shortcuts() {
    let source = owned_source(false, 0.);
    let (text, occurrence) = file(&source, false);
    let decoded =
        crate::step_read::read_periodic_phillips_material(&text, Default::default(), None).unwrap();
    let good = source
        .match_material(decoded.material(), occurrence, Default::default(), None)
        .unwrap();
    let used = good.fixed_work_used();
    let cells = good.cavity_cells_used();
    for n in [0, 1, 8192, 8193, 28674, used - 1] {
        assert!(matches!(
            source.match_material(
                decoded.material(),
                occurrence,
                PeriodicMaterialMatchOptions {
                    maximum_fixed_work: n,
                    ..Default::default()
                },
                None
            ),
            Err(PhillipsReadError::Limit)
        ));
    }
    assert!(
        source
            .match_material(
                decoded.material(),
                occurrence,
                PeriodicMaterialMatchOptions {
                    maximum_fixed_work: used,
                    ..Default::default()
                },
                None
            )
            .is_ok()
    );
    assert!(matches!(
        source.match_material(
            decoded.material(),
            occurrence,
            PeriodicMaterialMatchOptions {
                maximum_cavity_cells: cells - 1,
                ..Default::default()
            },
            None
        ),
        Err(PhillipsReadError::Limit)
    ));
    assert!(
        source
            .match_material(
                decoded.material(),
                occurrence,
                PeriodicMaterialMatchOptions {
                    maximum_cavity_cells: cells,
                    ..Default::default()
                },
                None
            )
            .is_ok()
    );
    for o in [
        PeriodicMaterialMatchOptions {
            maximum_fixed_work: 16_000_001,
            ..Default::default()
        },
        PeriodicMaterialMatchOptions {
            maximum_distance_mm: f64::NAN,
            ..Default::default()
        },
        PeriodicMaterialMatchOptions {
            maximum_derivative_mm: 1.,
            ..Default::default()
        },
        PeriodicMaterialMatchOptions {
            maximum_volume_discrepancy_mm3: 0.,
            ..Default::default()
        },
        PeriodicMaterialMatchOptions {
            maximum_cavity_cells: 0,
            ..Default::default()
        },
        PeriodicMaterialMatchOptions {
            maximum_cavity_cells: 262145,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            source.match_material(decoded.material(), occurrence, o, None),
            Err(PhillipsReadError::Invalid)
        ));
    }
    let mut moved = occurrence;
    moved[0][0] = 1.;
    assert!(matches!(
        source.match_material(decoded.material(), moved, Default::default(), None),
        Err(PhillipsReadError::Unsupported)
    ));
    let other = owned_source(false, 0.4);
    let (other_text, _) = file(&other, false);
    let changed =
        crate::step_read::read_periodic_phillips_material(&other_text, Default::default(), None)
            .unwrap();
    assert!(matches!(
        source.match_material(changed.material(), occurrence, Default::default(), None),
        Err(PhillipsReadError::Unsupported)
    ));
    let opposite = owned_source(true, 0.);
    let expected = source.world_volume_mm3();
    let mirror = opposite.world_volume_mm3();
    assert!(expected[0] <= mirror[1] && expected[1] >= mirror[0]);
    let (opposite_text, _) = file(&opposite, false);
    let opposite_read =
        crate::step_read::read_periodic_phillips_material(&opposite_text, Default::default(), None)
            .unwrap();
    assert!(matches!(
        source.match_material(
            opposite_read.material(),
            occurrence,
            Default::default(),
            None
        ),
        Err(PhillipsReadError::Unsupported)
    ));
    let mut bad = occurrence;
    bad[1] = [0.; 3];
    assert!(matches!(
        source.match_material(decoded.material(), bad, Default::default(), None),
        Err(PhillipsReadError::Invalid)
    ));
    bad = occurrence;
    bad[0][0] = f64::NAN;
    assert!(matches!(
        source.match_material(decoded.material(), bad, Default::default(), None),
        Err(PhillipsReadError::Invalid)
    ));
    bad = occurrence;
    bad[3][2] = -1.;
    assert!(matches!(
        source.match_material(decoded.material(), bad, Default::default(), None),
        Err(PhillipsReadError::Invalid)
    ));
    assert!(matches!(
        source.match_material(
            decoded.material(),
            occurrence,
            PeriodicMaterialMatchOptions {
                maximum_volume_discrepancy_mm3: 1e-30,
                ..Default::default()
            },
            None
        ),
        Err(PhillipsReadError::Unsupported)
    ));
    assert!(matches!(
        source.match_material(
            decoded.material(),
            occurrence,
            Default::default(),
            Some(&|| true)
        ),
        Err(PhillipsReadError::Cancelled)
    ));
}

#[test]
fn actual_authored_sources_exercise_identity_caps_changed_radii_and_inventory() {
    use crate::head_profile::{Head, Profile};
    let original = degree11::tests::fixture::actual_source();
    let source = owned_source(false, 0.);
    let occurrence = [[0.; 3], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
    source
        .match_material(&source, occurrence, Default::default(), None)
        .unwrap();
    for (radius, plain) in [(2.01, 0.5), (2., 0.)] {
        let head = Head::Crowned {
            radius_mm: radius,
            height_mm: 1.6,
            crown_radius_mm: 3.2,
            top_land_radius_mm: 1.235645633015642,
            neck_radius_mm: 0.1,
        };
        let profile = Profile::new_for_thread_union(head, 1., plain).unwrap();
        let changed = OwnedPeriodicMaterial::from_authored_parts(
            original.thread(),
            &profile,
            original.cavity().source(),
            crate::profiled_screw::SourceSpec {
                head: Some(head),
                shank_length_mm: plain,
                ..original.spec()
            },
            Default::default(),
            None,
        )
        .unwrap();
        assert!(matches!(
            source.match_material(&changed, occurrence, Default::default(), None),
            Err(PhillipsReadError::Unsupported)
        ));
    }
}
