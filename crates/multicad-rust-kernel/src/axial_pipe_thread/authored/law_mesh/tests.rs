use super::super::tests::bspt_spec;
use super::*;
fn print() -> Tessellation {
    Tessellation {
        deflection: 0.003,
        angle: 5f64.to_radians(),
        triangle_limit: 2_000_000,
    }
}
fn npt(k: f64, limit: Option<f64>) -> AuthoredThread {
    let p = 25.4 / 27.;
    let m = 3_f64.sqrt();
    let s = 1. / (m - k) + 1. / (m + k);
    let depth = p / s - 2. * 0.0645 * p;
    AuthoredThread::new(
        AuthoredSpec {
            major_radius_mm: (0.37360 * 25.4 - if k > 0. { 0.1615 * 25.4 / 16. } else { 0. }
                + depth)
                / 2.,
            pitch_mm: p,
            length_mm: 0.2639 * 25.4,
            radial_taper: k,
            section: AuthoredSection::NptFlat {
                depth_mm: depth,
                flat_width_mm: (p - depth * s) / 2.,
            },
            crest_radius_limit_mm: limit,
        },
        1e-7,
        1e-6,
    )
    .unwrap()
}
fn area(mesh: &DisplayMesh) -> f64 {
    mesh.indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let p = t.map(|i| {
                let i = 3 * i as usize;
                [mesh.vertices[i], mesh.vertices[i + 1], mesh.vertices[i + 2]]
            });
            0.5 * norm(cross(sub(p[1], p[0]), sub(p[2], p[0])))
        })
        .sum()
}
fn retained_volume(thread: &AuthoredThread) -> [f64; 2] {
    thread
        .raw_boundary_ownership(1e-7, None)
        .unwrap()
        .ownership
        .bind_surfaces(1e-7, None)
        .unwrap()
        .certify_axial(1e-7, 1e-6, None)
        .unwrap()
        .certify_caps(None)
        .unwrap()
        .certify_ordering(None)
        .unwrap()
        .certify_material(1e-7, None)
        .unwrap()
        .volume_with_interrupt(1e-6, None)
        .unwrap()
        .volume_mm3
}

#[test]
fn print_quality_law_solids_meet_chord_and_angle_bounds_and_enclose_the_exact_volume() {
    // 1/8-27 NPT male clipped by its 0.405" pipe OD, and a rounded BSPT R.
    let clipped = npt(1. / 32., Some(0.405 * 25.4 / 2.));
    let plan = clipped.stock_clip_plan(Default::default(), None).unwrap();
    let exact = plan
        .volume_bounds(1e-6, Default::default(), None)
        .unwrap()
        .source_volume_mm3();
    let rounded = AuthoredThread::new(bspt_spec(1. / 32.), 1e-7, 1e-6).unwrap();
    for (thread, exact, clips) in [
        (clipped, exact, true),
        (rounded, retained_volume(&rounded), false),
    ] {
        let law = thread
            .tessellate_law(print(), LawClosure::Solid, None)
            .unwrap();
        assert!(law.deviation_mm <= 0.003, "{}", law.deviation_mm);
        assert!(law.normal_angle_rad <= 5f64.to_radians());
        assert!(law.mesh.face_indices.len() <= 2_000_000);
        let proof = law.proof.unwrap();
        // Every facet lies within the deviation of the law surface.
        let reserve = law.deviation_mm * area(&law.mesh) + (exact[1] - exact[0]);
        let mid = (exact[0] + exact[1]) / 2.;
        assert!(
            (proof.signed_volume_mm3 - mid).abs() <= reserve,
            "{} {mid} {reserve}",
            proof.signed_volume_mm3
        );
        assert_eq!(
            law.mesh
                .face_indices
                .iter()
                .any(|&f| law.role(f) == Some(LawRole::Clip)),
            clips
        );
        for role in [LawRole::StartCap, LawRole::EndCap, LawRole::Band(1)] {
            let index = law.role_index(role);
            assert!(law.mesh.face_indices.contains(&index));
        }
        // Rims are counterclockwise and lie exactly on the cap planes.
        for (end, rim) in law.rims.iter().enumerate() {
            assert!(rim.len() >= law.rows);
            let z = [0., thread.spec().length_mm][end];
            let mut turned = 0.;
            for i in 0..rim.len() {
                let a = 3 * rim[i] as usize;
                let b = 3 * rim[(i + 1) % rim.len()] as usize;
                assert_eq!(law.mesh.vertices[a + 2], z);
                let ta = law.mesh.vertices[a + 1].atan2(law.mesh.vertices[a]);
                let tb = law.mesh.vertices[b + 1].atan2(law.mesh.vertices[b]);
                let step = (tb - ta).rem_euclid(std::f64::consts::TAU);
                assert!(step > 0. && step < 0.1);
                turned += step;
            }
            assert!((turned - std::f64::consts::TAU).abs() < 1e-9);
        }
    }
}

#[test]
fn bores_are_inward_lateral_surfaces_and_bad_requests_refuse() {
    for thread in [
        npt(-1. / 32., None),
        AuthoredThread::new(bspt_spec(-1. / 32.), 1e-7, 1e-6).unwrap(),
    ] {
        let law = thread
            .tessellate_law(print(), LawClosure::Bore, None)
            .unwrap();
        assert!(law.proof.is_none());
        assert!(law.deviation_mm <= 0.003);
        // Lateral faces only; caps are left to the receiving stock.
        assert!(
            law.mesh
                .face_indices
                .iter()
                .all(|&f| !matches!(law.role(f), Some(LawRole::StartCap | LawRole::EndCap)))
        );
        assert!(law.rims.iter().all(|r| r.len() >= law.rows));
        // Facet normals point towards the axis: material is outside.
        let t = &law.mesh.indices[..3];
        let p = t.iter().map(|&i| 3 * i as usize).collect::<Vec<_>>();
        let v = |i: usize| {
            [
                law.mesh.vertices[i],
                law.mesh.vertices[i + 1],
                law.mesh.vertices[i + 2],
            ]
        };
        let n = cross(sub(v(p[1]), v(p[0])), sub(v(p[2]), v(p[0])));
        let c = v(p[0]);
        assert!(n[0] * c[0] + n[1] * c[1] < 0.);
    }
    let clipped = npt(1. / 32., Some(0.405 * 25.4 / 2.));
    assert_eq!(
        clipped
            .tessellate_law(print(), LawClosure::Bore, None)
            .err(),
        Some(Error::Unsupported)
    );
    assert_eq!(
        clipped
            .tessellate_law(
                Tessellation {
                    triangle_limit: 1000,
                    ..print()
                },
                LawClosure::Solid,
                None
            )
            .err(),
        Some(Error::Limit)
    );
    assert_eq!(
        clipped
            .tessellate_law(print(), LawClosure::Solid, Some(&|| true))
            .err(),
        Some(Error::Cancelled)
    );
}

#[test]
fn law_bores_join_print_quality_cylindrical_stock_into_one_closed_shell() {
    use crate::axial_pipe_thread::authored::hardware::mesh::{MeshOptions, join_sleeve_law};
    for thread in [
        npt(-1. / 32., None),
        AuthoredThread::new(bspt_spec(-1. / 32.), 1e-7, 1e-6).unwrap(),
    ] {
        let law = thread
            .tessellate_law(print(), LawClosure::Bore, None)
            .unwrap();
        let start = law.role_index(LawRole::StartCap);
        let end = law.role_index(LawRole::EndCap);
        let radius = thread.spec().major_radius_mm + 2.;
        let length = thread.spec().length_mm;
        let deviation = law.deviation_mm;
        let lateral_area = area(&law.mesh);
        let (mesh, proof, distance) = join_sleeve_law(
            law,
            [0., length],
            [start, end, end + 1],
            radius,
            MeshOptions {
                tessellation: print(),
                maximum_work: 512_000_000,
                maximum_bytes: 64 * 1024 * 1024,
                additional_resident_bytes: 0,
            },
            None,
        )
        .unwrap();
        assert!(distance <= 0.003);
        assert_eq!(proof, crate::manufacturing::certify(&mesh).unwrap());
        let cavity = retained_volume(&thread);
        let n = crate::axial_pipe_thread::authored::hardware::mesh::sleeve_segments(radius, print())
            .unwrap() as f64;
        let polygon = 0.5 * n * radius * radius * (std::f64::consts::TAU / n).sin();
        let expected = polygon * length - (cavity[0] + cavity[1]) / 2.;
        assert!(
            (proof.signed_volume_mm3 - expected).abs()
                <= deviation * lateral_area + (cavity[1] - cavity[0]),
            "{} {expected}",
            proof.signed_volume_mm3
        );
    }
}

#[test]
fn sourced_npt_pipe_ends_at_both_sizes_close_at_print_quality() {
    // The 1/8-27 and 1/4-18 NBS plain-face pipe ends, 0.005 mm radial
    // allowance, clipped by their published pipe OD.
    for (r, p, l, d, w, od) in [
        (
            4.958071685363371,
            0.9407407407407407,
            6.70306,
            0.6930846207267441,
            0.07008747760177175,
            5.1435,
        ),
        (
            6.583342659989502,
            1.411111111111111,
            10.20572,
            1.0509158199790052,
            0.0986114510443532,
            6.858,
        ),
    ] {
        let thread = AuthoredThread::new(
            AuthoredSpec {
                major_radius_mm: r,
                pitch_mm: p,
                length_mm: l,
                radial_taper: 1. / 32.,
                section: AuthoredSection::NptFlat {
                    depth_mm: d,
                    flat_width_mm: w,
                },
                crest_radius_limit_mm: Some(od),
            },
            1e-7,
            1e-6,
        )
        .unwrap();
        let law = thread
            .tessellate_law(print(), LawClosure::Solid, None)
            .unwrap();
        assert!(law.deviation_mm <= 0.003 && law.normal_angle_rad <= 5f64.to_radians());
        let plan = thread.stock_clip_plan(Default::default(), None).unwrap();
        let exact = plan
            .volume_bounds(1e-6, Default::default(), None)
            .unwrap()
            .source_volume_mm3();
        let mid = (exact[0] + exact[1]) / 2.;
        let reserve = law.deviation_mm * area(&law.mesh) + (exact[1] - exact[0]);
        assert!((law.proof.unwrap().signed_volume_mm3 - mid).abs() <= reserve);
    }
}

#[test]
fn sourced_rc_receivers_join_stock_within_the_join_frontiers_at_print_quality() {
    use crate::axial_pipe_thread::authored::hardware::mesh::{MeshOptions, join_sleeve_law};
    // ISO 7 Rc 1/8 and Rc 1/4 receivers with 0.005 mm radial allowance.
    for (r, p, l, d) in [
        (
            4.883174107142857,
            0.9071428571428571,
            7.4,
            0.5808684050603288,
        ),
        (
            6.6043881578947365,
            1.3368421052631578,
            11.,
            0.8560165969310107,
        ),
    ] {
        let thread = AuthoredThread::new(
            AuthoredSpec {
                major_radius_mm: r,
                pitch_mm: p,
                length_mm: l,
                radial_taper: -1. / 32.,
                section: AuthoredSection::BsptRounded { depth_mm: d },
                crest_radius_limit_mm: None,
            },
            1e-7,
            1e-6,
        )
        .unwrap();
        let law = thread
            .tessellate_law(print(), LawClosure::Bore, None)
            .unwrap();
        assert!(law.deviation_mm <= 0.003 && law.normal_angle_rad <= 5f64.to_radians());
        let s = law.role_index(LawRole::StartCap);
        let e = law.role_index(LawRole::EndCap);
        let (mesh, proof, distance) = join_sleeve_law(
            law,
            [0., l],
            [s, e, e + 1],
            r + 2.,
            MeshOptions {
                tessellation: print(),
                maximum_work: 512_000_000,
                maximum_bytes: 64 * 1024 * 1024,
                additional_resident_bytes: 0,
            },
            None,
        )
        .unwrap();
        assert!(distance <= 0.003);
        assert!(mesh.face_indices.len() <= 2_000_000 && proof.signed_volume_mm3 > 0.);
    }
}
