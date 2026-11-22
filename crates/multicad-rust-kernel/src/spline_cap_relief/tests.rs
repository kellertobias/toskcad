use super::*;
use crate::{display::Tessellation, normal_sweep, spline_sweep::Cubic};
use multicad_precision::PortableMath as _;
fn square(x: f64, y: f64, s: f64) -> Vec<Edge> {
    let p = [[x, y], [x + s, y], [x + s, y + s], [x, y + s]];
    (0..4)
        .map(|i| {
            Edge::Analytic(normal_sweep::Edge::Line {
                a: p[i],
                b: p[(i + 1) % 4],
            })
        })
        .collect()
}
pub(crate) fn fixture(kind: Kind) -> Relief {
    let base = Sweep::new(square(0., 0., 20.), None, 5., 1e-7).unwrap();
    let points = [[3., 4.], [4., 5.], [5., 4.], [4., 3.]];
    let controls = [[3., 5.], [5., 5.], [5., 3.], [3., 3.]];
    let hole = (0..4)
        .map(|i| {
            let a = points[i];
            let b = points[(i + 1) % 4];
            let c = controls[i];
            Edge::Cubic(
                Cubic::new(
                    vec![
                        a,
                        std::array::from_fn(|j| a[j] + (c[j] - a[j]) * 2. / 3.),
                        std::array::from_fn(|j| b[j] + (c[j] - b[j]) * 2. / 3.),
                        b,
                    ],
                    vec![0., 1.],
                    false,
                )
                .unwrap(),
            )
        })
        .collect();
    let motif = Sweep::new_rings(
        square(2., 2., 5.),
        vec![
            hole,
            vec![Edge::Analytic(normal_sweep::Edge::Circle {
                center: [6., 6.],
                radius_mm: 0.25,
                clockwise: true,
            })],
        ],
        1.,
        1e-7,
    )
    .unwrap();
    Relief::new(base, vec![motif], kind, 1e-7, 1e-4).unwrap()
}
fn options() -> Tessellation {
    Tessellation {
        deflection: 0.002,
        angle: 0.08,
        triangle_limit: 100000,
    }
}
fn mesh_volume(mesh: &crate::display::DisplayMesh) -> f64 {
    mesh.indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let p = [t[0], t[1], t[2]].map(|i| {
                let i = i as usize;
                [
                    mesh.vertices[3 * i],
                    mesh.vertices[3 * i + 1],
                    mesh.vertices[3 * i + 2],
                ]
            });
            (p[0][0] * (p[1][1] * p[2][2] - p[1][2] * p[2][1])
                + p[0][1] * (p[1][2] * p[2][0] - p[1][0] * p[2][2])
                + p[0][2] * (p[1][0] * p[2][1] - p[1][1] * p[2][0]))
                / 6.
        })
        .sum()
}
#[test]
fn actual_counter_islands_share_all_curved_wires_and_keep_material_volume() {
    for kind in [Kind::Raised, Kind::Recessed] {
        let relief = fixture(kind);
        let graph = relief.boundary();
        assert_eq!(
            graph
                .faces
                .iter()
                .filter(|f| matches!(f.role, FaceRole::CounterCap { .. }))
                .count(),
            2
        );
        assert_eq!(
            graph
                .faces
                .iter()
                .find(|f| matches!(f.role, FaceRole::MotifCap { .. }))
                .unwrap()
                .loops
                .len(),
            3
        );
        let mut uses = vec![[0; 2]; graph.edges.len()];
        for face in &graph.faces {
            for wire in &face.loops {
                for (a, b) in wire.iter().zip(wire.iter().cycle().skip(1)) {
                    uses[a.edge][usize::from(a.forward)] += 1;
                    let ends = |e: &EdgeUse| {
                        let p = graph.edges[e.edge].vertices;
                        if e.forward { p } else { [p[1], p[0]] }
                    };
                    assert_eq!(ends(a)[1], ends(b)[0]);
                }
            }
        }
        assert!(uses.iter().all(|&u| u == [1, 1]));
        let sign = if kind == Kind::Raised { 1. } else { -1. };
        assert_eq!(
            relief.volume_mm3(),
            2000. + sign * relief.motifs()[0].volume_mm3()
        );
        let mesh = relief.mesh(options()).unwrap();
        assert!((mesh_volume(&mesh) - relief.volume_mm3()).abs() < 0.02);
        let mut edges = std::collections::BTreeMap::new();
        for t in mesh.indices.as_chunks::<3>().0.iter() {
            for i in 0..3 {
                let a = t[i];
                let b = t[(i + 1) % 3];
                let entry = edges.entry([a.min(b), a.max(b)]).or_insert([0; 2]);
                entry[usize::from(a < b)] += 1;
            }
        }
        assert!(edges.values().all(|&u| u == [1, 1]));
        assert_eq!(relief.base().height_mm(), 5.);
        assert_eq!(
            relief.motif_level_mm(0),
            Some(if kind == Kind::Raised { 6. } else { 4. })
        );
        assert_eq!(relief.motif_level_mm(1), None);
        assert_eq!(relief.minimum_mm(), 1e-4);
        assert_eq!(relief.tolerance_mm(), 1e-7);
        assert!(relief.centroid_mm().into_iter().all(f64::is_finite));
        assert_eq!(
            relief.bounds()[1][2],
            if kind == Kind::Raised { 6. } else { 5. }
        );
    }
}

#[test]
fn full_shared_step_graph_is_proved_in_both_directions_and_proper_pose() {
    for kind in [Kind::Raised, Kind::Recessed] {
        let relief = fixture(kind);
        for [origin, u, v] in [
            [[0.; 3], [1., 0., 0.], [0., 1., 0.]],
            [[80., 60., 20.], [0.8, 0.6, 0.], [-0.36, 0.48, 0.8]],
        ] {
            let text = crate::step::write_spline_cap_relief(&relief, origin, u, v).unwrap();
            assert_eq!(
                crate::step_read::verify_spline_cap_relief(&text, &relief, origin, u, v, 1e-7)
                    .unwrap(),
                relief.volume_mm3()
            );
            let (proof, imported) = crate::step_read::read_spline_cap_relief(&text, 1e-7).unwrap();
            assert!((proof.volume_mm3 - relief.volume_mm3()).abs() < 1e-8);
            assert_eq!(imported.relief.kind(), kind);
            assert_eq!(imported.relief.motifs().len(), 1);
            assert_eq!(imported.relief.motifs()[0].holes().len(), 2);
            let again = crate::step::write_spline_cap_relief(
                &imported.relief,
                imported.origin_mm,
                imported.u,
                imported.v,
            )
            .unwrap();
            assert!(crate::step_read::read_spline_cap_relief(&again, 1e-7).is_ok());
        }
    }
}

#[test]
fn literal_material_queries_caps_and_invalid_depths_keep_both_counters() {
    for kind in [Kind::Raised, Kind::Recessed] {
        let r = fixture(kind);
        let z = if kind == Kind::Raised { 5.5 } else { 4.5 };
        for (p, occupied) in [
            ([2.5, 2.5, z], kind == Kind::Raised),
            ([4., 4., z], kind == Kind::Recessed),
            ([6., 6., z], kind == Kind::Recessed),
            ([10., 10., 2.], true),
            ([-1., 0., 0.], false),
            ([2., 2., -1.], false),
            ([2., 2., 8.], false),
        ] {
            assert_eq!(r.contains(p, 1e-7).unwrap(), occupied);
        }
        assert!(r.contains([f64::NAN, 0., 0.], 1e-7).is_err());
        assert!(r.contains([1., 1., 1.], 0.).is_err());
        assert!(r.contains([3., 4., z], 1e-7).is_err());
        let cap_area = r
            .boundary
            .faces
            .iter()
            .enumerate()
            .filter_map(|(i, _)| r.cap_properties(i).map(|(a, _)| a))
            .sum::<f64>();
        assert!((cap_area - 800.).abs() < 1e-10);
        assert_eq!(r.cap_properties(usize::MAX), None);
        assert_eq!(
            r.support([0., 0., 1.]),
            if kind == Kind::Raised { 6. } else { 5. }
        );
        assert_eq!(r.support([1., 0., 0.]), 20.);
        assert_eq!(r.boundary_error_mm(), 0.);
    }
    let source = fixture(Kind::Raised);
    let base = source.base().clone();
    let motif = source.motifs()[0].clone();
    for (tol, min) in [
        (0., 1e-4),
        (f64::NAN, 1e-4),
        (1e-7, f64::INFINITY),
        (1e-7, 1e-8),
    ] {
        assert_eq!(
            Relief::new(base.clone(), vec![motif.clone()], Kind::Raised, tol, min).err(),
            Some(Error::Invalid)
        );
    }
    assert_eq!(
        Relief::new(base.clone(), vec![], Kind::Raised, 1e-7, 1e-4).err(),
        Some(Error::Invalid)
    );
    assert_eq!(
        Relief::new(
            base.clone(),
            vec![motif.clone(); MAX_MOTIFS + 1],
            Kind::Raised,
            1e-7,
            1e-4
        )
        .err(),
        Some(Error::Limit)
    );
    let short = Sweep::new(square(2., 2., 5.), None, 5e-5, 1e-7).unwrap();
    assert_eq!(
        Relief::new(base.clone(), vec![short], Kind::Raised, 1e-7, 1e-4).err(),
        Some(Error::Invalid)
    );
    for h in [5., 6.] {
        let deep = Sweep::new(square(2., 2., 5.), None, h, 1e-7).unwrap();
        assert_eq!(
            Relief::new(base.clone(), vec![deep], Kind::Recessed, 1e-7, 1e-4).err(),
            Some(Error::Invalid)
        );
    }
    let outside = Sweep::new(square(30., 2., 2.), None, 1., 1e-7).unwrap();
    assert!(Relief::new(base.clone(), vec![outside], Kind::Raised, 1e-7, 1e-4).is_err());
    assert!(
        Relief::new(
            base.clone(),
            vec![motif.clone(), motif.clone()],
            Kind::Raised,
            1e-7,
            1e-4
        )
        .is_err()
    );
    let calls = std::cell::Cell::new(0);
    let observe = || {
        calls.set(calls.get() + 1);
        false
    };
    Relief::new_with_interrupt(
        base.clone(),
        vec![motif.clone()],
        Kind::Raised,
        1e-7,
        1e-4,
        Some(&observe),
    )
    .unwrap();
    for at in [
        0,
        1,
        calls.get() / 2,
        calls.get() - 1,
        calls.get() - motif.edges().count() - 1,
    ] {
        let n = std::cell::Cell::new(0);
        let cancel = || {
            let i = n.get();
            n.set(i + 1);
            i == at
        };
        assert_eq!(
            Relief::new_with_interrupt(
                base.clone(),
                vec![motif.clone()],
                Kind::Raised,
                1e-7,
                1e-4,
                Some(&cancel)
            )
            .err(),
            Some(Error::Cancelled)
        );
    }
    assert_eq!(
        Relief::new(base, vec![motif], Kind::Raised, 1e-7, 1e-4).unwrap(),
        source
    );
}

#[test]
fn detached_mesh_limits_fail_before_publication_and_retry_the_same_geometry() {
    let r = fixture(Kind::Raised);
    let original = r.clone();
    let mut request = options();
    request.triangle_limit = 1;
    assert_eq!(
        r.mesh(request).err(),
        Some(crate::display::MeshError::BudgetExceeded)
    );
    request = options();
    request.deflection = 0.;
    assert_eq!(
        r.mesh(request).err(),
        Some(crate::display::MeshError::InvalidOptions)
    );
    request = options();
    request.deflection = 1e-19;
    assert!(r.mesh(request).is_err());
    assert_eq!(r, original);
    r.mesh(options()).unwrap();
}

#[test]
fn finite_mass_interpolation_and_shared_join_precision_refuse_before_mesh_publication() {
    assert_eq!(finish_metrics(0., [0.; 3], [0.; 3]), Err(Error::Unresolved));
    assert_eq!(
        finish_metrics(1., [f64::INFINITY, 0., 0.], [0.; 3]),
        Err(Error::Unresolved)
    );
    assert_eq!(
        finish_metrics(2., [2., 4., 6.], [1.; 3]).unwrap(),
        [2., 3., 4.]
    );
    let line = Curve::Line {
        a_mm: [1., 0., 2.],
        b_mm: [0., 1., 2.],
    };
    assert_eq!(line.point(0.25), [0.75, 0.25, 2.]);
    assert_eq!(line.moments([0.; 2]), [0.5, 1. / 6., 1. / 6.]);
    let source = fixture(Kind::Raised);
    let huge = Sweep::new(square(0., 0., 20.), None, 1e8, 1e-7).unwrap();
    assert_eq!(
        Relief::new(huge, source.motifs().to_vec(), Kind::Raised, 1e-7, 1e-4).err(),
        Some(Error::Limit)
    );
    let mut edges = square(0., 0., 20.);
    let Edge::Analytic(normal_sweep::Edge::Line { a, .. }) = &mut edges[1] else {
        panic!()
    };
    a[0] += 1e-12;
    let base = Sweep::new(edges, None, 5., 1e-7).unwrap();
    let r = Relief::new(base, source.motifs().to_vec(), Kind::Raised, 1e-7, 1e-4).unwrap();
    assert!(r.boundary_error_mm() > 0.);
    let mut options = options();
    options.deflection = r.boundary_error_mm() / 2.;
    assert_eq!(
        r.mesh(options).err(),
        Some(crate::display::MeshError::UnresolvedBoundary)
    );
    for limit in [10, 100, 200] {
        let mut o = super::tests::options();
        o.triangle_limit = limit;
        assert_eq!(
            source.mesh(o).err(),
            Some(crate::display::MeshError::BudgetExceeded)
        );
    }
}

#[test]
fn simultaneous_xy_join_and_nonbinary_height_errors_have_one_outward_reserve() {
    let source = fixture(Kind::Raised);
    let mut edges = square(0., 0., 20.);
    let Edge::Analytic(normal_sweep::Edge::Line { a, .. }) = &mut edges[1] else {
        panic!()
    };
    a[0] += 1e-12;
    let base = Sweep::new(edges, None, 0.1, 1e-7).unwrap();
    let join = base.join_adjustment_mm();
    let m = &source.motifs()[0];
    let motif = Sweep::new_rings(m.outline().to_vec(), m.holes().to_vec(), 0.2, 1e-7).unwrap();
    let relief = Relief::new(base, vec![motif], Kind::Raised, 1e-7, 1e-4).unwrap();
    let (_, z_error) = sum_error(0.1, 0.2);
    assert!(join > 0. && z_error != 0.);
    assert!(relief.boundary_error_mm() > join.max(z_error.abs()));
    assert!(relief.boundary_error_mm() >= join + z_error.abs());
    assert_eq!(
        relief.contains([2.5, 2.5, 0.1 + 0.2], 1e-7),
        Err(Error::Unresolved)
    );
    assert_eq!(relief.contains([2.5, 2.5, 0.25], 1e-7), Ok(true));
    assert_eq!(relief.contains([2.5, 2.5, 0.35], 1e-7), Ok(false));
    let mut request = options();
    request.deflection = relief.boundary_error_mm() / 2.;
    assert_eq!(
        relief.mesh(request).err(),
        Some(crate::display::MeshError::UnresolvedBoundary)
    );
    relief.mesh(options()).unwrap();
    let base = Sweep::new(square(0., 0., 20.), None, 0.4, 1e-7).unwrap();
    let m = &source.motifs()[0];
    let motif = Sweep::new_rings(m.outline().to_vec(), m.holes().to_vec(), 0.1, 1e-7).unwrap();
    let cut = Relief::new(base, vec![motif], Kind::Recessed, 1e-7, 1e-4).unwrap();
    let level = cut.motif_level_mm(0).unwrap();
    assert_eq!(
        cut.contains([2.5, 2.5, level], 1e-7),
        Err(Error::Unresolved)
    );
    assert_eq!(cut.contains([2.5, 2.5, 0.25], 1e-7), Ok(true));
    assert_eq!(cut.contains([2.5, 2.5, 0.35], 1e-7), Ok(false));
}

#[test]
fn reversing_every_actual_edge_preserves_its_image_and_restores_literal_source_data() {
    let source = fixture(Kind::Raised);
    for edge in source
        .base()
        .edges()
        .chain(source.motifs().iter().flat_map(Sweep::edges))
    {
        let reversed = edge.reversed();
        assert_eq!(reversed.reversed(), *edge);
        for t in [0., 0.125, 0.5, 0.875, 1.] {
            let a = edge.point(1. - t);
            let b = reversed.point(t);
            assert!((a[0] - b[0]).portable_hypot(a[1] - b[1]) < 1e-12);
        }
    }
}

#[test]
fn remaining_chord_and_angle_budgets_round_outward_down_without_filling_precision_gaps() {
    assert_eq!(remaining_budget(1., 0.1), 0.9_f64.next_down());
    assert_eq!(remaining_budget(1., 0.25), 0.75);
    assert_eq!(remaining_budget(1., 0.), 1.);
    let source = fixture(Kind::Raised);
    let mut requested = options();
    requested.angle = 5e-7;
    assert_eq!(
        source.mesh(requested).err(),
        Some(crate::display::MeshError::InvalidOptions)
    );
    source.mesh(options()).unwrap();
}
