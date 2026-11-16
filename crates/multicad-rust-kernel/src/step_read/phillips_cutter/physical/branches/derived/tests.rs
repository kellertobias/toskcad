use super::*;
fn prerequisites(b: &Boundary) -> (FiniteLifts<'_>, CubicSupportBounds<'_>) {
    let endpoints = b
        .endpoint_coupling_bounds(
            EndpointCouplingOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 1_000_000,
            },
            None,
        )
        .unwrap();
    let analytic = b
        .analytic_support_bounds(
            &endpoints,
            AnalyticSupportOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 10_000_000,
            },
            None,
        )
        .unwrap();
    let lifts = analytic
        .finite_edge_lifts(endpoints, LiftOptions::default(), None)
        .unwrap();
    let cubics = b
        .cubic_support_bounds(
            CubicSupportOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 10_000_000,
            },
            None,
        )
        .unwrap();
    (lifts, cubics)
}
#[test]
fn actual_body_derived_receipt_retains_complete_roles_and_measured_transfers() {
    for scale in [1., 1.137] {
        for posed in [false, true] {
            let text = super::super::representative::scaled_text(scale, posed);
            let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
            let (lifts, cubics) = prerequisites(body.cavity());
            let sewn = lifts
                .sewn_edge_sources(&body, &cubics, SewnOptions::default(), None)
                .unwrap();
            let receipt = body
                .derived_cavity_representative(&sewn, DerivedCavityOptions::default(), None)
                .unwrap();
            assert!(std::ptr::eq(receipt.body(), &body));
            assert_eq!(
                (
                    receipt.vertices().len(),
                    receipt.edges().len(),
                    receipt.faces().len()
                ),
                (49, 76, 29)
            );
            let topology = receipt.derived_cavity_source().boundary_topology().unwrap();
            for (i, v) in receipt.vertices().iter().enumerate() {
                assert_eq!(v.source_role, topology.vertices()[i].source_identity());
                let raw = body
                    .cavity()
                    .vertices()
                    .iter()
                    .find(|r| r.entity == v.raw_vertex_entity)
                    .unwrap();
                assert_eq!(raw.point_entity, v.raw_point_entity);
            }
            for (i, e) in receipt.edges().iter().enumerate() {
                assert_eq!(e.source_role, topology.edges()[i].source_identity());
                let raw = body
                    .cavity()
                    .edges()
                    .iter()
                    .find(|r| r.entity == e.raw_edge_entity)
                    .unwrap();
                assert_eq!(raw.geometry, e.raw_curve_entity);
                assert_eq!(raw.same_sense, e.raw_curve_forward);
                let source = topology.edges()[i]
                    .vertices()
                    .map(|v| receipt.vertices()[v].raw_vertex_entity);
                let mut actual = raw.vertices.map(|v| body.cavity().vertices()[v].entity);
                if e.raw_edge_reversed {
                    actual.reverse()
                };
                assert_eq!(source, actual);
            }
            for (i, f) in receipt.faces().iter().enumerate() {
                assert_eq!(f.source_role, topology.faces()[i].source_identity());
                let raw = body
                    .cavity()
                    .faces()
                    .iter()
                    .find(|r| r.material_entity == f.raw_material_entity)
                    .unwrap();
                assert_eq!(raw.surface_entity, f.raw_surface_entity);
                assert_eq!(raw.material_forward, f.raw_material_forward);
                assert_eq!(raw.cavity_forward, f.raw_cavity_forward);
            }
            let p = receipt.transfer_bounds();
            assert_eq!(
                (p.vertex_count(), p.edge_count(), p.face_count()),
                (49, 76, 29)
            );
            assert!(p.certified_span_count() > 1000);
            assert!(
                p.maximum_vertex_error_mm() > 0. && p.maximum_vertex_error_mm() < DEFAULT_TRANSFER
            );
            assert!(
                p.maximum_analytic_support_error_mm() > 0.
                    && p.maximum_analytic_support_error_mm() < DEFAULT_TRANSFER
            );
            assert!(
                p.maximum_polynomial_error_mm() > 0.
                    && p.maximum_polynomial_support_error_mm() < DEFAULT_TRANSFER
            );
            assert!(
                receipt
                    .derived_frame()
                    .into_iter()
                    .flatten()
                    .all(f64::is_finite)
            );
            assert!(
                receipt.work_used() > sewn.work_used() + SOURCE_WORK && receipt.work_used() < WORK
            );
            assert!(
                receipt.retained_budget_bytes()
                    >= body.retained_budget_bytes() + std::mem::size_of_val(&receipt)
            );
            assert_eq!(receipt.body().literal_input(), text);
        }
    }
}
#[test]
fn limits_precision_and_cancellation_never_publish_a_partial_derived_source() {
    let text = super::super::representative::scaled_text(1., false);
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let (lifts, cubics) = prerequisites(body.cavity());
    let sewn = lifts
        .sewn_edge_sources(&body, &cubics, SewnOptions::default(), None)
        .unwrap();
    for options in [
        DerivedCavityOptions {
            maximum_transfer_mm: f64::NAN,
            ..Default::default()
        },
        DerivedCavityOptions {
            maximum_transfer_mm: 0.,
            ..Default::default()
        },
        DerivedCavityOptions {
            maximum_transfer_mm: 1e-7,
            ..Default::default()
        },
        DerivedCavityOptions {
            maximum_work: WORK + 1,
            ..Default::default()
        },
        DerivedCavityOptions {
            maximum_bytes: 0,
            ..Default::default()
        },
        DerivedCavityOptions {
            maximum_bytes: BYTES + 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            body.derived_cavity_representative(&sewn, options, None),
            Err(PhillipsReadError::Invalid)
        ));
    }
    for options in [
        DerivedCavityOptions {
            maximum_work: sewn.work_used(),
            ..Default::default()
        },
        DerivedCavityOptions {
            maximum_bytes: sewn.retained_capacity_bytes(),
            ..Default::default()
        },
    ] {
        assert!(matches!(
            body.derived_cavity_representative(&sewn, options, None),
            Err(PhillipsReadError::Limit)
        ));
    }
    assert!(matches!(
        body.derived_cavity_representative(
            &sewn,
            DerivedCavityOptions {
                maximum_transfer_mm: 1e-15,
                ..Default::default()
            },
            None
        ),
        Err(PhillipsReadError::Unsupported)
    ));
    let calls = std::cell::Cell::new(0);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    let success = body
        .derived_cavity_representative(&sewn, Default::default(), Some(&count))
        .unwrap();
    for prefix in [
        0,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        calls.get() / 4,
        calls.get() / 2,
        calls.get() - 1,
    ] {
        let n = std::cell::Cell::new(0);
        let cancel = || {
            let i = n.get();
            n.set(i + 1);
            i >= prefix
        };
        assert!(
            matches!(
                body.derived_cavity_representative(&sewn, Default::default(), Some(&cancel)),
                Err(PhillipsReadError::Cancelled)
            ),
            "prefix {prefix}/{}",
            calls.get()
        );
    }
    let retry = body
        .derived_cavity_representative(&sewn, Default::default(), None)
        .unwrap();
    assert_eq!(retry.transfer_bounds(), success.transfer_bounds());
}

#[test]
fn complete_original_body_charge_is_included_before_source_allocation() {
    let text = super::super::representative::scaled_text(1., false);
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let (lifts, cubics) = prerequisites(body.cavity());
    let sewn = lifts
        .sewn_edge_sources(&body, &cubics, SewnOptions::default(), None)
        .unwrap();
    let old_omission = sewn.retained_capacity_bytes()
        + std::mem::size_of::<DerivedCavityRepresentative<'_, '_>>()
        + SCRATCH_BYTES;
    let complete = old_omission + body.retained_budget_bytes();
    assert!(complete < BYTES);
    // This budget used to pass while omitting most of the live original arena.
    assert!(matches!(
        body.derived_cavity_representative(
            &sewn,
            DerivedCavityOptions {
                maximum_bytes: complete - 1,
                ..Default::default()
            },
            None
        ),
        Err(PhillipsReadError::Limit)
    ));
    assert!(complete - 1 > old_omission);
    let receipt = body
        .derived_cavity_representative(
            &sewn,
            DerivedCavityOptions {
                maximum_bytes: complete,
                ..Default::default()
            },
            None,
        )
        .unwrap();
    assert_eq!(
        receipt.retained_budget_bytes(),
        body.retained_budget_bytes() + std::mem::size_of_val(&receipt)
    );
    const { assert!(RAW_SPLINE_SCRATCH < 4 * 1024 * 1024) };
}
#[test]
fn actual_shallow_v_planes_earn_complete_independent_derived_cavity() {
    let text = super::super::tests::text(1, false);
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let (lifts, cubics) = prerequisites(body.cavity());
    let sewn = lifts
        .sewn_edge_sources(&body, &cubics, SewnOptions::default(), None)
        .unwrap();
    let receipt = body
        .derived_cavity_representative(&sewn, Default::default(), None)
        .unwrap();
    assert!(matches!(
        receipt.derived_cavity_source().spec().inner,
        crate::phillips_drive::Inner::ShallowV { .. }
    ));
    let transfer = receipt.transfer_bounds();
    assert_eq!(
        (
            transfer.vertex_count(),
            transfer.edge_count(),
            transfer.face_count()
        ),
        (49, 76, 29)
    );
    assert!(
        transfer.maximum_vertex_error_mm() < 1e-7
            && transfer.maximum_analytic_support_error_mm() < 1e-7
            && transfer.maximum_polynomial_error_mm() < 1e-7
            && transfer.maximum_polynomial_support_error_mm() < 1e-7
    );
}

// Refine one actual imported carrier to the largest admitted cubic schema.
// All endpoints keep original POINT identity; new interior coefficients are
// explicit decimal POINT leaves, not a private fabricated Boundary instance.
fn maximum_span_text(text: &str) -> String {
    let body = read_phillips_body(text, 1e-7, 1_000_000).unwrap();
    let (spline, ids, spans) = body
        .cavity()
        .edges()
        .iter()
        .find_map(|e| match &e.curve {
            Curve::CubicIntersection {
                spline_entity,
                pole_entities,
                spans,
                ..
            } => Some((*spline_entity, pole_entities, spans)),
            _ => None,
        })
        .unwrap();
    let entities = part21::parse(text).unwrap();
    let mut next = entities.keys().last().unwrap() + 1;
    let mut point_records = String::new();
    let mut poles = vec![ids[0]];
    let mut knots = vec![0.];
    let point = |p: [Point; 4], t: f64| {
        let mut p = p;
        for level in (1..4).rev() {
            for i in 0..level {
                p[i] = std::array::from_fn(|k| (1. - t) * p[i][k] + t * p[i + 1][k]);
            }
        }
        p[0]
    };
    let derivative = |p: [Point; 4], t: f64| {
        let q: [Point; 3] =
            std::array::from_fn(|i| std::array::from_fn(|k| 3. * (p[i + 1][k] - p[i][k])));
        let a: Point = std::array::from_fn(|k| (1. - t) * q[0][k] + t * q[1][k]);
        let b: Point = std::array::from_fn(|k| (1. - t) * q[1][k] + t * q[2][k]);
        std::array::from_fn::<_, 3, _>(|k| (1. - t) * a[k] + t * b[k])
    };
    for (index, span) in spans.iter().enumerate() {
        let count = 4096 / spans.len() + usize::from(index < 4096 % spans.len());
        for part in 0..count {
            let a = part as f64 / count as f64;
            let b = (part + 1) as f64 / count as f64;
            let delta = (b - a) / 3.;
            let pa = point(span.poles_mm, a);
            let pb = point(span.poles_mm, b);
            let da = derivative(span.poles_mm, a);
            let db = derivative(span.poles_mm, b);
            for (which, p) in [
                std::array::from_fn(|k| pa[k] + delta * da[k]),
                std::array::from_fn(|k| pb[k] - delta * db[k]),
                pb,
            ]
            .into_iter()
            .enumerate()
            {
                if which == 2 && part + 1 == count {
                    poles.push(ids[3 * index + 3]);
                } else {
                    let id = next;
                    next += 1;
                    poles.push(id);
                    point_records.push_str(&format!(
                        "#{id}=CARTESIAN_POINT('',({:.17E},{:.17E},{:.17E}));\n",
                        p[0], p[1], p[2]
                    ));
                }
            }
            knots.push(if part + 1 == count {
                span.parameters[1]
            } else {
                span.parameters[0] + b * (span.parameters[1] - span.parameters[0])
            });
        }
    }
    assert_eq!(knots.len(), 4097);
    assert_eq!(poles.len(), 12289);
    let ids = poles
        .iter()
        .map(|id| format!("#{id}"))
        .collect::<Vec<_>>()
        .join(",");
    let multiplicities = (0..4097)
        .map(|i| if i == 0 || i == 4096 { "4" } else { "3" })
        .collect::<Vec<_>>()
        .join(",");
    let knots = knots
        .iter()
        .map(|k| format!("{k:.17E}"))
        .collect::<Vec<_>>()
        .join(",");
    let old = text
        .lines()
        .find(|l| l.starts_with(&format!("#{spline}=")))
        .unwrap();
    let new = format!(
        "#{spline}=B_SPLINE_CURVE_WITH_KNOTS('',3,({ids}),.UNSPECIFIED.,.F.,.F.,({multiplicities}),({knots}),.UNSPECIFIED.);"
    );
    let mut changed = text.replacen(old, &new, 1);
    let end = changed.rfind("ENDSEC;").unwrap();
    changed.insert_str(end, &point_records);
    changed
}
#[test]
fn maximum_actual_cubic_schema_fits_the_declared_workspace_and_global_work() {
    let text = maximum_span_text(&super::super::representative::scaled_text(1., false));
    let body = read_phillips_body(&text, 1e-7, WORK).unwrap();
    let (lifts, cubics) = prerequisites(body.cavity());
    let sewn = lifts
        .sewn_edge_sources(&body, &cubics, SewnOptions::default(), None)
        .unwrap();
    let receipt = body
        .derived_cavity_representative(&sewn, Default::default(), None)
        .unwrap();
    assert_eq!(
        (
            receipt.vertices().len(),
            receipt.edges().len(),
            receipt.faces().len()
        ),
        (49, 76, 29)
    );
    let b = body.cavity();
    let raw = b
        .edges()
        .iter()
        .find(|e| matches!(&e.curve,Curve::CubicIntersection{spans,..}if spans.len()==4096))
        .unwrap();
    let mut work = Work {
        remaining: WORK,
        interrupt: None,
        cancelled: std::cell::Cell::new(false),
    };
    let spans = geometry::spline(
        &Shell::new(body.source_entities(), 1e-7, true),
        raw.geometry,
        raw.owners.map(|i| b.faces()[i].surface_entity),
        if raw.same_sense {
            raw.vertices
        } else {
            [raw.vertices[1], raw.vertices[0]]
        }
        .map(|i| b.vertices()[i].point_entity),
        &mut work,
    )
    .unwrap();
    assert_eq!(spans.len(), 4096);
    assert!(spans.capacity() <= 8192);
    let closed = receipt
        .closed_boundary_certificate(
            crate::phillips_drive::ClosedBoundaryOptions::default(),
            None,
        )
        .unwrap();
    assert!(closed.work_used() < WORK);
    assert_eq!(closed.local_source_certificate().patch_domains().len(), 88);
    assert!(receipt.work_used() < WORK);
    assert!(receipt.retained_budget_bytes() + SCRATCH_BYTES < BYTES);
    eprintln!(
        "max-schema derived source work={} retained={} spanCapacity={}",
        receipt.work_used(),
        receipt.retained_budget_bytes(),
        spans.capacity()
    );
}

#[test]
fn actual_torus_radius_is_mandatory_in_final_whole_support_comparison() {
    let text = super::super::representative::scaled_text(1., false);
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let entities = body.source_entities();
    let torus = body
        .cavity()
        .faces()
        .iter()
        .find(|f| matches!(f.support, Surface::Torus { .. }))
        .unwrap();
    let words = part21::numeric_spellings(&entities[&torus.surface_entity][0]).unwrap();
    let record = text
        .lines()
        .find(|l| l.starts_with(&format!("#{}=", torus.surface_entity)))
        .unwrap();
    let changed_record = record.replacen(
        words[0],
        &format!("{:.17E}", words[0].parse::<f64>().unwrap() + 0.001),
        1,
    );
    let changed = text.replacen(record, &changed_record, 1);
    assert_ne!(changed, text);
    let altered = read_phillips_body(&changed, 1e-7, 1_000_000).unwrap();
    let mut work = Work {
        remaining: WORK,
        interrupt: None,
        cancelled: std::cell::Cell::new(false),
    };
    // Proposal ignores the changed major radius, so this reaches the final proof.
    let candidate = fit::candidate(altered.cavity(), &mut work).unwrap();
    let topology = candidate.source.boundary_topology().unwrap();
    let boundary = altered.cavity();
    let views = boundary
        .faces()
        .iter()
        .map(|f| SourceFaceView {
            surface: f.surface_entity,
            forward: f.cavity_forward,
            uses: f
                .uses
                .iter()
                .map(|&(e, s)| (boundary.edges()[e].entity, s))
                .collect(),
        })
        .collect::<Vec<_>>();
    let mut shell = Shell::new(altered.source_entities(), 1e-7, true);
    let graph = graph::read_views(
        &mut shell,
        &views,
        &topology,
        &candidate.pose,
        1e-7,
        &mut work,
    )
    .unwrap();
    assert!(geometry::verify(&shell, &graph, &topology, &candidate.pose, 1e-7, &mut work).is_err());
    assert_eq!(body.literal_input(), text);
}

#[test]
fn result_borrows_original_body_without_retaining_prerequisite_buffers() {
    let text = super::super::representative::scaled_text(1., false);
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let receipt = {
        let (lifts, cubics) = prerequisites(body.cavity());
        let sewn = lifts
            .sewn_edge_sources(&body, &cubics, SewnOptions::default(), None)
            .unwrap();
        let foreign = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
        assert!(matches!(
            foreign.derived_cavity_representative(&sewn, Default::default(), None),
            Err(PhillipsReadError::Invalid)
        ));
        body.derived_cavity_representative(&sewn, Default::default(), None)
            .unwrap()
    };
    assert!(std::ptr::eq(receipt.body(), &body));
    assert_eq!(receipt.edges().len(), 76);
}

#[test]
fn reordered_actual_cavity_faces_preserve_geometric_roles_and_raw_entities() {
    let text = super::super::representative::scaled_text(1., false);
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let inner = body
        .cavity()
        .faces()
        .iter()
        .filter_map(|f| {
            f.material_entity
                .map(|id| (id, matches!(f.support, Surface::Torus { .. })))
        })
        .collect::<BTreeMap<_, _>>();
    let (&shell, record) = body
        .source_entities()
        .iter()
        .find_map(|(id, v)| (v.len() == 1 && v[0].name == "CLOSED_SHELL").then_some((id, &v[0])))
        .unwrap();
    let Value::List(original) = &record.args[1] else {
        panic!()
    };
    let original = original
        .iter()
        .map(|v| reference(v).unwrap())
        .collect::<Vec<_>>();
    let mut cavity = inner.keys().copied().collect::<Vec<_>>();
    cavity.sort_by_key(|id| (!inner[id], *id));
    let reordered = original
        .iter()
        .copied()
        .filter(|id| !inner.contains_key(id))
        .chain(cavity)
        .map(|id| format!("#{id}"))
        .collect::<Vec<_>>()
        .join(",");
    let old = text
        .lines()
        .find(|l| l.starts_with(&format!("#{shell}=")))
        .unwrap();
    let changed = text.replacen(old, &format!("#{shell}=CLOSED_SHELL('',({reordered}));"), 1);
    assert_ne!(text, changed);
    let other = read_phillips_body(&changed, 1e-7, 1_000_000).unwrap();
    assert!(matches!(
        other.cavity().faces()[1].support,
        Surface::Torus { .. }
    ));
    let (lifts, cubics) = prerequisites(body.cavity());
    let sewn = lifts
        .sewn_edge_sources(&body, &cubics, SewnOptions::default(), None)
        .unwrap();
    let before = body
        .derived_cavity_representative(&sewn, Default::default(), None)
        .unwrap();
    let (lifts, cubics) = prerequisites(other.cavity());
    let sewn = lifts
        .sewn_edge_sources(&other, &cubics, SewnOptions::default(), None)
        .unwrap();
    let after = other
        .derived_cavity_representative(&sewn, Default::default(), None)
        .unwrap();
    for (a, b) in before.vertices().iter().zip(after.vertices()) {
        assert_eq!(a.source_role, b.source_role);
        assert_eq!(a.raw_vertex_entity, b.raw_vertex_entity);
        assert_eq!(a.raw_point_entity, b.raw_point_entity);
    }
    for (a, b) in before.edges().iter().zip(after.edges()) {
        assert_eq!(a.source_role, b.source_role);
        assert_eq!(a.raw_edge_entity, b.raw_edge_entity);
        assert_eq!(a.raw_curve_entity, b.raw_curve_entity);
        assert_eq!(a.raw_edge_reversed, b.raw_edge_reversed);
        assert_eq!(a.raw_curve_forward, b.raw_curve_forward);
    }
    for (a, b) in before.faces().iter().zip(after.faces()) {
        assert_eq!(a.source_role, b.source_role);
        assert_eq!(a.raw_material_entity, b.raw_material_entity);
        assert_eq!(a.raw_surface_entity, b.raw_surface_entity);
        assert_eq!(a.raw_material_forward, b.raw_material_forward);
        assert_eq!(a.raw_cavity_forward, b.raw_cavity_forward);
    }
    assert_eq!(before.transfer_bounds(), after.transfer_bounds());
}

#[test]
fn actual_derived_cavity_closed_source_and_raw_affine_image_are_bound_to_original_body() {
    use crate::phillips_drive::ClosedBoundaryOptions;
    use num_rational::BigRational as Q;
    for posed in [false, true] {
        let text = super::super::representative::scaled_text(1.137, posed);
        let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
        let (lifts, cubics) = prerequisites(body.cavity());
        let sewn = lifts
            .sewn_edge_sources(&body, &cubics, SewnOptions::default(), None)
            .unwrap();
        let representative = body
            .derived_cavity_representative(&sewn, Default::default(), None)
            .unwrap();
        let closed = representative
            .closed_boundary_certificate(ClosedBoundaryOptions::default(), None)
            .unwrap();
        assert!(std::ptr::eq(closed.representative(), &representative));
        assert!(std::ptr::eq(
            closed.local_source_certificate().source(),
            representative.derived_cavity_source()
        ));
        assert_eq!(closed.local_source_certificate().patch_domains().len(), 88);
        assert_eq!(
            closed
                .local_source_certificate()
                .chart_topology()
                .euler_characteristic,
            2
        );
        let columns = closed.actual_affine_columns();
        assert_eq!(columns, representative.derived_affine_columns());
        assert_eq!(
            representative.derived_frame(),
            [columns[0], columns[1], columns[2]]
        );
        let [_, u, v, n] = columns.map(|p| p.map(|x| Q::from_float(x).unwrap()));
        let cross: [Q; 3] = std::array::from_fn(|i| {
            &v[(i + 1) % 3] * &n[(i + 2) % 3] - &v[(i + 2) % 3] * &n[(i + 1) % 3]
        });
        let determinant = (0..3).fold(Q::from_integer(0.into()), |s, i| s + &u[i] * &cross[i]);
        let bounds = closed
            .determinant_bounds()
            .map(|x| Q::from_float(x).unwrap());
        assert!(bounds[0] > Q::from_integer(0.into()));
        assert!(bounds[0] <= determinant && determinant <= bounds[1]);
        if posed {
            assert_ne!(determinant, Q::from_integer(1.into()));
        }
        assert!(closed.work_used() > representative.work_used() && closed.work_used() < WORK);
        assert_eq!(
            closed.retained_budget_bytes(),
            representative.retained_budget_bytes() + std::mem::size_of_val(&closed)
        );
        assert_eq!(closed.representative().body().literal_input(), text);
    }
}
#[test]
fn closed_derived_source_resource_and_cancellation_failures_leave_original_receipt_reusable() {
    use crate::phillips_drive::{ClosedBoundaryCertificate, ClosedBoundaryOptions};
    let text = super::super::representative::scaled_text(1., false);
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let (lifts, cubics) = prerequisites(body.cavity());
    let sewn = lifts
        .sewn_edge_sources(&body, &cubics, SewnOptions::default(), None)
        .unwrap();
    let representative = body
        .derived_cavity_representative(&sewn, Default::default(), None)
        .unwrap();
    for options in [
        ClosedBoundaryOptions {
            maximum_work: 0,
            ..Default::default()
        },
        ClosedBoundaryOptions {
            maximum_work: WORK + 1,
            ..Default::default()
        },
        ClosedBoundaryOptions {
            maximum_bytes: 0,
            ..Default::default()
        },
        ClosedBoundaryOptions {
            maximum_bytes: BYTES + 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            representative.closed_boundary_certificate(options, None),
            Err(PhillipsReadError::Invalid)
        ));
    }
    let peak = representative.retained_budget_bytes()
        + std::mem::size_of::<ClosedDerivedCavityCertificate<'_, '_, '_>>()
        + ClosedBoundaryCertificate::scratch_budget_bytes();
    assert!(matches!(
        representative.closed_boundary_certificate(
            ClosedBoundaryOptions {
                maximum_bytes: peak - 1,
                ..Default::default()
            },
            None
        ),
        Err(PhillipsReadError::Limit)
    ));
    let good = representative
        .closed_boundary_certificate(
            ClosedBoundaryOptions {
                maximum_bytes: peak,
                ..Default::default()
            },
            None,
        )
        .unwrap();
    let work = good.work_used();
    for maximum_work in [
        representative.work_used(),
        representative.work_used() + 257,
        work - 1,
    ] {
        assert!(matches!(
            representative.closed_boundary_certificate(
                ClosedBoundaryOptions {
                    maximum_work,
                    ..Default::default()
                },
                None
            ),
            Err(PhillipsReadError::Limit)
        ));
    }
    let calls = std::cell::Cell::new(0);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    representative
        .closed_boundary_certificate(Default::default(), Some(&count))
        .unwrap();
    for prefix in [0, 1, 2, 3, 4, 5, calls.get() / 2, calls.get() - 1] {
        let n = std::cell::Cell::new(0);
        let cancel = || {
            let old = n.get();
            n.set(old + 1);
            old >= prefix
        };
        assert!(matches!(
            representative.closed_boundary_certificate(Default::default(), Some(&cancel)),
            Err(PhillipsReadError::Cancelled)
        ));
    }
    let retry = representative
        .closed_boundary_certificate(Default::default(), None)
        .unwrap();
    assert_eq!(retry.determinant_bounds(), good.determinant_bounds());
    assert_eq!(
        retry.local_source_certificate().reserves(),
        good.local_source_certificate().reserves()
    );
}

#[test]
fn actual_unfitted_shallow_v_plane_remains_mandatory_in_whole_source_comparison() {
    let text = super::super::tests::text(1, false);
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let mut work = Work {
        remaining: WORK,
        interrupt: None,
        cancelled: std::cell::Cell::new(false),
    };
    let baseline = fit::candidate(body.cavity(), &mut work).unwrap();
    let dot = |a: Point, b: Point| a.into_iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
    let changed_face = body
        .cavity()
        .faces()
        .iter()
        .find(|face| match face.support {
            Surface::Plane { axis } => {
                let n = baseline.pose.axes.map(|a| dot(a, axis.normal));
                n[0] > 0. && n[1] > 0. && n[2] > 0.
            }
            _ => false,
        })
        .unwrap();
    // The proposal uses negative-local-Y upper planes. Shift a real lower
    // plane, preserving every original vertex, edge and other support. A fit
    // alone is consequently unchanged and cannot admit this altered file.
    let shell = Shell::new(body.source_entities(), 1e-7, true);
    let surface = changed_face.surface_entity;
    let axis_id = reference(&shell.args(surface, "PLANE").unwrap()[1]).unwrap();
    let axis = shell.args(axis_id, "AXIS2_PLACEMENT_3D").unwrap();
    let old_point = reference(&axis[1]).unwrap();
    let p = shell.point(old_point).unwrap();
    let normal = shell.direction(reference(&axis[2]).unwrap()).unwrap();
    let p: Point = std::array::from_fn(|i| p[i] + 0.001 * normal[i]);
    let new_point = body.source_entities().keys().max().unwrap() + 1;
    let new_axis = new_point + 1;
    let axis_record = body.source_entities()[&axis_id][0]
        .source
        .unwrap()
        .replacen(&format!("#{old_point},"), &format!("#{new_point},"), 1);
    let old_line = text
        .lines()
        .find(|line| line.starts_with(&format!("#{surface}=")))
        .unwrap();
    let new_line = old_line.replacen(&format!("#{axis_id}"), &format!("#{new_axis}"), 1);
    assert_ne!(old_line, new_line);
    let mut changed = text.replacen(old_line, &new_line, 1);
    let end = changed.rfind("ENDSEC;").unwrap();
    changed.insert_str(
        end,
        &format!(
            "#{new_point}=CARTESIAN_POINT('',({:.17E},{:.17E},{:.17E}));\n#{new_axis}={};\n",
            p[0],
            p[1],
            p[2],
            axis_record.trim_end_matches(';')
        ),
    );
    let altered = read_phillips_body(&changed, 1e-7, 1_000_000).unwrap();
    let candidate = fit::candidate(altered.cavity(), &mut work).unwrap();
    assert_eq!(candidate.source.spec(), baseline.source.spec());
    let topology = candidate.source.boundary_topology().unwrap();
    let b = altered.cavity();
    let views = b
        .faces()
        .iter()
        .map(|f| SourceFaceView {
            surface: f.surface_entity,
            forward: f.cavity_forward,
            uses: f
                .uses
                .iter()
                .map(|&(e, s)| (b.edges()[e].entity, s))
                .collect(),
        })
        .collect::<Vec<_>>();
    let mut shell = Shell::new(altered.source_entities(), 1e-7, true);
    let graph = graph::read_views(
        &mut shell,
        &views,
        &topology,
        &candidate.pose,
        1e-7,
        &mut work,
    )
    .unwrap();
    assert!(geometry::verify(&shell, &graph, &topology, &candidate.pose, 1e-7, &mut work).is_err());
    assert_eq!(body.literal_input(), text);
}

#[test]
fn actual_sewing_collar_numeric_composition_encloses_and_refuses_invalid_or_overflow() {
    let collar = sewn_edge_collar(0.25, [0.125, 0.5]).unwrap();
    assert!((0.75..0.750000000000001).contains(&collar));
    for (correction, owners) in [
        (-1., [0., 0.]),
        (f64::NAN, [0., 0.]),
        (0., [f64::INFINITY, 0.]),
        (0., [0., -1.]),
        (f64::MAX, [f64::MAX, 0.]),
    ] {
        assert!(sewn_edge_collar(correction, owners).is_err());
    }
}
