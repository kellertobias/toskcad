use super::*;
fn solid() -> Solid {
    Solid::rust_wood_insert(
        super::super::step_coupled_candidate_tests::actual_source(4_f64.next_up(), false),
        EntityId::new(),
        Default::default(),
        Default::default(),
        None,
    )
    .unwrap()
}
fn sphere_writer() -> step::Writer {
    let mut writer = step::Writer::new();
    writer
        .append_sphere(
            multicad_rust_kernel::sphere::Sphere::new(1.).unwrap(),
            [20., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
        )
        .unwrap();
    writer
}
fn changed_mapping(text: &str, reflected: bool, overflow: bool) -> String {
    let document = read::CompoundDocument::parse(text, Default::default(), None).unwrap();
    let operator = document.represented_members()[0]
        .mapping_link()
        .unwrap()
        .operator_entity();
    let line = text
        .lines()
        .find(|l| l.starts_with(&format!("#{operator}=")))
        .unwrap();
    let args = line
        .split("CARTESIAN_TRANSFORMATION_OPERATOR(")
        .nth(1)
        .unwrap()
        .split(')')
        .next()
        .unwrap();
    let id = args.split(',').nth(if reflected { 1 } else { 2 }).unwrap();
    let old = text
        .lines()
        .find(|l| l.starts_with(&format!("{}=", id.trim())))
        .unwrap();
    let new = if reflected {
        format!("{}=DIRECTION('',(0.,-1.,0.));", id.trim())
    } else {
        format!(
            "{}=CARTESIAN_POINT('',({},0.,0.));",
            id.trim(),
            if overflow { "1.E308" } else { "1." }
        )
    };
    text.replacen(old, &new, 1)
}
#[test]
fn real_topology_selected_arena_frontiers_and_plan_error_contracts() {
    let solid = solid();
    let (body, frame) = solid.wood_insert_body().unwrap();
    let policy = solid.policy;
    let reserve = SCRATCH + size_of::<WoodInsertStepVerification>();
    body.with_coupled_plan(policy, None, |topology| {
        let text = topology
            .write_mapped_candidate(frame.origin, frame.u, frame.v, WORK, BYTES, None)
            .unwrap();
        let slack = text.capacity() - text.len();
        let fixed = SCRATCH + size_of::<WoodInsertStepVerification>() + slack;
        let frontier = topology.retained_bytes() + fixed;
        for limit in [frontier - 1, frontier] {
            assert_eq!(
                selected_scope(topology, slack, reserve, limit, frame, policy)
                    .err()
                    .unwrap()
                    .kind,
                GeometryErrorKind::Limit
            );
        }
        let scope = selected_scope(topology, slack, reserve, frontier + 1, frame, policy).unwrap();
        assert_eq!(scope.available, 1);
        assert_eq!(scope.fixed, fixed);
        assert_eq!(
            selected_scope(topology, usize::MAX, reserve, BYTES, frame, policy)
                .err()
                .unwrap()
                .kind,
            GeometryErrorKind::Limit
        );
        assert_eq!(
            selected_scope(topology, slack, reserve, BYTES + 1, frame, policy)
                .err()
                .unwrap()
                .kind,
            GeometryErrorKind::Input
        );
        Ok(())
    })
    .unwrap();
    assert_eq!(
        body.with_coupled_plan(policy, Some(&|| true), |_| Ok(()))
            .unwrap_err()
            .kind,
        GeometryErrorKind::Cancelled
    );
    for (cause, kind) in [
        (insert::Error::Invalid, GeometryErrorKind::Input),
        (insert::Error::Limit, GeometryErrorKind::Limit),
        (insert::Error::Cancelled, GeometryErrorKind::Cancelled),
        (insert::Error::Unresolved, GeometryErrorKind::Unsupported),
    ] {
        let error = plan_error(cause);
        assert_eq!(error.kind, kind);
        assert!(error.message.contains(&format!("{cause:?}")));
    }
}
#[test]
fn owned_actual_candidate_refusals_leave_prior_sphere_intact() {
    let solid = solid();
    let (body, frame) = solid.wood_insert_body().unwrap();
    let policy = solid.policy;
    body.with_coupled_plan(policy, None, |topology| {
        let text = topology
            .write_mapped_candidate(frame.origin, frame.u, frame.v, WORK, BYTES, None)
            .unwrap();
        // Actual retained String capacity is live even when the STEP text is
        // short. Reject its selected-arena peak before parsing or append.
        let mut oversized = text.clone();
        oversized.reserve_exact(BYTES);
        assert!(oversized.capacity() - oversized.len() >= BYTES);
        let mut writer = sphere_writer();
        assert_eq!(
            body.append_verified_candidate(topology, oversized, frame, policy, &mut writer)
                .unwrap_err()
                .kind,
            GeometryErrorKind::Limit,
        );
        assert_eq!(writer.finish().unwrap(), sphere_writer().finish().unwrap());
        for (reflected, overflow, kind) in [
            (false, false, GeometryErrorKind::InvalidResult),
            (false, true, GeometryErrorKind::Unsupported),
            (true, false, GeometryErrorKind::Unsupported),
        ] {
            let changed = changed_mapping(&text, reflected, overflow);
            let mut writer = sphere_writer();
            let error = body
                .append_verified_candidate(topology, changed, frame, policy, &mut writer)
                .unwrap_err();
            assert_eq!(
                error.kind, kind,
                "reflection{reflected},overflow{overflow}: {error:?}"
            );
            assert_eq!(writer.finish().unwrap(), sphere_writer().finish().unwrap());
        }
        Ok(())
    })
    .unwrap();
}
#[test]
fn qualified_raw_local_pose_cancels_every_original_curve_and_final_publication() {
    let solid = solid();
    let (body, frame) = solid.wood_insert_body().unwrap();
    let policy = solid.policy;
    let reserve = SCRATCH + size_of::<WoodInsertStepVerification>();
    body.with_coupled_plan(policy, None, |topology| {
        let text = topology
            .write_mapped_candidate(frame.origin, frame.u, frame.v, WORK, BYTES, None)
            .unwrap();
        let document = read::CompoundDocument::parse(&text, Default::default(), None).unwrap();
        let member = &document.represented_members()[0];
        let packet = document
            .member_packet(member.item_entity(), Default::default(), None)
            .unwrap();
        let parsed = packet.parsed();
        let mapping = parsed
            .physical_mapping(
                member.mapping_link().unwrap().operator_entity(),
                Default::default(),
                None,
            )
            .unwrap();
        let scope = selected_scope(
            topology,
            text.capacity() - text.len(),
            reserve,
            BYTES,
            frame,
            policy,
        )
        .unwrap();
        let local = body
            .coupled_local_proof(topology, parsed, &mapping, scope, None)
            .unwrap();
        let polls = std::cell::Cell::new(0usize);
        let observe = || {
            polls.set(polls.get() + 1);
            false
        };
        let baseline =
            super::super::step_coupled_pose::transport(frame, local.clone(), Some(&observe))
                .unwrap();
        let curves = body
            .boundary()
            .source_boundary(source::Side::Outer)
            .source_curves()
            .count()
            + body
                .boundary()
                .source_boundary(source::Side::Inner)
                .source_curves()
                .count();
        assert_eq!(polls.get(), curves + 1);
        assert!(curves > 0);
        for stop in 0..polls.get() {
            let counter = std::cell::Cell::new(0usize);
            let cancel = || {
                let n = counter.get();
                counter.set(n + 1);
                n == stop
            };
            assert_eq!(
                super::super::step_coupled_pose::transport(frame, local.clone(), Some(&cancel))
                    .unwrap_err()
                    .kind,
                GeometryErrorKind::Cancelled
            );
        }
        assert_eq!(
            super::super::step_coupled_pose::transport(frame, local, None).unwrap(),
            baseline
        );
        Ok(())
    })
    .unwrap();
}
