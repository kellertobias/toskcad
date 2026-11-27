use multicad_core::{
    commands::{Actor, Channel, Session},
    document::{Document, Entity, EntityData, EntityId, Mutation},
    drawings::*,
    project_file::{self, ProjectFile},
    references::ReferenceFrame,
};
fn actor() -> Actor {
    Actor::local(Channel::Ui)
}
fn fixture() -> (Session, EntityId, EntityId) {
    let mut doc = Document::new("Views").unwrap();
    let part = Entity::new(
        "Part",
        EntityData::Part {
            features: vec![],
            material: None,
            tags: Default::default(),
        },
    );
    let source = part.id;
    doc.insert(part).unwrap();
    let mut page = Page::new("A1");
    page.views.push(view(source));
    let book = DrawingBook { pages: vec![page] }.entity("Drawing").unwrap();
    let id = book.id;
    doc.insert(book).unwrap();
    (Session::new(doc).unwrap(), id, source)
}
fn view(source: EntityId) -> DrawingView {
    DrawingView {
        id: EntityId::new(),
        name: "Front".into(),
        source: ViewSource::Model {
            entities: vec![source],
        },
        frame: ReferenceFrame {
            origin_mm: [0.; 3],
            u: [1., 0., 0.],
            v: [0., 1., 0.],
        },
        kind: ViewKind::Projection,
        position_mm: [100., 80.],
        scale: Scale {
            numerator: 1,
            denominator: 2,
        },
        hidden_lines: true,
        tangent_edges: false,
        paper_chord_tolerance_mm: 0.01,
        guide: None,
        equipment_categories: None,
    }
}
#[test]
fn linked_view_roundtrips_and_invalidates_from_source_without_serializing_vectors() {
    let (mut session, id, source) = fixture();
    let book = read(session.document(), id).unwrap();
    let project = ProjectFile::new(session.document().clone()).unwrap();
    assert!(project.capabilities().contains(VIEW_CAPABILITY));
    let mut bytes = vec![];
    project_file::write(&project, &mut bytes, Default::default()).unwrap();
    let file = project_file::read(
        &mut bytes.as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap();
    assert_eq!(read(&file.project.document, id).unwrap(), book);
    let mut part = session.document().get(source).unwrap().clone();
    part.description = "Changed stock".into();
    let prepared = session
        .prepare(multicad_core::commands::Request {
            expected_revision: 0,
            actor: actor(),
            label: "Source edit".into(),
            edits: vec![Mutation::Replace(part)],
        })
        .unwrap();
    let receipt = session.commit(prepared).unwrap();
    assert!(receipt.invalidated.contains(&id));
    assert!(!serde_json::to_string(&book).unwrap().contains("preview_mm"));
}
#[test]
fn duplicate_allocates_view_identity_and_removal_keeps_capability_version_sticky() {
    let (mut session, id, _) = fixture();
    let original = read(session.document(), id).unwrap();
    let page = original.pages[0].id;
    let new_page = EntityId::new();
    let edit = prepare(
        &session,
        0,
        actor(),
        id,
        PageAction::Duplicate { page, new_page },
    )
    .unwrap();
    session.commit(edit).unwrap();
    let duplicate = read(session.document(), id).unwrap();
    assert_ne!(
        duplicate.pages[0].views[0].id,
        duplicate.pages[1].views[0].id
    );
    assert_eq!(
        duplicate.pages[0].views[0].source,
        duplicate.pages[1].views[0].source
    );
    let mut a = duplicate.pages[0].clone();
    a.views.clear();
    let e = prepare(&session, 1, actor(), id, PageAction::Replace(a)).unwrap();
    session.commit(e).unwrap();
    let e = prepare(&session, 2, actor(), id, PageAction::Remove(new_page)).unwrap();
    session.commit(e).unwrap();
    let EntityData::Drawing { definition, inputs } = &session.document().get(id).unwrap().data
    else {
        panic!()
    };
    assert_eq!(definition.schema_version, 2);
    assert!(inputs.is_empty());
    session.undo(3, actor()).unwrap();
    assert_eq!(read(session.document(), id).unwrap().pages.len(), 2);
}
#[test]
fn bad_view_frames_and_duplicate_source_ids_fail_without_source_changes() {
    let (session, id, source) = fixture();
    let mut book = read(session.document(), id).unwrap();
    book.pages[0].views[0].frame.v = [1., 0., 0.];
    assert!(replace_edits(session.document(), id, &book).is_err());
    book.pages[0].views[0] = view(source);
    book.pages[0].views[0].source = ViewSource::Model {
        entities: vec![source, source],
    };
    assert!(book.validate().is_err());
    assert_eq!(session.revision(), 0);
}
#[test]
fn stale_view_replacement_cannot_commit_and_paper_coordinates_preserve_independent_scales() {
    let (mut session, id, source) = fixture();
    let mut book = read(session.document(), id).unwrap();
    let mut b = view(source);
    b.scale = Scale {
        numerator: 2,
        denominator: 1,
    };
    b.position_mm = [150., 80.];
    assert_eq!(b.paper_point([10., 10.]), [170., 60.]);
    assert_eq!(book.pages[0].views[0].paper_point([10., 10.]), [105., 75.]);
    book.pages[0].views.push(b);
    let edit = prepare(
        &session,
        0,
        actor(),
        id,
        PageAction::Replace(book.pages[0].clone()),
    )
    .unwrap();
    let mut competing_page = read(session.document(), id).unwrap().pages[0].clone();
    competing_page.info.title = "Concurrent page edit".into();
    let other = prepare(
        &session,
        0,
        actor(),
        id,
        PageAction::Replace(competing_page),
    )
    .unwrap();
    session.commit(other).unwrap();
    assert!(session.commit(edit).is_err());
}

#[test]
fn section_guides_appear_edge_on_in_parent_views_trimmed_to_the_window() {
    use multicad_core::drawings::{SectionRange, SectionWindow, section_guide};
    let base = |u: [f64; 3], v: [f64; 3], kind: ViewKind| DrawingView {
        id: EntityId::new(),
        name: "v".into(),
        source: ViewSource::Model {
            entities: vec![EntityId::new()],
        },
        frame: multicad_core::references::ReferenceFrame {
            origin_mm: [0.; 3],
            u,
            v,
        },
        kind,
        position_mm: [0.; 2],
        scale: Scale {
            numerator: 1,
            denominator: 1,
        },
        hidden_lines: false,
        tangent_edges: false,
        paper_chord_tolerance_mm: 0.01,
        guide: None,
        equipment_categories: None,
    };
    // Top view (looking down Z) and a front section looking along +Y, cut at
    // y = 5 with a 20 mm depth and a window x in [10, 30].
    let top = base([1., 0., 0.], [0., 1., 0.], ViewKind::Projection);
    let section = base(
        [1., 0., 0.],
        [0., 0., 1.],
        ViewKind::Section {
            offset_mm: -5.,
            range: SectionRange::Depth { depth_mm: 20. },
            window: SectionWindow {
                u_mm: Some([10., 30.]),
                v_mm: None,
            },
        },
    );
    let guide = section_guide(&top, &section, ([0., 0.], [100., 50.])).unwrap();
    let mut xs = [guide.cut[0][0], guide.cut[1][0]];
    xs.sort_by(f64::total_cmp);
    let near = |a: f64, b: f64| assert!((a - b).abs() < 1e-9, "{a} != {b}");
    near(xs[0], 10.);
    near(xs[1], 30.);
    near(guide.cut[0][1], 5.);
    near(guide.cut[1][1], 5.);
    let far = guide.far.unwrap();
    near(far[0][1], 25.);
    // Without a window the guide spans the parent's bounds plus a margin.
    let unbounded = base(
        [1., 0., 0.],
        [0., 0., 1.],
        ViewKind::Section {
            offset_mm: -5.,
            range: SectionRange::Unlimited,
            window: SectionWindow::default(),
        },
    );
    let guide = section_guide(&top, &unbounded, ([0., 0.], [100., 50.])).unwrap();
    let mut xs = [guide.cut[0][0], guide.cut[1][0]];
    xs.sort_by(f64::total_cmp);
    near(xs[0], -5.);
    near(xs[1], 105.);
    assert!(guide.far.is_none());
    // A parent looking along the section normal sees no line.
    let front = base([1., 0., 0.], [0., 0., 1.], ViewKind::Projection);
    assert!(section_guide(&front, &section, ([0., 0.], [1., 1.])).is_none());
    assert!(section_guide(&top, &front, ([0., 0.], [1., 1.])).is_none());
}

#[test]
fn section_limits_keep_old_wire_form_and_validate() {
    use multicad_core::drawings::{SectionRange, SectionWindow};
    let legacy: ViewKind = serde_json::from_str(r#"{"kind":"section","offset_mm":4.0}"#).unwrap();
    assert_eq!(
        legacy,
        ViewKind::Section {
            offset_mm: 4.,
            range: SectionRange::Unlimited,
            window: SectionWindow::default()
        }
    );
    assert_eq!(
        serde_json::to_string(&legacy).unwrap(),
        r#"{"kind":"section","offset_mm":4.0}"#
    );
    let limited = ViewKind::Section {
        offset_mm: 4.,
        range: SectionRange::Depth { depth_mm: 12.5 },
        window: SectionWindow {
            u_mm: Some([-1., 2.]),
            v_mm: None,
        },
    };
    let wire = serde_json::to_string(&limited).unwrap();
    assert_eq!(serde_json::from_str::<ViewKind>(&wire).unwrap(), limited);
    assert!(
        serde_json::from_str::<ViewKind>(
            r#"{"kind":"section","offset_mm":1,"range":{"kind":"depth","depth_mm":1},"extra":1}"#
        )
        .is_err()
    );
}

#[test]
fn a_picked_parent_line_defines_a_finite_section_and_paper_points_invert() {
    use multicad_core::drawings::{SectionRange, SectionWindow, section_from_line, section_guide};
    let mut top = DrawingView {
        id: EntityId::new(),
        name: "Plan".into(),
        source: ViewSource::Model {
            entities: vec![EntityId::new()],
        },
        frame: ReferenceFrame {
            origin_mm: [10., 20., 0.],
            u: [1., 0., 0.],
            v: [0., 1., 0.],
        },
        kind: ViewKind::Projection,
        position_mm: [100., 80.],
        scale: Scale {
            numerator: 1,
            denominator: 2,
        },
        hidden_lines: false,
        tangent_edges: false,
        paper_chord_tolerance_mm: 0.01,
        guide: None,
        equipment_categories: None,
    };
    let paper = top.paper_point([12., -7.]);
    let back = top.model_point(paper);
    assert!((back[0] - 12.).abs() < 1e-12 && (back[1] + 7.).abs() < 1e-12);
    // A line along plan +X at plan y = 5, 40 mm long.
    let (frame, window) = section_from_line(&top, [0., 5.], [40., 5.]).unwrap();
    assert_eq!(frame.origin_mm, [10., 25., 0.]);
    assert_eq!(frame.u, [1., 0., 0.]);
    assert_eq!(frame.v, [0., 0., 1.]);
    assert_eq!(
        window,
        SectionWindow {
            u_mm: Some([0., 40.]),
            v_mm: None
        }
    );
    // Its guide in the plan lies back on the picked segment.
    let section = DrawingView {
        frame,
        kind: ViewKind::Section {
            offset_mm: 0.,
            range: SectionRange::Unlimited,
            window,
        },
        ..top.clone()
    };
    let guide = section_guide(&top, &section, ([-100., -100.], [100., 100.])).unwrap();
    for (end, expected) in guide.cut.iter().zip([[0., 5.], [40., 5.]]) {
        assert!(
            (end[0] - expected[0]).abs() < 1e-9 && (end[1] - expected[1]).abs() < 1e-9,
            "{end:?}"
        );
    }
    assert!(section_from_line(&top, [1., 1.], [1., 1.]).is_none());
    top.kind = ViewKind::Detail {
        center_mm: [3., 4.],
        radius_mm: 10.,
    };
    let back = top.model_point(top.paper_point([5., 6.]));
    assert!((back[0] - 5.).abs() < 1e-12 && (back[1] - 6.).abs() < 1e-12);
}

/// CAD-437: front, top and right views of the same frame origin.
fn three_views(source: EntityId) -> (Page, [EntityId; 3]) {
    let mut page = Page::new("P1");
    let mut front = view(source);
    front.frame.v = [0., 0., 1.];
    front.position_mm = [100., 150.];
    let mut top = view(source);
    top.name = "Top".into();
    top.position_mm = [100., 60.];
    let mut right = view(source);
    right.name = "Right".into();
    right.frame.u = [0., 1., 0.];
    right.frame.v = [0., 0., 1.];
    right.position_mm = [200., 150.];
    let ids = [front.id, top.id, right.id];
    page.views = vec![front, top, right];
    page.projections = vec![
        ProjectedView {
            view: ids[1],
            base: ids[0],
        },
        ProjectedView {
            view: ids[2],
            base: ids[0],
        },
    ];
    (page, ids)
}
fn at(page: &Page, id: EntityId) -> [f64; 2] {
    page.views.iter().find(|v| v.id == id).unwrap().position_mm
}

#[test]
fn moving_a_base_view_carries_its_projections_and_projections_slide_on_their_axis() {
    let (mut page, [front, top, right]) = three_views(EntityId::new());
    page.settings.paper = PaperSize::A3;
    page.settings.orientation = Orientation::Landscape;
    DrawingBook {
        pages: vec![page.clone()],
    }
    .validate()
    .unwrap();
    // The base carries both projections by the same delta.
    assert_eq!(move_view(&mut page, front, [10., -5.]).unwrap(), [10., -5.]);
    assert_eq!(at(&page, front), [110., 145.]);
    assert_eq!(at(&page, top), [110., 55.]);
    assert_eq!(at(&page, right), [210., 145.]);
    // A top view keeps the front's column; a side view keeps its row.
    assert_eq!(move_view(&mut page, top, [30., -20.]).unwrap(), [0., -20.]);
    assert_eq!(at(&page, top), [110., 35.]);
    assert_eq!(move_view(&mut page, right, [40., 25.]).unwrap(), [40., 0.]);
    assert_eq!(at(&page, right), [250., 145.]);
    // The whole group stops at the printable area instead of leaving it.
    let margin = page.settings.margin_mm;
    move_view(&mut page, front, [-1000., 0.]).unwrap();
    assert!((at(&page, front)[0] - margin).abs() < 1e-9);
    assert_eq!(at(&page, top)[0], at(&page, front)[0]);
    assert!(move_view(&mut page, front, [f64::NAN, 0.]).is_err());
    // Scale and origin changes on the base re-align the projections.
    page.views[0].scale = Scale {
        numerator: 1,
        denominator: 5,
    };
    page.views[1].frame.origin_mm = [50., 0., 0.];
    align_projections(&mut page);
    assert_eq!(page.views[1].scale, page.views[0].scale);
    assert!((at(&page, top)[0] - (at(&page, front)[0] + 10.)).abs() < 1e-9);
    // Removing the base frees its projections.
    remove_view(&mut page, front);
    assert!(page.projections.is_empty());
}

#[test]
fn projection_links_validate_roundtrip_with_their_capability_and_survive_duplicate() {
    let mut doc = Document::new("Projected").unwrap();
    let part = Entity::new(
        "Part",
        EntityData::Part {
            features: vec![],
            material: None,
            tags: Default::default(),
        },
    );
    let source = part.id;
    doc.insert(part).unwrap();
    let (page, [front, top, _]) = three_views(source);
    // Invalid links never enter the document.
    for bad in [
        ProjectedView {
            view: top,
            base: EntityId::new(),
        },
        ProjectedView {
            view: front,
            base: top,
        },
    ] {
        let mut p = page.clone();
        p.projections.push(bad);
        assert!(DrawingBook { pages: vec![p] }.validate().is_err());
    }
    let mut same_direction = page.clone();
    same_direction.views[1].frame = same_direction.views[0].frame;
    assert!(
        DrawingBook {
            pages: vec![same_direction]
        }
        .validate()
        .is_err()
    );
    let book = DrawingBook { pages: vec![page] }.entity("Drawing").unwrap();
    let id = book.id;
    doc.insert(book).unwrap();
    let mut session = Session::new(doc).unwrap();
    let project = ProjectFile::new(session.document().clone()).unwrap();
    assert!(project.capabilities().contains(PROJECTION_CAPABILITY));
    let mut bytes = vec![];
    project_file::write(&project, &mut bytes, Default::default()).unwrap();
    let mut older = project_file::supported_capabilities();
    older.remove(PROJECTION_CAPABILITY);
    assert!(project_file::read(&mut bytes.as_slice(), Default::default(), &older).is_err());
    let reopened = project_file::read(
        &mut bytes.as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap();
    assert_eq!(
        read(&reopened.project.document, id).unwrap(),
        read(session.document(), id).unwrap()
    );
    // A duplicated page links its own copies.
    let first = read(session.document(), id).unwrap().pages[0].id;
    let copy = EntityId::new();
    let prepared = prepare(
        &session,
        session.revision(),
        actor(),
        id,
        PageAction::Duplicate {
            page: first,
            new_page: copy,
        },
    )
    .unwrap();
    session.commit(prepared).unwrap();
    let book = read(session.document(), id).unwrap();
    let ids: Vec<_> = book.pages[1].views.iter().map(|v| v.id).collect();
    assert_eq!(book.pages[1].projections.len(), 2);
    assert!(
        book.pages[1]
            .projections
            .iter()
            .all(|l| ids.contains(&l.view) && ids.contains(&l.base))
    );
    // Without projections the wire form is unchanged.
    let wire = serde_json::to_string(&Page::new("x")).unwrap();
    assert!(!wire.contains("projections"));
}

/// CAD-436: the standard layout picks the largest standard scale that fits
/// the sheet, aligns top with front's column and right with its row, and
/// follows the page's projection angle.
#[test]
fn standard_layout_fits_the_sheet_and_aligns_projections() {
    // A 400 × 200 × 100 mm block: front 400 × 100, top 400 × 200, right 200 × 100.
    let extents = [[400., 100.], [400., 200.], [200., 100.], [420., 300.]];
    let mut chosen = vec![];
    for paper in [PaperSize::A4, PaperSize::A3, PaperSize::A1] {
        let mut page = Page::new("1");
        page.settings.paper = paper;
        for projection in [Projection::ThirdAngle, Projection::FirstAngle] {
            page.settings.projection = projection;
            let layout = layout_page(&page, 1, 1).unwrap();
            let arranged = standard_layout(&layout, projection, extents).unwrap();
            let [front, top, right, iso] = arranged.positions;
            assert_eq!(front[0], top[0], "top shares the front's column");
            assert_eq!(front[1], right[1], "right shares the front's row");
            assert!(iso[0] == right[0] && iso[1] == top[1]);
            let third = projection == Projection::ThirdAngle;
            assert_eq!(top[1] < front[1], third, "top above front in third angle");
            assert_eq!(right[0] > front[0], third, "right of front in third angle");
            // The chosen scale fits; the next larger standard scale does not.
            let fits = |scale: Scale| {
                let s = scale_factor(scale);
                let w = (400. + 420.) * s + VIEW_GAP_MM;
                let h = (100. + 300.) * s + VIEW_GAP_MM;
                fit_regions(&layout)
                    .iter()
                    .any(|r| w <= r[2] - r[0] + 1e-9 && h <= r[3] - r[1] + 1e-9)
            };
            assert!(fits(arranged.scale));
            let scales = standard_scales();
            let index = scales.iter().position(|s| *s == arranged.scale).unwrap();
            assert!(index == 0 || !fits(scales[index - 1]));
            // Everything stays on the printable paper.
            let [x0, y0, x1, y1] = layout.printable_mm;
            for p in arranged.positions {
                assert!(p[0] > x0 && p[0] < x1 && p[1] > y0 && p[1] < y1);
            }
        }
        chosen.push(arranged_scale(paper, extents));
    }
    assert_eq!(
        chosen,
        [
            Scale {
                numerator: 1,
                denominator: 5
            },
            Scale {
                numerator: 1,
                denominator: 5
            },
            Scale {
                numerator: 1,
                denominator: 2
            },
        ]
    );
    assert!(
        standard_layout(
            &layout_page(&Page::new("1"), 1, 1).unwrap(),
            Projection::FirstAngle,
            [[f64::NAN, 1.], [1., 1.], [1., 1.], [1., 1.]]
        )
        .is_err()
    );
}
fn arranged_scale(paper: PaperSize, extents: [[f64; 2]; 4]) -> Scale {
    let mut page = Page::new("1");
    page.settings.paper = paper;
    standard_layout(
        &layout_page(&page, 1, 1).unwrap(),
        Projection::FirstAngle,
        extents,
    )
    .unwrap()
    .scale
}

/// CAD-436: a single guided view takes the largest standard scale that fits
/// the free sheet with margins, so a larger sheet never gets a smaller scale.
#[test]
fn single_view_default_scale_is_the_largest_standard_scale_fitting_the_sheet() {
    let source = EntityId::new();
    let mut previous = 0.;
    for (paper, expected) in [
        (PaperSize::A4, (1, 2)),
        (PaperSize::A3, (1, 2)),
        (PaperSize::A2, (1, 1)),
        (PaperSize::A0, (2, 1)),
    ] {
        let mut page = Page::new("1");
        page.settings.paper = paper;
        let layout = layout_page(&page, 1, 1).unwrap();
        let fit = fit_view(
            &view(source),
            &layout,
            [-200., -50.],
            [200., 50.],
            FitMode::Standard,
        )
        .unwrap();
        assert_eq!(
            (fit.scale.numerator, fit.scale.denominator),
            expected,
            "{paper:?}"
        );
        let f = scale_factor(fit.scale);
        assert!(f >= previous);
        previous = f;
        assert!(
            fit_regions(&layout)
                .iter()
                .any(|r| 400. * f <= r[2] - r[0] && 100. * f <= r[3] - r[1])
        );
    }
}

/// CAD-438: a view shows every part until its page stores a filter; tags
/// and per-part toggles decide, invalid filters are refused, a duplicated
/// page keeps its filters and pages without filters keep their wire form.
#[test]
fn view_part_filters_follow_tags_and_toggles_and_validate() {
    let mut doc = Document::new("Filters").unwrap();
    let part = |name: &str, tags: &[&str]| {
        Entity::new(
            name,
            EntityData::Part {
                features: vec![],
                material: None,
                tags: tags.iter().map(|t| t.to_string()).collect(),
            },
        )
    };
    let (frame, bolt) = (part("Frame", &["frame"]), part("Bolt", &["hardware"]));
    let (f, b) = (frame.id, bolt.id);
    doc.insert(frame).unwrap();
    doc.insert(bolt).unwrap();
    let mut page = Page::new("1");
    let mut v = view(f);
    v.source = ViewSource::Model {
        entities: vec![f, b],
    };
    let id = v.id;
    page.views.push(v);
    assert_eq!(view_parts(&doc, &page.views[0]), [f, b]);
    assert!(page.part_filter(id).is_none(), "everything by default");
    let shows = |page: &Page, part| page.part_filter(id).is_none_or(|x| x.shows(&doc, part));
    page.set_tag_filter(id, ["frame".into()].into(), Default::default())
        .unwrap();
    assert!(shows(&page, f) && !shows(&page, b));
    page.set_tag_filter(id, Default::default(), ["frame".into()].into())
        .unwrap();
    assert!(!shows(&page, f) && shows(&page, b));
    page.set_part_shown(&doc, id, f, true).unwrap();
    assert!(shows(&page, f), "a toggle overrides the tag rule");
    page.set_part_shown(&doc, id, b, false).unwrap();
    assert!(!shows(&page, b));
    DrawingBook {
        pages: vec![page.clone()],
    }
    .validate()
    .unwrap();
    // Clearing everything removes the stored filter.
    page.set_tag_filter(id, Default::default(), Default::default())
        .unwrap();
    page.set_part_shown(&doc, id, f, true).unwrap();
    page.set_part_shown(&doc, id, b, true).unwrap();
    assert!(page.part_filters.is_empty());
    assert!(
        page.set_part_shown(&doc, EntityId::new(), f, false)
            .is_err()
    );
    assert!(
        page.set_tag_filter(id, [" padded".into()].into(), Default::default())
            .is_err()
    );
    for bad in [
        ViewPartFilter {
            view: EntityId::new(),
            hidden: [f].into(),
            ..Default::default()
        },
        ViewPartFilter {
            view: id,
            ..Default::default()
        },
        ViewPartFilter {
            view: id,
            hidden: [f].into(),
            shown: [f].into(),
            ..Default::default()
        },
    ] {
        let mut p = page.clone();
        p.part_filters = vec![bad];
        assert!(DrawingBook { pages: vec![p] }.validate().is_err());
    }
    // Stored, duplicated with the page and required by readers.
    page.set_part_shown(&doc, id, b, false).unwrap();
    let book = DrawingBook {
        pages: vec![page.clone()],
    }
    .entity("Drawing")
    .unwrap();
    let book_id = book.id;
    doc.insert(book).unwrap();
    let mut session = Session::new(doc).unwrap();
    let copy = EntityId::new();
    let prepared = prepare(
        &session,
        session.revision(),
        actor(),
        book_id,
        PageAction::Duplicate {
            page: page.id,
            new_page: copy,
        },
    )
    .unwrap();
    session.commit(prepared).unwrap();
    let book = read(session.document(), book_id).unwrap();
    let duplicate = &book.pages[1];
    assert_eq!(duplicate.part_filters.len(), 1);
    assert_eq!(duplicate.part_filters[0].view, duplicate.views[0].id);
    let project = ProjectFile::new(session.document().clone()).unwrap();
    assert!(project.capabilities().contains(PART_FILTER_CAPABILITY));
    let mut bytes = vec![];
    project_file::write(&project, &mut bytes, Default::default()).unwrap();
    let mut older = project_file::supported_capabilities();
    older.remove(PART_FILTER_CAPABILITY);
    assert!(project_file::read(&mut bytes.as_slice(), Default::default(), &older).is_err());
    assert!(
        !serde_json::to_string(&Page::new("x"))
            .unwrap()
            .contains("part_filters")
    );
}
