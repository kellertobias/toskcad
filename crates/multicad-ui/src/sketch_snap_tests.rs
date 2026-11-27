use super::*;
use multicad_core::constraint_sketch::{Curve, Point};
use multicad_core::references::TopologyKind;

/// Screen pixels per sketch millimeter in these tests.
const SCALE: f32 = 4.;

fn to_screen(p: [f64; 2]) -> Option<[f32; 2]> {
    Some([p[0] as f32 * SCALE, -(p[1] as f32) * SCALE])
}
fn query(pointer_mm: [f64; 2]) -> Query<'static> {
    Query {
        pointer_mm,
        pointer_px: to_screen(pointer_mm).unwrap(),
        to_screen: &to_screen,
        from: None,
    }
}
/// Pointer offset by whole screen pixels from a sketch position.
fn query_px(at: [f64; 2], dx: f32, dy: f32) -> Query<'static> {
    let mut pointer = to_screen(at).unwrap();
    pointer[0] += dx;
    pointer[1] += dy;
    Query {
        pointer_mm: [(pointer[0] / SCALE) as f64, (-pointer[1] / SCALE) as f64],
        pointer_px: pointer,
        to_screen: &to_screen,
        from: None,
    }
}

#[derive(Default)]
struct Builder {
    sketch: Sketch,
}
impl Builder {
    fn point(&mut self, p: [f64; 2]) -> EntityId {
        let id = EntityId::new();
        self.sketch.points.push(Point { id, position_mm: p });
        id
    }
    fn curve(&mut self, kind: CurveKind) -> EntityId {
        let id = EntityId::new();
        self.sketch.curves.push(Curve {
            id,
            construction: false,
            kind,
        });
        id
    }
    fn line(&mut self, start: EntityId, end: EntityId) -> EntityId {
        self.curve(CurveKind::Line { start, end })
    }
    fn polygon(&mut self, vertices: &[[f64; 2]]) -> (Vec<EntityId>, Vec<EntityId>) {
        let points: Vec<_> = vertices.iter().map(|v| self.point(*v)).collect();
        let lines = (0..points.len())
            .map(|k| self.line(points[k], points[(k + 1) % points.len()]))
            .collect();
        (points, lines)
    }
    fn geometry(&self) -> SnapGeometry {
        SnapGeometry::from_sketch(&self.sketch)
    }
}

/// Hovers at `pointer` at each listed time and returns the last outcome.
fn hover(assist: &mut SnapAssist, g: &SnapGeometry, q: &Query, times: &[f64]) -> Outcome {
    let mut out = Outcome::default();
    for t in times {
        out = assist.update(g, q, *t);
    }
    out
}

#[test]
fn dwell_acquires_only_after_the_documented_duration() {
    let mut b = Builder::default();
    let a = b.point([10., 10.]);
    let g = b.geometry();
    let mut assist = SnapAssist::default();
    let q = query_px([10., 10.], 2., -1.);
    let out = hover(&mut assist, &g, &q, &[0., 0.2, DWELL_SECONDS - 0.01]);
    assert!(out.event.is_none());
    assert!(assist.references().is_empty());
    let remaining = out.repaint_after.unwrap();
    assert!((remaining - 0.01).abs() < 1e-9);
    let out = assist.update(&g, &q, DWELL_SECONDS);
    assert_eq!(out.event, Some(DwellEvent::Acquired(0)));
    let r = &assist.references()[0];
    assert_eq!(r.position, [10., 10.], "cross sits at the exact point");
    assert_eq!(r.source, TargetKey::Point(a));
    // Continuing to hover does not toggle it off again.
    assert!(assist.update(&g, &q, 3.).event.is_none());
    assert_eq!(assist.references().len(), 1);
}

#[test]
fn briefly_passing_over_or_drifting_does_not_acquire() {
    let mut b = Builder::default();
    b.point([0., 0.]);
    let g = b.geometry();
    let mut assist = SnapAssist::default();
    assist.update(&g, &query_px([0., 0.], 0., 0.), 0.);
    assist.update(&g, &query([40., 40.]), 0.3);
    let out = assist.update(&g, &query_px([0., 0.], 0., 0.), 0.6);
    assert!(out.event.is_none(), "timer restarts after leaving");
    // Drift beyond the stability tolerance restarts the timer as well.
    assist.update(&g, &query([40., 40.]), 1.9);
    assist.update(&g, &query_px([0., 0.], -3., 0.), 2.);
    assist.update(&g, &query_px([0., 0.], 3., 0.), 2.3);
    let out = assist.update(&g, &query_px([0., 0.], 3., 0.), 2.6);
    assert!(
        out.event.is_none(),
        "6 px drift restarted the dwell at 2.3 s"
    );
    assert!(assist.references().is_empty());
    let out = assist.update(&g, &query_px([0., 0.], 2., 0.), 2.3 + DWELL_SECONDS);
    assert_eq!(out.event, Some(DwellEvent::Acquired(0)));
}

fn acquire(assist: &mut SnapAssist, g: &SnapGeometry, at: [f64; 2], start: f64) {
    let q = query_px(at, 1., 1.);
    let out = hover(assist, g, &q, &[start, start + DWELL_SECONDS]);
    assert!(
        matches!(out.event, Some(DwellEvent::Acquired(_))),
        "{at:?} was not acquired: {out:?}"
    );
}

#[test]
fn two_references_expose_their_exact_midpoint() {
    let mut b = Builder::default();
    b.point([0., 0.]);
    b.point([20., 10.]);
    let g = b.geometry();
    let mut assist = SnapAssist::default();
    acquire(&mut assist, &g, [0., 0.], 0.);
    acquire(&mut assist, &g, [20., 10.], 1.);
    assert_eq!(assist.references().len(), 2, "first reference retained");
    let out = assist.update(&g, &query_px([10., 5.], 3., -2.), 3.);
    let target = out.target.unwrap();
    assert_eq!(target.kind, TargetKind::ReferenceMidpoint);
    assert_eq!(target.position, [10., 5.]);
    assert_eq!(
        placement(Some(&target), [0.; 2], Some(10.)),
        ([10., 5.], None)
    );
}

#[test]
fn guides_snap_in_all_four_directions() {
    let mut b = Builder::default();
    b.point([5., 5.]);
    let g = b.geometry();
    let mut assist = SnapAssist::default();
    acquire(&mut assist, &g, [5., 5.], 0.);
    for (pointer, axis, expected) in [
        ([35., 5.6], Axis::Horizontal, [35., 5.]),
        ([-25., 4.4], Axis::Horizontal, [-25., 5.]),
        ([5.7, 40.], Axis::Vertical, [5., 40.]),
        ([4.5, -30.], Axis::Vertical, [5., -30.]),
    ] {
        let target = assist.update(&g, &query(pointer), 2.).target.unwrap();
        assert_eq!(target.kind, TargetKind::Guide(axis), "{pointer:?}");
        assert!((target.position[0] - expected[0]).abs() < 1e-6);
        assert!((target.position[1] - expected[1]).abs() < 1e-6);
        assert_eq!(
            target.position[if axis == Axis::Horizontal { 1 } else { 0 }],
            5.
        );
        assert_eq!(
            target.guides,
            vec![GuideLine {
                from: [5., 5.],
                axis
            }]
        );
    }
    // Beyond the guide tolerance (6 px = 1.5 mm here) nothing snaps.
    assert!(assist.update(&g, &query([35., 7.]), 2.).target.is_none());
}

#[test]
fn guides_from_two_references_intersect_and_meet_lines() {
    let mut b = Builder::default();
    b.point([0., 0.]);
    b.point([30., 20.]);
    let s = b.point([50., -10.]);
    let e = b.point([50., 40.]);
    b.line(s, e);
    let g = b.geometry();
    let mut assist = SnapAssist::default();
    acquire(&mut assist, &g, [0., 0.], 0.);
    acquire(&mut assist, &g, [30., 20.], 1.);
    let target = assist.update(&g, &query([30.5, 0.4]), 2.).target.unwrap();
    assert_eq!(target.kind, TargetKind::GuideIntersection);
    assert_eq!(target.position, [30., 0.]);
    assert_eq!(target.guides.len(), 2);
    let target = assist.update(&g, &query([49.6, 20.3]), 2.).target.unwrap();
    assert_eq!(target.kind, TargetKind::GuideOnLine);
    assert_eq!(target.position, [50., 20.]);
}

#[test]
fn line_midpoints_points_on_lines_and_competing_targets() {
    let mut b = Builder::default();
    let a = b.point([0., 0.]);
    let c = b.point([40., 0.]);
    let line = b.line(a, c);
    let g = b.geometry();
    let target = resolve(&g, &[], &query_px([20., 0.], 2., 1.), None).unwrap();
    assert_eq!(target.key, TargetKey::LineMidpoint(line));
    assert_eq!(target.position, [20., 0.]);
    let target = resolve(&g, &[], &query([31., 0.8]), None).unwrap();
    assert_eq!(target.kind, TargetKind::OnCurve);
    assert!((target.position[0] - 31.).abs() < 1e-6 && target.position[1] == 0.);
    // A vertex beats the nearer point on its own line and midpoints.
    let target = resolve(&g, &[], &query_px([40., 0.], -6., 0.), None).unwrap();
    assert_eq!(target.key, TargetKey::Point(c));
    assert_eq!(target.anchor, Some(c));
    // Outside every tolerance there is no geometry target; the grid applies.
    assert!(resolve(&g, &[], &query([25., 10.]), None).is_none());
    assert_eq!(placement(None, [24.6, 10.2], Some(10.)), ([20., 10.], None));
    assert_eq!(placement(None, [24.6, 10.2], None), ([24.6, 10.2], None));
}

#[test]
fn a_nearby_grid_point_never_displaces_a_geometry_target() {
    let mut b = Builder::default();
    let a = b.point([0., 1.3]);
    let c = b.point([40., 1.3]);
    b.line(a, c);
    let g = b.geometry();
    // The pointer sits on the grid point (10, 0) but within 6 px of the line.
    let q = query([10., 0.1]);
    let target = resolve(&g, &[], &q, None).unwrap();
    for grid in [Some(10.), None] {
        let (p, anchor) = placement(Some(&target), q.pointer_mm, grid);
        assert!((p[0] - 10.).abs() < 1e-6);
        assert_eq!(p[1], 1.3);
        assert_eq!(anchor, None);
    }
}

#[test]
fn hysteresis_keeps_the_previous_target_until_the_release_distance() {
    let mut b = Builder::default();
    let a = b.point([0., 0.]);
    let g = b.geometry();
    let key = TargetKey::Point(a);
    let q = query_px([0., 0.], 10., 0.);
    assert!(resolve(&g, &[], &q, None).is_none());
    assert_eq!(resolve(&g, &[], &q, Some(&key)).unwrap().key, key);
    let q = query_px([0., 0.], POINT_SNAP_PX * RELEASE_FACTOR + 0.5, 0.);
    assert!(resolve(&g, &[], &q, Some(&key)).is_none());
}

#[test]
fn centers_are_defined_only_for_circles_arcs_rectangles_and_regular_polygons() {
    let mut b = Builder::default();
    b.polygon(&[[0., 0.], [40., 0.], [40., 20.], [0., 20.]]);
    let hexagon: Vec<[f64; 2]> = (0..6)
        .map(|k| {
            let a = k as f64 * TAU / 6.;
            [100. + 10. * a.portable_cos(), 10. * a.portable_sin()]
        })
        .collect();
    b.polygon(&hexagon);
    // Irregular quadrilateral and scalene triangle: no center target.
    b.polygon(&[[200., 0.], [240., 0.], [230., 20.], [205., 25.]]);
    b.polygon(&[[300., 0.], [340., 0.], [310., 30.]]);
    let center = b.point([0., 100.]);
    b.curve(CurveKind::Circle {
        center,
        radius_mm: 5.,
    });
    let g = b.geometry();
    let mut centers: Vec<_> = g.shape_centers().collect();
    centers.sort_by(|a, b| a.0[0].total_cmp(&b.0[0]));
    assert_eq!(centers.len(), 2, "{centers:?}");
    assert_eq!(centers[0], ([20., 10.], ShapeKind::Rectangle));
    assert!((centers[1].0[0] - 100.).abs() < 1e-9 && centers[1].0[1].abs() < 1e-9);
    assert_eq!(centers[1].1, ShapeKind::RegularPolygon(6));
    let target = resolve(&g, &[], &query_px([20., 10.], 2., 2.), None).unwrap();
    assert_eq!(target.kind, TargetKind::ShapeCenter(ShapeKind::Rectangle));
    assert_eq!(target.position, [20., 10.]);
    let target = resolve(&g, &[], &query_px([0., 100.], 1., 0.), None).unwrap();
    assert_eq!(target.kind, TargetKind::Center);
    assert_eq!(target.anchor, Some(center));
    assert!(resolve(&g, &[], &query([220., 12.]), None).is_none());
}

#[test]
fn rectangle_with_an_attached_line_keeps_its_center() {
    let mut b = Builder::default();
    let (points, _) = b.polygon(&[[0., 0.], [10., 0.], [10., 10.], [0., 10.]]);
    let tail = b.point([30., 30.]);
    b.line(points[2], tail);
    let centers: Vec<_> = b.geometry().shape_centers().collect();
    assert_eq!(centers, vec![([5., 5.], ShapeKind::Rectangle)]);
}

#[test]
fn arcs_snap_only_within_their_sweep() {
    let mut b = Builder::default();
    let center = b.point([0., 0.]);
    let start = b.point([10., 0.]);
    let end = b.point([0., 10.]);
    b.curve(CurveKind::Arc {
        center,
        start,
        end,
        clockwise: false,
    });
    let g = b.geometry();
    let s = std::f64::consts::FRAC_1_SQRT_2 * 10.;
    let target = resolve(&g, &[], &query([s + 0.5, s + 0.5]), None).unwrap();
    assert_eq!(target.kind, TargetKind::OnCurve);
    assert!(((target.position[0]).portable_hypot(target.position[1]) - 10.).abs() < 1e-9);
    assert!(resolve(&g, &[], &query([-s - 0.5, -s - 0.5]), None).is_none());
}

#[test]
fn line_midpoints_and_centers_can_be_acquired_and_aligned() {
    let mut b = Builder::default();
    b.polygon(&[[0., 0.], [40., 0.], [40., 20.], [0., 20.]]);
    let g = b.geometry();
    let mut assist = SnapAssist::default();
    acquire(&mut assist, &g, [20., 10.], 0.);
    assert_eq!(
        assist.references()[0].kind,
        TargetKind::ShapeCenter(ShapeKind::Rectangle)
    );
    acquire(&mut assist, &g, [40., 10.], 1.);
    assert_eq!(assist.references()[1].kind, TargetKind::LineMidpoint);
    let target = assist.update(&g, &query([80., 10.6]), 2.).target.unwrap();
    assert_eq!(target.kind, TargetKind::Guide(Axis::Horizontal));
    assert_eq!(target.position[1], 10.);
    let target = assist.update(&g, &query([20.6, 60.]), 2.).target.unwrap();
    assert_eq!(target.kind, TargetKind::Guide(Axis::Vertical));
    assert_eq!(target.position[0], 20.);
}

#[test]
fn on_curve_references_persist_and_redwelling_releases() {
    let mut b = Builder::default();
    let a = b.point([0., 0.]);
    let c = b.point([40., 0.]);
    b.line(a, c);
    let g = b.geometry();
    let mut assist = SnapAssist::default();
    let q = query([13., 0.5]);
    let out = hover(&mut assist, &g, &q, &[0., DWELL_SECONDS, 2., 4.]);
    assert!(out.event.is_none(), "acquired point must not flip off");
    assert_eq!(assist.references().len(), 1);
    assert_eq!(assist.references()[0].position, [13., 0.]);
    assert_eq!(out.target.unwrap().kind, TargetKind::Reference);
    // Leave and dwell again on the reference: it is released.
    assist.update(&g, &query([20., 30.]), 5.);
    let out = hover(&mut assist, &g, &q, &[6., 6. + DWELL_SECONDS]);
    assert_eq!(out.event, Some(DwellEvent::Released(0)));
    assert!(assist.references().is_empty());
}

#[test]
fn a_third_reference_drops_the_oldest() {
    let mut b = Builder::default();
    for x in [0., 30., 60.] {
        b.point([x, x / 3.]);
    }
    let g = b.geometry();
    let mut assist = SnapAssist::default();
    for (k, x) in [0., 30., 60.].into_iter().enumerate() {
        acquire(&mut assist, &g, [x, x / 3.], k as f64);
    }
    let positions: Vec<_> = assist.references().iter().map(|r| r.position).collect();
    assert_eq!(positions, vec![[30., 10.], [60., 20.]]);
}

#[test]
fn references_follow_their_source_and_vanish_with_it() {
    let mut b = Builder::default();
    let a = b.point([0., 0.]);
    let mut assist = SnapAssist::default();
    acquire(&mut assist, &b.geometry(), [0., 0.], 0.);
    b.sketch.points[0].position_mm = [5., 5.];
    assist.revalidate(&b.geometry());
    assert_eq!(assist.references()[0].position, [5., 5.]);
    b.sketch.points.retain(|p| p.id != a);
    assist.revalidate(&b.geometry());
    assert!(assist.references().is_empty());
}

#[test]
fn clearing_and_switching_sketches_remove_transient_aids() {
    let mut b = Builder::default();
    b.point([0., 0.]);
    let g = b.geometry();
    let mut assist = SnapAssist::default();
    let sketch = EntityId::new();
    assist.bind(sketch);
    acquire(&mut assist, &g, [0., 0.], 0.);
    assist.bind(sketch);
    assert_eq!(assist.references().len(), 1);
    assist.bind(EntityId::new());
    assert!(assist.references().is_empty());
    acquire(&mut assist, &g, [0., 0.], 1.);
    assist.clear();
    assert!(assist.references().is_empty() && assist.shown().is_none());
}

#[test]
fn tolerances_are_screen_pixels_and_follow_zoom() {
    let mut b = Builder::default();
    b.point([0., 0.]);
    let g = b.geometry();
    // 1 mm away: 4 px at the test scale snaps; zoomed out 10x it is 0.4 px;
    // zoomed in so 1 mm is 12 px it no longer snaps.
    for (scale, snaps) in [(4f32, true), (0.4, true), (12., false)] {
        let screen = move |p: [f64; 2]| Some([p[0] as f32 * scale, -(p[1] as f32) * scale]);
        let q = Query {
            pointer_mm: [1., 0.],
            pointer_px: [scale, 0.],
            to_screen: &screen,
            from: None,
        };
        assert_eq!(resolve(&g, &[], &q, None).is_some(), snaps, "scale {scale}");
    }
}

/// `query_px` with the start of the line being drawn.
fn from_query(at: [f64; 2], dx: f32, dy: f32, from: [f64; 2]) -> Query<'static> {
    Query {
        from: Some(from),
        ..query_px(at, dx, dy)
    }
}

#[test]
fn perpendicular_foot_on_a_line_needs_a_line_start_and_lies_on_the_segment() {
    let mut b = Builder::default();
    let (p, q) = (b.point([0., 0.]), b.point([40., 30.]));
    b.line(p, q);
    let g = b.geometry();
    let from = [5., 20.];
    // Foot of the perpendicular from (5,20) onto the line through (0,0) and
    // (40,30): t = (5*40 + 20*30) / 2500 = 0.32 -> (12.8, 9.6).
    let foot = [12.8, 9.6];
    let target = resolve(&g, &[], &from_query(foot, 3., -2., from), None).unwrap();
    assert_eq!(target.kind, TargetKind::Perpendicular);
    let ab = [40., 30.];
    let drawn = sub(target.position, from);
    assert!(
        dot(drawn, ab).abs() < 1e-9,
        "not perpendicular: {:?}",
        target.position
    );
    assert!(len(sub(target.position, foot)) < 1e-9);
    assert_eq!(target.anchor, None);
    // It beats the on-curve point, but an endpoint still wins over it.
    let near = resolve(&g, &[], &query_px(foot, 3., -2.), None).unwrap();
    assert_eq!(near.kind, TargetKind::OnCurve, "no line start, no foot");
    // A foot beyond the segment is not offered.
    let target = resolve(&g, &[], &from_query([44., 33.], 2., 0., [60., 10.]), None);
    assert!(target.is_none_or(|t| t.kind != TargetKind::Perpendicular));
}

#[test]
fn circles_and_arcs_offer_radial_feet_and_tangent_points_from_the_line_start() {
    let mut b = Builder::default();
    let c = b.point([0., 0.]);
    let circle = b.curve(CurveKind::Circle {
        center: c,
        radius_mm: 10.,
    });
    let g = b.geometry();
    let from = [30., 0.];
    let near = resolve(&g, &[], &from_query([10., 0.], 2., 3., from), None).unwrap();
    assert_eq!(near.kind, TargetKind::Perpendicular);
    assert_eq!(near.key, TargetKey::Perpendicular(circle, 0));
    assert!(len(sub(near.position, [10., 0.])) < 1e-12);
    let far = resolve(&g, &[], &from_query([-10., 0.], 2., 3., from), None).unwrap();
    assert_eq!(far.key, TargetKey::Perpendicular(circle, 1));
    // Tangent from (30,0): the radius to the tangent point is perpendicular
    // to the drawn line.
    let angle = (10f64 / 30.).portable_acos();
    let tangent = [10. * angle.portable_cos(), 10. * angle.portable_sin()];
    let target = resolve(&g, &[], &from_query(tangent, -2., 2., from), None).unwrap();
    assert_eq!(target.kind, TargetKind::Tangent);
    assert!(dot(target.position, sub(target.position, from)).abs() < 1e-9);
    // A line start inside the circle has no tangent.
    assert!(
        relative_targets(&g, [1., 0.])
            .iter()
            .all(|t| t.1 != TargetKind::Tangent)
    );
    // A quarter arc only offers targets inside its sweep.
    let mut b = Builder::default();
    let (c, s, e) = (b.point([0., 0.]), b.point([10., 0.]), b.point([0., 10.]));
    b.curve(CurveKind::Arc {
        center: c,
        start: s,
        end: e,
        clockwise: false,
    });
    let g = b.geometry();
    let targets = relative_targets(&g, [20., 20.]);
    let feet: Vec<_> = targets
        .iter()
        .filter(|t| t.1 == TargetKind::Perpendicular)
        .collect();
    assert_eq!(feet.len(), 1, "only the near foot is on the arc");
    let h = 50f64.sqrt();
    assert!(len(sub(feet[0].2, [h, h])) < 1e-9);
    // From the opposite side only the far foot lies on the arc; both
    // tangent points fall outside its sweep.
    let behind = relative_targets(&g, [-20., -20.]);
    assert_eq!(behind.len(), 1, "{behind:?}");
    assert!(matches!(behind[0].0, TargetKey::Perpendicular(_, 1)));
    assert!(len(sub(behind[0].2, [h, h])) < 1e-9);
}

fn support_key(kind: TopologyKind, ordinal: u32) -> TopologyKey {
    TopologyKey { kind, ordinal }
}
#[test]
fn support_face_corners_and_edges_snap_exactly_and_existing_sketch_geometry_wins() {
    let mut b = Builder::default();
    let existing = b.point([30., 20.]);
    let corner = support_key(TopologyKind::Vertex, 0);
    let edge = support_key(TopologyKind::Edge, 1);
    let arc = support_key(TopologyKind::Edge, 2);
    let geometry = b
        .geometry()
        .with_support(vec![
            (corner, [0., 0.]),
            // Already in the sketch: reused by identity, not offered again.
            (support_key(TopologyKind::Vertex, 3), [30., 20.]),
        ])
        .with_support_edges(vec![
            (edge, SupportEdge::Line([0., -10.], [30., -10.])),
            (
                arc,
                SupportEdge::Arc {
                    center: [60., 0.],
                    radius: 5.,
                    start: 0.,
                    sweep: std::f64::consts::PI,
                },
            ),
        ]);

    let target = resolve(&geometry, &[], &query_px([0., 0.], 3., 2.), None).unwrap();
    assert_eq!(target.key, TargetKey::Support(corner));
    assert_eq!(target.kind, TargetKind::FaceCorner);
    assert_eq!(
        target.position,
        [0., 0.],
        "the exact face corner, not the pointer"
    );
    assert_eq!(
        target.anchor, None,
        "the projected point does not exist yet"
    );
    assert_eq!(target.kind.marker(), Marker::Vertex);

    let target = resolve(&geometry, &[], &query_px([12., -10.], 0., 3.), None).unwrap();
    assert_eq!(target.key, TargetKey::Support(edge));
    assert_eq!(target.kind, TargetKind::FaceEdge);
    assert_eq!(target.kind.label(), "Face edge");
    assert_eq!(target.position[1], -10., "exactly on the face edge");
    assert_eq!(target.kind.marker(), Marker::Edge);

    let on_arc = resolve(&geometry, &[], &query([60., 5.2]), None).unwrap();
    assert_eq!(on_arc.key, TargetKey::Support(arc));
    assert!((on_arc.position[1] - 5.).abs() < 1e-12);
    // Outside the arc's sweep there is no target.
    assert!(resolve(&geometry, &[], &query([60., -5.2]), None).is_none());

    // The existing point keeps its identity.
    let target = resolve(&geometry, &[], &query([30., 20.]), None).unwrap();
    assert_eq!(target.anchor, Some(existing));
    // A face corner beats an edge of the face nearby (point class first).
    let geometry = geometry.with_support_edges(vec![(
        support_key(TopologyKind::Edge, 4),
        SupportEdge::Line([-10., 0.], [10., 0.]),
    )]);
    let target = resolve(&geometry, &[], &query_px([0., 0.], 2., 0.), None).unwrap();
    assert_eq!(target.key, TargetKey::Support(corner));
    // An edge already drawn in the sketch is not offered again.
    let mut b = Builder::default();
    let (p, q) = (b.point([0., -10.]), b.point([30., -10.]));
    b.line(p, q);
    let drawn = b
        .geometry()
        .with_support_edges(vec![(edge, SupportEdge::Line([30., -10.], [0., -10.]))]);
    let target = resolve(&drawn, &[], &query_px([12., -10.], 0., 3.), None).unwrap();
    assert!(matches!(target.key, TargetKey::OnCurve(_)), "{target:?}");
    // A face corner is a dwell-acquirable reference (CAD-423, approved by
    // the product owner 2026-11-27); a point on a face edge is not.
    assert!(TargetKey::Support(corner).acquirable());
    assert!(!TargetKey::Support(edge).acquirable());
}

#[test]
fn a_short_concentric_sketch_arc_does_not_hide_the_rest_of_a_face_circle() {
    use std::f64::consts::{FRAC_PI_2, TAU};
    // A quarter arc of radius 10 about the origin, counter-clockwise from +x.
    let mut b = Builder::default();
    let center = b.point([0., 0.]);
    let start = b.point([10., 0.]);
    let end = b.point([0., 10.]);
    b.curve(CurveKind::Arc {
        center,
        start,
        end,
        clockwise: false,
    });
    let hole = support_key(TopologyKind::Edge, 7);
    let same_span = support_key(TopologyKind::Edge, 8);
    let face_circle = SupportEdge::Arc {
        center: [0., 0.],
        radius: 10.,
        start: 0.,
        sweep: TAU,
    };
    let quarter = SupportEdge::Arc {
        center: [0., 0.],
        radius: 10.,
        start: 0.,
        sweep: FRAC_PI_2,
    };
    let geometry = b
        .geometry()
        .with_support_edges(vec![(hole, face_circle), (same_span, quarter)]);
    // Opposite the sketch arc, the face circle is still a snap target.
    let target = resolve(&geometry, &[], &query([-10.4, 0.]), None).unwrap();
    assert_eq!(target.key, TargetKey::Support(hole), "{target:?}");
    assert_eq!(target.kind, TargetKind::FaceEdge);
    assert!((target.position[0] + 10.).abs() < 1e-12);
    // The face arc spanning exactly the sketch arc is left to the sketch arc.
    let s = std::f64::consts::FRAC_1_SQRT_2 * 10.;
    let target = resolve(&geometry, &[], &query([s + 0.3, s + 0.3]), None).unwrap();
    assert_eq!(target.kind, TargetKind::OnCurve, "{target:?}");
    // A sketch circle covers every face arc on it.
    let mut b = Builder::default();
    let center = b.point([0., 0.]);
    b.curve(CurveKind::Circle {
        center,
        radius_mm: 10.,
    });
    let covered = b.geometry().with_support_edges(vec![(hole, face_circle)]);
    let target = resolve(&covered, &[], &query([-10.4, 0.]), None).unwrap();
    assert_eq!(target.kind, TargetKind::OnCurve, "{target:?}");
}

#[test]
fn a_face_corner_offered_beside_a_far_sketch_point_uses_a_relative_tolerance() {
    // Far from the origin an existing point within rounding of the corner
    // keeps it, so the corner is reused rather than projected twice.
    let mut b = Builder::default();
    let far = 5.0e4;
    let existing = b.point([far + 1e-6, far]);
    let geometry = b
        .geometry()
        .with_support(vec![(support_key(TopologyKind::Vertex, 0), [far, far])]);
    let target = resolve(&geometry, &[], &query([far, far]), None).unwrap();
    assert_eq!(target.anchor, Some(existing), "{target:?}");
    // A clearly separate corner is still offered.
    let geometry = b.geometry().with_support(vec![(
        support_key(TopologyKind::Vertex, 0),
        [far + 0.5, far],
    )]);
    let target = resolve(&geometry, &[], &query([far + 0.5, far]), None).unwrap();
    assert_eq!(target.kind, TargetKind::FaceCorner, "{target:?}");
}

/// The top face outline of a 30 x 20 box whose vertical edges are filleted
/// with radius 4, plus a through hole of radius 3 split at a seam into two
/// half arcs, as the face's edges in sketch coordinates.
fn rounded_face(r: f64) -> Vec<(TopologyKey, SupportEdge)> {
    use std::f64::consts::{FRAC_PI_2, PI};
    let (w, h) = (30., 20.);
    let edge = |ordinal| support_key(TopologyKind::Edge, ordinal);
    let arc = |center, start| SupportEdge::Arc {
        center,
        radius: r,
        start,
        sweep: FRAC_PI_2,
    };
    vec![
        (edge(0), SupportEdge::Line([r, 0.], [w - r, 0.])),
        (edge(1), SupportEdge::Line([w, r], [w, h - r])),
        (edge(2), SupportEdge::Line([w - r, h], [r, h])),
        (edge(3), SupportEdge::Line([0., h - r], [0., r])),
        (edge(4), arc([w - r, r], -FRAC_PI_2)),
        (edge(5), arc([w - r, h - r], 0.)),
        (edge(6), arc([r, h - r], FRAC_PI_2)),
        (edge(7), arc([r, r], PI)),
        (
            edge(8),
            SupportEdge::Arc {
                center: [15., 10.],
                radius: 3.,
                start: 0.,
                sweep: PI,
            },
        ),
        (
            edge(9),
            SupportEdge::Arc {
                center: [15., 10.],
                radius: 3.,
                start: PI,
                sweep: PI,
            },
        ),
    ]
}

#[test]
fn every_circular_face_edge_offers_its_exact_centre_as_a_face_centre_target() {
    let r = 4.;
    let geometry = SnapGeometry::default().with_support_edges(rounded_face(r));
    let corners = [[26., 4.], [26., 16.], [4., 16.], [4., 4.]];
    for (i, centre) in corners.into_iter().enumerate() {
        let target = resolve(&geometry, &[], &query_px(centre, 3., -2.), None).unwrap();
        assert_eq!(
            target.key,
            TargetKey::SupportCenter(support_key(TopologyKind::Edge, 4 + i as u32)),
            "{target:?}"
        );
        assert_eq!(target.kind, TargetKind::FaceCenter);
        assert_eq!(target.kind.label(), "Face centre");
        assert_eq!(target.kind.marker(), Marker::Center);
        assert_eq!(target.position, centre, "exactly the arc centre");
        assert_eq!(
            target.anchor, None,
            "the projected centre does not exist yet"
        );
    }
    // The hole's two half arcs share one centre target.
    let centres: Vec<_> = candidates(&geometry, &[], &query([0., 0.]))
        .into_iter()
        .filter(|c| c.target.kind == TargetKind::FaceCenter)
        .map(|c| c.target.position)
        .collect();
    assert_eq!(
        centres.len(),
        5,
        "four corner arcs and one hole: {centres:?}"
    );
    let hole = resolve(&geometry, &[], &query_px([15., 10.], -2., 2.), None).unwrap();
    assert_eq!(
        hole.key,
        TargetKey::SupportCenter(support_key(TopologyKind::Edge, 8))
    );
    assert_eq!(hole.position, [15., 10.]);
    // Straight edges have no centre: away from every centre there is none.
    let target = resolve(&geometry, &[], &query([15., 0.3]), None).unwrap();
    assert_eq!(target.kind, TargetKind::FaceEdge, "{target:?}");
}

#[test]
fn a_sketch_point_on_a_face_centre_keeps_it_and_a_centre_beats_a_nearby_edge() {
    let mut b = Builder::default();
    let existing = b.point([4., 4.]);
    let geometry = b.geometry().with_support_edges(rounded_face(4.));
    let target = resolve(&geometry, &[], &query([4., 4.]), None).unwrap();
    assert_eq!(target.anchor, Some(existing), "{target:?}");
    assert!(
        candidates(&geometry, &[], &query([0., 0.]))
            .iter()
            .all(|c| c.target.position != [4., 4.] || c.target.kind != TargetKind::FaceCenter),
        "the occupied centre is not offered twice"
    );
    // A small hole: its centre lies within the edge snap distance of the
    // circle, and the centre (a point class) wins.
    let small = support_key(TopologyKind::Edge, 1);
    let geometry = SnapGeometry::default().with_support_edges(vec![(
        small,
        SupportEdge::Arc {
            center: [0., 0.],
            radius: 1.,
            start: 0.,
            sweep: TAU,
        },
    )]);
    let target = resolve(&geometry, &[], &query_px([0., 0.], 2., 0.), None).unwrap();
    assert_eq!(target.key, TargetKey::SupportCenter(small), "{target:?}");
}

fn parallel(a: [f64; 2], b: [f64; 2]) -> bool {
    (a[0] * b[1] - a[1] * b[0]).abs() < 1e-12 * (1. + len(a) * len(b))
}
fn unit(d: [f64; 2]) -> [f64; 2] {
    let l = len(d);
    [d[0] / l, d[1] / l]
}

#[test]
fn a_triangle_corner_offers_exactly_x_y_and_its_two_extended_sides() {
    let mut b = Builder::default();
    let (points, _) = b.polygon(&[[0., 0.], [20., 5.], [7., 15.]]);
    let g = b.geometry();
    let mut assist = SnapAssist::default();
    acquire(&mut assist, &g, [0., 0.], 0.);
    let r = assist.references()[0].clone();
    assert_eq!(r.source, TargetKey::Point(points[0]));
    let guides = r.guides();
    assert_eq!(guides.len(), 4, "X, Y and two sides: {guides:?}");
    assert_eq!(guides[0], (Axis::Horizontal, [1., 0.]));
    assert_eq!(guides[1], (Axis::Vertical, [0., 1.]));
    // Each side continues past the corner, away from its other end.
    let expected = [unit([-20., -5.]), unit([-7., -15.])];
    for (axis, d) in &guides[2..] {
        assert!(matches!(axis, Axis::Edge(_)));
        let e = expected.iter().find(|e| parallel(**e, *d)).expect("a side");
        assert!(dot(*e, *d) > 0.999_999, "points past the corner: {d:?}");
    }
    // Moving along an extension snaps onto it exactly.
    let along = unit([-20., -5.]);
    let normal = [-along[1], along[0]];
    let pointer = [
        -30. * along[0] + 0.4 * normal[0],
        -30. * along[1] + 0.4 * normal[1],
    ];
    let target = assist.update(&g, &query(pointer), 2.).target.unwrap();
    assert!(
        matches!(target.kind, TargetKind::Guide(Axis::Edge(_))),
        "{target:?}"
    );
    assert_eq!(target.kind.label(), "Edge extension");
    assert_eq!(target.kind.marker(), Marker::Alignment);
    assert!(
        parallel(target.position, along) && dot(target.position, along) < 0.,
        "on the extension past the corner: {:?}",
        target.position
    );
    assert_eq!(target.guides[0].from, [0., 0.]);
    assert_eq!(
        placement(Some(&target), pointer, Some(5.)).0,
        target.position
    );
    // Beyond the guide tolerance (6 px = 1.5 mm here) it does not snap.
    let away = [pointer[0] + 2. * normal[0], pointer[1] + 2. * normal[1]];
    let target = assist.update(&g, &query(away), 3.).target;
    assert!(
        !target.is_some_and(|t| matches!(t.kind, TargetKind::Guide(Axis::Edge(_)))),
        "no extension snap 2.4 mm off it"
    );
}

/// A 30 x 20 rectangle with corner radius 4, rotated by 30 degrees about the
/// origin: four sides and four counter-clockwise corner arcs sharing their
/// tangent points. Returns the builder and the rotation.
fn rotated_rounded_rectangle() -> (Builder, impl Fn([f64; 2]) -> [f64; 2]) {
    let (w, h, r) = (30., 20., 4.);
    let (s, c) = (0.5f64, 0.75f64.sqrt());
    let rot = move |p: [f64; 2]| [c * p[0] - s * p[1], s * p[0] + c * p[1]];
    let mut b = Builder::default();
    let tangent = [
        [r, 0.],
        [w - r, 0.],
        [w, r],
        [w, h - r],
        [w - r, h],
        [r, h],
        [0., h - r],
        [0., r],
    ]
    .map(|p| b.point(rot(p)));
    let centers = [[w - r, r], [w - r, h - r], [r, h - r], [r, r]].map(|p| b.point(rot(p)));
    for k in 0..4 {
        b.line(tangent[2 * k], tangent[2 * k + 1]);
        b.curve(CurveKind::Arc {
            center: centers[k],
            start: tangent[2 * k + 1],
            end: tangent[(2 * k + 2) % 8],
            clockwise: false,
        });
    }
    (b, rot)
}

#[test]
fn both_ends_of_a_rounded_corner_extend_their_sides_and_cross_at_the_virtual_corner() {
    let (b, rot) = rotated_rounded_rectangle();
    let g = b.geometry();
    let mut assist = SnapAssist::default();
    let (a, c) = (rot([30., 16.]), rot([26., 20.]));
    acquire(&mut assist, &g, a, 0.);
    acquire(&mut assist, &g, c, 1.);
    let [ra, rc] = [&assist.references()[0], &assist.references()[1]];
    // The side and the tangent arc meet in one line: one extension each.
    assert_eq!(ra.directions.len(), 1, "{:?}", ra.directions);
    assert!(parallel(ra.directions[0], rot([0., 1.])));
    assert_eq!(rc.directions.len(), 1, "{:?}", rc.directions);
    assert!(parallel(rc.directions[0], rot([1., 0.])));
    let (sa, sc) = (ra.serial, rc.serial);
    let corner = rot([30., 20.]);
    let target = assist
        .update(&g, &query_px(corner, 1., -1.), 2.)
        .target
        .unwrap();
    assert_eq!(
        target.key,
        TargetKey::GuideCrossing {
            first: (sa, Axis::Edge(0)),
            second: (sc, Axis::Edge(0)),
        }
    );
    assert_eq!(target.kind, TargetKind::GuideIntersection);
    assert_eq!(target.kind.marker(), Marker::Construction);
    assert!(
        len(sub(target.position, corner)) < 1e-12,
        "the exact virtual corner: {:?} vs {corner:?}",
        target.position
    );
    assert_eq!(target.guides.len(), 2);
    assert_eq!(
        placement(Some(&target), [0.; 2], Some(10.)).0,
        target.position
    );
    // X/Y guides of one reference also cross the other's extension.
    let crossings = candidates(&g, assist.references(), &query(corner))
        .into_iter()
        .filter(|c| matches!(c.target.key, TargetKey::GuideCrossing { .. }))
        .count();
    assert_eq!(
        crossings, 5,
        "edge x edge plus edge x X/Y in both directions"
    );
}

#[test]
fn face_corners_and_face_centres_are_acquired_with_their_face_edge_guides() {
    use std::f64::consts::FRAC_PI_3;
    let corner = support_key(TopologyKind::Vertex, 0);
    let arc_end = support_key(TopologyKind::Vertex, 1);
    let arc = support_key(TopologyKind::Edge, 3);
    // The arc's end, at 60 degrees about its centre (40, 0).
    let end = [40. + 10. * FRAC_PI_3.cos(), 10. * FRAC_PI_3.sin()];
    let g = SnapGeometry::default()
        .with_support(vec![(corner, [0., 0.]), (arc_end, end)])
        .with_support_edges(vec![
            (
                support_key(TopologyKind::Edge, 0),
                SupportEdge::Line([20., 5.], [0., 0.]),
            ),
            (
                support_key(TopologyKind::Edge, 1),
                SupportEdge::Line([0., 0.], [7., 15.]),
            ),
            (
                arc,
                SupportEdge::Arc {
                    center: [40., 0.],
                    radius: 10.,
                    start: 0.,
                    sweep: FRAC_PI_3,
                },
            ),
        ]);
    let mut assist = SnapAssist::default();
    acquire(&mut assist, &g, [0., 0.], 0.);
    let r = &assist.references()[0];
    assert_eq!(r.source, TargetKey::Support(corner));
    assert_eq!(r.kind, TargetKind::FaceCorner);
    assert_eq!(r.guides().len(), 4, "X, Y and two face edges");
    // An arc's end extends along its tangent there.
    acquire(&mut assist, &g, end, 1.);
    let r = &assist.references()[1];
    assert_eq!(r.source, TargetKey::Support(arc_end));
    assert_eq!(r.directions.len(), 1);
    let tangent = [-FRAC_PI_3.sin(), FRAC_PI_3.cos()];
    assert!(
        dot(r.directions[0], tangent) > 1. - 1e-12,
        "{:?}",
        r.directions
    );
    // A face centre is acquired too; it has X and Y guides only.
    acquire(&mut assist, &g, [40., 0.], 2.);
    let r = assist.references().last().unwrap();
    assert_eq!(r.source, TargetKey::SupportCenter(arc));
    assert_eq!(r.kind, TargetKind::FaceCenter);
    assert_eq!(r.guides().len(), 2);
    // Leaving and dwelling again releases it.
    assist.update(&g, &query([80., 80.]), 3.5);
    let q = query_px([40., 0.], 1., 1.);
    let out = hover(&mut assist, &g, &q, &[4., 4. + DWELL_SECONDS]);
    assert!(
        matches!(out.event, Some(DwellEvent::Released(_))),
        "{out:?}"
    );
    // The centre took the oldest (corner) reference's place; the arc end
    // is held. Once projected, a sketch point holds that corner and the
    // reference stays on it; when nothing is there any more it is dropped.
    assert_eq!(assist.references().len(), 1);
    let mut b = Builder::default();
    b.point(end);
    assist.revalidate(&b.geometry());
    assert_eq!(assist.references()[0].position, end);
    assist.revalidate(&SnapGeometry::default());
    assert!(
        assist.references().is_empty(),
        "a vanished corner is dropped"
    );
}
