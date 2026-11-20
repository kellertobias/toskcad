use super::*;
use multicad_core::manufacturing::PlanarEdge;
use multicad_core::references::TopologyElement;
use multicad_rust_kernel::{
    normal_relief::Kind,
    spline_sweep::{Cubic, Sweep},
};
struct Profile {
    id: EntityId,
    sweep: Sweep,
    frame: PlaneFrame,
    walls: Vec<Vec<Option<TopologyElement>>>,
    caps: [Option<TopologyElement>; 2],
    source_enclosure: source_precision::Receipt,
}
fn cardinal(v: [f64; 3]) -> bool {
    v.iter().filter(|&&x| x.abs() == 1.).count() == 1
        && v.into_iter().all(|x| x == 0. || x.abs() == 1.)
}
fn exact_frame(f: PlaneFrame) -> bool {
    cardinal(f.u) && cardinal(f.v) && cardinal(f.normal())
}
fn selected_frame(stock: &Solid, face: u32) -> Result<Option<PlaneFrame>, GeometryError> {
    Ok(stock.retained_cap_source(face)?.map(|s| {
        let f = place_frame(s.placement, s.frame);
        PlaneFrame {
            origin: f.origin_mm,
            u: f.u,
            v: f.v,
        }
    }))
}
struct SourceWork<'a> {
    remaining: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl SourceWork<'_> {
    fn charge(&mut self, n: usize) -> Result<(), GeometryError> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(fail(
                GeometryErrorKind::Cancelled,
                "Spline cap Boolean was superseded",
            ));
        }
        self.remaining = self.remaining.checked_sub(n).ok_or_else(|| {
            fail(
                GeometryErrorKind::Limit,
                "Source ancestry exceeds the cap Boolean work budget",
            )
        })?;
        Ok(())
    }
}
struct SourceIndex<'a> {
    elements:
        std::collections::BTreeMap<multicad_core::references::TopologyKey, &'a TopologyElement>,
    faces: std::collections::BTreeMap<u32, &'a TopologyElement>,
    unique: std::collections::BTreeMap<u32, BTreeSet<Origin>>,
    by_origin: std::collections::BTreeMap<&'a Origin, BTreeSet<u32>>,
}
impl<'a> SourceIndex<'a> {
    fn new(
        snapshot: &'a TopologySnapshot,
        work: &mut SourceWork<'_>,
    ) -> Result<Self, GeometryError> {
        let mut elements = std::collections::BTreeMap::new();
        let mut faces = std::collections::BTreeMap::new();
        for element in &snapshot.elements {
            work.charge(1 + element.origins.len() + element.boundary.len() + element.frames.len())?;
            if element.origins.len() > multicad_core::references::MAX_GEOMETRY_ORIGINS
                || element.boundary.len() > multicad_core::references::MAX_GEOMETRY_ORIGINS
            {
                return Err(fail(
                    GeometryErrorKind::Limit,
                    "Source ancestry exceeds its element budget",
                ));
            }
            if elements.insert(element.key, element).is_some() {
                return Err(fail(
                    GeometryErrorKind::InvalidResult,
                    "Source topology repeats an element key",
                ));
            }
            if element.key.kind == TopologyKind::Face {
                faces.insert(element.key.ordinal, element);
            }
        }
        let common = faces
            .values()
            .next()
            .map(|first| {
                faces
                    .values()
                    .skip(1)
                    .fold(first.origins.clone(), |common, f| {
                        common.intersection(&f.origins).cloned().collect()
                    })
            })
            .unwrap_or_default();
        let mut unique = std::collections::BTreeMap::new();
        let mut by_origin: std::collections::BTreeMap<&Origin, BTreeSet<u32>> =
            std::collections::BTreeMap::new();
        for (&ordinal, face) in &faces {
            work.charge(1 + face.origins.len())?;
            let own = face.origins.difference(&common).cloned().collect();
            unique.insert(ordinal, own);
            for origin in face
                .origins
                .iter()
                .filter(|origin| !common.contains(*origin))
            {
                by_origin.entry(origin).or_default().insert(ordinal);
            }
        }
        Ok(Self {
            elements,
            faces,
            unique,
            by_origin,
        })
    }
    fn side(
        &self,
        key: multicad_core::references::TopologyKey,
        cap: u32,
        opposite: u32,
        work: &mut SourceWork<'_>,
    ) -> Result<Option<TopologyElement>, GeometryError> {
        let Some(edge) = self.elements.get(&key) else {
            return Ok(None);
        };
        let mut candidates = BTreeSet::new();
        for origin in &edge.boundary {
            work.charge(1)?;
            if let Some(faces) = self.by_origin.get(origin) {
                work.charge(faces.len())?;
                candidates.extend(faces);
            }
        }
        let mut matches = Vec::new();
        for ordinal in candidates {
            if ordinal == cap || ordinal == opposite {
                continue;
            }
            let own = &self.unique[&ordinal];
            work.charge(1 + own.len())?;
            if !own.is_empty() && own.is_subset(&edge.boundary) {
                matches.push(self.faces[&ordinal]);
            }
        }
        Ok((matches.len() == 1).then(|| matches[0].clone()))
    }
}
fn construction_error(e: multicad_rust_kernel::spline_sweep::ConstructionError) -> GeometryError {
    use multicad_rust_kernel::spline_sweep::ConstructionError as E;
    fail(
        match e {
            E::Cancelled => GeometryErrorKind::Cancelled,
            E::Limit => GeometryErrorKind::Limit,
            _ => GeometryErrorKind::Unsupported,
        },
        if e == E::Cancelled {
            "Spline cap Boolean was superseded".into()
        } else {
            format!("The retained profile certificate refused: {e:?}")
        },
    )
}
fn analytic(
    stock: &Solid,
    cap: u32,
    _feature: EntityId,
    policy: GeometryPolicy,
    work: &mut SourceWork<'_>,
) -> Result<Option<Profile>, GeometryError> {
    let Some(frame) = selected_frame(stock, cap)? else {
        return Ok(None);
    };
    let Some(source) = normal_relief_body::source::retained_source(stock, cap, frame)? else {
        return Ok(None);
    };
    let snapshot = stock.rust_topology().ok_or_else(|| {
        fail(
            GeometryErrorKind::Unsupported,
            "Missing retained Rust source",
        )
    })??;
    let index = SourceIndex::new(&snapshot, work)?;
    let mut loops = source
        .contours
        .contours
        .loops
        .iter()
        .zip(&source.contours.edge_keys)
        .collect::<Vec<_>>();
    loops.sort_by_key(|(r, _)| !r.outer);
    if loops.iter().filter(|(r, _)| r.outer).count() != 1 {
        return Ok(None);
    }
    let mut rings = vec![];
    let mut walls = vec![];
    for (i, (ring, keys)) in loops.into_iter().enumerate() {
        let mut edges = ring
            .edges
            .iter()
            .map(|e| {
                Edge::Analytic(match *e {
                    PlanarEdge::Line { a, b } => AnalyticEdge::Line { a, b },
                    PlanarEdge::Arc {
                        a,
                        b,
                        center,
                        sweep_rad,
                    } if a == b && sweep_rad.abs() == std::f64::consts::TAU => {
                        AnalyticEdge::Circle {
                            center,
                            radius_mm: (a[0] - center[0]).hypot(a[1] - center[1]),
                            clockwise: sweep_rad < 0.,
                        }
                    }
                    PlanarEdge::Arc {
                        a,
                        b,
                        center,
                        sweep_rad,
                    } => AnalyticEdge::Arc {
                        a,
                        b,
                        center,
                        clockwise: sweep_rad < 0.,
                    },
                })
            })
            .collect::<Vec<_>>();
        let mut labels = keys
            .iter()
            .map(|&key| index.side(key, cap, source.opposite, work))
            .collect::<Result<Vec<_>, _>>()?;
        let area = multicad_rust_kernel::normal_sweep::signed_area(
            &edges
                .iter()
                .filter_map(|e| {
                    if let Edge::Analytic(a) = e {
                        Some(*a)
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>(),
        );
        if (area > 0.) != (i == 0) {
            edges = edges.iter().rev().map(Edge::reversed).collect();
            labels.reverse();
        }
        rings.push(edges);
        walls.push(labels);
    }
    let outline = rings.remove(0);
    let sweep = Sweep::new_rings_with_interrupt(
        outline,
        rings,
        source.depth_mm,
        policy.linear_tolerance_mm,
        work.interrupt,
    )
    .map_err(construction_error)?;
    let caps =
        [source.opposite, cap].map(|ordinal| index.faces.get(&ordinal).map(|face| (*face).clone()));
    Ok(Some(Profile {
        id: snapshot.owner,
        sweep,
        frame: source.base_frame,
        walls,
        caps,
        source_enclosure: source_precision::Receipt::default(),
    }))
}
fn tool(
    source: &Solid,
    feature: EntityId,
    policy: GeometryPolicy,
    work: &mut SourceWork<'_>,
) -> Result<Option<Profile>, GeometryError> {
    if let Some((body, frame)) = source.spline_sweep_body() {
        let snapshot = source.rust_topology().ok_or_else(|| {
            fail(
                GeometryErrorKind::Unsupported,
                "Missing retained Rust source",
            )
        })??;
        let mut offset = 0;
        let all = SourceIndex::new(&snapshot, work)?;
        let walls = body
            .sweep
            .rings()
            .map(|ring| {
                ring.iter()
                    .map(|_| {
                        let i = offset;
                        offset += 1;
                        all.faces.get(&(i + 2)).map(|f| (*f).clone())
                    })
                    .collect()
            })
            .collect();
        let caps = [0, 1].map(|ordinal| all.faces.get(&ordinal).map(|f| (*f).clone()));
        return Ok(Some(Profile {
            id: snapshot.owner,
            sweep: body.sweep.clone(),
            frame,
            walls,
            caps,
            source_enclosure: {
                let (d, a) = body.source_enclosure();
                source_precision::Receipt {
                    distance_mm: d,
                    angle_rad: a,
                }
            },
        }));
    }
    analytic(source, 1, feature, policy, work)
}
fn add(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let p = s - a;
    (s, ((a - (s - p)) + (b - p)).abs())
}
fn mapped(
    edge: &Edge,
    source: PlaneFrame,
    target: PlaneFrame,
    tol: f64,
    receipt: &mut source_precision::Receipt,
) -> Result<Edge, GeometryError> {
    let mut max_error = 0_f64;
    let mut map = |p: [f64; 2]| -> Result<[f64; 2], GeometryError> {
        let mut error = 0.;
        let mut world = source.origin;
        for i in 0..3 {
            for (axis, value) in [source.u, source.v].into_iter().zip(p) {
                let (s, e) = add(world[i], axis[i] * value);
                world[i] = s;
                error = source_precision::add_up(error, e);
            }
        }
        let mut result = [0.; 2];
        for (i, axis) in [target.u, target.v].into_iter().enumerate() {
            let at = axis.iter().position(|x| *x != 0.).expect("cardinal source");
            let (s, e) = add(world[at], -target.origin[at]);
            result[i] = axis[at] * s;
            error = source_precision::add_up(error, e);
        }
        if !error.is_finite() || error > tol / 64. {
            return Err(fail(
                GeometryErrorKind::Unsupported,
                "The source UV translation exceeds its finite precision budget",
            ));
        }
        max_error = max_error.max(error);
        Ok(result)
    };
    let reflected = dot(source.u, target.u) * dot(source.v, target.v)
        - dot(source.u, target.v) * dot(source.v, target.u)
        < 0.;
    let result = match edge {
        Edge::Analytic(AnalyticEdge::Line { a, b }) => Edge::Analytic(AnalyticEdge::Line {
            a: map(*a)?,
            b: map(*b)?,
        }),
        Edge::Analytic(AnalyticEdge::Arc {
            a,
            b,
            center,
            clockwise,
        }) => Edge::Analytic(AnalyticEdge::Arc {
            a: map(*a)?,
            b: map(*b)?,
            center: map(*center)?,
            clockwise: *clockwise != reflected,
        }),
        Edge::Analytic(AnalyticEdge::Circle {
            center,
            radius_mm,
            clockwise,
        }) => Edge::Analytic(AnalyticEdge::Circle {
            center: map(*center)?,
            radius_mm: *radius_mm,
            clockwise: *clockwise != reflected,
        }),
        Edge::Rational(_) => {
            return Err(fail(
                GeometryErrorKind::Unsupported,
                "A cap relief cannot yet map a rational spline source",
            ));
        }
        Edge::Cubic(c) => Edge::Cubic(
            Cubic::new(
                c.poles()
                    .iter()
                    .map(|&p| map(p))
                    .collect::<Result<Vec<_>, _>>()?,
                c.knots().to_vec(),
                c.reversed(),
            )
            .map_err(|_| {
                fail(
                    GeometryErrorKind::Unsupported,
                    "The mapped spline source is unqualified",
                )
            })?,
        ),
    };
    let matrix = [target.u, target.v].map(|axis| [dot(axis, source.u), dot(axis, source.v)]);
    let measured = source_precision::certify(edge, &result, matrix, max_error)?;
    receipt.distance_mm = receipt.distance_mm.max(measured.distance_mm);
    receipt.angle_rad = receipt.angle_rad.max(measured.angle_rad);
    Ok(result)
}
fn local_labels(
    original: Option<&TopologyElement>,
    feature: EntityId,
    role: String,
    frame: PlaneFrame,
    policy: GeometryPolicy,
) -> Result<Labels, GeometryError> {
    let mut origins = original.map(|e| e.origins.clone()).unwrap_or_default();
    origins.insert(Origin {
        source: feature,
        role,
    });
    let anchors = original
        .into_iter()
        .flat_map(|e| &e.frames)
        .map(|a| {
            Ok(AnchoredFrame {
                origin: a.origin.clone(),
                frame: revolution_body::local_attachment(frame, a.frame, policy)?,
            })
        })
        .collect::<Result<Vec<_>, GeometryError>>()?;
    Ok(Labels { origins, anchors })
}
pub(super) fn build(
    operands: &[&Solid],
    operation: BooleanOperation,
    feature: EntityId,
    policy: GeometryPolicy,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Option<Result<Solid, GeometryError>> {
    if operands.len() < 2
        || !operands[1..]
            .iter()
            .any(|s| s.spline_sweep_body().is_some())
        || operation == BooleanOperation::Intersect
    {
        return None;
    }
    Some((|| {
        policy.validate()?;
        if feature.is_nil()
            || operands.iter().any(|s| s.policy != policy)
            || policy.boolean_fuzzy_mm != 0.
        {
            return Err(fail(
                GeometryErrorKind::Input,
                "Exact cap Boolean operands need one policy and feature identity",
            ));
        }
        if operands.len() > multicad_rust_kernel::spline_cap_relief::MAX_MOTIFS + 1 {
            return Err(fail(
                GeometryErrorKind::Limit,
                "Too many retained cap profiles",
            ));
        }
        let kind = if operation == BooleanOperation::Union {
            Kind::Raised
        } else {
            Kind::Recessed
        };
        let checkpoint = || {
            if interrupt.is_some_and(|f| f()) {
                Err(fail(
                    GeometryErrorKind::Cancelled,
                    "Spline cap Boolean was superseded",
                ))
            } else {
                Ok(())
            }
        };
        checkpoint()?;
        let mut work = SourceWork {
            remaining: 4_000_000,
            interrupt,
        };
        let mut edge_count = 0_usize;
        let mut tools = operands[1..]
            .iter()
            .map(|s| {
                checkpoint()?;
                tool(s, feature, policy, &mut work)?
                    .ok_or_else(|| {
                        fail(
                            GeometryErrorKind::Unsupported,
                            "A Boolean tool has no retained straight profile",
                        )
                    })
                    .and_then(|profile| {
                        edge_count += profile.sweep.edges().count();
                        if edge_count > multicad_rust_kernel::spline_cap_relief::MAX_EDGES {
                            Err(fail(
                                GeometryErrorKind::Limit,
                                "Too many retained profile edges",
                            ))
                        } else {
                            Ok(profile)
                        }
                    })
            })
            .collect::<Result<Vec<_>, _>>()?;
        let first = &tools[0];
        if !exact_frame(first.frame) {
            return Err(fail(
                GeometryErrorKind::Unsupported,
                "Spline cap Booleans currently require signed-cardinal source charts",
            ));
        }
        let expected = first
            .frame
            .normal()
            .map(|n| if kind == Kind::Raised { n } else { -n });
        let mut stock = None;
        for face in operands[0].faces()? {
            checkpoint()?;
            if let Some((origin, normal)) = face.plane
                && normal == expected
            {
                let axis = normal
                    .iter()
                    .position(|x| x.abs() == 1.)
                    .expect("signed cardinal");
                if origin[axis] == first.frame.origin[axis] {
                    stock = analytic(operands[0], face.index, feature, policy, &mut work)?;
                    if stock.is_some() {
                        break;
                    }
                }
            }
        }
        let stock = stock.ok_or_else(|| {
            fail(
                GeometryErrorKind::Unsupported,
                "The tools do not start exactly on a qualified stock cap",
            )
        })?;
        let frame = stock.frame;
        if !exact_frame(frame) {
            return Err(fail(
                GeometryErrorKind::Unsupported,
                "The stock source chart is not signed-cardinal",
            ));
        }
        let axis = frame
            .normal()
            .iter()
            .position(|x| x.abs() == 1.)
            .expect("signed cardinal");
        let entry = position(frame, [0., 0., stock.sweep.height_mm()]);
        let mut motifs = Vec::new();
        let mut source_enclosure = source_precision::Receipt::default();
        for t in &mut tools {
            checkpoint()?;
            let mut mapped_enclosure = source_precision::Receipt::default();
            if !exact_frame(t.frame)
                || t.frame.origin[axis] != entry[axis]
                || t.frame.normal() != first_normal(kind, frame.normal())
            {
                return Err(fail(
                    GeometryErrorKind::Unsupported,
                    "Every tool must have the exact source plane and axial direction",
                ));
            }
            let reflection = dot(t.frame.normal(), frame.normal()) < 0.;
            let mut rings = Vec::new();
            for (r, ring) in t.sweep.rings().enumerate() {
                let mut edges = ring
                    .iter()
                    .map(|e| {
                        work.charge(match e {
                            Edge::Cubic(c) => (c.knots().len() - 1) * 64,
                            Edge::Rational(_) => 64,
                            Edge::Analytic(_) => 32,
                        })?;
                        mapped(
                            e,
                            t.frame,
                            frame,
                            policy.linear_tolerance_mm,
                            &mut mapped_enclosure,
                        )
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                if reflection {
                    edges = edges.iter().rev().map(Edge::reversed).collect();
                    t.walls[r].reverse();
                }
                rings.push(edges);
            }
            source_enclosure.distance_mm =
                source_enclosure.distance_mm.max(source_precision::add_up(
                    t.source_enclosure.distance_mm,
                    mapped_enclosure.distance_mm,
                ));
            source_enclosure.angle_rad = source_enclosure.angle_rad.max(source_precision::add_up(
                t.source_enclosure.angle_rad,
                mapped_enclosure.angle_rad,
            ));
            let outline = rings.remove(0);
            motifs.push(
                Sweep::new_rings_with_interrupt(
                    outline,
                    rings,
                    t.sweep.height_mm(),
                    policy.linear_tolerance_mm,
                    interrupt,
                )
                .map_err(construction_error)?,
            );
        }
        let geometry = Relief::new_with_interrupt(
            stock.sweep.clone(),
            motifs,
            kind,
            policy.linear_tolerance_mm,
            policy.minimum_feature_mm,
            interrupt,
        )
        .map_err(|e| {
            fail(
                match e {
                    multicad_rust_kernel::spline_sweep::ConstructionError::Cancelled => {
                        GeometryErrorKind::Cancelled
                    }
                    multicad_rust_kernel::spline_sweep::ConstructionError::Limit => {
                        GeometryErrorKind::Limit
                    }
                    _ => GeometryErrorKind::Unsupported,
                },
                format!(
                    "Retained spline cap material could not be certified: {e:?}{}",
                    if e == multicad_rust_kernel::spline_sweep::ConstructionError::Cancelled {
                        " superseded"
                    } else {
                        ""
                    }
                ),
            )
        })?;
        let labels = geometry
            .boundary()
            .faces
            .iter()
            .map(|face| {
                let original = match face.role {
                    FaceRole::BaseCap { top } => stock.caps[usize::from(top)].as_ref(),
                    FaceRole::BaseWall { ring, edge } => stock.walls[ring][edge].as_ref(),
                    FaceRole::MotifCap { motif } => tools[motif].caps[1].as_ref(),
                    FaceRole::MotifWall { motif, ring, edge } => {
                        tools[motif].walls[ring][edge].as_ref()
                    }
                    FaceRole::CounterCap { .. } => None,
                };
                let role = match face.role {
                    FaceRole::BaseCap { top } => format!("stock/cap/{top}"),
                    FaceRole::BaseWall { ring, edge } => format!(
                        "stock/wall/{}",
                        ContentHash::of(
                            format!(
                                "{:?}",
                                stock.sweep.rings().nth(ring).expect("source ring")[edge]
                            )
                            .as_bytes()
                        )
                    ),
                    FaceRole::MotifCap { motif } => {
                        format!("motif/{}/cap/{kind:?}", tools[motif].id)
                    }
                    FaceRole::CounterCap { motif, hole } => format!(
                        "motif/{}/counter/{}",
                        tools[motif].id,
                        ContentHash::of(
                            format!(
                                "{:?}",
                                tools[motif].walls[hole + 1]
                                    .iter()
                                    .flatten()
                                    .map(|f| &f.origins)
                                    .collect::<Vec<_>>()
                            )
                            .as_bytes()
                        )
                    ),
                    FaceRole::MotifWall { motif, ring, edge } => format!(
                        "motif/{}/wall/{kind:?}/{}",
                        tools[motif].id,
                        ContentHash::of(
                            format!(
                                "{:?}",
                                geometry.motifs()[motif]
                                    .rings()
                                    .nth(ring)
                                    .expect("source ring")[edge]
                            )
                            .as_bytes()
                        )
                    ),
                };
                let role = if original.is_some_and(|face| !face.origins.is_empty()) {
                    match face.role {
                        FaceRole::BaseWall { .. } => "stock/retained-surface".into(),
                        FaceRole::MotifWall { motif, .. } => {
                            format!("motif/{}/retained-surface/{kind:?}", tools[motif].id)
                        }
                        _ => role,
                    }
                } else {
                    role
                };
                local_labels(original, feature, role, frame, policy)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let body = SplineCapReliefBody::new(geometry, feature, labels, policy, source_enclosure)?;
        frame.checked_metric_bounds(body.bounds(PlaneFrame::default()), policy)?;
        let bounds = body.bounds(frame);
        policy.point(bounds.min)?;
        policy.point(bounds.max)?;
        Ok(Solid {
            body: SolidBody::SplineCapRelief(Rc::new(body), frame),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        })
    })())
}
fn first_normal(kind: Kind, n: [f64; 3]) -> [f64; 3] {
    n.map(|v| if kind == Kind::Raised { v } else { -v })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn indexed_ancestry_budget_cancellation_and_recovery_use_real_source_topology() {
        let stock =
            Solid::box_named([20., 20., 5.], EntityId::new(), GeometryPolicy::default()).unwrap();
        let snapshot = stock.rust_topology().unwrap().unwrap();
        let before = snapshot.clone();
        let mut work = SourceWork {
            remaining: 1,
            interrupt: None,
        };
        assert_eq!(
            SourceIndex::new(&snapshot, &mut work).err().unwrap().kind,
            GeometryErrorKind::Limit
        );
        let cancel = || true;
        let mut work = SourceWork {
            remaining: 4_000_000,
            interrupt: Some(&cancel),
        };
        assert_eq!(
            SourceIndex::new(&snapshot, &mut work).err().unwrap().kind,
            GeometryErrorKind::Cancelled
        );
        let mut work = SourceWork {
            remaining: 4_000_000,
            interrupt: None,
        };
        let index = SourceIndex::new(&snapshot, &mut work).unwrap();
        let selected = snapshot
            .elements
            .iter()
            .filter(|e| e.key.kind == TopologyKind::Edge)
            .find_map(|e| index.side(e.key, 0, 1, &mut work).unwrap());
        assert!(selected.is_some());
        assert_eq!(snapshot, before);
    }
}
