//! Owned local coupled insert candidate. No independent import/mass authority.
use super::*;
mod mapped;
use crate::step::{Coordinates, Entities, Entity, StepError, StepFile, Writer};
pub use mapped::write_mapped_candidate;
impl From<StepError> for Error {
    fn from(e: StepError) -> Self {
        match e {
            StepError::Limit => Self::Limit,
            StepError::InvalidGeometry | StepError::InvalidPlacement => Self::Invalid,
            StepError::Unsupported => Self::Unresolved,
        }
    }
}
/// Emit the actual local coupled candidate after all local proof gates. This
/// source-known output does not open Native export or foreign-import admission.
pub fn write_local_candidate(
    topology: &TopologyPlan<'_, '_, '_>,
    maximum_work: usize,
    maximum_bytes: usize,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<String, Error> {
    let mut writer = Writer::new();
    writer.append_owned_insert_local(topology, maximum_work, maximum_bytes, interrupt)?;
    writer.finish().map_err(Error::from)
}
impl Writer {
    /// Append the qualified local coupled BRep through the same atomic
    /// compound transaction as other members. Geometry remains at its actual
    /// local datum; this method does not imply placed Native admission.
    ///
    /// `maximum_bytes` covers the retained insert proof arenas, complete output
    /// String, and encoding scratch. Pre-existing generic Writer point-cache
    /// and solid-index metadata remain governed by its separate compound limits;
    /// this is not a whole-process allocator-peak claim. Failed appends restore
    /// published records/IDs while retaining allocation capacity, as Writer does.
    pub fn append_owned_insert_local(
        &mut self,
        topology: &TopologyPlan<'_, '_, '_>,
        maximum_work: usize,
        maximum_bytes: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<(), Error> {
        self.append_typed(|file| {
            encode_local(file, topology, maximum_work, maximum_bytes, interrupt, 0)
                .map(|(solid, _, _)| solid)
        })
    }
}
fn encode_local<'cancel>(
    file: &mut StepFile,
    topology: &TopologyPlan<'_, '_, '_>,
    maximum_work: usize,
    maximum_bytes: usize,
    interrupt: Option<&'cancel dyn Fn() -> bool>,
    extra_output: usize,
) -> Result<(Entity, Work<'cancel>, Memory), Error> {
    let mut w = Work::new(maximum_work, interrupt)?;
    w.charge(topology.work_used())?;
    if maximum_bytes == 0 || maximum_bytes > 64 * 1024 * 1024 {
        return Err(Error::Invalid);
    }
    let surface = topology.surface();
    let mut memory = Memory::new(topology.retained_bytes(), maximum_bytes)?;
    let prior_capacity = file.data.capacity();
    // The entire live String is counted once, including earlier compound members.
    let fixed = 4096usize.checked_add(prior_capacity).ok_or(Error::Limit)?;
    memory.add(fixed)?;
    let projection_bytes = maximum_bytes
        .checked_sub(topology.retained_bytes() - surface.retained_bytes() + fixed)
        .ok_or(Error::Limit)?;
    let p = surface.projection_bounds(
        w.remaining() + surface.work_used(),
        projection_bytes,
        interrupt,
    )?;
    w.charge(p.work_used - surface.work_used())?;
    if p.component_radius_mm[0][0] <= p.component_radius_mm[1][1] {
        return Err(Error::Unresolved);
    }
    let caps = topology.cap_regularity_bounds(
        w.remaining() + topology.work_used(),
        maximum_bytes - fixed,
        interrupt,
    )?;
    w.charge(caps.work_used - topology.work_used())?;
    let controls = surface
        .rows()
        .iter()
        .map(|r| r.poles().len())
        .sum::<usize>();
    let estimate = output_estimate(controls, topology.edges().len(), topology.faces().len())?
        .saturating_add(extra_output);
    let base = memory.used;
    memory.add(estimate)?;
    memory.used = base;
    file.data
        .try_reserve_exact(estimate)
        .map_err(|_| Error::Limit)?;
    memory.add(file.data.capacity() - prior_capacity)?;
    let solid = emit(file, topology, &mut w, &mut memory)?;
    w.charge(1)?;
    Ok((solid, w, memory))
}

fn point(file: &mut StepFile, p: [f64; 3], w: &mut Work<'_>) -> Result<Entity, Error> {
    w.charge(128)?;
    Ok(file.entity(format_args!("CARTESIAN_POINT('',{})", Coordinates(p))))
}
fn spline(
    file: &mut StepFile,
    degree: usize,
    poles: &[Entity],
    knots: &[f64],
    w: &mut Work<'_>,
) -> Result<Entity, Error> {
    w.charge(512 + 24 * poles.len() + 64 * knots.len())?;
    let mult = Multiplicities(knots.len(), degree);
    let values = KnotValues(knots);
    Ok(file.entity(format_args!("B_SPLINE_CURVE_WITH_KNOTS('',{degree},{},.UNSPECIFIED.,.F.,.F.,({mult}),({values}),.UNSPECIFIED.)",Entities(poles))))
}
#[allow(clippy::too_many_arguments)]
fn pcurve(
    file: &mut StepFile,
    support: Entity,
    degree: usize,
    points: &[[f64; 2]],
    knots: &[f64],
    context: Entity,
    w: &mut Work<'_>,
    memory: &mut Memory,
) -> Result<Entity, Error> {
    let mut poles = memory.allocate(points.len(), w)?;
    for &[u, v] in points {
        w.charge(128)?;
        poles.push(file.entity(format_args!("CARTESIAN_POINT('',({u:.16E},{v:.16E}))")));
    }
    let curve = spline(file, degree, &poles, knots, w)?;
    memory.release(&poles);
    drop(poles);
    let rep = file.entity(format_args!(
        "DEFINITIONAL_REPRESENTATION('',({curve}),{context})"
    ));
    Ok(file.entity(format_args!("PCURVE('',{support},{rep})")))
}
fn emit(
    file: &mut StepFile,
    t: &TopologyPlan<'_, '_, '_>,
    w: &mut Work<'_>,
    memory: &mut Memory,
) -> Result<Entity, Error> {
    let s = t.surface();
    let mut xyz = memory.allocate(s.vertices().len(), w)?;
    let mut vertices = memory.allocate(s.vertices().len(), w)?;
    for v in s.vertices() {
        let p = point(file, v.point_mm(), w)?;
        xyz.push(p);
        vertices.push(file.entity(format_args!("VERTEX_POINT('',{p})")));
    }
    let mut row_points = memory.allocate(s.rows().len(), w)?;
    let mut curves = memory.allocate(t.edges().len(), w)?;
    for r in s.rows() {
        let ends = r.endpoints();
        let mut poles = memory.allocate(r.poles().len(), w)?;
        for (i, &p) in r.poles().iter().enumerate() {
            poles.push(if i == 0 {
                xyz[ends[0]]
            } else if i + 1 == r.poles().len() {
                xyz[ends[1]]
            } else {
                point(file, p, w)?
            });
        }
        curves.push(spline(file, 11, &poles, r.knots(), w)?);
        row_points.push(poles);
    }
    for e in &t.edges()[s.rows().len()..] {
        curves.push(spline(file, 1, &e.ends.map(|i| xyz[i]), &[0., 1.], w)?);
    }
    let mut supports = memory.allocate(t.faces().len(), w)?;
    for f in t.faces() {
        w.charge(512)?;
        supports.push(match f.geometry{
            FaceGeometry::Tensor(i)=>{let patch=&s.patches()[i];let [a,b]=patch.rows.map(|r|&row_points[r]);let knots=s.rows()[patch.rows[0]].knots();
                let mult=Multiplicities(knots.len(),11);
                let values=KnotValues(knots);
                file.entity(format_args!("B_SPLINE_SURFACE_WITH_KNOTS('',1,11,({},{}),.UNSPECIFIED.,.F.,.F.,.F.,(2,2),({mult}),(0.,1.),({values}),.UNSPECIFIED.)",Entities(a),Entities(b)))},
            FaceGeometry::Cap(end)=>{let h=if end==0{0.}else{s.domain().boundary().source().full_length_mm()};let p=point(file,[0.,0.,h],w)?;let axis=file.axis(p,[0.,0.,1.],[1.,0.,0.]);file.entity(format_args!("PLANE('',{axis})"))},
        });
    }
    let context=file.entity(format_args!("(GEOMETRIC_REPRESENTATION_CONTEXT(2)PARAMETRIC_REPRESENTATION_CONTEXT()REPRESENTATION_CONTEXT('2D SPACE',''))"));
    let mut edge_pcurves = memory.allocate(t.edges().len(), w)?;
    for _ in t.edges() {
        edge_pcurves.push(memory.allocate(2, w)?);
    }
    for (face, f) in t.faces().iter().enumerate() {
        for wire in &f.loops {
            for c in wire {
                let e = &t.edges()[c.edge];
                let support = supports[face];
                let pc = match f.geometry {
                    FaceGeometry::Cap(_) => {
                        let r = cap_row(e, s)?;
                        let mut uv = memory.allocate(r.poles().len(), w)?;
                        uv.extend(r.poles().iter().map(|p| [p[0], p[1]]));
                        let pc = pcurve(file, support, 11, &uv, r.knots(), context, w, memory)?;
                        memory.release(&uv);
                        drop(uv);
                        pc
                    }
                    FaceGeometry::Tensor(i) => {
                        let p = &s.patches()[i];
                        let [a, b] = p.rows.map(|r| &s.rows()[r]);
                        let uv = tensor_uv(p.rows, [a.endpoints(), b.endpoints()], e)?;
                        pcurve(file, support, 1, &uv, &[0., 1.], context, w, memory)?
                    }
                };
                edge_pcurves[c.edge].push(pc);
            }
        }
    }
    let mut edges = memory.allocate(t.edges().len(), w)?;
    for (i, e) in t.edges().iter().enumerate() {
        w.charge(512)?;
        let curve = file.entity(format_args!(
            "SURFACE_CURVE('',{},{},.CURVE_3D.)",
            curves[i],
            Entities(&edge_pcurves[i])
        ));
        let [a, b] = e.ends.map(|j| vertices[j]);
        edges.push(file.entity(format_args!("EDGE_CURVE('',{a},{b},{curve},.T.)")));
    }
    let mut faces = memory.allocate(t.faces().len(), w)?;
    for (i, f) in t.faces().iter().enumerate() {
        let mut bounds = memory.allocate(f.loops.len(), w)?;
        for (j, wire) in f.loops.iter().enumerate() {
            let mut oriented = memory.allocate(wire.len(), w)?;
            for c in wire {
                w.charge(64)?;
                oriented.push(file.oriented_edge(edges[c.edge], c.forward));
            }
            bounds.push(file.bound(&oriented, j == 0));
            memory.release(&oriented);
            drop(oriented);
        }
        faces.push(file.face_on(&bounds, supports[i], f.same_sense));
        memory.release(&bounds);
        drop(bounds);
    }
    w.charge(1)?;
    Ok(file.solid(&faces))
}

fn output_estimate(controls: usize, edges: usize, faces: usize) -> Result<usize, Error> {
    controls
        .checked_mul(1024)
        .and_then(|n| n.checked_add(edges.checked_mul(4096)?))
        .and_then(|n| n.checked_add(faces.checked_mul(4096)?))
        .ok_or(Error::Limit)
}
fn cap_row<'s>(
    edge: &ShellEdge,
    surface: &'s SurfacePlan<'_, '_>,
) -> Result<&'s SurfaceRow, Error> {
    let row = edge.geometry.row_index().ok_or(Error::Unresolved)?;
    Ok(&surface.rows()[row])
}
fn tensor_uv(
    rows: [usize; 2],
    ends: [[usize; 2]; 2],
    edge: &ShellEdge,
) -> Result<[[f64; 2]; 2], Error> {
    match edge.geometry {
        EdgeGeometry::Row(row) => {
            let u = if row == rows[0] {
                0.
            } else if row == rows[1] {
                1.
            } else {
                return Err(Error::Unresolved);
            };
            Ok([[u, 0.], [u, 1.]])
        }
        EdgeGeometry::Ruling => {
            let at = |vertex| -> Result<[f64; 2], Error> {
                for (v, (&a, &b)) in ends[0].iter().zip(&ends[1]).enumerate() {
                    if vertex == a {
                        return Ok([0., v as f64]);
                    }
                    if vertex == b {
                        return Ok([1., v as f64]);
                    }
                }
                Err(Error::Unresolved)
            };
            Ok([at(edge.ends[0])?, at(edge.ends[1])?])
        }
    }
}
struct KnotValues<'a>(&'a [f64]);
impl std::fmt::Display for KnotValues<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (i, v) in self.0.iter().enumerate() {
            if i > 0 {
                f.write_str(",")?;
            }
            write!(f, "{v:.16E}")?;
        }
        Ok(())
    }
}
struct Multiplicities(usize, usize);
impl std::fmt::Display for Multiplicities {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in 0..self.0 {
            if i > 0 {
                f.write_str(",")?;
            }
            write!(
                f,
                "{}",
                if i == 0 || i + 1 == self.0 {
                    self.1 + 1
                } else {
                    self.1
                }
            )?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn source() -> Boundary {
        super::super::tests::metric_graph()
    }
    #[test]
    fn actual_writer_limits_cancellation_and_error_mapping_never_publish_partial_files() {
        let body = source();
        let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
        let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
        let t = TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, None).unwrap();
        for bytes in [0, 64 * 1024 * 1024 + 1] {
            assert!(matches!(
                write_local_candidate(&t, 16_000_000, bytes, None),
                Err(Error::Invalid)
            ));
        }
        for bytes in [
            t.retained_bytes(),
            t.retained_bytes() + 4096,
            // Parent arenas/output fit; the actual projection child's proof
            // workspace is the next allocation and must refuse atomically.
            t.retained_bytes() + 4096 + StepFile::new().data.capacity() + 1,
        ] {
            assert!(matches!(
                write_local_candidate(&t, 16_000_000, bytes, None),
                Err(Error::Limit)
            ));
        }
        assert!(matches!(
            write_local_candidate(&t, 1, 64 * 1024 * 1024, None),
            Err(Error::Limit)
        ));
        assert!(matches!(
            write_local_candidate(&t, 16_000_000, 64 * 1024 * 1024, Some(&|| true)),
            Err(Error::Cancelled)
        ));
        // Projection fits, but its earned work receipt leaves insufficient
        // budget for the following independent cap/Jordan proof. The parent
        // must propagate that actual child refusal before encoding anything.
        let projection = s
            .projection_bounds(16_000_000, 64 * 1024 * 1024, None)
            .unwrap();
        let budget = t.work_used() + projection.work_used - s.work_used() + 10;
        assert!(matches!(
            write_local_candidate(&t, budget, 64 * 1024 * 1024, None),
            Err(Error::Limit)
        ));
        // Both geometric children fit, then the real encoder exhausts its
        // remaining budget AFTER writing part of the new candidate transaction.
        let caps = t
            .cap_regularity_bounds(16_000_000, 64 * 1024 * 1024, None)
            .unwrap();
        let budget = projection.work_used - s.work_used() + caps.work_used + 2000;
        assert!(matches!(
            write_local_candidate(&t, budget, 64 * 1024 * 1024, None),
            Err(Error::Limit)
        ));
        assert_eq!(Error::from(StepError::InvalidGeometry), Error::Invalid);
        assert_eq!(Error::from(StepError::InvalidPlacement), Error::Invalid);
        assert_eq!(Error::from(StepError::Unsupported), Error::Unresolved);
        assert_eq!(Error::from(StepError::Limit), Error::Limit);
    }
    #[test]
    fn actual_partial_insert_append_rolls_back_to_previous_real_sphere() {
        let body = source();
        let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
        let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
        let t = TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, None).unwrap();
        let seed = || {
            let mut writer = Writer::new();
            writer
                .append_sphere(
                    crate::sphere::Sphere::new(1.).unwrap(),
                    [0.; 3],
                    [1., 0., 0.],
                    [0., 1., 0.],
                )
                .unwrap();
            writer
        };
        let expected = seed().finish().unwrap();
        for cancel in [false, true] {
            let mut writer = seed();
            let polls = std::cell::Cell::new(0);
            let stop = || {
                polls.set(polls.get() + 1);
                cancel && polls.get() == 20
            };
            let mut w = Work::new(if cancel { 16_000_000 } else { 2000 }, Some(&stop)).unwrap();
            let mut memory = Memory::new(t.retained_bytes(), 64 * 1024 * 1024).unwrap();
            let failed = writer.append_typed(|file| emit(file, &t, &mut w, &mut memory));
            assert_eq!(
                failed,
                Err(if cancel {
                    Error::Cancelled
                } else {
                    Error::Limit
                })
            );
            assert_eq!(writer.finish().unwrap(), expected);
        }
    }
    #[test]
    fn actual_edge_roles_cannot_be_encoded_in_a_foreign_face_context() {
        let body = source();
        let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
        let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
        let t = TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, None).unwrap();
        let p = &s.patches()[0];
        let ends = p.rows.map(|i| s.rows()[i].endpoints());
        let ruling = t
            .edges()
            .iter()
            .find(|e| e.geometry == EdgeGeometry::Ruling)
            .unwrap();
        assert!(matches!(cap_row(ruling, &s), Err(Error::Unresolved)));
        let foreign = t
            .edges()
            .iter()
            .find(|e| matches!(e.geometry,EdgeGeometry::Row(r) if !p.rows.contains(&r)))
            .unwrap();
        assert!(matches!(
            tensor_uv(p.rows, ends, foreign),
            Err(Error::Unresolved)
        ));
        let foreign = t
            .edges()
            .iter()
            .find(|e| {
                e.geometry == EdgeGeometry::Ruling
                    && e.ends
                        .iter()
                        .any(|v| !ends.iter().flatten().any(|x| x == v))
            })
            .unwrap();
        assert!(matches!(
            tensor_uv(p.rows, ends, foreign),
            Err(Error::Unresolved)
        ));
    }
    #[test]
    fn output_accounting_and_streaming_numeric_fields_preserve_failure() {
        assert_eq!(output_estimate(1, 2, 3).unwrap(), 1024 + 5 * 4096);
        for (c, e, f) in [
            (usize::MAX, 0, 0),
            (0, usize::MAX, 0),
            (0, 0, usize::MAX),
            (usize::MAX / 1024, 1, 0),
            (usize::MAX / 1024, 0, 1),
        ] {
            assert_eq!(output_estimate(c, e, f), Err(Error::Limit));
        }
        assert_eq!(Multiplicities(3, 11).to_string(), "12,11,12");
        let values = KnotValues(&[0., 0.5, 1.]);
        assert_eq!(
            values.to_string(),
            "0.0000000000000000E0,5.0000000000000000E-1,1.0000000000000000E0"
        );
        struct Stop(usize);
        impl std::fmt::Write for Stop {
            fn write_str(&mut self, _: &str) -> std::fmt::Result {
                if self.0 == 0 {
                    Err(std::fmt::Error)
                } else {
                    self.0 -= 1;
                    Ok(())
                }
            }
        }
        use std::fmt::Write;
        for n in 0..5 {
            assert!(write!(&mut Stop(n), "{values}").is_err());
            assert!(write!(&mut Stop(n), "{}", Multiplicities(8, 11)).is_err());
        }
    }
    #[test]
    fn a_real_positive_source_wall_is_not_an_unearned_literal_sheet_separation() {
        use crate::{
            authored_flat_thread::{AuthoredThread, SourceOptions},
            wood_insert::{BoundaryOptions, Insert},
        };
        let original = source();
        let old = original.source();
        let mut spec = old.inner().spec();
        spec.major_radius_mm = old.outer().minimum_radius_mm().next_down();
        let inner = AuthoredThread::new(
            spec,
            SourceOptions {
                position_limit_mm: 1e-7,
                normal_angle_limit_rad: 1e-7,
            },
        )
        .unwrap();
        let original = Insert::new(*old.outer(), inner, f64::EPSILON * 0.5, None).unwrap();
        let body = original
            .boundary_with_interrupt(BoundaryOptions::default(), None)
            .unwrap();
        let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
        let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
        let t = TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, None).unwrap();
        let proof = s
            .projection_bounds(16_000_000, 64 * 1024 * 1024, None)
            .unwrap();
        assert!(proof.component_radius_mm[0][0] <= proof.component_radius_mm[1][1]);
        assert!(matches!(
            write_local_candidate(&t, 16_000_000, 64 * 1024 * 1024, None),
            Err(Error::Unresolved)
        ));
    }
    #[test]
    fn actual_local_insert_file_uses_shared_tensor_boundaries_and_two_pcurves_per_edge() {
        for kind in 0..3 {
            let body = if kind == 0 {
                super::super::tests::metric_graph()
            } else {
                super::super::tests::graph(kind == 2, 1)
            };
            let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
            let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
            let t = TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, None).unwrap();
            let text = write_local_candidate(&t, 16_000_000, 64 * 1024 * 1024, None).unwrap();
            assert_eq!(text.matches("=MANIFOLD_SOLID_BREP(").count(), 1);
            assert_eq!(text.matches("=CLOSED_SHELL(").count(), 1);
            assert_eq!(text.matches("=ADVANCED_FACE(").count(), t.faces().len());
            assert_eq!(
                text.matches("=B_SPLINE_SURFACE_WITH_KNOTS(").count(),
                s.patches().len()
            );
            assert_eq!(text.matches("=PCURVE(").count(), 2 * t.edges().len());
            assert_eq!(text.matches("=PLANE(").count(), 2);
            assert!(!text.contains("=LINE("));
            eprintln!("LOCAL_STEP {kind}: {}bytes", text.len());
        }
    }
}
#[cfg(test)]
mod append_tests {
    use super::*;
    fn seed() -> Writer {
        let mut writer = Writer::new();
        writer
            .append_sphere(
                crate::sphere::Sphere::new(1.).unwrap(),
                [20., 0., 0.],
                [1., 0., 0.],
                [0., 1., 0.],
            )
            .unwrap();
        writer
    }
    #[test]
    fn actual_public_append_preserves_mixed_members_and_retries_after_partial_failures() {
        let body = super::super::tests::metric_graph();
        let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
        let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
        let t = TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, None).unwrap();
        let polls = std::cell::Cell::new(0);
        let probe = || {
            polls.set(polls.get() + 1);
            false
        };
        let mut expected = seed();
        expected
            .append_owned_insert_local(&t, 16_000_000, 64 * 1024 * 1024, Some(&probe))
            .unwrap();
        let total = polls.get();
        assert!(total > 100);
        let expected = expected.finish().unwrap();
        let original = seed().finish().unwrap();
        let projection = s
            .projection_bounds(16_000_000, 64 * 1024 * 1024, None)
            .unwrap();
        let caps = t
            .cap_regularity_bounds(16_000_000, 64 * 1024 * 1024, None)
            .unwrap();
        let partial_budget = projection.work_used - s.work_used() + caps.work_used + 2000;
        for cancel in [false, true] {
            let mut writer = seed();
            let polls = std::cell::Cell::new(0);
            let stop = || {
                polls.set(polls.get() + 1);
                cancel && polls.get() == total - 20
            };
            assert_eq!(
                writer.append_owned_insert_local(
                    &t,
                    if cancel { 16_000_000 } else { partial_budget },
                    64 * 1024 * 1024,
                    Some(&stop)
                ),
                Err(if cancel {
                    Error::Cancelled
                } else {
                    Error::Limit
                })
            );
            // Rebuild the same seed for a direct publication check. A second
            // retry additionally checks restored IDs and retained-capacity paths.
            if !cancel {
                assert_eq!(writer.finish().unwrap(), original);
                continue;
            }
            writer
                .append_owned_insert_local(&t, 16_000_000, 64 * 1024 * 1024, None)
                .unwrap();
            assert_eq!(writer.finish().unwrap(), expected);
        }
        let mut mixed = seed();
        mixed
            .append_owned_insert_local(&t, 16_000_000, 64 * 1024 * 1024, None)
            .unwrap();
        mixed
            .append_sphere(
                crate::sphere::Sphere::new(2.).unwrap(),
                [-20., 0., 0.],
                [1., 0., 0.],
                [0., 1., 0.],
            )
            .unwrap();
        let text = mixed.finish().unwrap();
        let faces =
            crate::step_read::owned_insert::represented_member_faces_fixture(&text).unwrap();
        assert_eq!(faces.len(), 3);
        assert_eq!(
            faces.iter().filter(|f| f.len() == t.faces().len()).count(),
            1
        );
        let all = faces
            .iter()
            .flatten()
            .copied()
            .collect::<std::collections::BTreeSet<_>>();
        assert_eq!(all.len(), faces.iter().map(Vec::len).sum::<usize>());
        // The current independent single-body packet intentionally does not
        // confuse a mixed compound with one occupied insert.
        assert!(matches!(
            crate::step_read::owned_insert::parse(
                &text,
                crate::step_read::owned_insert::Options::default(),
                None
            ),
            Err(crate::step_read::owned_insert::Error::Unsupported)
        ));
    }
}

#[cfg(test)]
mod tiny_trim_tests;
