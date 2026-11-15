//! Serialization of earned immutable literal geometry. This candidate has NOT
//! passed independent geometric readback and cannot authorize Native admission.
use super::*;
use crate::step::{Entities, Entity, StepFile};
use std::fmt::{self, Write};
mod numbers;
mod records;
use records::Records;
mod sealed {
    pub trait Sealed {}
}
pub trait CapLiteral: sealed::Sealed {
    fn events(&self) -> &CapEvents<'_, '_, '_, '_, '_>;
    fn segment(&self) -> usize;
    fn xyz(&self) -> &[[EncodedScalar; 3]];
    fn weights(&self) -> &[EncodedScalar];
    fn uv(&self) -> &[[EncodedScalar; 3]];
    fn uv_weights(&self) -> &[EncodedScalar];
    fn range(&self) -> &[EncodedScalar; 2];
    fn reversed(&self) -> bool;
    fn source_bounds(&self) -> [f64; 3];
    fn work_used(&self) -> usize;
    fn retained_bytes(&self) -> usize;
}
macro_rules! cap_impl {
    ($t:ident,$reverse:expr) => {
        impl sealed::Sealed for $t<'_, '_, '_, '_, '_, '_> {}
        impl CapLiteral for $t<'_, '_, '_, '_, '_, '_> {
            fn events(&self) -> &CapEvents<'_, '_, '_, '_, '_> {
                self.source()
            }
            fn segment(&self) -> usize {
                self.segment_id()
            }
            fn xyz(&self) -> &[[EncodedScalar; 3]] {
                self.points()
            }
            fn weights(&self) -> &[EncodedScalar] {
                self.weights()
            }
            fn uv(&self) -> &[[EncodedScalar; 3]] {
                self.surface_points()
            }
            fn uv_weights(&self) -> &[EncodedScalar] {
                self.surface_weights()
            }
            fn range(&self) -> &[EncodedScalar; 2] {
                self.parameter_range()
            }
            fn reversed(&self) -> bool {
                $reverse
            }
            fn source_bounds(&self) -> [f64; 3] {
                self.source_bounds()
            }
            fn work_used(&self) -> usize {
                self.work_used()
            }
            fn retained_bytes(&self) -> usize {
                self.retained_bytes()
            }
        }
    };
}
cap_impl!(BoundAffineCap, true);
cap_impl!(BoundRootCap, false);
/// Whole-file budget is distinct from the existing per-query budgets.
#[derive(Clone, Copy)]
pub struct FileOptions {
    pub maximum_encoding_work: usize,
    pub maximum_stage_work: usize,
    pub maximum_bytes: usize,
}
impl Default for FileOptions {
    fn default() -> Self {
        Self {
            maximum_encoding_work: 128_000_000,
            maximum_stage_work: 16_000_000,
            maximum_bytes: 64 * 1024 * 1024,
        }
    }
}
pub struct StepCandidate {
    text: String,
    encoding_work: usize,
    maximum_stage_work: usize,
    stages: usize,
    known_work: usize,
    bytes: usize,
    peak: usize,
}
impl StepCandidate {
    pub(crate) fn into_string(self) -> String {
        self.text
    }
    pub fn as_str(&self) -> &str {
        &self.text
    }
    pub fn maximum_stage_work_used(&self) -> usize {
        self.maximum_stage_work
    }
    pub fn encoding_stages(&self) -> usize {
        self.stages
    }
    pub fn encoding_work_used(&self) -> usize {
        self.encoding_work
    }
    /// Shared original prerequisites counted once; includes this encoding stage.
    pub fn aggregate_known_work_used(&self) -> usize {
        self.known_work
    }
    pub fn peak_accounted_bytes(&self) -> usize {
        self.peak
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
}
/// Separately budgeted encoding stages, one per face/edge plus small shared
/// header, vertex and shell stages. `maximum_stage_work` applies to EACH stage and `maximum_encoding_work`
/// bounds their sum;
/// reported total work covers ALL stages, without pretending the file is one
///16M operation. Pole blocks and buffer growth also have explicit stages. The byte ledger and numeric cache span the complete file. All supports must be earned beforehand
/// and belong to these EXACT immutable arenas. No source recipe is serialized
/// as a substitute for geometric coefficients.
pub fn encode_candidate(
    topology: &CapTopology<'_, '_, '_, '_, '_, '_, '_>,
    caps: &[&dyn CapLiteral],
    helices: &[&CompactHelix<'_, '_, '_, '_, '_, '_>],
    plains: &[&PlainWall<'_, '_, '_, '_, '_>],
    options: FileOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<StepCandidate, Error> {
    let compact = topology.source();
    let events = topology.events();
    let inherited = topology.retained_bytes()
        + compact.retained_bytes()
        + compact.surfaces().retained_bytes()
        + events.retained_bytes();
    let mut records = Records::new(options, inherited, interrupt)?;
    // Typed lookup arenas and all small index/Entity lists are reserved before
    // emission. Numeric/spline temporary arrays remain inside numeric scratch.
    records.reserve_indices(
        topology.vertices().len()
            + topology.edges().len() * 8
            + topology.faces().len() * 8
            + events.segments().len() * 2,
    )?;
    let mut known = topology.work_used();
    let mut caps_by_id = vec![None; events.segments().len()];
    for &cap in caps {
        records.w.charge(1)?;
        if !std::ptr::eq(cap.events(), events)
            || cap.segment() >= caps_by_id.len()
            || caps_by_id[cap.segment()].is_some()
        {
            return Err(Error::Invalid);
        }
        caps_by_id[cap.segment()] = Some(cap);
        records.hold_bytes(cap.retained_bytes())?;
        known = known
            .checked_add(cap.work_used() - events.work_used())
            .ok_or(Error::Limit)?;
    }
    records.w.charge(caps_by_id.len())?;
    if caps_by_id.iter().any(Option::is_none) {
        return Err(Error::Invalid);
    }
    let mut helix_by_id = vec![None; compact.edges().len()];
    for &helix in helices {
        records.w.charge(1)?;
        if !std::ptr::eq(helix.vertices_source(), events.source())
            || !std::ptr::eq(helix.source(), compact)
            || helix.edge() >= helix_by_id.len()
            || helix_by_id[helix.edge()].is_some()
        {
            return Err(Error::Invalid);
        }
        helix_by_id[helix.edge()] = Some(helix);
        records.hold_bytes(helix.retained_bytes())?;
        let prefix =
            compact.work_used() + events.source().work_used() - compact.source().work_used();
        known = known
            .checked_add(helix.work_used() - prefix)
            .ok_or(Error::Limit)?;
    }
    let mut plain_by_face = vec![None; compact.faces().len()];
    for &plain in plains {
        records.w.charge(1)?;
        if !std::ptr::eq(plain.source(), compact)
            || plain.face() >= plain_by_face.len()
            || plain_by_face[plain.face()].is_some()
        {
            return Err(Error::Invalid);
        }
        plain_by_face[plain.face()] = Some(plain);
        records.hold_bytes(plain.retained_bytes())?;
        known = known
            .checked_add(plain.work_used() - compact.work_used())
            .ok_or(Error::Limit)?;
    }
    let context=records.entity(format_args!("(GEOMETRIC_REPRESENTATION_CONTEXT(2)PARAMETRIC_REPRESENTATION_CONTEXT()REPRESENTATION_CONTEXT('2D SPACE',''))"))?;
    records.finish_stage()?;
    let mut surfaces = Vec::with_capacity(compact.faces().len());
    for (face, owner) in compact.faces().iter().enumerate() {
        surfaces.push(match owner.kind {
            ChartFaceKind::Thread(s) => {
                let s = &compact.surfaces().surfaces()[s];
                records.surface(
                    s.points(),
                    s.dimensions(),
                    s.degrees(),
                    s.row_weights(),
                    s.knots(),
                )?
            }
            ChartFaceKind::PlainWall { .. } => {
                let s = plain_by_face[face].ok_or(Error::Invalid)?;
                records.surface(s.points(), s.dimensions(), s.degrees(), &[], s.knots())?
            }
            ChartFaceKind::Plane(_) => {
                let c = &topology.faces()[face].loops()[0][0];
                let v = topology.edges()[c.edge].vertices[0];
                let z =
                    topology.coordinates(v).ok_or(Error::Unresolved)?[2].number(&mut records.w)?;
                let origin = records.point(&[N::integer(0), N::integer(0), z], 3)?;
                let zaxis = records.entity(format_args!("DIRECTION('',(0.,0.,1.))"))?;
                let xaxis = records.entity(format_args!("DIRECTION('',(1.,0.,0.))"))?;
                let axis = records.entity(format_args!(
                    "AXIS2_PLACEMENT_3D('',{origin},{zaxis},{xaxis})"
                ))?;
                records.entity(format_args!("PLANE('',{axis})"))?
            }
        });
        records.finish_stage()?;
    }
    let mut vertices = vec![None; topology.vertices().len()];
    for edge in topology.edges() {
        for &v in &edge.vertices {
            if vertices[v].is_none() {
                let p = topology.coordinates(v).ok_or(Error::Unresolved)?;
                let point = records.encoded_point(p, 3, None)?;
                vertices[v] = Some(records.entity(format_args!("VERTEX_POINT('',{point})"))?);
            }
        }
    }
    records.finish_stage()?;
    let mut edges = Vec::with_capacity(topology.edges().len());
    for (id, edge) in topology.edges().iter().enumerate() {
        let mut peers = Vec::new();
        for (face, owner) in topology.faces().iter().enumerate() {
            for c in owner.loops().iter().flatten() {
                records.w.charge(1)?;
                if c.edge == id {
                    peers.push(face)
                }
            }
        }
        if peers.len() != 2 {
            return Err(Error::Unresolved);
        }
        let (curve, pcurves, same_sense) = match edge.kind {
            TrimEdgeKind::CapSegment(segment) => {
                let cap = caps_by_id[segment].ok_or(Error::Invalid)?;
                let xyz = records.curve(
                    cap.xyz(),
                    cap.weights(),
                    &[N::integer(0), N::integer(1)],
                    cap.xyz().len() - 1,
                    3,
                    None,
                )?;
                let range = [
                    cap.range()[0].number(&mut records.w)?,
                    cap.range()[1].number(&mut records.w)?,
                ];
                let owner = &events.segments()[segment];
                let chart = &compact.surfaces().cells()[owner.source_cell()];
                let mut pcs = Vec::new();
                for face in peers {
                    let uv = match topology.faces()[face].kind {
                        ChartFaceKind::Thread(s) if s == chart.surface() => records.curve(
                            cap.uv(),
                            cap.uv_weights(),
                            &range,
                            cap.uv().len() - 1,
                            2,
                            Some(N::integer(chart.turn_shift())),
                        )?,
                        ChartFaceKind::Plane(_) => records.curve(
                            cap.xyz(),
                            cap.weights(),
                            &[N::integer(0), N::integer(1)],
                            cap.xyz().len() - 1,
                            2,
                            None,
                        )?,
                        _ => return Err(Error::Unresolved),
                    };
                    pcs.push(records.pcurve(surfaces[face], uv, context)?);
                }
                (xyz, pcs, !cap.reversed())
            }
            TrimEdgeKind::Original(original) => match compact.edges()[original].kind {
                ChartEdgeKind::Helix { .. } => {
                    let helix = helix_by_id[original].ok_or(Error::Invalid)?;
                    let knots = records.binary_knots(helix.knots())?;
                    let xyz = records.curve(
                        helix.points(),
                        std::slice::from_ref(helix.weight()),
                        &knots,
                        10,
                        3,
                        None,
                    )?;
                    let mut pcs = Vec::new();
                    for pc in helix.pcurves() {
                        if !peers.contains(&pc.face()) {
                            return Err(Error::Unresolved);
                        }
                        let a = knots[0].clone();
                        let b = knots.last().ok_or(Error::Unresolved)?.clone();
                        let shift = pc.v_shift().number(&mut records.w)?;
                        let u = binary(pc.u(), &mut records.w)?;
                        let shifted_a = a.add(&shift, &mut records.w)?;
                        let shifted_b = b.add(&shift, &mut records.w)?;
                        let uv = records.line([[u.clone(), shifted_a], [u, shifted_b]], &[a, b])?;
                        pcs.push(records.pcurve(surfaces[pc.face()], uv, context)?);
                    }
                    (xyz, pcs, helix.same_sense())
                }
                ChartEdgeKind::Original(old) => match compact.source().edges()[old].curve {
                    ClosedCurve::PlainCircle { end, quarter } => {
                        let face=*peers.iter().find(|&&f|matches!(topology.faces()[f].kind,ChartFaceKind::PlainWall {quarter:q} if q==quarter)).ok_or(Error::Unresolved)?;
                        let wall = plain_by_face[face].ok_or(Error::Invalid)?;
                        let width = wall.dimensions()[1];
                        let row = usize::from(end == End::End);
                        let points = &wall.points()[row * width..(row + 1) * width];
                        let knots = records.binary_knots(wall.knots())?;
                        let xyz = records.curve(points, &[], &knots, 10, 3, None)?;
                        let mut pcs = Vec::new();
                        for peer in peers {
                            let uv = if peer == face {
                                records.line(
                                    [
                                        [N::integer(row as i64), knots[0].clone()],
                                        [N::integer(row as i64), knots.last().unwrap().clone()],
                                    ],
                                    &[knots[0].clone(), knots.last().unwrap().clone()],
                                )?
                            } else {
                                records.curve(points, &[], &knots, 10, 2, None)?
                            };
                            pcs.push(records.pcurve(surfaces[peer], uv, context)?);
                        }
                        (xyz, pcs, true)
                    }
                    ClosedCurve::PlainRail { quarter } => {
                        let points = [
                            topology
                                .coordinates(edge.vertices[0])
                                .ok_or(Error::Unresolved)?
                                .clone(),
                            topology
                                .coordinates(edge.vertices[1])
                                .ok_or(Error::Unresolved)?
                                .clone(),
                        ];
                        let knots = [N::integer(0), N::integer(1)];
                        let xyz = records.curve(&points, &[], &knots, 1, 3, None)?;
                        let mut pcs = Vec::new();
                        for peer in peers {
                            let ChartFaceKind::PlainWall { quarter: q } =
                                topology.faces()[peer].kind
                            else {
                                return Err(Error::Unresolved);
                            };
                            let v = binary(
                                if q == quarter {
                                    f64::from(q) / 4.
                                } else {
                                    f64::from(q + 1) / 4.
                                },
                                &mut records.w,
                            )?;
                            let uv = records
                                .line([[N::integer(0), v.clone()], [N::integer(1), v]], &knots)?;
                            pcs.push(records.pcurve(surfaces[peer], uv, context)?);
                        }
                        (xyz, pcs, true)
                    }
                    _ => return Err(Error::Unresolved),
                },
            },
        };
        let surface_curve = records.entity(format_args!(
            "SURFACE_CURVE('',{curve},{},.CURVE_3D.)",
            Entities(&pcurves)
        ))?;
        let a = vertices[edge.vertices[0]].ok_or(Error::Unresolved)?;
        let b = vertices[edge.vertices[1]].ok_or(Error::Unresolved)?;
        edges.push(records.entity(format_args!(
            "EDGE_CURVE('',{a},{b},{surface_curve},{})",
            if same_sense { ".T." } else { ".F." }
        ))?);
        records.finish_stage()?;
    }
    let mut faces = Vec::new();
    for (face, owner) in topology.faces().iter().enumerate() {
        let mut bounds = Vec::new();
        for (ring, coedges) in owner.loops().iter().enumerate() {
            let mut uses = Vec::new();
            for c in coedges {
                uses.push(records.entity(format_args!(
                    "ORIENTED_EDGE('',*,*,{},{})",
                    edges[c.edge],
                    if c.sense { ".T." } else { ".F." }
                ))?)
            }
            let wire = records.entity(format_args!("EDGE_LOOP('',{})", Entities(&uses)))?;
            bounds.push(records.entity(format_args!(
                "{}('',{wire},.T.)",
                if ring == 0 {
                    "FACE_OUTER_BOUND"
                } else {
                    "FACE_BOUND"
                }
            ))?);
        }
        faces.push(records.entity(format_args!(
            "ADVANCED_FACE('',{},{},{})",
            Entities(&bounds),
            surfaces[face],
            if owner.same_sense { ".T." } else { ".F." }
        ))?);
    }
    let shell = records.entity(format_args!("CLOSED_SHELL('',{})", Entities(&faces)))?;
    let solid = records.entity(format_args!("MANIFOLD_SOLID_BREP('',{shell})"))?;
    records.finish(solid, known)
}

#[cfg(test)]
mod tests;
