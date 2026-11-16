//! Literal replay of an independently occupied source with physical head topology.
//! Candidate encoding alone never grants decoded material or source-law authority.
use super::*;
use crate::step::pipe_hardware::numbers;
use std::fmt::{self, Write};
const EDIT_BYTES: usize = 1024 * 1024;
const EDITS: usize = 256;
const HEADER: &str = concat!(
    "ISO-10303-21;\nHEADER;\n",
    "FILE_DESCRIPTION(('MultiCAD clipped D11 physical head candidate'),'2;1');\n",
    "FILE_NAME('MultiCAD','1970-01-01T00:00:00',(''),(''),'MultiCAD','MultiCAD','');\n",
    "FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\nENDSEC;\nDATA;\n"
);
const TAIL: &str = "ENDSEC;\nEND-ISO-10303-21;\n";
#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub maximum_work: usize,
    pub maximum_bytes: usize,
    pub additional_resident_bytes: usize,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            maximum_work: 96_000_000,
            maximum_bytes: MAX_BYTES,
            additional_resident_bytes: 0,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Receipt {
    pub inherited_material_work: usize,
    pub fresh_encoding_work: usize,
    pub aggregate_work: usize,
    /// Includes the preflighted 1MiB head-edit/formatter arena; not RSS.
    pub accounted_peak_bytes: usize,
}
pub struct CandidateFile {
    text: String,
    receipt: Receipt,
}
impl CandidateFile {
    pub fn as_str(&self) -> &str {
        &self.text
    }
    pub fn receipt(&self) -> Receipt {
        self.receipt
    }
}
struct Count(usize);
impl Write for Count {
    fn write_str(&mut self, s: &str) -> fmt::Result {
        self.0 = self.0.checked_add(s.len()).ok_or(fmt::Error)?;
        Ok(())
    }
}
struct Edit {
    id: u32,
    body: String,
}
struct Encoder<'a> {
    edits: Vec<Edit>,
    next: u32,
    bytes: usize,
    work: Work<'a>,
}
impl Encoder<'_> {
    fn record(&mut self, id: u32, args: fmt::Arguments<'_>) -> Result<u32> {
        self.work.charge(32 + EDITS)?;
        if self.edits.len() == EDITS || self.edits.iter().any(|e| e.id == id) {
            return Err(Error::Limit);
        }
        let mut count = Count(0);
        count.write_fmt(args).map_err(|_| Error::Limit)?;
        self.work
            .charge(count.0.checked_mul(2).ok_or(Error::Limit)?)?;
        if count.0 > EDIT_BYTES - self.bytes {
            return Err(Error::Limit);
        }
        let mut body = String::new();
        body.try_reserve_exact(count.0).map_err(|_| Error::Limit)?;
        self.bytes = self
            .bytes
            .checked_add(body.capacity())
            .ok_or(Error::Limit)?;
        if self.bytes > EDIT_BYTES {
            return Err(Error::Limit);
        }
        body.write_fmt(args).map_err(|_| Error::Limit)?;
        self.edits.push(Edit { id, body });
        Ok(id)
    }
    fn add(&mut self, args: fmt::Arguments<'_>) -> Result<u32> {
        let id = self.next;
        self.next = self.next.checked_add(1).ok_or(Error::Limit)?;
        self.record(id, args)
    }
    fn xyz(&mut self, name: &str, xyz: [f64; 3]) -> Result<u32> {
        let mut strings = [String::new(), String::new(), String::new()];
        for (s, x) in strings.iter_mut().zip(xyz) {
            let (work, bytes) = numbers::cost(x).map_err(|_| Error::Invalid)?;
            self.work
                .charge(work.checked_mul(64).ok_or(Error::Limit)?)?;
            // The edit arena preflights all small scalar temporaries too. Even
            // three maximum-length finite binary64 decimal leaves fit 16KiB.
            if bytes > 4096 || self.bytes > EDIT_BYTES - 16 * 1024 {
                return Err(Error::Limit);
            }
            *s = numbers::literal(x).map_err(|_| Error::Invalid)?;
        }
        self.add(format_args!(
            "{name}('',({},{},{}))",
            strings[0], strings[1], strings[2]
        ))
    }
    fn plane(&mut self, p: Plane) -> Result<u32> {
        let point = self.xyz("CARTESIAN_POINT", p.origin)?;
        let normal = self.xyz("DIRECTION", p.normal)?;
        let reference = self.xyz("DIRECTION", p.reference)?;
        let axis = self.add(format_args!(
            "AXIS2_PLACEMENT_3D('',#{point},#{normal},#{reference})"
        ))?;
        self.add(format_args!("PLANE('',#{axis})"))
    }
    fn bound(&mut self, uses: &[(usize, bool)], edges: &[u32]) -> Result<u32> {
        self.work.charge(uses.len() + 1)?;
        let mut ids = Vec::new();
        ids.try_reserve_exact(uses.len())
            .map_err(|_| Error::Limit)?;
        if ids.capacity() * size_of::<u32>() > 1024 {
            return Err(Error::Limit);
        }
        for &(edge, forward) in uses {
            let edge = *edges.get(edge).ok_or(Error::Invalid)?;
            ids.push(self.add(format_args!(
                "ORIENTED_EDGE('',*,*,#{edge},{})",
                sense(forward)
            ))?);
        }
        let loop_id = self.add(format_args!("EDGE_LOOP('',({}))", Refs(&ids)))?;
        self.add(format_args!("FACE_OUTER_BOUND('',#{loop_id},.T.)"))
    }
}
fn sense(b: bool) -> &'static str {
    if b { ".T." } else { ".F." }
}
struct Refs<'a>(&'a [u32]);
impl fmt::Display for Refs<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, id) in self.0.iter().enumerate() {
            if i != 0 {
                f.write_char(',')?;
            }
            write!(f, "#{id}")?;
        }
        Ok(())
    }
}
fn one<'a, 'i>(
    work: &mut Work<'_>,
    records: &'a BTreeMap<u32, Vec<Entity<'i>>>,
    id: u32,
    name: &str,
) -> Result<&'a Entity<'i>> {
    work.charge(128)?;
    match records.get(&id).map(Vec::as_slice) {
        Some([e]) if e.name == name => Ok(e),
        _ => Err(Error::Invalid),
    }
}
fn reference(v: &Value<'_>) -> Result<u32> {
    super::super::records::reference(v)
}
/// Encode a single identity-frame candidate from the actual prerequisite.
/// All original 3D nets/curves and planar XY pcurves are replayed literally.
/// The former verification-cap topology becomes the physical annulus; its
/// physical support origin is an exact head descriptor, explicitly sewn to
/// the unchanged actual terminal 3D rim by the prerequisite's collar.
pub fn propose(
    proof: &OccupiedClippedHead<'_, '_, '_>,
    options: Options,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<CandidateFile> {
    if options.maximum_work == 0
        || options.maximum_work > 96_000_000
        || options.maximum_bytes == 0
        || options.maximum_bytes > MAX_BYTES
    {
        return Err(Error::Invalid);
    }
    let inherited = proof.receipt.aggregate_work;
    let resident = proof
        .receipt
        .accounted_peak_bytes
        .checked_add(options.additional_resident_bytes)
        .ok_or(Error::Limit)?;
    let fixed = resident.checked_add(EDIT_BYTES).ok_or(Error::Limit)?;
    if fixed > options.maximum_bytes {
        return Err(Error::Limit);
    }
    let parsed = proof.thread.projection().parsed();
    let members = parsed.document.represented_members();
    if members.len() != 1 || members[0].mapping_link().is_some() {
        return Err(Error::Unsupported);
    }
    let records = parsed.document.records();
    let mut work = Work {
        used: 0,
        maximum: options.maximum_work,
        interrupt,
    };
    work.charge(inherited)?;
    work.charge(128)?;
    let maximum = *records.keys().next_back().ok_or(Error::Invalid)?;
    let mut edits = Vec::new();
    edits.try_reserve_exact(EDITS).map_err(|_| Error::Limit)?;
    let edit_bytes = edits.capacity() * size_of::<Edit>();
    if edit_bytes > EDIT_BYTES - 16 * 1024 {
        return Err(Error::Limit);
    }
    let mut enc = Encoder {
        edits,
        next: maximum.checked_add(1).ok_or(Error::Limit)?,
        bytes: edit_bytes,
        work,
    };
    let head = &proof.head;
    if head.vertices.len() != 12 || head.edges.len() != 18 || head.faces.len() != 8 {
        return Err(Error::Unsupported);
    }
    let cap = one(
        &mut enc.work,
        records,
        head.auxiliary_cap_entity,
        "ADVANCED_FACE",
    )?;
    let [_, Value::List(bounds), support, _] = cap.args.as_slice() else {
        return Err(Error::Invalid);
    };
    let [cap_bound] = bounds.as_slice() else {
        return Err(Error::Invalid);
    };
    let cap_bound = reference(cap_bound)?;
    let old_bound = one(&mut enc.work, records, cap_bound, "FACE_OUTER_BOUND")?;
    let [_, wire, Value::Enum("T")] = old_bound.args.as_slice() else {
        return Err(Error::Invalid);
    };
    let wire = reference(wire)?;
    enc.record(cap_bound, format_args!("FACE_BOUND('',#{wire},.T.)"))?;
    let support = reference(support)?;
    let plane = one(&mut enc.work, records, support, "PLANE")?;
    let [_, axis] = plane.args.as_slice() else {
        return Err(Error::Invalid);
    };
    let axis = one(
        &mut enc.work,
        records,
        reference(axis)?,
        "AXIS2_PLACEMENT_3D",
    )?;
    let [_, origin, _, _] = axis.args.as_slice() else {
        return Err(Error::Invalid);
    };
    let origin = reference(origin)?;
    // This one verification-support leaf is replaced by the physical head
    // descriptor. It must not alias any unchanged original 3D polynomial leaf.
    for index in parsed
        .curves
        .iter()
        .flat_map(|c| c.poles.iter())
        .chain(parsed.tensors.iter().flat_map(|t| t.poles.iter()))
    {
        enc.work.charge(4)?;
        if parsed.points[*index].id == origin {
            return Err(Error::Unsupported);
        }
    }
    for face in &parsed.faces {
        enc.work.charge(4)?;
        if face.id != head.auxiliary_cap_entity
            && matches!(face.support, Support::Plane { height } if parsed.points[height].id == origin)
        {
            return Err(Error::Unsupported);
        }
    }
    enc.work.charge(head.faces.len())?;
    let low = head
        .faces
        .iter()
        .find(|f| !f.inner.is_empty())
        .ok_or(Error::Invalid)?;
    enc.work.charge(8192)?;
    let z = numbers::literal(low.plane.origin[2]).map_err(|_| Error::Invalid)?;
    enc.record(origin, format_args!("CARTESIAN_POINT('',(0.,0.,{z}))"))?;
    let mut vertices = [0u32; 12];
    for (out, p) in vertices.iter_mut().zip(&head.vertices) {
        let p = enc.xyz("CARTESIAN_POINT", *p)?;
        *out = enc.add(format_args!("VERTEX_POINT('',#{p})"))?;
    }
    let mut edges = [0u32; 18];
    for (out, [a, b]) in edges.iter_mut().zip(&head.edges) {
        let delta = std::array::from_fn(|j| head.vertices[*b][j] - head.vertices[*a][j]);
        for (j, d) in delta.iter().enumerate() {
            enc.work.charge(8192)?;
            if !numbers::exact_difference(head.vertices[*b][j], head.vertices[*a][j], *d)
                .map_err(|_| Error::Invalid)?
            {
                return Err(Error::Unsupported);
            }
        }
        let p = enc.xyz("CARTESIAN_POINT", head.vertices[*a])?;
        let d = enc.xyz("DIRECTION", delta)?;
        let vector = enc.add(format_args!("VECTOR('',#{d},1.)"))?;
        let line = enc.add(format_args!("LINE('',#{p},#{vector})"))?;
        *out = enc.add(format_args!(
            "EDGE_CURVE('',#{},#{},#{line},.T.)",
            vertices[*a], vertices[*b]
        ))?;
    }
    enc.work.charge(17)?;
    let mut faces = Vec::with_capacity(17);
    for f in &parsed.faces {
        if f.id != head.auxiliary_cap_entity {
            faces.push(f.id);
        }
    }
    for f in &head.faces {
        let bound = enc.bound(&f.outer, &edges)?;
        if f.inner.is_empty() {
            let support = enc.plane(f.plane)?;
            faces.push(enc.add(format_args!(
                "ADVANCED_FACE('',(#{bound}),#{support},{})",
                sense(f.forward)
            ))?);
        } else {
            enc.record(
                head.auxiliary_cap_entity,
                format_args!(
                    "ADVANCED_FACE('',(#{bound},#{cap_bound}),#{support},{})",
                    sense(f.forward)
                ),
            )?;
            faces.push(head.auxiliary_cap_entity);
        }
    }
    let solid = one(&mut enc.work, records, parsed.solid, "MANIFOLD_SOLID_BREP")?;
    let [_, shell] = solid.args.as_slice() else {
        return Err(Error::Invalid);
    };
    let shell = reference(shell)?;
    enc.record(shell, format_args!("CLOSED_SHELL('',({}))", Refs(&faces)))?;
    // Each ordinary record costs one original lookup and one bounded edit
    // lookup (at most256 entries), plus its actual encoded byte visits.
    let mut length = HEADER.len() + TAIL.len();
    for (id, components) in records {
        enc.work.charge(32 + EDITS)?;
        let edit = enc.edits.iter().find(|e| e.id == *id);
        let body = edit.map(|e| e.body.len()).unwrap_or_else(|| {
            components
                .iter()
                .map(|e| e.source.expect("actual neutral component").len())
                .sum::<usize>()
                + usize::from(components.len() > 1) * 2
        });
        length = length
            .checked_add(body + id.to_string().len() + 4)
            .ok_or(Error::Limit)?;
    }
    for e in &enc.edits {
        if e.id > maximum {
            length = length
                .checked_add(e.body.len() + e.id.to_string().len() + 4)
                .ok_or(Error::Limit)?;
        }
    }
    let peak = fixed.checked_add(length).ok_or(Error::Limit)?;
    if peak > options.maximum_bytes {
        return Err(Error::Limit);
    }
    enc.work.charge(length)?;
    let mut text = String::new();
    text.try_reserve_exact(length).map_err(|_| Error::Limit)?;
    let peak = fixed.checked_add(text.capacity()).ok_or(Error::Limit)?;
    if peak > options.maximum_bytes {
        return Err(Error::Limit);
    }
    text.push_str(HEADER);
    for (id, components) in records {
        enc.work.charge(32 + EDITS)?;
        if let Some(e) = enc.edits.iter().find(|e| e.id == *id) {
            writeln!(text, "#{id}={};", e.body).map_err(|_| Error::Limit)?;
        } else {
            write!(text, "#{id}=").map_err(|_| Error::Limit)?;
            if components.len() > 1 {
                text.push('(');
            }
            for e in components {
                let s = e.source.expect("actual neutral component");
                enc.work.charge(s.len())?;
                text.push_str(s);
            }
            if components.len() > 1 {
                text.push(')');
            }
            text.push_str(";\n");
        }
    }
    for e in &enc.edits {
        if e.id > maximum {
            enc.work.charge(e.body.len() + 32)?;
            writeln!(text, "#{}={};", e.id, e.body).map_err(|_| Error::Limit)?;
        }
    }
    text.push_str(TAIL);
    enc.work.charge(0)?;
    if text.len() != length {
        return Err(Error::Invalid);
    }
    Ok(CandidateFile {
        text,
        receipt: Receipt {
            inherited_material_work: inherited,
            fresh_encoding_work: enc.work.used - inherited,
            aggregate_work: enc.work.used,
            accounted_peak_bytes: peak,
        },
    })
}
