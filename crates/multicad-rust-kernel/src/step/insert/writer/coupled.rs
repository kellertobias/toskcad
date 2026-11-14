//! Exact-token SOURCE candidate emission. Incidence/support binding, not material.
use super::*;
use crate::step::insert::coupled::{Coefficient, topology::CoupledTopology};
use num_bigint::{BigInt, Sign};
use std::fmt::{self, Write as _};
struct Counter(usize);
impl fmt::Write for Counter {
    fn write_str(&mut self, s: &str) -> fmt::Result {
        self.0 = self.0.checked_add(s.len()).ok_or(fmt::Error)?;
        Ok(())
    }
}
fn entity(
    file: &mut StepFile,
    args: fmt::Arguments<'_>,
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<Entity, Error> {
    // An encoded-byte unit covers two bounded fmt visits. The first allocates
    // nothing and earns the exact byte charge BEFORE mutating the transaction.
    let mut count = Counter(0);
    writeln!(&mut count, "{}={args};", Entity(file.next)).map_err(|_| Error::Limit)?;
    w.charge(count.0 + 1)?;
    reserve(file, count.0, memory, w)?;
    Ok(file.entity(args))
}
fn reserve(
    file: &mut StepFile,
    additional: usize,
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<(), Error> {
    w.charge(1)?;
    let required = file
        .data
        .len()
        .checked_add(additional)
        .ok_or(Error::Limit)?;
    let old = file.data.capacity();
    if required <= old {
        return Ok(());
    }
    let desired = old.saturating_add((old / 2).max(64 * 1024));
    // Current memory includes the OLD buffer. A realloc may temporarily retain
    // BOTH buffers, so admit the complete new capacity, not just its growth.
    let target = desired.max(required).min(memory.remaining());
    if target < required {
        return Err(Error::Limit);
    }
    let before = memory.used;
    memory.add(target)?;
    // Bulk copy is metered in machine words; allocation in pages. ASCII
    // formatting visits are separately metered byte-for-byte by entity().
    w.charge(1 + file.data.len().div_ceil(size_of::<usize>()) + target.div_ceil(4096))?;
    file.data
        .try_reserve_exact(target - file.data.len())
        .map_err(|_| Error::Limit)?;
    // Keep the actual returned capacity. Pinned RawVec exact growth stores the
    // requested capacity; this postcheck does not depend on allocator RSS.
    memory.used = before - old;
    memory.add(file.data.capacity())?;
    Ok(())
}
struct Scalar<'a> {
    source: &'a Coefficient,
    text: String,
    next: u32,
}
struct Scalars<'a> {
    buckets: Vec<u32>,
    entries: Vec<Scalar<'a>>,
}
impl<'a> Scalars<'a> {
    fn prepare(count: usize, memory: &mut Memory, w: &mut Work<'_>) -> Result<Self, Error> {
        // CoupledPlan's default source decomposition has at most2*32768rows,
        // each12controls/3coordinates: every possible entry fits a u32 index.
        let n = (count * 2).next_power_of_two();
        let mut buckets = memory.allocate(n, w)?;
        buckets.resize(n, u32::MAX);
        let entries = memory.allocate(count, w)?;
        Ok(Self { buckets, entries })
    }
    fn index(
        &mut self,
        c: &'a Coefficient,
        memory: &mut Memory,
        w: &mut Work<'_>,
    ) -> Result<usize, Error> {
        let bytes = c.numerator_magnitude_le();
        w.charge(bytes.len() + 2)?;
        let mut hash = 0xcbf29ce484222325_u64;
        for b in std::iter::once(c.signum() as u8).chain(bytes.iter().copied()) {
            hash ^= u64::from(b);
            hash = hash.wrapping_mul(0x100000001b3);
        }
        let bucket = hash as usize & (self.buckets.len() - 1);
        let mut candidate = self.buckets[bucket];
        while candidate != u32::MAX {
            let e = &self.entries[candidate as usize];
            w.charge(1 + bytes.len() + e.source.numerator_magnitude_le().len())?;
            if e.source.signum() == c.signum() && e.source.numerator_magnitude_le() == bytes {
                return Ok(candidate as usize);
            }
            candidate = e.next;
        }
        let text = scalar_text(c, memory, w)?;
        let index = self.entries.len();
        self.entries.push(Scalar {
            source: c,
            text,
            next: self.buckets[bucket],
        });
        self.buckets[bucket] = index as u32;
        Ok(index)
    }
    fn point(
        &mut self,
        file: &mut StepFile,
        p: &'a [Coefficient; 3],
        memory: &mut Memory,
        w: &mut Work<'_>,
    ) -> Result<Entity, Error> {
        let a = self.index(&p[0], memory, w)?;
        let b = self.index(&p[1], memory, w)?;
        let c = self.index(&p[2], memory, w)?;
        entity(
            file,
            format_args!(
                "CARTESIAN_POINT('',({},{},{}))",
                self.entries[a].text, self.entries[b].text, self.entries[c].text
            ),
            memory,
            w,
        )
    }
}
fn scalar_text(c: &Coefficient, memory: &mut Memory, w: &mut Work<'_>) -> Result<String, Error> {
    // Finite encoded source coefficients have at most1120 numerator bits:
    // a finite binary64 coordinate cache plus its bounded dyadic denominator.
    // On the pinned64-bit bigint build, radix10 uses the<32-limb division
    // path. Charge worst half-word divisions and actual integer decimal bytes.
    memory.add(64 * 1024)?;
    w.charge(c.numerator_magnitude_le().len() + 1)?;
    let value = BigInt::from_bytes_le(
        if c.signum() < 0 {
            Sign::Minus
        } else if c.signum() > 0 {
            Sign::Plus
        } else {
            Sign::NoSign
        },
        c.numerator_magnitude_le(),
    );
    let words = (value.bits() as usize).div_ceil(64);
    w.charge(1 + words * (words + 1) + 2 * words)?;
    let text = format!("{value}.");
    w.charge(text.len())?;
    memory.add(text.capacity())?;
    memory.used -= 64 * 1024;
    Ok(text)
}
fn spline(
    file: &mut StepFile,
    degree: usize,
    points: &[Entity],
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<Entity, Error> {
    let m = degree + 1;
    entity(
        file,
        format_args!(
            "B_SPLINE_CURVE_WITH_KNOTS('',{degree},{},.UNSPECIFIED.,.F.,.F.,({m},{m}),(0.,1.),.UNSPECIFIED.)",
            Entities(points)
        ),
        memory,
        w,
    )
}
fn definition(
    file: &mut StepFile,
    support: Entity,
    curve: Entity,
    context: Entity,
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<Entity, Error> {
    let d = entity(
        file,
        format_args!("DEFINITIONAL_REPRESENTATION('',({curve}),{context})"),
        memory,
        w,
    )?;
    entity(file, format_args!("PCURVE('',{support},{d})"), memory, w)
}
fn uv_curve(
    file: &mut StepFile,
    support: Entity,
    uv: [[f64; 2]; 2],
    context: Entity,
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<Entity, Error> {
    let mut points = [Entity(0); 2];
    for (i, p) in uv.into_iter().enumerate() {
        points[i] = entity(
            file,
            format_args!("CARTESIAN_POINT('',({:.16E},{:.16E}))", p[0], p[1]),
            memory,
            w,
        )?;
    }
    let curve = spline(file, 1, &points, memory, w)?;
    definition(file, support, curve, context, memory, w)
}
fn cap_curve<'a>(
    file: &mut StepFile,
    support: Entity,
    row: &'a crate::step::insert::coupled::CoupledRow,
    context: Entity,
    scalars: &mut Scalars<'a>,
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<Entity, Error> {
    let mut points = memory.allocate(12, w)?;
    for p in row.coefficients() {
        let x = scalars.index(&p[0], memory, w)?;
        let y = scalars.index(&p[1], memory, w)?;
        points.push(entity(
            file,
            format_args!(
                "CARTESIAN_POINT('',({},{}))",
                scalars.entries[x].text, scalars.entries[y].text
            ),
            memory,
            w,
        )?);
    }
    let curve = spline(file, 11, &points, memory, w)?;
    memory.release(&points);
    drop(points);
    definition(file, support, curve, context, memory, w)
}
impl CoupledTopology<'_, '_> {
    /// Source candidate in exact numerator coordinates, with an ACTUAL STEP
    /// representation map converting those local coordinates to physical mm.
    #[allow(clippy::too_many_arguments)]
    pub fn write_mapped_candidate(
        &self,
        origin: [f64; 3],
        u: [f64; 3],
        v: [f64; 3],
        maximum_work: usize,
        maximum_bytes: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<String, Error> {
        let mut writer = Writer::new();
        writer.append_coupled_insert_candidate(
            self,
            origin,
            u,
            v,
            maximum_work,
            maximum_bytes,
            interrupt,
        )?;
        writer.finish().map_err(Error::from)
    }
}
impl Writer {
    #[allow(clippy::too_many_arguments)]
    pub fn append_coupled_insert_candidate(
        &mut self,
        t: &CoupledTopology<'_, '_>,
        origin: [f64; 3],
        u: [f64; 3],
        v: [f64; 3],
        maximum_work: usize,
        maximum_bytes: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<(), Error> {
        let n = super::mapped::pose(origin, u, v)?;
        self.append_typed(|file| {
            encode(
                file,
                t,
                maximum_work,
                maximum_bytes,
                interrupt,
                [origin, u, v, n],
            )
        })
    }
}
fn encode(
    file: &mut StepFile,
    t: &CoupledTopology<'_, '_>,
    maximum_work: usize,
    maximum_bytes: usize,
    interrupt: Option<&dyn Fn() -> bool>,
    map: [[f64; 3]; 4],
) -> Result<Entity, Error> {
    let mut w = Work::new(maximum_work, interrupt)?;
    w.charge(t.work_used() + 1)?;
    if maximum_bytes == 0 || maximum_bytes > 64 * 1024 * 1024 {
        return Err(Error::Invalid);
    }
    let mut memory = Memory::new(t.retained_bytes(), maximum_bytes)?;
    let controls = t.plan().rows().len() * 12;
    memory.add(file.data.capacity() + 4096)?;
    let mut scalars = Scalars::prepare(controls * 3, &mut memory, &mut w)?;
    let body = emit(file, t, &mut scalars, &mut memory, &mut w)?;
    let mapped = map_numerators(
        file,
        body,
        t.plan().rows()[0].coefficients()[0][0].fractional_bits(),
        map,
        &mut memory,
        &mut w,
    )?;
    // Standalone finish has fixed one-member protocol records. Compound
    // membership metadata remains governed by the generic Writer contract.
    reserve(file, 4096, &mut memory, &mut w)?;
    w.charge(4096)?;
    Ok(mapped)
}
fn emit<'a>(
    file: &mut StepFile,
    t: &'a CoupledTopology<'_, '_>,
    scalars: &mut Scalars<'a>,
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<Entity, Error> {
    let s = t.plan();
    let mut xyz = memory.allocate(t.vertices().len(), w)?;
    let mut vertices = memory.allocate(t.vertices().len(), w)?;
    for vertex in t.vertices() {
        let [row, p] = vertex.coefficient_reference();
        let point = scalars.point(file, &s.rows()[row].coefficients()[p], memory, w)?;
        xyz.push(point);
        vertices.push(entity(
            file,
            format_args!("VERTEX_POINT('',{point})"),
            memory,
            w,
        )?);
    }
    let mut row_points = memory.allocate(s.rows().len(), w)?;
    let mut curves = memory.allocate(t.edges().len(), w)?;
    for (index, row) in s.rows().iter().enumerate() {
        let ends = t.row_endpoints(index);
        let mut points = memory.allocate(12, w)?;
        for (i, p) in row.coefficients().iter().enumerate() {
            points.push(if i == 0 {
                xyz[ends[0]]
            } else if i == 11 {
                xyz[ends[1]]
            } else {
                scalars.point(file, p, memory, w)?
            });
        }
        curves.push(spline(file, 11, &points, memory, w)?);
        row_points.push(points);
    }
    for edge in &t.edges()[s.rows().len()..] {
        curves.push(spline(file, 1, &edge.ends.map(|i| xyz[i]), memory, w)?);
    }
    let mut supports = memory.allocate(t.faces().len(), w)?;
    for face in t.faces() {
        supports.push(match face.geometry{
        FaceGeometry::Tensor(i)=>{let[a,b]=t.patches()[i].rows.map(|r|&row_points[r]);entity(file,format_args!("B_SPLINE_SURFACE_WITH_KNOTS('',1,11,({},{}),.UNSPECIFIED.,.F.,.F.,.F.,(2,2),(12,12),(0.,1.),(0.,1.),.UNSPECIFIED.)",Entities(a),Entities(b)),memory, w)?},
        FaceGeometry::Cap(_)=>{let edge=&t.edges()[face.loops[0][0].edge];let row=edge.geometry.row_index().expect("cap incidence consists of actual source rows");let z=scalars.index(&s.rows()[row].coefficients()[0][2],memory,w)?;let origin=entity(file,format_args!("CARTESIAN_POINT('',(0.,0.,{}))",scalars.entries[z].text),memory, w)?;let n=entity(file,format_args!("DIRECTION('',(0.,0.,1.))"),memory, w)?;let x=entity(file,format_args!("DIRECTION('',(1.,0.,0.))"),memory, w)?;let axis=entity(file,format_args!("AXIS2_PLACEMENT_3D('',{origin},{n},{x})"),memory, w)?;entity(file,format_args!("PLANE('',{axis})"),memory, w)?},
    });
    }
    let context = entity(
        file,
        format_args!(
            "(GEOMETRIC_REPRESENTATION_CONTEXT(2)PARAMETRIC_REPRESENTATION_CONTEXT()REPRESENTATION_CONTEXT('2D SPACE',''))"
        ),
        memory,
        w,
    )?;
    let mut edge_pcurves = memory.allocate(t.edges().len(), w)?;
    for _ in t.edges() {
        edge_pcurves.push(memory.allocate(2, w)?);
    }
    for (i, face) in t.faces().iter().enumerate() {
        for wire in &face.loops {
            for c in wire {
                let edge = &t.edges()[c.edge];
                let pc = match face.geometry {
                    FaceGeometry::Cap(_) => {
                        let row = edge.geometry.row_index().expect("actual cap row incidence");
                        cap_curve(
                            file,
                            supports[i],
                            &s.rows()[row],
                            context,
                            scalars,
                            memory,
                            w,
                        )?
                    }
                    FaceGeometry::Tensor(p) => {
                        let patch = t.patches()[p];
                        let uv =
                            tensor_uv(patch.rows, patch.rows.map(|r| t.row_endpoints(r)), edge)?;
                        uv_curve(file, supports[i], uv, context, memory, w)?
                    }
                };
                edge_pcurves[c.edge].push(pc);
            }
        }
    }
    let mut edges = memory.allocate(t.edges().len(), w)?;
    for (i, edge) in t.edges().iter().enumerate() {
        let curve = entity(
            file,
            format_args!(
                "SURFACE_CURVE('',{},{},.CURVE_3D.)",
                curves[i],
                Entities(&edge_pcurves[i])
            ),
            memory,
            w,
        )?;
        let [a, b] = edge.ends.map(|v| vertices[v]);
        edges.push(entity(
            file,
            format_args!("EDGE_CURVE('',{a},{b},{curve},.T.)"),
            memory,
            w,
        )?);
    }
    let mut faces = memory.allocate(t.faces().len(), w)?;
    for (i, face) in t.faces().iter().enumerate() {
        let mut bounds = memory.allocate(face.loops.len(), w)?;
        for (j, wire) in face.loops.iter().enumerate() {
            let mut oriented = memory.allocate(wire.len(), w)?;
            for c in wire {
                let sense = if c.forward { ".T." } else { ".F." };
                oriented.push(entity(
                    file,
                    format_args!("ORIENTED_EDGE('',*,*,{},{sense})", edges[c.edge]),
                    memory,
                    w,
                )?);
            }
            let loop_id = entity(
                file,
                format_args!("EDGE_LOOP('',{})", Entities(&oriented)),
                memory,
                w,
            )?;
            let kind = if j == 0 {
                "FACE_OUTER_BOUND"
            } else {
                "FACE_BOUND"
            };
            bounds.push(entity(
                file,
                format_args!("{kind}('',{loop_id},.T.)"),
                memory,
                w,
            )?);
            memory.release(&oriented);
            drop(oriented);
        }
        let sense = if face.same_sense { ".T." } else { ".F." };
        faces.push(entity(
            file,
            format_args!(
                "ADVANCED_FACE('',{},{},{sense})",
                Entities(&bounds),
                supports[i]
            ),
            memory,
            w,
        )?);
        memory.release(&bounds);
        drop(bounds);
    }
    let shell = entity(
        file,
        format_args!("CLOSED_SHELL('',{})", Entities(&faces)),
        memory,
        w,
    )?;
    entity(
        file,
        format_args!("MANIFOLD_SOLID_BREP('',{shell})"),
        memory,
        w,
    )
}
fn map_numerators(
    file: &mut StepFile,
    body: Entity,
    bits: usize,
    [origin, u, v, n]: [[f64; 3]; 4],
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<Entity, Error> {
    let z = entity(
        file,
        format_args!("CARTESIAN_POINT('',(0.,0.,0.))"),
        memory,
        w,
    )?;
    let a = entity(file, format_args!("DIRECTION('',(0.,0.,1.))"), memory, w)?;
    let b = entity(file, format_args!("DIRECTION('',(1.,0.,0.))"), memory, w)?;
    let source_axis = entity(
        file,
        format_args!("AXIS2_PLACEMENT_3D('',{z},{a},{b})"),
        memory,
        w,
    )?;
    let length = entity(
        file,
        format_args!("(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.))"),
        memory,
        w,
    )?;
    let angle = entity(
        file,
        format_args!("(NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.))"),
        memory,
        w,
    )?;
    let solid_angle = entity(
        file,
        format_args!("(NAMED_UNIT(*) SI_UNIT($,.STERADIAN.) SOLID_ANGLE_UNIT())"),
        memory,
        w,
    )?;
    let context = entity(
        file,
        format_args!(
            "(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT(({length},{angle},{solid_angle})) REPRESENTATION_CONTEXT('',''))"
        ),
        memory,
        w,
    )?;
    let local = entity(
        file,
        format_args!("SHAPE_REPRESENTATION('',({body},{source_axis}),{context})"),
        memory,
        w,
    )?;
    let source = entity(
        file,
        format_args!("REPRESENTATION_MAP({source_axis},{local})"),
        memory,
        w,
    )?;
    let o = entity(
        file,
        format_args!(
            "CARTESIAN_POINT('',({:.16E},{:.16E},{:.16E}))",
            origin[0], origin[1], origin[2]
        ),
        memory,
        w,
    )?;
    let mut directions = [Entity(0); 3];
    for (i, p) in [u, v, n].into_iter().enumerate() {
        directions[i] = entity(
            file,
            format_args!("DIRECTION('',({:.16E},{:.16E},{:.16E}))", p[0], p[1], p[2]),
            memory,
            w,
        )?;
    }
    // 2^-bits =5^bits/10^bits, exactly, not a17-digit scale proposal.
    w.charge(1 + (3 * bits).div_ceil(64).pow(2))?;
    let digits = BigInt::from(5).pow(bits as u32).to_str_radix(10);
    let scale = format!("0.{}{digits}", "0".repeat(bits - digits.len()));
    w.charge(scale.len())?;
    let [x, y, z] = directions;
    let target = entity(
        file,
        format_args!(
            "(CARTESIAN_TRANSFORMATION_OPERATOR({x},{y},{o},{scale})CARTESIAN_TRANSFORMATION_OPERATOR_3D({z})FUNCTIONALLY_DEFINED_TRANSFORMATION('','')GEOMETRIC_REPRESENTATION_ITEM()REPRESENTATION_ITEM(''))"
        ),
        memory,
        w,
    )?;
    entity(
        file,
        format_args!("MAPPED_ITEM('',{source},{target})"),
        memory,
        w,
    )
}
#[cfg(test)]
mod tests;
