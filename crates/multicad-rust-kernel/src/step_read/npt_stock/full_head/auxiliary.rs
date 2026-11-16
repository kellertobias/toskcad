use super::*;
use std::fmt::{self, Write};
const HEADER: &str = concat!(
    "ISO-10303-21;\nHEADER;\n",
    "FILE_DESCRIPTION(('DERIVED clipped D11 verification auxiliary; not physical material'),'2;1');\n",
    "FILE_NAME('DERIVED','1970-01-01T00:00:00',(''),(''),'MultiCAD','MultiCAD','');\n",
    "FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\nENDSEC;\nDATA;\n"
);
const TAIL: &str = "ENDSEC;\nEND-ISO-10303-21;\n";
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
fn edit(
    edits: &mut Vec<Edit>,
    id: u32,
    args: fmt::Arguments<'_>,
    b: &mut Budget<'_>,
) -> Result<()> {
    b.step(128)?;
    if edits.len() == 16 || edits.iter().any(|e| e.id == id) {
        return Err(Error::Limit);
    }
    let mut count = Count(0);
    count.write_fmt(args).map_err(|_| Error::Limit)?;
    b.step(count.0.checked_mul(2).ok_or(Error::Limit)?)?;
    if count.0 > b.limit - b.bytes {
        return Err(Error::Limit);
    }
    let mut body = String::new();
    body.try_reserve_exact(count.0).map_err(|_| Error::Limit)?;
    b.bytes = b.bytes.checked_add(body.capacity()).ok_or(Error::Limit)?;
    b.peak = b.peak.max(b.bytes);
    if b.bytes > b.limit {
        return Err(Error::Limit);
    }
    body.write_fmt(args).map_err(|_| Error::Limit)?;
    edits.push(Edit { id, body });
    Ok(())
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
pub(super) struct Derived {
    pub text: String,
    pub head: HeadInput,
    pub bindings: Vec<ComponentBinding>,
    pub chart: CapChart,
    pub retained_without_derived: usize,
}
fn span(input: &str, s: &str) -> Result<Range<usize>> {
    let start = (s.as_ptr() as usize)
        .checked_sub(input.as_ptr() as usize)
        .ok_or(Error::Invalid)?;
    let end = start.checked_add(s.len()).ok_or(Error::Limit)?;
    if end > input.len() || input.get(start..end) != Some(s) {
        return Err(Error::Invalid);
    }
    Ok(start..end)
}
pub(super) fn derive(
    doc: &CompoundDocument<'_>,
    plan: graph::Plan<'_>,
    b: &mut Budget<'_>,
) -> Result<Derived> {
    let records = doc.records();
    b.step(128)?;
    let old_max = *records.keys().next_back().ok_or(Error::Invalid)?;
    let ids = std::array::from_fn::<_, 7, _>(|i| old_max.checked_add(i as u32 + 1));
    let [
        Some(point),
        Some(normal),
        Some(reference),
        Some(axis),
        Some(plane),
        Some(bound),
        Some(face),
    ] = ids
    else {
        return Err(Error::Limit);
    };
    let mut edits = b.allocate(16)?;
    edit(
        &mut edits,
        point,
        format_args!("CARTESIAN_POINT('',(0.,0.,{}))", plan.cap_height.spelling()),
        b,
    )?;
    edit(
        &mut edits,
        normal,
        format_args!("DIRECTION('',(0.,0.,1.))"),
        b,
    )?;
    edit(
        &mut edits,
        reference,
        format_args!("DIRECTION('',(1.,0.,0.))"),
        b,
    )?;
    edit(
        &mut edits,
        axis,
        format_args!("AXIS2_PLACEMENT_3D('',#{point},#{normal},#{reference})"),
        b,
    )?;
    edit(&mut edits, plane, format_args!("PLANE('',#{axis})"), b)?;
    edit(
        &mut edits,
        bound,
        format_args!("FACE_OUTER_BOUND('',#{},.T.)", plan.inner_wire),
        b,
    )?;
    edit(
        &mut edits,
        face,
        format_args!("ADVANCED_FACE('',(#{bound}),#{plane},.T.)"),
        b,
    )?;
    for &(id, representation) in &plan.rehomes {
        edit(
            &mut edits,
            id,
            format_args!("PCURVE('',#{plane},#{representation})"),
            b,
        )?;
    }
    let mut selected = b.allocate(10)?;
    selected.extend_from_slice(&plan.thread_faces);
    selected.push(face);
    edit(
        &mut edits,
        plan.shell,
        format_args!("CLOSED_SHELL('',({}))", Refs(&selected)),
        b,
    )?;
    let mut components = 0usize;
    let mut length = HEADER.len() + TAIL.len();
    for (&id, parts) in records {
        b.step(128)?;
        components = components.checked_add(parts.len()).ok_or(Error::Limit)?;
        let body = if let Some(e) = edits.iter().find(|e| e.id == id) {
            e.body.len()
        } else {
            parts
                .iter()
                .map(|e| e.source.expect("actual physical component").len())
                .sum::<usize>()
                + usize::from(parts.len() > 1) * 2
        };
        length = length
            .checked_add(body + id.to_string().len() + 4)
            .ok_or(Error::Limit)?;
    }
    for e in &edits {
        if e.id > old_max {
            length = length
                .checked_add(e.body.len() + e.id.to_string().len() + 4)
                .ok_or(Error::Limit)?;
        }
    }
    let mut bindings = b.allocate(components)?;
    b.step(length)?;
    if length > b.limit - b.bytes {
        return Err(Error::Limit);
    }
    let mut text = String::new();
    text.try_reserve_exact(length).map_err(|_| Error::Limit)?;
    b.bytes = b.bytes.checked_add(text.capacity()).ok_or(Error::Limit)?;
    b.peak = b.peak.max(b.bytes);
    if b.bytes > b.limit {
        return Err(Error::Limit);
    }
    text.push_str(HEADER);
    for (&id, parts) in records {
        b.step(128)?;
        if let Some(e) = edits.iter().find(|e| e.id == id) {
            b.step(e.body.len())?;
            writeln!(text, "#{id}={};", e.body).map_err(|_| Error::Limit)?;
        } else {
            write!(text, "#{id}=").map_err(|_| Error::Limit)?;
            if parts.len() > 1 {
                text.push('(');
            }
            for (component, e) in parts.iter().enumerate() {
                let source = e.source.expect("actual physical component");
                b.step(source.len().checked_mul(2).ok_or(Error::Limit)?)?;
                let actual = span(doc.literal_input(), source)?;
                let start = text.len();
                text.push_str(source);
                bindings.push(ComponentBinding {
                    entity: id,
                    component,
                    actual,
                    derived: start..text.len(),
                });
            }
            if parts.len() > 1 {
                text.push(')');
            }
            text.push_str(";\n");
        }
    }
    for e in &edits {
        if e.id > old_max {
            b.step(e.body.len() + 32)?;
            writeln!(text, "#{}={};", e.id, e.body).map_err(|_| Error::Limit)?;
        }
    }
    text.push_str(TAIL);
    b.step(0)?;
    if text.len() != length {
        return Err(Error::Invalid);
    }
    let mut head = plan.head;
    head.auxiliary_cap_entity = face;
    let physical_height = head
        .faces
        .iter()
        .find(|f| !f.inner.is_empty())
        .ok_or(Error::Invalid)?
        .plane
        .origin[2];
    let [lo, hi] = plan.cap_height.bounds();
    let displacement = crate::tapered_thread::interval::I::new(lo, hi)
        - crate::tapered_thread::interval::I::point(physical_height);
    let chart = CapChart {
        physical_support_entity: plan.physical_plane,
        derived_support_entity: plane,
        physical_height_mm: physical_height,
        derived_height_bounds_mm: [lo, hi],
        z_displacement_bounds_mm: [displacement.lo, displacement.hi],
    };
    let retained_without_derived = head::storage(&head)?
        .checked_add(bindings.capacity() * size_of::<ComponentBinding>())
        .ok_or(Error::Limit)?;
    Ok(Derived {
        text,
        head,
        bindings,
        chart,
        retained_without_derived,
    })
}
