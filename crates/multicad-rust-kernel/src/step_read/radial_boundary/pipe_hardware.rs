//! Independent actual-file finite-head decoding. Only the explicitly sewn
//! reference earns occupied admission; original file operators remain retained.
use super::*;
use crate::radial_tensor::Surface;
use crate::tapered_thread::interval::I;
use crate::{
    axial_pipe_thread::{
        Error as GeometryError,
        authored::hardware::shell::{self, occupied},
    },
    step::pipe_hardware::numbers as exact_numbers,
};
use std::{cell::Cell, sync::Arc};
pub(crate) mod budget;
mod records;
pub(crate) use records::Records;
type GeometryResult<T> = std::result::Result<T, GeometryError>;
// Composite query: cold parser64M, decoder64M, retained material96M.
// Child ceilings and the75M delegated material reservation remain unchanged.
const MAX_WORK: usize = 64_000_000 + 64_000_000 + 96_000_000;
#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub maximum_work: usize,
    pub maximum_bytes: usize,
    pub additional_resident_bytes: usize,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            maximum_work: MAX_WORK,
            maximum_bytes: 64 * 1024 * 1024,
            additional_resident_bytes: 0,
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub enum ParameterMapping {
    /// Actual knot domain[0,n], reference domain[0,1]. STEP=n*reference;
    /// original STEP derivatives pull back with n, separately for each axis.
    Spline {
        spans: usize,
    },
    Tensor {
        spans: [usize; 2],
    },
    /// Actual LINE has unit-length parameter; common finite segment parameter
    /// is[0,1]. The exact endpoint vector is retained, not a sampled tangent.
    Line {
        delta: [f64; 3],
        length_bounds_mm: [f64; 2],
    },
    AngularIdentity,
    /// STEP cone v is slant distance, reference z is along its canonical axis.
    /// STEP_v=z/cos(angle); the derivative pullback includes this scale.
    Cone {
        axial_to_slant_scale: [f64; 2],
    },
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Target {
    Curve(usize),
    Surface(usize),
}
#[derive(Clone, Debug)]
pub struct ParameterMap {
    target: Target,
    entity: u32,
    mapping: ParameterMapping,
}
impl ParameterMap {
    pub fn target(&self) -> Target {
        self.target
    }
    pub fn actual_entity(&self) -> u32 {
        self.entity
    }
    pub fn mapping(&self) -> &ParameterMapping {
        &self.mapping
    }
    /// Pull back an independently enclosed actual STEP tangent into the
    /// reference chart. Axis0 is curve/angular/u; axis1 is tensor-v/cone-axial.
    /// This applies the earned chain factor; it does not certify caller input.
    pub fn reference_derivative_bounds(
        &self,
        step: [[f64; 2]; 3],
        axis: usize,
    ) -> GeometryResult<[[f64; 2]; 3]> {
        let scale = match self.mapping {
            ParameterMapping::Spline { spans } if axis == 0 => I::point(spans as f64),
            ParameterMapping::Tensor { spans } if axis < 2 => I::point(spans[axis] as f64),
            ParameterMapping::Line {
                length_bounds_mm, ..
            } if axis == 0 => I::new(length_bounds_mm[0], length_bounds_mm[1]),
            ParameterMapping::AngularIdentity if axis == 0 => I::point(1.),
            ParameterMapping::Cone { .. } if axis == 0 => I::point(1.),
            ParameterMapping::Cone {
                axial_to_slant_scale,
            } if axis == 1 => I::new(axial_to_slant_scale[0], axial_to_slant_scale[1]),
            _ => return Err(GeometryError::Invalid),
        }
        .positive()?;
        let mut output = [[0.; 2]; 3];
        for k in 0..3 {
            let bound = (I::new(step[k][0], step[k][1]).checked()? * scale).checked()?;
            output[k] = [bound.lo, bound.hi];
        }
        Ok(output)
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Receipt {
    pub original_file_bytes: usize,
    pub parsed_record_count: usize,
    pub parsed_value_count: usize,
    pub maximum_argument_depth: usize,
    pub parser_work_reserved: usize,
    pub parser_scan_work_used: usize,
    pub delegated_parser_work_reserved: usize,
    pub decoder_work_reserved: usize,
    pub occupied_work_reserved: usize,
    pub aggregate_work_reserved: usize,
    pub parser_peak_bytes_reserved: usize,
    pub decode_peak_bytes_reserved: usize,
    pub occupied_peak_bytes_reserved: usize,
}
pub struct ReadHead {
    pub(crate) records: Records,
    material: occupied::OccupiedHead,
    maps: Vec<ParameterMap>,
    receipt: Receipt,
    retained_bytes: usize,
}
impl ReadHead {
    pub fn original_input(&self) -> &str {
        self.records.input()
    }
    pub fn material(&self) -> &occupied::OccupiedHead {
        &self.material
    }
    pub fn parameter_maps(&self) -> &[ParameterMap] {
        &self.maps
    }
    pub fn receipt(&self) -> Receipt {
        self.receipt
    }
    pub fn retained_bytes_upper_bound(&self) -> usize {
        self.retained_bytes
    }
}
fn geometry_error(e: ReadError) -> GeometryError {
    match e {
        ReadError::Invalid => GeometryError::Invalid,
        ReadError::Unsupported => GeometryError::Unsupported,
        ReadError::Limit => GeometryError::Limit,
    }
}
fn checkpoint(interrupt: Option<&dyn Fn() -> bool>) -> GeometryResult<()> {
    if interrupt.is_some_and(|f| f()) {
        Err(GeometryError::Cancelled)
    } else {
        Ok(())
    }
}
pub fn read(
    input: &str,
    tolerance: f64,
    minimum: f64,
    options: Options,
    interrupt: Option<&dyn Fn() -> bool>,
) -> GeometryResult<ReadHead> {
    if options.maximum_work == 0
        || options.maximum_work > MAX_WORK
        || options.maximum_bytes == 0
        || options.maximum_bytes > 64 * 1024 * 1024
    {
        return Err(GeometryError::Invalid);
    }
    checkpoint(interrupt)?;
    let cold = budget::reserve(
        input,
        options.maximum_work.min(64_000_000),
        options.maximum_bytes,
        interrupt,
    )?;
    let record_storage = records::storage(cold.represented_record_count, cold.components)?;
    let parser_peak = cold
        .parser_peak_bytes_reserved
        .checked_add(options.additional_resident_bytes)
        .ok_or(GeometryError::Limit)?;
    // Actual decoded pole copies arise from reference leaves, not every scalar
    // in the arena.128B/reference encloses net/curve copying, Vec growth and
    // transient Vec-to-Arc conversion. Tables/fans are separately bounded.
    let decode_peak = cold
        .retained_arena_bytes_reserved
        .checked_add(
            cold.references
                .checked_mul(128)
                .ok_or(GeometryError::Limit)?,
        )
        .and_then(|n| n.checked_add(input.len()))
        .and_then(|n| n.checked_add(record_storage))
        .and_then(|n| n.checked_add(2 * 1024 * 1024))
        .and_then(|n| n.checked_add(options.additional_resident_bytes))
        .ok_or(GeometryError::Limit)?;
    let decoder_work = input
        .len()
        .checked_mul(4)
        .and_then(|n| n.checked_add(cold.values.checked_mul(64)?))
        .and_then(|n| n.checked_add(cold.represented_record_count.checked_mul(128)?))
        .ok_or(GeometryError::Limit)?;
    let preceding = cold
        .aggregate_work_reserved
        .checked_add(decoder_work)
        .ok_or(GeometryError::Limit)?;
    let occupied_allowance = options
        .maximum_work
        .checked_sub(preceding)
        .ok_or(GeometryError::Limit)?
        .min(96_000_000);
    #[cfg(test)]
    eprintln!(
        "READ_RESERVATION parser_peak={parser_peak} decode_peak={decode_peak} preceding_work={preceding} decoder_work={decoder_work} occupied_allowance={occupied_allowance} caller_resident={} maxbytes={} maxwork={}",
        options.additional_resident_bytes, options.maximum_bytes, options.maximum_work
    );
    if parser_peak > options.maximum_bytes
        || decode_peak > options.maximum_bytes
        || decoder_work > 64_000_000
        || occupied_allowance < 75_000_000
    {
        return Err(GeometryError::Limit);
    }
    checkpoint(interrupt)?;
    let failure = Cell::new(None);
    let parsed =
        compound::read_members_retaining_arena(input, tolerance, true, true, |s, id, _| {
            let result = decode(s, id, interrupt);
            result.map_err(|e| {
                failure.set(Some(e));
                ReadError::Unsupported
            })
        });
    let (entities, mut members) =
        parsed.map_err(|e| failure.get().unwrap_or_else(|| geometry_error(e)))?;
    checkpoint(interrupt)?;
    validate_numbers(&entities, interrupt)?;
    if members.len() != 1 {
        return Err(GeometryError::Unsupported);
    }
    let Decoded {
        ownership,
        maps,
        member,
    } = members.pop().ok_or(GeometryError::Invalid)?;
    let records = Records::from_actual(input, &entities, member, interrupt)?;
    // Copy authentic operators/spans/bytes while the arena is charged, then
    // release that whole AST before the retained material construction.
    drop(entities);
    let additional = records
        .bytes()?
        .checked_add(input.len())
        .and_then(|n| n.checked_add(options.additional_resident_bytes))
        .and_then(|n| n.checked_add(2 * 1024 * 1024))
        .ok_or(GeometryError::Limit)?;
    let material = occupied::certify_with_options(
        ownership,
        tolerance,
        minimum,
        occupied::Options {
            maximum_work: occupied_allowance,
            maximum_bytes: options.maximum_bytes,
            additional_resident_bytes: additional,
        },
        interrupt,
    )?;
    let occupied_receipt = material.receipt();
    let aggregate = preceding
        .checked_add(occupied_receipt.aggregate_work_reserved)
        .ok_or(GeometryError::Limit)?;
    if aggregate > options.maximum_work {
        return Err(GeometryError::Limit);
    }
    let retained_bytes = occupied_receipt
        .retained_result_bytes_upper_bound
        .checked_sub(additional)
        .and_then(|n| n.checked_add(records.bytes().ok()?))
        .and_then(|n| n.checked_add(maps.capacity() * std::mem::size_of::<ParameterMap>()))
        .and_then(|n| n.checked_add(std::mem::size_of::<ReadHead>()))
        .ok_or(GeometryError::Limit)?;
    Ok(ReadHead {
        records,
        material,
        maps,
        receipt: Receipt {
            original_file_bytes: cold.input_bytes,
            parsed_record_count: cold.represented_record_count,
            parsed_value_count: cold.values,
            maximum_argument_depth: cold.maximum_argument_depth,
            parser_work_reserved: cold.aggregate_work_reserved,
            parser_scan_work_used: cold.scan_work_used,
            delegated_parser_work_reserved: cold.delegated_parser_work_reserved,
            decoder_work_reserved: decoder_work,
            occupied_work_reserved: occupied_receipt.aggregate_work_reserved,
            aggregate_work_reserved: aggregate,
            parser_peak_bytes_reserved: parser_peak,
            decode_peak_bytes_reserved: decode_peak,
            occupied_peak_bytes_reserved: occupied_receipt.phase_peak_bytes_reserved,
        },
        retained_bytes,
    })
}
struct Decoded {
    ownership: shell::Ownership,
    maps: Vec<ParameterMap>,
    member: u32,
}
fn canonical_axis(a: radial::Axis) -> GeometryResult<()> {
    if a.origin[0] != 0.
        || a.origin[1] != 0.
        || (a.normal != [0., 0., 1.] && a.normal != [0., 0., -1.])
        || a.reference != [1., 0., 0.]
    {
        Err(GeometryError::Unsupported)
    } else {
        Ok(())
    }
}
fn decode(
    s: &mut Shell<'_, '_>,
    member: u32,
    interrupt: Option<&dyn Fn() -> bool>,
) -> GeometryResult<Decoded> {
    let face_ids = s.member_faces(member).map_err(geometry_error)?;
    if face_ids.len() != 13 {
        return Err(GeometryError::Unsupported);
    }
    let mut input = shell::Input {
        vertices: vec![],
        curves: vec![],
        edges: vec![],
        surfaces: vec![],
        faces: vec![],
    };
    let (mut vertices, mut curves, mut edges, mut surfaces) = (
        BTreeMap::new(),
        BTreeMap::new(),
        BTreeMap::new(),
        BTreeMap::new(),
    );
    let mut maps = Vec::new();
    let mut controls = 0usize;
    for id in face_ids {
        checkpoint(interrupt)?;
        let (surface, forward, bounds) = face(s, id).map_err(geometry_error)?;
        if bounds.is_empty() || bounds.len() > 2 {
            return Err(GeometryError::Unsupported);
        }
        let surface_index = if let Some(&i) = surfaces.get(&surface) {
            i
        } else {
            if input.surfaces.len() >= 48 {
                return Err(GeometryError::Limit);
            }
            let i = input.surfaces.len();
            let value = if s.simple(surface).map_err(geometry_error)?.name
                == "B_SPLINE_SURFACE_WITH_KNOTS"
            {
                let (t, spans) = tensor(s, surface, &mut controls, interrupt)?;
                maps.push(ParameterMap {
                    target: Target::Surface(i),
                    entity: surface,
                    mapping: ParameterMapping::Tensor { spans },
                });
                Surface::Tensor(t)
            } else {
                let value = super::surface(s, surface).map_err(geometry_error)?;
                if let Surface::Cone {
                    axis,
                    semi_angle_rad,
                    ..
                } = value
                {
                    canonical_axis(axis)?;
                    let scale =
                        (I::point(1.) / I::point(semi_angle_rad).cos().positive()?).checked()?;
                    maps.push(ParameterMap {
                        target: Target::Surface(i),
                        entity: surface,
                        mapping: ParameterMapping::Cone {
                            axial_to_slant_scale: [scale.lo, scale.hi],
                        },
                    });
                }
                value
            };
            input.surfaces.push(value);
            surfaces.insert(surface, i);
            i
        };
        let mut outer = None;
        let mut inner = Vec::new();
        for bound in bounds {
            let (outside, _, _) = s.bound(bound, false).map_err(geometry_error)?;
            let cycle = ring(s, bound, outside).map_err(geometry_error)?;
            if cycle.len() > 80 {
                return Err(GeometryError::Limit);
            }
            let mut coedges = Vec::new();
            for (edge_id, direction) in cycle {
                checkpoint(interrupt)?;
                let edge_index = if let Some(&i) = edges.get(&edge_id) {
                    i
                } else {
                    if input.edges.len() >= 128 {
                        return Err(GeometryError::Limit);
                    }
                    let [Value::String(_), a, b, c, fwd] =
                        s.args(edge_id, "EDGE_CURVE").map_err(geometry_error)?
                    else {
                        return Err(GeometryError::Invalid);
                    };
                    let ids = [
                        reference(a).map_err(geometry_error)?,
                        reference(b).map_err(geometry_error)?,
                    ];
                    let curve_id = reference(c).map_err(geometry_error)?;
                    let fwd = sense(fwd).map_err(geometry_error)?;
                    let mut ends = [0; 2];
                    for (k, id) in ids.into_iter().enumerate() {
                        ends[k] = if let Some(&i) = vertices.get(&id) {
                            i
                        } else {
                            if input.vertices.len() >= 96 {
                                return Err(GeometryError::Limit);
                            }
                            let p = s.vertex(id).map_err(geometry_error)?;
                            let i = input.vertices.len();
                            input.vertices.push(p);
                            vertices.insert(id, i);
                            i
                        };
                    }
                    let curve = if let Some(&i) = curves.get(&curve_id) {
                        i
                    } else {
                        if input.curves.len() >= 128 {
                            return Err(GeometryError::Limit);
                        }
                        let i = input.curves.len();
                        let raw = s.simple(curve_id).map_err(geometry_error)?.name;
                        let (value, mapping) = match raw {
                            "LINE" => line(s, curve_id, ends, fwd, &input.vertices)?,
                            "B_SPLINE_CURVE_WITH_KNOTS" => {
                                let (c, n) = spline(s, curve_id, &mut controls, interrupt)?;
                                (
                                    shell::Curve::Original(c),
                                    ParameterMapping::Spline { spans: n },
                                )
                            }
                            "CIRCLE" => {
                                let c = super::curve(s, curve_id).map_err(geometry_error)?;
                                let radial::Curve::Circle { axis, .. } = c else {
                                    return Err(GeometryError::Invalid);
                                };
                                canonical_axis(axis)?;
                                (shell::Curve::Original(c), ParameterMapping::AngularIdentity)
                            }
                            _ => return Err(GeometryError::Unsupported),
                        };
                        input.curves.push(value);
                        maps.push(ParameterMap {
                            target: Target::Curve(i),
                            entity: curve_id,
                            mapping,
                        });
                        curves.insert(curve_id, i);
                        i
                    };
                    let i = input.edges.len();
                    input.edges.push(radial::Edge {
                        vertices: ends,
                        curve,
                        forward: fwd,
                    });
                    edges.insert(edge_id, i);
                    i
                };
                coedges.push((edge_index, direction));
            }
            if outside {
                if outer.replace(coedges).is_some() {
                    return Err(GeometryError::Invalid);
                }
            } else {
                inner.push(coedges);
            }
        }
        input.faces.push(shell::Face {
            surface: surface_index,
            outer: outer.ok_or(GeometryError::Invalid)?,
            inner,
            forward,
        });
    }
    Ok(Decoded {
        ownership: shell::Ownership::new(input, interrupt)?,
        maps,
        member,
    })
}
fn knot_data(mult: &Value<'_>, knots: &Value<'_>, spans: usize) -> GeometryResult<()> {
    let (Value::List(mult), Value::List(knots)) = (mult, knots) else {
        return Err(GeometryError::Invalid);
    };
    if spans == 0 || spans > 524_288 || mult.len() != spans + 1 || knots.len() != spans + 1 {
        return Err(GeometryError::Unsupported);
    }
    for i in 0..=spans {
        if mult[i] != Value::Number(if i == 0 || i == spans { 10. } else { 9. })
            || knots[i] != Value::Number(i as f64)
        {
            return Err(GeometryError::Unsupported);
        }
    }
    Ok(())
}
fn add_controls(used: &mut usize, n: usize) -> GeometryResult<()> {
    *used = used.checked_add(n).ok_or(GeometryError::Limit)?;
    if *used > 524_288 {
        Err(GeometryError::Limit)
    } else {
        Ok(())
    }
}
fn spline(
    s: &Shell<'_, '_>,
    id: u32,
    controls: &mut usize,
    interrupt: Option<&dyn Fn() -> bool>,
) -> GeometryResult<(radial::Curve, usize)> {
    let [
        Value::String(_),
        Value::Number(degree),
        poles,
        Value::Enum("UNSPECIFIED"),
        Value::Enum("F"),
        Value::Enum("F"),
        mult,
        knots,
        Value::Enum("UNSPECIFIED"),
    ] = s
        .args(id, "B_SPLINE_CURVE_WITH_KNOTS")
        .map_err(geometry_error)?
    else {
        return Err(GeometryError::Invalid);
    };
    let Value::List(poles) = poles else {
        return Err(GeometryError::Invalid);
    };
    if *degree != 9. || poles.len() < 10 || !(poles.len() - 1).is_multiple_of(9) {
        return Err(GeometryError::Unsupported);
    }
    let spans = (poles.len() - 1) / 9;
    knot_data(mult, knots, spans)?;
    add_controls(controls, poles.len())?;
    let mut points = Vec::new();
    points
        .try_reserve_exact(poles.len())
        .map_err(|_| GeometryError::Limit)?;
    for p in poles {
        checkpoint(interrupt)?;
        points.push(
            s.point(reference(p).map_err(geometry_error)?)
                .map_err(geometry_error)?,
        );
    }
    Ok((
        radial::Curve::Spline {
            spans,
            poles: points,
        },
        spans,
    ))
}
fn tensor(
    s: &Shell<'_, '_>,
    id: u32,
    controls: &mut usize,
    interrupt: Option<&dyn Fn() -> bool>,
) -> GeometryResult<(RetainedTensor, [usize; 2])> {
    let [
        Value::String(_),
        Value::Number(du),
        Value::Number(dv),
        rows,
        Value::Enum("UNSPECIFIED"),
        Value::Enum("F"),
        Value::Enum("F"),
        Value::Enum("F"),
        um,
        vm,
        uk,
        vk,
        Value::Enum("UNSPECIFIED"),
    ] = s
        .args(id, "B_SPLINE_SURFACE_WITH_KNOTS")
        .map_err(geometry_error)?
    else {
        return Err(GeometryError::Invalid);
    };
    let Value::List(rows) = rows else {
        return Err(GeometryError::Invalid);
    };
    if *du != 9. || *dv != 9. || rows.len() < 10 || !(rows.len() - 1).is_multiple_of(9) {
        return Err(GeometryError::Unsupported);
    }
    let Value::List(first) = &rows[0] else {
        return Err(GeometryError::Invalid);
    };
    if first.len() < 10 || !(first.len() - 1).is_multiple_of(9) {
        return Err(GeometryError::Unsupported);
    }
    let dimensions = [rows.len(), first.len()];
    let spans = dimensions.map(|n| (n - 1) / 9);
    knot_data(um, uk, spans[0])?;
    knot_data(vm, vk, spans[1])?;
    let count = dimensions[0]
        .checked_mul(dimensions[1])
        .ok_or(GeometryError::Limit)?;
    add_controls(controls, count)?;
    let mut points = Vec::new();
    points
        .try_reserve_exact(count)
        .map_err(|_| GeometryError::Limit)?;
    for row in rows {
        let Value::List(row) = row else {
            return Err(GeometryError::Invalid);
        };
        if row.len() != dimensions[1] {
            return Err(GeometryError::Invalid);
        }
        for p in row {
            checkpoint(interrupt)?;
            points.push(
                s.point(reference(p).map_err(geometry_error)?)
                    .map_err(geometry_error)?,
            );
        }
    }
    RetainedTensor::new(spans, points).map(|t| (t, spans))
}
fn line(
    s: &Shell<'_, '_>,
    id: u32,
    ends: [usize; 2],
    forward: bool,
    vertices: &[[f64; 3]],
) -> GeometryResult<(shell::Curve, ParameterMapping)> {
    let [Value::String(_), p, v] = s.args(id, "LINE").map_err(geometry_error)? else {
        return Err(GeometryError::Invalid);
    };
    let origin = s
        .point(reference(p).map_err(geometry_error)?)
        .map_err(geometry_error)?;
    let [Value::String(_), d, Value::Number(magnitude)] = s
        .args(reference(v).map_err(geometry_error)?, "VECTOR")
        .map_err(geometry_error)?
    else {
        return Err(GeometryError::Invalid);
    };
    let [Value::String(_), value] = s
        .args(reference(d).map_err(geometry_error)?, "DIRECTION")
        .map_err(geometry_error)?
    else {
        return Err(GeometryError::Invalid);
    };
    let delta = numbers(value).map_err(geometry_error)?;
    let ends = if forward { ends } else { [ends[1], ends[0]] };
    if *magnitude != 1. || origin != vertices[ends[0]] {
        return Err(GeometryError::Unsupported);
    }
    for k in 0..3 {
        if !exact_numbers::exact_difference(vertices[ends[1]][k], origin[k], delta[k])? {
            return Err(GeometryError::Unsupported);
        }
    }
    let length = delta
        .map(I::point)
        .into_iter()
        .fold(I::point(0.), |s, d| s + d.square())
        .sqrt()
        .positive()?;
    Ok((
        shell::Curve::Line(ends),
        ParameterMapping::Line {
            delta,
            length_bounds_mm: [length.lo, length.hi],
        },
    ))
}
fn next_numeric<'a>(text: &mut &'a str) -> Option<&'a str> {
    while let Some(token) = crate::step_read::literal_member_binding::token(text) {
        if let crate::step_read::literal_member_binding::Token::Literal(t) = token
            && t.as_bytes()
                .first()
                .is_some_and(|b| matches!(b, b'+' | b'-' | b'0'..=b'9'))
        {
            return Some(t);
        }
    }
    None
}
fn validate_numbers(
    entities: &BTreeMap<u32, Vec<part21::Entity<'_>>>,
    interrupt: Option<&dyn Fn() -> bool>,
) -> GeometryResult<()> {
    fn validate_value(
        value: &Value<'_>,
        text: &mut &str,
        integer: bool,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> GeometryResult<()> {
        checkpoint(interrupt)?;
        match value {
            Value::Number(n) => {
                let token = next_numeric(text).ok_or(GeometryError::Invalid)?;
                if integer && !token.bytes().all(|b| b.is_ascii_digit()) {
                    return Err(GeometryError::Invalid);
                }
                if !exact_numbers::exact_leaf(token, *n)? {
                    return Err(GeometryError::Unsupported);
                }
            }
            Value::List(v) => {
                for value in v {
                    validate_value(value, text, integer, interrupt)?;
                }
            }
            Value::Typed(e) => {
                for value in &e.args {
                    validate_value(value, text, integer, interrupt)?;
                }
            }
            _ => {}
        }
        Ok(())
    }
    for records in entities.values() {
        for e in records {
            checkpoint(interrupt)?;
            if e.name == "GEOMETRIC_REPRESENTATION_CONTEXT" {
                let mut text = e.source.ok_or(GeometryError::Invalid)?;
                if !matches!(e.args.as_slice(),[Value::Number(n)]if *n==3.)
                    || next_numeric(&mut text) != Some("3")
                {
                    return Err(GeometryError::Invalid);
                }
                continue;
            }
            if matches!(
                e.name,
                "CARTESIAN_POINT"
                    | "DIRECTION"
                    | "VECTOR"
                    | "CIRCLE"
                    | "CONICAL_SURFACE"
                    | "CYLINDRICAL_SURFACE"
                    | "B_SPLINE_CURVE_WITH_KNOTS"
                    | "B_SPLINE_SURFACE_WITH_KNOTS"
            ) {
                let mut text = e.source.ok_or(GeometryError::Invalid)?;
                for (i, value) in e.args.iter().enumerate() {
                    let integer = match e.name {
                        "B_SPLINE_CURVE_WITH_KNOTS" => matches!(i, 1 | 6),
                        "B_SPLINE_SURFACE_WITH_KNOTS" => matches!(i, 1 | 2 | 8 | 9),
                        _ => false,
                    };
                    validate_value(value, &mut text, integer, interrupt)?;
                }
                if next_numeric(&mut text).is_some() {
                    return Err(GeometryError::Invalid);
                }
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
