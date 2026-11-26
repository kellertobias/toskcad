//! A moved or instanced occurrence of retained degree-eleven material. The
//! original local records are replayed literally (only entity identifiers
//! change) inside one declared map, placed by one `MAPPED_ITEM`. Qualification
//! decodes that member's LOCAL records independently, encloses the file's
//! exact EXPRESS operator, and transports every unchanged correspondence
//! bound through it to the source occurrence.
use super::*;
use crate::step::periodic_phillips;
use crate::step_read::owned_insert as document;
use crate::tapered_thread::interval::I;
fn document_error(e: document::Error) -> Error {
    match e {
        document::Error::Invalid => Error::Invalid,
        document::Error::Limit => Error::Limit,
        document::Error::Cancelled => Error::Cancelled,
        document::Error::Unsupported => Error::Unsupported,
    }
}
/// `occurrence` = origin, then the unit images of the local x, y, z axes.
/// The written operator states these axes, so they must already be a
/// right-handed orthonormal basis to a few rounding units.
fn proper(occurrence: [[f64; 3]; 4]) -> Result<(), Error> {
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let [origin, u, v, n] = occurrence;
    let cross = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    const EPSILON: f64 = 1e-12;
    if occurrence.iter().flatten().any(|x| !x.is_finite())
        || origin.iter().any(|x| x.abs() > 1e7)
        || (dot(u, u) - 1.).abs() > EPSILON
        || (dot(v, v) - 1.).abs() > EPSILON
        || (dot(n, n) - 1.).abs() > EPSILON
        || dot(u, v).abs() > EPSILON
        || (0..3).any(|i| (cross[i] - n[i]).abs() > 4. * EPSILON)
    {
        return Err(Error::Invalid);
    }
    Ok(())
}
/// Fresh qualification of ACTUAL bytes holding one mapped member.
pub fn qualify_placed_actual_file(
    source: &Material,
    text: &str,
    occurrence: [[f64; 3]; 4],
    options: Options,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<QualifiedMaterial, Error> {
    checkpoint(interrupt)?;
    proper(occurrence)?;
    if !valid(options) {
        return Err(Error::Invalid);
    }
    if text.len() > options.maximum_output_bytes {
        return Err(Error::Limit);
    }
    let actual = document::CompoundDocument::parse(text, Default::default(), interrupt)
        .map_err(document_error)?;
    let [member] = actual.represented_members() else {
        return Err(Error::Unsupported);
    };
    if member.mapping_link().is_none() {
        return Err(Error::Unsupported);
    }
    let item = member.item_entity();
    let operator = actual
        .member_operator(item, Default::default(), interrupt)
        .map_err(document_error)?;
    // The decimal scale is enclosed, not assumed exact: the map's linear
    // columns are that scale times the normalized EXPRESS basis.
    let [lo, hi] = operator.scale;
    if operator.handedness != 1 || !(lo > 0. && hi.is_finite()) {
        return Err(Error::Unsupported);
    }
    let scale = I::new(lo, hi);
    let column = |axis: [[f64; 2]; 3]| {
        axis.map(|[a, b]| {
            let v = I::new(a, b) * scale;
            [v.lo, v.hi]
        })
    };
    let local = actual
        .local_source_document(item, Default::default(), interrupt)
        .map_err(document_error)?;
    let [x, y, z] = operator.axes.map(column);
    let mut qualified = qualification_through(
        source,
        local.as_str(),
        local.as_str().len(),
        0,
        occurrence,
        [operator.origin_mm, x, y, z],
        options,
        interrupt,
    )?;
    let extra = actual
        .retained_budget_bytes()
        .checked_add(local.query_budget_bytes())
        .ok_or(Error::Limit)?;
    qualified.verification.known_resident_peak_bytes = qualified
        .verification
        .known_resident_peak_bytes
        .checked_add(extra)
        .filter(|x| *x <= options.maximum_known_resident_bytes)
        .ok_or(Error::Limit)?;
    checkpoint(interrupt)?;
    Ok(qualified)
}
/// Replay the original retained local records inside one mapped item at
/// `occurrence`, then qualify the written bytes independently.
pub fn export_retained_placed(
    source: &Material,
    occurrence: [[f64; 3]; 4],
    options: Options,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<QualifiedStep, Error> {
    checkpoint(interrupt)?;
    proper(occurrence)?;
    if !valid(options) {
        return Err(Error::Invalid);
    }
    let input = source.original_input();
    let mut writer = crate::step::Writer::new();
    writer
        .append_qualified_singleton(
            input,
            source.retained_budget_bytes() + input.len(),
            periodic_phillips::AppendOptions::default(),
            Some(occurrence),
            interrupt,
        )
        .map_err(|e| match e {
            periodic_phillips::Error::Limit => Error::Limit,
            periodic_phillips::Error::Cancelled => Error::Cancelled,
            periodic_phillips::Error::Unsupported => Error::Unsupported,
            _ => Error::Invalid,
        })?;
    let text = writer.finish().map_err(|e| match e {
        crate::step::StepError::Limit => Error::Limit,
        crate::step::StepError::Unsupported => Error::Unsupported,
        _ => Error::Invalid,
    })?;
    if text.len() > options.maximum_output_bytes {
        return Err(Error::Limit);
    }
    let qualification = qualify_placed_actual_file(source, &text, occurrence, options, interrupt)?;
    Ok(QualifiedStep {
        text,
        qualification,
    })
}
