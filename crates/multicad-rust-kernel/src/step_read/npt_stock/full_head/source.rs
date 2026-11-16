//! Source matching follows, and never replaces, independent actual material.
use super::*;
use crate::{
    step::npt_stock::{CorrespondenceOptions, CorrespondenceReceipt, SurfacePlan},
    tapered_thread::interval::I,
};
#[derive(Clone, Copy, Debug)]
pub struct HexHeadIntent {
    pub across_flats_mm: f64,
    pub height_mm: f64,
}
#[derive(Clone, Copy, Debug)]
pub struct SourceOptions {
    pub resources: Options,
    pub position_error_mm: f64,
    pub original_thread_derivative_error_mm: f64,
    pub material_tolerance_mm: f64,
    pub minimum_wall_mm: f64,
}
impl Default for SourceOptions {
    fn default() -> Self {
        Self {
            resources: Options::default(),
            position_error_mm: 1e-7,
            original_thread_derivative_error_mm: 1e-6,
            material_tolerance_mm: 1e-7,
            minimum_wall_mm: 1e-6,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct SourceReceipt {
    pub material: Receipt,
    /// C0/D1 for ORIGINAL thread law versus RAW retained STEP coefficients.
    /// This does not claim D1 after axial sewing, head UV-D1 or world transport.
    pub original_thread: CorrespondenceReceipt,
    pub head_author_corner_error_mm: f64,
    pub declared_head_planar_sewing_mm: f64,
    pub sewn_reference_position_error_mm: f64,
    pub head_matching_work: usize,
    pub aggregate_work: usize,
    pub accounted_peak_bytes: usize,
}
pub struct QualifiedFile<'input, 'plan, 'source> {
    actual: &'input str,
    source: &'plan SurfacePlan<'source>,
    receipt: SourceReceipt,
}
impl<'input, 'plan, 'source> QualifiedFile<'input, 'plan, 'source> {
    pub fn as_str(&self) -> &'input str {
        self.actual
    }
    pub fn source_plan(&self) -> &'plan SurfacePlan<'source> {
        self.source
    }
    pub fn receipt(&self) -> SourceReceipt {
        self.receipt
    }
}
fn source_error(e: crate::step::npt_stock::Error) -> Error {
    match e {
        crate::step::npt_stock::Error::Invalid => Error::Invalid,
        crate::step::npt_stock::Error::Unresolved => Error::Unsupported,
        crate::step::npt_stock::Error::Limit => Error::Limit,
        crate::step::npt_stock::Error::Cancelled => Error::Cancelled,
    }
}
fn head_error(
    material: &OccupiedClippedHead<'_, '_, '_>,
    source: &SurfacePlan<'_>,
    intent: HexHeadIntent,
    tolerance: f64,
    b: &mut Budget<'_>,
) -> Result<f64> {
    b.step(256)?;
    // Author definitions are interval expressions, not rounded corner caches.
    let across = I::point(intent.across_flats_mm);
    let radius = across / I::point(3.).sqrt();
    let half = across / I::point(2.);
    let r2 = radius / I::point(2.);
    let zero = I::point(0.);
    let xy = [
        [radius, zero],
        [r2, half],
        [-r2, half],
        [-radius, zero],
        [-r2, -half],
        [r2, -half],
    ];
    let low = I::point(source.curves().source().source().spec().length_mm);
    let high = low + I::point(intent.height_mm);
    let expected = std::array::from_fn::<_, 12, _>(|i| {
        [xy[i % 6][0], xy[i % 6][1], if i < 6 { low } else { high }]
    });
    let mut assigned = [false; 12];
    let mut mapping = [usize::MAX; 12];
    let mut error = 0f64;
    for (index, actual) in material.head().vertices.iter().enumerate() {
        let mut matched = None;
        for (i, p) in expected.iter().enumerate() {
            b.step(256)?;
            let delta = (0..3).fold(I::point(0.), |sum, k| {
                sum + I::point((I::point(actual[k]) - p[k]).magnitude())
            });
            if !delta.valid() {
                return Err(Error::Unsupported);
            }
            if delta.hi <= tolerance {
                if matched.is_some() {
                    return Err(Error::Unsupported);
                }
                matched = Some((i, delta.hi));
            }
        }
        let (i, d) = matched.ok_or(Error::Unsupported)?;
        if assigned[i] {
            return Err(Error::Unsupported);
        }
        assigned[i] = true;
        mapping[index] = i;
        error = error.max(d);
    }
    // Global convexity has already been independently proved. Explicit cyclic
    // corner correspondence gives the piecewise-affine whole-head C0 map.
    let mut upper = [0usize; 6];
    for (i, p) in material.head_polygon_mm().iter().enumerate() {
        b.step(128)?;
        let id = material
            .head()
            .vertices
            .iter()
            .position(|v| *v == [p[0], p[1], material.head_stations_mm()[1]])
            .ok_or(Error::Unsupported)?;
        upper[i] = mapping[id];
        if upper[i] < 6 {
            return Err(Error::Unsupported);
        }
    }
    if !(0..6).all(|i| upper[(i + 1) % 6] - 6 == (upper[i] - 6 + 1) % 6) {
        return Err(Error::Unsupported);
    }
    Ok(error)
}
/// Qualify exact supplied bytes against an ORIGINAL clipped thread plan and
/// a canonical fixed-frame hex-head author definition. Independent whole-head
/// occupied admission always runs first. No Native/export flag is granted.
pub fn qualify_actual_file<'input, 'plan, 'source>(
    actual: &'input str,
    source: &'plan SurfacePlan<'source>,
    head: HexHeadIntent,
    options: SourceOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<QualifiedFile<'input, 'plan, 'source>> {
    let r = options.resources;
    if r.maximum_work == 0
        || r.maximum_work > MAX_WORK
        || r.maximum_bytes == 0
        || r.maximum_bytes > MAX_BYTES
        || !options.position_error_mm.is_finite()
        || options.position_error_mm <= 0.
        || !options.original_thread_derivative_error_mm.is_finite()
        || options.original_thread_derivative_error_mm <= 0.
        || !options.material_tolerance_mm.is_finite()
        || options.material_tolerance_mm <= 0.
        || options.material_tolerance_mm > 1e-4
        || !options.minimum_wall_mm.is_finite()
        || options.minimum_wall_mm <= 0.
        || !head.across_flats_mm.is_finite()
        || head.across_flats_mm <= 0.
        || !head.height_mm.is_finite()
        || head.height_mm <= 0.
    {
        return Err(Error::Invalid);
    }
    let source_work = source.work_used();
    let source_bytes = source.retained_budget_bytes();
    let maximum_work = r
        .maximum_work
        .checked_sub(source_work)
        .filter(|n| *n > 0)
        .ok_or(Error::Limit)?;
    let additional = r
        .additional_resident_bytes
        .checked_add(source_bytes)
        .ok_or(Error::Limit)?;
    let packet = DecodedFile::parse(
        actual,
        Options {
            maximum_work,
            additional_resident_bytes: additional,
            ..r
        },
        interrupt,
    )?;
    // Packet already reserved this exact second cold parse, including source
    // plan residency. Its logical document prefix is not debited a second time.
    let doc = CompoundDocument::parse(
        packet.derived_auxiliary_text(),
        DocumentOptions {
            maximum_work: packet.derived_parser_work_reserved,
            retained_budget_bytes: r
                .maximum_bytes
                .checked_sub(packet.retained_beyond_auxiliary_input)
                .filter(|n| *n > 0)
                .ok_or(Error::Limit)?,
        },
        interrupt,
    )
    .map_err(neutral_error)?;
    let parse_bytes = r
        .maximum_bytes
        .checked_sub(packet.retained_beyond_auxiliary_input)
        .and_then(|n| n.checked_sub(doc.retained_budget_bytes()))
        .filter(|n| *n > 0)
        .ok_or(Error::Limit)?;
    let parse_work = maximum_work
        .checked_sub(packet.receipt().aggregate_work)
        .filter(|n| *n > 0)
        .ok_or(Error::Limit)?;
    let parsed = Parsed::parse(
        &doc,
        doc.represented_members()[0].item_entity(),
        super::super::Options {
            maximum_work: parse_work.min(16_000_000),
            maximum_bytes: parse_bytes,
        },
        interrupt,
    )?;
    let occupied = packet.certify(
        &parsed,
        options.material_tolerance_mm,
        options.minimum_wall_mm,
        interrupt,
    )?;
    let old_head_bytes = head::storage(occupied.material().head())?;
    let outside = doc
        .retained_budget_bytes()
        .checked_add(parsed.query_budget_bytes())
        .and_then(|n| n.checked_add(packet.retained_beyond_auxiliary_input))
        .and_then(|n| n.checked_add(old_head_bytes))
        .ok_or(Error::Limit)?;
    let matching_resident = outside.checked_add(16 * 1024).ok_or(Error::Limit)?;
    if matching_resident > r.maximum_bytes {
        return Err(Error::Limit);
    }
    let remaining = r
        .maximum_work
        .checked_sub(source_work)
        .and_then(|n| n.checked_sub(occupied.receipt().aggregate_work))
        .filter(|n| *n > 0)
        .ok_or(Error::Limit)?;
    let mut b = Budget {
        used: 0,
        remaining: remaining.min(16_000_000),
        bytes: matching_resident,
        peak: matching_resident,
        limit: r.maximum_bytes,
        interrupt,
    };
    let corner = head_error(
        occupied.material(),
        source,
        head,
        options.position_error_mm,
        &mut b,
    )?;
    let head_error = (I::point(corner) + I::point(occupied.material().planar_sewing_mm())).hi;
    if !head_error.is_finite() || head_error > options.position_error_mm {
        return Err(Error::Unsupported);
    }
    let matching_work = b.used;
    let matching_peak = b.peak;
    let available = r
        .maximum_work
        .checked_sub(source_work)
        .and_then(|n| n.checked_sub(occupied.receipt().aggregate_work))
        .and_then(|n| n.checked_sub(matching_work))
        .filter(|n| *n > 0)
        .ok_or(Error::Limit)?;
    let child_bytes = r
        .maximum_bytes
        .checked_sub(outside)
        .filter(|n| *n > 0)
        .ok_or(Error::Limit)?;
    let selected = Selected::parse(
        &doc,
        doc.represented_members()[0].item_entity(),
        super::super::Options {
            maximum_work: available.min(16_000_000),
            maximum_bytes: child_bytes,
        },
        interrupt,
    )?;
    let retained = selected.retain(
        MaterialOptions {
            maximum_work: available.min(64_000_000),
            maximum_bytes: child_bytes,
            ..Default::default()
        },
        interrupt,
    )?;
    let retained_peak = outside
        .checked_add(retained.construction_query_peak_bytes())
        .ok_or(Error::Limit)?;
    if retained_peak > r.maximum_bytes {
        return Err(Error::Limit);
    }
    drop(selected);
    // Correspondence includes source+owned prefixes once. The outside scope
    // already includes source residency, so remove only that known shared plan.
    let compare_outside = outside.checked_sub(source_bytes).ok_or(Error::Limit)?;
    let compare_bytes = r
        .maximum_bytes
        .checked_sub(compare_outside)
        .filter(|n| *n > 0)
        .ok_or(Error::Limit)?;
    let compare_work = r
        .maximum_work
        .checked_sub(occupied.receipt().aggregate_work)
        .and_then(|n| n.checked_sub(matching_work))
        .filter(|n| *n > 0)
        .ok_or(Error::Limit)?;
    let raw_position = options.position_error_mm - occupied.material().join_sewing_mm();
    if !raw_position.is_finite() || raw_position <= 0. {
        return Err(Error::Unsupported);
    }
    let thread = source
        .correspond_literal(
            &retained,
            CorrespondenceOptions {
                position_error_mm: raw_position,
                original_derivative_error_mm: options.original_thread_derivative_error_mm,
                maximum_work: compare_work.min(96_000_000),
                maximum_bytes: compare_bytes,
            },
            interrupt,
        )
        .map_err(source_error)?;
    let sewn_thread =
        (I::point(thread.position_error_mm()) + I::point(occupied.material().join_sewing_mm())).hi;
    let sewn = sewn_thread.max(head_error);
    if !sewn.is_finite() || sewn > options.position_error_mm {
        return Err(Error::Unsupported);
    }
    let aggregate = occupied
        .receipt()
        .aggregate_work
        .checked_add(matching_work)
        .and_then(|n| n.checked_add(thread.aggregate_work_used()))
        .ok_or(Error::Limit)?;
    if aggregate > r.maximum_work {
        return Err(Error::Limit);
    }
    let peak = occupied
        .receipt()
        .accounted_peak_bytes
        .max(matching_peak)
        .max(retained_peak)
        .max(
            compare_outside
                .checked_add(thread.accounted_peak_bytes())
                .ok_or(Error::Limit)?,
        );
    if peak > r.maximum_bytes {
        return Err(Error::Limit);
    }
    if interrupt.is_some_and(|f| f()) {
        return Err(Error::Cancelled);
    }
    Ok(QualifiedFile {
        actual,
        source,
        receipt: SourceReceipt {
            material: occupied.receipt(),
            original_thread: thread,
            head_author_corner_error_mm: corner,
            declared_head_planar_sewing_mm: occupied.material().planar_sewing_mm(),
            sewn_reference_position_error_mm: sewn,
            head_matching_work: matching_work,
            aggregate_work: aggregate,
            accounted_peak_bytes: peak,
        },
    })
}
