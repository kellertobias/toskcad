//! Re-export of SAME actual imported local geometry under actual B, then A.
//! Exact local records/topology/UV transfer independently; only directed world
//! C0 is asserted for the output operator. No authored law, mass equality or
//! sewn/chart D1 claim is inferred from this position comparison.
use super::*;
use crate::step_read::npt_stock::LiteralOccurrencePlacement;
#[cfg(test)]
mod tests;
const COMPARE_WORK: usize = 2048;
const PROPOSAL_WORK: usize = 2048;
#[derive(Clone, Copy, Debug)]
pub struct ReplayOptions {
    pub resources: Options,
    pub maximum_position_error_mm: f64,
}
impl Default for ReplayOptions {
    fn default() -> Self {
        Self {
            resources: Default::default(),
            maximum_position_error_mm: 1e-7,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct ReplayReceipt {
    pub inherited_source_work: usize,
    pub inherited_occurrence_work: usize,
    pub fresh_encoding_work: usize,
    pub fresh_output_binding_work: usize,
    pub operator_comparison_work: usize,
    pub aggregate_work: usize,
    pub accounted_peak_bytes: usize,
    /// Complete A(B(globalT(auxiliary) union actual head)) enclosure comparison.
    /// Local source/UV/IDs are exact; this is not raw mass or a D1 assertion.
    pub maximum_world_position_error_mm: f64,
}
pub struct VerifiedReplay<'source, 'pose> {
    actual: MappedActualHead,
    source: &'source MappedActualHead,
    occurrence: Option<&'pose LiteralOccurrencePlacement>,
    receipt: ReplayReceipt,
}
impl<'source, 'pose> VerifiedReplay<'source, 'pose> {
    pub fn actual(&self) -> &MappedActualHead {
        &self.actual
    }
    pub fn original(&self) -> &'source MappedActualHead {
        self.source
    }
    pub fn occurrence(&self) -> Option<&'pose LiteralOccurrencePlacement> {
        self.occurrence
    }
    pub fn as_str(&self) -> &str {
        self.actual.actual_file()
    }
    pub fn receipt(&self) -> ReplayReceipt {
        self.receipt
    }
}
fn validate(options: ReplayOptions, interrupt: Option<&dyn Fn() -> bool>) -> Result<()> {
    if interrupt.is_some_and(|f| f()) {
        return Err(Error::Cancelled);
    }
    let r = options.resources;
    if r.maximum_work == 0
        || r.maximum_work > MAX_WORK
        || r.maximum_bytes == 0
        || r.maximum_bytes > MAX_BYTES
        || !options.maximum_position_error_mm.is_finite()
        || options.maximum_position_error_mm <= 0.
        || options.maximum_position_error_mm > 1e-4
    {
        return Err(Error::Invalid);
    }
    Ok(())
}
/// B is the actual independently parsed file operator, NEVER a guessed
/// intrinsic author factor. A is the separately retained complete Native chain.
fn operator(
    source: &MappedActualHead,
    a: Option<&LiteralOccurrencePlacement>,
) -> ([[I; 3]; 3], [I; 3]) {
    let f = source.actual_frame();
    let scale = I::new(f.scale[0], f.scale[1]);
    let b = std::array::from_fn::<_, 3, _>(|i| {
        std::array::from_fn::<_, 3, _>(|j| scale * I::new(f.axes[j][i][0], f.axes[j][i][1]))
    });
    let origin = f.origin.map(|x| I::new(x[0], x[1]));
    match a {
        None => (b, origin),
        Some(a) => {
            let rows = a
                .linear_rows_bounds()
                .map(|r| r.map(|x| I::new(x[0], x[1])));
            let offset = a.translation_bounds_mm().map(|x| I::new(x[0], x[1]));
            let out = std::array::from_fn(|i| {
                std::array::from_fn(|j| (0..3).fold(I::point(0.), |v, k| v + rows[i][k] * b[k][j]))
            });
            let origin =
                std::array::from_fn(|i| (0..3).fold(offset[i], |v, k| v + rows[i][k] * origin[k]));
            (out, origin)
        }
    }
}
fn pose_work(a: Option<&LiteralOccurrencePlacement>) -> usize {
    a.map_or(0, LiteralOccurrencePlacement::work_used)
}
fn pose_bytes(a: Option<&LiteralOccurrencePlacement>) -> usize {
    a.map_or(0, LiteralOccurrencePlacement::retained_budget_bytes)
}
fn outside(
    source: &MappedActualHead,
    a: Option<&LiteralOccurrencePlacement>,
    r: Options,
) -> Result<usize> {
    r.additional_resident_bytes
        .checked_add(source.receipt.retained_bytes)
        .and_then(|n| n.checked_sub(source.local.receipt().retained_bytes))
        .and_then(|n| n.checked_add(pose_bytes(a)))
        .and_then(|n| n.checked_add(size_of::<VerifiedReplay<'_, '_>>()))
        .ok_or(Error::Limit)
}
/// Independently parse and bind an ACTUAL supplied output. Existing full local
/// regularity, finite graph, occupied material and topology transfer ONLY by
/// exact all-record/ID/UV equality. Its positive uniform proper output map is
/// freshly earned, preserving embedding by an invertible map; C0 alone never
/// authorizes topology. Cold public import still earns material independently.
pub fn qualify<'source, 'pose>(
    source: &'source MappedActualHead,
    a: Option<&'pose LiteralOccurrencePlacement>,
    input: Arc<str>,
    options: ReplayOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<VerifiedReplay<'source, 'pose>> {
    validate(options, interrupt)?;
    let r = options.resources;
    let inherited = source
        .receipt
        .aggregate_work
        .checked_add(pose_work(a))
        .ok_or(Error::Limit)?;
    let before = inherited
        .checked_add(COMPARE_WORK)
        .filter(|n| *n <= r.maximum_work)
        .ok_or(Error::Limit)?;
    let extra_work = before
        .checked_sub(source.local.receipt().aggregate_work)
        .ok_or(Error::Limit)?;
    let external = outside(source, a, r)?;
    let actual = MappedActualHead::read_reusing_local(
        input,
        Arc::clone(&source.local),
        OwnedHeadOptions {
            resources: Options {
                maximum_work: r.maximum_work.checked_sub(extra_work).ok_or(Error::Limit)?,
                additional_resident_bytes: external,
                ..r
            },
            material_tolerance_mm: source.local.material_tolerance_mm(),
            minimum_wall_mm: source.local.minimum_wall_mm(),
        },
        interrupt,
    )?;
    if interrupt.is_some_and(|f| f()) {
        return Err(Error::Cancelled);
    }
    // COMPARE_WORK is precharged BEFORE both source/decoded matrix products
    // and complete-box evaluation. There is no sampled-point approximation.
    let (target, origin) = operator(source, a);
    let (emitted, offset) = operator(&actual, None);
    let bound = source
        .local
        .local_sewn_enclosure_mm()
        .map(|x| I::new(x[0], x[1]));
    let delta = std::array::from_fn::<_, 3, _>(|i| {
        (0..3).fold(origin[i] - offset[i], |v, j| {
            v + (target[i][j] - emitted[i][j]) * bound[j]
        })
    });
    if delta.iter().any(|x| !x.valid()) {
        return Err(Error::Unsupported);
    }
    let error = delta
        .into_iter()
        .fold(I::point(0.), |v, x| v + I::point(x.magnitude()))
        .hi;
    if !error.is_finite() || error > options.maximum_position_error_mm {
        return Err(Error::Unsupported);
    }
    let fresh = actual
        .receipt
        .aggregate_work
        .checked_sub(source.local.receipt().aggregate_work)
        .ok_or(Error::Limit)?;
    let aggregate = before
        .checked_add(fresh)
        .filter(|n| *n <= r.maximum_work)
        .ok_or(Error::Limit)?;
    // Only this shared local payload is resident once; original mapped owner
    // (including original bytes), output bytes and all A buffers coexist.
    let inherited_peak = source
        .receipt
        .accounted_peak_bytes
        .checked_sub(source.receipt.additional_resident_bytes)
        .and_then(|n| n.checked_add(pose_bytes(a)))
        .and_then(|n| n.checked_add(r.additional_resident_bytes))
        .and_then(|n| n.checked_add(size_of::<VerifiedReplay<'_, '_>>()))
        .ok_or(Error::Limit)?;
    let pose_peak = match a {
        Some(a) => a
            .accounted_peak_bytes()
            .checked_add(source.receipt.retained_bytes)
            .and_then(|n| n.checked_add(r.additional_resident_bytes))
            .and_then(|n| n.checked_add(size_of::<VerifiedReplay<'_, '_>>()))
            .ok_or(Error::Limit)?,
        None => 0,
    };
    let peak = actual
        .receipt
        .accounted_peak_bytes
        .max(inherited_peak)
        .max(pose_peak);
    if peak > r.maximum_bytes {
        return Err(Error::Limit);
    }
    if interrupt.is_some_and(|f| f()) {
        return Err(Error::Cancelled);
    }
    Ok(VerifiedReplay {
        actual,
        source,
        occurrence: a,
        receipt: ReplayReceipt {
            inherited_source_work: source.receipt.aggregate_work,
            inherited_occurrence_work: pose_work(a),
            fresh_encoding_work: 0,
            fresh_output_binding_work: fresh,
            operator_comparison_work: COMPARE_WORK,
            aggregate_work: aggregate,
            accounted_peak_bytes: peak,
            maximum_world_position_error_mm: error,
        },
    })
}
/// NON-admitting output proposal. No STEP/material capability is conferred.
#[derive(Clone, Copy, Debug)]
pub struct PreparationReceipt {
    pub inherited_source_work: usize,
    pub inherited_occurrence_work: usize,
    pub fresh_encoding_work: usize,
    pub proposal_work: usize,
    pub aggregate_work: usize,
    pub accounted_peak_bytes: usize,
}
pub struct ReplayCandidate {
    bytes: Arc<str>,
    receipt: PreparationReceipt,
}
impl ReplayCandidate {
    pub fn as_str(&self) -> &str {
        &self.bytes
    }
    pub fn receipt(&self) -> PreparationReceipt {
        self.receipt
    }
    pub fn into_bytes(self) -> Arc<str> {
        self.bytes
    }
}
/// Produce ONLY a numerical proposal, for staging then actual file decoding.
/// Source/occurrence prefixes are counted once; no local proof is fabricated.
pub fn propose(
    source: &MappedActualHead,
    a: Option<&LiteralOccurrencePlacement>,
    options: ReplayOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<ReplayCandidate> {
    validate(options, interrupt)?;
    let r = options.resources;
    let prefix = source
        .receipt
        .aggregate_work
        .checked_add(pose_work(a))
        .and_then(|n| n.checked_add(PROPOSAL_WORK))
        .filter(|n| *n <= r.maximum_work)
        .ok_or(Error::Limit)?;
    let extra = prefix
        .checked_sub(source.local.receipt().aggregate_work)
        .ok_or(Error::Limit)?;
    let external = outside(source, a, r)?;
    external
        .checked_add(source.local.receipt().retained_bytes)
        .and_then(|n| n.checked_add(16 * 1024))
        .filter(|n| *n <= r.maximum_bytes)
        .ok_or(Error::Limit)?;
    let (rows, origin) = operator(source, a);
    let candidate = encoding::propose_operator(
        &source.local,
        encoding::OperatorProposal {
            rows,
            origin,
            work: 0,
            retained: 0,
            peak: 0,
        },
        encoding::EncodingOptions {
            maximum_work: r.maximum_work.checked_sub(extra).ok_or(Error::Limit)?,
            maximum_bytes: r.maximum_bytes,
            additional_resident_bytes: external,
        },
        interrupt,
    )?;
    let encoded = candidate.receipt();
    let aggregate = prefix
        .checked_add(encoded.fresh_encoding_work)
        .filter(|n| *n <= r.maximum_work)
        .ok_or(Error::Limit)?;
    Ok(ReplayCandidate {
        bytes: candidate.actual_bytes(),
        receipt: PreparationReceipt {
            inherited_source_work: source.receipt.aggregate_work,
            inherited_occurrence_work: pose_work(a),
            fresh_encoding_work: encoded.fresh_encoding_work,
            proposal_work: PROPOSAL_WORK,
            aggregate_work: aggregate,
            accounted_peak_bytes: encoded.accounted_peak_bytes,
        },
    })
}
/// In-memory convenience: candidate followed by independent output decoding.
pub fn encode<'source, 'pose>(
    source: &'source MappedActualHead,
    a: Option<&'pose LiteralOccurrencePlacement>,
    options: ReplayOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<VerifiedReplay<'source, 'pose>> {
    let candidate = propose(source, a, options, interrupt)?;
    let prepared = candidate.receipt();
    let bytes = candidate.into_bytes();
    let removed = prepared
        .fresh_encoding_work
        .checked_add(prepared.proposal_work)
        .ok_or(Error::Limit)?;
    let mut verified = qualify(
        source,
        a,
        bytes,
        ReplayOptions {
            resources: Options {
                maximum_work: options
                    .resources
                    .maximum_work
                    .checked_sub(removed)
                    .filter(|n| *n > 0)
                    .ok_or(Error::Limit)?,
                ..options.resources
            },
            ..options
        },
        interrupt,
    )?;
    verified.receipt.fresh_encoding_work = prepared.fresh_encoding_work;
    verified.receipt.operator_comparison_work = verified
        .receipt
        .operator_comparison_work
        .checked_add(prepared.proposal_work)
        .ok_or(Error::Limit)?;
    verified.receipt.aggregate_work = verified
        .receipt
        .aggregate_work
        .checked_add(removed)
        .filter(|n| *n <= options.resources.maximum_work)
        .ok_or(Error::Limit)?;
    verified.receipt.accounted_peak_bytes = verified
        .receipt
        .accounted_peak_bytes
        .max(prepared.accounted_peak_bytes);
    if interrupt.is_some_and(|f| f()) {
        return Err(Error::Cancelled);
    }
    Ok(verified)
}
