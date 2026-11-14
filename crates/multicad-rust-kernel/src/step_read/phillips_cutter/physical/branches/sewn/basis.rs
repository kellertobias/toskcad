//! Exact discrete basis metadata; f64 equality is not INTEGER authority.
use super::*;
fn integer(spelling: &str, expected: u32) -> bool {
    let digits = spelling.strip_prefix('+').unwrap_or(spelling);
    !digits.is_empty()
        && digits.bytes().all(|c| c.is_ascii_digit())
        && digits.parse::<u32>() == Ok(expected)
}
pub(super) fn prove(
    body: &Body<'_>,
    entity: u32,
    count: usize,
    work: &mut Budget<'_>,
) -> BranchResult<()> {
    work.charge(1)?;
    let original = body.original_numeric_leaves(entity, work.remaining, work.interrupt)?;
    work.charge(original.work_used())?;
    footprint(
        original.retained_capacity_bytes(),
        0,
        0,
        QUERY_SCRATCH,
        work,
    )?;
    let leaves = original.leaves();
    // The sealed same-Body physical parser proved a simple non-rational cubic
    // B_SPLINE_CURVE_WITH_KNOTS with exactly 4/3 endpoint/interior multiplicity
    // counts. Here the ORIGINAL tokens, not their rounded f64s, earn those
    // discrete values. Ref IDs and quoted labels are not numeric leaves.
    // Its exact nine-field schema contains one numeric degree, two equally
    // sized numeric knot/multiplicity lists, and a reference-only pole list.
    // The same-Body collector preserves exactly those numeric leaves, so the
    // count is already 1+2*(span_count+1); it is not an independent fallible
    // condition or permission to discard any original coefficient.
    let bytes = leaves[..count + 2]
        .iter()
        .map(|v| v.spelling().len())
        .fold(0usize, usize::saturating_add);
    work.charge(bytes.saturating_mul(2))?;
    if !integer(leaves[0].spelling(), 3) {
        return Err(PhillipsReadError::Unsupported);
    }
    for i in 0..=count {
        work.charge(1)?;
        if !integer(
            leaves[i + 1].spelling(),
            if i == 0 || i == count { 4 } else { 3 },
        ) {
            return Err(PhillipsReadError::Unsupported);
        }
    }
    // Strict decoded knot order implies strict original numeric order under
    // the validated decimal->f64 rounding contract. Unit-span Bezier shape is
    // independent of the actual retained knot widths; no global derivative
    // parameter or substituted knot coefficient is claimed by this receipt.
    Ok(())
}
