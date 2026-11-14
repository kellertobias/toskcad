//! Measured correspondence to the ORIGINAL finite indexed/evaluator paths.
//! No source rate, station or support is normalized or replaced.
use super::*;
#[derive(Clone, Copy)]
pub(super) struct Allowance {
    pub indexed_mm: f64,
    pub evaluator_mm: f64,
    pub total_mm: f64,
}
fn ulp_bound(magnitude: I) -> R<I> {
    let bound = magnitude.magnitude();
    if !bound.is_finite() {
        return Err(Error::Unsupported);
    }
    Ok(ip(finite(ip(bound.next_up()) - ip(bound))?.hi))
}
pub(super) fn measure(thread: &Thread, interrupt: &dyn Fn() -> bool) -> R<Allowance> {
    checkpoint(interrupt)?;
    let spec = thread.spec();
    let pitch = ip(spec.pitch_mm);
    let mut indexed = 0_f64;
    // Same finite and neighbor-period domain as the original classifier and
    // PhasePrecision property receipt. Thread admission bounds this to <=265
    // periods; each has six source joins and both literal evaluation paths.
    let margin = 2 * i64::from(spec.starts) + 4;
    let last = (spec.length_mm / spec.pitch_mm).ceil() as i64 + margin;
    for k in -margin..=last {
        checkpoint(interrupt)?;
        for knot in thread.profile_breaks() {
            let exact = ip(k as f64) * pitch + ip(knot);
            for literal in [
                (k as f64).mul_add(spec.pitch_mm, knot),
                k as f64 * spec.pitch_mm + knot,
            ] {
                indexed = indexed.max((ip(literal) - exact).magnitude());
            }
        }
    }
    let true_tau = ip(2.)
        * I::new(
            std::f64::consts::PI.next_down(),
            std::f64::consts::PI.next_up(),
        );
    let literal_tau = ip(std::f64::consts::TAU);
    let rate = ip(thread.helical_rate().abs());
    let lead = pitch * ip(f64::from(spec.starts));
    let defect = (rate * true_tau - lead)
        .magnitude()
        .max((rate * literal_tau - lead).magnitude());
    let turns = (ip(spec.length_mm) / lead).hi.ceil() + 2.;
    // Rust rem_euclid uses x%rhs followed by an optional positive-modulus add:
    // https://doc.rust-lang.org/src/core/num/f64.rs.html#2181-2184
    // LLVM frem has libm fmod semantics under the default FP environment:
    // https://llvm.org/docs/LangRef.html#frem-instruction
    // The remainder is rounded to its result format; rem_euclid may also add
    // the positive modulus. Charge a full maximum ULP per operation: one for
    // nonnegative z remainder, two for angle remainder, two for final phase
    // remainder, and one each for multiply/subtract. This also covers ties,
    // binade transitions and subnormal results without a relative-error model.
    let pitch_ulp = ulp_bound(pitch)?;
    let angle_ulp = ulp_bound(literal_tau)?;
    // The existing source classifier separately encloses atan2 evaluation by
    // 16eps radians; retain that same allowance rather than assuming exact atan.
    let angle_error = ip(16. * f64::EPSILON) + ip(2.) * angle_ulp;
    let product = rate * (true_tau + angle_error);
    let arithmetic =
        ip(3.) * pitch_ulp + rate * angle_error + ulp_bound(product)? + ulp_bound(pitch + product)?;
    // Admitted ranges give |k*p+break|<2e9, lead<=8e7, rate<1.4e7,
    // turns<=130. Full ULPs on these ranges, including directed I operations,
    // put the complete allowance below .001mm; this proves caller finiteness
    // independently of a requested manufacturing/clearance tolerance.
    let evaluator = finite(arithmetic + ip(turns) * ip(defect))?;
    let total = finite(ip(indexed) + evaluator)?;
    checkpoint(interrupt)?;
    Ok(Allowance {
        indexed_mm: indexed,
        evaluator_mm: evaluator.hi,
        total_mm: total.hi,
    })
}
#[cfg(test)]
mod tests;
