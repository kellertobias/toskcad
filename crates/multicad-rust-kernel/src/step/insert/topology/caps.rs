//! Actual cap rims: strict angular regularity plus earned one-turn winding.
use super::*;
use crate::step::insert::projection::polynomial::P;
#[cfg(test)]
mod tests;
const DEPTH: usize = 32;
const SCRATCH: usize = 64 * 1024;
#[derive(Clone, Copy, Debug)]
pub struct CapRegularity {
    pub winding: [[i32; 2]; 2],
    pub minimum_signed_angular_derivative: f64,
    pub certified_segments: usize,
    pub work_used: usize,
}
impl TopologyPlan<'_, '_, '_> {
    /// Each actual rim has one strict angular direction on EVERY complete
    /// polynomial span. Together with its independently earned winding +/-1,
    /// this proves a simple radial Jordan rim; winding alone does not do so.
    pub fn cap_regularity_bounds(
        &self,
        maximum_work: usize,
        maximum_bytes: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<CapRegularity, Error> {
        if maximum_bytes == 0 || maximum_bytes > 64 * 1024 * 1024 {
            return Err(Error::Invalid);
        }
        let mut memory = Memory::new(self.retained_bytes(), maximum_bytes)?;
        memory.add(SCRATCH)?;
        // Parent fixed polynomial scratch remains live during the child query.
        let winding = self.winding_bounds(maximum_work, maximum_bytes - SCRATCH, interrupt)?;
        let mut w = Work::new(maximum_work, interrupt)?;
        w.charge(winding.work_used)?;
        let mut stack = memory.allocate(DEPTH + 1, &mut w)?;
        let mut minimum = f64::INFINITY;
        let mut count = 0;
        for face in self.faces() {
            let FaceGeometry::Cap(end) = face.geometry else {
                continue;
            };
            let height = if end == 0 {
                0.
            } else {
                self.surface().domain().boundary().source().full_length_mm()
            };
            for (index, wire) in face.loops.iter().enumerate() {
                for c in wire {
                    // The preceding winding query independently earned that
                    // every cap coedge owns a row; this is not a new guard.
                    let row = self.edges()[c.edge]
                        .geometry
                        .row_index()
                        .expect("cap row ownership earned by winding");
                    let row = &self.surface().rows()[row];
                    coplanar(row.poles(), height, &mut w)?;
                    let sign = f64::from(winding.winding[end as usize][index])
                        * if c.forward { 1. } else { -1. };
                    for poles in row.poles().windows(12).step_by(11) {
                        let x =
                            P::controls(&std::array::from_fn::<_, 12, _>(|i| literal(poles[i][0])));
                        let y =
                            P::controls(&std::array::from_fn::<_, 12, _>(|i| literal(poles[i][1])));
                        let p = x
                            .product(y.derivative(&mut w)?, &mut w)?
                            .sub(y.product(x.derivative(&mut w)?, &mut w)?, &mut w)?
                            .scale(sign, &mut w)?;
                        minimum = minimum.min(positive(p, &mut stack, &mut w)?);
                        count += 1;
                    }
                }
            }
        }
        w.charge(1)?;
        Ok(CapRegularity {
            winding: winding.winding,
            minimum_signed_angular_derivative: minimum,
            certified_segments: count,
            work_used: w.used,
        })
    }
}
fn literal(x: f64) -> I {
    if x == 0. || (x.fract() == 0. && x.abs() <= 9_007_199_254_740_992.) {
        I::point(x)
    } else {
        I::new(x.next_down(), x.next_up())
    }
}
fn coplanar(poles: &[[f64; 3]], height: f64, w: &mut Work<'_>) -> Result<(), Error> {
    // This is the OWNED writer's identical binary64 -> identical decimal spelling
    // coupling. A foreign decoder must first earn equality of its actual numeric
    // leaves; equality of parsed binary64 values is not that proof.
    for p in poles {
        w.charge(1)?;
        if p[2] != height {
            return Err(Error::Unresolved);
        }
    }
    Ok(())
}
fn positive(p: P, stack: &mut Vec<(P, usize)>, w: &mut Work<'_>) -> Result<f64, Error> {
    stack.push((p, 0));
    let mut minimum = f64::INFINITY;
    while let Some((p, depth)) = stack.pop() {
        w.charge(p.degree + 1)?;
        let h = p.hull();
        if !h.valid() {
            return Err(Error::Unresolved);
        }
        if h.lo > 0. {
            minimum = minimum.min(h.lo);
            continue;
        }
        if h.hi <= 0. || depth == DEPTH {
            return Err(Error::Unresolved);
        }
        let [a, b] = p.split(w)?;
        stack.push((b, depth + 1));
        stack.push((a, depth + 1));
    }
    Ok(minimum)
}
