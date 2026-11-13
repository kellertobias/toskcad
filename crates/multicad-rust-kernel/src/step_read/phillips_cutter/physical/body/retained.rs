//! Original numeric leaves remain the symbolic authority. Adjacent binary64
//! enclosures are evaluation caches, not replacement exact source coefficients.
use super::*;
mod memory;
pub(super) use memory::budget;
#[cfg(test)]
mod tests;
#[derive(Clone, Copy, Debug)]
pub struct OriginalNumericLeaf<'input> {
    spelling: &'input str,
    bounds: [f64; 2],
}
impl<'input> OriginalNumericLeaf<'input> {
    pub fn spelling(&self) -> &'input str {
        self.spelling
    }
    pub fn bounds(&self) -> [f64; 2] {
        self.bounds
    }
}
#[derive(Debug)]
pub struct OriginalNumericLeaves<'a, 'input> {
    body: &'a Body<'input>,
    entity: u32,
    leaves: Vec<OriginalNumericLeaf<'input>>,
    work_used: usize,
}
impl<'a, 'input> OriginalNumericLeaves<'a, 'input> {
    pub fn body(&self) -> &'a Body<'input> {
        self.body
    }
    pub fn entity(&self) -> u32 {
        self.entity
    }
    pub fn leaves(&self) -> &[OriginalNumericLeaf<'input>] {
        &self.leaves
    }
    pub fn work_used(&self) -> usize {
        self.work_used
    }
    pub fn retained_capacity_bytes(&self) -> usize {
        std::mem::size_of::<Self>().saturating_add(
            self.leaves
                .capacity()
                .saturating_mul(std::mem::size_of::<OriginalNumericLeaf<'input>>()),
        )
    }
}
#[derive(Clone, Copy, Debug)]
pub struct ExactPointLeaf<'a, 'input> {
    body: &'a Body<'input>,
    entity: u32,
    coordinates: [OriginalNumericLeaf<'input>; 3],
    work_used: usize,
}
impl<'a, 'input> ExactPointLeaf<'a, 'input> {
    pub fn body(&self) -> &'a Body<'input> {
        self.body
    }
    pub fn entity(&self) -> u32 {
        self.entity
    }
    pub fn coordinate_spellings(&self) -> [&'input str; 3] {
        self.coordinates.map(|c| c.spelling)
    }
    pub fn coordinate_bounds_mm(&self) -> [[f64; 2]; 3] {
        self.coordinates.map(|c| c.bounds)
    }
    pub fn work_used(&self) -> usize {
        self.work_used
    }
}
fn numbers_in(values: &[Value<'_>], work: &mut Work<'_>, numbers: &mut Vec<f64>) -> Result<()> {
    for value in values {
        work.step(1)?;
        match value {
            Value::Number(number) => numbers.push(*number),
            Value::List(values) => numbers_in(values, work, numbers)?,
            Value::Typed(entity) => numbers_in(&entity.args, work, numbers)?,
            _ => {}
        }
    }
    Ok(())
}
impl<'input> Body<'input> {
    /// Conservative source-retention charge, including the borrowed input span
    /// exactly once. Cached under reader Work before Body publication; this is
    /// not an exact allocator-footprint measurement or manufacturing receipt.
    pub fn retained_budget_bytes(&self) -> usize {
        self.retained_budget_bytes
    }
    /// Direct numeric leaves of one simple literal record, in source order.
    /// References stay references: no inferred coordinate/gauge is substituted.
    pub fn original_numeric_leaves<'a>(
        &'a self,
        entity: u32,
        work_limit: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<OriginalNumericLeaves<'a, 'input>, PhillipsReadError> {
        let mut work = Work {
            remaining: work_limit,
            interrupt,
            cancelled: Cell::new(false),
        };
        let result: Result<OriginalNumericLeaves<'a, 'input>> = (|| {
            work.step(1)?;
            let shell = Shell::new(&self.entities, self.tolerance, true);
            let record = shell.simple(entity)?;
            let source = record.source.ok_or(ReadError::Unsupported)?;
            work.step(source.len())?;
            let tokens = crate::step_read::part21::numeric_spellings(record)
                .ok_or(ReadError::Unsupported)?;
            let mut numbers = Vec::new();
            numbers_in(&record.args, &mut work, &mut numbers)?;
            require(tokens.len() == numbers.len())?;
            let mut leaves = Vec::with_capacity(tokens.len());
            for (spelling, value) in tokens.into_iter().zip(numbers) {
                work.step(1)?;
                let bounds = [value.next_down(), value.next_up()];
                require(bounds.into_iter().all(f64::is_finite))?;
                leaves.push(OriginalNumericLeaf { spelling, bounds });
            }
            Ok(OriginalNumericLeaves {
                body: self,
                entity,
                leaves,
                work_used: work_limit - work.remaining,
            })
        })();
        if work.cancelled.get() {
            Err(PhillipsReadError::Cancelled)
        } else {
            result.map_err(PhillipsReadError::from)
        }
    }
    pub fn exact_point_leaf<'a>(
        &'a self,
        entity: u32,
        work_limit: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<ExactPointLeaf<'a, 'input>, PhillipsReadError> {
        let leaves = self.original_numeric_leaves(entity, work_limit, interrupt)?;
        let shell = Shell::new(&self.entities, self.tolerance, true);
        let [Value::String(_), values] = shell
            .args(entity, "CARTESIAN_POINT")
            .map_err(PhillipsReadError::from)?
        else {
            return Err(PhillipsReadError::Invalid);
        };
        // Preserve the exact three-dimensional POINT schema and original tokens.
        numbers(values).map_err(PhillipsReadError::from)?;
        // The immutable record above is exactly String + List[Number;3]. The
        // generic collector therefore returned exactly three matching tokens;
        // there is no additional fallible conversion or fabricated point here.
        let coordinates = [leaves.leaves[0], leaves.leaves[1], leaves.leaves[2]];
        let mut finish = Work {
            remaining: work_limit - leaves.work_used,
            interrupt,
            cancelled: Cell::new(false),
        };
        let ready = finish.step(1);
        if finish.cancelled.get() {
            return Err(PhillipsReadError::Cancelled);
        }
        ready.map_err(PhillipsReadError::from)?;
        Ok(ExactPointLeaf {
            body: self,
            entity,
            coordinates,
            work_used: work_limit - finish.remaining,
        })
    }
}
