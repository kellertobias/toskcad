//! The STEP-defined transform, not a fitted Native matrix or a source recipe.
use super::*;
pub struct Placement<'selected, 'document, 'input> {
    selected: &'selected Selected<'document, 'input>,
    linear: [[I; 3]; 3],
    inverse: [[I; 3]; 3],
    origin: [I; 3],
    inverse_origin: [I; 3],
    determinant: I,
    work: usize,
    bytes: usize,
}
fn interval([lo, hi]: [f64; 2]) -> I {
    I::new(lo, hi)
}
fn bounds(i: I) -> [f64; 2] {
    [i.lo, i.hi]
}
fn mapping_error(error: crate::step_read::owned_insert::Error) -> Error {
    match error {
        crate::step_read::owned_insert::Error::Invalid => Error::Invalid,
        crate::step_read::owned_insert::Error::Unsupported => Error::Unsupported,
        crate::step_read::owned_insert::Error::Limit => Error::Limit,
        crate::step_read::owned_insert::Error::Cancelled => Error::Cancelled,
    }
}
pub(super) fn transport(
    columns: [[I; 3]; 3],
    mut origin: [I; 3],
    input: [[f64; 2]; 3],
) -> Result<[[f64; 2]; 3]> {
    let input = input.map(interval);
    if input.iter().any(|i| !i.valid()) {
        return Err(Error::Invalid);
    }
    for (axis, value) in origin.iter_mut().enumerate() {
        for (column, input) in columns.iter().zip(input) {
            *value = *value + column[axis] * input;
        }
        if !value.valid() {
            return Err(Error::Unsupported);
        }
    }
    Ok(origin.map(bounds))
}
impl<'selected, 'document, 'input> Placement<'selected, 'document, 'input> {
    pub fn selected(&self) -> &'selected Selected<'document, 'input> {
        self.selected
    }
    pub fn mapping_link(&self) -> Option<MappingLink> {
        self.selected.member.mapping_link()
    }
    /// Bounds of the exact EXPRESS normalized/projected/scaled operator.
    /// These are correlated enclosures, not a new defining binary64 matrix.
    pub fn linear_columns_bounds(&self) -> [[[f64; 2]; 3]; 3] {
        self.linear.map(|v| v.map(bounds))
    }
    pub fn inverse_columns_bounds(&self) -> [[[f64; 2]; 3]; 3] {
        self.inverse.map(|v| v.map(bounds))
    }
    pub fn origin_mm(&self) -> [[f64; 2]; 3] {
        self.origin.map(bounds)
    }
    pub fn inverse_origin_mm(&self) -> [[f64; 2]; 3] {
        self.inverse_origin.map(bounds)
    }
    pub fn determinant_bounds(&self) -> [f64; 2] {
        bounds(self.determinant)
    }
    pub fn map_box_mm(&self, input: [[f64; 2]; 3]) -> Result<[[f64; 2]; 3]> {
        transport(self.linear, self.origin, input)
    }
    pub fn inverse_box_mm(&self, input: [[f64; 2]; 3]) -> Result<[[f64; 2]; 3]> {
        transport(self.inverse, self.inverse_origin, input)
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn query_budget_bytes(&self) -> usize {
        self.bytes
    }
}
impl<'document, 'input> Selected<'document, 'input> {
    pub fn placement<'selected>(
        &'selected self,
        options: Options,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Placement<'selected, 'document, 'input>> {
        let mut budget = Budget::new(options, interrupt)?;
        budget.step(self.parsed.work)?;
        let fixed = size_of::<Placement<'_, '_, '_>>();
        if self.parsed.bytes + fixed > budget.limit {
            return Err(Error::Limit);
        }
        budget.bytes = self.parsed.bytes + fixed;
        budget.peak = budget.bytes;
        let mut axes =
            std::array::from_fn(|j| std::array::from_fn(|i| I::point(f64::from(i == j))));
        let mut origin = [I::point(0.); 3];
        let mut scale = I::point(1.);
        if let Some(link) = self.member.mapping_link() {
            // Remaining work belongs to this SAME selected-member ledger.
            if budget.remaining == 0 || budget.limit == budget.bytes {
                return Err(Error::Limit);
            }
            let operator = mapping::evaluate_neutral_operator(
                self.parsed.document.records(),
                link.operator_entity(),
                MappingOptions {
                    maximum_work: budget.remaining,
                    retained_budget_bytes: budget.limit - budget.bytes,
                },
                interrupt,
            )
            .map_err(mapping_error)?;
            budget.step(operator.work)?;
            budget.peak = budget.peak.max(budget.bytes + operator.bytes);
            if operator.handedness != 1 {
                return Err(Error::Unsupported);
            }
            origin = operator.origin.map(interval);
            axes = operator.axes.map(|v| v.map(interval));
            scale = interval(operator.scale);
        }
        budget.step(256)?;
        let linear = axes.map(|v| v.map(|x| x * scale));
        // EXPRESS construction earns an EXACT orthonormal basis. This
        // transpose/scale inverse follows that schema identity; it is NEVER
        // the inverse of a rounded approximately-proper Native matrix.
        let inverse = std::array::from_fn(|j| std::array::from_fn(|i| axes[i][j] / scale));
        let inverse_origin = std::array::from_fn(|i| {
            -(0..3).fold(I::point(0.), |sum, j| sum + inverse[j][i] * origin[j])
        });
        let determinant = scale.square() * scale;
        if linear
            .iter()
            .chain(&inverse)
            .flatten()
            .chain(&inverse_origin)
            .any(|i| !i.valid())
            || !determinant.valid()
            || determinant.lo <= 0.
        {
            return Err(Error::Unsupported);
        }
        Ok(Placement {
            selected: self,
            linear,
            inverse,
            origin,
            inverse_origin,
            determinant,
            work: budget.used,
            bytes: budget.peak,
        })
    }
}
