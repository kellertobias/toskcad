//! The literal flux becomes occupied volume only after the same packet earns
//! global embedding, cap Jordan degree, material nesting and outward orientation.
use super::*;
pub struct OccupiedVolume<'a, 'b, 'c, 'input> {
    material: &'a MaterialBounds<'b, 'c, 'input>,
    flux: OrientedFlux<'b, 'c, 'input>,
    work: usize,
    bytes: usize,
}
impl<'a, 'b, 'c, 'input> OccupiedVolume<'a, 'b, 'c, 'input> {
    pub fn material(&self) -> &'a MaterialBounds<'b, 'c, 'input> {
        self.material
    }
    pub fn volume_mm3(&self) -> [f64; 2] {
        self.flux.oriented_mm3()
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    /// Conservative retained/query storage, not process resident memory.
    pub fn query_budget_bytes(&self) -> usize {
        self.bytes
    }
}
impl<'b, 'c, 'input> MaterialBounds<'b, 'c, 'input> {
    pub fn occupied_volume<'a>(
        &'a self,
        options: FluxOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<OccupiedVolume<'a, 'b, 'c, 'input>> {
        options.validate()?;
        let bindings = self.binding();
        // Both receipts include the identical parse/C0 prerequisite once.
        // Their remaining charges accumulate, without renewing either gate.
        let extra_work = self.work_used() - bindings.work_used() + 1;
        let extra_bytes = self.retained_budget_bytes() - bindings.retained_budget_bytes()
            + size_of::<OccupiedVolume<'_, '_, '_, '_>>();
        if extra_work >= options.maximum_work || extra_bytes >= options.maximum_bytes {
            return Err(Error::Limit);
        }
        let flux = bindings.oriented_flux(
            FluxOptions {
                maximum_work: options.maximum_work - extra_work,
                maximum_bytes: options.maximum_bytes - extra_bytes,
                ..options
            },
            interrupt,
        )?;
        // The global material proof already establishes outwardness. Refuse
        // any numerical result that cannot certify strictly positive volume.
        if flux.oriented_mm3()[0] <= 0. {
            return Err(Error::Unsupported);
        }
        if interrupt.is_some_and(|cancel| cancel()) {
            return Err(Error::Cancelled);
        }
        Ok(OccupiedVolume {
            material: self,
            work: extra_work + flux.work_used(),
            bytes: extra_bytes + flux.query_budget_bytes(),
            flux,
        })
    }
}
