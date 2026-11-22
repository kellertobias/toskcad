//! Source-owned finite Cartesian side support, including its varying Z.
//! Finite edge trims, shared junctions and closed-file admission are separate.
use super::super::literal::{CartesianCap, Precision, encode_spatial_at};
use super::*;
pub struct LiteralSideSupport<'q, 'f, 's> {
    source: &'q SideSupport<'f, 's>,
    geometry: CartesianCap,
    bounds: [f64; 2],
    work: usize,
    bytes: usize,
}
impl LiteralSideSupport<'_, '_, '_> {
    pub fn source(&self) -> &SideSupport<'_, '_> {
        self.source
    }
    pub fn points(&self) -> &[[EncodedScalar; 3]] {
        &self.geometry.poles
    }
    pub fn weights(&self) -> &[EncodedScalar] {
        &self.geometry.weights
    }
    pub fn source_bounds(&self) -> [f64; 2] {
        self.bounds
    }
    pub fn field_bounds(&self) -> [f64; 2] {
        self.geometry.error
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
}
impl<'f, 's> SideSupport<'f, 's> {
    pub fn literal<'q>(
        &'q self,
        options: CurveOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<LiteralSideSupport<'q, 'f, 's>, Error> {
        options.validate()?;
        let (mut w, mut m) = crate::step::hanger::domain::ledgers(
            DomainOptions {
                maximum_work: options.work,
                maximum_bytes: options.bytes,
            },
            self.field.retained_bytes() + self.retained_bytes(),
            interrupt,
        )?;
        w.charge(self.field.work_used() + self.work_used())?;
        m.add(SCRATCH)?;
        let base = m.used;
        m.add(size_of::<LiteralSideSupport<'_, '_, '_>>())?;
        let mut h = m.allocate(self.poles.len(), &mut w)?;
        for pole in &self.poles {
            let weight = pole.weight.number(&mut w)?;
            h.push([
                pole.point[0].number(&mut w)?.mul(&weight, &mut w)?,
                pole.point[1].number(&mut w)?.mul(&weight, &mut w)?,
                pole.point[2].number(&mut w)?.mul(&weight, &mut w)?,
                weight,
            ]);
        }
        let geometry = encode_spatial_at(&h, options, &mut m, &mut w, Precision::Bits64)?;
        m.release(&h);
        drop(h);
        let parameter = self.parameter();
        let span = binary(parameter[1], &mut w)?
            .sub(&binary(parameter[0], &mut w)?, &mut w)?
            .enclosure(&mut w)?;
        let position = p(geometry.error[0]) + p(self.position_error_mm());
        let derivative = p(geometry.error[1]) / span + p(self.derivative_error_mm());
        if !position.valid()
            || !derivative.valid()
            || position.hi > options.position_error_mm
            || derivative.hi > options.derivative_error_mm
        {
            return Err(Error::Unresolved);
        }
        w.charge(1)?;
        Ok(LiteralSideSupport {
            source: self,
            geometry,
            bounds: [position.hi, derivative.hi],
            work: w.used,
            bytes: m.used - base,
        })
    }
}
#[cfg(test)]
mod tests;
