//! Actual authored Phillips screw exchange; independent import remains gated.
use super::*;
use crate::{phillips_drive::PhillipsScrew, step::phillips_cutter::Error};
pub(super) fn source_error(e: crate::phillips_drive::Error) -> StepError {
    match e {
        crate::phillips_drive::Error::Invalid => StepError::InvalidGeometry,
        crate::phillips_drive::Error::Unresolved
        | crate::phillips_drive::Error::Limit
        | crate::phillips_drive::Error::Cancelled => StepError::Limit,
    }
}
pub(super) fn plan_error(e: Error) -> StepError {
    match e {
        Error::InvalidGeometry => StepError::InvalidGeometry,
        Error::InvalidPlacement => StepError::InvalidPlacement,
        Error::Unsupported => StepError::Unsupported,
        Error::Limit | Error::Cancelled => StepError::Limit,
    }
}
pub fn write_phillips_screw(
    source: &PhillipsScrew,
    origin: Point,
    u: Point,
    v: Point,
    tolerance_mm: f64,
) -> Result<String, Error> {
    write_phillips_screw_with_interrupt(source, origin, u, v, tolerance_mm, None)
}
pub fn write_phillips_screw_with_interrupt(
    source: &PhillipsScrew,
    origin: Point,
    u: Point,
    v: Point,
    tolerance_mm: f64,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<String, Error> {
    let mut writer = Writer::new();
    writer.append_phillips_screw_with_interrupt(source, origin, u, v, tolerance_mm, interrupt)?;
    writer.finish().map_err(Error::from)
}
impl Writer {
    pub fn append_phillips_screw(
        &mut self,
        source: &PhillipsScrew,
        origin: Point,
        u: Point,
        v: Point,
        tolerance_mm: f64,
    ) -> Result<(), Error> {
        self.append_phillips_screw_with_interrupt(source, origin, u, v, tolerance_mm, None)
    }
    pub fn append_phillips_screw_with_interrupt(
        &mut self,
        screw: &PhillipsScrew,
        origin: Point,
        u: Point,
        v: Point,
        tolerance_mm: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<(), Error> {
        let p = Placement::new(origin, u, v)?;
        let aborted = std::cell::Cell::new(false);
        let stop = || {
            aborted.set(aborted.get() || interrupt.is_some_and(|f| f()));
            aborted.get()
        };
        self.append_typed(|file| {
            let source = HeadSource::phillips(screw);
            crate::step::helical::encode_with_interrupt(
                file,
                screw.thread(),
                Closure::Profiled(source),
                p,
                tolerance_mm,
                Some(&stop),
            )
            .map_err(|error| {
                if aborted.get() {
                    Error::Cancelled
                } else {
                    Error::from(error)
                }
            })
        })
    }
}

#[cfg(test)]
mod tests;
