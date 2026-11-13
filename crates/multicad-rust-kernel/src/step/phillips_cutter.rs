//! Authored Phillips cutter export. Readback authority remains source-known;
//! this does not admit a screw recess or reconstruct a foreign cutter recipe.
use super::*;
use crate::{
    phillips_drive::{AnalyticCurve, BoundaryTopology, CubicTrimSpan, CurveKey, Cutter, Surface},
    tapered_thread::interval::I,
};
pub(crate) mod curves;
mod plan;
pub(crate) mod support;
#[cfg(test)]
mod tests;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    InvalidGeometry,
    InvalidPlacement,
    Unsupported,
    Limit,
    Cancelled,
}
impl From<StepError> for Error {
    fn from(e: StepError) -> Self {
        match e {
            StepError::InvalidGeometry => Self::InvalidGeometry,
            StepError::InvalidPlacement => Self::InvalidPlacement,
            StepError::Unsupported => Self::Unsupported,
            StepError::Limit => Self::Limit,
        }
    }
}
impl From<crate::phillips_drive::Error> for Error {
    fn from(e: crate::phillips_drive::Error) -> Self {
        match e {
            crate::phillips_drive::Error::Invalid => Self::InvalidGeometry,
            crate::phillips_drive::Error::Unresolved => Self::Unsupported,
            crate::phillips_drive::Error::Limit => Self::Limit,
            crate::phillips_drive::Error::Cancelled => Self::Cancelled,
        }
    }
}
fn checkpoint(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), Error> {
    if interrupt.is_some_and(|f| f()) {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
pub fn write_phillips_cutter(
    source: &Cutter,
    origin: Point,
    u: Point,
    v: Point,
    tolerance_mm: f64,
) -> Result<String, Error> {
    write_phillips_cutter_with_interrupt(source, origin, u, v, tolerance_mm, None)
}
pub fn write_phillips_cutter_with_interrupt(
    source: &Cutter,
    origin: Point,
    u: Point,
    v: Point,
    tolerance_mm: f64,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<String, Error> {
    let mut writer = Writer::new();
    writer.append_phillips_cutter_with_interrupt(source, origin, u, v, tolerance_mm, interrupt)?;
    writer.finish().map_err(Error::from)
}
impl Writer {
    pub fn append_phillips_cutter(
        &mut self,
        source: &Cutter,
        origin: Point,
        u: Point,
        v: Point,
        tolerance_mm: f64,
    ) -> Result<(), Error> {
        self.append_phillips_cutter_with_interrupt(source, origin, u, v, tolerance_mm, None)
    }
    pub fn append_phillips_cutter_with_interrupt(
        &mut self,
        source: &Cutter,
        origin: Point,
        u: Point,
        v: Point,
        tolerance_mm: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<(), Error> {
        let pose = support::Pose::new(origin, u, v)?;
        self.append_typed(|file| {
            let plan = plan::Plan::new(source, &pose, tolerance_mm, interrupt)?;
            plan.emit(file, interrupt)
        })
    }
}
