//! Blind finite rounded stock, preserving the actual cutter-tip trim graph.
use super::*;
use crate::{rounded_thread::BlindBore, tapered_thread::interval::I};

pub fn write_rounded_thread_blind_bore(
    bore: &BlindBore,
    origin: Point,
    u: Point,
    v: Point,
    tolerance: f64,
) -> Result<String, StepError> {
    write_one(|writer| writer.append_rounded_thread_blind_bore(bore, origin, u, v, tolerance))
}

pub(crate) fn dimensions(bore: &BlindBore, tolerance: f64) -> Result<(), StepError> {
    if !tolerance.is_finite() || tolerance <= 0. {
        return Err(StepError::InvalidGeometry);
    }
    let radius = bore.stock_radius_mm().as_f64();
    let height = bore.stock_height_mm().as_f64();
    let [rlo, rhi] = bore.stock_radius_bounds_mm();
    let [hlo, hhi] = bore.stock_height_bounds_mm();
    // Actual exact stock -> literal STEP scalar transfer, kept separate from
    // the original thread's already qualified finite curve/tensor fitting.
    let transfer = I::point(2.) * I::point((I::new(rlo, rhi) - I::point(radius)).magnitude())
        + I::point((I::new(hlo, hhi) - I::point(height)).magnitude());
    // Placement and independent reader both require the literal two-axis
    // Gram entries within1e-10 of an orthonormal frame. Including their cross
    // column gives operator norm<2; reserve this actual world amplification.
    let world_transfer = I::point(2.) * transfer;
    if !world_transfer.valid()
        || world_transfer.hi > tolerance.min(1e-7) / 32.
        || radius - bore.cutter_radial_bounds_mm()[1] <= 2. * tolerance.min(1e-7)
        || bore.floor_thickness_bounds_mm()[0] <= 2. * tolerance.min(1e-7)
    {
        return Err(StepError::Limit);
    }
    Ok(())
}

pub(crate) fn precision(bore: &BlindBore, origin: Point, tolerance: f64) -> Result<(), StepError> {
    dimensions(bore, tolerance)?;
    if origin.iter().any(|v| !v.is_finite()) {
        return Err(StepError::InvalidGeometry);
    }
    // The neutral reader caches REAL leaves as binary64. Cache equality is
    // not exact decimal equality. Reserve whole-stock/world conversion and
    // analytic-axis arithmetic, separately from the original thread fitting
    // floor. A large exact dyadic stock radius must not bypass this floor
    // merely because its authored transfer is zero. Unsupported sub-ULP
    // requests conservatively refuse, even for particular exact input leaves.
    let origin = origin.into_iter().map(f64::abs).fold(0., f64::max);
    let scale = I::point(bore.stock_radius_bounds_mm()[1])
        + I::point(bore.stock_height_bounds_mm()[1])
        + I::point(origin)
        + I::point(1.);
    let conversion = I::point(256. * f64::EPSILON) * scale;
    if !conversion.valid() || conversion.hi > tolerance.min(1e-7) / 32. {
        return Err(StepError::Limit);
    }
    Ok(())
}

pub(super) fn encode(
    file: &mut StepFile,
    bore: &BlindBore,
    placement: Placement,
    tolerance: f64,
) -> Result<Entity, StepError> {
    precision(bore, placement.origin, tolerance)?;
    let mut candidate = StepFile::new();
    let solid = rounded_thread::encode_blind_shell(&mut candidate, bore, &placement, tolerance)?;
    candidate.finish(vec![solid]);
    crate::step_read::verify_rounded_thread_blind_bore(
        &candidate.data,
        bore,
        placement.origin,
        placement.u,
        placement.v,
        tolerance,
    )
    .map_err(|e| {
        if e == crate::step_read::ReadError::Limit {
            StepError::Limit
        } else {
            StepError::InvalidGeometry
        }
    })?;
    rounded_thread::encode_blind_shell(file, bore, &placement, tolerance)
}
