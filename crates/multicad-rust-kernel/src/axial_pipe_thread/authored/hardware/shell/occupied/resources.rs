//! Bounded-phase reservations for the unchanged retained-material constructors.
//! These are conservative upper bounds, never claimed as measured phase work.
use super::*;
use std::mem::size_of;
// Ownership2M, trims2M, axial2M plus nested radialcharts1M, caps2M,
// ordering2M, material64M. Each unchanged constructor owns its fixed ceiling.
pub(super) const DELEGATED_WORK: usize = 75_000_000;
pub(super) struct Reservation {
    pub source: usize,
    pub retained: usize,
    pub peak: usize,
}
fn add(total: &mut usize, n: usize) -> Result<(), Error> {
    *total = total.checked_add(n).ok_or(Error::Limit)?;
    Ok(())
}
fn bytes<T>(n: usize) -> Result<usize, Error> {
    n.checked_mul(size_of::<T>()).ok_or(Error::Limit)
}
pub(super) fn reserve(s: &Input, options: Options) -> Result<Reservation, Error> {
    if options.maximum_work == 0
        || options.maximum_work > 96_000_000
        || options.maximum_bytes == 0
        || options.maximum_bytes > 64 * 1024 * 1024
    {
        return Err(Error::Invalid);
    }
    if options.maximum_work < DELEGATED_WORK {
        return Err(Error::Limit);
    }
    let mut source = size_of::<Input>();
    add(&mut source, bytes::<[f64; 3]>(s.vertices.capacity())?)?;
    add(&mut source, bytes::<Curve>(s.curves.capacity())?)?;
    add(&mut source, bytes::<radial::Edge>(s.edges.capacity())?)?;
    add(&mut source, bytes::<Surface>(s.surfaces.capacity())?)?;
    add(&mut source, bytes::<Face>(s.faces.capacity())?)?;
    let mut tensor_controls = 0usize;
    let mut spline_controls = 0usize;
    for curve in &s.curves {
        if let Curve::Original(radial::Curve::Spline { poles, .. }) = curve {
            add(&mut source, bytes::<[f64; 3]>(poles.capacity())?)?;
            add(&mut spline_controls, poles.len())?;
        }
    }
    for surface in &s.surfaces {
        if let Surface::Tensor(t) = surface {
            // RetainedTensor holds Arc<[Point]>, with no spare slice capacity.
            // Duplicate Arc references are conservatively counted repeatedly.
            add(&mut source, bytes::<[f64; 3]>(t.poles().len())?)?;
            add(&mut source, 2 * size_of::<usize>())?;
            add(&mut tensor_controls, t.poles().len())?;
        }
    }
    for face in &s.faces {
        add(&mut source, bytes::<(usize, bool)>(face.outer.capacity())?)?;
        add(
            &mut source,
            bytes::<Vec<(usize, bool)>>(face.inner.capacity())?,
        )?;
        for ring in &face.inner {
            add(&mut source, bytes::<(usize, bool)>(ring.capacity())?)?;
        }
    }
    let mut base = source;
    add(&mut base, options.additional_resident_bytes)?;
    // Capped-thread spline clones are admitted before allocating/copying poles.
    add(&mut base, bytes::<[f64; 3]>(spline_controls)?)?;
    let mut retained = base;
    // Projected chart controls48B each, possible Vec growth/reversal slack,
    // both retained projected rail collections, their headers and table maps.
    add(
        &mut retained,
        tensor_controls.checked_mul(128).ok_or(Error::Limit)?,
    )?;
    add(
        &mut retained,
        spline_controls.checked_mul(256).ok_or(Error::Limit)?,
    )?;
    add(&mut retained, 2 * 1024 * 1024)?;
    let mut peak = base;
    // Doubled live chart/reversal arrays and rail construction temporaries.
    // Source tensor Arc clones share storage already counted above.
    add(
        &mut peak,
        tensor_controls.checked_mul(256).ok_or(Error::Limit)?,
    )?;
    add(
        &mut peak,
        spline_controls.checked_mul(512).ok_or(Error::Limit)?,
    )?;
    // All tables are capped at96vertices/128edges/48faces. Material DFS depth32
    // retains at most33 degree18^2 interval arrays; old+two split outputs,
    // 100point cell, weights and pending capacity need less than256KiB.
    // The4MiB phase allowance also encloses trims/caps/order DFS, metadata,
    // bounded fan maps and temporary cap/reindex table allocations.
    add(&mut peak, 4 * 1024 * 1024)?;
    if peak > options.maximum_bytes {
        return Err(Error::Limit);
    }
    Ok(Reservation {
        source,
        retained,
        peak,
    })
}
