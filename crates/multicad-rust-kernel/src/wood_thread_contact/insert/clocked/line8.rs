//! Source-aware Line8 entry point sharing the established clocked projection,
//! sealed source boxes, measured alignment and whole-period radius proof.
use super::*;
#[derive(Clone, Copy)]
pub struct ClockedLine8Placement<'a> {
    pub source: &'a crate::line8_nut::Nut,
    pub stock_frame: Frame,
    pub outer_frame: Frame,
    pub minimum_wall_mm: f64,
}
pub struct ClockedLine8Mating<'a> {
    nut: ClockedLine8Placement<'a>,
    external: MatingExternal<'a>,
    pub mating: ClockedMatingBounds,
    pub radial_map_error_mm: f64,
    pub axial_map_error_mm: f64,
    pub angle_error_rad: f64,
    pub cells_used: usize,
    pub relief_clearance_lower_mm: f64,
    pub stock_wall_lower_mm: f64,
    pub complete_turns_lower: f64,
    pub accounted_peak_bytes: usize,
}
impl<'a> ClockedLine8Mating<'a> {
    pub fn nut(&self) -> ClockedLine8Placement<'a> {
        self.nut
    }
    pub fn external(&self) -> MatingExternal<'a> {
        self.external
    }
}
pub fn certify_clocked_line8_mating<'a>(
    nut: ClockedLine8Placement<'a>,
    external: MatingExternal<'a>,
    options: ClockedMatingOptions,
    interrupt: &dyn Fn() -> bool,
) -> R<ClockedLine8Mating<'a>> {
    checkpoint(interrupt)?;
    let MatingScrew::Profiled(screw) = external.source else {
        return Err(Error::Unsupported);
    };
    if screw.recess().flat_hex_geometry().is_none() {
        return Err(Error::Unsupported);
    }
    if !nut.minimum_wall_mm.is_finite()
        || nut.minimum_wall_mm <= 0.
        || nut.minimum_wall_mm > 1e7
        || options.maximum_cells == 0
        || options.maximum_cells > 4096
    {
        return Err(Error::Invalid);
    }
    // The established16384-bit scalar guard remains unchanged. Compositions
    // normalize each bounded scalar; raw products/sums have a finite <=65536-bit
    // temporary bound.8MiB covers128live exact values, pinned GCD scratch and
    // interval caches; the existing phase frontier adds at most4097cells.
    let bytes = nut.source.retained_bytes()
        + std::mem::size_of::<crate::profiled_screw::Screw>()
        + 8 * 1024 * 1024
        + (options.maximum_cells + 1) * std::mem::size_of::<[f64; 2]>();
    if bytes > 64 * 1024 * 1024 {
        return Err(Error::Limit);
    }
    let stock = SourceFrame::literal(nut.stock_frame)?
        .translated(0, nut.source.bore_from_start_mm())?
        .placed(nut.outer_frame)?;
    let tap = stock.copied().phase(nut.source.phase_rotation())?;
    let wall =
        finite(ip(nut.source.minimum_thread_wall_lower_bound_mm()) * stock.minimum_scale())?.lo;
    if wall < nut.minimum_wall_mm {
        return Err(Error::Unsupported);
    }
    let receiver = Receiver {
        source: nut.source.thread(),
        length: nut.source.relief().start_mm,
        frame: tap,
    };
    let shared = certify_clocked_receiver(&receiver, external, options, interrupt)?;
    let turns = finite(shared.external_span / ip(screw.thread().spec().pitch_mm))?;
    if turns.lo < 2. {
        return Err(Error::Unsupported);
    }
    let projection =
        Projection::expressions(&SourceFrame::literal(external.frame)?, &stock, interrupt)?;
    // Full actual head/neck/plain support versus the ORIGINAL stock entry
    // halfspace is a separate obligation from the tap radial clearance.
    let prefix = projection.enclose(SourceBox::prefix_source(external.source)?);
    if prefix[2].hi >= 0. {
        return Err(Error::Unsupported);
    }
    let prefix_gap =
        finite(-prefix[2] * I::new(stock.minimum_scale().lo, stock.maximum_scale().hi))?;
    let relief_map = align_source(
        external.source,
        SourceBox::thread_source(external.source)?,
        projection,
    )?;
    let radius = screw.thread().source_coordinate_bounds_mm()[0][1];
    let clearance = finite(
        (ip(nut.source.relief().radius_mm)
            - relief_map.scale * ip(radius)
            - ip(relief_map.radial_error))
            * stock.minimum_scale(),
    )?;
    if clearance.lo < options.requirements.minimum_clearance_mm {
        return Err(Error::Unsupported);
    }
    // Above the effective4.099 prefix the actual relief contains this whole
    // original thread envelope wherever original stock exists; below entry the
    // original stock is absent. No thread continuation is discarded.
    let mut mating = shared.bounds;
    let engagement = receiver.frame.normal_plane_distance(receiver.length)?;
    if engagement.lo < options.requirements.minimum_engagement_mm {
        return Err(Error::Unsupported);
    }
    mating.engagement_mm = [engagement.lo, engagement.hi];
    mating.prefix_entry_clearance_mm = [prefix_gap.lo, prefix_gap.hi];
    checkpoint(interrupt)?;
    Ok(ClockedLine8Mating {
        nut,
        external,
        mating,
        radial_map_error_mm: shared.radial_map_error_mm,
        axial_map_error_mm: shared.axial_map_error_mm,
        angle_error_rad: shared.angle_error_rad,
        cells_used: shared.cells_used,
        relief_clearance_lower_mm: clearance.lo,
        stock_wall_lower_mm: wall,
        complete_turns_lower: turns.lo,
        accounted_peak_bytes: bytes,
    })
}
#[cfg(test)]
mod tests;
