//! Analytic support for independent trim correspondence. This does not
//! grant generic B-rep admission to a self-intersecting complete torus.
use super::*;

pub(super) fn read(
    source: &Source<'_>,
    definition: &uv::Definition,
    budget: &mut Budget<'_>,
) -> Result<Surface> {
    match source.surface(definition.surface) {
        Ok(surface) => return Ok(surface),
        Err(ReadError::Unsupported) => {}
        Err(other) => return Err(other),
    }
    budget.charge(4096)?;
    let Some([_, placement, major, minor]) =
        source.component(definition.surface, "TOROIDAL_SURFACE")?
    else {
        return Err(ReadError::Unsupported);
    };
    check(super::si_radians(source, budget)?, ReadError::Unsupported)?;
    let raw_major = precision::scalar(source, definition.surface, major, budget)?;
    let raw_minor = precision::scalar(source, definition.surface, minor, budget)?;
    check(
        raw_major.valid() && raw_major.lo > 0. && raw_minor.valid() && raw_minor.lo > 0.,
        ReadError::Unsupported,
    )?;
    // Only the analytic support is made available to the
    // independent trim certificates. Source::surface keeps its original
    // major>=minor requirement for complete retained-solid construction.
    Ok(Surface::Torus {
        frame: source.frame(reference(placement)?)?,
        major: radius(major, false)?,
        minor: radius(minor, false)?,
    })
}

#[cfg(test)]
mod tests;
