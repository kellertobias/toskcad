//! Area and area-weighted moments in the immutable authored local chart.
//! These are source integrals; neither mesh points nor a world-area invariant
//! define them. A later affine placement needs its own normal-dependent proof.
use super::*;
mod hanger;
mod integrate;
mod metric;
mod series;
mod source;
#[cfg(test)]
mod tests;
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum ThreadBandRole {
    Lower = 0,
    Crest = 1,
    Upper = 2,
    Root = 3,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum FaceKey {
    StartCap,
    EndCap,
    Thread { start: u32, band: ThreadBandRole },
    EntryCone,
    TerminalShoulder,
    TerminalCylinder,
    MaterialOuterWall,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FacePropertyOptions {
    pub area_tolerance_mm2: f64,
    pub first_moment_tolerance_mm3: f64,
    pub centroid_tolerance_mm: f64,
    /// Maximum number of bounded interval cell proofs, including refinements.
    pub maximum_cells: usize,
}
impl Default for FacePropertyOptions {
    fn default() -> Self {
        Self {
            area_tolerance_mm2: 1e-7,
            first_moment_tolerance_mm3: 1e-7,
            centroid_tolerance_mm: 1e-7,
            maximum_cells: 4096,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FaceProperties {
    pub face: FaceKey,
    pub area_mm2: Enclosure,
    pub first_moment_mm3: [Enclosure; 3],
    pub centroid_mm: [Enclosure; 3],
}
impl Cavity {
    /// Semantic role count. Property queries separately qualify the supported
    /// boundary stage; this helper alone is not a topology certificate.
    pub fn face_count(&self) -> u32 {
        2 + 4 * self.thread.spec().starts
            + u32::from(self.entry.is_some())
            + 2 * u32::from(self.tail.is_some())
    }
    pub fn face_key(&self, index: u32) -> Option<FaceKey> {
        match index {
            0 => Some(FaceKey::StartCap),
            1 => Some(FaceKey::EndCap),
            _ => {
                let n = self.thread.spec().starts;
                if index < 2 + 4 * n {
                    let band = match (index - 2) % 4 {
                        0 => ThreadBandRole::Lower,
                        1 => ThreadBandRole::Crest,
                        2 => ThreadBandRole::Upper,
                        _ => ThreadBandRole::Root,
                    };
                    Some(FaceKey::Thread {
                        start: (index - 2) / 4,
                        band,
                    })
                } else {
                    let index = index - (2 + 4 * n);
                    if self.entry.is_some() && index == 0 {
                        Some(FaceKey::EntryCone)
                    } else if self.tail.is_some() {
                        match index - u32::from(self.entry.is_some()) {
                            0 => Some(FaceKey::TerminalShoulder),
                            1 => Some(FaceKey::TerminalCylinder),
                            _ => None,
                        }
                    } else {
                        None
                    }
                }
            }
        }
    }
    pub fn face_index(&self, key: FaceKey) -> Option<u32> {
        (0..self.face_count()).find(|&i| self.face_key(i) == Some(key))
    }
    /// Period is the actual source q/pitch period, including negative q/v cells.
    pub fn thread_face_key(&self, period: i64, band: usize) -> Result<FaceKey, Error> {
        thread_face_key(&self.thread, period, band)
    }

    pub fn local_face_properties(
        &self,
        face: FaceKey,
        options: FacePropertyOptions,
    ) -> Result<FaceProperties, Error> {
        self.local_face_properties_with_interrupt(face, options, None)
    }
    pub fn local_face_properties_with_interrupt(
        &self,
        face: FaceKey,
        options: FacePropertyOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<FaceProperties, Error> {
        query(self, None, face, options, interrupt)
    }
}
impl Material {
    pub fn face_count(&self) -> u32 {
        self.cavity().face_count() + 1
    }
    pub fn face_key(&self, index: u32) -> Option<FaceKey> {
        if index == self.cavity().face_count() {
            Some(FaceKey::MaterialOuterWall)
        } else {
            self.cavity().face_key(index)
        }
    }
    pub fn face_index(&self, key: FaceKey) -> Option<u32> {
        if key == FaceKey::MaterialOuterWall {
            Some(self.cavity().face_count())
        } else {
            self.cavity().face_index(key)
        }
    }
    /// The exact stock center/bottom is removed, as in stock_chart_origin().
    pub fn local_face_properties(
        &self,
        face: FaceKey,
        options: FacePropertyOptions,
    ) -> Result<FaceProperties, Error> {
        self.local_face_properties_with_interrupt(face, options, None)
    }
    pub fn local_face_properties_with_interrupt(
        &self,
        face: FaceKey,
        options: FacePropertyOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<FaceProperties, Error> {
        query(self.cavity(), Some(self), face, options, interrupt)
    }
}
fn query(
    body: &Cavity,
    material: Option<&Material>,
    face: FaceKey,
    options: FacePropertyOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<FaceProperties, Error> {
    let mut work = integrate::Work::new(options, interrupt)?;
    if material
        .map_or_else(|| body.face_index(face), |m| m.face_index(face))
        .is_none()
    {
        return Err(Error::Invalid);
    }
    // This first graph/property stage retains the same explicit positive-entry
    // domain. Small valid entries need source threshold splitting, not g<0.
    if body
        .entry
        .is_some_and(|e| e.radius <= body.thread.maximum_radius_mm())
    {
        return Err(Error::Unresolved);
    }
    let values = source::evaluate(
        source::Domain::treated(body),
        material.map(Material::radius_interval),
        face,
        options,
        &mut work,
    )?;
    let result = finish(face, values, options)?;
    work.checkpoint()?;
    Ok(result)
}
fn enclosure(i: I) -> Enclosure {
    Enclosure {
        lower: i.lo,
        upper: i.hi,
    }
}
struct Measures {
    area_mm2: Enclosure,
    first_moment_mm3: [Enclosure; 3],
    centroid_mm: [Enclosure; 3],
}
fn finish(
    face: FaceKey,
    values: [I; 4],
    options: FacePropertyOptions,
) -> Result<FaceProperties, Error> {
    let values = finish_values(values, options)?;
    Ok(FaceProperties {
        face,
        area_mm2: values.area_mm2,
        first_moment_mm3: values.first_moment_mm3,
        centroid_mm: values.centroid_mm,
    })
}
fn finish_values(values: [I; 4], options: FacePropertyOptions) -> Result<Measures, Error> {
    let values = checked_values(values)?;
    if values[0].lo <= 0. {
        return Err(Error::Unresolved);
    }
    let centroid = std::array::from_fn::<_, 3, _>(|i| values[i + 1] / values[0]);
    if values[0].width() > options.area_tolerance_mm2
        || values[1..]
            .iter()
            .any(|v| v.width() > options.first_moment_tolerance_mm3)
        || centroid
            .iter()
            .any(|v| !v.valid() || v.width() > options.centroid_tolerance_mm)
    {
        return Err(Error::Unresolved);
    }
    Ok(Measures {
        area_mm2: enclosure(values[0]),
        first_moment_mm3: std::array::from_fn(|i| enclosure(values[i + 1])),
        centroid_mm: centroid.map(enclosure),
    })
}

fn checked_values(values: [I; 4]) -> Result<[I; 4], Error> {
    if values.into_iter().all(I::valid) {
        Ok(values)
    } else {
        Err(Error::Unresolved)
    }
}

fn thread_face_key(thread: &AuthoredThread, period: i64, band: usize) -> Result<FaceKey, Error> {
    if band >= 5 {
        return Err(Error::Invalid);
    }
    let n = i64::from(thread.spec().starts);
    let group = (period.rem_euclid(n) + [0, 0, 0, 1, 1][band]) % n;
    let start = if thread.spec().left_hand {
        group
    } else {
        (n - group) % n
    } as u32;
    let role = [
        ThreadBandRole::Crest,
        ThreadBandRole::Upper,
        ThreadBandRole::Root,
        ThreadBandRole::Lower,
        ThreadBandRole::Crest,
    ][band];
    Ok(FaceKey::Thread { start, band: role })
}

impl AuthoredThread {
    /// Measures the actual shared-affine finite source, without a treatment.
    pub fn local_face_properties(
        &self,
        face: FaceKey,
        options: FacePropertyOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<FaceProperties, Error> {
        bare_properties(self, None, face, options, interrupt)
    }
}
/// Optional subtraction is only the paired source's actual planar cap domain.
/// Both integrals share one work ledger and retain their own phase/start law.
pub(crate) fn bare_properties(
    source: &AuthoredThread,
    subtract: Option<&AuthoredThread>,
    face: FaceKey,
    options: FacePropertyOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<FaceProperties, Error> {
    let mut work = integrate::Work::new(options, interrupt)?;
    let valid = match face {
        FaceKey::StartCap | FaceKey::EndCap => true,
        FaceKey::Thread { start, .. } => start < source.spec().starts && subtract.is_none(),
        _ => false,
    };
    if !valid || subtract.is_some_and(|s| s.full_length_mm() != source.full_length_mm()) {
        return Err(Error::Invalid);
    }
    // Independent operands each receive half the requested integral widths.
    // The final quotient/centroid width is checked after directed subtraction.
    let integration_options = if subtract.is_some() {
        FacePropertyOptions {
            area_tolerance_mm2: options.area_tolerance_mm2 * 0.5,
            first_moment_tolerance_mm3: options.first_moment_tolerance_mm3 * 0.5,
            centroid_tolerance_mm: options.centroid_tolerance_mm * 0.5,
            ..options
        }
    } else {
        options
    };
    let mut values = source::evaluate(
        source::Domain::bare(source),
        None,
        face,
        integration_options,
        &mut work,
    )?;
    if let Some(inner) = subtract {
        let other = source::evaluate(
            source::Domain::bare(inner),
            None,
            face,
            integration_options,
            &mut work,
        )?;
        for axis in 0..4 {
            values[axis] = values[axis] - other[axis];
        }
    }
    let result = finish(face, values, options)?;
    work.checkpoint()?;
    Ok(result)
}
