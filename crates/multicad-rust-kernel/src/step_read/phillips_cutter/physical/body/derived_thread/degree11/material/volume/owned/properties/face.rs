//! Per-face actual material measures. Queries never measure display triangles.
use super::*;
use crate::{
    head_profile::RepairedHeadSurfaceOptions,
    phillips_drive::{FaceKey, FacePropertyOptions},
    tapered_thread::SurfaceOptions11,
};
#[derive(Clone, Copy, Debug)]
pub struct Degree11FacePropertyOptions {
    pub area_width_mm2: f64,
    pub first_moment_width_mm3: f64,
    pub centroid_width_mm: f64,
    pub maximum_cells: usize,
    /// Fresh per-face work; inherited immutable material is not reconstructed.
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for Degree11FacePropertyOptions {
    fn default() -> Self {
        Self {
            area_width_mm2: 1e-5,
            first_moment_width_mm3: 1e-5,
            centroid_width_mm: 1e-5,
            maximum_cells: 8192,
            maximum_work: END_WORK,
            maximum_bytes: MAX_BYTES,
        }
    }
}
pub struct Degree11FaceMeasurement<'s> {
    source: &'s OwnedDegree11Material,
    pub source_face: usize,
    pub area_mm2: [f64; 2],
    pub first_moment_mm3: [[f64; 2]; 3],
    pub centroid_mm: [[f64; 2]; 3],
    work: usize,
    peak: usize,
}
impl Degree11FaceMeasurement<'_> {
    pub fn source(&self) -> &OwnedDegree11Material {
        self.source
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn accounted_peak_bytes(&self) -> usize {
        self.peak
    }
}
fn iv(b: [f64; 2]) -> I {
    I::new(b[0], b[1])
}
fn pi() -> I {
    I::new(
        std::f64::consts::PI.next_down(),
        std::f64::consts::PI.next_up(),
    )
}
fn fail(e: crate::tapered_thread::Error) -> PhillipsReadError {
    use crate::tapered_thread::Error as E;
    match e {
        E::Invalid => PhillipsReadError::Invalid,
        E::Limit => PhillipsReadError::Limit,
        E::Cancelled => PhillipsReadError::Cancelled,
        _ => PhillipsReadError::Unsupported,
    }
}
fn head_error(e: crate::head_profile::RepairedHeadError) -> PhillipsReadError {
    use crate::head_profile::RepairedHeadError as E;
    match e {
        E::Invalid => PhillipsReadError::Invalid,
        E::Limit => PhillipsReadError::Limit,
        E::Cancelled => PhillipsReadError::Cancelled,
        E::Unresolved => PhillipsReadError::Unsupported,
    }
}
fn cavity_error(e: crate::phillips_drive::Error) -> PhillipsReadError {
    use crate::phillips_drive::Error as E;
    match e {
        E::Invalid => PhillipsReadError::Invalid,
        E::Limit => PhillipsReadError::Limit,
        E::Cancelled => PhillipsReadError::Cancelled,
        E::Unresolved => PhillipsReadError::Unsupported,
    }
}
impl OwnedDegree11Material {
    /// Measure one immutable face inventory index. This keeps an inspector's
    /// selected-face query independent of unrelated curved-face integrations.
    pub fn face_properties(
        &self,
        index: usize,
        o: Degree11FacePropertyOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<Degree11FaceMeasurement<'_>, PhillipsReadError> {
        let source_face = self.face_source(index).ok_or(PhillipsReadError::Invalid)?;
        if [
            o.area_width_mm2,
            o.first_moment_width_mm3,
            o.centroid_width_mm,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.)
            || o.maximum_cells == 0
            || o.maximum_cells > 65536
            || o.maximum_work > END_WORK
            || o.maximum_bytes > MAX_BYTES
        {
            return Err(PhillipsReadError::Invalid);
        }
        let mut work = EndWork {
            remaining: o.maximum_work,
            interrupt,
        };
        work.step(1024)?;
        let resident = self
            .retained_budget_bytes()
            .checked_add(std::mem::size_of::<Degree11FaceMeasurement<'_>>() + 16 * 1024)
            .ok_or(PhillipsReadError::Limit)?;
        if resident > o.maximum_bytes {
            return Err(PhillipsReadError::Limit);
        }
        let mut peak = resident;
        let accuracy = FacePropertyOptions {
            area_tolerance_mm2: o.area_width_mm2 / 16.,
            first_moment_tolerance_mm3: o.first_moment_width_mm3 / 16.,
            centroid_tolerance_mm: o.centroid_width_mm / 16.,
            maximum_cells: o.maximum_cells,
        };
        let mut cavity = |key: FaceKey,
                          work: &mut EndWork<'_>|
         -> std::result::Result<[I; 4], PhillipsReadError> {
            // Existing Cutter API counts cumulative integration cells. Reserve
            // a conservative fixed scalar-operation/workspace envelope per cell.
            let cells = o.maximum_cells.min(8192);
            let scratch = cells
                .checked_mul(4096)
                .and_then(|n| n.checked_add(1024 * 1024))
                .ok_or(PhillipsReadError::Limit)?;
            let need = resident
                .checked_add(scratch)
                .ok_or(PhillipsReadError::Limit)?;
            if need > o.maximum_bytes {
                return Err(PhillipsReadError::Limit);
            }
            peak = peak.max(need);
            work.step(cells.checked_mul(16384).ok_or(PhillipsReadError::Limit)?)?;
            let values = self
                .cavity()
                .local_face_properties_with_cancel(
                    key,
                    FacePropertyOptions {
                        maximum_cells: cells,
                        ..accuracy
                    },
                    &|| interrupt.is_some_and(|f| f()),
                )
                .map_err(cavity_error)?;
            let area = iv(values.area_mm2);
            let z = I::point(self.chart().outer_top_z);
            Ok([
                area,
                iv(values.first_moment_mm3[0]),
                iv(values.first_moment_mm3[1]),
                iv(values.first_moment_mm3[2]) + z * area,
            ])
        };
        let mut certified_centroid = None;
        let values = match source_face {
            Degree11MaterialFace::ThreadSide { index: face, .. } => {
                let tensor = self
                    .thread_faces()
                    .get(face)
                    .ok_or(PhillipsReadError::Unsupported)?
                    .source();
                let shared = std::mem::size_of::<crate::tapered_thread::RetainedTensor11>()
                    + std::mem::size_of_val(tensor.poles())
                    + 64;
                let outside = resident
                    .checked_sub(shared)
                    .ok_or(PhillipsReadError::Limit)?;
                if work.remaining == 0 {
                    return Err(PhillipsReadError::Limit);
                }
                let q = tensor
                    .surface_properties(
                        SurfaceOptions11 {
                            area_width_mm2: o.area_width_mm2,
                            first_moment_width_mm3: o.first_moment_width_mm3,
                            centroid_width_mm: o.centroid_width_mm,
                            max_work: work.remaining,
                            max_bytes: o
                                .maximum_bytes
                                .checked_sub(outside)
                                .ok_or(PhillipsReadError::Limit)?,
                            max_cells: o.maximum_cells.min(8192),
                            max_depth: 48,
                        },
                        interrupt,
                    )
                    .map_err(fail)?;
                work.step(q.work_used)?;
                peak = peak.max(
                    outside
                        .checked_add(q.accounted_peak_bytes)
                        .ok_or(PhillipsReadError::Limit)?,
                );
                certified_centroid = Some(q.centroid_mm.map(iv));
                [
                    iv(q.area_mm2),
                    iv(q.first_moment_mm3[0]),
                    iv(q.first_moment_mm3[1]),
                    iv(q.first_moment_mm3[2]),
                ]
            }
            Degree11MaterialFace::HeadMeridian { segment, .. } => {
                let head = self.head().view();
                let shared = std::mem::size_of_val(head.profile_proposal());
                let outside = resident
                    .checked_sub(shared)
                    .ok_or(PhillipsReadError::Limit)?;
                let q = head
                    .surface_properties(
                        RepairedHeadSurfaceOptions {
                            area_tolerance_mm2: accuracy.area_tolerance_mm2,
                            first_moment_tolerance_mm3: accuracy.first_moment_tolerance_mm3,
                            centroid_tolerance_mm: accuracy.centroid_tolerance_mm,
                            maximum_cells: o.maximum_cells,
                            maximum_work: work.remaining.min(16_000_000),
                            maximum_bytes: o
                                .maximum_bytes
                                .checked_sub(outside)
                                .ok_or(PhillipsReadError::Limit)?,
                        },
                        interrupt,
                    )
                    .map_err(head_error)?;
                work.step(q.work_used())?;
                peak = peak.max(
                    outside
                        .checked_add(q.accounted_peak_bytes())
                        .ok_or(PhillipsReadError::Limit)?,
                );
                let q = q
                    .surfaces()
                    .get(segment)
                    .ok_or(PhillipsReadError::Unsupported)?;
                certified_centroid = Some(q.centroid_mm.map(iv));
                [
                    iv(q.area_mm2),
                    iv(q.first_moment_mm3[0]),
                    iv(q.first_moment_mm3[1]),
                    iv(q.first_moment_mm3[2]),
                ]
            }
            Degree11MaterialFace::ReversedCavity { face, .. } => cavity(face, &mut work)?,
            Degree11MaterialFace::TopAnnulus { .. } => {
                let head = self.head().view();
                let station = head.station(0).ok_or(PhillipsReadError::Unsupported)?;
                let disk = pi() * I::point(station[0]).square();
                let hole = cavity(FaceKey::Mouth, &mut work)?;
                [
                    disk - hole[0],
                    -hole[1],
                    -hole[2],
                    I::point(station[1]) * disk - hole[3],
                ]
            }
            Degree11MaterialFace::ThreadTip { .. } | Degree11MaterialFace::JoinAnnulus { .. } => {
                let q = self.end_disk_properties(
                    Degree11EndDiskOptions {
                        maximum_area_width_mm2: accuracy.area_tolerance_mm2,
                        maximum_first_moment_width_mm3: accuracy.first_moment_tolerance_mm3,
                        maximum_centroid_width_mm: accuracy.centroid_tolerance_mm,
                        maximum_work: work.remaining,
                        maximum_bytes: o
                            .maximum_bytes
                            .checked_sub(resident - self.retained_budget_bytes())
                            .ok_or(PhillipsReadError::Limit)?,
                    },
                    interrupt,
                )?;
                work.step(q.work_used())?;
                peak = peak.max(
                    q.accounted_peak_bytes()
                        .checked_add(resident - self.retained_budget_bytes())
                        .ok_or(PhillipsReadError::Limit)?,
                );
                if matches!(source_face, Degree11MaterialFace::ThreadTip { .. }) {
                    let q = q.ends()[1];
                    [
                        iv(q.area_mm2),
                        iv(q.first_moment_mm3[0]),
                        iv(q.first_moment_mm3[1]),
                        iv(q.first_moment_mm3[2]),
                    ]
                } else {
                    let head = self.head().view();
                    let station = head
                        .station(head.segment_count())
                        .ok_or(PhillipsReadError::Unsupported)?;
                    let inner = q.ends()[0];
                    let area = pi() * I::point(station[0]).square() - iv(inner.area_mm2);
                    [
                        area,
                        -iv(inner.first_moment_mm3[0]),
                        -iv(inner.first_moment_mm3[1]),
                        I::point(station[1]) * area,
                    ]
                }
            }
        };
        let area = values[0];
        let moments = [values[1], values[2], values[3]];
        let ratio = moments.map(|m| m / area);
        // Reuse the same child's centred integral certificate. Re-forming only
        // raw M/A would lose its shared-origin/area correlation unnecessarily.
        let centroid = certified_centroid.map_or(ratio, |c| {
            std::array::from_fn(|i| I::new(ratio[i].lo.max(c[i].lo), ratio[i].hi.min(c[i].hi)))
        });
        if peak > o.maximum_bytes {
            return Err(PhillipsReadError::Limit);
        }
        if !area.valid()
            || area.lo <= 0.
            || width(area) > o.area_width_mm2
            || moments
                .iter()
                .any(|&v| !v.valid() || width(v) > o.first_moment_width_mm3)
            || centroid
                .iter()
                .any(|&v| !v.valid() || width(v) > o.centroid_width_mm)
        {
            return Err(PhillipsReadError::Unsupported);
        }
        work.step(1)?;
        Ok(Degree11FaceMeasurement {
            source: self,
            source_face: index,
            area_mm2: pair(area),
            first_moment_mm3: moments.map(pair),
            centroid_mm: centroid.map(pair),
            work: o.maximum_work - work.remaining,
            peak,
        })
    }
}
