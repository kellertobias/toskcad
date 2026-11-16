//! Explicit common-chart carrier transfer and proposed head-envelope containment.
//! Retained meridian support/station gaps remain; no closed outer/material/mass.
use super::*;
use crate::{
    phillips_drive::{HeadEnvelopeClearance, HeadEnvelopeOptions},
    step::helical::prepare_retained_reference,
    tapered_thread::interval::I,
};
#[cfg(test)]
pub(in crate::step_read::phillips_cutter::physical::body) mod tests;
const WORK: usize = 16_000_000;
const BYTES: usize = 64 * 1024 * 1024;
const CELLS: usize = 65536;
#[derive(Clone, Copy)]
pub struct CommonMouthOptions {
    pub maximum_carrier_distance_mm: f64,
    pub minimum_head_clearance_mm: f64,
    pub maximum_cells: usize,
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for CommonMouthOptions {
    fn default() -> Self {
        Self {
            maximum_carrier_distance_mm: 1e-7,
            minimum_head_clearance_mm: 1e-7,
            maximum_cells: CELLS,
            maximum_work: WORK,
            maximum_bytes: BYTES,
        }
    }
}
/// NEW coupled chart, not the original exact plane/outer placement. Source
/// outer points map as O+U*x+V*y+N*(z-topZ). The subtraction is symbolic; no
/// independently rounded translated origin is stored or asserted coincident.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CommonMouthChart {
    pub columns: [Point; 4],
    pub outer_top_z: f64,
}
#[derive(Clone, Copy, Debug, Default)]
pub struct CommonTensorTransfer {
    pub face_entity: u32,
    pub surface_entity: u32,
    pub control_rows: usize,
    pub control_columns: usize,
    pub long_rails: [u32; 2],
    pub entry_vertex_entity: u32,
    pub entry_control_column: usize,
    pub distance_mm: f64,
}
/// Same-Body-bound source/carrier/envelope receipt. The old qualified finite
/// head support/trim matching and intrinsic Profile reconstruction reserve stay
/// explicitly separate from newly measured common-chart residuals. This is NOT
/// exact raw sewing, a closed outer source, material closure, flux or file mass.
pub struct CommonMouthHeadEnvelope<'q, 'r, 'a, 'input> {
    closed: &'q ClosedDerivedCavityCertificate<'r, 'a, 'input>,
    chart: CommonMouthChart,
    tensors: [CommonTensorTransfer; 32],
    count: usize,
    mouth_distance: f64,
    head_frame_distance: f64,
    head: HeadEnvelopeClearance<'r, 'a>,
    head_reserve: f64,
    work: usize,
    retained: usize,
}
impl<'q, 'r, 'a, 'input> CommonMouthHeadEnvelope<'q, 'r, 'a, 'input> {
    pub fn closed_cavity(&self) -> &'q ClosedDerivedCavityCertificate<'r, 'a, 'input> {
        self.closed
    }
    pub fn chart(&self) -> CommonMouthChart {
        self.chart
    }
    pub fn tensor_transfers(&self) -> &[CommonTensorTransfer] {
        &self.tensors[..self.count]
    }
    /// Actual original decimal plane distance over the ENTIRE cavity mouth disk.
    pub fn raw_mouth_plane_distance_mm(&self) -> f64 {
        self.mouth_distance
    }
    /// Whole finite line/arc-support extent transfer between proposal frames.
    pub fn whole_head_frame_transfer_mm(&self) -> f64 {
        self.head_frame_distance
    }
    /// Inherited finite-matching arithmetic reserve, NOT a new measured residual.
    pub fn inherited_head_reconstruction_reserve_mm(&self) -> f64 {
        self.head_reserve
    }
    pub fn proposed_head_envelope(&self) -> &HeadEnvelopeClearance<'r, 'a> {
        &self.head
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.retained
    }
}
fn norm(v: [I; 3]) -> Result<f64> {
    let d = v
        .into_iter()
        .fold(I::point(0.), |s, v| s + v.square())
        .sqrt();
    require(d.valid() && d.hi.is_finite())?;
    Ok(d.hi)
}
fn cross(a: Point, b: Point) -> Point {
    std::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
}
fn mapped(chart: CommonMouthChart, xy: [I; 2], z: I) -> [I; 3] {
    std::array::from_fn(|i| {
        I::point(chart.columns[0][i])
            + I::point(chart.columns[1][i]) * xy[0]
            + I::point(chart.columns[2][i]) * xy[1]
            + I::point(chart.columns[3][i]) * (z - I::point(chart.outer_top_z))
    })
}
fn raw_plane_distance(body: &Body<'_>, chart: CommonMouthChart, radius: f64) -> Result<f64> {
    let shell = Shell::new(&body.entities, body.tolerance, true);
    let f = body
        .cavity
        .faces()
        .iter()
        .find(|f| f.material_entity.is_none())
        .ok_or(ReadError::Invalid)?;
    let raw = geometry::support(&shell, f.surface_entity, f.cavity_forward)?;
    let frame = support::Frame::literal(&raw.axis).map_err(source_error)?;
    let normal = frame.axes[2];
    let delta = std::array::from_fn::<_, 3, _>(|i| I::point(chart.columns[0][i]) - frame.origin[i]);
    let dot = |v: [I; 3]| (0..3).fold(I::point(0.), |s, i| s + normal[i] * v[i]);
    let center = dot(delta).magnitude();
    let du = dot(chart.columns[1].map(I::point)).magnitude();
    let dv = dot(chart.columns[2].map(I::point)).magnitude();
    let bound = I::point(center) + I::point(radius) * (I::point(du) + I::point(dv));
    require(bound.valid())?;
    Ok(bound.hi)
}
fn head_frame_distance(body: &Body<'_>, chart: CommonMouthChart) -> Result<f64> {
    let [o, u, v] = body.frame;
    let n = cross(u, v);
    let mut radius = 0_f64;
    let mut axial = 0_f64;
    for segment in body.profile.segments() {
        radius = radius.max(segment.start[0].abs()).max(segment.end[0].abs());
        axial = axial.max(segment.start[1].abs()).max(segment.end[1].abs());
        if let crate::head_profile::Curve::Arc {
            center, radius_mm, ..
        } = segment.curve
        {
            radius = radius.max((I::point(center[0].abs()) + I::point(radius_mm)).hi);
            axial = axial.max((I::point(center[1].abs()) + I::point(radius_mm)).hi);
        }
    }
    let r = I::new(-radius, radius);
    let z = I::new(-axial, axial);
    let delta = std::array::from_fn(|i| {
        I::point(chart.columns[0][i])
            - I::point(o[i])
            - I::point(chart.columns[3][i]) * I::point(chart.outer_top_z)
            + (I::point(chart.columns[1][i]) - I::point(u[i])) * r
            + (I::point(chart.columns[2][i]) - I::point(v[i])) * r
            + (I::point(chart.columns[3][i]) - I::point(n[i])) * z
    });
    norm(delta)
}
fn pure_error(e: crate::phillips_drive::Error) -> PhillipsReadError {
    if e == crate::phillips_drive::Error::Cancelled {
        PhillipsReadError::Cancelled
    } else {
        cutter_error(e).into()
    }
}
impl<'r, 'a, 'input> ClosedDerivedCavityCertificate<'r, 'a, 'input> {
    pub fn common_mouth_head_envelope<'q>(
        &'q self,
        outer: &OuterCarrierBounds<'_, 'input>,
        options: CommonMouthOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<CommonMouthHeadEnvelope<'q, 'r, 'a, 'input>, PhillipsReadError> {
        let body = self.representative().body();
        let mut work = Work {
            remaining: options.maximum_work,
            interrupt,
            cancelled: Cell::new(false),
        };
        let result = (|| {
            work.step(1)?;
            if !std::ptr::eq(body, outer.body())
                || !options.maximum_carrier_distance_mm.is_finite()
                || options.maximum_carrier_distance_mm <= 0.
                || options.maximum_carrier_distance_mm > 1e-7
                || !options.minimum_head_clearance_mm.is_finite()
                || options.minimum_head_clearance_mm <= 0.
                || options.maximum_cells == 0
                || options.maximum_cells > CELLS
                || options.maximum_work > WORK
                || options.maximum_bytes == 0
                || options.maximum_bytes > BYTES
            {
                return Err(ReadError::Invalid);
            }
            work.step(self.work_used().saturating_add(outer.work_used()))?;
            let points = outer
                .tensors()
                .iter()
                .map(|t| t.control_rows.saturating_mul(t.control_columns))
                .fold(0usize, usize::saturating_add);
            let max_points = outer
                .tensors()
                .iter()
                .map(|t| t.control_rows.saturating_mul(t.control_columns))
                .max()
                .unwrap_or(0);
            // Actual Body once; existing owned outer receipt payload once. Fresh Vec
            // controls/reference factors include fourfold growth/reallocation reserve.
            let scratch = max_points
                .saturating_mul(4 * std::mem::size_of::<Point>())
                .saturating_add(points.saturating_mul(128))
                .saturating_add(HeadEnvelopeClearance::scratch_budget_bytes(
                    options.maximum_cells,
                ))
                .saturating_add(1024 * 1024);
            let retained = self
                .retained_budget_bytes()
                .saturating_add(std::mem::size_of::<CommonMouthHeadEnvelope<'_, '_, '_, '_>>());
            if retained
                .saturating_add(outer.retained_capacity_bytes())
                .saturating_add(scratch)
                > options.maximum_bytes
            {
                return Err(ReadError::Limit);
            }
            work.step(points.saturating_mul(256).saturating_add(4096))?;
            let basis = prepare_retained_reference(&body.thread, body.visible[0], body.tolerance)
                .map_err(helical::step_error)?;
            work.step(1)?;
            let chart = CommonMouthChart {
                columns: self.actual_affine_columns(),
                outer_top_z: body.profile.top()[1],
            };
            let radius = self
                .local_source_certificate()
                .source()
                .envelope()
                .entry_radius_mm_bounds()[1];
            let mouth = raw_plane_distance(body, chart, radius)?;
            let head_frame = head_frame_distance(body, chart)?;
            require(
                mouth <= options.maximum_carrier_distance_mm
                    && head_frame <= options.maximum_carrier_distance_mm,
            )?;
            let mut transfers = [CommonTensorTransfer::default(); 32];
            // Exact immutable OuterCarrierBounds/Body binding and identical deterministic
            // reference inputs already earned full sizes/knots/owner/finite-trim matching.
            for (i, t) in outer.tensors().iter().enumerate() {
                work.step(1)?;
                let actual = body.literal_tensor_control_points_raw(i)?;
                let mut distance = 0_f64;
                for (row, r) in actual.iter().zip(&basis.bands[i].poles) {
                    for ((p, a), z) in row.iter().zip(&basis.axial.poles).zip(&basis.heights) {
                        work.step(1)?;
                        let r = r.map(I::point);
                        let a = a.map(I::point);
                        let xy = [r[0] * a[0] - r[1] * a[1], r[0] * a[1] + r[1] * a[0]];
                        let expected = mapped(chart, xy, I::point(*z));
                        let residual =
                            std::array::from_fn(|k| support::adjacent(p[k]) - expected[k]);
                        let deviation = norm(residual)?;
                        distance = distance.max(deviation);
                    }
                }
                require(distance <= options.maximum_carrier_distance_mm)?;
                transfers[i] = CommonTensorTransfer {
                    face_entity: t.face_entity,
                    surface_entity: t.surface_entity,
                    control_rows: t.control_rows,
                    control_columns: t.control_columns,
                    long_rails: t.long_rails,
                    entry_vertex_entity: t.entry_vertex_entity,
                    entry_control_column: t.entry_control_column,
                    distance_mm: distance,
                };
            }
            work.step(1)?;
            Ok((chart, transfers, mouth, head_frame, retained))
        })();
        let (chart, tensors, mouth_distance, head_frame_distance, retained) =
            result.map_err(|e| {
                if work.cancelled.get() {
                    PhillipsReadError::Cancelled
                } else {
                    e.into()
                }
            })?;
        let inherited = self.local_source_certificate().work_used();
        let head = self
            .local_source_certificate()
            .certify_head_envelope(
                &body.profile,
                HeadEnvelopeOptions {
                    minimum_clearance_mm: options.minimum_head_clearance_mm,
                    maximum_cells: options.maximum_cells,
                    maximum_work: work.remaining.saturating_add(inherited),
                    maximum_bytes: options.maximum_bytes,
                },
                interrupt,
            )
            .map_err(pure_error)?;
        let used = options.maximum_work - work.remaining + head.work_used() - inherited;
        if interrupt.is_some_and(|f| f()) {
            return Err(PhillipsReadError::Cancelled);
        }
        Ok(CommonMouthHeadEnvelope {
            closed: self,
            chart,
            tensors,
            count: outer.tensors().len(),
            mouth_distance,
            head_frame_distance,
            head,
            head_reserve: body.profile.rounding_error_mm(),
            work: used,
            retained,
        })
    }
}
