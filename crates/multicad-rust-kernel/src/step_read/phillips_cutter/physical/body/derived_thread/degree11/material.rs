//! Closed DISTINCT head/thread-minus-cavity material; no mass or file admission.
//!
//! The occupied D11 thread starts on the exact head-bottom station. The head's
//! monotone atlas and thread's strict axial graph put their interiors in
//! opposite halfspaces. A strict full-rim gap makes the radial homotopy from
//! the SAME D11 entry curve to the head circle an embedded planar annulus.
//! Removing both virtual join disks and retaining the actual tip disk closes
//! their union. The SAME head's certified cavity meets that union only at its
//! shaped mouth: replace the top disk by the actual outer-minus-mouth annulus
//! and reverse every cavity side. The complete finite cavity atlas supplies
//! all faces and mouth curves; no circular substitute or fitted thread is used.
use super::*;
use crate::{
    head_profile::RepairedHeadMeridian,
    phillips_drive::{
        ClosedBoundaryCertificate, Cutter, FaceKey, RepairedHeadContainment,
        RepairedHeadContainmentOptions,
    },
};
const MAX_WORK: usize = 64_000_000;
const MAX_BYTES: usize = 64 * 1024 * 1024;
const MAX_CELLS: usize = 65_536;
const MAX_DEPTH: usize = 24;
const RIM_SCRATCH: usize = 64 * 1024;

#[derive(Clone, Copy)]
pub struct Degree11MaterialOptions {
    /// Local radial/axial margin, not affine-world Euclidean wall thickness.
    pub minimum_clearance_mm: f64,
    /// Whole request, including the occupied carrier aggregate exactly once.
    pub maximum_work: usize,
    pub maximum_bytes: usize,
    pub maximum_rim_cells: usize,
    pub maximum_rim_depth: usize,
    pub maximum_containment_cells: usize,
}
impl Default for Degree11MaterialOptions {
    fn default() -> Self {
        Self {
            minimum_clearance_mm: 1e-7,
            maximum_work: MAX_WORK,
            maximum_bytes: MAX_BYTES,
            maximum_rim_cells: MAX_CELLS,
            maximum_rim_depth: MAX_DEPTH,
            maximum_containment_cells: MAX_CELLS,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Degree11JoinCollar {
    pub inherited_sewing_world_mm: f64,
    pub head_world_mm: f64,
    pub cavity_world_mm: f64,
    pub thread_world_mm: f64,
    pub head_intrinsic_local_mm: f64,
    /// Actual inverse-chart entrywise absolute sum: world componentwise error
    /// to local L1 error. Also bounds local radial and axial displacement.
    pub inverse_component_to_local_l1: f64,
    pub local_total_mm: f64,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Degree11MaterialFace {
    HeadMeridian {
        segment: usize,
        raw_material_face: u32,
    },
    ThreadSide {
        index: usize,
        raw_material_face: u32,
    },
    ThreadTip {
        raw_material_face: u32,
    },
    JoinAnnulus {
        raw_shoulder_face: u32,
    },
    TopAnnulus {
        raw_material_face: u32,
    },
    ReversedCavity {
        face: FaceKey,
        raw_material_face: u32,
    },
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Degree11RimUse {
    pub face_index: usize,
    /// Reverse this face's complete original U curve on the entry station.
    pub reversed: bool,
}
pub struct Degree11ClosedMaterial<'s, 'a, 'input> {
    head: &'s RepairedHeadMeridian<'s>,
    thread: &'s OccupiedDegree11Thread<'s, 'a, 'input>,
    cavity: &'s ClosedBoundaryCertificate<'s>,
    determinant: [f64; 2],
    containment: RepairedHeadContainment<'s, 's, 's>,
    collar: Degree11JoinCollar,
    entry_radius: f64,
    join_margin: f64,
    cavity_faces: [(FaceKey, u32); 88],
    cavity_count: usize,
    new_work: usize,
    retained: usize,
    peak: usize,
    rim_cells: usize,
}
impl<'s, 'a, 'input> Degree11ClosedMaterial<'s, 'a, 'input> {
    pub fn body(&self) -> &'a Body<'input> {
        self.thread.body()
    }
    pub fn head(&self) -> &'s RepairedHeadMeridian<'s> {
        self.head
    }
    pub fn thread(&self) -> &'s OccupiedDegree11Thread<'s, 'a, 'input> {
        self.thread
    }
    pub fn cavity(&self) -> &'s Cutter {
        self.cavity.source()
    }
    pub fn cavity_certificate(&self) -> &'s ClosedBoundaryCertificate<'s> {
        self.cavity
    }
    pub fn chart(&self) -> CommonMouthChart {
        self.thread.chart()
    }
    pub fn determinant_bounds(&self) -> [f64; 2] {
        self.determinant
    }
    pub fn containment(&self) -> &RepairedHeadContainment<'s, 's, 's> {
        &self.containment
    }
    pub fn transfer_collar(&self) -> Degree11JoinCollar {
        self.collar
    }
    pub fn entry_radius_upper_mm(&self) -> f64 {
        self.entry_radius
    }
    pub fn join_radial_margin_mm(&self) -> f64 {
        self.join_margin
    }
    pub fn rim_cells(&self) -> usize {
        self.rim_cells
    }
    pub fn new_work_used(&self) -> usize {
        self.new_work
    }
    pub fn work_used(&self) -> usize {
        self.thread.total_work_used() + self.new_work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.retained
    }
    pub fn accounted_peak_bytes(&self) -> usize {
        self.peak
    }
    pub fn face_count(&self) -> usize {
        self.head.segment_count() + self.thread.carrier().faces().len() + 3 + self.cavity_count
    }
    pub fn face_source(&self, mut index: usize) -> Option<Degree11MaterialFace> {
        let body = self.body();
        let h = self.head.segment_count();
        if index < h {
            return Some(Degree11MaterialFace::HeadMeridian {
                segment: index,
                raw_material_face: body.profile_faces[index],
            });
        }
        index -= h;
        let faces = self.thread.carrier().faces();
        if index < faces.len() {
            return Some(Degree11MaterialFace::ThreadSide {
                index,
                raw_material_face: faces[index].transfer().face_entity,
            });
        }
        index -= faces.len();
        match index {
            0 => Some(Degree11MaterialFace::ThreadTip {
                raw_material_face: body.material_faces[body.tensor_face_count + 1],
            }),
            1 => Some(Degree11MaterialFace::JoinAnnulus {
                raw_shoulder_face: body.material_faces[body.tensor_face_count],
            }),
            2 => Some(Degree11MaterialFace::TopAnnulus {
                raw_material_face: body.cap.material_face,
            }),
            _ => self.cavity_faces[..self.cavity_count].get(index - 3).map(
                |&(face, raw_material_face)| Degree11MaterialFace::ReversedCavity {
                    face,
                    raw_material_face,
                },
            ),
        }
    }
    /// Complete inner wire, in the removed entry disk's actual winding.
    /// These uses refer to the same immutable D11 face U curves, not new edges.
    pub fn join_inner_uses(&self) -> impl ExactSizeIterator<Item = Degree11RimUse> + '_ {
        let faces = self.thread.carrier().faces();
        let count = faces.len();
        let reverse = faces[0].forward();
        (0..count).map(move |i| Degree11RimUse {
            face_index: if reverse { count - 1 - i } else { i },
            reversed: reverse,
        })
    }
    pub fn entry_rim_poles(
        &self,
        index: usize,
    ) -> Option<impl ExactSizeIterator<Item = Point> + '_> {
        let source = self.thread.carrier().faces().get(index)?.source();
        let [width, height] = source.dimensions();
        Some((0..width).map(move |u| source.poles()[u * height]))
    }
}
fn valid(o: Degree11MaterialOptions) -> bool {
    o.minimum_clearance_mm.is_finite()
        && o.minimum_clearance_mm > 0.
        && o.maximum_work > 0
        && o.maximum_work <= MAX_WORK
        && o.maximum_bytes > 0
        && o.maximum_bytes <= MAX_BYTES
        && o.maximum_rim_cells > 0
        && o.maximum_rim_cells <= MAX_CELLS
        && o.maximum_rim_depth <= MAX_DEPTH
        && o.maximum_containment_cells > 0
        && o.maximum_containment_cells <= MAX_CELLS
}
fn cross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    std::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
}
fn dot(a: [I; 3], b: [I; 3]) -> I {
    a.into_iter()
        .zip(b)
        .fold(I::point(0.), |s, (a, b)| s + a * b)
}
fn inverse_bound(chart: CommonMouthChart) -> Result<f64> {
    let [_, u, v, n] = chart.columns.map(|p| p.map(I::point));
    let determinant = dot(u, cross(v, n));
    if !determinant.valid() || determinant.lo <= 0. {
        return Err(ReadError::Unsupported);
    }
    let rows = [cross(v, n), cross(n, u), cross(u, v)].map(|r| r.map(|x| x / determinant));
    let norm = rows
        .into_iter()
        .flatten()
        .fold(I::point(0.), |s, x| s + I::point(x.magnitude()));
    require(norm.valid() && norm.hi.is_finite() && norm.hi > 0.)?;
    Ok(norm.hi)
}
fn collar(
    head: &ClosedRepairedHead<'_, '_, '_, '_, '_>,
    thread: &OccupiedDegree11Thread<'_, '_, '_>,
    work: &mut Work<'_>,
) -> Result<Degree11JoinCollar> {
    work.step(4096)?;
    let common = head.common_mouth();
    let rep = common.closed_cavity().representative();
    let t = rep.transfer_bounds();
    let cavity = t
        .maximum_vertex_error_mm()
        .max(t.maximum_analytic_support_error_mm())
        .max(t.maximum_polynomial_error_mm())
        .max(t.maximum_polynomial_support_error_mm());
    let support = head
        .support_transfers()
        .map(|x| x.distance_mm)
        .fold(0f64, f64::max);
    let h = (I::point(support)
        + I::point(common.whole_head_frame_transfer_mm())
        + I::point(common.raw_mouth_plane_distance_mm()))
    .hi;
    let p = thread
        .carrier()
        .faces()
        .iter()
        .map(|f| f.transfer().distance_mm)
        .fold(0f64, f64::max);
    let sewing = rep.inherited_sewing_collar_mm();
    let repaired = (0..head.source().segment_count())
        .map(|i| {
            let c = head
                .source()
                .segment_certificate(i)
                .expect("qualified repaired head segment");
            (I::point(c.support_transfer_mm) + I::point(c.finite_trim_transfer_mm)).hi
        })
        .fold(0f64, f64::max);
    let intrinsic =
        (I::point(head.inherited_finite_matching_admission_mm()) + I::point(repaired)).hi;
    let inverse = inverse_bound(head.chart())?;
    let total = compose_collar([sewing, h, cavity, p, intrinsic], inverse)?;
    Ok(Degree11JoinCollar {
        inherited_sewing_world_mm: sewing,
        head_world_mm: h,
        cavity_world_mm: cavity,
        thread_world_mm: p,
        head_intrinsic_local_mm: intrinsic,
        inverse_component_to_local_l1: inverse,
        local_total_mm: total,
    })
}
fn compose_collar(bounds: [f64; 5], inverse: f64) -> Result<f64> {
    require(bounds.iter().all(|v| v.is_finite() && *v >= 0.))?;
    require(inverse.is_finite() && inverse > 0.)?;
    let [sewing, head, cavity, thread, intrinsic] = bounds;
    let total = (I::point(intrinsic)
        + I::point(inverse)
            * (I::point(sewing) + I::point(head) + I::point(cavity) + I::point(thread)))
    .hi;
    require(total.is_finite() && total >= 0.)?;
    Ok(total)
}
type Curve = [[I; 2]; 12];
#[derive(Clone, Copy)]
struct Node {
    poles: Curve,
    depth: usize,
}
fn split(mut p: Curve) -> [Curve; 2] {
    let mut a = [[I::point(0.); 2]; 12];
    let mut b = a;
    a[0] = p[0];
    b[11] = p[11];
    for n in (1..12).rev() {
        for i in 0..n {
            p[i] = std::array::from_fn(|j| (p[i][j] + p[i + 1][j]) * I::point(0.5));
        }
        a[12 - n] = p[0];
        b[n - 1] = p[n - 1];
    }
    [a, b]
}
fn radius(p: [I; 2]) -> I {
    (p[0].square() + p[1].square()).sqrt()
}
fn rim_bound(
    thread: &OccupiedDegree11Thread<'_, '_, '_>,
    head_radius: f64,
    required: f64,
    o: Degree11MaterialOptions,
    work: &mut Work<'_>,
) -> Result<(f64, usize)> {
    // Fixed26-node stack plus all curve/split/output temporaries <16KiB;
    // caller reserves64KiB before invoking this query. No heap data is built.
    work.step(MAX_DEPTH + 2)?;
    let mut stack: [Option<Node>; MAX_DEPTH + 2] = [None; MAX_DEPTH + 2];
    let mut upper = 0f64;
    let mut cells = 0;
    for face in thread.carrier().faces() {
        let source = face.source();
        let height = source.dimensions()[1];
        for u in 0..source.spans()[0] {
            work.step(128)?;
            let poles = std::array::from_fn(|i| {
                let p = source.poles()[(11 * u + i) * height];
                [I::point(p[0]), I::point(p[1])]
            });
            stack[0] = Some(Node { poles, depth: 0 });
            let mut count = 1;
            while count > 0 {
                count -= 1;
                let node = stack[count].take().expect("bounded D11 rim DFS");
                work.step(512)?;
                cells += 1;
                if cells > o.maximum_rim_cells {
                    return Err(ReadError::Limit);
                }
                let r = node.poles.map(radius);
                let bound = r.iter().map(|r| r.hi).fold(0f64, f64::max);
                if (I::point(head_radius) - I::point(bound)).lo > required {
                    upper = upper.max(bound);
                    continue;
                }
                // Actual endpoint radii give an honest no-clearance refusal.
                if (I::point(head_radius) - I::point(r[0].lo.max(r[11].lo))).hi <= required {
                    return Err(ReadError::Unsupported);
                }
                if node.depth >= o.maximum_rim_depth {
                    return Err(ReadError::Unsupported);
                }
                work.step(2048)?;
                let [a, b] = split(node.poles);
                stack[count] = Some(Node {
                    poles: b,
                    depth: node.depth + 1,
                });
                stack[count + 1] = Some(Node {
                    poles: a,
                    depth: node.depth + 1,
                });
                count += 2;
            }
        }
    }
    Ok((upper, cells))
}
fn cavity_error(e: crate::phillips_drive::Error) -> PhillipsReadError {
    match e {
        crate::phillips_drive::Error::Cancelled => PhillipsReadError::Cancelled,
        crate::phillips_drive::Error::Limit => PhillipsReadError::Limit,
        crate::phillips_drive::Error::Invalid | crate::phillips_drive::Error::Unresolved => {
            PhillipsReadError::Unsupported
        }
    }
}
fn compatible(
    origin: &ClosedRepairedHead<'_, '_, '_, '_, '_>,
    body: &Body<'_>,
    profile: &crate::head_profile::Profile,
    chart: CommonMouthChart,
    station: f64,
    thread: &OccupiedDegree11Thread<'_, '_, '_>,
) -> bool {
    std::ptr::eq(origin, thread.carrier().originating_head())
        && std::ptr::eq(body, thread.body())
        && chart == thread.chart()
        && std::ptr::eq(profile, thread.carrier().original_head_profile())
        && station == thread.cap_stations_mm()[0]
}
impl<'h, 'q, 'r, 'a, 'input> ClosedRepairedHead<'h, 'q, 'r, 'a, 'input> {
    pub fn closed_degree11_material<'s>(
        &'s self,
        thread: &'s OccupiedDegree11Thread<'s, 'a, 'input>,
        o: Degree11MaterialOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<Degree11ClosedMaterial<'s, 'a, 'input>, PhillipsReadError>
    where
        'a: 's,
        'r: 's,
    {
        if !valid(o) {
            return Err(PhillipsReadError::Invalid);
        }
        let available = o
            .maximum_work
            .checked_sub(thread.total_work_used())
            .ok_or(PhillipsReadError::Limit)?;
        let mut work = Work {
            remaining: available,
            interrupt,
            cancelled: Cell::new(false),
        };
        let result = (|| {
            work.step(1)?;
            if !compatible(
                self,
                self.body(),
                self.source().profile_proposal(),
                self.chart(),
                self.virtual_bottom_station()[1],
                thread,
            ) {
                return Err(ReadError::Unsupported);
            }
            let retained = thread
                .retained_budget_bytes()
                .checked_add(std::mem::size_of::<Degree11ClosedMaterial<'_, '_, '_>>())
                .ok_or(ReadError::Limit)?;
            let scratch =
                RepairedHeadContainment::scratch_budget_bytes(o.maximum_containment_cells);
            let peak = retained
                .checked_add(scratch)
                .and_then(|x| x.checked_add(RIM_SCRATCH))
                .ok_or(ReadError::Limit)?;
            if peak > o.maximum_bytes {
                return Err(ReadError::Limit);
            }
            let collar = collar(self, thread, &mut work)?;
            let required = (I::point(o.minimum_clearance_mm) + I::point(collar.local_total_mm)).hi;
            let (entry_radius, rim_cells) = rim_bound(
                thread,
                self.virtual_bottom_station()[0],
                required,
                o,
                &mut work,
            )?;
            let join_margin = (I::point(self.virtual_bottom_station()[0])
                - I::point(entry_radius)
                - I::point(collar.local_total_mm))
            .lo;
            require(join_margin > o.minimum_clearance_mm)?;
            Ok((
                retained,
                peak,
                collar,
                required,
                entry_radius,
                rim_cells,
                join_margin,
            ))
        })();
        let (retained, peak, collar, required, entry_radius, rim_cells, join_margin) =
            result.map_err(|e| bound_error(e, &work))?;
        let closed = self.common_mouth().closed_cavity();
        let local = closed.local_source_certificate();
        let head: &'s RepairedHeadMeridian<'s> = self.source();
        // Existing containment reports its SAME local-cavity and head-source
        // prefix again. That prefix is already in occupied.total_work_used.
        // Pass only prefix+available fresh budget, then debit the checked
        // difference. No old stage is renewed or hidden in this aggregate.
        let prefix = local
            .work_used()
            .checked_add(head.work_used())
            .ok_or(PhillipsReadError::Limit)?;
        let limit = prefix
            .checked_add(work.remaining)
            .ok_or(PhillipsReadError::Limit)?
            .min(16_000_000);
        let containment = local
            .certify_repaired_head(
                head,
                RepairedHeadContainmentOptions {
                    minimum_clearance_mm: required,
                    maximum_cells: o.maximum_containment_cells,
                    maximum_work: limit,
                    maximum_bytes: o.maximum_bytes,
                },
                interrupt,
            )
            .map_err(cavity_error)?;
        let fresh = containment
            .work_used()
            .checked_sub(prefix)
            .ok_or(PhillipsReadError::Unsupported)?;
        work.step(fresh).map_err(|e| bound_error(e, &work))?;
        // Fixed88 source-key/entity pairs, charged before initialization.
        work.step(88 * 2).map_err(|e| bound_error(e, &work))?;
        let mut cavity_faces = [(FaceKey::Mouth, 0); 88];
        let mut count = 0;
        for domain in local.patch_domains() {
            work.step(128).map_err(|e| bound_error(e, &work))?;
            let face = domain.source_face;
            if face != FaceKey::Mouth && !cavity_faces[..count].iter().any(|&(f, _)| f == face) {
                let raw = closed
                    .representative()
                    .faces()
                    .iter()
                    .find(|f| f.source_role == face)
                    .and_then(|f| f.raw_material_entity)
                    .ok_or(PhillipsReadError::Unsupported)?;
                cavity_faces[count] = (face, raw);
                count += 1;
            }
        }
        work.step(1).map_err(|e| bound_error(e, &work))?;
        Ok(Degree11ClosedMaterial {
            head,
            thread,
            cavity: local,
            determinant: closed.determinant_bounds(),
            containment,
            collar,
            entry_radius,
            join_margin,
            cavity_faces,
            cavity_count: count,
            new_work: available - work.remaining,
            retained,
            peak,
            rim_cells,
        })
    }
}
#[cfg(test)]
mod tests;

mod volume;
pub use volume::{
    Degree11EndDiskOptions, Degree11EndDiskProperties, Degree11EndDisks, Degree11FaceMeasurement,
    Degree11FacePropertyOptions, Degree11MaterialQuery, Degree11MaterialQueryError,
    Degree11MaterialQueryOptions, Degree11MaterialVolume, Degree11MaterialVolumeOptions,
    Degree11OwnedMaterialOptions, Degree11ReplayMatch, Degree11ReplayMatchOptions,
    OwnedDegree11Material, SuppliedDegree11MaterialMatch, SuppliedDegree11MaterialOptions,
};
