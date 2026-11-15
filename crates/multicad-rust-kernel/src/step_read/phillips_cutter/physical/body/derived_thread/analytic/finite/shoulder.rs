//! Original nominal shared-crest shoulder to a distinct radial annulus.
//! The collapsed crest interval is deliberately omitted from that annulus.
//! This is a join receipt, not a complete material or source-file admission.
use super::*;
use crate::head_profile::RepairedHeadMeridian;
mod domain;
mod material;
mod plane;
pub use material::{
    AuthoredPeriodicOptions, OwnedPeriodicMaterial, PeriodicClosedMaterial, PeriodicCurveOwner,
    PeriodicMaterialFace, PeriodicMaterialOptions, PeriodicMaterialProvenance,
    PeriodicMaterialVolume, PeriodicOwnedOptions, PeriodicVolumeOptions,
};
#[cfg(test)]
mod tests;
const TOTAL: usize = 64_000_000;
const SCRATCH: usize = 128 * 1024;
#[derive(Clone, Copy)]
pub struct PeriodicShoulderOptions {
    pub maximum_distance_mm: f64,
    pub maximum_cell_derivative_mm: f64,
    pub maximum_work: usize,
    pub maximum_total_work: usize,
    pub maximum_bytes: usize,
}
impl Default for PeriodicShoulderOptions {
    fn default() -> Self {
        Self {
            maximum_distance_mm: 1e-7,
            maximum_cell_derivative_mm: 1e-7,
            maximum_work: WORK,
            maximum_total_work: TOTAL,
            maximum_bytes: BYTES,
        }
    }
}
pub struct PeriodicShoulderRepair<'s, 't, 'a, 'input> {
    finite: &'s PeriodicFiniteCurves<'t, 'a, 'input>,
    head: &'s RepairedHeadMeridian<'a>,
    complement: PeriodicFiniteCurve,
    plane_entity: u32,
    plane_distance: f64,
    domain: domain::Domain,
    new_work: usize,
    retained: usize,
}
impl<'s, 't, 'a, 'input> PeriodicShoulderRepair<'s, 't, 'a, 'input> {
    pub fn finite_curves(&self) -> &'s PeriodicFiniteCurves<'t, 'a, 'input> {
        self.finite
    }
    pub fn head_source(&self) -> &'s RepairedHeadMeridian<'a> {
        self.head
    }
    pub fn original_complement(&self) -> &PeriodicFiniteCurve {
        &self.complement
    }
    pub fn original_plane_entity(&self) -> u32 {
        self.plane_entity
    }
    pub fn plane_distance_mm(&self) -> f64 {
        self.plane_distance
    }
    pub fn head_bottom_station_mm(&self) -> [f64; 2] {
        self.domain.bottom
    }
    /// Local radial coordinate of the three noncollapsed annulus bands.
    /// Angular position is the same real phase expression as periodic_source.
    pub fn periodic_source(&self) -> &PeriodicThreadSource {
        &self.domain.source
    }
    pub fn radial_bounds_mm(
        &self,
        band: usize,
        u: [f64; 2],
        radial: [f64; 2],
        maximum_work: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> FResult<[f64; 2]> {
        let mut work = Work {
            remaining: maximum_work,
            interrupt,
            cancelled: Cell::new(false),
        };
        if maximum_work > WORK {
            return Err(PhillipsReadError::Invalid);
        }
        let result = self
            .domain
            .radial(
                band,
                I::new(u[0], u[1]),
                I::new(radial[0], radial[1]),
                &mut work,
            )
            .map(|r| [r.lo, r.hi]);
        if work.cancelled.get() {
            Err(PhillipsReadError::Cancelled)
        } else {
            result
        }
    }
    pub fn new_work_used(&self) -> usize {
        self.new_work
    }
    pub fn inherited_work_used(&self) -> usize {
        self.finite.total_work_used()
    }
    pub fn total_work_used(&self) -> usize {
        self.new_work + self.inherited_work_used()
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.retained
    }
}
impl<'h, 'q, 'r, 'a, 'input> ClosedRepairedHead<'h, 'q, 'r, 'a, 'input> {
    pub fn periodic_shared_crest_shoulder<'s, 't>(
        &'s self,
        finite: &'s PeriodicFiniteCurves<'t, 'a, 'input>,
        o: PeriodicShoulderOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> FResult<PeriodicShoulderRepair<'s, 't, 'a, 'input>> {
        let mut work = Work {
            remaining: o.maximum_work,
            interrupt,
            cancelled: Cell::new(false),
        };
        let result = (|| {
            work.step(1)?;
            if !o.maximum_distance_mm.is_finite()
                || o.maximum_distance_mm <= 0.
                || o.maximum_distance_mm > 1e-7
                || !o.maximum_cell_derivative_mm.is_finite()
                || o.maximum_cell_derivative_mm <= 0.
                || o.maximum_cell_derivative_mm > 1e-7
                || o.maximum_work > WORK
                || o.maximum_total_work > TOTAL
                || o.maximum_bytes > BYTES
            {
                return Err(PhillipsReadError::Invalid);
            }
            if !std::ptr::eq(self.body(), finite.thread.body())
                || self.chart() != finite.thread.chart()
            {
                return Err(PhillipsReadError::Invalid);
            }
            let allowance = o.maximum_work.min(
                o.maximum_total_work
                    .checked_sub(finite.total_work_used())
                    .ok_or(PhillipsReadError::Limit)?,
            );
            work.remaining = allowance.checked_sub(1).ok_or(PhillipsReadError::Limit)?;
            let retained = finite
                .retained_budget_bytes()
                .saturating_add(std::mem::size_of::<PeriodicShoulderRepair<'_, '_, '_, '_>>());
            if retained.saturating_add(SCRATCH) > o.maximum_bytes {
                return Err(PhillipsReadError::Limit);
            }
            work.step(1024)?;
            let bottom = self.virtual_bottom_station();
            let source = finite.thread.source();
            let body = self.body();
            let shell = Shell::new(&body.entities, body.tolerance, true);
            let face = body.material_faces[body.tensor_face_count];
            let (surface, forward, bounds) = helical::face(&shell, face)?;
            if bounds.len() != 1 {
                return Err(PhillipsReadError::Unsupported);
            }
            let domain = domain::Domain::new(*source, bottom, finite.stations, &mut work)?;
            let uses = helical::orient(
                &helical::ring(&shell, bounds[0], true)?,
                !body.thread.spec().left_hand,
            );
            // profiled::junction earned this exact face, plane, single loop and four
            // original edges before the opaque head was constructed. The first three
            // are the same owning entry edges retained by the supplied finite packet.
            assert!(forward);
            assert_eq!(uses.len(), 4);
            assert_eq!(
                uses[..3],
                [
                    (finite.curves[0].edge_entity, true),
                    (finite.curves[4].edge_entity, true),
                    (finite.curves[8].edge_entity, true)
                ]
            );
            assert!(!uses[3].1);
            let edge = uses[3].0;
            let a = shell.args(edge, "EDGE_CURVE")?;
            let vertices = [reference(&a[1])?, reference(&a[2])?];
            assert_eq!(
                vertices,
                [
                    finite.curves[0].vertex_entities[0],
                    finite.curves[8].vertex_entities[1]
                ]
            );
            let vertex_points = vertices
                .map(|v| shell.args(v, "VERTEX_POINT").and_then(|a| reference(&a[1])))
                .into_iter()
                .collect::<Result<Vec<_>>>()?;
            let mut complement = PeriodicFiniteCurve {
                face_entity: face,
                edge_entity: edge,
                curve_entity: reference(&a[3])?,
                role: 3,
                forward_in_face: false,
                curve_same_sense: sense(&a[4])?,
                vertex_entities: vertices,
                vertex_point_entities: [vertex_points[0], vertex_points[1]],
                kind: PeriodicFiniteKind::Circle,
                original_curve_controls: Vec::new(),
                original_surface_controls: Vec::new(),
                trim_control_indices: [0, 1],
                distance_mm: 0.,
                cell_derivative_mm: 0.,
                endpoint_incidence_reserve_mm: 0.,
            };
            let target =
                circle::target::ArcTarget::new(*source, [0, 3], bottom[0], bottom[1], &mut work)?;
            circle::qualify_target(
                body,
                &shell,
                self.chart(),
                target,
                &mut complement,
                &mut work,
            )?;
            let a = shell.args(surface, "PLANE")?;
            let axis = shell.args(reference(&a[1])?, "AXIS2_PLACEMENT_3D")?;
            let origin = point(body, reference(&axis[1])?, &mut work)?;
            let normal = circle::direction(body, &shell, reference(&axis[2])?, &mut work)?;
            let reference = circle::direction(body, &shell, reference(&axis[3])?, &mut work)?;
            let plane_distance =
                plane::disk_transfer(origin, normal, reference, self.chart(), bottom, &mut work)?;
            if plane_distance > o.maximum_distance_mm
                || complement.distance_mm > o.maximum_distance_mm
                || complement.cell_derivative_mm > o.maximum_cell_derivative_mm
            {
                return Err(PhillipsReadError::Unsupported);
            }
            work.step(1)?;
            Ok(PeriodicShoulderRepair {
                finite,
                head: self.source(),
                complement,
                plane_entity: surface,
                plane_distance,
                domain,
                new_work: allowance - work.remaining,
                retained,
            })
        })();
        if work.cancelled.get() {
            Err(PhillipsReadError::Cancelled)
        } else {
            result
        }
    }
}
