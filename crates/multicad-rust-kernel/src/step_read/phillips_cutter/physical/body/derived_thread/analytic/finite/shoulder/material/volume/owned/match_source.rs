//! Actual closed-source boundary charts, including cap annuli. Both operands
//! independently own material/void/embedding/mass proofs. This compares those
//! distinct repaired sources; it does not upgrade raw STEP to exact material.
use super::*;
mod world;
#[derive(Clone, Copy)]
pub struct PeriodicMaterialMatchOptions {
    pub maximum_distance_mm: f64,
    pub maximum_derivative_mm: f64,
    pub maximum_volume_discrepancy_mm3: f64,
    /// Fixed stages, including delegated head/thread proof work and the two
    /// bounded canonical atlas builders. This EXCLUDES adaptive Taylor cells;
    /// those have their independently enforced cap below. No aggregate
    /// sixteen-million primitive-operation claim is made by this API.
    pub maximum_fixed_work: usize,
    /// Paired order-eight Taylor attempts, including rejected parent cells.
    /// Each attempt has four fixed Series<9> evaluations, a fixed-size lift,
    /// and a cancellation checkpoint. Depth-first scratch holds <=64 cells.
    /// This is a cell quota, never converted into guessed primitive work.
    pub maximum_cavity_cells: usize,
}
impl Default for PeriodicMaterialMatchOptions {
    fn default() -> Self {
        Self {
            maximum_distance_mm: 1e-7,
            maximum_derivative_mm: 1e-7,
            maximum_volume_discrepancy_mm3: 1e-5,
            maximum_fixed_work: 16_000_000,
            maximum_cavity_cells: 262144,
        }
    }
}
pub struct PeriodicMaterialMatch<'s, 't> {
    source: &'s OwnedPeriodicMaterial,
    target: &'t OwnedPeriodicMaterial,
    distance: f64,
    derivative: f64,
    volume_discrepancy: f64,
    work: usize,
    cavity_cells: usize,
}
impl PeriodicMaterialMatch<'_, '_> {
    pub fn source(&self) -> &OwnedPeriodicMaterial {
        self.source
    }
    pub fn target(&self) -> &OwnedPeriodicMaterial {
        self.target
    }
    pub fn maximum_world_distance_mm(&self) -> f64 {
        self.distance
    }
    /// Same normalized chart coordinates; NOT a unit-normal angle bound.
    pub fn maximum_world_chart_derivative_mm(&self) -> f64 {
        self.derivative
    }
    pub fn independent_volume_discrepancy_mm3(&self) -> f64 {
        self.volume_discrepancy
    }
    /// Consumed fixed-stage work; adaptive cell attempts are reported separately.
    pub fn fixed_work_used(&self) -> usize {
        self.work
    }
    pub fn cavity_cells_used(&self) -> usize {
        self.cavity_cells
    }
}
fn compatible(ok: bool) -> FResult<()> {
    if ok {
        Ok(())
    } else {
        Err(PhillipsReadError::Unsupported)
    }
}
fn abs(x: I) -> I {
    I::point(x.magnitude())
}
fn maximum(x: impl IntoIterator<Item = I>) -> I {
    I::point(x.into_iter().map(|x| x.hi).fold(0., f64::max))
}
fn scalar_delta(a: f64, b: f64) -> I {
    if a == b {
        I::point(0.)
    } else {
        abs(I::point(a) - I::point(b))
    }
}
fn expression_delta(a: [f64; 2], b: [f64; 2]) -> I {
    if a == b {
        I::point(0.)
    } else {
        abs(expression(a) - expression(b))
    }
}
impl OwnedPeriodicMaterial {
    pub fn match_material<'s, 't>(
        &'s self,
        target: &'t Self,
        occurrence: [[f64; 3]; 4],
        options: PeriodicMaterialMatchOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> FResult<PeriodicMaterialMatch<'s, 't>> {
        let mut work = Work {
            remaining: options.maximum_fixed_work,
            interrupt,
            cancelled: Cell::new(false),
        };
        let result = (|| {
            work.step(1)?;
            if [options.maximum_distance_mm, options.maximum_derivative_mm]
                .into_iter()
                .any(|v| !v.is_finite() || v <= 0. || v > 1e-7)
                || !options.maximum_volume_discrepancy_mm3.is_finite()
                || options.maximum_volume_discrepancy_mm3 <= 0.
                || options.maximum_fixed_work > 16_000_000
                || options.maximum_cavity_cells == 0
                || options.maximum_cavity_cells > 262144
            {
                return Err(PhillipsReadError::Invalid);
            }
            compatible(
                self.face_count() == target.face_count()
                    && (0..self.face_count()).all(|i| self.face_source(i) == target.face_source(i)),
            )?;
            compatible(self.thread.left_hand() == target.thread.left_hand())?;
            work.step(8192)?;
            let world = world::compare(self.chart, target.chart, occurrence)?;
            let head = self.head.view();
            let other_head = target.head.view();
            let h = head
                .translated_correspondence(
                    &other_head,
                    self.head_shift,
                    target.head_shift,
                    work.remaining.min(16_000_000),
                    interrupt,
                )
                .map_err(
                    crate::step_read::phillips_cutter::physical::body::closed_head::repaired_error,
                )?;
            work.step(h.work_used)?;
            let t = self.thread.finite_correspondence(
                &target.thread,
                self.caps,
                target.caps,
                work.remaining.min(16_000_000),
                interrupt,
            )?;
            work.step(t.work_used)?;
            // Fixed canonical boundary builders: <=128 nodes,256 edges,88 charts;
            // their primitive construction/validation fits below one million steps.
            work.step(1_000_000)?;
            let a = self
                .cavity
                .surface_patches_with_cancel(&|| interrupt.is_some_and(|f| f()))
                .map_err(cavity_error)?;
            work.step(1_000_000)?;
            let b = target
                .cavity
                .surface_patches_with_cancel(&|| interrupt.is_some_and(|f| f()))
                .map_err(cavity_error)?;
            // Tighten the trim derivative proof before lifting normalized mouth rays.
            // The final complete cap/world bounds below still enforce caller limits.
            let c = a
                .source_correspondence(
                    &b,
                    crate::phillips_drive::PatchCorrespondenceOptions {
                        maximum_position_mm: 1e-10,
                        maximum_derivative_mm: 1e-10,
                        maximum_cells: options.maximum_cavity_cells,
                    },
                    interrupt,
                )
                .map_err(cavity_error)?;
            work.step(8192)?;
            let mouth_shift = expression_delta(self.mouth, target.mouth);
            let source_top = head.station(0).expect("qualified top station");
            let target_top = other_head.station(0).expect("qualified top station");
            let source_bottom = self.bottom;
            let target_bottom = target.bottom;
            let minimum =
                I::point(a.radial_margin_mm_bounds()[0].min(b.radial_margin_mm_bounds()[0]));
            compatible(minimum.lo > 0.)?;
            // q=p/|p|; projectors have norm<=1 and their difference<=2|dq|.
            let mouth_error = I::point(c.position_mm);
            let mouth_derivative = I::point(c.du_mm);
            let source_speed = I::point(c.maximum_source_du_mm);
            let target_speed = I::point(c.maximum_target_du_mm);
            let direction = I::point(2.) * mouth_error / minimum;
            let direction_derivative = mouth_derivative / minimum
                + target_speed * mouth_error / minimum.square()
                + I::point(2.) * direction * target_speed / minimum;
            let top_radius_error = scalar_delta(source_top[0], target_top[0]);
            let top_outer = top_radius_error + I::point(target_top[0]) * direction;
            let top_outer_derivative = top_radius_error * source_speed / minimum
                + I::point(target_top[0]) * direction_derivative;
            // Each annulus chart is the affine radial blend between its actual mouth
            // point and the outer head circle on THAT SAME mouth ray. Positive source
            // clearances make these the complete physical top disk minus the void.
            let top_distance = maximum([mouth_error, top_outer]) + mouth_shift;
            let top_du = maximum([mouth_derivative, top_outer_derivative]);
            let top_dv = mouth_error + top_outer;
            let join_radius_delta = scalar_delta(source_bottom[0], target_bottom[0]);
            let radial = maximum([I::point(t.maximum_radius_error_mm), join_radius_delta]);
            let target_outer = I::point(target_bottom[0].max(target.thread.radius_nodes_mm()[0]));
            let source_outer = I::point(source_bottom[0].max(self.thread.radius_nodes_mm()[0]));
            let angle = I::point(t.maximum_angle_error_radians);
            let join_z = expression_delta(self.caps[0], target.caps[0]);
            let join_distance = radial + target_outer * angle + join_z;
            let join_du = I::point(t.maximum_du_radius_error_mm)
                + I::point(t.maximum_target_du_radius_mm) * angle
                + radial * I::point(t.maximum_source_du_angle_radians)
                + target_outer
                    * (I::point(t.maximum_du_angle_error_radians)
                        + I::point(t.maximum_target_du_angle_radians) * angle);
            let join_dv = I::point(t.maximum_radius_error_mm)
                + join_radius_delta
                + (I::point(target_bottom[0]) + I::point(target.thread.radius_nodes_mm()[0]))
                    * angle;
            let tau = I::point(2.) * rotation::Rotation::new(&mut work)?.pi(&mut work)?;
            // Disk radial maps contract their corresponding rim differences; revolution
            // angular derivatives contribute 2pi*radial C0, independently of meridian C1.
            let local_distance = maximum([
                I::point(h.position_mm),
                I::point(t.position_mm),
                I::point(c.position_mm) + mouth_shift,
                top_distance,
                join_distance,
            ]);
            let local_derivative = maximum([
                I::point(h.derivative_mm),
                tau * I::point(h.position_mm),
                I::point(h.position_mm),
                I::point(t.du_mm),
                I::point(t.dv_mm),
                I::point(t.position_mm),
                I::point(c.du_mm),
                I::point(c.dv_mm),
                top_du,
                top_dv,
                join_du,
                join_dv,
            ]);
            let mut source_derivative = maximum([
                I::point(self.thread.radius_nodes_mm()[0]),
                I::point(t.maximum_source_du_mm),
                I::point(t.maximum_source_dv_mm),
                I::point(c.maximum_source_du_mm),
                I::point(c.maximum_source_dv_mm),
                source_speed * I::point(source_top[0]) / minimum,
                I::point(2.) * I::point(source_top[0]),
                I::point(t.maximum_source_du_radius_mm)
                    + source_outer * I::point(t.maximum_source_du_angle_radians),
                source_outer + I::point(self.thread.radius_nodes_mm()[0]),
            ]);
            let mut radius = 0_f64;
            for i in 0..head.segment_count() {
                work.step(4096)?;
                let query = head.certified_interval_bounds(i, [0., 1.]);
                radius = radius.max(query.position_mm[0][1]);
                source_derivative = maximum([
                    source_derivative,
                    abs(interval(query.derivative_mm[0])) + abs(interval(query.derivative_mm[1])),
                ]);
            }
            source_derivative = maximum([
                source_derivative,
                tau * I::point(radius),
                I::point(radius),
                source_speed,
            ]);
            let square = self.enclosure.into_iter().fold(I::point(0.), |s, x| {
                s + I::point(x[0].abs().max(x[1].abs())).square()
            });
            let extent = I::new(0., square.hi).sqrt();
            let distance =
                world.translation + world.linear * extent + world.target_norm * local_distance;
            let derivative =
                world.linear * source_derivative + world.target_norm * local_derivative;
            let mass = abs(
                interval(self.world_volume_mm3()) * world.occurrence_determinant
                    - interval(target.world_volume_mm3()),
            );
            compatible(
                [distance, derivative, mass]
                    .into_iter()
                    .all(|x| x.valid() && x.hi.is_finite())
                    && distance.hi <= options.maximum_distance_mm
                    && derivative.hi <= options.maximum_derivative_mm
                    && mass.hi <= options.maximum_volume_discrepancy_mm3,
            )?;
            work.step(1)?;
            Ok(PeriodicMaterialMatch {
                source: self,
                target,
                distance: distance.hi,
                derivative: derivative.hi,
                volume_discrepancy: mass.hi,
                work: options.maximum_fixed_work - work.remaining,
                cavity_cells: c.cells_used,
            })
        })();
        result.map_err(|e| {
            if work.cancelled.get() {
                PhillipsReadError::Cancelled
            } else {
                e
            }
        })
    }
}

#[cfg(test)]
mod tests;
