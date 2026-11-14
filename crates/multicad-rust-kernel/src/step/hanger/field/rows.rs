//! Complete source-owned meridian control rows at true periodic seams.
//! Coefficients are not point samples; no sewn/literal body is admitted here.
use super::*;
use crate::{
    authored_hanger_bolt::{CurveKey, End, FaceKey},
    authored_internal_treatment::ThreadBandRole,
};
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BandPiece {
    Whole,
    CrestLowPhase,
    CrestHighPhase,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CellSeam {
    pub period: i32,
    pub upper_turn: bool,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MeridianRowKey {
    pub face: FaceKey,
    pub total_period: i32,
    pub piece: BandPiece,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HomogeneousControl {
    /// [w*x,w*y,w*z,w] of the original actual source curve.
    pub source_bounds: [[f64; 2]; 4],
    /// One deterministic approximation; literal token transfer remains later.
    pub literal: [f64; 4],
}
pub struct MeridianRow {
    key: MeridianRowKey,
    degree: usize,
    controls: [HomogeneousControl; 3],
    work: usize,
}
impl MeridianRow {
    pub fn key(&self) -> MeridianRowKey {
        self.key
    }
    pub fn degree(&self) -> usize {
        self.degree
    }
    pub fn controls(&self) -> &[HomogeneousControl] {
        &self.controls[..=self.degree]
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
}
fn point_control(point: SeamPoint) -> HomogeneousControl {
    HomogeneousControl {
        source_bounds: [
            point.source_bounds_mm[0],
            point.source_bounds_mm[1],
            point.source_bounds_mm[2],
            [1., 1.],
        ],
        literal: [point.point_mm[0], point.point_mm[1], point.point_mm[2], 1.],
    }
}
impl SourceField<'_> {
    pub fn meridian_row(
        &self,
        face: FaceKey,
        piece: BandPiece,
        seam: CellSeam,
        options: DomainOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<MeridianRow, Error> {
        let mut w = super::super::domain::work(options, self.bytes, interrupt)?;
        // Both outer row workspace and nested point-query temporaries are
        // included; source arena is borrowed/counts once. No heap allocation.
        let mut memory = Memory::new(self.bytes, options.maximum_bytes)?;
        memory.add(2 * SCRATCH)?;
        let (component, band) = match face {
            FaceKey::MetricThread { band, .. } => (Component::Metric, band),
            FaceKey::WoodThread { band, .. } => (Component::Wood, band),
            _ => return Err(Error::Invalid),
        };
        let c = self.covers[index(component)];
        // q>=0 and v>=0: j>=ceil(L/p) has no positive-area source cell.
        // Its possible boundary points remain available via seam_point.
        if seam.period < -(c.starts as i32) || seam.period >= c.last_cell {
            return Err(Error::Invalid);
        }
        if (band == ThreadBandRole::Crest) == (piece == BandPiece::Whole) {
            return Err(Error::Invalid);
        }
        let mut source_period = None;
        for curve in self.boundary.curves() {
            w.charge(1)?;
            if curve.faces.contains(&face)
                && let CurveKey::Cap {
                    component: owner,
                    end: End::Start,
                    period,
                    band: owner_band,
                } = curve.key
                && owner == component
                && owner_band == band
            {
                source_period = Some(period);
                break;
            }
        }
        let actual = source_period.ok_or(Error::Invalid)? as i32;
        let reference = seam.period
            - if piece == BandPiece::CrestLowPhase {
                1
            } else {
                0
            };
        if reference.rem_euclid(c.starts as i32) != actual {
            return Err(Error::Invalid);
        }
        let (first, last, carry) = match band {
            ThreadBandRole::Upper => (
                PhaseNode::Knot(Knot::CrestEnd),
                PhaseNode::Knot(Knot::UpperRoot),
                0,
            ),
            ThreadBandRole::Root => (
                PhaseNode::Knot(Knot::UpperRoot),
                PhaseNode::Knot(Knot::LowerRoot),
                0,
            ),
            ThreadBandRole::Lower => (
                PhaseNode::Knot(Knot::LowerRoot),
                PhaseNode::Knot(Knot::CrestStart),
                0,
            ),
            ThreadBandRole::Crest => {
                if piece == BandPiece::CrestLowPhase {
                    (PhaseNode::CrestSplit, PhaseNode::Knot(Knot::CrestEnd), 0)
                } else {
                    (PhaseNode::Knot(Knot::CrestStart), PhaseNode::CrestSplit, 1)
                }
            }
        };
        let mut point = |cell, node| -> Result<SeamPoint, Error> {
            if w.remaining() == 0 {
                return Err(Error::Limit);
            }
            let result = self.seam_point(
                component,
                cell,
                seam.upper_turn,
                node,
                DomainOptions {
                    maximum_work: w.remaining(),
                    ..options
                },
                interrupt,
            )?;
            w.charge(result.work_used)?;
            Ok(result)
        };
        let a = point(seam.period, first)?;
        let b = point(seam.period + carry, last)?;
        let empty = HomogeneousControl {
            source_bounds: [[0.; 2]; 4],
            literal: [0.; 4],
        };
        let mut controls = [empty; 3];
        controls[0] = point_control(a);
        let degree = if component == Component::Metric && band == ThreadBandRole::Root {
            w.charge(512)?;
            let s = self.boundary.source();
            let metric = s.metric();
            let h = metric.root_meridian().homogeneous_controls()[1].map(hull);
            let offset = p(metric.spec().pitch_mm) * p(f64::from(a.key.total_period))
                - p(metric.full_length_mm())
                - p(s.spec().plain_length_mm);
            let value = [h[1], p(0.), h[0] + h[2] * offset, h[2]];
            controls[1] = HomogeneousControl {
                source_bounds: value.map(|v| [v.lo, v.hi]),
                literal: value.map(I::midpoint),
            };
            2
        } else {
            1
        };
        controls[degree] = point_control(b);
        w.charge(1)?;
        Ok(MeridianRow {
            key: MeridianRowKey {
                face,
                total_period: a.key.total_period,
                piece,
            },
            degree,
            controls,
            work: w.used,
        })
    }
}
#[cfg(test)]
mod tests;
