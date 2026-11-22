//! Certified scalar volume without surface-area and inertia integration.
use super::*;
#[cfg(test)]
mod tests;
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VolumeProperties {
    pub volume_mm3: [f64; 2],
    pub characteristic_length_mm: f64,
}
impl RadialBoundary {
    /// Enclose the actual retained finite boundary volume with the same
    /// sewing/closure reserve and normalized width as the full mass query.
    /// No face-area or inertia accuracy is required for a scalar query.
    pub fn volume_with_interrupt(
        &self,
        relative_accuracy: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<VolumeProperties, Error> {
        self.volume_with_budget(relative_accuracy, interrupt, 512_000_000)
    }
    fn volume_with_budget(
        &self,
        accuracy: f64,
        interrupt: Option<&dyn Fn() -> bool>,
        limit: usize,
    ) -> Result<VolumeProperties, Error> {
        if !accuracy.is_finite() || !(0. ..1.).contains(&accuracy) || accuracy == 0. {
            return Err(Error::Invalid);
        }
        let mut work = Work {
            remaining: limit,
            interrupt,
        };
        work.step()?;
        let sources = faces::build(self, &mut work)?;
        let bounds = bounds::prove(self, &sources, &mut work)?;
        let target = I::point(accuracy) * measure::power_i(I::point(bounds.length), 3);
        (I::point(2.) * bounds.volume_error).bounded_by(target.lo)?;
        let center = self.ordered.caps.axial.chart.nominal_height * 0.5;
        let epsilon = accuracy / (256. * sources.len() as f64);
        let mut raw = I::point(0.);
        for source in sources {
            let values = integrate::face(
                source.sources,
                source.density,
                &self.ordered.caps.axial.chart,
                center,
                integrate::Accuracy {
                    scales: [bounds.volume_scale; QUANTITIES],
                    epsilon,
                    selection: measure::Selection::Volume,
                },
                &mut work,
            )?;
            raw = raw + values[0];
        }
        let volume = (raw + symmetric(bounds.volume_error)).positive()?;
        require_width(volume, target.lo)?;
        work.step()?;
        Ok(VolumeProperties {
            volume_mm3: pair(volume),
            characteristic_length_mm: bounds.length,
        })
    }
}
