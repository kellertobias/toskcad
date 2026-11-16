//! Actual retained-source face area and first moments. Mass/inertia integrands
//! are absent; closure collars and world/centroid width gates are unchanged.
use super::*;
#[cfg(test)]
mod tests;
impl RadialBoundary {
    /// Actual face areas and centroids at the same centered characteristic
    /// scale as properties_with_interrupt. Every face area width is bounded
    /// by accuracy*L² and each world centroid coordinate by accuracy*L. L is
    /// the larger of the collar-inflated source radius and axial span.
    /// Bounds are per face: summing returned area intervals retains the actual
    /// aggregate width, which is not promised to meet one face's target.
    /// This API does not substitute a nominal thread or inferred support.
    pub fn face_properties_with_interrupt(
        &self,
        relative_accuracy: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Vec<FaceProperties>, Error> {
        self.face_properties_with_budget(relative_accuracy, interrupt, 512_000_000)
    }
    fn face_properties_with_budget(
        &self,
        accuracy: f64,
        interrupt: Option<&dyn Fn() -> bool>,
        limit: usize,
    ) -> Result<Vec<FaceProperties>, Error> {
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
        let length = I::point(bounds.length);
        let area_target = I::point(accuracy) * length.square();
        for &error in &bounds.area_errors {
            (I::point(2.) * error).bounded_by(area_target.lo)?;
        }
        let center = self.ordered.caps.axial.chart.nominal_height * 0.5;
        let scales = std::array::from_fn(|i| match i {
            10 => length.square(),
            11..=13 => measure::power_i(length, 3),
            _ => bounds.volume_scale,
        });
        for scale in scales {
            scale.positive()?;
        }
        let mut output = Vec::with_capacity(sources.len());
        for (source, &error) in sources.into_iter().zip(&bounds.area_errors) {
            let measurement = |value: &measure::Values| {
                let area = (value[10] + symmetric(error)).positive()?;
                let local =
                    std::array::from_fn(|i| (value[11 + i] + symmetric(error * length)) / area);
                let centroid = world(self, local, center)?;
                require_width(area, area_target.lo)?;
                for value in centroid {
                    require_width(value, (I::point(accuracy) * length).lo)?;
                }
                Ok::<_, Error>((area, centroid))
            };
            let acceptance = |value: &measure::Values| measurement(value).is_ok();
            let value = integrate::surface_face(
                source.sources,
                source.density,
                &self.ordered.caps.axial.chart,
                center,
                scales,
                &acceptance,
                &mut work,
            )?;
            let (area, centroid) = measurement(&value)?;
            output.push(FaceProperties {
                face_index: source.index,
                area_mm2: pair(area),
                centroid_mm: centroid.map(pair),
                planar_support: source.plane,
            });
        }
        work.step()?;
        Ok(output)
    }
}
