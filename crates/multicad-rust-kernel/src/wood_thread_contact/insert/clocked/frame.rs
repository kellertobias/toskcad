//! Original literal frame expressions, retaining symbolic translations and
//! compositions. No rounded replacement Frame can enter a source projection.
use super::*;
pub(super) struct SourceFrame {
    origin: [Q; 3],
    axes: [[Q; 3]; 3],
    metric: Metric,
}
fn copied(m: &Metric) -> Metric {
    Metric {
        minimum_scale: m.minimum_scale,
        maximum_scale: m.maximum_scale,
        determinant: m.determinant,
    }
}
impl SourceFrame {
    pub fn copied(&self) -> Self {
        Self {
            origin: self.origin.clone(),
            axes: self.axes.clone(),
            metric: copied(&self.metric),
        }
    }
    pub fn minimum_scale(&self) -> I {
        self.metric.minimum_scale
    }
    pub fn maximum_scale(&self) -> I {
        self.metric.maximum_scale
    }

    pub fn literal(frame: Frame) -> R<Self> {
        let metric = scale(frame)?;
        Ok(Self {
            origin: vector(frame.origin_mm)?,
            axes: [
                vector(frame.axes[0])?,
                vector(frame.axes[1])?,
                vector(frame.axes[2])?,
            ],
            metric,
        })
    }
    pub fn translated(self, axis: usize, distance: f64) -> R<Self> {
        let d = q(distance)?;
        let mut origin = self.origin.clone();
        for (i, p) in origin.iter_mut().enumerate() {
            *p = bounded(&self.origin[i] + &self.axes[axis][i] * &d)?;
        }
        Ok(Self { origin, ..self })
    }
    pub fn placed(self, outer: Frame) -> R<Self> {
        let outer = Self::literal(outer)?;
        let mut axes = self.axes.clone();
        let mut origin = self.origin.clone();
        for i in 0..3 {
            origin[i] = bounded(
                &outer.origin[i]
                    + (0..3)
                        .map(|k| &outer.axes[k][i] * &self.origin[k])
                        .sum::<Q>(),
            )?;
            for (j, axis) in axes.iter_mut().enumerate() {
                axis[i] = bounded((0..3).map(|k| &outer.axes[k][i] * &self.axes[j][k]).sum())?;
            }
        }
        Ok(Self {
            origin,
            axes,
            metric: Metric {
                minimum_scale: self.metric.minimum_scale * outer.metric.minimum_scale,
                maximum_scale: self.metric.maximum_scale * outer.metric.maximum_scale,
                determinant: self.metric.determinant * outer.metric.determinant,
            },
        })
    }
    pub fn phase(self, phase: [[f64; 2]; 2]) -> R<Self> {
        let p = phase.map(|r| r.map(ip));
        let mut delta = ip(0.);
        for a in 0..2 {
            for b in 0..2 {
                delta = delta
                    + ip((p[0][a] * p[0][b] + p[1][a] * p[1][b]
                        - ip(if a == b { 1. } else { 0. }))
                    .magnitude());
            }
        }
        let determinant = p[0][0] * p[1][1] - p[0][1] * p[1][0];
        if delta.hi > 1e-6 || determinant.lo <= 0. {
            return Err(Error::Unsupported);
        }
        let mut axes = self.axes.clone();
        for (j, axis) in axes.iter_mut().enumerate().take(2) {
            for (i, value) in axis.iter_mut().enumerate() {
                *value = bounded(
                    &self.axes[0][i] * q(phase[0][j])? + &self.axes[1][i] * q(phase[1][j])?,
                )?;
            }
        }
        Ok(Self {
            origin: self.origin,
            axes,
            metric: Metric {
                minimum_scale: self.metric.minimum_scale * positive(ip(1.) - delta)?.sqrt(),
                maximum_scale: self.metric.maximum_scale * (ip(1.) + delta).sqrt(),
                determinant: self.metric.determinant * determinant,
            },
        })
    }
    pub fn axis_values(&self, index: usize) -> &[Q; 3] {
        &self.axes[index]
    }
    pub fn axis(&self, index: usize) -> R<[I; 3]> {
        Ok([
            exact(&self.axes[index][0])?,
            exact(&self.axes[index][1])?,
            exact(&self.axes[index][2])?,
        ])
    }
    pub fn normal_plane_distance(&self, length: f64) -> R<I> {
        let c = cross(&self.axes[0], &self.axes[1])?;
        let d = dot(&self.axes[2], &c)?;
        let d = positive(exact(&d)?)?;
        let c = [exact(&c[0])?, exact(&c[1])?, exact(&c[2])?];
        let norm = positive(c.into_iter().fold(ip(0.), |s, v| s + v.square()))?.sqrt();
        finite(ip(length) * d / norm)
    }
}
impl Projection {
    pub(super) fn expressions(
        from: &SourceFrame,
        to: &SourceFrame,
        interrupt: &dyn Fn() -> bool,
    ) -> R<Self> {
        let c = [
            cross(&to.axes[1], &to.axes[2])?,
            cross(&to.axes[2], &to.axes[0])?,
            cross(&to.axes[0], &to.axes[1])?,
        ];
        let d = dot(&to.axes[0], &c[0])?;
        if d <= Q::zero() {
            return Err(Error::Unsupported);
        }
        let mut offset = from.origin.clone();
        for (i, value) in offset.iter_mut().enumerate() {
            *value = bounded(&from.origin[i] - &to.origin[i])?;
        }
        let origin = [
            bounded(dot(&c[0], &offset)? / &d)?,
            bounded(dot(&c[1], &offset)? / &d)?,
            bounded(dot(&c[2], &offset)? / &d)?,
        ];
        let mut columns: [[Q; 3]; 3] = std::array::from_fn(|_| std::array::from_fn(|_| Q::zero()));
        for (j, column) in columns.iter_mut().enumerate() {
            checkpoint(interrupt)?;
            for i in 0..3 {
                column[i] = bounded(dot(&c[i], &from.axes[j])? / &d)?;
            }
        }
        let mut origin_i = [ip(0.); 3];
        let mut columns_i = [[ip(0.); 3]; 3];
        for i in 0..3 {
            origin_i[i] = exact(&origin[i])?;
            for j in 0..3 {
                columns_i[j][i] = exact(&columns[j][i])?;
            }
        }
        Ok(Self {
            origin,
            columns,
            origin_i,
            columns_i,
            metric: copied(&to.metric),
        })
    }
}
