//! Derived literal STEP coordinate semantics; no requested tolerance is an error.
use super::super::decode::{ChartEnclosure, icross, idot, references};
use super::super::*;
use crate::{head_profile::TransportError, tapered_thread::interval::I};
#[cfg(test)]
mod tests;
#[derive(Clone, Copy)]
pub(super) struct Axis {
    pub origin: [I; 3],
    pub axes: [[I; 3]; 3],
}
fn number(value: f64, conversion: bool) -> I {
    if conversion {
        I::new(value.next_down(), value.next_up())
    } else {
        I::point(value)
    }
}
fn checked(value: I) -> Result<I> {
    value.checked().map_err(|_| ReadError::Unsupported)
}
fn unit(vector: [I; 3]) -> Result<[I; 3]> {
    let scale = vector.into_iter().map(I::magnitude).fold(0_f64, f64::max);
    if !scale.is_finite() || scale == 0. {
        return Err(ReadError::Invalid);
    }
    let scaled = vector.map(|v| I::new((v.lo / scale).next_down(), (v.hi / scale).next_up()));
    let square = checked(idot(scaled, scaled))?;
    if square.lo <= 0. {
        return Err(ReadError::Unsupported);
    }
    let length = square.sqrt();
    let out = scaled.map(|v| v / length);
    for v in out {
        checked(v)?;
    }
    Ok(out)
}
fn raw_direction(shell: &Shell<'_, '_>, id: u32, conversion: bool) -> Result<[I; 3]> {
    let [Value::String(_), values] = shell.args(id, "DIRECTION")? else {
        return Err(ReadError::Invalid);
    };
    Ok(numbers(values)?.map(|v| number(v, conversion)))
}
impl Axis {
    pub fn read(shell: &Shell<'_, '_>, id: u32, conversion: bool) -> Result<Self> {
        let [Value::String(_), origin, normal, reference_axis] =
            shell.args(id, "AXIS2_PLACEMENT_3D")?
        else {
            return Err(ReadError::Invalid);
        };
        let n = unit(raw_direction(shell, reference(normal)?, conversion)?)?;
        let v = unit(icross(
            n,
            unit(raw_direction(
                shell,
                reference(reference_axis)?,
                conversion,
            )?)?,
        ))?;
        Ok(Self {
            origin: shell
                .point(reference(origin)?)?
                .map(|v| number(v, conversion)),
            axes: [icross(v, n), v, n],
        })
    }
}
#[derive(Clone, Copy)]
pub(super) struct Canonical {
    chosen: Frame,
    literal: ChartEnclosure,
    numeric: ChartEnclosure,
    thread_axis: Axis,
    origin: [I; 3],
}
fn zero(values: [I; 3]) -> Result<()> {
    if values
        .into_iter()
        .all(|v| v.valid() && v.lo <= 0. && v.hi >= 0.)
    {
        Ok(())
    } else {
        Err(ReadError::Unsupported)
    }
}
fn distance(actual: Point, bounds: [I; 3]) -> f64 {
    actual.into_iter().zip(bounds).fold(0_f64, |sum, (a, v)| {
        (sum + (a - v.lo).abs().max((a - v.hi).abs())).next_up()
    })
}
impl Canonical {
    pub fn new(shell: &Shell<'_, '_>, faces: &[u32], frame: Frame, count: usize) -> Result<Self> {
        let first = *faces.first().ok_or(ReadError::Invalid)?;
        let cap = *faces.get(count).ok_or(ReadError::Invalid)?;
        let (surface, _, _) = face(shell, cap)?;
        let [Value::String(_), axis] = shell.args(surface, "PLANE")? else {
            return Err(ReadError::Invalid);
        };
        let axis = reference(axis)?;
        let (tensor, _, _) = face(shell, first)?;
        let args = shell.args(tensor, "B_SPLINE_SURFACE_WITH_KNOTS")?;
        let Value::List(rows) = args.get(3).ok_or(ReadError::Invalid)? else {
            return Err(ReadError::Invalid);
        };
        let ids = refs(rows.first().ok_or(ReadError::Invalid)?)?;
        let pole = shell.point(*ids.first().ok_or(ReadError::Invalid)?)?;
        let numeric_axis = Axis::read(shell, axis, true)?;
        let matrix = |conversion| -> Result<(ChartEnclosure, [I; 3])> {
            let a = Axis::read(shell, axis, conversion)?;
            let p = pole.map(|v| number(v, conversion));
            let along = idot(std::array::from_fn(|i| p[i] - a.origin[i]), a.axes[2]);
            let origin = std::array::from_fn(|i| a.origin[i] + a.axes[2][i] * along);
            Ok((ChartEnclosure::from_axes(origin, a.axes)?, origin))
        };
        let (literal, _) = matrix(false)?;
        let (numeric, origin) = matrix(true)?;
        Ok(Self {
            chosen: frame,
            literal,
            numeric,
            thread_axis: numeric_axis,
            origin,
        })
    }
    pub fn local_point(&self, point: Point, conversion: bool) -> Result<[I; 3]> {
        if conversion {
            self.numeric.local(point, true)
        } else {
            self.literal.local(point, false)
        }
    }
    pub fn project_point(&self, point: [I; 3], conversion: bool) -> Result<[I; 3]> {
        if conversion {
            self.numeric.point(point)
        } else {
            self.literal.point(point)
        }
    }
    pub fn local_vector(&self, vector: Point, conversion: bool) -> Result<[I; 3]> {
        self.project_vector(vector.map(|v| number(v, conversion)), conversion)
    }
    pub fn project_vector(&self, vector: [I; 3], conversion: bool) -> Result<[I; 3]> {
        if conversion {
            self.numeric.vector(vector)
        } else {
            self.literal.vector(vector)
        }
    }
    pub fn local_axis(&self, shell: &Shell<'_, '_>, id: u32, conversion: bool) -> Result<Axis> {
        let a = Axis::read(shell, id, conversion)?;
        Ok(Axis {
            origin: self.local_point(shell.placement(id)?.origin, conversion)?,
            axes: [
                self.project_vector(a.axes[0], conversion)?,
                self.project_vector(a.axes[1], conversion)?,
                self.project_vector(a.axes[2], conversion)?,
            ],
        })
    }
    pub fn parallel_axis(&self, shell: &Shell<'_, '_>, id: u32) -> Result<()> {
        zero(icross(
            Axis::read(shell, id, true)?.axes[2],
            self.thread_axis.axes[2],
        ))
    }
    pub fn coplanar_origins(&self, shell: &Shell<'_, '_>, circle: u32, plane: u32) -> Result<()> {
        self.parallel_axis(shell, circle)?;
        self.parallel_axis(shell, plane)?;
        let a = Axis::read(shell, circle, true)?;
        let b = Axis::read(shell, plane, true)?;
        zero([
            idot(
                std::array::from_fn(|i| a.origin[i] - b.origin[i]),
                b.axes[2],
            ),
            I::point(0.),
            I::point(0.),
        ])
    }
    pub fn same_axis_line(&self, shell: &Shell<'_, '_>, a: u32, b: u32) -> Result<()> {
        self.parallel_axis(shell, a)?;
        self.parallel_axis(shell, b)?;
        let a = Axis::read(shell, a, true)?;
        let b = Axis::read(shell, b, true)?;
        zero(icross(
            std::array::from_fn(|i| a.origin[i] - b.origin[i]),
            b.axes[2],
        ))
    }
    pub fn coaxial_with_thread(&self, shell: &Shell<'_, '_>, id: u32) -> Result<()> {
        self.parallel_axis(shell, id)?;
        let a = Axis::read(shell, id, true)?;
        zero(icross(
            std::array::from_fn(|i| a.origin[i] - self.origin[i]),
            self.thread_axis.axes[2],
        ))
    }
    pub fn transport(&self, shell: &Shell<'_, '_>, head: &[u32]) -> Result<TransportError> {
        let mut pending = head.to_vec();
        let mut visited = BTreeSet::new();
        let mut points = Vec::new();
        let mut axes = Vec::new();
        let mut extent = 1_f64;
        while let Some(id) = pending.pop() {
            if !visited.insert(id) {
                continue;
            }
            for entity in shell.components(id)? {
                if entity.name == "CARTESIAN_POINT"
                    && matches!(entity.args.as_slice(),[Value::String(_),Value::List(v)] if v.len()==3)
                {
                    let p = shell.point(id)?;
                    for conversion in [false, true] {
                        let bound = self
                            .local_point(p, conversion)?
                            .into_iter()
                            .fold(0_f64, |sum, v| (sum + v.magnitude()).next_up());
                        extent = extent.max(bound);
                    }
                    points.push(p);
                }
                if entity.name == "AXIS2_PLACEMENT_3D" {
                    axes.push(id);
                }
                for value in &entity.args {
                    references(value, &mut pending);
                }
            }
        }
        // Coordinates of a finite support are relative to its own origin.
        // Both that origin and the trim/control points were included above;
        // their difference is bounded by twice the absolute chart bound.
        extent = (2. * extent).next_up();
        // Add literal radii for circular/revolved interiors between trims.
        for &id in &visited {
            for entity in shell.components(id)? {
                if matches!(
                    entity.name,
                    "CIRCLE"
                        | "CYLINDRICAL_SURFACE"
                        | "CONICAL_SURFACE"
                        | "SPHERICAL_SURFACE"
                        | "TOROIDAL_SURFACE"
                ) {
                    for v in &entity.args {
                        if let Value::Number(r) = v {
                            extent = (extent + r.abs().next_up()).next_up();
                        }
                    }
                }
            }
        }
        let mut out = TransportError::ZERO;
        for p in points {
            let actual = self.chosen.local(p);
            out.arithmetic_mm = out
                .arithmetic_mm
                .max(distance(actual, self.local_point(p, false)?));
            out.numeric_conversion_mm = out
                .numeric_conversion_mm
                .max(distance(actual, self.local_point(p, true)?));
        }
        let mut axis_error = TransportError::ZERO;
        for id in axes {
            let actual = shell.placement(id)?;
            let vectors = [actual.u, actual.v, actual.normal];
            for conversion in [false, true] {
                let bounds = self.local_axis(shell, id, conversion)?;
                let mut error = 0_f64;
                for (p, bounds) in vectors.into_iter().zip(bounds.axes) {
                    let local = self
                        .chosen
                        .rows
                        .map(|r| dot(r, p) / self.chosen.determinant);
                    error = (error + (distance(local, bounds) * extent).next_up()).next_up();
                }
                if conversion {
                    axis_error.numeric_conversion_mm = axis_error.numeric_conversion_mm.max(error);
                } else {
                    axis_error.arithmetic_mm = axis_error.arithmetic_mm.max(error);
                }
            }
        }
        out.arithmetic_mm = (out.arithmetic_mm + axis_error.arithmetic_mm).next_up();
        out.numeric_conversion_mm =
            (out.numeric_conversion_mm + axis_error.numeric_conversion_mm).next_up();
        require(out.total_mm().is_finite())?;
        Ok(out)
    }
}
