//! Finite coaxial blind pockets retain the original cylindrical stock chart.
use crate::{
    cylinder::AxialCylinder,
    revolution::Revolution,
    tapered_thread::{Error, interval::I},
};
#[derive(Clone, Copy, Debug, PartialEq)]
/// Literal material-chart columns; never normalized by this provider.
pub struct Frame {
    pub origin_mm: [f64; 3],
    pub axes: [[f64; 3]; 3],
}
#[derive(Clone, Copy, Debug, PartialEq)]
/// Existing cylinder material placement, distinct from an authored reference frame.
pub enum Placement {
    /// The exact cylinder center/bottom are world coordinates.
    Axial,
    /// Maps local radial coordinates and travel from the first rim center.
    Framed(Frame),
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Seat {
    pub radius_mm: f64,
    pub depth_mm: f64,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Spec {
    pub radius_mm: f64,
    pub depth_mm: f64,
    pub seat: Option<Seat>,
    pub minimum_wall_mm: f64,
    pub minimum_floor_mm: f64,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Role {
    Back,
    Outer,
    Mouth,
    SeatWall,
    Shoulder,
    Bore,
    Floor,
    Axis,
}
#[derive(Clone, Debug)]
pub struct Pocket {
    stock: AxialCylinder,
    geometry: Revolution,
    frame: Frame,
    spec: Spec,
    roles: Vec<Role>,
    end: bool,
    source_error_mm: f64,
    placement_error_mm: f64,
}
fn dot(a: [I; 3], b: [I; 3]) -> I {
    a.into_iter()
        .zip(b)
        .fold(I::point(0.), |s, (a, b)| s + a * b)
}
fn cross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn length(a: [I; 3]) -> I {
    a.into_iter()
        .fold(I::point(0.), |s, x| s + x.square())
        .sqrt()
}
fn ratio(x: crate::Rational) -> I {
    let v = x.as_f64();
    I::new(v.next_down(), v.next_up())
}
fn metric(frame: Frame) -> Result<(f64, f64, [[I; 3]; 3]), Error> {
    if frame
        .origin_mm
        .into_iter()
        .chain(frame.axes.into_iter().flatten())
        .any(|x| !x.is_finite() || x.abs() > 1e7)
    {
        return Err(Error::Invalid);
    }
    let a = frame.axes.map(|v| v.map(I::point));
    let cof = [cross(a[1], a[2]), cross(a[2], a[0]), cross(a[0], a[1])];
    let det = dot(a[0], cof[0]).positive()?;
    let mut error = I::point(0.);
    for i in 0..3 {
        for j in 0..3 {
            error = error
                + I::point((dot(a[i], a[j]) - I::point(if i == j { 1. } else { 0. })).magnitude());
        }
    }
    if error.checked()?.hi >= 0.01 {
        return Err(Error::Unsupported);
    }
    let lo = (I::point(1.) - error).sqrt().positive()?.lo;
    let hi = (I::point(1.) + error).sqrt().checked()?.hi;
    Ok((lo, hi, cof.map(|row| row.map(|v| v / det))))
}
impl Pocket {
    pub fn new(
        stock: AxialCylinder,
        placement: Placement,
        entry: Frame,
        spec: Spec,
        tolerance_mm: f64,
        minimum_mm: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        let check = || {
            if interrupt.is_some_and(|f| f()) {
                Err(Error::Cancelled)
            } else {
                Ok(())
            }
        };
        check()?;
        let values = [
            spec.radius_mm,
            spec.depth_mm,
            spec.minimum_wall_mm,
            spec.minimum_floor_mm,
            tolerance_mm,
            minimum_mm,
        ];
        if values.into_iter().any(|v| !v.is_finite() || v <= 0.)
            || values[..4].iter().any(|&v| v < minimum_mm)
            || tolerance_mm >= minimum_mm
        {
            return Err(Error::Invalid);
        }
        let radius = ratio(stock.circle().radius());
        let height = ratio(stock.height());
        let r = stock.circle().radius().as_f64();
        let h = stock.height().as_f64();
        let (frame, origin_error) = match placement {
            Placement::Axial => {
                let [x, y] = stock.circle().center();
                let coords = [x, y, stock.bottom()];
                let origin = coords.map(|v| v.as_f64());
                (
                    Frame {
                        origin_mm: origin,
                        axes: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
                    },
                    coords
                        .into_iter()
                        .zip(origin)
                        .fold(I::point(0.), |sum, (v, f)| {
                            sum + I::point((ratio(v) - I::point(f)).magnitude())
                        }),
                )
            }
            Placement::Framed(frame) => (frame, I::point(0.)),
        };
        let (scale, upper, inverse) = metric(frame)?;
        metric(entry)?;
        let f = frame.axes;
        let literal = [
            f[0][1] * f[1][2] - f[0][2] * f[1][1],
            f[0][2] * f[1][0] - f[0][0] * f[1][2],
            f[0][0] * f[1][1] - f[0][1] * f[1][0],
        ];
        if literal != f[2] {
            return Err(Error::Unsupported);
        }
        let max_radius = spec.seat.map_or(spec.radius_mm, |s| s.radius_mm);
        if let Some(seat) = spec.seat
            && (!seat.radius_mm.is_finite()
                || !seat.depth_mm.is_finite()
                || seat.depth_mm < minimum_mm
                || (I::point(seat.radius_mm) - I::point(spec.radius_mm)).lo < minimum_mm
                || (I::point(spec.depth_mm) - I::point(seat.depth_mm)).lo < minimum_mm)
        {
            return Err(Error::Invalid);
        }
        let floor = h - spec.depth_mm;
        let seat_floor = spec.seat.map_or(h, |s| h - s.depth_mm);
        let subtraction =
            I::point((I::point(h) - I::point(spec.depth_mm) - I::point(floor)).magnitude())
                + I::point((I::point(h) - I::point(floor) - I::point(spec.depth_mm)).magnitude())
                + I::point(spec.seat.map_or(0., |s| {
                    (I::point(h) - I::point(s.depth_mm) - I::point(seat_floor)).magnitude()
                }))
                + I::point(spec.seat.map_or(0., |s| {
                    (I::point(h) - I::point(seat_floor) - I::point(s.depth_mm)).magnitude()
                }));
        let source_error = (origin_error
            + I::point(upper)
                * (I::point((radius - I::point(r)).magnitude())
                    + I::point((height - I::point(h)).magnitude())
                    + subtraction))
            .checked()?
            .hi;
        if source_error > tolerance_mm / 8. {
            return Err(Error::Unresolved);
        }
        if ((radius - I::point(max_radius)) * I::point(scale) - I::point(2. * source_error)).lo
            < spec.minimum_wall_mm
            || ((height - I::point(spec.depth_mm)) * I::point(scale) - I::point(2. * source_error))
                .lo
                < spec.minimum_floor_mm
        {
            return Err(Error::Unsupported);
        }
        let delta =
            std::array::from_fn(|i| I::point(entry.origin_mm[i]) - I::point(frame.origin_mm[i]));
        let local = inverse.map(|row| dot(row, delta));
        let axes = entry
            .axes
            .map(|column| inverse.map(|row| dot(row, column.map(I::point))));
        let end = axes[2][2].hi < 0.;
        let sign = if end { -1. } else { 1. };
        let cap = if end { h } else { 0. };
        let xy = |v: [I; 3]| [v[0], v[1], I::point(0.)];
        let gram = I::point((dot(xy(axes[0]), xy(axes[0])) - I::point(1.)).magnitude())
            + I::point((dot(xy(axes[1]), xy(axes[1])) - I::point(1.)).magnitude())
            + I::point(2. * dot(xy(axes[0]), xy(axes[1])).magnitude());
        if gram.checked()?.hi >= 0.01 {
            return Err(Error::Unsupported);
        }
        let tilt = I::point(axes[0][2].magnitude()) + I::point(axes[1][2].magnitude());
        let drift = I::point(axes[2][0].magnitude())
            + I::point(axes[2][1].magnitude())
            + I::point((axes[2][2] - I::point(sign)).magnitude());
        let offset = [local[0], local[1], local[2] - I::point(cap)];
        let placement_error = (I::point(upper)
            * (length(offset)
                + I::point(max_radius) * (gram + tilt)
                + I::point(spec.depth_mm) * drift)
            + I::point(source_error))
        .checked()?
        .hi;
        if placement_error > tolerance_mm / 4. {
            return Err(Error::Unsupported);
        }
        let mut points = vec![[0., 0.], [r, 0.], [r, h], [max_radius, h]];
        let mut roles = vec![Role::Back, Role::Outer, Role::Mouth];
        if let Some(seat) = spec.seat {
            points.extend([[seat.radius_mm, seat_floor], [spec.radius_mm, seat_floor]]);
            roles.extend([Role::SeatWall, Role::Shoulder]);
        }
        points.extend([[spec.radius_mm, floor], [0., floor]]);
        roles.extend([Role::Bore, Role::Floor, Role::Axis]);
        if !end {
            for p in &mut points {
                p[1] = h - p[1];
            }
            points.reverse();
            roles.reverse();
            roles.rotate_left(1);
        }
        for (i, role) in roles.iter().enumerate() {
            if *role != Role::Axis {
                let a = points[i];
                let b = points[(i + 1) % points.len()];
                let span = ((I::point(a[0]) - I::point(b[0])).square()
                    + (I::point(a[1]) - I::point(b[1])).square())
                .sqrt()
                    * I::point(scale);
                if span.checked()?.lo < minimum_mm {
                    return Err(Error::Unsupported);
                }
            }
        }
        let geometry = Revolution::new(vec![points], std::f64::consts::TAU).map_err(|cause| {
            if cause == crate::revolution::Error::Limit {
                Error::Limit
            } else {
                Error::Unresolved
            }
        })?;
        check()?;
        Ok(Self {
            stock,
            geometry,
            frame,
            spec,
            roles,
            end,
            source_error_mm: source_error,
            placement_error_mm: placement_error,
        })
    }
    pub fn stock(&self) -> &AxialCylinder {
        &self.stock
    }
    pub fn geometry(&self) -> &Revolution {
        &self.geometry
    }
    pub fn frame(&self) -> Frame {
        self.frame
    }
    pub fn spec(&self) -> Spec {
        self.spec
    }
    pub fn roles(&self) -> &[Role] {
        &self.roles
    }
    pub fn end(&self) -> bool {
        self.end
    }
    /// Euclidean numeric transfer bound for retained stock supports and authored
    /// cut stations; excludes the intentional removal of material.
    pub fn source_error_mm(&self) -> f64 {
        self.source_error_mm
    }
    /// Whole requested-tool to actual pocket bound, including the full seat radius.
    pub fn placement_error_mm(&self) -> f64 {
        self.placement_error_mm
    }
}
#[cfg(test)]
mod tests;
