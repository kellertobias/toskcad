//! Area, centroid and enclosed-volume integrals over a trimmed face.
//!
//! With `F` the integrand on the parameter plane and `G(u, v)` its integral
//! along u from a fixed base, Green's theorem turns the area integral into
//! the boundary integral of `G dv`. Each loop that winds about a periodic
//! surface contributes the full-period integral up to its start, and a face
//! that closes at a pole contributes it up to the pole.
use super::domain::{Domain, UvLoop};
use super::surface::PreparedInverse;
use super::{
    Boundary, Error, Face, GAUSS, Point, Result, Surface, cross, distance, dot, ensure, norm,
    scale, sub,
};

/// Volume flux, area and the three first moments of area.
pub(crate) type Moments = [f64; 5];

fn accumulate(into: &mut Moments, value: Moments, weight: f64) {
    for (slot, value) in into.iter_mut().zip(value) {
        *slot += value * weight;
    }
}

fn integrand(surface: &Surface, uv: [f64; 2], anchor: Point) -> Moments {
    let (s, [su, sv]) = surface.evaluate(uv);
    let normal = cross(su, sv);
    let area = norm(normal);
    [
        dot(sub(s, anchor), normal) / 3.,
        area,
        s[0] * area,
        s[1] * area,
        s[2] * area,
    ]
}

/// Gauss-Legendre integral of `f` over panels, each cut `refine` times.
fn quadrature(
    panels: &[f64],
    refine: usize,
    mut f: impl FnMut(f64) -> Result<Moments>,
) -> Result<Moments> {
    let mut total = [0.; 5];
    for pair in panels.windows(2) {
        let width = (pair[1] - pair[0]) / refine as f64;
        for piece in 0..refine {
            let middle = pair[0] + width * (piece as f64 + 0.5);
            for (node, weight) in GAUSS {
                accumulate(
                    &mut total,
                    f(middle + 0.5 * width * node)?,
                    0.5 * width * weight,
                );
            }
        }
    }
    Ok(total)
}

struct Integral<'a> {
    surface: &'a Surface,
    inverse: PreparedInverse<'a>,
    anchor: Point,
    base: [f64; 2],
    refine: usize,
}

impl Integral<'_> {
    /// `G(u, v)`: the integrand accumulated along u from the base.
    fn strip(&self, [u, v]: [f64; 2]) -> Moments {
        let (low, high, sign) = if u < self.base[0] {
            (u, self.base[0], -1.)
        } else {
            (self.base[0], u, 1.)
        };
        let panels = self.surface.breaks(0, low, high);
        quadrature(&panels, self.refine, |s| {
            Ok(integrand(self.surface, [s, v], self.anchor))
        })
        .unwrap_or_default()
        .map(|g| g * sign)
    }

    /// `H(v)`: the full-period strip accumulated along v from the base.
    fn band(&self, v: f64) -> Moments {
        let period = self.surface.periods()[0].unwrap_or(0.);
        let (low, high, sign) = if v < self.base[1] {
            (v, self.base[1], -1.)
        } else {
            (self.base[1], v, 1.)
        };
        let panels = self.surface.breaks(1, low, high);
        quadrature(&panels, self.refine, |w| {
            Ok(self.strip([self.base[0] + period, w]))
        })
        .unwrap_or_default()
        .map(|h| h * sign)
    }

    fn boundary(&self, boundary: &Boundary, trim: &UvLoop, tolerance: f64) -> Result<Moments> {
        let periods = self.surface.periods();
        let mut total = [0.; 5];
        for (i, segment) in trim.segments.iter().enumerate() {
            let Some(segment) = segment else { continue };
            let curve = &boundary.edges[segment.edge].curve;
            let ends = [trim.at(i, periods), trim.at(i + 1, periods)];
            let piece = quadrature(&segment.t, self.refine, |t| {
                let (point, tangent) = curve.evaluate(t);
                let along = (t - segment.t[0]) / (segment.t[1] - segment.t[0]);
                let seed =
                    [0, 1].map(|axis| ends[0][axis] + (ends[1][axis] - ends[0][axis]) * along);
                let uv = self.inverse.invert(point, Some(seed));
                let (s, [su, sv]) = self.surface.evaluate(uv);
                ensure(distance(s, point) <= tolerance, Error::Invalid)?;
                let (a, b, c) = (dot(su, su), dot(su, sv), dot(sv, sv));
                let (ru, rv) = (dot(tangent, su), dot(tangent, sv));
                let determinant = a * c - b * b;
                let rate = if determinant > 1e-12 * a * c {
                    (a * rv - b * ru) / determinant
                } else {
                    // Where the u derivative stalls, as at a spline end whose
                    // parametrization has zero speed, any motion along u
                    // moves the point by less than the tolerance the volume
                    // error already counts, and a tangent along the v
                    // derivative moves at the rate T.Sv/|Sv|^2. Parallel
                    // derivatives that do not stall leave it undefined.
                    let rate = rv / c;
                    let residual = norm(sub(tangent, scale(sv, rate)));
                    ensure(
                        a.sqrt() <= tolerance && c > 0. && residual <= 1e-9 * norm(tangent),
                        Error::Unsupported,
                    )?;
                    rate
                };
                Ok(self.strip(uv).map(|g| g * rate))
            })?;
            accumulate(&mut total, piece, 1.);
        }
        if trim.winding[0] != 0 {
            accumulate(
                &mut total,
                self.band(trim.points[0].uv[1]),
                -f64::from(trim.winding[0]),
            );
        }
        Ok(total)
    }
}

/// Integrals over the face region, taken counter-clockwise about the
/// surface normal: `[flux, area, moment x, moment y, moment z]`.
pub(crate) fn integrate(
    boundary: &Boundary,
    face: &Face,
    domain: &Domain,
    anchor: Point,
    refine: usize,
    tolerance: f64,
) -> Result<Moments> {
    let surface = &face.surface;
    let mut integral = Integral {
        surface,
        inverse: PreparedInverse::new(surface),
        anchor,
        base: [0.; 2],
        refine,
    };
    match domain {
        Domain::Complete([low, high]) => {
            integral.base = *low;
            Ok(integral.band(high[1]))
        }
        Domain::Trimmed { loops, pole } => {
            let start = loops[0].points[0].uv;
            let low = loops
                .iter()
                .flat_map(|trim| &trim.points)
                .fold(f64::INFINITY, |low, point| low.min(point.uv[0]));
            // Winding loops share the first loop's period window.
            integral.base = [
                if loops[0].winding[0] != 0 {
                    start[0]
                } else {
                    low
                },
                start[1],
            ];
            let mut total = [0.; 5];
            for trim in loops {
                accumulate(
                    &mut total,
                    integral.boundary(boundary, trim, tolerance)?,
                    1.,
                );
            }
            if let Some((v, _)) = pole {
                accumulate(
                    &mut total,
                    integral.band(*v),
                    f64::from(loops[0].winding[0]),
                );
            }
            Ok(total)
        }
    }
}
