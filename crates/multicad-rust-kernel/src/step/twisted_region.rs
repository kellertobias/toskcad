//! Closed twisted extrusions of planar loops of affine-radius spirals.
//!
//! A loop lies in the section plane at the lower height; every section at a
//! greater height is that loop turned about the axis by the twist rate times
//! the rise. Each loop piece `r(θ) = r₀ + (r₁ - r₀)(θ - θ₀)/(θ₁ - θ₀)` is the
//! certified degree-eleven spiral of the helical thread exchange, and each
//! side face is the tensor product of that piece with the certified axial
//! rotation of the helical gear exchange, at an affine height. Caps stay
//! exactly planar. A constant radius gives an arc of the axis circle.
use super::{
    helical::{Spline, product},
    herringbone::{Chart, axial_rotation},
    *,
};

const MAX_ERROR_MM: f64 = 1e-7;
const MAX_PIECES: usize = 256;

/// One loop piece from its first angle and radius to its second.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Spiral {
    pub angle_rad: [f64; 2],
    pub radius_mm: [f64; 2],
}

impl Spiral {
    fn reversed(self) -> Self {
        Self {
            angle_rad: [self.angle_rad[1], self.angle_rad[0]],
            radius_mm: [self.radius_mm[1], self.radius_mm[0]],
        }
    }

    fn point(self, end: usize) -> [f64; 2] {
        let (sin, cos) = self.angle_rad[end].portable_sin_cos();
        [self.radius_mm[end] * cos, self.radius_mm[end] * sin]
    }

    /// Exact `½∮ r² dθ` of this piece.
    fn signed_area(self) -> f64 {
        let span = self.angle_rad[1] - self.angle_rad[0];
        let [r0, r1] = self.radius_mm;
        span * (r0 * r0 + r0 * r1 + r1 * r1) / 6.
    }

    fn chart(self, ends: [[f64; 2]; 2], tolerance: f64) -> Result<Chart, StepError> {
        let [a, b] = self.angle_rad;
        // The spline turns by `sense * t`, so a decreasing angle runs on
        // the increasing parameter `-θ`.
        let spline = if b > a {
            Spline::new(a, b, self.radius_mm, 1., tolerance)?
        } else {
            Spline::new(-a, -b, self.radius_mm, -1., tolerance)?
        };
        let mut chart = Chart::from_spin(spline);
        let last = chart.poles.len() - 1;
        for (index, expected) in [0, last].into_iter().zip(ends) {
            let [x, y] = std::array::from_fn(|i| chart.poles[index][i] - expected[i]);
            if x.portable_hypot(y) > tolerance {
                return Err(StepError::Limit);
            }
            chart.poles[index] = expected;
        }
        Ok(chart)
    }
}

/// Whether every loop is a closed, simple-looking chain this encoder takes.
fn validate(loops: &[Vec<Spiral>], [z0, z1]: [f64; 2], twist: f64) -> Result<(), StepError> {
    let finite = |x: f64| x.is_finite();
    if loops.is_empty() || ![z0, z1, twist].into_iter().all(finite) || z1 <= z0 || twist == 0. {
        return Err(StepError::InvalidGeometry);
    }
    for pieces in loops {
        if pieces.len() < 2 || pieces.len() > MAX_PIECES {
            return Err(StepError::InvalidGeometry);
        }
        for piece in pieces {
            if !piece
                .angle_rad
                .into_iter()
                .chain(piece.radius_mm)
                .all(finite)
                || piece.radius_mm.into_iter().any(|r| r <= 0.)
                || piece.angle_rad[0] == piece.angle_rad[1]
            {
                return Err(StepError::InvalidGeometry);
            }
        }
    }
    Ok(())
}

pub fn write_twisted_spirals(
    loops: &[Vec<Spiral>],
    heights_mm: [f64; 2],
    twist_rad_per_mm: f64,
    origin: Point,
    u: Point,
    v: Point,
    tolerance_mm: f64,
) -> Result<String, StepError> {
    validate(loops, heights_mm, twist_rad_per_mm)?;
    write_one(|writer| {
        for pieces in loops {
            let placement = Placement::new(origin, u, v)?;
            writer.append(|file| {
                encode(
                    file,
                    pieces,
                    heights_mm,
                    twist_rad_per_mm,
                    &placement,
                    tolerance_mm,
                )
            })?;
        }
        Ok(())
    })
}

fn encode(
    file: &mut StepFile,
    pieces: &[Spiral],
    [z0, z1]: [f64; 2],
    twist: f64,
    placement: &Placement,
    tolerance: f64,
) -> Result<Entity, StepError> {
    if !tolerance.is_finite() || tolerance <= 0. {
        return Err(StepError::InvalidGeometry);
    }
    let tolerance = tolerance.min(MAX_ERROR_MM);
    let height = z1 - z0;
    let radius = pieces
        .iter()
        .flat_map(|piece| piece.radius_mm)
        .fold(0., f64::max);
    let angle = twist * height;
    let scale = placement
        .origin
        .into_iter()
        .map(f64::abs)
        .fold(0., f64::max)
        + radius * (1. + angle.abs())
        + z1.abs()
        + 1.;
    if 128. * f64::EPSILON * scale > tolerance / 4. {
        return Err(StepError::Limit);
    }
    // Material lies left of a counterclockwise outer loop.
    let area: f64 = pieces.iter().map(|piece| piece.signed_area()).sum();
    if !area.is_finite() || area.abs() <= tolerance * radius {
        return Err(StepError::InvalidGeometry);
    }
    let pieces: Vec<Spiral> = if area > 0. {
        pieces.to_vec()
    } else {
        pieces.iter().rev().map(|piece| piece.reversed()).collect()
    };
    let count = pieces.len();
    let ends: Vec<[f64; 2]> = pieces.iter().map(|piece| piece.point(0)).collect();
    for (index, piece) in pieces.iter().enumerate() {
        let [x, y] = piece.point(1);
        let [a, b] = ends[(index + 1) % count];
        if (x - a).portable_hypot(y - b) > tolerance / 16. {
            return Err(StepError::InvalidGeometry);
        }
    }
    let transverse = pieces
        .iter()
        .enumerate()
        .map(|(index, piece)| {
            piece.chart([ends[index], ends[(index + 1) % count]], tolerance / 16.)
        })
        .collect::<Result<Vec<_>, _>>()?;
    let (axial, rise) = axial_rotation(angle, height, tolerance / (16. * radius))?;
    let heights: Vec<f64> = rise.iter().map(|z| z0 + z).collect();
    let phases = [[1., 0.], *axial.poles.last().expect("rotation endpoints")];
    let levels = [z0, z1];

    let stations: Vec<Vec<Vertex>> = (0..2)
        .map(|level| {
            transverse
                .iter()
                .map(|chart| {
                    file.vertex(
                        placement.point_f64(product(chart.poles[0], phases[level]), levels[level]),
                    )
                })
                .collect()
        })
        .collect();
    let rims: Vec<Vec<Entity>> = (0..2)
        .map(|level| {
            transverse
                .iter()
                .enumerate()
                .map(|(index, chart)| {
                    let [a, b] = [stations[level][index], stations[level][(index + 1) % count]];
                    let points = chart.poles.iter().map(|&point| {
                        placement.point_f64(product(point, phases[level]), levels[level])
                    });
                    chart.edge(file, [a, b], points)
                })
                .collect()
        })
        .collect();
    let rails: Vec<Entity> = transverse
        .iter()
        .enumerate()
        .map(|(index, chart)| {
            let points = axial
                .poles
                .iter()
                .zip(&heights)
                .map(|(&spin, &z)| placement.point_f64(product(chart.poles[0], spin), z));
            axial.edge(file, [stations[0][index], stations[1][index]], points)
        })
        .collect();

    let mut faces = Vec::with_capacity(count + 2);
    for (index, chart) in transverse.iter().enumerate() {
        let next = (index + 1) % count;
        let rows = chart
            .poles
            .iter()
            .map(|&point| {
                let row = axial
                    .poles
                    .iter()
                    .zip(&heights)
                    .map(|(&spin, &z)| file.point(placement.point_f64(product(point, spin), z)))
                    .collect::<Vec<_>>();
                Entities(&row).to_string()
            })
            .collect::<Vec<_>>()
            .join(",");
        let [um, vm] = [chart.multiplicities(), axial.multiplicities()];
        let [uk, vk] = [chart.knot_values(), axial.knot_values()];
        let surface = file.entity(format_args!(
            "B_SPLINE_SURFACE_WITH_KNOTS('',{},{},({rows}),.UNSPECIFIED.,.F.,.F.,.F.,({um}),({vm}),({uk}),({vk}),.UNSPECIFIED.)",
            chart.degree, axial.degree,
        ));
        let wire = [
            (rims[0][index], true),
            (rails[next], true),
            (rims[1][index], false),
            (rails[index], false),
        ]
        .map(|(edge, sense)| file.oriented_edge(edge, sense));
        let bound = file.bound(&wire, true);
        faces.push(file.face_on(&[bound], surface, true));
    }
    for (level, top) in [(0, false), (1, true)] {
        let mut edges = rims[level]
            .iter()
            .map(|&edge| file.oriented_edge(edge, top))
            .collect::<Vec<_>>();
        if !top {
            edges.reverse();
        }
        let bound = file.bound(&edges, true);
        let point = file.point(placement.point_f64([0.; 2], levels[level]));
        let normal = placement.normal.map(|n| if top { n } else { -n });
        faces.push(file.face(&[bound], point, normal, placement.u));
    }
    Ok(file.solid(&faces))
}

#[cfg(test)]
mod tests;
