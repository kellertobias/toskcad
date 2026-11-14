//! A single closed rectangular stock shell with a two-radius blind recess.
//! The bore center is the local origin and depth increases along positive Z.
use super::*;
use std::collections::BTreeMap;
use std::f64::consts::FRAC_PI_2;

#[derive(Clone, Copy, Debug)]
pub struct RectangularBlindCounterbore {
    pub min_mm: Point2,
    pub max_mm: Point2,
    pub height_mm: f64,
    pub bore_radius_mm: f64,
    pub bore_depth_mm: f64,
    pub seat_radius_mm: f64,
    pub seat_depth_mm: f64,
    pub minimum_wall_mm: f64,
    pub minimum_floor_mm: f64,
}
struct Builder {
    vertices: Vec<Point>,
    vertex_ids: BTreeMap<[u64; 3], usize>,
    edges: Vec<Edge>,
    edge_ids: BTreeMap<[usize; 2], usize>,
    faces: Vec<Face>,
}
impl Builder {
    fn new() -> Self {
        Self {
            vertices: Vec::new(),
            vertex_ids: BTreeMap::new(),
            edges: Vec::new(),
            edge_ids: BTreeMap::new(),
            faces: Vec::new(),
        }
    }
    fn vertex(&mut self, p: Point) -> usize {
        let key = p.map(|x| if x == 0. { 0 } else { x.to_bits() });
        *self.vertex_ids.entry(key).or_insert_with(|| {
            let id = self.vertices.len();
            self.vertices.push(p);
            id
        })
    }
    fn face(
        &mut self,
        points: &[Point],
        curves: &[Curve],
        support: Support,
        center: Point,
        area: f64,
        source: usize,
    ) {
        let ids: Vec<_> = points.iter().map(|&p| self.vertex(p)).collect();
        let mut boundary = Vec::with_capacity(ids.len());
        for (i, (&a, &b)) in ids
            .iter()
            .zip(ids.iter().cycle().skip(1))
            .enumerate()
            .take(ids.len())
        {
            let mut key = [a, b];
            key.sort();
            let edge = if let Some(&edge) = self.edge_ids.get(&key) {
                let old = &self.edges[edge];
                assert_eq!(old.vertices, [b, a]);
                let expected = match old.curve {
                    Curve::Line => Curve::Line,
                    Curve::Arc {
                        center_mm,
                        radius_mm,
                        sweep_rad,
                    } => Curve::Arc {
                        center_mm,
                        radius_mm,
                        sweep_rad: -sweep_rad,
                    },
                };
                // Both uses are authored by this sealed factory. A mismatch
                // is a programming error, not a recoverable source input.
                assert_eq!(expected, curves[i]);
                edge
            } else {
                let edge = self.edges.len();
                self.edges.push(Edge {
                    vertices: [a, b],
                    curve: curves[i],
                });
                self.edge_ids.insert(key, edge);
                edge
            };
            boundary.push(OrientedEdge {
                edge,
                forward: self.edges[edge].vertices == [a, b],
            });
        }
        self.faces.push(Face {
            support,
            boundary,
            center_mm: center,
            area_mm2: area,
            sources: vec![source],
        });
    }
}
fn ring(radius: f64, z: f64) -> [Point; 4] {
    [
        [radius, 0., z],
        [0., radius, z],
        [-radius, 0., z],
        [0., -radius, z],
    ]
}
fn arc(radius: f64, z: f64, forward: bool) -> Curve {
    Curve::Arc {
        center_mm: [0., 0., z],
        radius_mm: radius,
        sweep_rad: if forward { FRAC_PI_2 } else { -FRAC_PI_2 },
    }
}
fn cap_properties(points: &[Point], curves: &[Curve]) -> (Point, f64) {
    let mut twice_area = 0.;
    let mut mx = 0.;
    let mut my = 0.;
    for (i, (a, b)) in points
        .iter()
        .zip(points.iter().cycle().skip(1))
        .enumerate()
        .take(points.len())
    {
        match curves[i] {
            Curve::Line => {
                twice_area += a[0] * b[1] - b[0] * a[1];
                mx += (b[1] - a[1]) * (a[0] * a[0] + a[0] * b[0] + b[0] * b[0]) / 6.;
                my -= (b[0] - a[0]) * (a[1] * a[1] + a[1] * b[1] + b[1] * b[1]) / 6.;
            }
            Curve::Arc {
                radius_mm: r,
                sweep_rad,
                ..
            } => {
                twice_area += r * r * sweep_rad;
                let primitive = |x: f64| x - x * x * x / 3.;
                mx += r * r * r / 2. * (primitive(b[1] / r) - primitive(a[1] / r));
                my -= r * r * r / 2. * (primitive(b[0] / r) - primitive(a[0] / r));
            }
        }
    }
    let area = twice_area / 2.;
    ([mx / area, my / area, points[0][2]], area.abs())
}
// Error-free TwoSum keeps declared wall/floor reserves at the original
// binary source values, including an exactly representable limiting station.
fn sum_leq(a: f64, b: f64, limit: f64) -> bool {
    let sum = a + b;
    let bb = sum - a;
    let residual = (a - (sum - bb)) + (b - bb);
    sum < limit || (sum == limit && residual <= 0.)
}
impl ProfileStack {
    /// Exact line/circle boundary construction. Source indices are four stock
    /// sides (0..4), back(4), entry(5), seat wall/shoulder(6), bore/floor(7).
    /// This does not approximate a different stock section by its bounds.
    pub fn blind_counterbore(
        s: RectangularBlindCounterbore,
        tolerance_mm: f64,
        minimum_mm: f64,
    ) -> Result<Self, Error> {
        let values = s.min_mm.into_iter().chain(s.max_mm).chain([
            s.height_mm,
            s.bore_radius_mm,
            s.bore_depth_mm,
            s.seat_radius_mm,
            s.seat_depth_mm,
            s.minimum_wall_mm,
            s.minimum_floor_mm,
        ]);
        if !tolerance_mm.is_finite()
            || tolerance_mm < 0.
            || !minimum_mm.is_finite()
            || minimum_mm <= tolerance_mm
            || values.into_iter().any(|x| !valid_coordinate(x))
        {
            return Err(Error::Invalid);
        }
        Profile::Rectangle {
            min: s.min_mm,
            max: s.max_mm,
        }
        .check(minimum_mm)?;
        if [
            s.height_mm,
            s.bore_radius_mm,
            s.bore_depth_mm,
            s.seat_depth_mm,
            s.minimum_wall_mm,
            s.minimum_floor_mm,
        ]
        .into_iter()
        .any(|x| x < minimum_mm)
            || !sum_leq(s.bore_radius_mm, minimum_mm, s.seat_radius_mm)
            || !sum_leq(s.seat_depth_mm, minimum_mm, s.bore_depth_mm)
        {
            return Err(Error::Invalid);
        }
        if [-s.min_mm[0], -s.min_mm[1], s.max_mm[0], s.max_mm[1]]
            .into_iter()
            .any(|bound| !sum_leq(s.seat_radius_mm, s.minimum_wall_mm, bound))
            || !sum_leq(s.bore_depth_mm, s.minimum_floor_mm, s.height_mm)
        {
            return Err(Error::Unsupported);
        }
        let [xmin, ymin] = s.min_mm;
        let [xmax, ymax] = s.max_mm;
        let corners = [[xmax, ymin], [xmax, ymax], [xmin, ymax], [xmin, ymin]];
        let at = |p: Point2, z| [p[0], p[1], z];
        let mut b = Builder::new();
        let normals = [[1., 0., 0.], [0., 1., 0.], [-1., 0., 0.], [0., -1., 0.]];
        for i in 0..4 {
            let j = (i + 1) % 4;
            let p = [
                at(corners[i], 0.),
                at(corners[j], 0.),
                at(corners[j], s.height_mm),
                at(corners[i], s.height_mm),
            ];
            let width = (corners[j][0] - corners[i][0]).hypot(corners[j][1] - corners[i][1]);
            b.face(
                &p,
                &[Curve::Line; 4],
                Support::Plane { normal: normals[i] },
                [
                    (p[0][0] + p[1][0]) / 2.,
                    (p[0][1] + p[1][1]) / 2.,
                    s.height_mm / 2.,
                ],
                width * s.height_mm,
                i,
            );
        }
        b.face(
            &corners.map(|p| at(p, s.height_mm)),
            &[Curve::Line; 4],
            Support::Plane {
                normal: [0., 0., 1.],
            },
            [(xmin + xmax) / 2., (ymin + ymax) / 2., s.height_mm],
            (xmax - xmin) * (ymax - ymin),
            4,
        );
        let mouth = ring(s.seat_radius_mm, 0.);
        let seat = ring(s.seat_radius_mm, s.seat_depth_mm);
        let shoulder = ring(s.bore_radius_mm, s.seat_depth_mm);
        let floor = ring(s.bore_radius_mm, s.bore_depth_mm);
        for i in 0..4 {
            let j = (i + 1) % 4;
            let p = [at(corners[i], 0.), mouth[i], mouth[j], at(corners[j], 0.)];
            let curves = [
                Curve::Line,
                arc(s.seat_radius_mm, 0., true),
                Curve::Line,
                Curve::Line,
            ];
            let (center, area) = cap_properties(&p, &curves);
            b.face(
                &p,
                &curves,
                Support::Plane {
                    normal: [0., 0., -1.],
                },
                center,
                area,
                5,
            );
        }
        for (a, z, radius, source) in [
            (&mouth, &seat, s.seat_radius_mm, 6),
            (&shoulder, &floor, s.bore_radius_mm, 7),
        ] {
            for i in 0..4 {
                let j = (i + 1) % 4;
                let p = [z[i], z[j], a[j], a[i]];
                b.face(
                    &p,
                    &[
                        arc(radius, z[0][2], true),
                        Curve::Line,
                        arc(radius, a[0][2], false),
                        Curve::Line,
                    ],
                    Support::InwardCylinder {
                        center_mm: [0.; 2],
                        radius_mm: radius,
                    },
                    [
                        (a[i][0] + a[j][0]) * 2. / PI,
                        (a[i][1] + a[j][1]) * 2. / PI,
                        (a[0][2] + z[0][2]) / 2.,
                    ],
                    FRAC_PI_2 * radius * (z[0][2] - a[0][2]),
                    source,
                );
            }
        }
        for i in 0..4 {
            let j = (i + 1) % 4;
            let p = [seat[i], shoulder[i], shoulder[j], seat[j]];
            let curves = [
                Curve::Line,
                arc(s.bore_radius_mm, s.seat_depth_mm, true),
                Curve::Line,
                arc(s.seat_radius_mm, s.seat_depth_mm, false),
            ];
            let (center, area) = cap_properties(&p, &curves);
            b.face(
                &p,
                &curves,
                Support::Plane {
                    normal: [0., 0., -1.],
                },
                center,
                area,
                6,
            );
        }
        let p = [floor[0], floor[3], floor[2], floor[1]];
        b.face(
            &p,
            &[arc(s.bore_radius_mm, s.bore_depth_mm, false); 4],
            Support::Plane {
                normal: [0., 0., -1.],
            },
            [0., 0., s.bore_depth_mm],
            PI * s.bore_radius_mm * s.bore_radius_mm,
            7,
        );
        let removed = PI
            * (s.bore_radius_mm * s.bore_radius_mm * s.bore_depth_mm
                + (s.seat_radius_mm - s.bore_radius_mm)
                    * (s.seat_radius_mm + s.bore_radius_mm)
                    * s.seat_depth_mm);
        Ok(Self {
            tolerance_mm,
            bounds: [[xmin, ymin, 0.], [xmax, ymax, s.height_mm]],
            volume: (xmax - xmin) * (ymax - ymin) * s.height_mm - removed,
            vertices: b.vertices,
            edges: b.edges,
            faces: b.faces,
        })
    }
}

#[cfg(test)]
mod tests;
