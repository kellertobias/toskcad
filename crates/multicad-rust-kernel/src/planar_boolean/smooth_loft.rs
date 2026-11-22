//! Smooth lofts through three or more corresponding convex sections on
//! parallel planes. Each section corner follows a C¹ cubic rail whose height
//! is its parameter, so every intermediate horizontal section is the polygon
//! of the rail points at one parameter. Each wall is the ruled surface
//! between two neighboring rails: an exact degree (1, 3) B-spline surface.
use super::*;
use multicad_precision::PortableMath as _;

const MAX_EDGES: usize = 64;
const MAX_STATIONS: usize = 16;
/// Subdivision depth of one Bernstein convexity certificate.
const MAX_DEPTH: u32 = 24;

// Eight point Gauss-Legendre quadrature integrates degree 15 exactly.
const NODES: [f64; 4] = [
    0.183_434_642_495_649_8,
    0.525_532_409_916_329,
    0.796_666_477_413_626_7,
    0.960_289_856_497_536_3,
];
const WEIGHTS: [f64; 4] = [
    0.362_683_783_378_362,
    0.313_706_645_877_887_3,
    0.222_381_034_453_374_5,
    0.101_228_536_290_376_3,
];

fn gauss() -> impl Iterator<Item = (f64, f64)> {
    // Nodes and weights on the unit interval.
    NODES
        .into_iter()
        .zip(WEIGHTS)
        .flat_map(|(x, w)| [((1. - x) * 0.5, w * 0.5), ((1. + x) * 0.5, w * 0.5)])
}

/// A closed solid with two planar caps and ruled walls between C¹ cubic
/// corner rails. Strict convexity of every intermediate section excludes
/// crossing walls and zero-area sections.
#[derive(Clone, Debug, PartialEq)]
pub struct SmoothLoft {
    /// Station heights above the first cap, starting at zero and increasing.
    heights: Vec<f64>,
    /// Per section corner, the in-plane rail poles of the clamped cubic
    /// B-spline with the station heights as knots: interior knots have
    /// multiplicity two, so each rail has two poles per station.
    rails: Vec<Vec<[f64; 2]>>,
    clockwise: bool,
}

impl SmoothLoft {
    /// Interpolate corresponding sections at increasing heights with C¹
    /// rails. Each station tangent is the derivative of the parabola through
    /// it and its neighbors, so three stations give one exact parabola per
    /// corner.
    pub fn interpolate(
        sections: &[Vec<[f64; 2]>],
        heights: &[f64],
        tolerance: f64,
    ) -> Result<Self, Error> {
        validate_tolerance(tolerance)?;
        let stations = sections.len();
        if !(3..=MAX_STATIONS).contains(&stations) || heights.len() != stations {
            return Err(Error::Invalid);
        }
        let count = sections[0].len();
        if count > MAX_EDGES {
            return Err(Error::Limit);
        }
        if count < 3 || sections.iter().any(|ring| ring.len() != count) {
            return Err(Error::Invalid);
        }
        validate_heights(heights, tolerance)?;
        let steps: Vec<_> = heights.windows(2).map(|pair| pair[1] - pair[0]).collect();
        let rails = (0..count)
            .map(|corner| {
                let points: Vec<_> = sections.iter().map(|ring| ring[corner]).collect();
                let slopes: Vec<[f64; 2]> = points
                    .windows(2)
                    .zip(&steps)
                    .map(|(pair, &step)| std::array::from_fn(|i| (pair[1][i] - pair[0][i]) / step))
                    .collect();
                let parabola = |left: usize| -> [[f64; 2]; 3] {
                    // Derivatives at the three stations of the parabola
                    // through stations left, left + 1 and left + 2.
                    let [a, b] = [steps[left], steps[left + 1]];
                    let [s, t] = [slopes[left], slopes[left + 1]];
                    let middle: [f64; 2] = std::array::from_fn(|i| (b * s[i] + a * t[i]) / (a + b));
                    [
                        std::array::from_fn(|i| 2. * s[i] - middle[i]),
                        middle,
                        std::array::from_fn(|i| 2. * t[i] - middle[i]),
                    ]
                };
                let tangents: Vec<_> = (0..stations)
                    .map(|k| {
                        if k == 0 {
                            parabola(0)[0]
                        } else if k + 1 == stations {
                            parabola(stations - 3)[2]
                        } else {
                            parabola(k - 1)[1]
                        }
                    })
                    .collect();
                let mut poles = Vec::with_capacity(2 * stations);
                poles.push(points[0]);
                for (k, &step) in steps.iter().enumerate() {
                    let [p, q] = [points[k], points[k + 1]];
                    let [m, n] = [tangents[k], tangents[k + 1]];
                    poles.push(std::array::from_fn(|i| p[i] + m[i] * step / 3.));
                    poles.push(std::array::from_fn(|i| q[i] - n[i] * step / 3.));
                }
                poles.push(points[stations - 1]);
                poles
            })
            .collect();
        Self::from_rails(rails, heights.to_vec(), tolerance)
    }

    /// Admit rails given by their poles, two per station height.
    pub fn from_rails(
        rails: Vec<Vec<[f64; 2]>>,
        heights: Vec<f64>,
        tolerance: f64,
    ) -> Result<Self, Error> {
        validate_tolerance(tolerance)?;
        let stations = heights.len();
        if rails.len() > MAX_EDGES {
            return Err(Error::Limit);
        }
        if !(2..=MAX_STATIONS).contains(&stations)
            || rails.len() < 3
            || rails.iter().any(|rail| rail.len() != 2 * stations)
            || rails
                .iter()
                .flatten()
                .flatten()
                .any(|v| !v.is_finite() || v.abs() > MAX_COORDINATE)
        {
            return Err(Error::Invalid);
        }
        validate_heights(&heights, tolerance)?;
        let first: Vec<_> = rails.iter().map(|rail| rail[0]).collect();
        let area = (1..first.len() - 1)
            .map(|i| cross2(sub2(first[i], first[0]), sub2(first[i + 1], first[0])))
            .sum::<f64>();
        let loft = Self {
            heights,
            rails,
            clockwise: area < 0.,
        };
        loft.prove_sections(tolerance)?;
        Ok(loft)
    }

    pub fn heights(&self) -> &[f64] {
        &self.heights
    }

    pub fn height(&self) -> f64 {
        self.heights[self.heights.len() - 1]
    }

    /// In-plane rail poles per corner, two per station.
    pub fn rails(&self) -> &[Vec<[f64; 2]>] {
        &self.rails
    }

    /// The local height of every rail pole: the Greville abscissae of the
    /// station knots, which make each rail's height its parameter.
    pub fn pole_heights(&self) -> Vec<f64> {
        let h = &self.heights;
        let mut z = Vec::with_capacity(2 * h.len());
        z.push(h[0]);
        for pair in h.windows(2) {
            let step = pair[1] - pair[0];
            z.push(pair[0] + step / 3.);
            z.push(pair[1] - step / 3.);
        }
        z.push(h[h.len() - 1]);
        z
    }

    pub fn clockwise(&self) -> bool {
        self.clockwise
    }

    pub fn count(&self) -> usize {
        self.rails.len()
    }

    /// The two cap rings, at heights zero and `height()`.
    pub fn rings(&self) -> [Vec<Point>; 2] {
        let top = self.height();
        [
            self.rails
                .iter()
                .map(|rail| [rail[0][0], rail[0][1], 0.])
                .collect(),
            self.rails
                .iter()
                .map(|rail| {
                    let [x, y] = rail[rail.len() - 1];
                    [x, y, top]
                })
                .collect(),
        ]
    }

    /// The in-plane Bézier controls of one corner rail over one segment.
    fn segment(&self, corner: usize, segment: usize) -> [[f64; 2]; 4] {
        let rail = &self.rails[corner];
        let h = &self.heights;
        let start = if segment == 0 {
            rail[0]
        } else {
            let [before, after] = [h[segment] - h[segment - 1], h[segment + 1] - h[segment]];
            let [b, a] = [rail[2 * segment], rail[2 * segment + 1]];
            std::array::from_fn(|i| (after * b[i] + before * a[i]) / (before + after))
        };
        let last = h.len() - 1;
        let end = if segment + 1 == last {
            rail[rail.len() - 1]
        } else {
            let next = segment + 1;
            let [before, after] = [h[next] - h[next - 1], h[next + 1] - h[next]];
            let [b, a] = [rail[2 * next], rail[2 * next + 1]];
            std::array::from_fn(|i| (after * b[i] + before * a[i]) / (before + after))
        };
        [start, rail[2 * segment + 1], rail[2 * segment + 2], end]
    }

    fn segments(&self) -> usize {
        self.heights.len() - 1
    }

    fn step(&self, segment: usize) -> f64 {
        self.heights[segment + 1] - self.heights[segment]
    }

    /// A rail point at local parameter `t` of one segment.
    fn rail_point(&self, corner: usize, segment: usize, t: f64) -> [f64; 2] {
        bezier(self.segment(corner, segment), t)
    }

    /// The point of wall `edge` at parameter `u` along the section edge.
    pub fn point(&self, edge: usize, u: f64, segment: usize, t: f64) -> Point {
        let next = (edge + 1) % self.count();
        let [a, b] = [edge, next].map(|corner| self.rail_point(corner, segment, t));
        [
            (1. - u) * a[0] + u * b[0],
            (1. - u) * a[1] + u * b[1],
            self.heights[segment] + t * self.step(segment),
        ]
    }

    fn prove_sections(&self, tolerance: f64) -> Result<(), Error> {
        let count = self.count();
        let sign = if self.clockwise { -1. } else { 1. };
        let coordinate_error = 64.
            * f64::EPSILON
            * self
                .rails
                .iter()
                .flatten()
                .flat_map(|p| p.iter().map(|x| x.abs()))
                .fold(1., f64::max);
        let mut budget = MAX_WORK;
        for segment in 0..self.segments() {
            let controls: Vec<_> = (0..count).map(|c| self.segment(c, segment)).collect();
            for i in 0..count {
                let j = (i + 1) % count;
                let edge: [[f64; 2]; 4] =
                    std::array::from_fn(|n| sub2(controls[j][n], controls[i][n]));
                for k in (0..count).filter(|&k| k != i && k != j) {
                    let offset: [[f64; 2]; 4] =
                        std::array::from_fn(|n| sub2(controls[k][n], controls[i][n]));
                    // The signed edge-to-vertex determinant is a degree six
                    // Bernstein polynomial in the segment parameter.
                    let mut coefficients = [0.; 7];
                    for (a, edge) in edge.iter().enumerate() {
                        for (b, offset) in offset.iter().enumerate() {
                            coefficients[a + b] +=
                                sign * BINOMIAL3[a] * BINOMIAL3[b] * cross2(*edge, *offset)
                                    / BINOMIAL6[a + b];
                        }
                    }
                    let length = edge
                        .iter()
                        .map(|v| v[0].portable_hypot(v[1]))
                        .fold(0., f64::max);
                    let reach = offset
                        .iter()
                        .map(|v| v[0].portable_hypot(v[1]))
                        .fold(0., f64::max);
                    let guard = tolerance * length
                        + 16. * coordinate_error * (length + reach + coordinate_error);
                    positive(coefficients, guard, &mut budget)?;
                }
            }
        }
        Ok(())
    }

    pub fn bounds(&self) -> [Point; 2] {
        let z = self.pole_heights();
        bounds(
            self.rails
                .iter()
                .flat_map(|rail| rail.iter().zip(&z).map(|(&[x, y], &z)| [x, y, z])),
        )
    }

    pub fn cap_areas_mm2(&self) -> [f64; 2] {
        self.rings().map(|ring| section_area(&ring))
    }

    pub fn volume_mm3(&self) -> f64 {
        let count = self.count();
        let mut volume = 0.;
        for segment in 0..self.segments() {
            // A section area is degree six in the segment parameter, which
            // eight Gauss points integrate exactly.
            let area: f64 = gauss()
                .map(|(t, w)| {
                    let ring: Vec<_> = (0..count)
                        .map(|corner| {
                            let [x, y] = self.rail_point(corner, segment, t);
                            [x, y, 0.]
                        })
                        .collect();
                    w * section_area(&ring)
                })
                .sum();
            volume += area * self.step(segment);
        }
        volume
    }

    /// Area of one wall by composite Gauss quadrature, refined until two
    /// levels agree to a relative 1e-12.
    pub fn wall_area_mm2(&self, edge: usize) -> Result<f64, Error> {
        let next = (edge + 1) % self.count();
        let mut total = 0.;
        for segment in 0..self.segments() {
            let [a, b] = [edge, next].map(|corner| self.segment(corner, segment));
            let step = self.step(segment);
            let integrate = |parts: usize| -> f64 {
                let mut sum = 0.;
                for part in 0..parts {
                    for (s, w) in gauss() {
                        let t = (part as f64 + s) / parts as f64;
                        let [p, q] = [bezier(a, t), bezier(b, t)];
                        let [dp, dq] = [derivative(a, t), derivative(b, t)];
                        let along = [q[0] - p[0], q[1] - p[1], 0.];
                        for (u, v) in gauss() {
                            let rail = [
                                (1. - u) * dp[0] + u * dq[0],
                                (1. - u) * dp[1] + u * dq[1],
                                step,
                            ];
                            sum += w * v * norm(cross(along, rail));
                        }
                    }
                }
                sum / parts as f64
            };
            let mut parts = 1;
            let mut previous = integrate(parts);
            loop {
                parts *= 2;
                if parts > 1024 {
                    return Err(Error::Limit);
                }
                let current = integrate(parts);
                if (current - previous).abs() <= 1e-12 * current.abs() {
                    total += current;
                    break;
                }
                previous = current;
            }
        }
        Ok(total)
    }

    pub fn mesh(&self, options: Tessellation) -> Result<DisplayMesh, MeshError> {
        options.check()?;
        let count = self.count();
        let mut rows = Vec::with_capacity(self.segments());
        for segment in 0..self.segments() {
            let controls: Vec<_> = (0..count).map(|c| self.segment(c, segment)).collect();
            // Second parameter derivatives of a rail and the change of each
            // section edge vector bound the chordal and twist errors.
            let curvature = controls
                .iter()
                .map(|c| {
                    (0..2)
                        .map(|n| {
                            let [x, y] =
                                std::array::from_fn(|i| c[n + 2][i] - 2. * c[n + 1][i] + c[n][i]);
                            6. * x.portable_hypot(y)
                        })
                        .fold(0., f64::max)
                })
                .fold(0., f64::max);
            let twist = (0..count)
                .map(|i| {
                    let j = (i + 1) % count;
                    (0..3)
                        .map(|n| {
                            let [x, y] = std::array::from_fn(|axis| {
                                (controls[j][n + 1][axis] - controls[i][n + 1][axis])
                                    - (controls[j][n][axis] - controls[i][n][axis])
                            });
                            3. * x.portable_hypot(y)
                        })
                        .fold(0., f64::max)
                })
                .fold(0., f64::max);
            let step = self.step(segment);
            let chord = (curvature / (4. * options.deflection)).sqrt();
            let skew = twist / (2. * options.deflection);
            let rotation = (curvature + twist) / (step * (options.angle * 0.5).portable_sin());
            let divisions = chord.max(skew).max(rotation).max(1.).ceil();
            if !divisions.is_finite()
                || 2. * count as f64 * divisions > options.triangle_limit as f64
            {
                return Err(MeshError::BudgetExceeded);
            }
            rows.push(divisions as usize);
        }
        let total: usize = rows.iter().sum();
        // Bound before allocating.
        if 2 * count * (total + 1) > options.triangle_limit {
            return Err(MeshError::BudgetExceeded);
        }
        let mut mesh = DisplayMesh {
            vertices: Vec::with_capacity(3 * (count * (total + 1) + 2)),
            indices: Vec::with_capacity(6 * count * (total + 1)),
            face_indices: Vec::with_capacity(2 * count * (total + 1)),
        };
        for (segment, &divisions) in rows.iter().enumerate() {
            let first = usize::from(segment > 0);
            for row in first..=divisions {
                let t = row as f64 / divisions as f64;
                for corner in 0..count {
                    mesh.vertices.extend(self.point(corner, 0., segment, t));
                }
            }
        }
        let triangle = |mesh: &mut DisplayMesh, mut indices: [usize; 3], face| {
            if self.clockwise {
                indices.swap(1, 2);
            }
            mesh.triangle(indices, face);
        };
        for row in 0..total {
            for i in 0..count {
                let j = (i + 1) % count;
                let [a, b] = [row * count + i, row * count + j];
                let [d, c] = [a + count, b + count];
                triangle(&mut mesh, [a, b, c], (2 + i) as u32);
                triangle(&mut mesh, [a, c, d], (2 + i) as u32);
            }
        }
        for (cap, ring) in self.rings().iter().enumerate() {
            let center = mesh.vertices.len() / 3;
            let average = scale(
                ring.iter()
                    .fold([0.; 3], |sum, &p| add(sum, sub(p, ring[0]))),
                1. / count as f64,
            );
            mesh.vertices.extend(add(ring[0], average));
            let offset = cap * total * count;
            for i in 0..count {
                let [a, b] = [offset + i, offset + (i + 1) % count];
                triangle(
                    &mut mesh,
                    if cap == 0 {
                        [center, b, a]
                    } else {
                        [center, a, b]
                    },
                    cap as u32,
                );
            }
        }
        Ok(mesh)
    }
}

const BINOMIAL3: [f64; 4] = [1., 3., 3., 1.];
const BINOMIAL6: [f64; 7] = [1., 6., 15., 20., 15., 6., 1.];

fn validate_heights(heights: &[f64], tolerance: f64) -> Result<(), Error> {
    if heights[0] != 0.
        || heights.iter().any(|height| !height.is_finite())
        || heights
            .windows(2)
            .any(|pair| pair[1] - pair[0] <= tolerance || !pair[1].is_finite())
        || heights[heights.len() - 1] > MAX_COORDINATE
    {
        return Err(Error::Invalid);
    }
    Ok(())
}

/// Prove a Bernstein polynomial exceeds `guard` on the whole unit interval.
/// Positive coefficients bound it below; a sampled endpoint at or below the
/// guard refutes strict convexity.
fn positive(coefficients: [f64; 7], guard: f64, budget: &mut usize) -> Result<(), Error> {
    let mut stack = vec![(coefficients, 0_u32)];
    while let Some((c, depth)) = stack.pop() {
        charge(budget, 1)?;
        if c[0] <= guard || c[6] <= guard {
            return Err(Error::Unresolved);
        }
        if c.iter().all(|&value| value > guard) {
            continue;
        }
        if depth == MAX_DEPTH {
            return Err(Error::Unresolved);
        }
        let mut left = [0.; 7];
        let mut right = [0.; 7];
        let mut work = c;
        for level in 0..7 {
            left[level] = work[0];
            right[6 - level] = work[6 - level];
            for n in 0..6 - level {
                work[n] = 0.5 * (work[n] + work[n + 1]);
            }
        }
        stack.push((left, depth + 1));
        stack.push((right, depth + 1));
    }
    Ok(())
}

fn bezier(c: [[f64; 2]; 4], t: f64) -> [f64; 2] {
    let s = 1. - t;
    let [a, b, d, e] = [s * s * s, 3. * s * s * t, 3. * s * t * t, t * t * t];
    std::array::from_fn(|i| a * c[0][i] + b * c[1][i] + d * c[2][i] + e * c[3][i])
}

fn derivative(c: [[f64; 2]; 4], t: f64) -> [f64; 2] {
    let s = 1. - t;
    let [a, b, d] = [3. * s * s, 6. * s * t, 3. * t * t];
    std::array::from_fn(|i| {
        a * (c[1][i] - c[0][i]) + b * (c[2][i] - c[1][i]) + d * (c[3][i] - c[2][i])
    })
}

fn sub2(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}

fn cross2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

fn section_area(ring: &[Point]) -> f64 {
    let origin = ring[0];
    ring.windows(2)
        .map(|pair| cross(sub(pair[0], origin), sub(pair[1], origin))[2])
        .sum::<f64>()
        .abs()
        * 0.5
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rect(w: f64, h: f64, [x, y]: [f64; 2]) -> Vec<[f64; 2]> {
        vec![[x, y], [x + w, y], [x + w, y + h], [x, y + h]]
    }

    #[test]
    fn three_stations_follow_one_parabola_per_corner_with_exact_volume() {
        let sections = [
            rect(10., 6., [0., 0.]),
            rect(14., 8., [-2., -1.]),
            rect(10., 6., [0., 0.]),
        ];
        let loft = SmoothLoft::interpolate(&sections, &[0., 15., 30.], 1e-7).unwrap();
        // Width 10 + 8s(1-s) and depth 6 + 8s(1-s)/2 over s = z/30.
        let width = |s: f64| 10. + 16. * s * (1. - s);
        let depth = |s: f64| 6. + 8. * s * (1. - s);
        for (segment, t) in [(0, 0.), (0, 0.4), (0, 1.), (1, 0.3), (1, 1.)] {
            let s = (15. * segment as f64 + 15. * t) / 30.;
            let a = loft.point(0, 0., segment, t);
            let b = loft.point(1, 0., segment, t);
            let c = loft.point(2, 0., segment, t);
            assert!((b[0] - a[0] - width(s)).abs() < 1e-12);
            assert!((c[1] - b[1] - depth(s)).abs() < 1e-12);
            assert!((a[2] - 30. * s).abs() < 1e-12);
        }
        // Independent Simpson integration of the degree four section area.
        let n = 2000;
        let simpson: f64 = (0..=n)
            .map(|i| {
                let s = i as f64 / n as f64;
                let weight = if i == 0 || i == n {
                    1.
                } else if i % 2 == 1 {
                    4.
                } else {
                    2.
                };
                weight * width(s) * depth(s)
            })
            .sum::<f64>()
            * 30.
            / (3. * n as f64);
        assert!((loft.volume_mm3() - simpson).abs() < 1e-9);
        assert_eq!(loft.cap_areas_mm2(), [60., 60.]);
        let bounds = loft.bounds();
        assert!(bounds[0][0] <= -2. && bounds[1][2] == 30.);
        // The x = const side walls are planar trapezoid-like curved strips;
        // their area exceeds the flat projection.
        let front = loft.wall_area_mm2(0).unwrap();
        assert!(front > 10. * 30.);
        let mesh = loft
            .mesh(Tessellation {
                deflection: 0.01,
                angle: 0.2,
                triangle_limit: 1_000_000,
            })
            .unwrap();
        assert_eq!(mesh.face_indices.len() * 3, mesh.indices.len());
        assert!(mesh.face_indices.iter().all(|&face| face < 6));
    }

    #[test]
    fn crossing_or_degenerate_sections_are_refused() {
        let square = rect(4., 4., [0., 0.]);
        let mut crossing = square.clone();
        crossing.swap(1, 2);
        assert_eq!(
            SmoothLoft::interpolate(
                &[square.clone(), crossing, square.clone()],
                &[0., 5., 10.],
                1e-7
            ),
            Err(Error::Unresolved)
        );
        // A middle station far outside both ends folds the rails through
        // the corner order of the end sections.
        let flipped: Vec<_> = square.iter().rev().copied().collect();
        assert!(
            SmoothLoft::interpolate(
                &[square.clone(), flipped, square.clone()],
                &[0., 5., 10.],
                1e-7
            )
            .is_err()
        );
        assert_eq!(
            SmoothLoft::interpolate(
                &[square.clone(), square.clone(), square.clone()],
                &[0., 5., 5.],
                1e-7
            ),
            Err(Error::Invalid)
        );
        assert_eq!(
            SmoothLoft::interpolate(&[square.clone(), square.clone()], &[0., 5.], 1e-7),
            Err(Error::Invalid)
        );
    }

    #[test]
    fn clockwise_sections_keep_their_orientation_and_volume() {
        let sections: Vec<Vec<_>> = [
            rect(10., 6., [0., 0.]),
            rect(14., 8., [-2., -1.]),
            rect(10., 6., [0., 0.]),
        ]
        .into_iter()
        .map(|ring| ring.into_iter().rev().collect())
        .collect();
        let clockwise = SmoothLoft::interpolate(&sections, &[0., 15., 30.], 1e-7).unwrap();
        assert!(clockwise.clockwise());
        let forward: Vec<_> = sections
            .iter()
            .map(|ring| ring.iter().rev().copied().collect())
            .collect();
        let counter = SmoothLoft::interpolate(&forward, &[0., 15., 30.], 1e-7).unwrap();
        assert!(!counter.clockwise());
        assert!((clockwise.volume_mm3() - counter.volume_mm3()).abs() < 1e-9);
    }
}
