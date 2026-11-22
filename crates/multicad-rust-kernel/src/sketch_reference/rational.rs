//! Exact finite rational Bézier edges expressed in the consumer frame.
//!
//! The edge keeps its source weights; an affine change of frame maps its
//! Euclidean poles. Plane contacts are isolated on the homogeneous depth
//! polynomial, whose Bernstein coefficients have the sign of the depth.
use super::{Contact, Method, validate};
use multicad_precision::PortableMath as _;

const MAX_DEPTH: usize = 60;
const MAX_PIECES: usize = 20_000;

/// A rational Bézier span with positive weights. Each control is
/// `[u, v, depth]` in the consumer sketch frame.
#[derive(Clone, Debug, PartialEq)]
pub struct Bezier {
    pub poles: Vec<[f64; 3]>,
    pub weights: Vec<f64>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum Outcome {
    Coincident,
    Projection {
        poles: Vec<[f64; 2]>,
        weights: Vec<f64>,
    },
    Contacts(Vec<([f64; 2], Contact)>),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Invalid,
    OffPlane,
    /// The projection degenerates to a point.
    Collapsed,
    /// The projection reverses along its own image or closes on itself.
    Retraced,
    Unresolved,
}

#[derive(Clone, Copy)]
struct Piece {
    controls: [[f64; 4]; 4],
    degree: usize,
    range: [f64; 2],
}

impl Piece {
    fn depths(&self) -> impl Iterator<Item = f64> + '_ {
        self.controls[..=self.degree].iter().map(|c| c[2] / c[3])
    }

    fn split(&self) -> [Self; 2] {
        let degree = self.degree;
        let mut left = [[0.; 4]; 4];
        let mut right = [[0.; 4]; 4];
        let mut work = self.controls;
        for level in 0..=degree {
            left[level] = work[0];
            right[degree - level] = work[degree - level];
            for i in 0..degree - level {
                work[i] = std::array::from_fn(|k| (work[i][k] + work[i + 1][k]) * 0.5);
            }
        }
        let middle = (self.range[0] + self.range[1]) * 0.5;
        [
            Self {
                controls: left,
                degree,
                range: [self.range[0], middle],
            },
            Self {
                controls: right,
                degree,
                range: [middle, self.range[1]],
            },
        ]
    }

    fn diameter(&self) -> f64 {
        let points: Vec<_> = self.controls[..=self.degree]
            .iter()
            .map(|c| [c[0] / c[3], c[1] / c[3], c[2] / c[3]])
            .collect();
        points
            .iter()
            .flat_map(|a| points.iter().map(move |b| distance(*a, *b)))
            .fold(0., f64::max)
    }
}

impl Bezier {
    fn checked(&self, tolerance: f64) -> Result<(), Error> {
        validate(&self.poles, tolerance).map_err(|_| Error::Invalid)?;
        if !(3..=4).contains(&self.poles.len())
            || self.weights.len() != self.poles.len()
            || self
                .weights
                .iter()
                .any(|w| !w.is_finite() || *w <= 0. || *w > 1e8)
        {
            return Err(Error::Invalid);
        }
        Ok(())
    }

    fn degree(&self) -> usize {
        self.poles.len() - 1
    }

    fn homogeneous(&self) -> [[f64; 4]; 4] {
        let mut controls = [[0.; 4]; 4];
        for (control, (&[u, v, d], &w)) in controls
            .iter_mut()
            .zip(self.poles.iter().zip(&self.weights))
        {
            *control = [u * w, v * w, d * w, w];
        }
        controls
    }

    fn at(&self, t: f64) -> ([f64; 3], [f64; 3]) {
        let degree = self.degree();
        let mut work = self.homogeneous();
        for level in (2..=degree).rev() {
            for i in 0..level {
                work[i] = lerp(work[i], work[i + 1], t);
            }
        }
        let [a, b] = [work[0], work[1]];
        let h = lerp(a, b, t);
        let dh: [f64; 4] = std::array::from_fn(|k| degree as f64 * (b[k] - a[k]));
        let point = [h[0] / h[3], h[1] / h[3], h[2] / h[3]];
        let derivative = std::array::from_fn(|k| (dh[k] - point[k] * dh[3]) / h[3]);
        (point, derivative)
    }

    pub fn query(&self, method: Method, tolerance: f64) -> Result<Outcome, Error> {
        self.checked(tolerance)?;
        // Positive weights keep the curve in the hull of its poles.
        let coplanar = self.poles.iter().all(|p| p[2].abs() <= tolerance);
        match method {
            Method::Intersect if coplanar => Ok(Outcome::Coincident),
            Method::Intersect => self.contacts(tolerance).map(Outcome::Contacts),
            Method::Coplanar if !coplanar => Err(Error::OffPlane),
            Method::Coplanar | Method::Project => self.project(tolerance),
        }
    }

    fn project(&self, tolerance: f64) -> Result<Outcome, Error> {
        let poles: Vec<[f64; 2]> = self.poles.iter().map(|&[u, v, _]| [u, v]).collect();
        let first = poles[0];
        let reach = poles
            .iter()
            .map(|p| (p[0] - first[0]).portable_hypot(p[1] - first[1]))
            .fold(0., f64::max);
        if reach <= tolerance {
            return Err(Error::Collapsed);
        }
        let last = poles[poles.len() - 1];
        let chord = [last[0] - first[0], last[1] - first[1]];
        let length = chord[0].portable_hypot(chord[1]);
        if length <= tolerance {
            return Err(Error::Retraced);
        }
        // A projected control polygon on one line is a segment image. It is
        // regular only when the poles advance monotonically along the chord;
        // otherwise the image doubles back over itself.
        let unit = [chord[0] / length, chord[1] / length];
        let collinear = poles.iter().all(|p| {
            ((p[0] - first[0]) * unit[1] - (p[1] - first[1]) * unit[0]).abs() <= tolerance
        });
        if collinear
            && poles.windows(2).any(|pair| {
                (pair[1][0] - pair[0][0]) * unit[0] + (pair[1][1] - pair[0][1]) * unit[1] < 0.
            })
        {
            return Err(Error::Retraced);
        }
        Ok(Outcome::Projection {
            poles,
            weights: self.weights.clone(),
        })
    }

    fn contacts(&self, tolerance: f64) -> Result<Vec<([f64; 2], Contact)>, Error> {
        // (parameter, crossing) candidates; `None` marks an unresolved touch.
        let mut found: Vec<(f64, Option<i8>)> = Vec::new();
        let mut stack = vec![(
            Piece {
                controls: self.homogeneous(),
                degree: self.degree(),
                range: [0., 1.],
            },
            0,
        )];
        // Depths within rounding of the plane are touches, as for conics.
        let scale = self
            .poles
            .iter()
            .flatten()
            .map(|v| v.abs())
            .fold(1., f64::max);
        let margin = 64. * f64::EPSILON * scale;
        let mut pieces = 0;
        while let Some((piece, depth)) = stack.pop() {
            pieces += 1;
            if pieces > MAX_PIECES {
                return Err(Error::Unresolved);
            }
            let depths: Vec<f64> = piece.depths().collect();
            let low = depths.iter().copied().fold(f64::INFINITY, f64::min);
            let high = depths.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            if low > margin || high < -margin {
                continue;
            }
            let rising = depths.windows(2).all(|pair| pair[1] >= pair[0]);
            let falling = depths.windows(2).all(|pair| pair[1] <= pair[0]);
            if (rising || falling) && low <= 0. && 0. <= high && low < high {
                // Variation diminishing: the depth is monotone on this piece
                // and changes sign once.
                found.push((
                    self.root(piece.range, rising),
                    Some(if rising { 1 } else { -1 }),
                ));
                continue;
            }
            if piece.diameter() <= tolerance * 0.25 || depth == MAX_DEPTH {
                if piece.diameter() > tolerance {
                    return Err(Error::Unresolved);
                }
                found.push(((piece.range[0] + piece.range[1]) * 0.5, None));
                continue;
            }
            for half in piece.split() {
                stack.push((half, depth + 1));
            }
        }
        found.sort_by(|a, b| a.0.total_cmp(&b.0));
        // Merge contacts that lie within the sketch precision of each other.
        // Opposite crossings or an unresolved touch inside one cluster are a
        // tangency within that precision; equal crossings are one root
        // reported on both sides of a subdivision boundary.
        let mut clusters: Vec<Vec<(f64, Option<i8>)>> = Vec::new();
        for entry in found {
            let point = self.at(entry.0).0;
            if let Some(cluster) = clusters.last_mut()
                && distance(self.at(cluster[cluster.len() - 1].0).0, point) <= tolerance
            {
                cluster.push(entry);
                continue;
            }
            clusters.push(vec![entry]);
        }
        let mut result = Vec::with_capacity(clusters.len());
        for cluster in clusters {
            let crossings: Vec<_> = cluster.iter().map(|entry| entry.1).collect();
            let crossing = match crossings[0] {
                Some(sign) if crossings.iter().all(|c| *c == Some(sign)) => Some(sign),
                _ => None,
            };
            let parameter = if crossing.is_some() {
                cluster[0].0
            } else {
                let (low, high) = (cluster[0].0, cluster[cluster.len() - 1].0);
                self.extremum(low, high)
            };
            let (point, derivative) = self.at(parameter);
            if point[2].abs() > tolerance {
                return Err(Error::Unresolved);
            }
            // A crossing at a vanishing depth derivative is not transverse.
            let speed = derivative[0]
                .portable_hypot(derivative[1])
                .portable_hypot(derivative[2]);
            let crossing = crossing.filter(|_| derivative[2].abs() > 1e-9 * speed);
            result.push((
                [point[0], point[1]],
                Contact {
                    parameter,
                    crossing,
                },
            ));
        }
        Ok(result)
    }

    /// Bisect a monotone sign change of the depth to rounding level.
    fn root(&self, [mut low, mut high]: [f64; 2], rising: bool) -> f64 {
        let depth = |t: f64| self.at(t).0[2];
        if depth(low) == 0. {
            return low;
        }
        if depth(high) == 0. {
            return high;
        }
        for _ in 0..200 {
            let middle = (low + high) * 0.5;
            if middle <= low || middle >= high {
                break;
            }
            if (depth(middle) < 0.) == rising {
                low = middle;
            } else {
                high = middle;
            }
        }
        (low + high) * 0.5
    }

    /// The parameter of smallest absolute depth on a touching cluster.
    fn extremum(&self, low: f64, high: f64) -> f64 {
        let depth = |t: f64| self.at(t).0[2].abs();
        let (mut a, mut b) = (low, high);
        for _ in 0..200 {
            if b - a <= f64::EPSILON * 4. {
                break;
            }
            let m1 = a + (b - a) / 3.;
            let m2 = b - (b - a) / 3.;
            if depth(m1) <= depth(m2) {
                b = m2;
            } else {
                a = m1;
            }
        }
        (a + b) * 0.5
    }
}

fn lerp(a: [f64; 4], b: [f64; 4], t: f64) -> [f64; 4] {
    std::array::from_fn(|i| (b[i] - a[i]).mul_add(t, a[i]))
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    (a[0] - b[0])
        .portable_hypot(a[1] - b[1])
        .portable_hypot(a[2] - b[2])
}

#[cfg(test)]
mod tests {
    use super::*;

    fn arch(depth: f64, weight: f64) -> Bezier {
        // An arch in the (u, depth) plane; a plane `depth = 0` cuts it.
        Bezier {
            poles: vec![[0., 0., -depth], [5., 0., 10. - depth], [10., 0., -depth]],
            weights: vec![1., weight, 1.],
        }
    }

    fn top(weight: f64) -> f64 {
        // Peak of the symmetric conic at t = 1/2.
        10. * weight / (1. + weight)
    }

    #[test]
    fn transverse_contacts_tangent_touch_and_miss() {
        for weight in [1., 0.7] {
            let Outcome::Contacts(hits) = arch(3., weight).query(Method::Intersect, 1e-6).unwrap()
            else {
                panic!()
            };
            assert_eq!(hits.len(), 2);
            assert_eq!(hits[0].1.crossing, Some(1));
            assert_eq!(hits[1].1.crossing, Some(-1));
            assert!((hits[0].1.parameter + hits[1].1.parameter - 1.).abs() < 1e-12);
            let Outcome::Contacts(touch) = arch(top(weight), weight)
                .query(Method::Intersect, 1e-6)
                .unwrap()
            else {
                panic!()
            };
            assert_eq!(touch.len(), 1, "{touch:?}");
            assert_eq!(touch[0].1.crossing, None);
            assert!((touch[0].1.parameter - 0.5).abs() < 1e-6);
            assert_eq!(
                arch(top(weight) + 1., weight).query(Method::Intersect, 1e-6),
                Ok(Outcome::Contacts(vec![]))
            );
        }
    }

    #[test]
    fn projections_keep_weights_and_refuse_collapse_and_retrace() {
        let curve = Bezier {
            poles: vec![[0., 0., 0.], [5., 10., 0.], [10., 0., 0.]],
            weights: vec![1., 0.7, 1.],
        };
        assert_eq!(
            curve.query(Method::Coplanar, 1e-6),
            Ok(Outcome::Projection {
                poles: vec![[0., 0.], [5., 10.], [10., 0.]],
                weights: vec![1., 0.7, 1.],
            })
        );
        assert_eq!(
            curve.query(Method::Intersect, 1e-6),
            Ok(Outcome::Coincident)
        );
        let edge_on = Bezier {
            poles: vec![[0., 0., 0.], [0., 0., 5.], [0., 0., 10.]],
            weights: vec![1., 1., 1.],
        };
        assert_eq!(edge_on.query(Method::Project, 1e-6), Err(Error::Collapsed));
        assert_eq!(edge_on.query(Method::Coplanar, 1e-6), Err(Error::OffPlane));
        let folded = Bezier {
            poles: vec![[0., 0., 0.], [10., 0., 5.], [0., 0., 10.]],
            weights: vec![1., 1., 1.],
        };
        assert_eq!(folded.query(Method::Project, 1e-6), Err(Error::Retraced));
        let bad = Bezier {
            poles: vec![[0., 0., 0.], [10., 0., 5.], [0., 0., 10.]],
            weights: vec![1., 0., 1.],
        };
        assert_eq!(bad.query(Method::Project, 1e-6), Err(Error::Invalid));
    }
}
