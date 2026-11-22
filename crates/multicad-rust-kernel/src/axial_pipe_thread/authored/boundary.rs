//! Source-owned raw boundary, deliberately separate from material admission.
use super::*;
use crate::radial_tensor::{Axis, BoundaryInput, BoundaryOwnership, Curve, Edge, Face, Surface};
use crate::tapered_thread::RetainedTensor;
use std::f64::consts::TAU;
mod rounded;
type S = Series<11>;
const DEGREE: usize = 9;
pub struct RawAuthoredBoundary {
    pub ownership: BoundaryOwnership,
    /// Euclidean positional source-to-raw-boundary bound, including cap trim
    /// endpoint projection. It is not a normal or decoded mass certificate.
    pub source_error_mm: f64,
    /// The stronger L1 Cartesian tensor positional bound.
    pub tensor_error_mm: f64,
    /// Euclidean analytic cone/circle support and coefficient reserve.
    pub analytic_error_mm: f64,
    pub cap_faces: [usize; 2],
    pub band_faces: Vec<usize>,
    /// Actual root support faces. Rounded roots are tensor bands, not cones.
    pub root_faces: Vec<usize>,
}
struct Work<'a> {
    remaining: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Work<'_> {
    fn step(&mut self) -> Result<(), Error> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        self.remaining = self.remaining.checked_sub(1).ok_or(Error::Limit)?;
        Ok(())
    }
}
struct Curves {
    spans: usize,
    poles: Vec<[f64; 4]>,
    error: [f64; 4],
    nodes: Vec<f64>,
}
fn choose(n: usize, k: usize) -> f64 {
    (0..k).fold(1., |v, i| v * (n - i) as f64 / (i + 1) as f64)
}
fn power(v: I, n: usize) -> I {
    (0..n).fold(I::point(1.), |x, _| x * v)
}
fn node(i: usize, n: usize, length: f64) -> f64 {
    if i == 0 {
        0.
    } else if i == n {
        length
    } else {
        length * i as f64 / n as f64
    }
}
fn curve_values(t: &AuthoredThread, band: Option<usize>, x: S) -> Result<[S; 4], Error> {
    let w = Series::interval(I::point(2.) * pi() / I::point(t.spec.pitch_mm));
    let (mut angle, scale) = if let Some(band) = band {
        let (q, a) = t.profile(band, x)?;
        (-w * q, a)
    } else {
        (w * x, x)
    };
    if matches!(t.spec.section, AuthoredSection::BsptRounded { .. })
        && ((band == Some(0) && x.c[0].lo == 0. && x.c[0].hi == 0.)
            || (band == Some(4) && x.c[0].lo == 1. && x.c[0].hi == 1.))
    {
        // The authored endpoints q=0 and q=p have exactly the same angular
        // phase. Retain this true periodic identity before rounding controls;
        // no literal imported support or positive junction is rewritten.
        angle.c[0] = I::point(0.);
    }
    let [sin, cos] = angle.sin_cos();
    let out = [cos, sin, scale * cos, scale * sin];
    for value in out {
        value.checked()?;
    }
    Ok(out)
}
fn fit(
    t: &AuthoredThread,
    band: Option<usize>,
    budget: f64,
    work: &mut Work<'_>,
) -> Result<Curves, Error> {
    let length = if band.is_some() { 1. } else { t.spec.length_mm };
    let whole = curve_values(t, band, S::variable(I::new(0., length)))?;
    let m = whole.map(|v| v.c[10].magnitude());
    let bound = m.into_iter().fold(0., f64::max);
    let estimate = (bound * length.powi(10) / (1024. * budget))
        .powf(0.1)
        .ceil()
        .max(1.);
    if !estimate.is_finite() || estimate > 256. {
        return Err(Error::Limit);
    }
    let mut spans = estimate as usize;
    loop {
        work.step()?;
        let count = DEGREE * spans + 1;
        let mut poles = vec![[0.; 4]; count];
        let mut error = [0_f64; 4];
        for cell in 0..spans {
            work.step()?;
            let lo = node(cell, spans, length);
            let hi = node(cell + 1, spans, length);
            let width = I::point(hi) - I::point(lo);
            let endpoints = [
                curve_values(t, band, S::variable(I::point(lo)))?,
                curve_values(t, band, S::variable(I::point(hi)))?,
            ];
            let derivative = curve_values(t, band, S::variable(I::new(lo, hi)))?;
            for i in 0..=DEGREE {
                let side = usize::from(i >= 5);
                let jet = if side == 0 { i } else { DEGREE - i };
                let delta = if side == 0 { width } else { -width };
                let mut value = [I::point(0.); 4];
                for j in 0..=jet {
                    let factor =
                        I::point(choose(jet, j)) / I::point(choose(DEGREE, j)) * power(delta, j);
                    for axis in 0..4 {
                        value[axis] = value[axis] + factor * endpoints[side][axis].c[j];
                    }
                }
                for axis in 0..4 {
                    let value = value[axis].checked()?;
                    let chosen = value.midpoint();
                    let round = (I::point(chosen) - value).magnitude();
                    let remainder = (I::point(derivative[axis].c[10].magnitude())
                        * power(width, 10)
                        / I::point(1024.))
                    .hi;
                    error[axis] = error[axis].max((I::point(round) + I::point(remainder)).hi);
                    let index = DEGREE * cell + i;
                    if cell > 0 && i == 0 {
                        identical(&poles[index][axis], &chosen)?;
                    }
                    poles[index][axis] = chosen;
                }
            }
        }
        if error.into_iter().all(|e| e <= budget) {
            return Ok(Curves {
                spans,
                poles,
                error,
                nodes: (0..=spans).map(|i| node(i, spans, length)).collect(),
            });
        }
        if spans >= 256 {
            return Err(Error::Unresolved);
        }
        spans = (spans * 2).min(256);
    }
}
struct Net {
    tensor: RetainedTensor,
    error: f64,
}
fn net(t: &AuthoredThread, u: &Curves, z: &Curves, work: &mut Work<'_>) -> Result<Net, Error> {
    let h = z.poles.len();
    let count = net_size(u.poles.len(), h)?;
    let mut poles = Vec::with_capacity(count);
    let mut rounding = 0_f64;
    let k = I::point(t.spec.radial_taper);
    for a in &u.poles {
        work.step()?;
        for (j, b) in z.poles.iter().enumerate() {
            let a = a.map(I::point);
            let b = b.map(I::point);
            let x = a[2] * b[0] - a[3] * b[1] + k * (a[0] * b[2] - a[1] * b[3]);
            let y = a[2] * b[1] + a[3] * b[0] + k * (a[0] * b[3] + a[1] * b[2]);
            let cell = (j / DEGREE).min(z.spans - 1);
            let local = j - cell * DEGREE;
            let zz = if local == 0 {
                I::point(z.nodes[cell])
            } else if local == DEGREE {
                I::point(z.nodes[cell + 1])
            } else {
                I::point(z.nodes[cell])
                    + (I::point(z.nodes[cell + 1]) - I::point(z.nodes[cell]))
                        * I::point(local as f64)
                        / I::point(DEGREE as f64)
            };
            let values = [x.checked()?, y.checked()?, zz.checked()?];
            let point = values.map(I::midpoint);
            rounding = rounding.max(
                values
                    .into_iter()
                    .zip(point)
                    .fold(I::point(0.), |e, (v, p)| {
                        e + I::point((I::point(p) - v).magnitude())
                    })
                    .checked()?
                    .hi,
            );
            poles.push(point);
        }
    }
    let eu0 = u.error[0].max(u.error[1]);
    let eu1 = u.error[2].max(u.error[3]);
    let ev0 = z.error[0].max(z.error[1]);
    let ev1 = z.error[2].max(z.error[3]);
    let source = I::point(4.)
        * (I::point(eu1) + (I::point(t.spec.major_radius_mm) + I::point(eu1)) * I::point(ev0))
        + I::point(4.)
            * I::point(t.spec.radial_taper.abs())
            * (I::point(t.spec.length_mm) * I::point(eu0)
                + (I::point(1.) + I::point(eu0)) * I::point(ev1));
    let error = (source + I::point(rounding)).checked()?.hi;
    Ok(Net {
        tensor: RetainedTensor::new([u.spans, z.spans], poles)?,
        error,
    })
}
fn identical<T: PartialEq + ?Sized>(a: &T, b: &T) -> Result<(), Error> {
    if a == b {
        Ok(())
    } else {
        Err(Error::Unresolved)
    }
}
fn net_size(width: usize, height: usize) -> Result<usize, Error> {
    let count = width.checked_mul(height).ok_or(Error::Limit)?;
    (count <= 200_000).then_some(count).ok_or(Error::Limit)
}
fn error_budget(error: f64, tolerance: f64) -> Result<(), Error> {
    if error > tolerance / 16. {
        Err(Error::Unresolved)
    } else {
        Ok(())
    }
}
struct Graph(BoundaryInput);
impl Graph {
    fn vertex(&mut self, p: [f64; 3]) -> usize {
        let i = self.0.vertices.len();
        self.0.vertices.push(p);
        i
    }
    fn surface(&mut self, s: Surface) -> usize {
        let i = self.0.surfaces.len();
        self.0.surfaces.push(s);
        i
    }
    fn edge(&mut self, vertices: [usize; 2], curve: Curve) -> usize {
        let forward = matches!(&curve, Curve::Spline { .. });
        let c = self.0.curves.len();
        self.0.curves.push(curve);
        let i = self.0.edges.len();
        self.0.edges.push(Edge {
            vertices,
            curve: c,
            forward,
        });
        i
    }
    fn face(
        &mut self,
        surface: usize,
        mut edges: Vec<(usize, bool)>,
        reverse: bool,
        forward: bool,
    ) {
        if reverse {
            edges.reverse();
            for e in &mut edges {
                e.1 = !e.1;
            }
        }
        self.0.faces.push(Face {
            surface,
            outer: edges,
            forward,
        });
    }
}
impl AuthoredThread {
    /// Raw local coefficients and abstract incidence only. Binding, finite
    /// trims, occupied material and literal embedding remain separate gates.
    pub fn raw_boundary_ownership(
        &self,
        tolerance: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<RawAuthoredBoundary, Error> {
        if !tolerance.is_finite() || tolerance <= 0. {
            return Err(Error::Invalid);
        }
        if let Some(limit) = self.spec.crest_radius_limit_mm {
            let maximum = I::point(self.spec.major_radius_mm)
                + max(
                    I::point(0.),
                    I::point(self.spec.radial_taper) * I::point(self.spec.length_mm),
                );
            if maximum.hi > limit {
                return Err(Error::Unsupported);
            }
        }
        let tolerance = tolerance.min(self.tolerance).min(1e-7);
        let budget = tolerance / (32. * (self.spec.major_radius_mm + self.spec.length_mm + 1.));
        let mut work = Work {
            remaining: 50_000,
            interrupt,
        };
        work.step()?;
        let axial = fit(self, None, budget, &mut work)?;
        let mut nets = Vec::new();
        let mut error = 0_f64;
        for band in 0..self.band_count() {
            let u = fit(self, Some(band), budget, &mut work)?;
            let n = net(self, &u, &axial, &mut work)?;
            error = error.max(n.error);
            nets.push(n);
        }
        error_budget(error, tolerance)?;
        if matches!(self.spec.section, AuthoredSection::BsptRounded { .. }) {
            return rounded::build(self, nets, error, tolerance, interrupt);
        }
        let mut graph = Graph(BoundaryInput {
            vertices: vec![],
            curves: vec![],
            edges: vec![],
            surfaces: vec![],
            faces: vec![],
        });
        let mut surfaces = Vec::new();
        let mut rails = Vec::new();
        let mut corners = Vec::<[usize; 2]>::new();
        let mut rims = Vec::<[usize; 2]>::new();
        let mut previous_rail = None::<Vec<[f64; 3]>>;
        for (band, n) in nets.into_iter().enumerate() {
            work.step()?;
            let dims = n.tensor.spans().map(|s| DEGREE * s + 1);
            let [w, h] = dims;
            let points = n.tensor.poles();
            let left = if band == 0 {
                [graph.vertex(points[0]), graph.vertex(points[h - 1])]
            } else {
                let ids = corners[corners.len() - 1];
                let previous = previous_rail.as_deref().expect("previous authored band");
                identical(previous, &points[..h])?;
                ids
            };
            let right = [
                graph.vertex(points[(w - 1) * h]),
                graph.vertex(points[w * h - 1]),
            ];
            if band == 0 {
                corners.push(left);
                rails.push(graph.edge(
                    left,
                    Curve::Spline {
                        spans: n.tensor.spans()[1],
                        poles: points[..h].to_vec(),
                    },
                ));
            }
            corners.push(right);
            previous_rail = Some(points[(w - 1) * h..].to_vec());
            rails.push(graph.edge(
                right,
                Curve::Spline {
                    spans: n.tensor.spans()[1],
                    poles: points[(w - 1) * h..].to_vec(),
                },
            ));
            rims.push(std::array::from_fn(|end| {
                graph.edge(
                    [left[end], right[end]],
                    Curve::Spline {
                        spans: n.tensor.spans()[0],
                        poles: (0..w).map(|i| points[i * h + end * (h - 1)]).collect(),
                    },
                )
            }));
            surfaces.push(graph.surface(Surface::Tensor(n.tensor)));
        }
        let root = I::point(self.spec.major_radius_mm) - I::point(self.depth());
        let mut gap = [Vec::new(), Vec::new()];
        for end in 0..2 {
            let z = if end == 0 { 0. } else { self.spec.length_mm };
            let radius = (root + I::point(self.spec.radial_taper) * I::point(z)).midpoint();
            let a = graph.0.vertices[corners[3][end]];
            let b = graph.0.vertices[corners[0][end]];
            let from = a[1].atan2(a[0]);
            let to = b[1].atan2(b[0]);
            let sweep = -(from - to).rem_euclid(TAU);
            let angle = from + sweep / 2.;
            let middle = graph.vertex([radius * angle.cos(), radius * angle.sin(), z]);
            let ids = [corners[3][end], middle, corners[0][end]];
            for pair in ids.windows(2) {
                gap[end].push(graph.edge(
                    [pair[0], pair[1]],
                    Curve::Circle {
                        axis: Axis {
                            origin: [0., 0., z],
                            normal: [0., 0., 1.],
                            reference: [1., 0., 0.],
                        },
                        radius_mm: radius,
                    },
                ));
            }
        }
        for end in 0..2 {
            let mut edges = (0..3).map(|b| (rims[b][end], true)).collect::<Vec<_>>();
            edges.extend(gap[end].iter().map(|&e| (e, true)));
            let z = if end == 0 { 0. } else { self.spec.length_mm };
            let s = graph.surface(Surface::Plane(Axis {
                origin: [0., 0., z],
                normal: [0., 0., if end == 0 { -1. } else { 1. }],
                reference: [1., 0., 0.],
            }));
            graph.face(s, edges, end == 1, true);
        }
        let tensor_error = error;
        let z = if self.spec.radial_taper < 0. {
            self.spec.length_mm
        } else {
            0.
        };
        let radius = (root + I::point(self.spec.radial_taper) * I::point(z)).midpoint();
        let semi = self.spec.radial_taper.abs().atan();
        let tangent = I::point(semi).sin() / I::point(semi).cos();
        let analytic_error =
            (I::point((tangent - I::point(self.spec.radial_taper.abs())).magnitude())
                * I::point(self.spec.length_mm)
                + I::point(64. * f64::EPSILON)
                    * (I::point(self.spec.major_radius_mm)
                        + I::point(self.spec.length_mm)
                        + I::point(1.)))
            .checked()?
            .hi;
        error = (I::point(2.) * I::point(tensor_error) + I::point(analytic_error))
            .checked()?
            .hi;
        let cone = graph.surface(Surface::Cone {
            axis: Axis {
                origin: [0., 0., z],
                normal: [0., 0., if self.spec.radial_taper < 0. { -1. } else { 1. }],
                reference: [1., 0., 0.],
            },
            radius_mm: radius,
            semi_angle_rad: semi,
        });
        for band in 0..3 {
            graph.face(
                surfaces[band],
                vec![
                    (rims[band][0], true),
                    (rails[band + 1], true),
                    (rims[band][1], false),
                    (rails[band], false),
                ],
                true,
                false,
            );
        }
        let mut edges = gap[0].iter().map(|&e| (e, true)).collect::<Vec<_>>();
        edges.push((rails[0], true));
        edges.extend(gap[1].iter().rev().map(|&e| (e, false)));
        edges.push((rails[3], false));
        graph.face(cone, edges, true, true);
        error_budget(error, tolerance)?;
        Ok(RawAuthoredBoundary {
            ownership: BoundaryOwnership::new(graph.0, tolerance, interrupt)?,
            source_error_mm: error,
            tensor_error_mm: tensor_error,
            analytic_error_mm: analytic_error,
            cap_faces: [0, 1],
            band_faces: vec![2, 3, 4],
            root_faces: vec![5],
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn actual_receipt_limits_and_shared_rail_tamper_are_checked_before_reuse() {
        let t = super::super::tests::source(1. / 32.);
        let raw = t.raw_boundary_ownership(1e-7, None).unwrap();
        let tensors = raw
            .ownership
            .input()
            .surfaces
            .iter()
            .filter_map(|s| {
                if let Surface::Tensor(t) = s {
                    Some(t)
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        assert_eq!(tensors.len(), 3);
        let tensor = tensors[0];
        let h = DEGREE * tensor.spans()[1] + 1;
        let rail = &tensor.poles()[..h];
        let mut tampered = rail.to_vec();
        tampered[h / 2][0] = tampered[h / 2][0].next_up();
        assert_eq!(identical(rail, &tampered[..]), Err(Error::Unresolved));
        assert!(identical(rail, rail).is_ok());
        assert!(error_budget(raw.tensor_error_mm, 1e-7).is_ok());
        assert_eq!(
            error_budget(raw.tensor_error_mm, raw.tensor_error_mm * 8.),
            Err(Error::Unresolved)
        );
        let width = DEGREE * tensor.spans()[0] + 1;
        assert_eq!(net_size(width, h), Ok(tensor.poles().len()));
        assert_eq!(net_size(200_001, 1), Err(Error::Limit));
        assert_eq!(net_size(usize::MAX, 2), Err(Error::Limit));
        // A high OD limit is physically inactive, so it must leave all finite
        // tensor controls and incidence unchanged.
        let inactive = AuthoredThread::new(
            AuthoredSpec {
                crest_radius_limit_mm: Some(10.),
                ..t.spec()
            },
            1e-7,
            1e-6,
        )
        .unwrap()
        .raw_boundary_ownership(1e-7, None)
        .unwrap();
        assert_eq!(raw.ownership.input(), inactive.ownership.input());
    }
}
