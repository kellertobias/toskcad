//! Measured C0 cubic/support envelopes only, not an exact intersection, finite
//! face filling, endpoint/analytic-curve, winding, embedding or mass certificate.
use super::*;
mod polynomial;
use polynomial::{add, hull, square};
#[cfg(test)]
mod tests;

/// Directed implicit-function data under literal AXIS2 semantics. Units of
/// value/gradient follow the polynomial (plane mm, quadric mm², torus mm⁴).
/// This is not a normalized distance, surface normal or finite-body property.
#[derive(Clone, Copy, Debug)]
pub struct ImplicitSupportEnclosure {
    pub value: [f64; 2],
    pub world_gradient: [[f64; 2]; 3],
    pub local_point: [[f64; 2]; 3],
    pub world_reconstruction: [[f64; 2]; 3],
}
impl Surface {
    pub fn implicit_enclosure(
        &self,
        world_bounds: [[f64; 2]; 3],
    ) -> std::result::Result<ImplicitSupportEnclosure, PhillipsReadError> {
        let point = world_bounds.map(|[a, b]| I::new(a, b));
        if point.into_iter().any(|i| !i.valid()) {
            return Err(PhillipsReadError::Invalid);
        }
        let chart = Chart::new(self).map_err(PhillipsReadError::from)?;
        let local = chart.project(point);
        let value = chart.implicit(local);
        let gradient = chart.world_gradient(local);
        let reconstructed = chart.world(local);
        if !value.valid()
            || local
                .into_iter()
                .chain(gradient)
                .chain(reconstructed)
                .any(|i| !i.valid())
        {
            return Err(PhillipsReadError::Limit);
        }
        let bounds = |i: I| [i.lo, i.hi];
        Ok(ImplicitSupportEnclosure {
            value: bounds(value),
            world_gradient: gradient.map(bounds),
            local_point: local.map(bounds),
            world_reconstruction: reconstructed.map(bounds),
        })
    }
}

#[derive(Clone, Copy, Debug)]
pub struct CubicSupportOptions {
    /// Admission limit only. The returned measured bound is derived from the
    /// actual literal controls/supports, not assigned from this requested value.
    pub maximum_distance_mm: f64,
    pub work_limit: usize,
}
#[derive(Clone, Copy, Debug)]
pub struct EdgeSupportBound {
    pub edge_index: usize,
    /// Same owning-face order as the immutable parsed edge's `owners` array.
    pub owner_distances_mm: [f64; 2],
}
/// Borrows the exact immutable arena it measures; cannot qualify another file.
#[derive(Debug)]
pub struct CubicSupportBounds<'a> {
    boundary: &'a Boundary,
    edges: Vec<EdgeSupportBound>,
    work_used: usize,
}
impl<'a> CubicSupportBounds<'a> {
    pub fn boundary(&self) -> &'a Boundary {
        self.boundary
    }
    pub fn edges(&self) -> &[EdgeSupportBound] {
        &self.edges
    }
    pub fn work_used(&self) -> usize {
        self.work_used
    }
    /// Receipt header and owned edge allocation; excludes the borrowed arena.
    /// Saturation remains a conservative charge for a caller's resource limit.
    pub fn retained_capacity_bytes(&self) -> usize {
        std::mem::size_of::<Self>().saturating_add(
            self.edges
                .capacity()
                .saturating_mul(std::mem::size_of::<EdgeSupportBound>()),
        )
    }
}
pub(super) mod chart;
use chart::{Chart, Kind};
impl Chart {
    fn distance(&self, p: [[I; 3]; 4]) -> Result<f64> {
        let [x, y, z] = [0, 1, 2].map(|i| p.map(|p| p[i]));
        let zbox = hull(&z);
        let radial = add(&square(&x), &square(&y));
        let rho = hull(&radial).sqrt();
        let residual = self.implicit_curve(p);
        let denominator = match self.kind {
            Kind::Plane => return finite_bound(hull(&residual).magnitude()),
            Kind::Cylinder(r) => rho + r,
            Kind::Cone(r, k) => {
                let q = z.map(|v| r + k * v);
                let lower = hull(&q);
                if lower.lo <= 0. {
                    return Err(ReadError::Unsupported);
                }
                rho + lower
            }
            Kind::Torus(a, b) => {
                let q = ((rho - a).square() + zbox.square()).sqrt();
                (q + b) * ((rho + a).square() + zbox.square() - b.square())
            }
        };
        if !denominator.valid() || denominator.lo <= 0. {
            return Err(ReadError::Unsupported);
        }
        finite_bound((I::point(hull(&residual).magnitude()) / denominator).hi)
    }
}

fn finite_bound(value: f64) -> Result<f64> {
    if value.is_finite() && value >= 0. {
        Ok(value)
    } else {
        Err(ReadError::Limit)
    }
}
impl Boundary {
    pub fn cubic_support_bounds(
        &self,
        options: CubicSupportOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<CubicSupportBounds<'_>, PhillipsReadError> {
        if !options.maximum_distance_mm.is_finite() || options.maximum_distance_mm <= 0. {
            return Err(PhillipsReadError::Invalid);
        }
        let mut work = Work {
            remaining: options.work_limit,
            interrupt,
            cancelled: Cell::new(false),
        };
        let result = (|| {
            let mut edges = Vec::new();
            for (edge_index, edge) in self.edges().iter().enumerate() {
                let Curve::CubicIntersection { spans, .. } = &edge.curve else {
                    continue;
                };
                let mut owner_distances_mm = [0_f64; 2];
                for (side, owner) in edge.owners.into_iter().enumerate() {
                    let chart = Chart::new(&self.faces()[owner].support)?;
                    for span in spans {
                        let mut pending = vec![(span.poles_mm.map(|p| chart.local(p)), 0)];
                        while let Some((p, depth)) = pending.pop() {
                            work.step(512)?;
                            let bound = chart.distance(p)?;
                            if bound <= options.maximum_distance_mm {
                                owner_distances_mm[side] = owner_distances_mm[side].max(bound);
                            } else {
                                if depth == 24 {
                                    return Err(ReadError::Limit);
                                }
                                pending.extend(polynomial::split(p).map(|p| (p, depth + 1)));
                            }
                        }
                    }
                }
                edges.push(EdgeSupportBound {
                    edge_index,
                    owner_distances_mm,
                });
            }
            Ok(CubicSupportBounds {
                boundary: self,
                edges,
                work_used: options.work_limit - work.remaining,
            })
        })();
        if work.cancelled.get() {
            Err(PhillipsReadError::Cancelled)
        } else {
            result.map_err(PhillipsReadError::from)
        }
    }
}
