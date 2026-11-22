//! Whole original-generator isocurves of extrusion and revolution supports.
//! UV LINE parameters remain free; bounded generator knot domains are retained.
use super::{Budget, Proof, Result, affine, analytic, vector};
use crate::brep::{Curve, Edge, Surface};
use crate::tapered_thread::interval::I;
type V = [I; 3];
const CONTROL_WORK: usize = 256;
const WORKSPACE: usize = 4 * 1024 * 1024;

fn add(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i] + b[i])
}
fn sub(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i] - b[i])
}
fn scale(a: V, b: I) -> V {
    a.map(|x| x * b)
}
fn dot(a: V, b: V) -> I {
    (0..3).fold(I::point(0.), |s, i| s + a[i] * b[i])
}
fn cross(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn valid(v: V) -> bool {
    v.into_iter().all(I::valid)
}

#[derive(Clone, Copy)]
struct Map {
    translation: V,
    columns: [V; 3],
}
impl Map {
    fn vector(self, p: V) -> V {
        (0..3).fold([I::point(0.); 3], |s, i| {
            add(s, scale(self.columns[i], p[i]))
        })
    }
    fn point(self, p: V) -> V {
        add(self.translation, self.vector(p))
    }
    fn translation(v: V) -> Self {
        Self {
            translation: v,
            columns: [
                vector([1., 0., 0.]),
                vector([0., 1., 0.]),
                vector([0., 0., 1.]),
            ],
        }
    }
    fn rotation(origin: V, axis: V, angle: f64) -> Self {
        let (cosine, sine) = (I::point(angle).cos(), I::point(angle).sin());
        // The actual retained Surface operator uses this complete Rodrigues
        // expression. No additional unit-axis or rigid-placement assumption.
        let turn = |p: V| {
            let along = scale(axis, dot(p, axis));
            let radial = sub(p, along);
            add(
                along,
                add(scale(radial, cosine), scale(cross(axis, radial), sine)),
            )
        };
        let columns = [
            vector([1., 0., 0.]),
            vector([0., 1., 0.]),
            vector([0., 0., 1.]),
        ]
        .map(turn);
        let mut result = Self {
            translation: [I::point(0.); 3],
            columns,
        };
        result.translation = sub(origin, result.vector(origin));
        result
    }
}

fn generator(
    curve: &Curve,
    map: Map,
    rate: f64,
    edge: &Edge,
    tol: f64,
    budget: &mut Budget<'_>,
) -> Result<Proof> {
    match curve {
        Curve::Line { origin, direction } => affine::line(
            map.point(vector(*origin)),
            scale(map.vector(vector(*direction)), I::point(rate)),
            edge,
            tol,
            budget,
        ),
        Curve::Nurbs { closed: true, .. } => Ok(Proof::Unsupported),
        Curve::Nurbs {
            curve,
            closed: false,
        } => {
            let count = curve.points().len();
            let bytes = count
                .checked_mul(std::mem::size_of::<V>())
                .and_then(|n| n.checked_add(std::mem::size_of::<Vec<V>>()));
            if count > crate::brep::nurbs::MAX_CONTROL_POINTS || bytes.is_none_or(|n| n > WORKSPACE)
            {
                return Err(super::ReadError::Limit);
            }
            budget.charge(64)?;
            let mut controls = Vec::new();
            controls
                .try_reserve_exact(count)
                .map_err(|_| super::ReadError::Limit)?;
            if controls
                .capacity()
                .saturating_mul(std::mem::size_of::<V>())
                .saturating_add(std::mem::size_of::<Vec<V>>())
                > WORKSPACE
            {
                return Err(super::ReadError::Limit);
            }
            for &point in curve.points() {
                // 9 interval products, 12 sums, point conversion, finite
                // validation, index/copy/loop work fit this bounded debit.
                budget.charge(CONTROL_WORK)?;
                let p = map.point(vector(point));
                if !valid(p) {
                    return Ok(Proof::Unsupported);
                }
                controls.push(p);
            }
            // Actual knots/positive weights are unchanged. The helper proves
            // finite coverage and affine/reversed parameter correspondence.
            // The selected map must increase the ACTUAL UV LINE parameter;
            // reversing its moving coordinate reverses the original generator.
            affine::nurbs(
                curve.degree(),
                curve.knots(),
                curve.weights(),
                &controls,
                rate > 0.,
                edge,
                tol,
                budget,
            )
        }
        Curve::Circle { frame, radius } => analytic::periodic(
            map.point(vector(frame.origin)),
            map.vector(scale(vector(frame.x), I::point(*radius))),
            map.vector(scale(vector(frame.y), I::point(*radius * rate.signum()))),
            edge,
            tol,
            budget,
        ),
        Curve::Ellipse {
            frame,
            semiaxes: [a, b],
        } => analytic::periodic(
            map.point(vector(frame.origin)),
            map.vector(scale(vector(frame.x), I::point(*a))),
            map.vector(scale(vector(frame.y), I::point(*b * rate.signum()))),
            edge,
            tol,
            budget,
        ),
        _ => Ok(Proof::Unsupported),
    }
}

pub(super) fn certify(
    surface: &Surface,
    uv: &Curve,
    edge: &Edge,
    tol: f64,
    budget: &mut Budget<'_>,
) -> Result<Proof> {
    // Fixed stack-only map setup: <=3 Rodrigues basis evaluations, complete
    // argument checks, trig pairs, vector copies, and branch bookkeeping.
    budget.charge(2048)?;
    if !tol.is_finite()
        || tol <= 0.
        || edge.range.into_iter().any(|x| !x.is_finite())
        || edge.range[0] >= edge.range[1]
    {
        return Ok(Proof::Unsupported);
    }
    let Curve::Line { origin, direction } = uv else {
        return Ok(Proof::Unsupported);
    };
    if origin.iter().chain(direction).any(|x| !x.is_finite())
        || origin[2] != 0.
        || direction[2] != 0.
    {
        return Ok(Proof::Unsupported);
    }
    let along_u = direction[0] != 0. && direction[1] == 0.;
    let along_v = direction[0] == 0. && direction[1] != 0.;
    if !along_u && !along_v {
        return Ok(Proof::Unsupported);
    }
    match surface {
        Surface::Extrusion {
            curve,
            direction: axis,
        } if along_u => generator(
            curve,
            Map::translation(scale(vector(*axis), I::point(origin[1]))),
            direction[0],
            edge,
            tol,
            budget,
        ),
        Surface::Extrusion {
            curve,
            direction: axis,
        } => {
            // A fixed generator parameter must be in its ORIGINAL natural
            // domain; at() never repairs it by clamping or extrapolating.
            let p = affine::at(curve, origin[0], budget)?;
            affine::line(
                add(p, scale(vector(*axis), I::point(origin[1]))),
                scale(vector(*axis), I::point(direction[1])),
                edge,
                tol,
                budget,
            )
        }
        Surface::Revolution {
            curve,
            origin: center,
            axis,
        } if along_v => generator(
            curve,
            Map::rotation(vector(*center), vector(*axis), origin[0]),
            direction[1],
            edge,
            tol,
            budget,
        ),
        Surface::Revolution {
            curve,
            origin: center,
            axis,
        } => {
            let p = affine::at(curve, origin[1], budget)?;
            let delta = sub(p, vector(*center));
            let along = scale(vector(*axis), dot(delta, vector(*axis)));
            let radial = sub(delta, along);
            // U ranges over all real angles, hence this SAME actual operator
            // covers its entire periodic image, not merely sampled points.
            analytic::periodic(
                add(vector(*center), along),
                radial,
                scale(
                    cross(vector(*axis), radial),
                    I::point(direction[0].signum()),
                ),
                edge,
                tol,
                budget,
            )
        }
        _ => Ok(Proof::Unsupported),
    }
}

#[cfg(test)]
mod tests;
