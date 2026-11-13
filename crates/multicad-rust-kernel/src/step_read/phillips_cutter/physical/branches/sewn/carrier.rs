//! Evaluation caches for expressions over the ORIGINAL retained decimal data.
//! Interval fields enclose that expression; they are not freely chosen source
//! parameters or replacement points. Fixed lift parameters are explicit new
//! source choices, with a separately measured transfer from the raw lift family.
use super::*;
#[derive(Clone, Copy)]
pub(super) enum Carrier {
    Line {
        origin: [I; 3],
        vector: [I; 3],
        parameters: [f64; 2],
        transfer: [I; 2],
    },
    Arc {
        frame: support::Frame,
        radii: [I; 2],
        start: f64,
        sweep: f64,
        transfer: [I; 2],
    },
}
fn interval(v: [f64; 2]) -> I {
    I::new(v[0], v[1])
}
fn choose(v: [f64; 2]) -> f64 {
    v[0].midpoint(v[1])
}
fn discrepancy(value: f64, range: [f64; 2]) -> I {
    I::point((I::point(value) - interval(range)).magnitude())
}
impl Carrier {
    pub(super) fn new(source: lifts::AnalyticCarrier) -> Self {
        match source {
            lifts::AnalyticCarrier::Line(lift) => {
                let chosen = lift.parameters.map(choose);
                let error = [0, 1].map(|i| discrepancy(chosen[i], lift.parameters[i]));
                Self::Line {
                    origin: lift.origin,
                    vector: lift.vector,
                    parameters: chosen,
                    transfer: [
                        lift.magnitude * I::point(error[0].hi.max(error[1].hi)),
                        lift.magnitude * (error[0] + error[1]),
                    ],
                }
            }
            lifts::AnalyticCarrier::Arc { lift, radii } => {
                let radius = I::point(radii[0].hi.max(radii[1].hi));
                let start = choose(lift.start);
                let sweep = choose(lift.sweep);
                let ds = discrepancy(start, lift.start);
                let dw = discrepancy(sweep, lift.sweep);
                let angle_error = ds + dw;
                // True AXIS2 frames are orthonormal. Angular speed and second
                // angular derivative of this ellipse are <=max(a,b). These
                // compare unit-edge parameterizations, not arc-length speeds.
                let transfer = [
                    radius * angle_error,
                    radius * (dw + I::point(interval(lift.sweep).magnitude()) * angle_error),
                ];
                Self::Arc {
                    frame: lift.frame,
                    radii,
                    start,
                    sweep,
                    transfer,
                }
            }
        }
    }
    pub(super) fn transfer(&self) -> [I; 2] {
        match self {
            Self::Line { transfer, .. } | Self::Arc { transfer, .. } => *transfer,
        }
    }
    pub(super) fn domain_envelope(&self) -> ([I; 3], [I; 3]) {
        match self {
            Self::Line { .. } => self.evaluate(I::new(0., 1.)),
            Self::Arc {
                frame,
                radii,
                sweep,
                ..
            } => {
                // All runtime sin/cos enclosures are clipped to[-1,1]. This
                // complete box bounds every validated unit subdomain using
                // the SAME add/multiply intermediate operation order.
                arc_fields(*frame, *radii, I::new(-1., 1.), I::new(-1., 1.), *sweep)
            }
        }
    }
    pub(super) fn evaluate(&self, s: I) -> ([I; 3], [I; 3]) {
        match self {
            Self::Line {
                origin,
                vector,
                parameters,
                ..
            } => {
                let width = I::point(parameters[1]) - I::point(parameters[0]);
                let parameter = I::point(parameters[0]) + width * s;
                (
                    std::array::from_fn(|i| origin[i] + vector[i] * parameter),
                    vector.map(|x| x * width),
                )
            }
            Self::Arc {
                frame,
                radii,
                start,
                sweep,
                ..
            } => {
                let theta = I::point(*start) + I::point(*sweep) * s;
                let (sin, cos) = (theta.sin(), theta.cos());
                arc_fields(*frame, *radii, sin, cos, *sweep)
            }
        }
    }
}
fn arc_fields(
    frame: support::Frame,
    radii: [I; 2],
    sin: I,
    cos: I,
    sweep: f64,
) -> ([I; 3], [I; 3]) {
    let local = [radii[0] * cos, radii[1] * sin];
    let derivative = [
        -radii[0] * sin * I::point(sweep),
        radii[1] * cos * I::point(sweep),
    ];
    (
        std::array::from_fn(|i| {
            frame.origin[i] + frame.axes[0][i] * local[0] + frame.axes[1][i] * local[1]
        }),
        std::array::from_fn(|i| {
            frame.axes[0][i] * derivative[0] + frame.axes[1][i] * derivative[1]
        }),
    )
}
