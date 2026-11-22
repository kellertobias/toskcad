//! Whole finite source coedges and point owners. Periodic crest coedges use
//! the original phase q4..pitch+q1, with exact one-turn endpoint ownership.
use super::*;
#[derive(Clone, Copy, Debug)]
pub struct PointBounds {
    position: [[f64; 2]; 3],
    work: usize,
}
impl PointBounds {
    pub fn position_mm(&self) -> [[f64; 2]; 3] {
        self.position
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
}
fn tau() -> I {
    I::new(
        std::f64::consts::TAU.next_down(),
        std::f64::consts::TAU.next_up(),
    )
}
impl Nut {
    fn angular_cuts(&self) -> [I; 5] {
        let r = point(self.relief.radius_mm);
        let a = point(self.width[3]);
        let root = (r * r - a * a).sqrt();
        let alpha = (a / root).atan();
        let pi = tau() / point(2.);
        [alpha, pi - alpha, pi + alpha, tau() - alpha, tau() + alpha]
    }
    fn thread_point(&self, q: I, r: I, z: I) -> [I; 3] {
        let gain = tau() / point(self.thread.spec().pitch_mm)
            * point(if self.thread.spec().left_hand {
                -1.
            } else {
                1.
            });
        let angle = (z - q) * gain;
        let xy = [r * angle.cos(), r * angle.sin()];
        let m = self.phase_rotation.map(|row| row.map(point));
        [
            m[0][0] * xy[0] + m[0][1] * xy[1],
            m[1][0] * xy[0] + m[1][1] * xy[1],
            z,
        ]
    }
    fn vertex(&self, v: Vertex) -> Result<[I; 3], Error> {
        Ok(match v {
            Vertex::Stock { end, corner } if end < 2 && corner < 8 => {
                let p = self.stock.base()[0][corner as usize];
                [
                    if end == 0 {
                        point(-self.bore_from_start_mm)
                    } else {
                        point(self.stock.height_mm()) - point(self.bore_from_start_mm)
                    },
                    point(p[0]),
                    point(p[1]),
                ]
            }
            Vertex::Crown(i) | Vertex::ReliefStart(i) if i < 4 => {
                let r = point(self.relief.radius_mm);
                let a = point(self.width[3]);
                let root = (r * r - a * a).sqrt();
                [
                    if i == 0 || i == 3 { root } else { -root },
                    if i < 2 { a } else { -a },
                    point(if matches!(v, Vertex::Crown(_)) {
                        self.height[3]
                    } else {
                        self.relief.start_mm
                    }),
                ]
            }
            Vertex::Thread { end, junction } if end < 2 && junction < 4 => {
                let b = self.thread.bands()[junction as usize];
                self.thread_point(
                    point(b.phase_mm[1]),
                    point(b.radius_mm[1]),
                    point(if end == 0 { 0. } else { self.relief.start_mm }),
                )
            }
            _ => return Err(Error::Invalid),
        })
    }
    pub fn vertex_bounds(
        &self,
        v: Vertex,
        maximum_work: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<PointBounds, Error> {
        checkpoint(interrupt)?;
        if maximum_work < 2048 {
            return Err(Error::Limit);
        }
        let position = self.vertex(v)?;
        for p in position {
            valid(p)?;
        }
        checkpoint(interrupt)?;
        Ok(PointBounds {
            position: position.map(|v| [v.lo, v.hi]),
            work: 2048,
        })
    }
    pub fn curve_bounds(
        &self,
        edge: Edge,
        range: [f64; 2],
        maximum_work: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<PointBounds, Error> {
        checkpoint(interrupt)?;
        if maximum_work < 2048 {
            return Err(Error::Limit);
        }
        if !range[0].is_finite()
            || !range[1].is_finite()
            || range[0] < 0.
            || range[1] > 1.
            || range[0] > range[1]
        {
            return Err(Error::Invalid);
        }
        let ends = edge.endpoints()?;
        let t = I::new(range[0], range[1]);
        let position = if range == [0., 0.] {
            self.vertex(ends[0])?
        } else if range == [1., 1.] {
            self.vertex(ends[1])?
        } else {
            match edge {
                Edge::StockEnd { .. } | Edge::Rail(_) | Edge::SplitRail { .. } => {
                    let a = self.vertex(ends[0])?;
                    let b = self.vertex(ends[1])?;
                    std::array::from_fn(|i| a[i] + (b[i] - a[i]) * t)
                }
                Edge::Crown(i) | Edge::ReliefStart(i) => {
                    let cuts = self.angular_cuts();
                    let angle = cuts[i as usize] + (cuts[i as usize + 1] - cuts[i as usize]) * t;
                    let r = point(self.relief.radius_mm);
                    let x = r * angle.cos();
                    let y = r * angle.sin();
                    let z = if matches!(edge, Edge::ReliefStart(_)) {
                        point(self.relief.start_mm)
                    } else if i == 0 || i == 2 {
                        let sign = if i == 0 { 1. } else { -1. };
                        let slope = (point(self.height[3]) - point(self.height[2]))
                            / (point(self.width[2]) - point(self.width[3]));
                        point(self.height[3]) - (point(sign) * y - point(self.width[3])) * slope
                    } else {
                        point(self.height[3])
                    };
                    [x, y, z]
                }
                Edge::ThreadRing { end, role } => {
                    let band = [3, 4, 1, 2][role as usize];
                    let b = self.thread.bands()[band];
                    let q0 = point(b.phase_mm[0]);
                    let q1 = if role == 1 {
                        point(self.thread.spec().pitch_mm)
                            + point(self.thread.bands()[0].phase_mm[1])
                    } else {
                        point(b.phase_mm[1])
                    };
                    let r0 = point(b.radius_mm[0]);
                    let r1 = point(b.radius_mm[1]);
                    self.thread_point(
                        q0 + (q1 - q0) * t,
                        r0 + (r1 - r0) * t,
                        point(if end == 0 { 0. } else { self.relief.start_mm }),
                    )
                }
                Edge::Helix(i) => {
                    let b = self.thread.bands()[i as usize];
                    self.thread_point(
                        point(b.phase_mm[1]),
                        point(b.radius_mm[1]),
                        t * point(self.relief.start_mm),
                    )
                }
            }
        };
        for p in position {
            valid(p)?;
        }
        checkpoint(interrupt)?;
        Ok(PointBounds {
            position: position.map(|v| [v.lo, v.hi]),
            work: 2048,
        })
    }
    /// Whole normalized-coedge second derivative. Together with the source point
    /// enclosures this gives the standard M/(8N²) piecewise-linear chord bound.
    pub fn curve_second_derivative_bound_mm(&self, edge: Edge) -> Result<f64, Error> {
        edge.endpoints()?;
        let scale = (self.maximum_thread / point(self.thread.maximum_radius_mm())).hi;
        let result = match edge {
            Edge::StockEnd { .. } | Edge::Rail(_) | Edge::SplitRail { .. } => point(0.),
            Edge::Crown(i) | Edge::ReliefStart(i) => {
                let cuts = self.angular_cuts();
                let span = cuts[i as usize + 1] - cuts[i as usize];
                let factor = if matches!(edge, Edge::Crown(0) | Edge::Crown(2)) {
                    point(1.)
                        + (point(self.height[3]) - point(self.height[2]))
                            / (point(self.width[2]) - point(self.width[3]))
                } else {
                    point(1.)
                };
                point(self.relief.radius_mm) * span * span * factor
            }
            Edge::Helix(i) => {
                let b = self.thread.bands()[i as usize];
                let gain = tau() * point(self.relief.start_mm) / point(self.thread.spec().pitch_mm);
                point(scale) * point(b.radius_mm[1]) * gain * gain
            }
            Edge::ThreadRing { role, .. } => {
                let band = [3, 4, 1, 2][role as usize];
                let b = self.thread.bands()[band];
                let span = if role == 1 {
                    point(self.thread.spec().pitch_mm) + point(self.thread.bands()[0].phase_mm[1])
                        - point(b.phase_mm[0])
                } else {
                    point(b.phase_mm[1]) - point(b.phase_mm[0])
                };
                let angular = tau() * span / point(self.thread.spec().pitch_mm);
                let dr = (point(b.radius_mm[1]) - point(b.radius_mm[0])).magnitude();
                point(scale)
                    * (point(self.thread.maximum_radius_mm()) * angular * angular
                        + point(2.) * point(dr) * angular)
            }
        };
        Ok(valid(result)?.hi)
    }
}
