//! Complete finite-band comparison of two exact-periodic sources. Shared
//! normalized u and v coordinates retain the real phase/axial expressions;
//! no sampled mesh, cached-rate Thread or endpoint fit defines either side.
use super::*;
#[derive(Clone, Copy, Debug)]
pub struct PeriodicThreadCorrespondence {
    pub position_mm: f64,
    pub du_mm: f64,
    pub dv_mm: f64,
    pub work_used: usize,
}
fn term(p: [f64; 2]) -> I {
    if p[1] == 0. {
        I::point(p[0])
    } else {
        I::point(p[0]) + I::point(p[1])
    }
}
fn difference(a: [f64; 2], b: [f64; 2]) -> I {
    if a == b {
        I::point(0.)
    } else {
        term(a) - term(b)
    }
}
fn abs(x: I) -> I {
    I::point(x.magnitude())
}
fn caps_valid(caps: [[f64; 2]; 2], length: f64) -> bool {
    use num_rational::BigRational as Q;
    if caps
        .into_iter()
        .flatten()
        .any(|x| !x.is_finite() || x.abs() > 1e7)
    {
        return false;
    }
    let value = caps.map(|p| Q::from_float(p[0]).unwrap() + Q::from_float(p[1]).unwrap());
    value[0] >= Q::from_integer(0.into())
        && value[0] < value[1]
        && value[1] <= Q::from_float(length).unwrap()
}
impl PeriodicThreadSource {
    /// C0/Du/Dv for each corresponding full finite band. v spans each actual
    /// pair of exact two-term cap stations, rather than nominal thread length.
    pub fn finite_correspondence(
        &self,
        other: &Self,
        source_caps_mm: [[f64; 2]; 2],
        target_caps_mm: [[f64; 2]; 2],
        maximum_work: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<PeriodicThreadCorrespondence, PhillipsReadError> {
        let mut work = Work {
            remaining: maximum_work,
            interrupt,
            cancelled: Cell::new(false),
        };
        let result = (|| {
            // Fixed four-band expressions, <128 interval operations/band,
            // each <=16 primitive steps. Include bounded exact cap arithmetic.
            work.step(16_384)?;
            if maximum_work > WORK
                || !caps_valid(source_caps_mm, self.length)
                || !caps_valid(target_caps_mm, other.length)
            {
                return Err(ReadError::Invalid);
            }
            require(self.hand == other.hand)?;
            let rotation = rotation::Rotation::new(&mut work)?;
            let tau = I::point(2.) * rotation.pi(&mut work)?;
            let rate = tau * I::point(self.hand) / I::point(self.pitch);
            let other_rate = tau * I::point(other.hand) / I::point(other.pitch);
            let delta_rate = if self.pitch == other.pitch {
                I::point(0.)
            } else {
                tau * I::point(self.hand)
                    * (I::point(1.) / I::point(self.pitch) - I::point(1.) / I::point(other.pitch))
            };
            let caps = source_caps_mm.map(term);
            let other_caps = target_caps_mm.map(term);
            let dz = I::point(
                abs(difference(source_caps_mm[0], target_caps_mm[0]))
                    .hi
                    .max(abs(difference(source_caps_mm[1], target_caps_mm[1])).hi),
            );
            let height = caps[1] - caps[0];
            let other_height = other_caps[1] - other_caps[0];
            let dh = if source_caps_mm == target_caps_mm {
                I::point(0.)
            } else {
                height - other_height
            };
            let phase = self.phase_node_expressions_mm();
            let other_phase = other.phase_node_expressions_mm();
            let r = self.radius_nodes_mm();
            let other_r = other.radius_nodes_mm();
            let mut out = PeriodicThreadCorrespondence {
                position_mm: 0.,
                du_mm: 0.,
                dv_mm: 0.,
                work_used: 0,
            };
            for band in 0..4 {
                work.step(1)?;
                let dp = I::point(
                    abs(difference(phase[band], other_phase[band]))
                        .hi
                        .max(abs(difference(phase[band + 1], other_phase[band + 1])).hi),
                );
                let zphase = I::new(other_caps[0].lo, other_caps[1].hi)
                    - I::new(term(other_phase[band]).lo, term(other_phase[band + 1]).hi);
                let angle = abs(rate) * (dz + dp) + abs(delta_rate) * abs(zphase);
                let dr = I::point(
                    abs(I::point(r[band]) - I::point(other_r[band]))
                        .hi
                        .max(abs(I::point(r[band + 1]) - I::point(other_r[band + 1])).hi),
                );
                let radius = I::point(other_r[band].max(other_r[band + 1]));
                let slope = I::point(r[band + 1]) - I::point(r[band]);
                let other_slope = I::point(other_r[band + 1]) - I::point(other_r[band]);
                let width = term(phase[band + 1]) - term(phase[band]);
                let other_width = term(other_phase[band + 1]) - term(other_phase[band]);
                let delta_width = if phase[band..=band + 1] == other_phase[band..=band + 1] {
                    I::point(0.)
                } else {
                    width - other_width
                };
                let u_angle = abs(rate * width);
                let other_u_angle = abs(other_rate * other_width);
                let du_angle = abs(rate) * abs(delta_width) + abs(delta_rate) * abs(other_width);
                let v_angle = abs(rate * height);
                let other_v_angle = abs(other_rate * other_height);
                let dv_angle = abs(rate) * abs(dh) + abs(delta_rate) * abs(other_height);
                let position = dr + radius * angle + dz;
                let du = abs(slope - other_slope)
                    + abs(other_slope) * angle
                    + dr * u_angle
                    + radius * (du_angle + other_u_angle * angle);
                let dv = dr * v_angle + radius * (dv_angle + other_v_angle * angle) + abs(dh);
                require(
                    [position, du, dv]
                        .into_iter()
                        .all(|x| x.valid() && x.hi.is_finite()),
                )?;
                out.position_mm = out.position_mm.max(position.hi);
                out.du_mm = out.du_mm.max(du.hi);
                out.dv_mm = out.dv_mm.max(dv.hi);
            }
            work.step(0)?;
            out.work_used = maximum_work - work.remaining;
            Ok(out)
        })();
        result.map_err(|e| bound_error(e, &work))
    }
}

#[cfg(test)]
mod tests;
