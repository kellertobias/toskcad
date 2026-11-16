//! Literal control-hull homotopy and actual endpoint polygon winding.
use super::*;
#[derive(Clone, Copy, Debug)]
pub struct WindingBounds {
    /// Start and end cap, each with its actual two directed loops.
    pub winding: [[i32; 2]; 2],
    pub certified_segments: usize,
    pub work_used: usize,
}
impl TopologyPlan<'_, '_, '_> {
    /// Each complete polynomial span and its chord occupy the SAME strict
    /// origin halfplane. Their straight homotopy therefore avoids the origin;
    /// winding is the actual literal endpoint polygon's ray count, not a
    /// source/header prediction. Shared closed endpoints were earned by prepare.
    pub fn winding_bounds(
        &self,
        maximum_work: usize,
        maximum_bytes: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<WindingBounds, Error> {
        let mut w = Work::new(maximum_work, interrupt)?;
        w.charge(self.work_used() + 1)?;
        if maximum_bytes == 0 || maximum_bytes > 64 * 1024 * 1024 {
            return Err(Error::Invalid);
        }
        let mut memory = Memory::new(self.retained_bytes(), maximum_bytes)?;
        // Fixed local coefficient intervals, halfplane candidates and receipt.
        memory.add(1024)?;
        let mut winding = [[0; 2]; 2];
        let mut segments = 0;
        for f in self.faces() {
            let FaceGeometry::Cap(end) = f.geometry else {
                continue;
            };
            for (loop_index, wire) in f.loops.iter().enumerate() {
                let mut value = 0;
                for c in wire {
                    let EdgeGeometry::Row(row) = self.edges()[c.edge].geometry else {
                        return Err(Error::Unresolved);
                    };
                    let row = &self.surface().rows()[row];
                    let count = row.knots().len() - 1;
                    for j in 0..count {
                        let span = if c.forward { j } else { count - 1 - j };
                        let poles = &row.poles()[11 * span..=11 * (span + 1)];
                        halfplane(poles, &mut w)?;
                        let [a, b] = if c.forward {
                            [poles[0], poles[11]]
                        } else {
                            [poles[11], poles[0]]
                        };
                        value += crossing(a, b, &mut w)?;
                        segments += 1;
                    }
                }
                winding[end as usize][loop_index] = value;
            }
            let values = winding[end as usize];
            if values[0].abs() != 1 || values[1] != -values[0] {
                return Err(Error::Unresolved);
            }
        }
        w.charge(1)?;
        Ok(WindingBounds {
            winding,
            certified_segments: segments,
            work_used: w.used,
        })
    }
}
fn literal(x: f64) -> I {
    if x == 0. || (x.fract() == 0. && x.abs() <= 9_007_199_254_740_992.) {
        I::point(x)
    } else {
        I::new(x.next_down(), x.next_up())
    }
}
fn halfplane(poles: &[[f64; 3]], w: &mut Work<'_>) -> Result<(), Error> {
    for [x, y] in [
        [1., 0.],
        [1., 1.],
        [0., 1.],
        [-1., 1.],
        [-1., 0.],
        [-1., -1.],
        [0., -1.],
        [1., -1.],
    ] {
        let mut fits = true;
        for p in poles {
            w.charge(1)?;
            let value = literal(p[0]) * I::point(x) + literal(p[1]) * I::point(y);
            if !value.valid() || value.lo <= 0. {
                fits = false;
                break;
            }
        }
        if fits {
            return Ok(());
        }
    }
    Err(Error::Unresolved)
}
pub(super) fn crossing(a: [f64; 3], b: [f64; 3], w: &mut Work<'_>) -> Result<i32, Error> {
    w.charge(1)?;
    let up = a[1] <= 0. && b[1] > 0.;
    let down = b[1] <= 0. && a[1] > 0.;
    if !up && !down {
        return Ok(0);
    }
    let cross = literal(a[0]) * literal(b[1]) - literal(a[1]) * literal(b[0]);
    if !cross.valid() || (cross.lo <= 0. && cross.hi >= 0.) {
        return Err(Error::Unresolved);
    }
    Ok(if up && cross.lo > 0. {
        1
    } else if down && cross.hi < 0. {
        -1
    } else {
        0
    })
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn literal_homotopy_and_half_open_ray_ownership_refuse_ambiguous_geometry() {
        let mut w = Work::new(16_000_000, None).unwrap();
        halfplane(&[[1., -0.5, 0.], [2., 0.5, 0.]], &mut w).unwrap();
        assert!(matches!(
            halfplane(&[[1., 0., 0.], [-1., 0., 0.]], &mut w),
            Err(Error::Unresolved)
        ));
        assert!(matches!(
            halfplane(&[[f64::NAN, 0., 0.]], &mut w),
            Err(Error::Unresolved)
        ));
        assert_eq!(crossing([1., 0., 0.], [1., 1., 0.], &mut w).unwrap(), 1);
        assert_eq!(crossing([1., 1., 0.], [1., 0., 0.], &mut w).unwrap(), -1);
        assert_eq!(crossing([-1., -1., 0.], [-1., 1., 0.], &mut w).unwrap(), 0);
        assert_eq!(crossing([-1., 1., 0.], [-1., -1., 0.], &mut w).unwrap(), 0);
        assert_eq!(crossing([1., 1., 0.], [2., 2., 0.], &mut w).unwrap(), 0);
        assert!(matches!(
            crossing([-1., -1., 0.], [1., 1., 0.], &mut w),
            Err(Error::Unresolved)
        ));
        assert!(matches!(
            crossing([f64::INFINITY, -1., 0.], [1., 1., 0.], &mut w),
            Err(Error::Unresolved)
        ));
    }
}
