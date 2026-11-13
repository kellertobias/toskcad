//! Conforming dyadic restriction of the original expanded shoulder only.
use super::*;
pub(super) const MAX_STRIPS: usize = 64;
pub(super) struct Spokes {
    strips: usize,
    ids: Vec<usize>,
}
pub(super) enum Verdict {
    Done,
    Radial,
    Refine,
}
impl Spokes {
    pub fn new(
        count: usize,
        strips: usize,
        maps: [usize; 2],
        w: &mut Work<'_>,
    ) -> Result<Self, Error> {
        let entries = count * (strips - 1);
        let check = |capacity| {
            storage::bytes(
                &[
                    (maps[0], std::mem::size_of::<usize>()),
                    (maps[1], std::mem::size_of::<usize>()),
                    (capacity, std::mem::size_of::<usize>()),
                ],
                0,
                4 * 1024 * 1024,
            )
            .map(|_| ())
        };
        check(entries)?;
        let mut ids = Vec::new();
        if entries > 0 {
            w.charge(entries)?;
            ids.try_reserve_exact(entries).map_err(|_| Error::Limit)?;
            check(ids.capacity())?;
            ids.resize(entries, usize::MAX);
            w.step()?;
        }
        Ok(Self { strips, ids })
    }
    fn point(
        &mut self,
        shell: &mut Shell,
        inner: usize,
        outer: usize,
        level: usize,
        errors: [f64; 2],
        w: &mut Work<'_>,
    ) -> Result<(usize, f64), Error> {
        if level == 0 {
            return Ok((inner, errors[0]));
        }
        if level == self.strips {
            return Ok((outer, errors[1]));
        }
        w.charge(4)?;
        let t = I::point(level as f64 / self.strips as f64);
        let one = I::point(1.) - t;
        let value: [I; 3] = std::array::from_fn(|k| {
            I::point(shell.mesh.vertices[3 * inner + k]) * one
                + I::point(shell.mesh.vertices[3 * outer + k]) * t
        });
        let slot = inner * (self.strips - 1) + level - 1;
        if self.ids[slot] == usize::MAX {
            w.step()?;
            self.ids[slot] = shell.mesh.vertices.len() / 3;
            shell
                .mesh
                .vertices
                .extend(value.map(|x| x.lo.midpoint(x.hi)));
        }
        let id = self.ids[slot];
        let bound = (one * I::point(errors[0])
            + t * I::point(errors[1])
            + I::point(error(&shell.mesh, id, value)))
        .hi;
        Ok((id, bound))
    }
    pub fn emit(
        &mut self,
        shell: &mut Shell,
        body: &Cavity,
        edge: Edge,
        outer: ([usize; 2], [f64; 2]),
        options: Tessellation,
        w: &mut Work<'_>,
    ) -> Result<Verdict, Error> {
        let (outer, errors) = outer;
        let face = 2 + 4 * body.thread.spec().starts + u32::from(body.entry.is_some());
        for strip in 0..self.strips {
            let mut ids = [0; 4];
            let mut e = [0.; 4];
            for (k, (end, level)) in [(0, strip), (1, strip), (1, strip + 1), (0, strip + 1)]
                .into_iter()
                .enumerate()
            {
                (ids[k], e[k]) = self.point(
                    shell,
                    edge.vertices[end],
                    outer[end],
                    level,
                    [edge.errors[end], errors[end]],
                    w,
                )?;
            }
            let field = if self.strips == 1 {
                field(body, edge, body.tail.unwrap().station, Surface::Shoulder)
            } else {
                field_at(
                    body,
                    edge,
                    body.tail.unwrap().station,
                    Surface::Shoulder,
                    Jet::constant(I::point(strip as f64 / self.strips as f64))
                        .add(Jet::variable(1).scale(I::point(1. / self.strips as f64))),
                )
            };
            if !quad(
                shell,
                source::Field {
                    point: field.point,
                    normal: field.normal,
                },
                ids,
                e,
                face,
                options,
                w,
            )? {
                w.step()?;
                let chord = certificate::chord(&field, [e[0], e[1], e[2].max(e[3])]);
                return Ok(if chord.is_finite() && chord > options.deflection {
                    Verdict::Radial
                } else {
                    Verdict::Refine
                });
            }
        }
        Ok(Verdict::Done)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn actual_radial_spoke_interpolation_cancels_before_any_new_publication() {
        let body = super::super::super::tests::source(false, 1, Some(0.4), true);
        let options = super::super::super::tests::options();
        let mut w = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        let mut shell = chart::shell(&body, options, 1, 64, &mut w)
            .unwrap()
            .unwrap();
        let edge = shell.top[0];
        let count = shell.mesh.vertices.len() / 3;
        let mut map = vec![usize::MAX; count];
        let outer = super::super::outer(
            &body,
            edge,
            body.tail.unwrap().station,
            &mut map,
            &mut shell.mesh,
            &mut w,
        )
        .unwrap();
        let mut spokes = Spokes::new(count, 2, [map.capacity(), map.capacity()], &mut w).unwrap();
        let before = shell.mesh.vertices.clone();
        let mut cancelled = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: Some(&|| true),
        };
        assert!(matches!(
            spokes.emit(&mut shell, &body, edge, outer, options, &mut cancelled),
            Err(Error::Cancelled)
        ));
        assert_eq!(shell.mesh.vertices, before);
    }
    #[test]
    fn bounded_source_vertex_plan_cannot_exceed_original_closure_workspace() {
        let mut w = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        // A legal upper source-vertex count cannot request an unaccounted
        // dyadic table: refuse its footprint before allocating that table.
        assert!(matches!(
            Spokes::new(131072, MAX_STRIPS, [131072; 2], &mut w),
            Err(Error::Limit)
        ));
    }
}
