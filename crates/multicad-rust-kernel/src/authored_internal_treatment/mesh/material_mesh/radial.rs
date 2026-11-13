//! Local dyadic restriction of the actual ruled annulus, never a radial fit.
use super::*;
pub(super) const MAX_STRIPS: usize = 64;
pub(super) enum Verdict {
    Done,
    Radial,
    Refine,
}
pub(super) struct Spokes {
    strips: usize,
    ids: Vec<usize>,
}
fn memory(entries: usize, outer: usize, errors: usize, w: &mut Work<'_>) -> Result<(), Error> {
    // A real allocation/publication checkpoint accompanies both the planned
    // footprint and the actual retained-capacity check.
    w.step()?;
    storage::bytes(
        &[
            (entries, std::mem::size_of::<usize>()),
            (outer, std::mem::size_of::<usize>()),
            (errors, std::mem::size_of::<f64>()),
        ],
        0,
        8 * 1024 * 1024,
    )
    .map(|_| ())
}
impl Spokes {
    pub fn new(
        count: usize,
        caps: usize,
        strips: usize,
        outer_capacity: usize,
        errors: &mut Vec<f64>,
        w: &mut Work<'_>,
    ) -> Result<Self, Error> {
        w.step()?;
        if strips == 1 {
            return Ok(Self {
                strips,
                ids: Vec::new(),
            });
        }
        // The sealed source chart retains <=131072 source vertices and
        // <=32768 cap receipts. Closure adds <=65538 vertices; strips is a
        // private power of two <=64. All products/additions here fit usize
        // even on 32-bit targets. No new ceiling affects the original path.
        let entries = count * (strips - 1);
        let extra = caps * (2 * strips);
        // Include both retained maps and the enlarged endpoint error table in
        // the original stock-workspace allowance before allocating anything.
        memory(
            entries,
            outer_capacity,
            errors.capacity().max(errors.len() + extra),
            w,
        )?;
        w.charge(entries)?;
        errors.try_reserve_exact(extra).map_err(|_| Error::Limit)?;
        let mut ids = Vec::new();
        ids.try_reserve_exact(entries).map_err(|_| Error::Limit)?;
        memory(ids.capacity(), outer_capacity, errors.capacity(), w)?;
        ids.resize(entries, usize::MAX);
        w.step()?;
        Ok(Self { strips, ids })
    }
    fn point(
        &mut self,
        shell: &mut Shell,
        inner: usize,
        outer: usize,
        level: usize,
        error: [f64; 2],
        w: &mut Work<'_>,
    ) -> Result<(usize, f64), Error> {
        if level == 0 {
            return Ok((inner, error[0]));
        }
        if level == self.strips {
            return Ok((outer, error[1]));
        }
        w.charge(4)?;
        // strips is a bounded power of two, so this is the exact source
        // radial parameter. The owner is the ORIGINAL inner vertex and this
        // dyadic fraction; adjacent cap cells share the same literal point.
        let t = I::point(level as f64 / self.strips as f64);
        let one = I::point(1.) - t;
        let value: [I; 3] = std::array::from_fn(|k| {
            I::point(shell.mesh.vertices[3 * inner + k]) * one
                + I::point(shell.mesh.vertices[3 * outer + k]) * t
        });
        let slot = inner * (self.strips - 1) + level - 1;
        if self.ids[slot] == usize::MAX {
            let id = shell.mesh.vertices.len() / 3;
            w.step()?;
            shell
                .mesh
                .vertices
                .extend(value.map(|x| x.lo.midpoint(x.hi)));
            self.ids[slot] = id;
        }
        let id = self.ids[slot];
        let literal = close::error(&shell.mesh, id, value);
        // Convex transfer of the two source endpoint errors plus the actual
        // directed interpolation/rounding discrepancy. No exact welding claim.
        let bound = (one * I::point(error[0]) + t * I::point(error[1]) + I::point(literal)).hi;
        Ok((id, bound))
    }
}
fn restricted(body: &Cavity, disk: &Disk, radius: I, low: f64, high: f64) -> source::Field {
    let z = if disk.end == 0 {
        0.
    } else {
        body.full_length_mm()
    };
    // Prepared cap fields are InnerCap/OuterCap: independent of radial t.
    let inner = close::field(body, disk.edge, z, disk.mode);
    let outer = curve(body, span(disk.edge), radius, z);
    let t =
        Jet::constant(I::point(low)).add(Jet::variable(1).scale(I::point(high) - I::point(low)));
    source::Field {
        point: std::array::from_fn(|k| {
            inner.point[k].add(outer.point[k].sub(inner.point[k]).mul(t))
        }),
        normal: inner.normal,
    }
}
pub(super) fn emit(
    spokes: &mut Spokes,
    shell: &mut Shell,
    source: StockSurface<'_>,
    disk: &Disk,
    outer: ([usize; 2], [f64; 2]),
    options: Tessellation,
    w: &mut Work<'_>,
) -> Result<Verdict, Error> {
    let StockSurface { body, radius, .. } = source;
    let (outer, error) = outer;
    for strip in 0..spokes.strips {
        w.step()?;
        let mut ids = [0; 4];
        let mut errors = [0.; 4];
        for (k, (end, level)) in [(0, strip), (1, strip), (1, strip + 1), (0, strip + 1)]
            .into_iter()
            .enumerate()
        {
            (ids[k], errors[k]) = spokes.point(
                shell,
                disk.ids[end],
                outer[end],
                level,
                [disk.errors[end], error[end]],
                w,
            )?;
        }
        // Keep the original one-strip field and facet order byte-identical.
        let field = if spokes.strips == 1 {
            annulus(body, disk, radius)
        } else {
            restricted(
                body,
                disk,
                radius,
                strip as f64 / spokes.strips as f64,
                (strip + 1) as f64 / spokes.strips as f64,
            )
        };
        w.step()?;
        let chord = certificate::chord(&field, [errors[0], errors[1], errors[2].max(errors[3])]);
        if !close::quad(shell, field, ids, errors, disk.end, options, w)? {
            return Ok(if chord.is_finite() && chord > options.deflection {
                Verdict::Radial
            } else {
                Verdict::Refine
            });
        }
    }
    Ok(Verdict::Done)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn actual_source_spoke_evaluation_cancels_before_shared_vertex_publication() {
        let material = super::super::partition_closure_tests::material(8);
        let body = material.cavity();
        let options = Tessellation {
            deflection: 0.1,
            angle: 0.3,
            triangle_limit: 1_000_000,
        };
        let mut w = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        let (original, trace) = chart::traced_shell(body, options, 1, 64, &mut w)
            .unwrap()
            .unwrap();
        let mut prepared = close::prepare(original, body, options, &mut w)
            .unwrap()
            .unwrap();
        let count = prepared.shell.mesh.vertices.len() / 3;
        let disk = prepared.caps[0];
        let source = StockSurface {
            body,
            radius: material.radius_interval(),
            trace: &trace,
        };
        let mut map = vec![usize::MAX; count];
        let mut errors = vec![0.; count];
        let a = boundary_point(
            source,
            &disk,
            0,
            &mut map,
            &mut errors,
            &mut prepared.shell.mesh,
            &mut w,
        )
        .unwrap();
        let b = boundary_point(
            source,
            &disk,
            1,
            &mut map,
            &mut errors,
            &mut prepared.shell.mesh,
            &mut w,
        )
        .unwrap();
        let original = prepared.shell.mesh.clone();
        for stop in [1, 2] {
            prepared.shell.mesh = original.clone();
            w.interrupt = None;
            let mut spokes = Spokes::new(
                count,
                prepared.caps.len(),
                2,
                map.capacity(),
                &mut errors,
                &mut w,
            )
            .unwrap();
            let calls = std::cell::Cell::new(0);
            let cancel = || {
                calls.set(calls.get() + 1);
                calls.get() > stop
            };
            let mut interrupted = Work {
                remaining: 16_000_000,
                used: 0,
                interrupt: Some(&cancel),
            };
            assert!(matches!(
                emit(
                    &mut spokes,
                    &mut prepared.shell,
                    source,
                    &disk,
                    ([a.0, b.0], [a.1, b.1]),
                    options,
                    &mut interrupted
                ),
                Err(Error::Cancelled)
            ));
        }
    }
    #[test]
    fn allocation_and_publication_checkpoint_cancellation_never_returns_a_partial_cache() {
        let calls = std::cell::Cell::new(0usize);
        for stop in 0..5 {
            let cancel = || {
                calls.set(calls.get() + 1);
                calls.get() > stop
            };
            calls.set(0);
            let mut errors = Vec::new();
            let mut w = Work {
                remaining: 16_000_000,
                used: 0,
                interrupt: Some(&cancel),
            };
            assert!(matches!(
                Spokes::new(1, 1, 2, 1, &mut errors, &mut w),
                Err(Error::Cancelled)
            ));
        }
    }
    #[test]
    fn retained_capacity_is_preflighted_before_another_spoke_allocation() {
        let mut errors = Vec::<f64>::new();
        errors.try_reserve_exact(2 * 1024 * 1024).unwrap();
        let retained = errors.capacity();
        let mut w = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        assert!(matches!(
            Spokes::new(1, 1, 2, 1, &mut errors, &mut w),
            Err(Error::Limit)
        ));
        assert_eq!(errors.capacity(), retained);
        assert!(errors.is_empty());
        let cancel = || true;
        w.interrupt = Some(&cancel);
        assert!(matches!(
            Spokes::new(1, 1, 2, 1, &mut errors, &mut w),
            Err(Error::Cancelled)
        ));
        assert_eq!(errors.capacity(), retained);
    }
}
