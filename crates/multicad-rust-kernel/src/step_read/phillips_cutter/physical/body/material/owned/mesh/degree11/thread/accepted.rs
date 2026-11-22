//! Invocation-local reuse of a fully proved ORIGINAL root partition. The
//! source/root ordinal and caller quality are fixed by the enclosing call;
//! only identical local levels AND literal vertex bits permit reuse.
use super::*;
pub(super) struct Accepted {
    levels: [usize; 2],
    points: Vec<P>,
    error: f64,
}
fn id(grid: &Grid, patch: Patch, u: usize, v: usize) -> usize {
    ((patch.offset[0] + u) % grid.shape[0]) * (grid.shape[1] + 1) + patch.offset[1] + v
}
impl Accepted {
    pub(super) fn bytes(&self) -> usize {
        self.points.capacity() * std::mem::size_of::<P>()
    }
    pub(super) fn capture(
        grid: &Grid,
        patch: Patch,
        levels: [usize; 2],
        work: &mut Work<'_>,
    ) -> R<Self> {
        let [nu, nv] = grid.divisions;
        let count = (nu + 1)
            .checked_mul(nv + 1)
            .ok_or(PhillipsReadError::Limit)?;
        work.storage(count.saturating_mul(std::mem::size_of::<P>()) + SCRATCH)?;
        let mut points = Vec::with_capacity(count);
        for u in 0..=nu {
            for v in 0..=nv {
                work.charge(32)?;
                points.push(
                    grid.vertices[id(grid, patch, u, v)].ok_or(PhillipsReadError::Unsupported)?,
                );
            }
        }
        let result = Self {
            levels,
            points,
            error: grid.error,
        };
        work.storage(result.bytes() + SCRATCH)?;
        Ok(result)
    }
    pub(super) fn reuse(
        &self,
        grid: &mut Grid,
        patch: Patch,
        levels: [usize; 2],
        work: &mut Work<'_>,
    ) -> R<bool> {
        work.charge(4)?;
        if levels != self.levels {
            return Ok(false);
        }
        let [nu, nv] = grid.divisions;
        // Validate ALL preexisting shared vertices before changing any slot.
        // Vacant slots can use these already-proved literal vertices directly.
        for u in 0..=nu {
            for v in 0..=nv {
                work.charge(64)?;
                let expected = self.points[u * (nv + 1) + v];
                if grid.vertices[id(grid, patch, u, v)]
                    .is_some_and(|actual| actual.map(f64::to_bits) != expected.map(f64::to_bits))
                {
                    return Ok(false);
                }
            }
        }
        for u in 0..=nu {
            for v in 0..=nv {
                work.charge(32)?;
                let index = id(grid, patch, u, v);
                grid.vertices[index] = Some(self.points[u * (nv + 1) + v]);
            }
        }
        for u in 0..nu {
            for v in 0..nv {
                work.charge(128)?;
                let corners = [[u, v], [u + 1, v], [u + 1, v + 1], [u, v + 1]]
                    .map(|[u, v]| id(grid, patch, u, v) as u32);
                for mut triangle in [[0, 1, 2], [0, 2, 3]] {
                    if !patch.forward {
                        triangle.swap(1, 2);
                    }
                    grid.indices.extend(triangle.map(|j| corners[j]));
                    grid.faces.push(patch.face);
                }
            }
        }
        grid.error = grid.error.max(self.error);
        Ok(true)
    }
}
