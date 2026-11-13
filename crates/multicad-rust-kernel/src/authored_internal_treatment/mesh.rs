//! Actual source-facet candidates. Native literal embedding/manifold proof is
//! mandatory before these detached arrays may become manufacturing output.
use super::*;
use crate::display::{DisplayMesh, Tessellation};
mod certificate;
mod chart;
mod clip;
mod close;
mod exact;
mod jet;
mod material_mesh;
mod owners;
mod regular;
mod source;
mod storage;
mod trace;
use source::Kind;
#[derive(Clone, Debug)]
pub struct Candidate {
    pub mesh: DisplayMesh,
    pub source_chord_bound_mm: f64,
    pub source_normal_angle_bound_rad: f64,
    pub work_used: usize,
}
mod sampling;
use sampling::{Attempt, Grid};
struct Work<'a> {
    remaining: usize,
    used: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Work<'_> {
    fn step(&mut self) -> Result<(), Error> {
        self.charge(1)
    }
    fn charge(&mut self, amount: usize) -> Result<(), Error> {
        if self.interrupt.is_some_and(|check| check()) {
            return Err(Error::Cancelled);
        }
        if amount > self.remaining {
            return Err(Error::Limit);
        }
        self.remaining -= amount;
        self.used += amount;
        Ok(())
    }
}
#[derive(Clone, Copy)]
struct Node {
    band: usize,
    start: u32,
}
/// Qualifies every actual facet against the immutable shared analytic source.
/// Small-entry threshold domains are currently an explicit unresolved case.
/// This function alone does not certify literal global embedding/manifoldness.
pub fn tessellate_candidate(
    body: &Cavity,
    options: Tessellation,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Candidate, Error> {
    options.check().map_err(|_| Error::Invalid)?;
    let mut work = Work {
        remaining: 16_000_000,
        used: 0,
        interrupt,
    };
    work.step()?;
    if body
        .entry
        .is_some_and(|entry| entry.radius <= body.thread.maximum_radius_mm())
    {
        return Err(Error::Unresolved);
    }
    let mut grid = Grid::new();
    loop {
        let repartition = match chart::grid_attempt(body, options, &grid, &mut work)? {
            Attempt::Complete(shell) => {
                if let Some(mut candidate) = close::close(shell, body, options, &mut work)? {
                    candidate.work_used = work.used;
                    return Ok(candidate);
                }
                false
            }
            Attempt::Refine => false,
            Attempt::Repartition => true,
        };
        grid.advance(repartition)?;
    }
}

fn wall_face(body: &Cavity, node: Node) -> u32 {
    let n = body.thread.spec().starts;
    let shift = [0, 0, 0, 1, 1][node.band];
    let group = (node.start + shift) % n;
    let group = if body.thread.spec().left_hand {
        group
    } else {
        (n - group) % n
    };
    2 + 4 * group + [1, 2, 3, 0, 1][node.band]
}
fn face(body: &Cavity, node: Node, kind: Kind) -> u32 {
    let base = 2 + 4 * body.thread.spec().starts;
    match kind {
        Kind::Wall => wall_face(body, node),
        Kind::Cone => base,
    }
}
#[derive(Clone, Copy)]
struct Qualified {
    ids: [usize; 3],
    reverse: bool,
    chord: f64,
    angle: f64,
}
fn qualify(
    field: &source::Field,
    points: &[[f64; 3]],
    errors: &[f64],
    ids: [usize; 3],
    options: Tessellation,
    work: &mut Work<'_>,
) -> Result<Option<Qualified>, Error> {
    Ok(
        match qualify_attempt(field, points, errors, ids, options, work)? {
            Attempt::Complete(proof) => Some(proof),
            Attempt::Refine | Attempt::Repartition => None,
        },
    )
}
fn qualify_attempt(
    field: &source::Field,
    points: &[[f64; 3]],
    errors: &[f64],
    ids: [usize; 3],
    options: Tessellation,
    work: &mut Work<'_>,
) -> Result<Attempt<Qualified>, Error> {
    work.step()?;
    let chord = certificate::chord(field, ids.map(|i| errors[i]));
    if !chord.is_finite() || chord > options.deflection {
        return Ok(Attempt::Refine);
    }
    let Some((reverse, angle)) =
        certificate::oriented_triangle(ids.map(|i| points[i]), field.normal, options.angle)
    else {
        work.charge(3)?;
        return Ok(
            if sampling::short_literal_edge(ids.map(|i| points[i]), ids.map(|i| errors[i])) {
                Attempt::Repartition
            } else {
                Attempt::Refine
            },
        );
    };
    Ok(Attempt::Complete(Qualified {
        ids,
        reverse,
        chord,
        angle,
    }))
}

fn append_qualified(
    mesh: &mut DisplayMesh,
    indices: &[usize],
    proof: Qualified,
    face: u32,
    options: Tessellation,
) -> Result<(), Error> {
    if mesh.face_indices.len() == options.triangle_limit {
        return Err(Error::Limit);
    }
    let ids = proof.ids.map(|i| indices[i]);
    mesh.triangle(
        if proof.reverse {
            [ids[0], ids[2], ids[1]]
        } else {
            ids
        },
        face,
    );
    Ok(())
}

#[cfg(test)]
mod tests;

impl Material {
    /// Candidate vertices are stock-relative. The original exact rational
    /// center/bottom is retained separately by stock_chart_origin(); no WCS
    /// origin is rounded or silently added to these emitted coordinates.
    pub fn tessellate_candidate(
        &self,
        options: Tessellation,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Candidate, Error> {
        material_mesh::tessellate(self, options, interrupt)
    }
}

pub(super) mod boundary;
