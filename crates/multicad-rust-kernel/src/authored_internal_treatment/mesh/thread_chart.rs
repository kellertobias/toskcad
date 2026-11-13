//! Borrowed actual source chart, shared by treated and bare thread walls.
//! No no-op Cavity is constructed and no nominal Thread proxy is consulted.
use super::*;
use crate::authored_flat_thread::Error as SourceError;
#[cfg(test)]
mod tests;
#[derive(Clone, Copy)]
pub(super) struct ThreadChart<'a> {
    pub thread: &'a AuthoredThread,
    pub entry: Option<Entry>,
    pub top: f64,
}
impl<'a> ThreadChart<'a> {
    pub fn treated(body: &'a Cavity) -> Self {
        Self {
            thread: &body.thread,
            entry: body.entry,
            top: body.tail.map_or(body.full_length_mm(), |t| t.station),
        }
    }
    fn bare(thread: &'a AuthoredThread) -> Self {
        Self {
            thread,
            entry: None,
            top: thread.full_length_mm(),
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ShellEnd {
    Start,
    End,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThreadRimEdge {
    pub vertex_indices: [usize; 2],
    pub angular_turns: [[f64; 2]; 2],
    /// Retained phase meridian rho=m*q+b, with directed source coefficients.
    pub phase_slope: [f64; 2],
    pub phase_intercept_mm: [f64; 2],
    pub point_error_mm: [f64; 2],
}
/// Immutable uncapped wall candidate bound to its actual authored source.
/// Whole-source facet receipts are retained. Closedness and literal embedding
/// are separate obligations before any manufacturing admission.
pub struct ThreadShell {
    thread: AuthoredThread,
    shell: chart::Shell,
    trace: trace::Trace,
    work_used: usize,
}
impl AuthoredThread {
    /// Actual wall facets only; the two finite rims remain open. At most16M
    /// aggregate weighted work and the existing64MiB exact ownership budget.
    pub fn tessellate_shell(
        &self,
        options: Tessellation,
        maximum_work: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<ThreadShell, SourceError> {
        shell(self, options, maximum_work, interrupt).map_err(failure)
    }
}
fn failure(e: Error) -> SourceError {
    match e {
        Error::Invalid => SourceError::Invalid,
        Error::Limit => SourceError::Limit,
        Error::Unresolved => SourceError::Unresolved,
        Error::Cancelled => SourceError::Cancelled,
    }
}
fn shell(
    thread: &AuthoredThread,
    options: Tessellation,
    maximum_work: usize,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<ThreadShell, Error> {
    options.check().map_err(|_| Error::Invalid)?;
    let mut w = Work {
        remaining: maximum_work.min(16_000_000),
        used: 0,
        interrupt,
    };
    w.step()?;
    let chart = ThreadChart::bare(thread);
    let mut grid = Grid::new();
    loop {
        let repartition = match chart::traced_chart(
            chart,
            options,
            (grid.subdivisions, grid.rows(), grid.shifted),
            &mut w,
        )? {
            Attempt::Complete((shell, trace)) => {
                w.step()?;
                return Ok(ThreadShell {
                    thread: *thread,
                    shell,
                    trace,
                    work_used: w.used,
                });
            }
            Attempt::Refine => false,
            Attempt::Repartition => true,
        };
        grid.advance(repartition)?;
    }
}
impl ThreadShell {
    pub fn source(&self) -> &AuthoredThread {
        &self.thread
    }
    pub fn mesh(&self) -> &DisplayMesh {
        &self.shell.mesh
    }
    pub fn source_chord_bound_mm(&self) -> f64 {
        self.shell.chord
    }
    pub fn source_normal_angle_bound_rad(&self) -> f64 {
        self.shell.normal
    }
    pub fn work_used(&self) -> usize {
        self.work_used
    }
    pub fn retained_bytes(&self) -> usize {
        use std::mem::size_of;
        size_of::<Self>()
            + self.trace.memory_bytes()
            + self.shell.mesh.vertices.capacity() * size_of::<f64>()
            + self.shell.mesh.indices.capacity() * size_of::<u32>()
            + self.shell.mesh.face_indices.capacity() * size_of::<u32>()
            + (self.shell.mouth.capacity() + self.shell.top.capacity()) * size_of::<chart::Edge>()
    }
    pub fn rim_count(&self, end: ShellEnd) -> usize {
        self.edges(end).len()
    }
    pub fn rim_edge(&self, end: ShellEnd, index: usize) -> Option<ThreadRimEdge> {
        let edge = self.edges(end).get(index)?;
        let pair = |v: I| [v.lo, v.hi];
        Some(ThreadRimEdge {
            vertex_indices: edge.vertices,
            angular_turns: edge.v.map(pair),
            phase_slope: pair(edge.slope),
            phase_intercept_mm: pair(edge.intercept),
            point_error_mm: edge.errors,
        })
    }
    fn edges(&self, end: ShellEnd) -> &[chart::Edge] {
        match end {
            ShellEnd::Start => &self.shell.mouth,
            ShellEnd::End => &self.shell.top,
        }
    }
}
