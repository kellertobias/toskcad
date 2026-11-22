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
    pub exact_top: Option<crate::Rational>,
}
impl<'a> ThreadChart<'a> {
    pub fn treated(body: &'a Cavity) -> Self {
        Self {
            thread: &body.thread,
            entry: body.entry,
            top: body.tail.map_or(body.full_length_mm(), |t| t.station),
            exact_top: if body.tail.is_none() {
                body.exact_receiving_extent()
                    .or(body.thread.exact_axial_extent())
            } else {
                None
            },
        }
    }
    fn bare(thread: &'a AuthoredThread) -> Self {
        Self {
            thread,
            entry: None,
            top: thread.full_length_mm(),
            exact_top: thread.exact_axial_extent(),
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
    /// Unsigned source v; physical angle is hand*2pi*v. A computational
    /// sampling cut may start above0 and end above1. These bounds do not prove
    /// equality or order of footpoints; the paired mesher uses exact tokens.
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
    end_station_mm: f64,
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
        shell(
            self,
            self.full_length_mm(),
            self.exact_axial_extent(),
            false,
            options,
            maximum_work,
            interrupt,
        )
        .map_err(failure)
    }
}
impl AuthoredThread {
    /// Finite prefix of this same immutable authored source. The original
    /// specification, phase coefficients and full-length metadata are retained;
    /// the returned open wall ends at the explicit binary64 station.
    pub fn tessellate_prefix_shell(
        &self,
        end_station_mm: f64,
        options: Tessellation,
        maximum_work: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<ThreadShell, SourceError> {
        if !end_station_mm.is_finite()
            || end_station_mm < self.spec().pitch_mm
            || !self.contains_axial_parameter(end_station_mm)
        {
            return Err(SourceError::Invalid);
        }
        shell(
            self,
            end_station_mm,
            None,
            false,
            options,
            maximum_work,
            interrupt,
        )
        .map_err(failure)
    }
}
impl AuthoredThread {
    /// Alternate computational angular seam. Original phase, coefficients and
    /// source endpoints remain unchanged; every finite field proof runs again.
    pub fn tessellate_prefix_shell_shifted(
        &self,
        end_station_mm: f64,
        options: Tessellation,
        maximum_work: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<ThreadShell, SourceError> {
        if !end_station_mm.is_finite()
            || end_station_mm < self.spec().pitch_mm
            || !self.contains_axial_parameter(end_station_mm)
        {
            return Err(SourceError::Invalid);
        }
        shell(
            self,
            end_station_mm,
            None,
            true,
            options,
            maximum_work,
            interrupt,
        )
        .map_err(failure)
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
    end_station_mm: f64,
    exact_top: Option<crate::Rational>,
    shifted: bool,
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
    let chart = ThreadChart {
        top: end_station_mm,
        exact_top,
        ..ThreadChart::bare(thread)
    };
    let mut grid = Grid::new();
    grid.shifted = shifted;
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
                    end_station_mm,
                });
            }
            Attempt::Refine => false,
            Attempt::Repartition => true,
        };
        grid.advance(repartition)?;
    }
}
impl ThreadShell {
    pub fn end_station_mm(&self) -> f64 {
        self.end_station_mm
    }
    pub(crate) fn exact_rim_turns(
        &self,
        vertex: usize,
        maximum_work: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<(num_rational::BigRational, usize), SourceError> {
        self.parameter(vertex, maximum_work, interrupt)
            .map_err(failure)
    }
    fn parameter(
        &self,
        vertex: usize,
        maximum_work: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<(num_rational::BigRational, usize), Error> {
        use super::exact::Q;
        let mut work = Work {
            remaining: maximum_work.min(16_000_000),
            used: 0,
            interrupt,
        };
        work.step()?;
        let original = self.trace.parameter(vertex, &mut work)?;
        let zero = Q::integer(0, &mut work)?;
        let one = Q::integer(1, &mut work)?;
        // Trace is owned by the sealed grid: v in[cut,cut+1), cut being0 or
        // the fixed positive dyadic seam shift<1. No arbitrary floor or float
        // comparison decides the physical angle; trueperiodicity is exact.
        let mut physical = if original.cmp(&one, &mut work)? == std::cmp::Ordering::Less {
            original
        } else {
            original.sub(&one, &mut work)?
        };
        if self.thread.spec().left_hand
            && physical.cmp(&zero, &mut work)? != std::cmp::Ordering::Equal
        {
            physical = one.sub(&physical, &mut work)?;
        }
        let value = physical.rational(&mut work)?;
        work.step()?;
        Ok((value, work.used))
    }
    pub fn into_mesh(self) -> DisplayMesh {
        self.shell.mesh
    }
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

#[cfg(test)]
mod consumption_tests {
    use super::*;
    #[test]
    fn actual_bare_flat_source_wall_can_be_consumed_without_a_nominal_rebuild() {
        let source = crate::authored_hanger_bolt::tests::fixture(false, false).wood;
        let shell = source
            .tessellate_shell(
                Tessellation {
                    deflection: 0.1,
                    angle: 0.3,
                    triangle_limit: 1_000_000,
                },
                16_000_000,
                None,
            )
            .unwrap();
        let count = shell.mesh().face_indices.len();
        assert_eq!(shell.into_mesh().face_indices.len(), count);
    }
}
