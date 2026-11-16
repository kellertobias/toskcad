//! Fixed-frame finite head over an independently requalified actual D11 thread.
//! The auxiliary original cap remains a verification object, never a physical
//! duplicate of the annulus. This stage admits only the explicitly sewn
//! reference and does not decode arbitrary head files or bind an authored law.
use super::*;
use crate::tapered_thread::interval::I;
pub(crate) mod encoding;
mod geometry;
#[cfg(test)]
mod tests;
mod topology;
const MAX_WORK: usize = 96_000_000;
const MAX_BYTES: usize = 64 * 1024 * 1024;
const SCRATCH: usize = 64 * 1024;
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Plane {
    pub origin: [f64; 3],
    pub normal: [f64; 3],
    pub reference: [f64; 3],
}
#[derive(Clone, Debug, PartialEq)]
pub struct HeadFace {
    pub plane: Plane,
    /// Indices of actual straight head edges and their coedge direction.
    pub outer: Vec<(usize, bool)>,
    /// Original EDGE_CURVE entities on the SAME auxiliary thread cap. Only
    /// the annulus has this loop; original IDs are never coordinate-merged.
    pub inner: Vec<(u32, bool)>,
    pub forward: bool,
}
#[derive(Clone, Debug, PartialEq)]
pub struct HeadInput {
    pub vertices: Vec<[f64; 3]>,
    pub edges: Vec<[usize; 2]>,
    pub faces: Vec<HeadFace>,
    pub auxiliary_cap_entity: u32,
}
#[derive(Clone, Copy, Debug)]
pub struct HeadOptions {
    pub maximum_work: usize,
    pub maximum_bytes: usize,
    /// Other caller state beyond the borrowed neutral document, actual Parsed
    /// query and moved HeadInput, all counted separately by this stage.
    pub additional_resident_bytes: usize,
}
impl Default for HeadOptions {
    fn default() -> Self {
        Self {
            maximum_work: MAX_WORK,
            maximum_bytes: MAX_BYTES,
            additional_resident_bytes: 0,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Receipt {
    pub inherited_document_work: usize,
    /// Includes the Parsed prefix once, projection and actual occupied stage.
    pub material_work_with_parsed_prefix: usize,
    pub new_head_work: usize,
    pub aggregate_work: usize,
    pub borrowed_document_bytes: usize,
    pub head_input_bytes: usize,
    pub accounted_peak_bytes: usize,
}
pub struct OccupiedClippedHead<'parsed, 'document, 'input> {
    thread: Occupied<'parsed, 'document, 'input>,
    head: HeadInput,
    polygon: [[f64; 2]; 6],
    stations: [f64; 2],
    planar_sewing: f64,
    join_sewing: f64,
    minimum_wall: f64,
    receipt: Receipt,
}
impl<'p, 'd, 'i> OccupiedClippedHead<'p, 'd, 'i> {
    pub fn capped_thread(&self) -> &Occupied<'p, 'd, 'i> {
        &self.thread
    }
    pub fn head(&self) -> &HeadInput {
        &self.head
    }
    pub fn head_polygon_mm(&self) -> &[[f64; 2]; 6] {
        &self.polygon
    }
    pub fn head_stations_mm(&self) -> [f64; 2] {
        self.stations
    }
    pub fn planar_sewing_mm(&self) -> f64 {
        self.planar_sewing
    }
    pub fn join_sewing_mm(&self) -> f64 {
        self.join_sewing
    }
    /// Conservative clearance of the complete actual control hull from the
    /// hex halfspaces; this is a lower bound, not an exact surface distance.
    pub fn minimum_wall_mm(&self) -> f64 {
        self.minimum_wall
    }
    pub fn receipt(&self) -> Receipt {
        self.receipt
    }
}
struct Work<'a> {
    used: usize,
    maximum: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Work<'_> {
    fn charge(&mut self, n: usize) -> Result<()> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        self.used = self
            .used
            .checked_add(n)
            .filter(|n| *n <= self.maximum)
            .ok_or(Error::Limit)?;
        Ok(())
    }
}
fn storage(head: &HeadInput) -> Result<usize> {
    let mut bytes = std::mem::size_of::<HeadInput>();
    for n in [
        head.vertices
            .capacity()
            .checked_mul(std::mem::size_of::<[f64; 3]>()),
        head.edges
            .capacity()
            .checked_mul(std::mem::size_of::<[usize; 2]>()),
        head.faces
            .capacity()
            .checked_mul(std::mem::size_of::<HeadFace>()),
    ] {
        bytes = bytes
            .checked_add(n.ok_or(Error::Limit)?)
            .ok_or(Error::Limit)?;
    }
    for face in &head.faces {
        bytes = bytes
            .checked_add(
                face.outer
                    .capacity()
                    .checked_mul(std::mem::size_of::<(usize, bool)>())
                    .ok_or(Error::Limit)?,
            )
            .and_then(|n| {
                n.checked_add(
                    face.inner
                        .capacity()
                        .checked_mul(std::mem::size_of::<(u32, bool)>())?,
                )
            })
            .ok_or(Error::Limit)?;
    }
    Ok(bytes)
}
fn bounded_peak(bytes: usize, maximum: usize) -> Result<usize> {
    if bytes > maximum {
        Err(Error::Limit)
    } else {
        Ok(bytes)
    }
}
fn interval([lo, hi]: [f64; 2]) -> I {
    I::new(lo, hi)
}
fn require(v: bool) -> Result<()> {
    if v { Ok(()) } else { Err(Error::Unsupported) }
}

/// Requalify the SAME immutable actual capped thread before occupied contact.
/// The physical extension removes precisely its upper auxiliary cap and puts
/// the original rim into the annulus, with opposite uses to the retained sides.
pub fn certify<'p, 'd, 'i>(
    parsed: &'p Parsed<'d, 'i>,
    head: HeadInput,
    tolerance: f64,
    minimum: f64,
    options: HeadOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<OccupiedClippedHead<'p, 'd, 'i>> {
    if !tolerance.is_finite()
        || tolerance <= 0.
        || tolerance > 1e-4
        || !minimum.is_finite()
        || minimum <= 0.
        || options.maximum_work == 0
        || options.maximum_work > MAX_WORK
        || options.maximum_bytes == 0
        || options.maximum_bytes > MAX_BYTES
    {
        return Err(Error::Invalid);
    }
    if head.vertices.len() != 12
        || head.edges.len() != 18
        || head.faces.len() != 8
        || head
            .faces
            .iter()
            .any(|f| f.outer.len() > 6 || f.inner.len() > 4)
    {
        return Err(Error::Unsupported);
    }
    let document_work = parsed.document().work_used();
    let document_bytes = parsed.document().retained_budget_bytes();
    let head_bytes = storage(&head)?;
    let other = document_bytes
        .checked_add(head_bytes)
        .and_then(|n| n.checked_add(options.additional_resident_bytes))
        .ok_or(Error::Limit)?;
    let local_peak = bounded_peak(
        other
            .checked_add(parsed.query_budget_bytes())
            .and_then(|n| n.checked_add(SCRATCH))
            .ok_or(Error::Limit)?,
        options.maximum_bytes,
    )?;
    let mut work = Work {
        used: 0,
        maximum: options
            .maximum_work
            .checked_sub(document_work)
            .and_then(|n| n.checked_sub(parsed.work_used()))
            .ok_or(Error::Limit)?,
        interrupt,
    };
    work.charge(1)?;
    let cap = topology::qualify(parsed, &head, &mut work)?;
    let child_work = options
        .maximum_work
        .checked_sub(document_work)
        .and_then(|n| n.checked_sub(work.used))
        .ok_or(Error::Limit)?
        .min(64_000_000);
    let child_bytes = options
        .maximum_bytes
        .checked_sub(other)
        .ok_or(Error::Limit)?;
    let material_options = MaterialOptions {
        maximum_work: child_work,
        maximum_bytes: child_bytes,
        ..Default::default()
    };
    let thread = parsed
        .projection(material_options, interrupt)?
        .occupied(material_options, interrupt)?;
    if thread.projection().cap_face_indices()[1] != cap {
        return Err(Error::Unsupported);
    }
    // The child receipt includes Parsed once. Debit only its fresh work into
    // the same parent meter, retaining the common work/cancellation guard.
    let fresh_child_work = thread
        .work_used()
        .checked_sub(parsed.work_used())
        .ok_or(Error::Limit)?;
    work.charge(fresh_child_work)?;
    let geometry = geometry::qualify(parsed, &thread, &head, tolerance, minimum, &mut work)?;
    work.charge(1)?;
    // Projection/Occupied retain only the immutable Parsed borrow and inline
    // scalars after their scratch drops. Their query byte receipt is a PAST
    // phase peak, not retained heap. The 64 KiB local reservation covers that
    // inline owner plus topology and the post-child geometry phase. If either
    // typed child gains heap fields, this phase-max argument must be updated.
    let peak = bounded_peak(
        local_peak.max(
            other
                .checked_add(thread.query_budget_bytes())
                .ok_or(Error::Limit)?,
        ),
        options.maximum_bytes,
    )?;
    let new_head_work = work
        .used
        .checked_sub(fresh_child_work)
        .ok_or(Error::Limit)?;
    let receipt = Receipt {
        inherited_document_work: document_work,
        material_work_with_parsed_prefix: thread.work_used(),
        new_head_work,
        aggregate_work: document_work
            .checked_add(parsed.work_used())
            .and_then(|n| n.checked_add(work.used))
            .ok_or(Error::Limit)?,
        borrowed_document_bytes: document_bytes,
        head_input_bytes: head_bytes,
        accounted_peak_bytes: peak,
    };
    Ok(OccupiedClippedHead {
        thread,
        head,
        polygon: geometry.polygon,
        stations: geometry.stations,
        planar_sewing: geometry.planar_sewing,
        join_sewing: geometry.join_sewing,
        minimum_wall: geometry.minimum_wall,
        receipt,
    })
}

impl HeadInput {
    /// Construct finite proposal values only. Occupied admission is exclusively
    /// `certify`; the actual scalar caches and original rim IDs are retained.
    pub fn regular_hex(
        parsed: &Parsed<'_, '_>,
        auxiliary_cap_entity: u32,
        across: f64,
        height: f64,
    ) -> Result<Self> {
        if !across.is_finite()
            || across <= 0.
            || !height.is_finite()
            || height <= 0.
            || across > 1e7
            || height > 1e7
        {
            return Err(Error::Invalid);
        }
        let cap = parsed
            .faces
            .iter()
            .find(|f| f.id == auxiliary_cap_entity)
            .ok_or(Error::Invalid)?;
        let Support::Plane { height: station } = cap.support else {
            return Err(Error::Unsupported);
        };
        let low = parsed.points[station].coordinates[2].value();
        let high = low + height;
        if !high.is_finite() || high.abs() > 1e7 {
            return Err(Error::Limit);
        }
        let r = across / 3_f64.sqrt();
        let h = across / 2.;
        let xy = [
            [r, 0.],
            [r / 2., h],
            [-r / 2., h],
            [-r, 0.],
            [-r / 2., -h],
            [r / 2., -h],
        ];
        let vertices = [low, high]
            .into_iter()
            .flat_map(|z| xy.map(|p| [p[0], p[1], z]))
            .collect::<Vec<_>>();
        let mut edges = Vec::with_capacity(18);
        for start in [0, 6] {
            for i in 0..6 {
                edges.push([start + i, start + (i + 1) % 6]);
            }
        }
        for i in 0..6 {
            edges.push([i, i + 6]);
        }
        let axis = |z| Plane {
            origin: [0., 0., z],
            normal: [0., 0., 1.],
            reference: [1., 0., 0.],
        };
        let mut faces = Vec::with_capacity(8);
        faces.push(HeadFace {
            plane: axis(low),
            outer: (0..6).rev().map(|i| (i, false)).collect(),
            inner: cap
                .coedges
                .iter()
                .map(|c| (parsed.edges[c.edge].id, c.forward))
                .collect(),
            forward: false,
        });
        faces.push(HeadFace {
            plane: axis(high),
            outer: (0..6).map(|i| (i + 6, true)).collect(),
            inner: vec![],
            forward: true,
        });
        for i in 0..6 {
            let j = (i + 1) % 6;
            let dx = xy[j][0] - xy[i][0];
            let dy = xy[j][1] - xy[i][1];
            faces.push(HeadFace {
                plane: Plane {
                    origin: vertices[i],
                    normal: [dy, -dx, 0.],
                    reference: [dx, dy, 0.],
                },
                outer: vec![(i, true), (12 + j, true), (6 + i, false), (12 + i, false)],
                inner: vec![],
                forward: true,
            });
        }
        Ok(Self {
            vertices,
            edges,
            faces,
            auxiliary_cap_entity,
        })
    }
}
