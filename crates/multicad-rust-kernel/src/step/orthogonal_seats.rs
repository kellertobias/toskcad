//! Retained perpendicular-seat exchange. Candidate encoding and independently
//! qualified publication are separate APIs.
use super::*;
use crate::{axial_pipe_thread::authored::StockClipPlacement, orthogonal_bore_seats::Seats};
mod qualified;
pub use qualified::{
    AppendOptions, QualificationError, QualificationOptions, QualifiedAppendReceipt,
    QualifiedOrthogonalSeatsStep, VerificationOptions, VerifiedMember, qualify,
};
#[cfg(test)]
mod tests;
#[derive(Clone, Copy)]
pub struct Options {
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            maximum_work: 16_000_000,
            maximum_bytes: 64 * 1024 * 1024,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Invalid,
    Limit,
    Cancelled,
    Encoding(StepError),
}
impl From<StepError> for Error {
    fn from(e: StepError) -> Self {
        Self::Encoding(e)
    }
}
pub struct Candidate<'a> {
    source: &'a Seats,
    placement: StockClipPlacement,
    text: String,
    faces: Vec<u32>,
    work: usize,
    bytes: usize,
}
impl<'a> Candidate<'a> {
    pub fn source(&self) -> &'a Seats {
        self.source
    }
    pub fn placement(&self) -> StockClipPlacement {
        self.placement
    }
    pub fn as_str(&self) -> &str {
        &self.text
    }
    pub fn face_entities(&self) -> &[u32] {
        &self.faces
    }
    pub fn encoding_work(&self) -> usize {
        self.work
    }
    pub fn accounted_peak_bytes(&self) -> usize {
        self.bytes
    }
}
fn check(interrupt: &dyn Fn() -> bool) -> Result<(), Error> {
    if interrupt() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
/// Original stock/cut receipts are prerequisites. This operation charges its
/// bounded encoding separately; it does not advertise material or STEP support.
pub fn candidate<'a>(
    source: &'a Seats,
    placement: StockClipPlacement,
    options: Options,
    interrupt: &dyn Fn() -> bool,
) -> Result<Candidate<'a>, Error> {
    check(interrupt)?;
    if options.maximum_work == 0
        || options.maximum_work > 16_000_000
        || options.maximum_bytes == 0
        || options.maximum_bytes > 64 * 1024 * 1024
    {
        return Err(Error::Invalid);
    }
    let boundary = source.boundary();
    let edges = boundary
        .faces
        .iter()
        .flat_map(|f| &f.rings)
        .map(Vec::len)
        .sum::<usize>();
    let work = edges * 1024 + source.bores().len() * 8192 + source.seats().len() * 4096;
    let bytes =
        source.retained_bytes() + edges * 16_384 + source.bores().len() * 32_768 + 16 * 1024;
    if work > options.maximum_work || bytes > options.maximum_bytes {
        return Err(Error::Limit);
    }
    let mut writer = Writer::new();
    let mut face_entities = Vec::new();
    writer.append_typed(|file| {
        let (solid, faces) = encode(file, source, placement, interrupt)?;
        face_entities = faces
            .into_iter()
            .map(|f| u32::try_from(f.0).map_err(|_| Error::Limit))
            .collect::<Result<_, _>>()?;
        Ok::<_, Error>(solid)
    })?;
    check(interrupt)?;
    let text = writer.finish()?;
    let work = work.checked_add(text.len()).ok_or(Error::Limit)?;
    let actual = source.retained_bytes()
        + text.capacity()
        + face_entities.capacity() * std::mem::size_of::<u32>();
    if work > options.maximum_work || actual > options.maximum_bytes {
        return Err(Error::Limit);
    }
    check(interrupt)?;
    Ok(Candidate {
        source,
        placement,
        text,
        faces: face_entities,
        work,
        bytes: bytes.max(actual),
    })
}
fn encode(
    file: &mut StepFile,
    source: &Seats,
    placement: StockClipPlacement,
    interrupt: &dyn Fn() -> bool,
) -> Result<(Entity, Vec<Entity>), Error> {
    let boundary = source.boundary();
    let matrix = placement.matrix();
    let translation = placement.translation_mm();
    let transform = |p: Point| {
        std::array::from_fn(|i| {
            matrix[i]
                .iter()
                .zip(p)
                .fold(translation[i], |s, (a, b)| s + a * b)
        })
    };
    let mut points = Vec::new();
    let mut indices = BTreeMap::new();
    for p in boundary.faces.iter().flat_map(|f| &f.rings).flatten() {
        check(interrupt)?;
        if let std::collections::btree_map::Entry::Vacant(e) = indices.entry(point_key(*p)) {
            e.insert(points.len());
            points.push(transform(*p));
        }
    }
    let mut brep = BRep::new(file, points);
    let supports = boundary
        .faces
        .iter()
        .map(|f| {
            let placed = |p| brep.vertices[indices[&point_key(p)]].point;
            let (n, u) = orthogonal_support(f, &f.rings[0], placed)?;
            let point = brep.vertices[indices[&point_key(f.rings[0][0])]].cartesian;
            Ok::<_, Error>((point, n, u))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut bounds = vec![Vec::new(); boundary.faces.len()];
    let mut walls = Vec::new();
    let mut cones = Vec::new();
    for (i, bore) in source.bores().iter().enumerate() {
        check(interrupt)?;
        let (_, normal, reference) = supports[bore.faces[1]];
        let from = bore.centers_mm[0];
        let to = bore.centers_mm[1];
        let axis = (0..3)
            .find(|&j| from[j] != to[j])
            .expect("qualified original bore");
        let sign = (to[axis] - from[axis]).signum();
        let height = (to[axis] - from[axis]).abs();
        let mut layers = vec![(0., from, bore.radius_mm), (height, to, bore.radius_mm)];
        for end in 0..2 {
            if let Some(seat) = source
                .seats()
                .iter()
                .find(|s| s.bore == i && usize::from(s.end) == end)
            {
                let at = end;
                layers[at].2 = seat.radius_mm;
                let mut p = bore.centers_mm[end];
                p[axis] += if end == 0 {
                    sign * seat.depth_mm
                } else {
                    -sign * seat.depth_mm
                };
                layers.push((
                    if end == 0 {
                        seat.depth_mm
                    } else {
                        height - seat.depth_mm
                    },
                    p,
                    bore.radius_mm,
                ));
            }
        }
        layers.sort_by(|a, b| a.0.total_cmp(&b.0));
        let rims = layers
            .iter()
            .map(|(h, p, r)| {
                CircularRim::new(brep.file, transform(*p), normal, reference, *r)
                    .map(|rim| (*h, rim))
            })
            .collect::<Result<Vec<_>, _>>()?;
        bounds[bore.faces[0]].push(circular_bound(brep.file, rims[0].1.edge, true, false));
        bounds[bore.faces[1]].push(circular_bound(
            brep.file,
            rims.last().expect("finite end rim").1.edge,
            false,
            false,
        ));
        for pair in rims.windows(2) {
            let face = round_wall(
                brep.file,
                &pair[0].1,
                &pair[1].1,
                pair[1].0 - pair[0].0,
                normal,
                reference,
                false,
            )?;
            if pair[0].1.radius == pair[1].1.radius {
                walls.push(face);
            } else {
                cones.push(face);
            }
        }
    }
    let mut faces = Vec::new();
    for (i, face) in boundary.faces.iter().enumerate() {
        check(interrupt)?;
        let (point, normal, reference) = supports[i];
        let mut rings = face
            .rings
            .iter()
            .enumerate()
            .map(|(i, ring)| {
                brep.bound(
                    &ring
                        .iter()
                        .map(|p| indices[&point_key(*p)])
                        .collect::<Vec<_>>(),
                    i == 0,
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        rings.append(&mut bounds[i]);
        faces.push(brep.file.face(&rings, point, normal, reference));
    }
    faces.extend(walls);
    faces.extend(cones);
    let solid = brep.file.solid(&faces);
    Ok((solid, faces))
}
