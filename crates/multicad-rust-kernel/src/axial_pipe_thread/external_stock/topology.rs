//! Physical ownership: the crest needs no zero-width shoulder. Auxiliary
//! shoulder radial seams and tail axial seams are absent from physical faces.
use super::*;
use std::collections::{BTreeMap, BTreeSet};
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Face {
    Thread(StockClipFace),
    Shoulder,
    TailCylinder,
    TailCap,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Vertex {
    Thread(StockClipVertex),
    ShoulderRoot(StockClipRail),
    TailTop(StockClipRail),
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Curve {
    Thread(StockClipCurve),
    ShoulderOuter(StockClipBand),
    TailTop(StockClipBand),
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Coedge {
    pub curve: Curve,
    pub forward: bool,
}
#[derive(Clone, Debug)]
pub struct FaceLoop {
    face: Face,
    loops: Vec<Vec<Coedge>>,
}
impl FaceLoop {
    pub fn face(&self) -> Face {
        self.face
    }
    pub fn loops(&self) -> &[Vec<Coedge>] {
        &self.loops
    }
}
#[derive(Clone, Debug)]
pub struct Boundary {
    source: ExternalStockCut,
    faces: Vec<FaceLoop>,
    curves: Vec<Curve>,
    vertices: Vec<Vertex>,
    work: usize,
    bytes: usize,
}
const BANDS: [StockClipBand; 4] = [
    StockClipBand::Lower,
    StockClipBand::Crest,
    StockClipBand::Upper,
    StockClipBand::Root,
];
const RAILS: [StockClipRail; 4] = [
    StockClipRail::RootLower,
    StockClipRail::LowerCrest,
    StockClipRail::CrestUpper,
    StockClipRail::UpperRoot,
];
fn end(rail: StockClipRail) -> Vertex {
    Vertex::Thread(StockClipVertex {
        station: StockClipStation::End,
        rail,
    })
}
fn outer(rail: StockClipRail) -> Vertex {
    if matches!(rail, StockClipRail::LowerCrest | StockClipRail::CrestUpper) {
        end(rail)
    } else {
        Vertex::ShoulderRoot(rail)
    }
}
fn rails(band: StockClipBand) -> [StockClipRail; 2] {
    let i = match band {
        StockClipBand::Lower => 0,
        StockClipBand::Crest => 1,
        StockClipBand::Upper => 2,
        StockClipBand::Root => 3,
    };
    [RAILS[i], RAILS[(i + 1) % 4]]
}
fn bottom(band: StockClipBand) -> Curve {
    if band == StockClipBand::Crest {
        Curve::Thread(StockClipCurve::Ring {
            station: StockClipStation::End,
            band,
        })
    } else {
        Curve::ShoulderOuter(band)
    }
}
impl Curve {
    pub fn endpoints(self) -> [Vertex; 2] {
        match self {
            Self::Thread(c) => c.endpoints().map(Vertex::Thread),
            Self::ShoulderOuter(b) => rails(b).map(outer),
            Self::TailTop(b) => rails(b).map(Vertex::TailTop),
        }
    }
}
fn c(curve: Curve, forward: bool) -> Coedge {
    Coedge { curve, forward }
}
impl ExternalStockCut {
    pub fn boundary(
        &self,
        options: StockClipOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Boundary, Error> {
        let mut w = Work::new(options, size_of::<Self>() + 128 * 1024, interrupt)?;
        w.charge(self.work)?;
        let mut faces = Vec::with_capacity(9);
        for face in self.thread.physical_faces() {
            w.charge(1)?;
            if face == StockClipFace::EndCap {
                continue;
            }
            faces.push(FaceLoop {
                face: Face::Thread(face),
                loops: vec![
                    self.thread
                        .face_loop(face)
                        .coedges()
                        .iter()
                        .map(|v| c(Curve::Thread(v.curve), v.forward))
                        .collect(),
                ],
            });
        }
        let inner = |band| {
            Curve::Thread(StockClipCurve::Ring {
                station: StockClipStation::End,
                band,
            })
        };
        faces.push(FaceLoop {
            face: Face::Shoulder,
            loops: vec![vec![
                c(inner(StockClipBand::Lower), false),
                c(inner(StockClipBand::Root), false),
                c(inner(StockClipBand::Upper), false),
                c(bottom(StockClipBand::Upper), true),
                c(bottom(StockClipBand::Root), true),
                c(bottom(StockClipBand::Lower), true),
            ]],
        });
        faces.push(FaceLoop {
            face: Face::TailCylinder,
            loops: vec![
                BANDS
                    .into_iter()
                    .rev()
                    .map(|b| c(bottom(b), false))
                    .collect(),
                BANDS
                    .into_iter()
                    .map(|b| c(Curve::TailTop(b), true))
                    .collect(),
            ],
        });
        faces.push(FaceLoop {
            face: Face::TailCap,
            loops: vec![
                BANDS
                    .into_iter()
                    .rev()
                    .map(|b| c(Curve::TailTop(b), false))
                    .collect(),
            ],
        });
        let mut edges = BTreeMap::<Curve, Vec<(usize, bool)>>::new();
        let mut vertices = BTreeSet::new();
        let mut fans = BTreeMap::<Vertex, BTreeMap<Curve, Curve>>::new();
        for (face_id, face) in faces.iter().enumerate() {
            for ring in &face.loops {
                if ring.len() < 3 {
                    return Err(Error::Unresolved);
                }
                for (i, use_) in ring.iter().enumerate() {
                    w.charge(128)?;
                    let mut ends = use_.curve.endpoints();
                    if !use_.forward {
                        ends.reverse()
                    }
                    let next = &ring[(i + 1) % ring.len()];
                    let mut next_ends = next.curve.endpoints();
                    if !next.forward {
                        next_ends.reverse()
                    }
                    if ends[0] == ends[1] || ends[1] != next_ends[0] {
                        return Err(Error::Unresolved);
                    }
                    vertices.extend(ends);
                    edges
                        .entry(use_.curve)
                        .or_default()
                        .push((face_id, use_.forward));
                    if fans
                        .entry(ends[1])
                        .or_default()
                        .insert(use_.curve, next.curve)
                        .is_some()
                    {
                        return Err(Error::Unresolved);
                    }
                }
            }
        }
        for uses in edges.values() {
            w.charge(16)?;
            if uses.len() != 2 || uses[0].1 == uses[1].1 {
                return Err(Error::Unresolved);
            }
        }
        for fan in fans.values() {
            let first = *fan.first_key_value().ok_or(Error::Unresolved)?.0;
            let mut current = first;
            let mut count = 0;
            loop {
                w.charge(32)?;
                current = *fan.get(&current).ok_or(Error::Unresolved)?;
                count += 1;
                if current == first {
                    break;
                }
                if count > fan.len() {
                    return Err(Error::Unresolved);
                }
            }
            if count != fan.len() || count < 2 {
                return Err(Error::Unresolved);
            }
        }
        let mut reached = BTreeSet::from([0]);
        loop {
            let n = reached.len();
            for uses in edges.values() {
                w.charge(16)?;
                if reached.contains(&uses[0].0) || reached.contains(&uses[1].0) {
                    reached.insert(uses[0].0);
                    reached.insert(uses[1].0);
                }
            }
            if reached.len() == n {
                break;
            }
        }
        let chi = faces
            .iter()
            .map(|f| 2 - f.loops.len() as isize)
            .sum::<isize>();
        if reached.len() != faces.len() || vertices.len() as isize - edges.len() as isize + chi != 2
        {
            return Err(Error::Unresolved);
        }
        let curves = edges.into_keys().collect::<Vec<_>>();
        let vertices = vertices.into_iter().collect::<Vec<_>>();
        let bytes = size_of::<Boundary>()
            + faces.capacity() * size_of::<FaceLoop>()
            + curves.capacity() * size_of::<Curve>()
            + vertices.capacity() * size_of::<Vertex>()
            + faces
                .iter()
                .map(|f| {
                    f.loops.capacity() * size_of::<Vec<Coedge>>()
                        + f.loops
                            .iter()
                            .map(|r| r.capacity() * size_of::<Coedge>())
                            .sum::<usize>()
                })
                .sum::<usize>();
        if bytes + 128 * 1024 > w.limit {
            return Err(Error::Limit);
        }
        w.charge(1)?;
        Ok(Boundary {
            source: *self,
            faces,
            curves,
            vertices,
            work: w.used,
            bytes,
        })
    }
}
impl Boundary {
    pub fn source(&self) -> &ExternalStockCut {
        &self.source
    }
    pub fn faces(&self) -> &[FaceLoop] {
        &self.faces
    }
    pub fn curves(&self) -> &[Curve] {
        &self.curves
    }
    pub fn vertices(&self) -> &[Vertex] {
        &self.vertices
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
}
