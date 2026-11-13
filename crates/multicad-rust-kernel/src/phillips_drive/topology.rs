//! Closed expression-owned source topology. Its incidence certificate is not
//! independent evidence about literal emitted STEP data or manufacturing mesh.
use super::*;
mod analytic;
mod geometry;
#[cfg(test)]
mod tests;
pub(in crate::phillips_drive) mod validate;
pub use analytic::AnalyticCurve;
pub use geometry::{LineSupport, Surface};
pub use validate::TopologyCertificate;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VertexKey {
    Mouth {
        quarter: u8,
        role: VertexRole,
    },
    Event {
        role: EventRole,
        quarter: u8,
        lower: bool,
    },
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MeridianStation {
    Entry,
    Floor,
}
impl MeridianStation {
    fn event(self) -> EventRole {
        match self {
            Self::Entry => EventRole::EntryJoin,
            Self::Floor => EventRole::FloorJoin,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CurveKey {
    Mouth {
        edge: usize,
    },
    Trim {
        envelope: EnvelopeRole,
        wall: WallKind,
        quarter: u8,
        lower: bool,
    },
    MeridianCircle {
        station: MeridianStation,
        quarter: u8,
    },
    WallJoin {
        quarter: u8,
        lower: bool,
    },
    Center {
        quarter: u8,
    },
    FloorGenerator {
        quarter: u8,
    },
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FaceKey {
    Mouth,
    Envelope {
        role: EnvelopeRole,
        quarter: u8,
    },
    Wall {
        kind: WallKind,
        quarter: u8,
        lower: bool,
    },
}
#[derive(Clone, Copy, Debug)]
pub struct BoundaryVertex {
    key: VertexKey,
    point: [I; 3],
}
impl BoundaryVertex {
    pub fn source_identity(&self) -> VertexKey {
        self.key
    }
    pub fn point_mm_bounds(&self) -> [Bounds; 3] {
        self.point.map(bounds)
    }
}
#[derive(Clone, Copy, Debug)]
pub struct BoundaryEdge {
    ends: [usize; 2],
    key: CurveKey,
}
impl BoundaryEdge {
    pub fn vertices(&self) -> [usize; 2] {
        self.ends
    }
    pub fn source_identity(&self) -> CurveKey {
        self.key
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EdgeUse {
    pub edge: usize,
    pub forward: bool,
}
#[derive(Clone, Debug)]
pub struct BoundaryFace {
    key: FaceKey,
    uses: Vec<EdgeUse>,
}
impl BoundaryFace {
    pub fn source_identity(&self) -> FaceKey {
        self.key
    }
    pub fn edge_uses(&self) -> &[EdgeUse] {
        &self.uses
    }
}
#[derive(Clone, Debug)]
pub struct BoundaryTopology<'a> {
    events: BoundaryEvents<'a>,
    vertices: Vec<BoundaryVertex>,
    edges: Vec<BoundaryEdge>,
    faces: Vec<BoundaryFace>,
    certificate: TopologyCertificate,
}
impl<'a> BoundaryTopology<'a> {
    pub fn events(&self) -> &BoundaryEvents<'a> {
        &self.events
    }
    pub fn vertices(&self) -> &[BoundaryVertex] {
        &self.vertices
    }
    pub fn edges(&self) -> &[BoundaryEdge] {
        &self.edges
    }
    pub fn faces(&self) -> &[BoundaryFace] {
        &self.faces
    }
    pub fn edge_index(&self, key: CurveKey) -> Option<usize> {
        self.edges.iter().position(|e| e.key == key)
    }
    pub fn certificate(&self) -> TopologyCertificate {
        self.certificate
    }
}
// Resolve the finite source schema once. Later assembly uses typed indices,
// so no edge can silently fall back to a different named mouth vertex.
struct MouthNodes {
    ids: [[usize; 5]; 4],
}
impl MouthNodes {
    fn read(vertices: &[SectionVertex]) -> Result<Self, Error> {
        let roles = [
            VertexRole::LowerEnvelope,
            VertexRole::LowerJoin,
            VertexRole::InnerTip,
            VertexRole::UpperJoin,
            VertexRole::UpperEnvelope,
        ];
        let mut ids = [[0; 5]; 4];
        for quarter in 0..4 {
            for (index, role) in roles.into_iter().enumerate() {
                let found: Vec<_> = vertices
                    .iter()
                    .enumerate()
                    .filter(|(_, v)| v.source_identity() == (quarter, role))
                    .map(|(i, _)| i)
                    .collect();
                if found.len() != 1 {
                    return Err(Error::Unresolved);
                }
                ids[quarter as usize][index] = found[0];
            }
        }
        Ok(Self { ids })
    }
    fn join(&self, quarter: u8, lower: bool) -> usize {
        self.ids[quarter as usize][if lower { 1 } else { 3 }]
    }
    fn tip(&self, quarter: u8) -> usize {
        self.ids[quarter as usize][2]
    }
}
// Publication consumes an assembled source graph only after the independent
// raw incidence checker succeeds. Failed assemblies cannot escape as bodies.
struct Assembly<'a> {
    events: BoundaryEvents<'a>,
    vertices: Vec<BoundaryVertex>,
    edges: Vec<BoundaryEdge>,
    faces: Vec<BoundaryFace>,
}
impl<'a> Assembly<'a> {
    fn finish(self) -> Result<BoundaryTopology<'a>, Error> {
        let Self {
            events,
            vertices,
            edges,
            faces,
        } = self;
        validate::certify(
            &vertices.iter().map(|v| v.point).collect::<Vec<_>>(),
            &edges.iter().map(|e| e.ends).collect::<Vec<_>>(),
            &faces.iter().map(|f| f.uses.clone()).collect::<Vec<_>>(),
        )
        .map(|certificate| BoundaryTopology {
            events,
            vertices,
            edges,
            faces,
            certificate,
        })
    }
}
fn i([lo, hi]: Bounds) -> I {
    I::new(lo, hi)
}
fn checkpoint(cancel: &dyn Fn() -> bool) -> Result<(), Error> {
    if cancel() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
pub(in crate::phillips_drive) fn event_key(role: EventRole, quarter: u8, lower: bool) -> VertexKey {
    match role {
        EventRole::MouthEnvelope => VertexKey::Mouth {
            quarter,
            role: if lower {
                VertexRole::LowerEnvelope
            } else {
                VertexRole::UpperEnvelope
            },
        },
        EventRole::FloorApex => VertexKey::Event {
            role,
            quarter: 0,
            lower: false,
        },
        EventRole::WallContact => VertexKey::Event {
            role,
            quarter,
            lower: false,
        },
        _ => VertexKey::Event {
            role,
            quarter,
            lower,
        },
    }
}
fn vertex_id(vertices: &[BoundaryVertex], key: VertexKey) -> Result<usize, Error> {
    vertices
        .iter()
        .position(|v| v.key == key)
        .ok_or(Error::Unresolved)
}
fn edge_id(edges: &[BoundaryEdge], key: CurveKey) -> Result<usize, Error> {
    edges
        .iter()
        .position(|e| e.key == key)
        .ok_or(Error::Unresolved)
}
fn mouth_edge(events: &BoundaryEvents<'_>, quarter: u8, role: EdgeRole) -> Result<CurveKey, Error> {
    events
        .mouth()
        .edges()
        .iter()
        .position(|e| e.source_identity() == (quarter, role))
        .map(|edge| CurveKey::Mouth { edge })
        .ok_or(Error::Unresolved)
}
fn uses(edges: &[BoundaryEdge], keys: Vec<(CurveKey, bool)>) -> Result<Vec<EdgeUse>, Error> {
    keys.into_iter()
        .map(|(key, forward)| edge_id(edges, key).map(|edge| EdgeUse { edge, forward }))
        .collect()
}
fn trim_uses(
    events: &BoundaryEvents<'_>,
    quarter: u8,
    lower: bool,
    envelope: Option<EnvelopeRole>,
    wall: Option<WallKind>,
    forward: bool,
) -> Vec<(CurveKey, bool)> {
    let mut indices: Vec<_> = events
        .trims()
        .iter()
        .enumerate()
        .filter(|(_, t)| {
            envelope.is_none_or(|e| t.envelope_role() == e)
                && wall.is_none_or(|w| t.wall_kind() == w)
        })
        .map(|(index, _)| index)
        .collect();
    if !forward {
        indices.reverse()
    };
    indices
        .into_iter()
        .map(|index| {
            (
                CurveKey::Trim {
                    envelope: events.trims()[index].envelope_role(),
                    wall: events.trims()[index].wall_kind(),
                    quarter,
                    lower,
                },
                forward,
            )
        })
        .collect()
}
pub(super) fn construct<'a>(
    source: &'a Cutter,
    cancel: &dyn Fn() -> bool,
) -> Result<BoundaryTopology<'a>, Error> {
    let events = source.boundary_events_with_cancel(cancel)?;
    checkpoint(cancel)?;
    let mouth_nodes = MouthNodes::read(events.mouth().vertices())?;
    let mut vertices = Vec::new();
    let mut edges = Vec::new();
    let mut faces = Vec::new();
    for v in events.mouth().vertices() {
        let (quarter, role) = v.source_identity();
        let [x, y] = v.point_mm_bounds().map(i);
        vertices.push(BoundaryVertex {
            key: VertexKey::Mouth { quarter, role },
            point: [x, y, I::point(0.)],
        })
    }
    for role in [
        EventRole::EntryJoin,
        EventRole::FloorJoin,
        EventRole::WallJoin,
        EventRole::WallContact,
        EventRole::FloorApex,
    ] {
        checkpoint(cancel)?;
        for quarter in 0..4 {
            for lower in [false, true] {
                let key = event_key(role, quarter, lower);
                if vertices.iter().any(|v| v.key == key) {
                    continue;
                }
                let point = events
                    .event_orbit_point_mm_bounds(role, quarter, lower)?
                    .map(i);
                vertices.push(BoundaryVertex { key, point })
            }
        }
    }
    for (index, e) in events.mouth().edges().iter().enumerate() {
        edges.push(BoundaryEdge {
            ends: e.vertices(),
            key: CurveKey::Mouth { edge: index },
        })
    }
    for station in [MeridianStation::Entry, MeridianStation::Floor] {
        for quarter in 0..4 {
            let role = station.event();
            edges.push(BoundaryEdge {
                ends: [
                    vertex_id(&vertices, event_key(role, quarter, false))?,
                    vertex_id(&vertices, event_key(role, (quarter + 1) % 4, true))?,
                ],
                key: CurveKey::MeridianCircle { station, quarter },
            })
        }
    }
    for t in events.trims() {
        for quarter in 0..4 {
            for lower in [false, true] {
                let [start, end] = t.end_events();
                edges.push(BoundaryEdge {
                    ends: [
                        vertex_id(&vertices, event_key(start, quarter, lower))?,
                        vertex_id(&vertices, event_key(end, quarter, lower))?,
                    ],
                    key: CurveKey::Trim {
                        envelope: t.envelope_role(),
                        wall: t.wall_kind(),
                        quarter,
                        lower,
                    },
                })
            }
        }
    }
    for quarter in 0..4 {
        checkpoint(cancel)?;
        for lower in [false, true] {
            edges.push(BoundaryEdge {
                ends: [
                    mouth_nodes.join(quarter, lower),
                    vertex_id(&vertices, event_key(EventRole::WallJoin, quarter, lower))?,
                ],
                key: CurveKey::WallJoin { quarter, lower },
            })
        }
        edges.push(BoundaryEdge {
            ends: [
                mouth_nodes.tip(quarter),
                vertex_id(&vertices, event_key(EventRole::WallContact, quarter, false))?,
            ],
            key: CurveKey::Center { quarter },
        });
        edges.push(BoundaryEdge {
            ends: [
                vertex_id(&vertices, event_key(EventRole::WallContact, quarter, false))?,
                vertex_id(&vertices, event_key(EventRole::FloorApex, 0, false))?,
            ],
            key: CurveKey::FloorGenerator { quarter },
        });
    }
    faces.push(BoundaryFace {
        key: FaceKey::Mouth,
        uses: (0..events.mouth().edges().len())
            .rev()
            .map(|edge| EdgeUse {
                edge,
                forward: false,
            })
            .collect(),
    });
    for quarter in 0..4 {
        checkpoint(cancel)?;
        let next = (quarter + 1) % 4;
        for kind in [WallKind::Outer, WallKind::Inner] {
            for lower in [true, false] {
                let role = match (kind, lower) {
                    (WallKind::Outer, true) => EdgeRole::OuterLower,
                    (WallKind::Outer, false) => EdgeRole::OuterUpper,
                    (WallKind::Inner, true) => EdgeRole::InnerLower,
                    (WallKind::Inner, false) => EdgeRole::InnerUpper,
                };
                let mut keys = vec![(mouth_edge(&events, quarter, role)?, true)];
                match (kind, lower) {
                    (WallKind::Outer, true) => {
                        keys.push((CurveKey::WallJoin { quarter, lower }, true));
                        keys.extend(trim_uses(&events, quarter, lower, None, Some(kind), false));
                    }
                    (WallKind::Outer, false) => {
                        keys.extend(trim_uses(&events, quarter, lower, None, Some(kind), true));
                        keys.push((CurveKey::WallJoin { quarter, lower }, false));
                    }
                    (WallKind::Inner, true) => {
                        keys.push((CurveKey::Center { quarter }, true));
                        keys.extend(trim_uses(&events, quarter, lower, None, Some(kind), false));
                        keys.push((CurveKey::WallJoin { quarter, lower }, false));
                    }
                    (WallKind::Inner, false) => {
                        keys.push((CurveKey::WallJoin { quarter, lower }, true));
                        keys.extend(trim_uses(&events, quarter, lower, None, Some(kind), true));
                        keys.push((CurveKey::Center { quarter }, false));
                    }
                }
                faces.push(BoundaryFace {
                    key: FaceKey::Wall {
                        kind,
                        quarter,
                        lower,
                    },
                    uses: uses(&edges, keys)?,
                });
            }
        }
        for role in [
            EnvelopeRole::EntryRoundover,
            EnvelopeRole::OuterFlank,
            EnvelopeRole::Floor,
        ] {
            let mut keys = vec![(
                match role {
                    EnvelopeRole::EntryRoundover => {
                        mouth_edge(&events, quarter, EdgeRole::Envelope)?
                    }
                    EnvelopeRole::OuterFlank => CurveKey::MeridianCircle {
                        station: MeridianStation::Entry,
                        quarter,
                    },
                    EnvelopeRole::Floor => CurveKey::MeridianCircle {
                        station: MeridianStation::Floor,
                        quarter,
                    },
                },
                true,
            )];
            keys.extend(trim_uses(&events, next, true, Some(role), None, true));
            if role == EnvelopeRole::Floor {
                keys.push((CurveKey::FloorGenerator { quarter: next }, true));
                keys.push((CurveKey::FloorGenerator { quarter }, false));
            } else {
                keys.push((
                    CurveKey::MeridianCircle {
                        station: if role == EnvelopeRole::EntryRoundover {
                            MeridianStation::Entry
                        } else {
                            MeridianStation::Floor
                        },
                        quarter,
                    },
                    false,
                ));
            }
            keys.extend(trim_uses(&events, quarter, false, Some(role), None, false));
            faces.push(BoundaryFace {
                key: FaceKey::Envelope { role, quarter },
                uses: uses(&edges, keys)?,
            });
        }
    }
    checkpoint(cancel)?;
    Assembly {
        events,
        vertices,
        edges,
        faces,
    }
    .finish()
}
