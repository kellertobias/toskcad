//! Source-owned chart refinement; collapsed sides are semantic identities.
use super::*;
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(super) enum Row {
    Mouth,
    Entry,
    Floor,
    Join,
    Contact,
}
impl Row {
    fn event(self) -> EventRole {
        match self {
            Self::Mouth => EventRole::MouthEnvelope,
            Self::Entry => EventRole::EntryJoin,
            Self::Floor => EventRole::FloorJoin,
            Self::Join => EventRole::WallJoin,
            Self::Contact => EventRole::WallContact,
        }
    }
    fn from_event(role: EventRole) -> Result<Self, Error> {
        match role {
            EventRole::MouthEnvelope => Ok(Self::Mouth),
            EventRole::EntryJoin => Ok(Self::Entry),
            EventRole::FloorJoin => Ok(Self::Floor),
            EventRole::WallJoin => Ok(Self::Join),
            EventRole::WallContact => Ok(Self::Contact),
            EventRole::FloorApex => Err(Error::Invalid),
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(super) enum Node {
    Original(usize),
    Projection {
        floor: bool,
        quarter: u8,
        lower: bool,
    },
    Center,
    Apex,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(super) enum Curve {
    Circle(Row, u8),
    Trim(usize, u8, bool),
    Trace(usize, u8, bool),
    Connector(Row, u8, bool),
    Ray(Row, u8, bool),
    Generator(u8),
}
#[derive(Clone, Copy, Debug)]
pub(super) enum Side {
    Curve(usize),
    Collapsed(usize),
}
#[derive(Clone, Debug)]
pub(super) struct Edge {
    pub key: Curve,
    pub ends: [usize; 2],
}
#[derive(Clone, Debug)]
pub(super) struct Complex {
    pub points: Vec<[I; 3]>,
    pub edges: Vec<Edge>,
    pub sides: Vec<[Side; 4]>,
    pub certificate: TopologyCertificate,
}
fn original(atlas: &PatchAtlas<'_>, key: VertexKey) -> Result<Node, Error> {
    atlas
        .topology
        .vertices()
        .iter()
        .position(|v| v.source_identity() == key)
        .map(Node::Original)
        .ok_or(Error::Unresolved)
}
fn actual(atlas: &PatchAtlas<'_>, row: Row, q: u8, lower: bool) -> Result<Node, Error> {
    original(
        atlas,
        crate::phillips_drive::topology::event_key(row.event(), q, lower),
    )
}
fn projected(atlas: &PatchAtlas<'_>, row: Row, q: u8, lower: bool) -> Result<Node, Error> {
    match row {
        Row::Mouth => actual(atlas, row, q, lower),
        Row::Contact => original(
            atlas,
            VertexKey::Mouth {
                quarter: q,
                role: VertexRole::InnerTip,
            },
        ),
        Row::Join => original(
            atlas,
            VertexKey::Mouth {
                quarter: q,
                role: if lower {
                    VertexRole::LowerJoin
                } else {
                    VertexRole::UpperJoin
                },
            },
        ),
        Row::Entry | Row::Floor => Ok(Node::Projection {
            floor: row == Row::Floor,
            quarter: q,
            lower,
        }),
    }
}
fn node_point(atlas: &PatchAtlas<'_>, node: Node) -> Result<[I; 3], Error> {
    match node {
        Node::Original(index) => atlas
            .topology
            .vertices()
            .get(index)
            .map(|p| p.point_mm_bounds().map(b))
            .ok_or(Error::Invalid),
        Node::Center => Ok([I::point(0.); 3]),
        Node::Apex => Ok(atlas.topology.events().event(EventRole::FloorApex).point),
        Node::Projection {
            floor,
            quarter,
            lower,
        } => {
            let e = atlas.topology.events();
            let [x, y, z] = e
                .event(if floor {
                    EventRole::FloorJoin
                } else {
                    EventRole::EntryJoin
                })
                .point;
            let [x, y] = rotate(
                [x + e.beta_tangent * z, if lower { -y } else { y }],
                quarter,
            );
            Ok([x, y, I::point(0.)])
        }
    }
}
fn curve_ends(atlas: &PatchAtlas<'_>, curve: Curve) -> Result<[Node; 2], Error> {
    match curve {
        Curve::Circle(row, q) => Ok([
            actual(atlas, row, q, false)?,
            actual(atlas, row, (q + 1) % 4, true)?,
        ]),
        Curve::Trim(i, q, l) => {
            let [a, z] = atlas.topology.events().trims()[i].end_events();
            Ok([
                actual(atlas, Row::from_event(a)?, q, l)?,
                actual(atlas, Row::from_event(z)?, q, l)?,
            ])
        }
        Curve::Trace(i, q, l) => {
            let [a, z] = atlas.topology.events().trims()[i].end_events();
            Ok([
                projected(atlas, Row::from_event(a)?, q, l)?,
                projected(atlas, Row::from_event(z)?, q, l)?,
            ])
        }
        Curve::Connector(row, q, l) => {
            Ok([projected(atlas, row, q, l)?, actual(atlas, row, q, l)?])
        }
        Curve::Ray(row, q, l) => Ok([Node::Center, projected(atlas, row, q, l)?]),
        Curve::Generator(q) => Ok([actual(atlas, Row::Contact, q, false)?, Node::Apex]),
    }
}
fn point_on_curve(atlas: &PatchAtlas<'_>, curve: Curve, t: f64) -> Result<[I; 3], Error> {
    let u = I::point(t);
    let e = atlas.topology.events();
    match curve {
        Curve::Trim(index, q, l) | Curve::Trace(index, q, l) => {
            let [x, y, z] = e.trim_point_mm_bounds(index, t)?.map(b);
            let depth = if matches!(curve, Curve::Trace(..)) {
                I::point(0.)
            } else {
                z
            };
            let x = if matches!(curve, Curve::Trace(..)) {
                x + e.beta_tangent * z
            } else {
                x
            };
            let [x, y] = rotate([x, if l { -y } else { y }], q);
            Ok([x, y, depth])
        }
        Curve::Circle(row, q) => {
            let p = e.event(row.event()).point;
            let phi = (p[1] / p[0]).atan();
            let theta = phi + u * (pi() / I::point(2.) - I::point(2.) * phi);
            let c = e.source();
            let r = match row {
                Row::Mouth => c.envelope.entry,
                Row::Entry => {
                    c.envelope.entry - I::point(c.spec.mouth_radius_mm) * c.envelope.cos_a
                }
                Row::Floor => I::point(c.spec.g_diameter_mm) / I::point(2.),
                Row::Join => (p[0].square() + p[1].square()).sqrt(),
                Row::Contact => c.corners[0].b_half,
            };
            let [x, y] = rotate([r * theta.cos(), r * theta.sin()], q);
            Ok([x, y, p[2]])
        }
        Curve::Connector(..) | Curve::Ray(..) | Curve::Generator(..) => {
            let [a, z] = curve_ends(atlas, curve)?.map(|n| node_point(atlas, n));
            let [a, z] = [a?, z?];
            Ok(std::array::from_fn(|i| a[i] + u * (z[i] - a[i])))
        }
    }
}
#[derive(Default)]
struct Draft {
    nodes: Vec<Node>,
    points: Vec<[I; 3]>,
    edges: Vec<Edge>,
    sides: Vec<[Side; 4]>,
}
#[derive(Clone, Copy)]
enum SideSpec {
    Curve(Curve),
    Collapsed(Node),
}
impl Draft {
    fn node(&mut self, atlas: &PatchAtlas<'_>, key: Node) -> Result<usize, Error> {
        if let Some(index) = self.nodes.iter().position(|n| *n == key) {
            return Ok(index);
        }
        let point = node_point(atlas, key)?;
        let id = self.nodes.len();
        self.nodes.push(key);
        self.points.push(point);
        Ok(id)
    }
    fn side(&mut self, atlas: &PatchAtlas<'_>, mut key: Curve) -> Result<Side, Error> {
        match key {
            Curve::Connector(Row::Mouth, q, l) => {
                return self
                    .node(atlas, actual(atlas, Row::Mouth, q, l)?)
                    .map(Side::Collapsed);
            }
            Curve::Connector(Row::Contact, q, _) => key = Curve::Connector(Row::Contact, q, false),
            Curve::Ray(Row::Contact, q, _) => key = Curve::Ray(Row::Contact, q, false),
            _ => {}
        }
        if let Some(index) = self.edges.iter().position(|e| e.key == key) {
            return Ok(Side::Curve(index));
        }
        let [a, z] = curve_ends(atlas, key)?;
        let ends = [self.node(atlas, a)?, self.node(atlas, z)?];
        let index = self.edges.len();
        self.edges.push(Edge { key, ends });
        Ok(Side::Curve(index))
    }
}
impl Complex {
    pub(super) fn build(atlas: &PatchAtlas<'_>, cancel: &dyn Fn() -> bool) -> Result<Self, Error> {
        let mut result = Draft::default();
        let mut faces = Vec::new();
        for (index, patch) in atlas.patches.iter().enumerate() {
            checkpoint(cancel)?;
            let spec = match patch.law {
                Law::Envelope(i, q) => {
                    let [a, z] = atlas.topology.events().trims()[i].end_events();
                    [
                        SideSpec::Curve(Curve::Circle(Row::from_event(a)?, q)),
                        SideSpec::Curve(Curve::Circle(Row::from_event(z)?, q)),
                        SideSpec::Curve(Curve::Trim(i, q, false)),
                        SideSpec::Curve(Curve::Trim(i, (q + 1) % 4, true)),
                    ]
                }
                Law::Wall(i, q, l) => {
                    let [a, z] = atlas.topology.events().trims()[i].end_events();
                    [
                        SideSpec::Curve(Curve::Connector(Row::from_event(a)?, q, l)),
                        SideSpec::Curve(Curve::Connector(Row::from_event(z)?, q, l)),
                        SideSpec::Curve(Curve::Trace(i, q, l)),
                        SideSpec::Curve(Curve::Trim(i, q, l)),
                    ]
                }
                Law::Trace(i, q, l) => {
                    let [a, z] = atlas.topology.events().trims()[i].end_events();
                    [
                        SideSpec::Curve(Curve::Ray(Row::from_event(a)?, q, l)),
                        SideSpec::Curve(Curve::Ray(Row::from_event(z)?, q, l)),
                        SideSpec::Collapsed(Node::Center),
                        SideSpec::Curve(Curve::Trace(i, q, l)),
                    ]
                }
                Law::Gap(q) => [
                    SideSpec::Curve(Curve::Ray(Row::Mouth, q, false)),
                    SideSpec::Curve(Curve::Ray(Row::Mouth, (q + 1) % 4, true)),
                    SideSpec::Collapsed(Node::Center),
                    SideSpec::Curve(Curve::Circle(Row::Mouth, q)),
                ],
                Law::Floor(q) => [
                    SideSpec::Curve(Curve::Circle(Row::Contact, q)),
                    SideSpec::Collapsed(Node::Apex),
                    SideSpec::Curve(Curve::Generator(q)),
                    SideSpec::Curve(Curve::Generator((q + 1) % 4)),
                ],
            };
            let [a, b, c, d] = spec.map(|side| match side {
                SideSpec::Curve(curve) => result.side(atlas, curve),
                SideSpec::Collapsed(node) => result.node(atlas, node).map(Side::Collapsed),
            });
            let sides = [a?, b?, c?, d?];
            let mut uses = Vec::new();
            for (j, forward) in [(2, true), (1, true), (3, false), (0, false)] {
                if let Side::Curve(edge) = sides[j] {
                    uses.push(EdgeUse { edge, forward })
                }
            }
            if atlas.patch_winding(index)? < 0 {
                uses.reverse();
                for u in &mut uses {
                    u.forward = !u.forward
                }
            }
            faces.push(uses);
            result.sides.push(sides);
        }
        let certificate = crate::phillips_drive::topology::validate::certify(
            &result.points,
            &result.edges.iter().map(|e| e.ends).collect::<Vec<_>>(),
            &faces,
        )?;
        Ok(Self {
            points: result.points,
            edges: result.edges,
            sides: result.sides,
            certificate,
        })
    }
    pub(super) fn edge_point(
        &self,
        atlas: &PatchAtlas<'_>,
        edge: usize,
        t: f64,
    ) -> Result<[I; 3], Error> {
        if !t.is_finite() || !(0. ..=1.).contains(&t) {
            return Err(Error::Invalid);
        }
        let e = self.edges.get(edge).ok_or(Error::Invalid)?;
        if t == 0. || t == 1. {
            return Ok(self.points[e.ends[t as usize]]);
        }
        point_on_curve(atlas, e.key, t)
    }
}
#[cfg(test)]
mod tests;
