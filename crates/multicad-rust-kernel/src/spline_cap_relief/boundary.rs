//! Shared finite curves and oriented face wires; counter islands are real caps.
use super::*;
type Point = [f64; 3];
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FaceRole {
    BaseCap {
        top: bool,
    },
    BaseWall {
        ring: usize,
        edge: usize,
    },
    MotifCap {
        motif: usize,
    },
    CounterCap {
        motif: usize,
        hole: usize,
    },
    MotifWall {
        motif: usize,
        ring: usize,
        edge: usize,
    },
}
#[derive(Clone, Debug, PartialEq)]
pub enum Curve {
    Line { a_mm: Point, b_mm: Point },
    Profile { edge: Edge, z_mm: f64 },
}
impl Curve {
    pub(crate) fn moments(&self, origin: [f64; 2]) -> [f64; 3] {
        match self {
            Self::Line { a_mm, b_mm } => crate::normal_sweep::Edge::Line {
                a: [a_mm[0], a_mm[1]],
                b: [b_mm[0], b_mm[1]],
            }
            .moments(origin),
            Self::Profile { edge, .. } => edge.moments(origin),
        }
    }
    pub fn point(&self, fraction: f64) -> Point {
        match self {
            Self::Line { a_mm, b_mm } => {
                if fraction == 0. {
                    *a_mm
                } else if fraction == 1. {
                    *b_mm
                } else {
                    std::array::from_fn(|i| a_mm[i] + fraction * (b_mm[i] - a_mm[i]))
                }
            }
            Self::Profile { edge, z_mm } => {
                let p = edge.point(fraction);
                [p[0], p[1], *z_mm]
            }
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct BoundaryEdge {
    pub vertices: [usize; 2],
    pub curve: Curve,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EdgeUse {
    pub edge: usize,
    pub forward: bool,
}
#[derive(Clone, Debug, PartialEq)]
pub enum Surface {
    Plane {
        z_mm: f64,
        upward: bool,
    },
    Extrusion {
        edge: Edge,
        span_mm: [f64; 2],
        inward: bool,
    },
}
#[derive(Clone, Debug, PartialEq)]
pub struct BoundaryFace {
    pub role: FaceRole,
    pub surface: Surface,
    pub loops: Vec<Vec<EdgeUse>>,
}
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Boundary {
    pub vertices: Vec<Point>,
    pub edges: Vec<BoundaryEdge>,
    pub faces: Vec<BoundaryFace>,
}
fn reverse(mut edges: Vec<EdgeUse>) -> Vec<EdgeUse> {
    edges.reverse();
    for edge in &mut edges {
        edge.forward = !edge.forward;
    }
    edges
}
struct Rims {
    low: Vec<EdgeUse>,
    high: Vec<EdgeUse>,
    walls: Vec<BoundaryFace>,
}
fn ring(
    graph: &mut Boundary,
    source: &[Edge],
    span: [f64; 2],
    inward: bool,
    role: impl Fn(usize) -> FaceRole,
    budget: &mut validation::Budget<'_>,
) -> Result<Rims, Error> {
    let vertices = span.map(|z| {
        source
            .iter()
            .map(|edge| {
                let p = edge.point(0.);
                let id = graph.vertices.len();
                graph.vertices.push([p[0], p[1], z]);
                id
            })
            .collect::<Vec<_>>()
    });
    let rims: [Vec<EdgeUse>; 2] = std::array::from_fn(|cap| {
        source
            .iter()
            .enumerate()
            .map(|(i, edge)| {
                let id = graph.edges.len();
                graph.edges.push(BoundaryEdge {
                    vertices: [vertices[cap][i], vertices[cap][(i + 1) % source.len()]],
                    curve: Curve::Profile {
                        edge: edge.clone(),
                        z_mm: span[cap],
                    },
                });
                EdgeUse {
                    edge: id,
                    forward: true,
                }
            })
            .collect()
    });
    let rails = (0..source.len())
        .map(|i| {
            let id = graph.edges.len();
            let ends = [vertices[0][i], vertices[1][i]];
            graph.edges.push(BoundaryEdge {
                vertices: ends,
                curve: Curve::Line {
                    a_mm: graph.vertices[ends[0]],
                    b_mm: graph.vertices[ends[1]],
                },
            });
            id
        })
        .collect::<Vec<_>>();
    let mut walls = Vec::with_capacity(source.len());
    for (i, edge) in source.iter().enumerate() {
        validation::consume(budget)?;
        let next = (i + 1) % source.len();
        let mut wire = vec![
            rims[0][i],
            EdgeUse {
                edge: rails[next],
                forward: true,
            },
            EdgeUse {
                edge: rims[1][i].edge,
                forward: false,
            },
            EdgeUse {
                edge: rails[i],
                forward: false,
            },
        ];
        if inward {
            wire = reverse(wire);
        }
        walls.push(BoundaryFace {
            role: role(i),
            surface: Surface::Extrusion {
                edge: edge.clone(),
                span_mm: span,
                inward,
            },
            loops: vec![wire],
        });
    }
    let [low, high] = rims;
    Ok(Rims { low, high, walls })
}
fn cap(role: FaceRole, z: f64, upward: bool, loops: Vec<Vec<EdgeUse>>) -> BoundaryFace {
    BoundaryFace {
        role,
        surface: Surface::Plane { z_mm: z, upward },
        loops,
    }
}
pub(super) fn build(
    relief: &Relief,
    budget: &mut validation::Budget<'_>,
) -> Result<Boundary, Error> {
    let mut graph = Boundary::default();
    let mut low = vec![];
    let mut high = vec![];
    let mut walls = vec![];
    for (r, source) in relief.base.rings().enumerate() {
        let rims = ring(
            &mut graph,
            source,
            [0., relief.base.height_mm()],
            false,
            |edge| FaceRole::BaseWall { ring: r, edge },
            budget,
        )?;
        low.push(reverse(rims.low));
        high.push(rims.high);
        walls.extend(rims.walls);
    }
    let mut faces = vec![];
    let raised = relief.kind == Kind::Raised;
    for (m, motif) in relief.motifs.iter().enumerate() {
        let level = relief.levels[m];
        let span = if raised {
            [relief.base.height_mm(), level]
        } else {
            [level, relief.base.height_mm()]
        };
        let mut far = vec![];
        for (r, source) in motif.rings().enumerate() {
            let rims = ring(
                &mut graph,
                source,
                span,
                !raised,
                |edge| FaceRole::MotifWall {
                    motif: m,
                    ring: r,
                    edge,
                },
                budget,
            )?;
            let (at_cap, away) = if raised {
                (rims.low, rims.high)
            } else {
                (rims.high, rims.low)
            };
            if r == 0 {
                high.push(reverse(at_cap));
            } else {
                faces.push(cap(
                    FaceRole::CounterCap {
                        motif: m,
                        hole: r - 1,
                    },
                    relief.base.height_mm(),
                    true,
                    vec![reverse(at_cap)],
                ));
            }
            far.push(away);
            walls.extend(rims.walls);
        }
        faces.push(cap(FaceRole::MotifCap { motif: m }, level, true, far));
    }
    graph
        .faces
        .push(cap(FaceRole::BaseCap { top: false }, 0., false, low));
    graph.faces.push(cap(
        FaceRole::BaseCap { top: true },
        relief.base.height_mm(),
        true,
        high,
    ));
    graph.faces.extend(faces);
    graph.faces.extend(walls);
    Ok(graph)
}
