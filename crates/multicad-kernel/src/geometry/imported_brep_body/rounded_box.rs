//! A box with every edge filleted: an exact retained boundary of six inset
//! planes, twelve quarter cylinders and eight spherical vertex blends.
//!
//! The solid is the Minkowski sum of the inner box with a ball, so each
//! corner blend is an octant of the sphere about the inner box corner and
//! meets its three edge cylinders along quarter circles. The boundary is
//! validated and integrated independently and must reproduce the closed
//! form volume.
use super::*;
use crate::TreatmentEdge;
use crate::geometry::authored_sweep::{AuthoredNaming, FaceName};
use brep::{Boundary, Edge, Face, Frame, Loop, Shell};
use multicad_rust_kernel::orthogonal::Box3;
use std::f64::consts::PI;

type Corner = [bool; 3];

fn axis(i: usize, sign: f64) -> [f64; 3] {
    std::array::from_fn(|k| if k == i { sign } else { 0. })
}

fn sign(high: bool) -> f64 {
    if high { 1. } else { -1. }
}

fn corner_index(corner: Corner) -> usize {
    (0..3).map(|i| usize::from(corner[i]) << i).sum()
}

fn invalid_box() -> GeometryError {
    error(
        "fillet",
        GeometryErrorKind::InvalidResult,
        "Rust rounded box boundary does not reproduce its closed form",
    )
}

impl Solid {
    /// Fillet all twelve edges of one world-aligned box cell. Returns `None`
    /// unless `edges` names each straight edge of that box exactly once.
    pub(crate) fn rust_fillet_all_box_edges(
        &self,
        edges: &[TreatmentEdge],
        radius: f64,
        feature: EntityId,
        policy: GeometryPolicy,
        topology: &TopologySnapshot,
    ) -> Option<Result<Self, GeometryError>> {
        if edges.len() != 12 || self.policy != policy || policy.boolean_fuzzy_mm != 0. {
            return None;
        }
        let cells = self.orthogonal_cells().ok()??;
        let [stock] = cells.as_slice() else {
            return None;
        };
        let stock = *stock;
        // Key each selected edge by its axis and the box sides it lies on.
        let mut origins = BTreeMap::new();
        for edge in edges {
            let Some(AnalyticGeometry::Line { a_mm, b_mm }) = topology
                .elements
                .iter()
                .find(|element| element.key == edge.key)
                .and_then(|element| element.geometry)
            else {
                return None;
            };
            let mut differing = (0..3).filter(|&i| a_mm[i] != b_mm[i]);
            let along = differing.next()?;
            if differing.next().is_some()
                || a_mm[along].min(b_mm[along]) != stock.min[along]
                || a_mm[along].max(b_mm[along]) != stock.max[along]
            {
                return None;
            }
            let mut sides = [false; 3];
            for i in (0..3).filter(|&i| i != along) {
                sides[i] = if a_mm[i] == stock.max[i] {
                    true
                } else if a_mm[i] == stock.min[i] {
                    false
                } else {
                    return None;
                };
            }
            if origins
                .insert((along, sides), edge.origin.clone())
                .is_some()
            {
                return None;
            }
        }
        if origins.len() != 12 {
            return None;
        }
        Some(self.rounded_box(stock, radius, &origins, feature, policy, topology))
    }

    fn rounded_box(
        &self,
        stock: Box3,
        radius: f64,
        edge_origins: &BTreeMap<(usize, [bool; 3]), Origin>,
        feature: EntityId,
        policy: GeometryPolicy,
        topology: &TopologySnapshot,
    ) -> Result<Self, GeometryError> {
        let size: [f64; 3] = std::array::from_fn(|i| stock.max[i] - stock.min[i]);
        policy.length(radius)?;
        for side in size {
            policy.length(side - 2. * radius).map_err(|cause| {
                error(
                    "fillet",
                    cause.kind,
                    "Edge treatment failed: the fillet radius leaves a straight box face below the minimum feature size",
                )
            })?;
        }
        let tolerance = policy.linear_tolerance_mm;
        // Inner box corner of each octant: the center of its sphere blend.
        let center = |corner: Corner| -> [f64; 3] {
            std::array::from_fn(|i| {
                if corner[i] {
                    stock.max[i] - radius
                } else {
                    stock.min[i] + radius
                }
            })
        };
        // The vertex where corner `c` touches the box face normal to `i`.
        let mut vertices = Vec::with_capacity(24);
        let mut vertex = BTreeMap::new();
        for c in 0..8 {
            let corner: Corner = std::array::from_fn(|i| c >> i & 1 == 1);
            for i in 0..3 {
                let mut point = center(corner);
                point[i] += sign(corner[i]) * radius;
                vertex.insert((c, i), vertices.len());
                vertices.push(point);
            }
        }
        let mut edges = Vec::with_capacity(48);
        let mut add = |curve, ends: [usize; 2], vertices: &[[f64; 3]]| {
            let points = ends.map(|v| vertices[v]);
            let edge = Edge::new(curve, ends, points, tolerance).map_err(failure)?;
            edges.push(edge);
            Ok::<_, GeometryError>(edges.len() - 1)
        };
        // Straight trims: along each box edge, one line on each adjacent face.
        let mut lines = BTreeMap::new();
        for along in 0..3 {
            for sides in 0..4 {
                let mut corner = [false; 3];
                let [i, j] = [(along + 1) % 3, (along + 2) % 3];
                corner[i] = sides & 1 == 1;
                corner[j] = sides & 2 == 2;
                let mut end = corner;
                end[along] = true;
                for face in [i, j] {
                    let ends = [
                        vertex[&(corner_index(corner), face)],
                        vertex[&(corner_index(end), face)],
                    ];
                    let [a, b] = ends.map(|v| vertices[v]);
                    let curve = Curve::Line {
                        origin: a,
                        direction: std::array::from_fn(|k| b[k] - a[k]),
                    };
                    lines.insert((along, corner, face), add(curve, ends, &vertices)?);
                }
            }
        }
        // Quarter circles: at each corner, between its touches on faces i, j.
        let mut arcs = BTreeMap::new();
        for c in 0..8 {
            let corner: Corner = std::array::from_fn(|i| c >> i & 1 == 1);
            for k in 0..3 {
                let [i, j] = [(k + 1) % 3, (k + 2) % 3];
                let x = axis(i, sign(corner[i]));
                let y = axis(j, sign(corner[j]));
                let frame = Frame {
                    origin: center(corner),
                    x,
                    y,
                    z: brep_cross(x, y),
                };
                let ends = [vertex[&(c, i)], vertex[&(c, j)]];
                arcs.insert(
                    (c, k),
                    add(Curve::Circle { frame, radius }, ends, &vertices)?,
                );
            }
        }
        let mut faces = Vec::with_capacity(26);
        let mut names = Vec::with_capacity(26);
        // Six planar faces, each bounded by its four inset lines.
        for normal_axis in 0..3 {
            for high in [false, true] {
                let [i, j] = [(normal_axis + 1) % 3, (normal_axis + 2) % 3];
                let z = axis(normal_axis, sign(high));
                let x = axis(i, 1.);
                let mut origin: [f64; 3] =
                    std::array::from_fn(|k| 0.5 * (stock.min[k] + stock.max[k]));
                origin[normal_axis] = if high {
                    stock.max[normal_axis]
                } else {
                    stock.min[normal_axis]
                };
                let frame = Frame {
                    origin,
                    x,
                    y: brep_cross(z, x),
                    z,
                };
                // Walk the face rectangle through its corners in order.
                let mut uses = Vec::with_capacity(4);
                let rectangle = [[false, false], [true, false], [true, true], [false, true]];
                for step in 0..4 {
                    let [from, to] = [rectangle[step], rectangle[(step + 1) % 4]];
                    let corner_of = |[a, b]: [bool; 2]| {
                        let mut corner = [false; 3];
                        corner[normal_axis] = high;
                        corner[i] = a;
                        corner[j] = b;
                        corner
                    };
                    let (along, fixed) = if from[0] != to[0] {
                        (i, from)
                    } else {
                        (j, from)
                    };
                    let mut low = corner_of(fixed);
                    low[along] = false;
                    let edge = lines[&(along, low, normal_axis)];
                    let forward = corner_of(from)[along] == low[along];
                    uses.push((edge, forward));
                }
                faces.push(Face {
                    surface: Surface::Plane(frame),
                    same_sense: true,
                    loops: vec![Loop { uses, vertex: None }],
                });
                names.push(self.box_face_name(topology, normal_axis, high, stock));
            }
        }
        // Twelve quarter cylinders about the inner box edges.
        for along in 0..3 {
            for sides in 0..4 {
                let [i, j] = [(along + 1) % 3, (along + 2) % 3];
                let mut corner = [false; 3];
                corner[i] = sides & 1 == 1;
                corner[j] = sides & 2 == 2;
                let mut end = corner;
                end[along] = true;
                let x = axis(i, sign(corner[i]));
                let z = axis(along, 1.);
                let frame = Frame {
                    origin: center(corner),
                    x,
                    y: brep_cross(z, x),
                    z,
                };
                let uses = vec![
                    (lines[&(along, corner, i)], true),
                    (arcs[&(corner_index(end), along)], true),
                    (lines[&(along, corner, j)], false),
                    (arcs[&(corner_index(corner), along)], false),
                ];
                faces.push(Face {
                    surface: Surface::Cylinder { frame, radius },
                    same_sense: true,
                    loops: vec![Loop { uses, vertex: None }],
                });
                let mut key = [false; 3];
                key[i] = corner[i];
                key[j] = corner[j];
                names.push(FaceName {
                    origins: vec![edge_origins[&(along, key)].clone()],
                    support: None,
                    anchor: None,
                });
            }
        }
        // Eight sphere octants, each bounded by its three quarter circles.
        for c in 0..8 {
            let corner: Corner = std::array::from_fn(|i| c >> i & 1 == 1);
            let s = corner.map(sign);
            // Both poles lie outside the closed octant of the blend.
            let z = brep_unit([s[0], -s[1], 0.]);
            let x = [0., 0., 1.];
            let frame = Frame {
                origin: center(corner),
                x,
                y: brep_cross(z, x),
                z,
            };
            faces.push(Face {
                surface: Surface::Sphere { frame, radius },
                same_sense: true,
                loops: vec![Loop {
                    uses: (0..3).map(|k| (arcs[&(c, k)], true)).collect(),
                    vertex: None,
                }],
            });
            names.push(FaceName {
                origins: (0..3)
                    .map(|along| {
                        let mut key = corner;
                        key[along] = false;
                        edge_origins[&(along, key)].clone()
                    })
                    .collect(),
                support: None,
                anchor: None,
            });
        }
        let boundary = Boundary {
            vertices,
            edges,
            shells: vec![Shell {
                faces: (0..faces.len()).collect(),
                reversed: false,
            }],
            faces,
        };
        let solid = brep::Solid::new(boundary, tolerance, None).map_err(failure)?;
        let inner = size.map(|side| side - 2. * radius);
        let expected = inner[0] * inner[1] * inner[2]
            + 2. * radius * (inner[0] * inner[1] + inner[1] * inner[2] + inner[0] * inner[2])
            + PI * radius * radius * (inner[0] + inner[1] + inner[2])
            + 4. / 3. * PI * radius.powi(3);
        if (solid.volume_mm3() - expected).abs() > solid.volume_error_mm3().max(1e-9 * expected) {
            return Err(invalid_box());
        }
        Self::from_retained_brep(
            solid,
            feature,
            Rc::new(AuthoredNaming { faces: names }),
            policy,
        )
    }

    /// The source box face on `normal_axis` at its high or low side.
    fn box_face_name(
        &self,
        topology: &TopologySnapshot,
        normal_axis: usize,
        high: bool,
        stock: Box3,
    ) -> FaceName {
        let level = if high {
            stock.max[normal_axis]
        } else {
            stock.min[normal_axis]
        };
        let element = topology.elements.iter().find(|element| {
            element.key.kind == TopologyKind::Face
                && matches!(element.geometry, Some(AnalyticGeometry::Plane { frame })
                    if (frame.normal()[normal_axis] - sign(high)).abs() < 1e-12
                        && frame.origin_mm[normal_axis] == level)
        });
        let Some(element) = element else {
            return FaceName {
                origins: vec![],
                support: None,
                anchor: None,
            };
        };
        let mut origins: Vec<_> = element.origins.iter().cloned().collect();
        let anchor = element.frames.first();
        if let Some(anchor) = anchor
            && let Some(position) = origins.iter().position(|o| *o == anchor.origin)
        {
            origins.swap(0, position);
        }
        FaceName {
            anchor: anchor
                .filter(|anchor| origins.first() == Some(&anchor.origin))
                .map(|anchor| anchor.frame),
            origins,
            support: None,
        }
    }
}

fn brep_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn brep_unit(v: [f64; 3]) -> [f64; 3] {
    let length = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    v.map(|x| x / length)
}
