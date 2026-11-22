//! Shared finite helical end-ring and prismatic wall topology.

use super::*;
use multicad_rust_kernel::headed_screw::hex_vertices;

fn at_z([x, y]: [f64; 2], z: f64) -> Point {
    [x, y, z]
}

pub(super) fn hex_edge(across_flats: f64, side: u32) -> [[f64; 2]; 2] {
    let corners = hex_vertices(across_flats);
    [
        corners[side as usize],
        corners[(side as usize + 1) % corners.len()],
    ]
}

pub(super) fn edge_midpoint(across_flats: f64, side: u32, z: f64) -> Point {
    let [[ax, ay], [bx, by]] = hex_edge(across_flats, side);
    [ax.midpoint(bx), ay.midpoint(by), z]
}

pub(super) fn placed_frame(frame: PlaneFrame, origin: Point, u: Point, v: Point) -> ReferenceFrame {
    let axes = [frame.u, frame.v, frame.normal()];
    let direction = |local: Point| std::array::from_fn(|i| dot(local, axes.map(|axis| axis[i])));
    ReferenceFrame {
        origin_mm: position(frame, origin),
        u: direction(u),
        v: direction(v),
    }
}

pub(super) fn hex_wall_plane(
    frame: PlaneFrame,
    af: f64,
    i: u32,
    z: f64,
    sign: f64,
) -> ReferenceFrame {
    let [[ax, ay], [bx, by]] = hex_edge(af, i);
    let length = (bx - ax).portable_hypot(by - ay);
    placed_frame(
        frame,
        [ax, ay, z],
        [sign * (bx - ax) / length, sign * (by - ay) / length, 0.],
        [0., 0., 1.],
    )
}

pub(super) struct Boundary<'a, Role> {
    body: &'a HelicalBody,
    frame: PlaneFrame,
    roles: &'a [Role],
    elements: &'a mut Vec<TopologyElement>,
    face_indices: Vec<usize>,
    next_edge: u32,
    next_vertex: u32,
    inherit_faces: bool,
}

impl<'a, Role: Copy + PartialEq> Boundary<'a, Role> {
    pub(super) fn new(
        body: &'a HelicalBody,
        frame: PlaneFrame,
        roles: &'a [Role],
        elements: &'a mut Vec<TopologyElement>,
        inherit_faces: bool,
    ) -> Self {
        let face_indices = elements
            .iter()
            .enumerate()
            .filter_map(|(index, element)| {
                (element.key.kind == TopologyKind::Face).then_some(index)
            })
            .collect();
        Self {
            body,
            frame,
            roles,
            elements,
            face_indices,
            next_edge: (3 * body.rails.len()) as u32,
            next_vertex: (2 * body.rails.len()) as u32,
            inherit_faces,
        }
    }

    pub(super) fn face(&mut self, role: Role) -> &mut TopologyElement {
        let ordinal = self
            .roles
            .iter()
            .position(|candidate| *candidate == role)
            .unwrap();
        &mut self.elements[self.face_indices[ordinal]]
    }

    fn geometry(&mut self, name: String, geometry: AnalyticGeometry, faces: &[Role]) {
        self.element(name, geometry.topology_kind(), Some(geometry), faces);
    }

    /// A finite edge retained by the owning physical source, whose curve is
    /// outside the canonical analytic descriptor vocabulary (as with helices).
    pub(super) fn surface_edge(&mut self, name: String, faces: &[Role]) {
        self.element(name, TopologyKind::Edge, None, faces);
    }

    fn element(
        &mut self,
        name: String,
        kind: TopologyKind,
        geometry: Option<AnalyticGeometry>,
        faces: &[Role],
    ) {
        let origin = self.body.origin(name);
        let mut origins = BTreeSet::from([origin.clone()]);
        for &role in faces {
            let inherit_faces = self.inherit_faces;
            let face = self.face(role);
            if kind == TopologyKind::Edge {
                face.boundary.insert(origin.clone());
            }
            if inherit_faces {
                origins.extend(face.origins.iter().cloned());
            }
        }
        let next = if kind == TopologyKind::Edge {
            &mut self.next_edge
        } else {
            &mut self.next_vertex
        };
        self.elements.push(TopologyElement {
            key: TopologyKey {
                kind,
                ordinal: *next,
            },
            origins,
            boundary: BTreeSet::new(),
            frames: vec![],
            geometry,
        });
        *next += 1;
    }

    pub(super) fn arc(
        &mut self,
        name: String,
        z: f64,
        radius_mm: f64,
        angles: [f64; 2],
        faces: &[Role],
    ) {
        self.geometry(
            name,
            AnalyticGeometry::CircularArc {
                frame: placed_frame(self.frame, [0., 0., z], [1., 0., 0.], [0., 1., 0.]),
                radius_mm,
                start_rad: angles[0],
                sweep_rad: angles[1],
            },
            faces,
        );
    }

    pub(super) fn circle_at(
        &mut self,
        name: String,
        center: Point,
        radius_mm: f64,
        angles: [f64; 2],
        faces: &[Role],
    ) {
        self.geometry(
            name,
            AnalyticGeometry::CircularArc {
                frame: placed_frame(self.frame, center, [1., 0., 0.], [0., 1., 0.]),
                radius_mm,
                start_rad: angles[0],
                sweep_rad: angles[1],
            },
            faces,
        );
    }

    pub(super) fn vertex(&mut self, name: String, point: Point, faces: &[Role]) {
        self.geometry(
            name,
            AnalyticGeometry::Point {
                position_mm: position(self.frame, point),
            },
            faces,
        );
    }

    pub(super) fn line(&mut self, name: String, a: Point, b: Point, faces: &[Role]) {
        self.geometry(
            name,
            AnalyticGeometry::Line {
                a_mm: position(self.frame, a),
                b_mm: position(self.frame, b),
            },
            faces,
        );
    }

    pub(super) fn hex_prism(
        &mut self,
        name: &str,
        across_flats: f64,
        depths: [f64; 2],
        caps: [Role; 2],
        wall: impl Fn(u32) -> Role,
    ) {
        let corners = hex_vertices(across_flats);
        self.polygon_prism(name, corners, depths, caps, wall);
    }

    pub(super) fn polygon_prism(
        &mut self,
        name: &str,
        corners: [[f64; 2]; 6],
        depths: [f64; 2],
        caps: [Role; 2],
        wall: impl Fn(u32) -> Role,
    ) {
        for (i, &p) in corners.iter().enumerate() {
            let j = (i + 1) % corners.len();
            let previous = (i + corners.len() - 1) % corners.len();
            let side = wall(i as u32);
            let adjacent = wall(previous as u32);
            let [a, b] = depths.map(|z| at_z(p, z));
            self.line(format!("{name}/post/{i}"), a, b, &[adjacent, side]);
            for (cap, z) in depths.into_iter().enumerate() {
                self.line(
                    format!("{name}/rim/{cap}/{i}"),
                    at_z(p, z),
                    at_z(corners[j], z),
                    &[caps[cap], side],
                );
                self.geometry(
                    format!("{name}/corner/{cap}/{i}"),
                    AnalyticGeometry::Point {
                        position_mm: position(self.frame, at_z(p, z)),
                    },
                    &[caps[cap], adjacent, side],
                );
            }
        }
    }
}
