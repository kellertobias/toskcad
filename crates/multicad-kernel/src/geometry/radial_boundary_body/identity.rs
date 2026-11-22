//! Finite source values and incidence, independent of source table ordinals.
use multicad_core::project_file::ContentHash;
use multicad_rust_kernel::radial_tensor::{Axis, BoundaryInput, Curve, Surface};

pub(super) struct Identity {
    pub member: ContentHash,
    pub faces: Vec<ContentHash>,
    pub edges: Vec<ContentHash>,
    pub vertices: Vec<ContentHash>,
}
struct Key(Vec<u8>);
impl Key {
    fn new(kind: &[u8]) -> Self {
        Self(kind.to_vec())
    }
    fn number(&mut self, value: f64) {
        self.0.extend(
            (if value == 0. { 0. } else { value })
                .to_bits()
                .to_le_bytes(),
        );
    }
    fn count(&mut self, count: usize) {
        self.0.extend((count as u64).to_le_bytes());
    }
    fn point(&mut self, point: [f64; 3]) {
        for value in point {
            self.number(value);
        }
    }
    fn axis(&mut self, axis: Axis) {
        self.point(axis.origin);
        self.point(axis.normal);
        self.point(axis.reference);
    }
    fn hash(&mut self, value: ContentHash) {
        self.0.extend(value.0);
    }
    fn flag(&mut self, value: bool) {
        self.0.push(u8::from(value));
    }
    fn finish(self) -> ContentHash {
        ContentHash::of(&self.0)
    }
    fn hashes(&mut self, mut values: Vec<ContentHash>) {
        values.sort_unstable();
        self.count(values.len());
        for value in values {
            self.hash(value);
        }
    }
}

pub(super) fn keys(input: &BoundaryInput) -> Identity {
    let surfaces = input
        .surfaces
        .iter()
        .map(|surface| {
            let mut key = Key::new(b"multicad.retained.surface/1");
            match surface {
                Surface::Plane(axis) => {
                    key.0.push(0);
                    key.axis(*axis);
                }
                Surface::Cone {
                    axis,
                    radius_mm,
                    semi_angle_rad,
                } => {
                    key.0.push(1);
                    key.axis(*axis);
                    key.number(*radius_mm);
                    key.number(*semi_angle_rad);
                }
                Surface::Tensor(tensor) => {
                    key.0.push(2);
                    for span in tensor.spans() {
                        key.count(span);
                    }
                    key.count(tensor.poles().len());
                    for &point in tensor.poles() {
                        key.point(point);
                    }
                }
            }
            key.finish()
        })
        .collect::<Vec<_>>();
    let curves = input
        .curves
        .iter()
        .map(|curve| {
            let mut key = Key::new(b"multicad.retained.curve/1");
            match curve {
                Curve::Circle { axis, radius_mm } => {
                    key.0.push(0);
                    key.axis(*axis);
                    key.number(*radius_mm);
                }
                Curve::Spline { spans, poles } => {
                    key.0.push(1);
                    key.count(*spans);
                    key.count(poles.len());
                    for &point in poles {
                        key.point(point);
                    }
                }
            }
            key.finish()
        })
        .collect::<Vec<_>>();
    let geometry = input
        .edges
        .iter()
        .map(|edge| {
            let mut key = Key::new(b"multicad.retained.edge-geometry/1");
            key.hash(curves[edge.curve]);
            for vertex in edge.vertices {
                key.point(input.vertices[vertex]);
            }
            key.flag(edge.forward);
            key.finish()
        })
        .collect::<Vec<_>>();
    let faces = input
        .faces
        .iter()
        .map(|face| {
            let mut key = Key::new(b"multicad.retained.face/1");
            key.hash(surfaces[face.surface]);
            key.flag(face.forward);
            let uses = face
                .outer
                .iter()
                .map(|&(edge, forward)| (geometry[edge], forward))
                .collect::<Vec<_>>();
            key.count(uses.len());
            let start = (0..uses.len())
                .min_by(|a, b| {
                    (0..uses.len())
                        .map(|i| uses[(a + i) % uses.len()])
                        .cmp((0..uses.len()).map(|i| uses[(b + i) % uses.len()]))
                })
                .unwrap_or(0);
            for i in 0..uses.len() {
                let (edge, forward) = uses[(start + i) % uses.len()];
                key.hash(edge);
                key.flag(forward);
            }
            key.finish()
        })
        .collect::<Vec<_>>();
    let mut incidence = vec![Vec::new(); input.edges.len()];
    for (i, face) in input.faces.iter().enumerate() {
        for &(edge, forward) in &face.outer {
            incidence[edge].push((faces[i], forward));
        }
    }
    let edges = geometry
        .iter()
        .zip(incidence)
        .map(|(&geometry, mut incident)| {
            incident.sort_unstable();
            let mut key = Key::new(b"multicad.retained.edge-incidence/1");
            key.hash(geometry);
            key.count(incident.len());
            for (face, forward) in incident {
                key.hash(face);
                key.flag(forward);
            }
            key.finish()
        })
        .collect::<Vec<_>>();
    let mut incidence = vec![Vec::new(); input.vertices.len()];
    for (i, edge) in input.edges.iter().enumerate() {
        for (end, &vertex) in edge.vertices.iter().enumerate() {
            incidence[vertex].push((edges[i], end));
        }
    }
    let vertices = input
        .vertices
        .iter()
        .zip(incidence)
        .map(|(&point, mut incident)| {
            incident.sort_unstable();
            let mut key = Key::new(b"multicad.retained.vertex-incidence/1");
            key.point(point);
            key.count(incident.len());
            for (edge, end) in incident {
                key.hash(edge);
                key.count(end);
            }
            key.finish()
        })
        .collect::<Vec<_>>();
    let mut member = Key::new(b"multicad.retained.member/1");
    member.hashes(faces.clone());
    member.hashes(edges.clone());
    member.hashes(vertices.clone());
    Identity {
        member: member.finish(),
        faces,
        edges,
        vertices,
    }
}
