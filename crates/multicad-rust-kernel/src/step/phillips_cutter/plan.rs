use super::*;
use support::{LiteralSupport, Pose, Support, adjacent, b, map_point, norm};
pub(super) struct Plan<'a> {
    topology: BoundaryTopology<'a>,
    points: Vec<Point>,
    surfaces: Vec<LiteralSupport>,
    curves: Vec<curves::Curve>,
    owners: Vec<[usize; 2]>,
}
fn owners(topology: &BoundaryTopology<'_>) -> Vec<[usize; 2]> {
    let mut uses = vec![Vec::new(); topology.edges().len()];
    for (face, source) in topology.faces().iter().enumerate() {
        for edge in source.edge_uses() {
            uses[edge.edge].push(face);
        }
    }
    // The sealed topology has already certified exactly two opposite uses.
    uses.into_iter().map(|faces| [faces[0], faces[1]]).collect()
}
impl<'a> Plan<'a> {
    pub fn new(
        source: &'a Cutter,
        pose: &Pose,
        tolerance: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        if !tolerance.is_finite() || tolerance <= 0. {
            return Err(Error::InvalidGeometry);
        }
        let tolerance = tolerance.min(1e-7);
        let topology = source.boundary_topology_with_cancel(&|| interrupt.is_some_and(|f| f()))?;
        let mut points = Vec::new();
        for vertex in topology.vertices() {
            checkpoint(interrupt)?;
            let expected = vertex.point_mm_bounds().map(b);
            let point = pose.point(expected.map(I::midpoint));
            if point.into_iter().any(|v| !v.is_finite() || v.abs() > 1e12) {
                return Err(Error::Limit);
            }
            if points.contains(&point) {
                return Err(Error::Limit);
            }
            let exact = map_point(pose, expected);
            let error = norm(std::array::from_fn(|i| adjacent(point[i]) - exact[i]));
            if !error.valid() || error.hi > tolerance / 16. {
                return Err(Error::Limit);
            }
            points.push(point);
        }
        for edge in topology.edges() {
            let [a, z] = edge.vertices();
            let distance = norm(std::array::from_fn(|i| {
                adjacent(points[a][i]) - adjacent(points[z][i])
            }));
            if !distance.valid() || distance.lo <= 16. * tolerance {
                return Err(Error::Limit);
            }
        }
        let mut surfaces = Vec::new();
        for face in 0..topology.faces().len() {
            checkpoint(interrupt)?;
            let expected = Support::new(topology.face_surface(face)?)?;
            let literal = expected.literal(pose);
            let error = expected.distance_bound(pose, &literal, source.bounding_box_mm())?;
            if error > tolerance / 16. {
                return Err(Error::Limit);
            }
            surfaces.push(literal);
        }
        // Four source trims each produce eight corner orbits. Charge their
        // shared capacity during construction, before allocating world spans.
        let mut remaining = 375;
        let mut generate = |i| {
            let spans = topology.events().trim_carrier_with_cancel(
                i,
                tolerance / 64.,
                remaining,
                &|| interrupt.is_some_and(|f| f()),
            )?;
            remaining -= spans.len();
            Ok::<_, crate::phillips_drive::Error>(spans)
        };
        let spans = [generate(0)?, generate(1)?, generate(2)?, generate(3)?];
        let curves = (0..topology.edges().len())
            .map(|index| {
                checkpoint(interrupt)?;
                curves::Curve::new(&topology, index, pose, &points, &spans, tolerance)
            })
            .collect::<Result<Vec<_>, Error>>()?;
        let owners = owners(&topology);
        Ok(Self {
            topology,
            points,
            surfaces,
            curves,
            owners,
        })
    }
    pub fn emit(
        &self,
        file: &mut StepFile,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Entity, Error> {
        let mut vertices = Vec::new();
        for &point in &self.points {
            checkpoint(interrupt)?;
            vertices.push(file.vertex(point));
        }
        let mut surfaces = Vec::new();
        for surface in &self.surfaces {
            checkpoint(interrupt)?;
            surfaces.push(surface.emit(file));
        }
        let mut edges = Vec::new();
        for (index, curve) in self.curves.iter().enumerate() {
            checkpoint(interrupt)?;
            let ends = self.topology.edges()[index].vertices().map(|v| vertices[v]);
            let supports = self.owners[index].map(|i| surfaces[i]);
            edges.push(curve.emit(file, ends, supports, interrupt)?);
        }
        let mut faces = Vec::new();
        for (index, face) in self.topology.faces().iter().enumerate() {
            checkpoint(interrupt)?;
            let uses = face
                .edge_uses()
                .iter()
                .map(|u| file.oriented_edge(edges[u.edge], u.forward))
                .collect::<Vec<_>>();
            let bound = file.bound(&uses, true);
            faces.push(file.face_on(&[bound], surfaces[index], self.surfaces[index].forward));
        }
        Ok(file.solid(&faces))
    }
}
