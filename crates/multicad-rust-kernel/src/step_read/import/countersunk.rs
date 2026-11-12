//! Full inward cone/cylinder chains, retaining their actual finite face keys.
use super::*;
use counterbored::Bore as SeatedBore;
impl Shell<'_, '_> {
    pub(super) fn countersunk_polygon(&self, faces: &[Face]) -> Option<ImportSolid> {
        let edges = self.indexed_edges(faces.len());
        faces.iter().enumerate().find_map(|(lower, face)| {
            let Face::Planar { origin, normal, .. } = face else {
                return None;
            };
            let direction = normal.map(|n| zero(-n));
            if !preferred_axis(direction) {
                return None;
            }
            self.match_countersunk_polygon(faces, &edges, lower, *origin, direction)
        })
    }
    fn match_countersunk_polygon(
        &self,
        faces: &[Face],
        edges: &[Vec<u32>],
        lower: usize,
        bottom: Point,
        direction: Point,
    ) -> Option<ImportSolid> {
        let pocketed::PolygonCaps {
            upper,
            height,
            mut covered,
            profile,
        } = self.polygon_caps(faces, edges, lower, bottom, direction)?;
        let mut bores = Vec::new();
        for &rim in &edges[lower] {
            if let Curve::Circle(circle) = self.edges[&rim].curve {
                bores.push(self.conical_seated_bore(
                    faces,
                    edges,
                    [lower, upper],
                    (rim, circle),
                    direction,
                    &mut covered,
                )?);
            }
        }
        if covered.len() != faces.len()
            || !bores.iter().any(|b| b.seats.iter().any(Option::is_some))
        {
            return None;
        }
        let keys = self.seated_source_keys(faces, edges)?;
        self.planar_profile(
            &profile,
            lower,
            direction,
            PlanarMatch::Countersunk(height, &bores, &keys, [lower, upper]),
        )
    }
    fn conical_seated_bore(
        &self,
        faces: &[Face],
        edges: &[Vec<u32>],
        caps: [usize; 2],
        mouth: (u32, CircleEdge),
        direction: Point,
        covered: &mut BTreeSet<usize>,
    ) -> Option<SeatedBore> {
        let (rim, mouth) = mouth;
        let mut current = (caps[0], rim, mouth);
        let mut seats = [None; 2];
        let mut core = None;
        let mut source_faces = Vec::new();
        loop {
            let (prior, rim, from) = current;
            let wall = self.adjacent_face(rim, prior)?;
            if !covered.insert(wall)
                || edges[wall].len() != 3
                || !edges[wall]
                    .iter()
                    .any(|id| matches!(self.edges[id].curve, Curve::Line))
            {
                return None;
            }
            source_faces.push(wall);
            // Exactly three edges, including a line and the known current
            // circle, leave at most one distinct far circular rim.
            let (far, to) = edges[wall].iter().find_map(|&id| {
                if id != rim {
                    if let Curve::Circle(c) = self.edges[&id].curve {
                        return Some((id, c));
                    }
                    None
                } else {
                    None
                }
            })?;
            let delta = sub(to.center_mm, from.center_mm);
            let travel = dot(delta, direction);
            if travel <= self.tolerance
                || length(sub(delta, direction.map(|n| n * travel))) > self.tolerance
            {
                return None;
            }
            let end = self.adjacent_face(far, wall)? == caps[1];
            match &faces[wall] {
                Face::Cylindrical { flux, .. }
                    if *flux < 0. && core.is_none() && from.radius_mm == to.radius_mm =>
                {
                    core = Some((from.radius_mm, source_faces.len() - 1));
                }
                Face::Conical(c) if !c.outward => {
                    if source_faces.len() == 1 && !end && from.radius_mm > to.radius_mm {
                        seats[0] = Some((from.radius_mm, travel));
                    } else if end && core.is_some() && to.radius_mm > from.radius_mm {
                        seats[1] = Some((to.radius_mm, travel));
                    } else {
                        return None;
                    }
                }
                _ => return None,
            }
            if end {
                return core.map(|(radius, core_face)| SeatedBore {
                    center: mouth.center_mm,
                    radius,
                    seats,
                    faces: source_faces,
                    core_face,
                });
            }
            current = (wall, far, to);
        }
    }
}

#[cfg(test)]
mod tests;
