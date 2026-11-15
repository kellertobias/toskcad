//! Merge faces that partition one plane, and the straight edges those
//! partitions split, into the maximal faces and edges of the same boundary.
//!
//! Exchange writers divide a planar boundary into coplanar faces to carry
//! colours, pads or markings. The partition does not change the solid, but
//! its edges join two faces of one plane, which no analytic family expects.
//! Every removed edge is validated by the ordinary edge reader first, and
//! every merged face and edge is read again by the family readers, so the
//! pass can only remove topology that lies inside one transmitted plane.
use super::*;

#[derive(Clone, Copy, PartialEq)]
struct Use {
    edge: u32,
    forward: bool,
}

enum Bound {
    /// A bound without edges keeps its transmitted record.
    Kept(u32),
    Ring {
        outer: bool,
        uses: Vec<Use>,
    },
}

struct FaceModel {
    id: u32,
    surface: u32,
    forward: bool,
    plane: Option<Plane>,
    bounds: Vec<Bound>,
    dirty: bool,
    live: bool,
}

impl FaceModel {
    fn uses(&self) -> impl Iterator<Item = Use> + '_ {
        self.bounds.iter().flat_map(|bound| match bound {
            Bound::Kept(_) => [].iter().copied(),
            Bound::Ring { uses, .. } => uses.iter().copied(),
        })
    }
}

/// A straight edge that replaces a chain of collinear transmitted edges.
struct Merged {
    id: u32,
    ends: [u32; 2],
    line: u32,
    forward: bool,
}

struct Model {
    shell: u32,
    faces: Vec<FaceModel>,
    ends: BTreeMap<u32, [u32; 2]>,
    merged: Vec<Merged>,
}

/// Rewrite the arena when at least one solid carries a planar partition.
/// Any irregularity leaves the solid exactly as transmitted.
pub(super) fn unify<'a>(
    entities: &mut BTreeMap<u32, Vec<Entity<'a>>>,
    solids: &[u32],
    tolerance: f64,
) -> bool {
    let mut next = entities.keys().next_back().copied().unwrap_or(0);
    let models = solids
        .iter()
        .filter_map(|&solid| Model::plan(entities, solid, tolerance, &mut next))
        .collect::<Vec<_>>();
    let Some(records) = models
        .iter()
        .map(|model| model.records(&mut next))
        .collect::<Option<Vec<_>>>()
        .filter(|records| !records.is_empty())
    else {
        return false;
    };
    let before = owned(entities, solids);
    for (model, records) in models.iter().zip(records) {
        for (id, name, args) in records {
            let source = None;
            entities.insert(id, vec![Entity { source, name, args }]);
        }
        let faces = model.faces.iter().filter(|face| face.live);
        let faces = faces.map(|face| Value::Reference(face.id)).collect();
        entities.insert(
            model.shell,
            vec![Entity {
                source: None,
                name: "CLOSED_SHELL",
                args: vec![Value::String(""), Value::List(faces)],
            }],
        );
    }
    let after = owned(entities, solids);
    for id in before.difference(&after) {
        entities.remove(id);
    }
    true
}

/// The exclusively owned topology below the solids.
fn owned(entities: &BTreeMap<u32, Vec<Entity<'_>>>, solids: &[u32]) -> BTreeSet<u32> {
    let mut seen = BTreeSet::new();
    let mut pending = solids.to_vec();
    while let Some(id) = pending.pop() {
        if let Some([entity]) = entities.get(&id).map(Vec::as_slice)
            && compound::topology(entity.name)
            && seen.insert(id)
        {
            for arg in &entity.args {
                match arg {
                    Value::Reference(id) => pending.push(*id),
                    Value::List(items) => {
                        pending.extend(items.iter().filter_map(|item| reference(item).ok()));
                    }
                    _ => {}
                }
            }
        }
    }
    seen
}

fn distance_to_line(point: Point, origin: Point, direction: Point) -> f64 {
    length(cross(sub(point, origin), direction))
}

impl Model {
    fn plan(
        entities: &BTreeMap<u32, Vec<Entity<'_>>>,
        solid: u32,
        tolerance: f64,
        next: &mut u32,
    ) -> Option<Self> {
        let mut shell = Shell::new(entities, tolerance, true);
        let mut model = Self::read(&shell, solid)?;
        let merged_faces = model.merge_faces(&mut shell)?;
        model.merge_edges(&mut shell, next)?;
        (merged_faces || !model.merged.is_empty()).then_some(model)
    }

    fn read(shell: &Shell<'_, '_>, solid: u32) -> Option<Self> {
        let closed = reference(shell.args(solid, "MANIFOLD_SOLID_BREP").ok()?.get(1)?).ok()?;
        let mut ends = BTreeMap::new();
        let faces = shell
            .member_faces(solid)
            .ok()?
            .into_iter()
            .map(|id| Self::face(shell, id, &mut ends))
            .collect::<Option<Vec<_>>>()?;
        Some(Self {
            shell: closed,
            faces,
            ends,
            merged: Vec::new(),
        })
    }

    fn face(
        shell: &Shell<'_, '_>,
        id: u32,
        ends: &mut BTreeMap<u32, [u32; 2]>,
    ) -> Option<FaceModel> {
        let args = shell.args(id, "ADVANCED_FACE").ok()?;
        let surface = reference(args.get(2)?).ok()?;
        let forward = sense(args.get(3)?).ok()?;
        let mut edged = true;
        let bounds = refs(args.get(1)?)
            .ok()?
            .into_iter()
            .map(|bound| {
                let (_, wire, direction) = shell.bound(bound, true).ok()?;
                let outer = shell.simple(bound).ok()?.name == "FACE_OUTER_BOUND";
                if shell.simple(wire).ok()?.name != "EDGE_LOOP" {
                    edged = false;
                    return Some(Bound::Kept(bound));
                }
                let mut oriented = refs(shell.args(wire, "EDGE_LOOP").ok()?.get(1)?).ok()?;
                if !direction {
                    oriented.reverse();
                }
                let uses = oriented
                    .into_iter()
                    .map(|oriented| {
                        let args = shell.args(oriented, "ORIENTED_EDGE").ok()?;
                        let edge = reference(args.get(3)?).ok()?;
                        let forward = sense(args.get(4)?).ok()? == direction;
                        let curve = shell.args(edge, "EDGE_CURVE").ok()?;
                        let vertices = [
                            reference(curve.get(1)?).ok()?,
                            reference(curve.get(2)?).ok()?,
                        ];
                        ends.insert(edge, vertices);
                        Some(Use { edge, forward })
                    })
                    .collect::<Option<Vec<_>>>()?;
                Some(Bound::Ring { outer, uses })
            })
            .collect::<Option<Vec<_>>>()?;
        let plane = shell
            .args(surface, "PLANE")
            .ok()
            .filter(|_| edged)
            .and_then(|args| shell.placement(reference(args.get(1)?).ok()?).ok())
            .map(|mut plane| {
                if !forward {
                    plane.normal = plane.normal.map(|axis| -axis);
                }
                plane
            });
        Some(FaceModel {
            id,
            surface,
            forward,
            plane,
            bounds,
            dirty: false,
            live: true,
        })
    }

    fn start(&self, item: Use) -> u32 {
        self.ends[&item.edge][usize::from(!item.forward)]
    }

    fn end(&self, item: Use) -> u32 {
        self.ends[&item.edge][usize::from(item.forward)]
    }

    /// Whether every vertex of `face` lies in the oriented plane of `seed`.
    fn coplanar(&self, shell: &mut Shell<'_, '_>, seed: usize, face: usize) -> bool {
        let (Some(plane), Some(other)) = (&self.faces[seed].plane, &self.faces[face].plane) else {
            return false;
        };
        dot(plane.normal, other.normal) > 0.
            && length(cross(plane.normal, other.normal)) <= 1e-10
            && self.faces[face].uses().all(|item| {
                shell.vertex(self.start(item)).is_ok_and(|point| {
                    dot(sub(point, plane.origin), plane.normal).abs() <= shell.tolerance
                })
            })
    }

    fn merge_faces(&mut self, shell: &mut Shell<'_, '_>) -> Option<bool> {
        let mut incidence = BTreeMap::<u32, Vec<usize>>::new();
        for (index, face) in self.faces.iter().enumerate() {
            for item in face.uses() {
                incidence.entry(item.edge).or_default().push(index);
            }
        }
        let mut group = vec![usize::MAX; self.faces.len()];
        let mut changed = false;
        for seed in 0..self.faces.len() {
            if group[seed] != usize::MAX {
                continue;
            }
            group[seed] = seed;
            let mut members = vec![seed];
            let mut pending = vec![seed];
            while let Some(face) = pending.pop() {
                for item in self.faces[face].uses().collect::<Vec<_>>() {
                    for &other in &incidence[&item.edge] {
                        if group[other] == usize::MAX && self.coplanar(shell, seed, other) {
                            group[other] = seed;
                            members.push(other);
                            pending.push(other);
                        }
                    }
                }
            }
            if members.len() > 1 {
                self.merge_group(shell, &members)?;
                changed = true;
            }
        }
        Some(changed)
    }

    fn merge_group(&mut self, shell: &mut Shell<'_, '_>, members: &[usize]) -> Option<()> {
        let uses = members
            .iter()
            .flat_map(|&face| self.faces[face].uses())
            .collect::<Vec<_>>();
        let mut directions = BTreeMap::<u32, [usize; 2]>::new();
        for item in &uses {
            directions.entry(item.edge).or_default()[usize::from(item.forward)] += 1;
        }
        let mut starts = BTreeMap::new();
        for &item in &uses {
            match directions[&item.edge] {
                // An edge inside the plane is validated before it is dropped.
                [1, 1] => shell.edge(item.edge).ok().map(|_| ())?,
                [0, 1] | [1, 0] => starts
                    .insert(self.start(item), item)
                    .is_none()
                    .then_some(())?,
                _ => return None,
            }
        }
        let mut bounds = Vec::new();
        while let Some((first, mut item)) = starts.pop_first() {
            let mut ring = vec![item];
            while self.end(item) != first {
                item = starts.remove(&self.end(item))?;
                ring.push(item);
            }
            let outer = false;
            bounds.push(Bound::Ring { outer, uses: ring });
        }
        (!bounds.is_empty()).then_some(())?;
        for &face in &members[1..] {
            self.faces[face].live = false;
        }
        let seed = &mut self.faces[members[0]];
        seed.bounds = bounds;
        seed.dirty = true;
        Some(())
    }

    /// The transmitted line of a straight edge that passes the edge reader.
    fn line(&self, shell: &mut Shell<'_, '_>, edge: u32) -> Option<(u32, Point, Point)> {
        if let Some(merged) = self.merged.iter().find(|merged| merged.id == edge) {
            let (origin, direction, _) = shell.line(merged.line, false).ok()?;
            return Some((merged.line, origin, direction));
        }
        let line = reference(shell.args(edge, "EDGE_CURVE").ok()?.get(3)?).ok()?;
        (shell.simple(line).ok()?.name == "LINE").then_some(())?;
        shell.edge(edge).ok()?;
        let (origin, direction, _) = shell.line(line, false).ok()?;
        Some((line, origin, direction))
    }

    fn merge_edges(&mut self, shell: &mut Shell<'_, '_>, next: &mut u32) -> Option<()> {
        let mut edges = BTreeMap::<u32, BTreeSet<u32>>::new();
        let mut faces = BTreeMap::<u32, Vec<usize>>::new();
        for (index, face) in self.faces.iter().enumerate().filter(|(_, face)| face.live) {
            for item in face.uses() {
                for vertex in self.ends[&item.edge] {
                    edges.entry(vertex).or_default().insert(item.edge);
                }
                faces.entry(item.edge).or_default().push(index);
            }
        }
        for vertex in edges.keys().copied().collect::<Vec<_>>() {
            let [a, b] = *edges[&vertex].iter().copied().collect::<Vec<_>>() else {
                continue;
            };
            let far = |edge: u32| {
                let [start, end] = self.ends[&edge];
                if start == vertex { end } else { start }
            };
            let [from, to] = [far(a), far(b)];
            // Both edges pass the edge reader, which also reads their vertices.
            let (Some((line, origin, direction)), Some(_)) =
                (self.line(shell, a), self.line(shell, b))
            else {
                continue;
            };
            let [start, middle, end] = [from, vertex, to].map(|id| shell.points[&id]);
            if distance_to_line(end, origin, direction) > shell.tolerance
                || dot(sub(end, middle), sub(middle, start)) <= 0.
            {
                continue;
            }
            *next = next.checked_add(1)?;
            let id = *next;
            let forward = dot(direction, sub(end, start)) > 0.;
            let touched = [faces[&a].as_slice(), faces[&b].as_slice()].concat();
            for &face in &touched {
                self.splice(face, [a, b], id, vertex)?;
            }
            self.ends.insert(id, [from, to]);
            self.merged
                .retain(|merged| merged.id != a && merged.id != b);
            self.merged.push(Merged {
                id,
                ends: [from, to],
                line,
                forward,
            });
            faces.insert(id, touched);
            edges.remove(&vertex);
            for (vertex, old) in [(from, a), (to, b)] {
                let incident = edges.get_mut(&vertex)?;
                incident.remove(&old);
                incident.insert(id);
            }
        }
        Some(())
    }

    /// Replace each adjacent pair of uses of `pair` around `vertex` in one face.
    fn splice(&mut self, face: usize, pair: [u32; 2], id: u32, vertex: u32) -> Option<()> {
        let ends = &self.ends;
        let face = &mut self.faces[face];
        for bound in &mut face.bounds {
            let Bound::Ring { uses, .. } = bound else {
                continue;
            };
            while let Some(first) = (0..uses.len()).find(|&index| {
                let item = uses[index];
                pair.contains(&item.edge) && ends[&item.edge][usize::from(item.forward)] == vertex
            }) {
                let second = (first + 1) % uses.len();
                (pair.contains(&uses[second].edge) && uses[second].edge != uses[first].edge)
                    .then_some(())?;
                uses[first] = Use {
                    edge: id,
                    forward: uses[first].edge == pair[0],
                };
                uses.remove(second);
            }
            uses.iter()
                .all(|item| !pair.contains(&item.edge))
                .then_some(())?;
        }
        face.dirty = true;
        Some(())
    }

    /// The records of the merged edges and of every rewritten face.
    fn records(&self, next: &mut u32) -> Option<Vec<(u32, &'static str, Vec<Value<'static>>)>> {
        let flag = |value: bool| Value::Enum(if value { "T" } else { "F" });
        let mut fresh = || {
            *next = next.checked_add(1)?;
            Some(*next)
        };
        let mut records = Vec::new();
        for merged in &self.merged {
            let [start, end] = merged.ends.map(Value::Reference);
            let line = Value::Reference(merged.line);
            let args = vec![Value::String(""), start, end, line, flag(merged.forward)];
            records.push((merged.id, "EDGE_CURVE", args));
        }
        for face in self.faces.iter().filter(|face| face.live && face.dirty) {
            let mut bounds = Vec::new();
            for bound in &face.bounds {
                let (outer, uses) = match bound {
                    Bound::Kept(id) => {
                        bounds.push(Value::Reference(*id));
                        continue;
                    }
                    Bound::Ring { outer, uses } => (*outer, uses),
                };
                let mut oriented = Vec::new();
                for item in uses {
                    let id = fresh()?;
                    let edge = Value::Reference(item.edge);
                    let args = vec![
                        Value::String(""),
                        Value::Derived,
                        Value::Derived,
                        edge,
                        flag(item.forward),
                    ];
                    records.push((id, "ORIENTED_EDGE", args));
                    oriented.push(Value::Reference(id));
                }
                let wire = fresh()?;
                let args = vec![Value::String(""), Value::List(oriented)];
                records.push((wire, "EDGE_LOOP", args));
                let id = fresh()?;
                let name = if outer {
                    "FACE_OUTER_BOUND"
                } else {
                    "FACE_BOUND"
                };
                let args = vec![Value::String(""), Value::Reference(wire), flag(true)];
                records.push((id, name, args));
                bounds.push(Value::Reference(id));
            }
            let args = vec![
                Value::String(""),
                Value::List(bounds),
                Value::Reference(face.surface),
                flag(face.forward),
            ];
            records.push((face.id, "ADVANCED_FACE", args));
        }
        Some(records)
    }
}

impl Shell<'_, '_> {
    /// Refuse, before any face is measured, a shell that no regular family
    /// can describe: one with another support kind, or one whose edge joins
    /// more than two faces in balanced directions. The second is a pinched,
    /// non-manifold contact that lenient readers accept. Both verdicts are
    /// then independent of the order in which the faces are listed.
    pub(super) fn outside_regular_families(&self, faces: &[u32]) -> Result<()> {
        let mut uses = BTreeMap::<u32, [usize; 2]>::new();
        for &face in faces {
            let Some(model) = Model::face(self, face, &mut BTreeMap::new()) else {
                continue;
            };
            let kind = self.simple(model.surface).map(|surface| surface.name);
            if self.allow_circles
                && !matches!(
                    kind,
                    Ok("PLANE" | "CYLINDRICAL_SURFACE" | "CONICAL_SURFACE") | Err(_)
                )
            {
                return Err(ReadError::Unsupported);
            }
            for item in model.uses() {
                uses.entry(item.edge).or_default()[usize::from(item.forward)] += 1;
            }
        }
        if uses
            .values()
            .any(|&[reverse, forward]| reverse == forward && forward > 1)
        {
            return Err(ReadError::Unsupported);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
