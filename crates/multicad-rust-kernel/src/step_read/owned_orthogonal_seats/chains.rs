//! Original finite cap/cone/cylinder chains, without an occupied-set claim.
use super::*;
#[derive(Clone)]
pub struct BoreChain {
    caps: [u32; 2],
    rims: Vec<u32>,
    walls: Vec<u32>,
    circles: Vec<CircleEdge>,
    cylinder: usize,
}
impl BoreChain {
    pub fn cap_face_entities(&self) -> [u32; 2] {
        self.caps
    }
    pub fn rim_edge_entities(&self) -> &[u32] {
        &self.rims
    }
    pub fn wall_face_entities(&self) -> &[u32] {
        &self.walls
    }
    /// Binary64 geometric caches only. Exact numeric authority remains in the
    /// original Boundary's retained records/spellings and rim owner IDs.
    pub fn rim_circles(&self) -> &[CircleEdge] {
        &self.circles
    }
    pub fn cylinder_index(&self) -> usize {
        self.cylinder
    }
}
pub struct Chains<'boundary, 'document, 'input> {
    boundary: &'boundary Boundary<'document, 'input>,
    chains: Vec<BoreChain>,
    work: usize,
    bytes: usize,
}
impl<'boundary, 'document, 'input> Chains<'boundary, 'document, 'input> {
    pub fn boundary(&self) -> &'boundary Boundary<'document, 'input> {
        self.boundary
    }
    pub fn bores(&self) -> &[BoreChain] {
        &self.chains
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn query_budget_bytes(&self) -> usize {
        self.bytes
    }
    pub fn aggregate_known_work(&self) -> usize {
        self.boundary.work_used() + self.work
    }
}
fn adjacent(edge: &Edge, prior: usize) -> std::result::Result<usize, Error> {
    edge.uses
        .iter()
        .find(|u| u.face != prior)
        .map(|u| u.face)
        .ok_or(Error::Invalid)
}
impl<'document, 'input> Boundary<'document, 'input> {
    /// Account this additional query separately from immutable boundary and
    /// shared document prerequisites. Every nonplanar face must be covered by
    /// exactly one through chain; no wall is inferred from a nominal recipe.
    pub fn through_chains<'boundary>(
        &'boundary self,
        options: Options,
        interrupt: &dyn Fn() -> bool,
    ) -> std::result::Result<Chains<'boundary, 'document, 'input>, Error> {
        let mut budget = Budget {
            work: 0,
            maximum: options.maximum_work,
            interrupt,
        };
        budget.charge(0)?;
        if options.maximum_work == 0
            || options.maximum_work > 16_000_000
            || options.maximum_bytes == 0
            || options.maximum_bytes > 64 * 1024 * 1024
        {
            return Err(Error::Invalid);
        }
        // <=64 faces,256 edges, at most16 chains of four rims/three walls.
        // Capacity doubling, original-owner sets and live traversal scratch.
        const BYTES: usize = 64 * 1024;
        if options.maximum_bytes < BYTES {
            return Err(Error::Limit);
        }
        let mut covered = BTreeSet::new();
        let mut chains = Vec::new();
        for (cap, face) in self.faces.iter().enumerate() {
            let Face::Planar { normal, .. } = face else {
                continue;
            };
            for (&rim, edge) in &self.shell.edges {
                budget.charge(32 + edge.uses.len())?;
                let Curve::Circle(circle) = edge.curve else {
                    continue;
                };
                if !edge.uses.iter().any(|u| u.face == cap) {
                    continue;
                }
                let wall = adjacent(edge, cap)?;
                if covered.contains(&wall) {
                    continue;
                }
                if matches!(self.faces[wall], Face::Planar { .. }) {
                    return Err(Error::Unsupported);
                }
                if chains.len() == 16 {
                    return Err(Error::Limit);
                }
                let mut chain = BoreChain {
                    caps: [self.face_entities[cap]; 2],
                    rims: vec![rim],
                    walls: Vec::new(),
                    circles: vec![circle],
                    cylinder: usize::MAX,
                };
                let mut prior = cap;
                let mut current = rim;
                let direction = normal.map(|n| -n);
                loop {
                    budget.charge(64)?;
                    let wall = adjacent(&self.shell.edges[&current], prior)?;
                    if !covered.insert(wall) || chain.walls.len() == 3 {
                        return Err(Error::Unsupported);
                    }
                    let mut far = None;
                    for (&id, edge) in &self.shell.edges {
                        budget.charge(32 + edge.uses.len())?;
                        if id == current || !edge.uses.iter().any(|u| u.face == wall) {
                            continue;
                        }
                        if let Curve::Circle(circle) = edge.curve
                            && far.replace((id, circle)).is_some()
                        {
                            return Err(Error::Unsupported);
                        }
                    }
                    let (next, circle) = far.ok_or(Error::Unsupported)?;
                    let from = *chain.circles.last().ok_or(Error::Invalid)?;
                    let delta = sub(circle.center_mm, from.center_mm);
                    let travel = dot(delta, direction);
                    if travel <= self.shell.tolerance
                        || length(sub(delta, direction.map(|n| n * travel))) > self.shell.tolerance
                    {
                        return Err(Error::Unsupported);
                    }
                    let end = adjacent(&self.shell.edges[&next], wall)?;
                    let terminal = matches!(self.faces[end], Face::Planar { .. });
                    match &self.faces[wall] {
                        Face::Cylindrical { flux, .. }
                            if *flux < 0. && chain.cylinder == usize::MAX =>
                        {
                            chain.cylinder = chain.walls.len()
                        }
                        Face::Conical(cone) if !cone.outward => {
                            let entry = chain.walls.is_empty()
                                && !terminal
                                && from.radius_mm > circle.radius_mm;
                            let exit = terminal
                                && chain.cylinder != usize::MAX
                                && circle.radius_mm > from.radius_mm;
                            if !entry && !exit {
                                return Err(Error::Unsupported);
                            }
                        }
                        _ => return Err(Error::Unsupported),
                    }
                    chain.walls.push(self.face_entities[wall]);
                    chain.rims.push(next);
                    chain.circles.push(circle);
                    if terminal {
                        if end == cap || chain.cylinder == usize::MAX {
                            return Err(Error::Unsupported);
                        }
                        let Face::Planar {
                            normal: opposite, ..
                        } = self.faces[end]
                        else {
                            return Err(Error::Invalid);
                        };
                        if length(sub(opposite, direction)) > 1e-10 {
                            return Err(Error::Unsupported);
                        }
                        // A full circular planar floor is a blind pocket, not
                        // an original stock cap. Both through caps keep lines.
                        for cap in [cap, end] {
                            let mut line = false;
                            for edge in self.shell.edges.values() {
                                budget.charge(32 + edge.uses.len())?;
                                line |= matches!(edge.curve, Curve::Line)
                                    && edge.uses.iter().any(|u| u.face == cap);
                            }
                            if !line {
                                return Err(Error::Unsupported);
                            }
                        }
                        chain.caps[1] = self.face_entities[end];
                        break;
                    }
                    prior = wall;
                    current = next;
                }
                chains.push(chain);
            }
        }
        for (index, face) in self.faces.iter().enumerate() {
            budget.charge(32)?;
            if !matches!(face, Face::Planar { .. }) && !covered.contains(&index) {
                return Err(Error::Unsupported);
            }
        }
        if chains.is_empty() {
            return Err(Error::Unsupported);
        }
        Ok(Chains {
            boundary: self,
            chains,
            work: budget.work,
            bytes: BYTES,
        })
    }
}
