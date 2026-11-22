use super::{
    Curve, Edge, Error, Face, OrientedEdge, Point, Shell, Support, cross, cylinder_patch::ANGULAR,
    dot, sub,
};
use crate::{
    normal_sweep::{CertifiedPiercing, Edge as ProfileEdge, PiercedSweep, Piercing, Sweep},
    trim_curve::{Ellipse, EllipseArc},
};
use multicad_precision::PortableMath as _;
use std::{collections::BTreeMap, f64::consts::FRAC_PI_2};

mod sweep_shell;
pub use sweep_shell::{normal_sweep_shell, pierced_sweep_shell};

#[derive(Clone, Copy)]
pub struct SheetWebEdge {
    pub ends: [[f64; 2]; 2],
    pub source: usize,
    pub flange: Option<usize>,
}

pub struct SheetFlange<'a> {
    /// Counterclockwise outward/thickness section, closed by a final root line.
    pub section: &'a [(ProfileEdge, usize)],
    pub ends: [usize; 2],
}

/// A rectangular or circular cut through a flange, in the sheet-local frame.
pub struct SheetFlangeOpening<'a> {
    /// Section edge ordinals of the entry and exit walls.
    pub walls: [usize; 2],
    pub origin_mm: Point,
    pub u: Point,
    pub v: Point,
    pub through_mm: Point,
    /// Counterclockwise outline and authored face provenance per edge.
    pub outline: &'a [(ProfileEdge, usize)],
}

/// Join rectangular stock and up to four separated flange sweeps at full
/// shared root faces. All edges and wall supports remain analytic.
pub fn formed_sheet(
    web: &[SheetWebEdge],
    holes: &[Vec<(ProfileEdge, usize)>],
    caps: [usize; 2],
    flanges: &[SheetFlange<'_>],
    signed_thickness: f64,
    tolerance: f64,
    minimum: f64,
) -> Result<Shell, Error> {
    formed_sheet_with_openings(
        web,
        holes,
        caps,
        flanges,
        &[],
        signed_thickness,
        tolerance,
        minimum,
    )
}

/// Join a web and separated flange sweeps with certified through openings.
/// A nonempty `openings` slice supplies one list per flange in flange order.
#[allow(clippy::too_many_arguments)]
pub fn formed_sheet_with_openings(
    web: &[SheetWebEdge],
    holes: &[Vec<(ProfileEdge, usize)>],
    caps: [usize; 2],
    flanges: &[SheetFlange<'_>],
    openings: &[&[SheetFlangeOpening<'_>]],
    signed_thickness: f64,
    tolerance: f64,
    minimum: f64,
) -> Result<Shell, Error> {
    formed_sheet_with_overhangs(
        web,
        holes,
        caps,
        flanges,
        openings,
        &[],
        signed_thickness,
        tolerance,
        minimum,
    )
}

/// As [`formed_sheet_with_openings`], where a flange root may continue past
/// its welded web edge. A nonempty `overhangs` slice gives, per flange, the
/// root length before the start and after the end of that edge; the web
/// must leave that part of the root free, as a corner relief does.
#[allow(clippy::too_many_arguments)]
pub fn formed_sheet_with_overhangs(
    web: &[SheetWebEdge],
    holes: &[Vec<(ProfileEdge, usize)>],
    caps: [usize; 2],
    flanges: &[SheetFlange<'_>],
    openings: &[&[SheetFlangeOpening<'_>]],
    overhangs: &[[f64; 2]],
    signed_thickness: f64,
    tolerance: f64,
    minimum: f64,
) -> Result<Shell, Error> {
    if (!openings.is_empty() && openings.len() != flanges.len())
        || (!overhangs.is_empty() && overhangs.len() != flanges.len())
        || overhangs
            .iter()
            .flatten()
            .any(|length| !length.is_finite() || *length < 0.)
    {
        return Err(Error::Invalid);
    }
    if flanges.is_empty() || flanges.len() > 4 || web.len() > 36 {
        return Err(Error::Unsupported);
    }
    let thickness = signed_thickness.abs();
    let base_ring = web
        .iter()
        .map(|edge| {
            let [a, b] = edge.ends;
            (ProfileEdge::Line { a, b }, edge.source)
        })
        .collect::<Vec<_>>();
    let rings: Vec<_> = std::iter::once(base_ring.as_slice())
        .chain(holes.iter().map(Vec::as_slice))
        .collect();
    let base = checked_sweep(&rings, thickness, tolerance, minimum)?;
    let up = [0., 0., signed_thickness.signum()];
    let web_frame = SweepFrame {
        origin: [0.; 3],
        u: [1., 0., 0.],
        v: [0., 1., 0.],
        direction: up,
    };
    let mut builder = Builder::new(tolerance);
    let mut volume = base.volume_mm3();
    let mut bounds = web_frame.bounds(&base);
    let mut prior_flanges = Vec::<PlacedFlange>::with_capacity(flanges.len());
    for (index, flange) in flanges.iter().enumerate() {
        let mut matches = web.iter().filter(|edge| edge.flange == Some(index));
        let hinge = matches.next().ok_or(Error::Invalid)?;
        if matches.next().is_some() {
            return Err(Error::Invalid);
        }
        let [a, b] = hinge.ends.map(|[x, y]| [x, y, 0.]);
        let along = sub(b, a);
        let welded = dot(along, along).sqrt();
        let along = scale(along, welded.recip());
        let outward = [along[1], -along[0], 0.];
        let [before, after] = overhangs.get(index).copied().unwrap_or([0., 0.]);
        let a = sub(a, scale(along, before));
        let width = before + welded + after;
        let section = checked_sweep(&[flange.section], width, tolerance, minimum)?;
        let (&root, walls) = section.rings()[0].split_last().ok_or(Error::Invalid)?;
        if root
            != (ProfileEdge::Line {
                a: [0., thickness],
                b: [0., 0.],
            })
        {
            return Err(Error::Unsupported);
        }
        let back = back_slab(walls, thickness, tolerance, minimum)?;
        let wall_count = walls.len();
        let reach = section.support([1., 0., 0.]);
        let free = |box_: [[f64; 2]; 2]| {
            web.iter()
                .all(|edge| !web_enters(edge.ends, a, along, outward, box_, tolerance))
        };
        // The front material and any exposed root face stay clear of the web.
        let mut clear = free([[tolerance, reach + minimum], [-minimum, width + minimum]]);
        let welded_span = [before, before + welded];
        if before > 0. {
            clear &= free([
                [-minimum, tolerance],
                [-minimum, welded_span[0] - tolerance],
            ]);
        }
        if after > 0. {
            clear &= free([
                [-minimum, tolerance],
                [welded_span[1] + tolerance, width + minimum],
            ]);
        }
        if !clear {
            return Err(Error::Unsupported);
        }
        let frame = SweepFrame {
            origin: a,
            u: outward,
            v: up,
            direction: along,
        };
        let placed = PlacedFlange {
            frame,
            section,
            width,
            back,
        };
        if prior_flanges
            .iter()
            .any(|other| !placed.separated_from(other, minimum))
        {
            return Err(Error::Unsupported);
        }
        let section = &placed.section;
        let omit = (0..flange.section.len()).map(|i| i == wall_count);
        let first_wall = builder.faces.len();
        builder.root_split = (before > 0. || after > 0.).then_some(RootSplit {
            points: [[0., thickness], [0., 0.]],
            welded: welded_span,
        });
        let swept = builder.sweep(&[flange.section], section, frame, flange.ends, omit);
        builder.root_split = None;
        swept?;
        volume += match openings.get(index).copied() {
            Some(openings) if !openings.is_empty() => {
                builder.pierce_flange(section, openings, frame, first_wall, thickness, minimum)?
            }
            _ => section.volume_mm3(),
        };
        let flange_bounds = frame.bounds(section);
        for i in 0..3 {
            bounds[0][i] = bounds[0][i].min(flange_bounds[0][i]);
            bounds[1][i] = bounds[1][i].max(flange_bounds[1][i]);
        }
        prior_flanges.push(placed);
    }
    if web
        .iter()
        .filter_map(|edge| edge.flange)
        .any(|index| index >= flanges.len())
    {
        return Err(Error::Invalid);
    }
    builder
        .sweep(
            &rings,
            &base,
            web_frame,
            caps,
            web.iter().map(|e| e.flange.is_some()),
        )
        .and_then(|()| builder.finish(bounds, volume))
}

/// The section walls behind the hinge plane, if any, as one height range
/// that lies wholly beyond one web face. Walls in front of the plane may
/// touch it only where they meet the root.
///
/// The section's material behind the plane is bounded by those walls and a
/// chord on the plane between their crossings. When every such wall lies
/// strictly below the web, or every one strictly above it, the material
/// does too, so a hem folded under or over the web cannot reach it.
fn back_slab(
    walls: &[ProfileEdge],
    thickness: f64,
    tolerance: f64,
    minimum: f64,
) -> Result<Option<[f64; 2]>, Error> {
    let (root_a, root_b) = ([0., thickness], [0., 0.]);
    let (mut below, mut above, mut range) = (true, true, None::<[f64; 2]>);
    for &edge in walls {
        let on_root = [edge.start(), edge.end()].map(|p| p == root_a || p == root_b);
        if on_root[0] && on_root[1] {
            return Err(Error::Unsupported);
        }
        let extent = edge.support([-1., 0.]);
        if extent <= tolerance {
            if !on_root[0] && !on_root[1] && extent > -minimum {
                return Err(Error::Unsupported);
            }
            continue;
        }
        let [low, high] = [-edge.support([0., -1.]), edge.support([0., 1.])];
        below &= high < -minimum;
        above &= low > thickness + minimum;
        range = Some(range.map_or([low, high], |[l, h]| [l.min(low), h.max(high)]));
    }
    match range {
        Some(_) if !below && !above => Err(Error::Unsupported),
        range => Ok(range),
    }
}

/// Whether a web boundary edge reaches a box in a hinge's outward/along
/// coordinates. A box in front of the hinge that no edge reaches lies
/// outside the web, as the web lies behind its own hinge edge there.
fn web_enters(
    ends: [[f64; 2]; 2],
    origin: Point,
    along: Point,
    outward: Point,
    [[out_low, out_high], [along_low, along_high]]: [[f64; 2]; 2],
    tolerance: f64,
) -> bool {
    let local = |[x, y]: [f64; 2]| {
        let offset = sub([x, y, 0.], origin);
        [dot(offset, outward), dot(offset, along)]
    };
    let [p, q] = ends.map(local);
    let lower = [out_low, along_low];
    let upper = [out_high, along_high];
    let (mut enter, mut exit) = (0f64, 1f64);
    for i in 0..2 {
        let delta = q[i] - p[i];
        if delta == 0. {
            if p[i] < lower[i] || p[i] > upper[i] {
                return false;
            }
            continue;
        }
        let [s, t] = [(lower[i] - p[i]) / delta, (upper[i] - p[i]) / delta];
        enter = enter.max(s.min(t));
        exit = exit.min(s.max(t));
    }
    if enter > exit {
        return false;
    }
    // A box straddling the hinge plane guards an exposed root face, where any
    // contact counts. In front of the hinge, contact only on the box's
    // hinge-side face stays at the tolerance already granted to the weld.
    out_low < 0. || p[0] + (q[0] - p[0]) * 0.5 * (enter + exit) > tolerance
}

/// A flange sweep placed in the web frame, with its height range behind
/// the hinge plane when it folds back over or under the web.
struct PlacedFlange {
    frame: SweepFrame,
    section: Sweep,
    width: f64,
    back: Option<[f64; 2]>,
}

impl PlacedFlange {
    /// Frame-aligned boxes that cover the flange material: the part in front
    /// of the hinge plane and, for a hem, the slab behind it.
    fn parts(&self) -> Vec<[[f64; 2]; 3]> {
        let s = |d| self.section.support(d);
        let height = [-s([0., -1., 0.]), s([0., 1., 0.])];
        let mut parts = vec![[[0., s([1., 0., 0.])], height, [0., self.width]]];
        if let Some(range) = self.back {
            parts.push([[-s([-1., 0., 0.]), 0.], range, [0., self.width]]);
        }
        parts
    }

    fn support(&self, part: [[f64; 2]; 3], direction: Point) -> f64 {
        let local = self.frame.coordinates(direction);
        let boxed: f64 = (0..3)
            .map(|i| (local[i] * part[i][0]).max(local[i] * part[i][1]))
            .sum();
        dot(self.frame.origin, direction) + boxed.min(self.section.support(local))
    }

    /// Whether a plane normal to an outward, along or thickness axis of either
    /// flange separates every pair of their material parts.
    fn separated_from(&self, other: &Self, minimum: f64) -> bool {
        let axes = [self.frame, other.frame]
            .into_iter()
            .flat_map(|frame| [frame.u, frame.direction])
            .chain([self.frame.v]);
        let directions: Vec<Point> = axes.flat_map(|axis| [axis, scale(axis, -1.)]).collect();
        self.parts().into_iter().all(|mine| {
            other.parts().into_iter().all(|theirs| {
                directions.iter().any(|&d| {
                    self.support(mine, d) + other.support(theirs, scale(d, -1.)) < -minimum
                })
            })
        })
    }
}

fn checked_sweep(
    rings: &[&[(ProfileEdge, usize)]],
    height: f64,
    tolerance: f64,
    minimum: f64,
) -> Result<Sweep, Error> {
    let edges: Vec<Vec<_>> = rings
        .iter()
        .map(|ring| ring.iter().map(|&(edge, _)| edge).collect())
        .collect();
    Sweep::new(edges, height, tolerance, minimum).map_err(sweep_error)
}

fn sweep_error(error: crate::normal_sweep::Error) -> Error {
    match error {
        crate::normal_sweep::Error::Limit => Error::Limit,
        _ => Error::Invalid,
    }
}

#[derive(Clone, Copy)]
struct SweepFrame {
    origin: Point,
    u: Point,
    v: Point,
    direction: Point,
}

impl SweepFrame {
    fn point(self, [u, v]: [f64; 2], level: f64) -> Point {
        add(
            add(add(self.origin, scale(self.u, u)), scale(self.v, v)),
            scale(self.direction, level),
        )
    }

    fn vector(self, [u, v, axial]: Point) -> Point {
        add(
            add(scale(self.u, u), scale(self.v, v)),
            scale(self.direction, axial),
        )
    }

    fn coordinates(self, vector: Point) -> Point {
        [self.u, self.v, self.direction].map(|axis| dot(axis, vector))
    }

    fn bounds(self, sweep: &Sweep) -> [Point; 2] {
        [-1., 1.].map(|sign| {
            std::array::from_fn(|i| {
                self.origin[i]
                    + sign
                        * sweep.support([
                            sign * self.u[i],
                            sign * self.v[i],
                            sign * self.direction[i],
                        ])
            })
        })
    }
}

struct Builder {
    vertices: Vec<Point>,
    edges: Vec<Edge>,
    by_ends: BTreeMap<[usize; 2], usize>,
    faces: Vec<Face>,
    tolerance: f64,
    split_long_arcs: bool,
    /// While a flange sweeps: its root ends and the welded level range.
    root_split: Option<RootSplit>,
}

/// A flange root welded to the web over part of its span. The rest of the
/// root stays an exposed planar face, for example beside a corner relief.
struct RootSplit {
    points: [[f64; 2]; 2],
    welded: [f64; 2],
}

impl RootSplit {
    fn exposed(&self, height: f64, tolerance: f64) -> Vec<[f64; 2]> {
        [[0., self.welded[0]], [self.welded[1], height]]
            .into_iter()
            .filter(|[low, high]| high - low > tolerance)
            .collect()
    }
}

impl Builder {
    fn new(tolerance: f64) -> Self {
        Self {
            vertices: Vec::new(),
            edges: Vec::new(),
            by_ends: BTreeMap::new(),
            faces: Vec::new(),
            tolerance,
            split_long_arcs: false,
            root_split: None,
        }
    }

    fn vertex(&mut self, point: Point) -> usize {
        if let Some(index) = self.vertices.iter().position(|&candidate| {
            let offset = sub(candidate, point);
            dot(offset, offset) <= self.tolerance.powi(2)
        }) {
            return index;
        }
        let index = self.vertices.len();
        self.vertices.push(point);
        index
    }

    fn edge(&mut self, points: [Point; 2], curve: Curve) -> Result<OrientedEdge, Error> {
        let vertices = points.map(|point| self.vertex(point));
        let [a, b] = vertices;
        if a == b {
            return Err(Error::Invalid);
        }
        let mut key = vertices;
        key.sort_unstable();
        if let Some(&edge) = self.by_ends.get(&key) {
            let prior = self.edges[edge];
            if curve != Curve::Line || prior.curve != Curve::Line {
                return Err(Error::Invalid);
            }
            return Ok(OrientedEdge {
                edge,
                forward: prior.vertices == vertices,
            });
        }
        let edge = self.edges.len();
        self.by_ends.insert(key, edge);
        self.edges.push(Edge { vertices, curve });
        Ok(OrientedEdge {
            edge,
            forward: true,
        })
    }

    fn sweep(
        &mut self,
        rings: &[&[(ProfileEdge, usize)]],
        sweep: &Sweep,
        frame: SweepFrame,
        caps: [usize; 2],
        omit: impl Iterator<Item = bool>,
    ) -> Result<(), Error> {
        let mut omit = omit.chain(std::iter::repeat(false));
        let height = sweep.height_mm();
        let mut cap_loops = [Vec::new(), Vec::new()];
        for (ring, authored) in sweep.rings().iter().zip(rings) {
            let edges = ring.iter().zip(*authored).map(|(&edge, &(_, source))| {
                (edge, source, omit.next().expect("unbounded omission flags"))
            });
            let ends = self.sweep_ring(edges, frame, height)?;
            for (cap, ring) in cap_loops.iter_mut().zip(ends) {
                cap.push(ring);
            }
        }
        for (i, mut rings) in cap_loops.into_iter().enumerate() {
            let direction = if i == 0 { -1. } else { 1. };
            self.faces.push(Face {
                boundary: rings.remove(0),
                holes: rings,
                support: Support::Plane {
                    normal: scale(frame.direction, direction),
                },
                point_mm: frame.point(sweep.centroid(), i as f64 * height),
                area_mm2: sweep.area_mm2(),
                source: 0,
                source_face: caps[i],
            });
        }
        Ok(())
    }

    fn sweep_ring(
        &mut self,
        edges: impl Iterator<Item = (ProfileEdge, usize, bool)>,
        frame: SweepFrame,
        height: f64,
    ) -> Result<[Vec<OrientedEdge>; 2], Error> {
        let reversed = dot(cross(frame.u, frame.v), frame.direction) < 0.;
        let mut end_loops = [Vec::new(), Vec::new()];
        for (edge, source, omit) in edges {
            let parts = if self.split_long_arcs && matches!(edge, ProfileEdge::Arc { .. }) {
                quarter_parts(edge)
            } else {
                1
            };
            for arc in segments(edge) {
                for i in 0..parts {
                    let fraction = [i, i + 1].map(|j| j as f64 / parts as f64);
                    let [a, b] = fraction.map(|t| arc.point(t));
                    let lifts = [0., height].map(|level| {
                        let points = [a, b].map(|p| frame.point(p, level));
                        let curve = curve(arc, fraction, frame, level)?;
                        self.edge(points, curve)
                    });
                    let [bottom, top] = lifts;
                    let (bottom, top) = (bottom?, top?);
                    end_loops[0].push(bottom);
                    end_loops[1].push(top);
                    let split = self
                        .root_split
                        .as_ref()
                        .filter(|split| split.points.contains(&a) && split.points.contains(&b));
                    let spans = match (omit, split) {
                        (false, _) => vec![[0., height]],
                        (true, None) => vec![],
                        (true, Some(split)) => split.exposed(height, self.tolerance),
                    };
                    for [low, high] in spans {
                        let mut level = |level: f64, lift: OrientedEdge| {
                            if level == 0. || level == height {
                                return Ok(lift);
                            }
                            let points = [a, b].map(|p| frame.point(p, level));
                            self.edge(points, curve(arc, fraction, frame, level)?)
                        };
                        let (bottom, top) = (level(low, bottom)?, level(high, top)?);
                        let end = self.rail(frame, b, [low, high])?;
                        let start = self.rail(frame, a, [low, high])?;
                        let mut boundary = vec![bottom];
                        boundary.extend(end);
                        boundary.push(reverse(top));
                        boundary.extend(start.into_iter().rev().map(reverse));
                        if reversed {
                            reverse_loop(&mut boundary);
                            boundary.rotate_right(1);
                        }
                        self.faces.push(Face {
                            boundary,
                            holes: vec![],
                            support: side_support(arc, frame, reversed),
                            point_mm: frame.point(a, low),
                            area_mm2: arc.length_mm() * (high - low) / parts as f64,
                            source: 0,
                            source_face: source,
                        });
                    }
                }
            }
        }
        reverse_loop(&mut end_loops[usize::from(reversed)]);
        Ok(end_loops)
    }

    /// The straight rail through one section point between two levels, split
    /// where a flange root's weld to the web begins or ends.
    fn rail(
        &mut self,
        frame: SweepFrame,
        point: [f64; 2],
        [low, high]: [f64; 2],
    ) -> Result<Vec<OrientedEdge>, Error> {
        let mut levels = vec![low];
        if let Some(split) = &self.root_split
            && split.points.contains(&point)
        {
            levels.extend(
                split
                    .welded
                    .into_iter()
                    .filter(|&level| level > low + self.tolerance && level < high - self.tolerance),
            );
        }
        levels.push(high);
        levels
            .windows(2)
            .map(|pair| {
                self.edge(
                    [pair[0], pair[1]].map(|level| frame.point(point, level)),
                    Curve::Line,
                )
            })
            .collect()
    }

    fn pierce_flange(
        &mut self,
        section: &Sweep,
        openings: &[SheetFlangeOpening<'_>],
        flange: SweepFrame,
        first_wall: usize,
        thickness: f64,
        minimum: f64,
    ) -> Result<f64, Error> {
        // A flange has one section ring, with the last edge welded to the web.
        // The other section edges each emitted exactly one face in ordinal order.
        let root = section.rings()[0].len() - 1;
        if openings.iter().any(|hole| hole.walls.contains(&root)) {
            return Err(Error::Unsupported);
        }
        let specs = openings
            .iter()
            .map(|hole| Piercing {
                sides: hole.walls,
                origin: flange.coordinates(sub(hole.origin_mm, flange.origin)),
                u: flange.coordinates(hole.u),
                v: flange.coordinates(hole.v),
                through: flange.coordinates(hole.through_mm),
                outline: hole.outline.iter().map(|&(edge, _)| edge).collect(),
            })
            .collect();
        let pierced =
            PiercedSweep::new(section, specs, self.tolerance, minimum).map_err(sweep_error)?;
        for (checked, authored) in pierced.openings().iter().zip(openings) {
            let height = checked.thickness_mm();
            if (height - thickness).abs() > self.tolerance {
                return Err(Error::Unsupported);
            }
            self.pierce(
                checked,
                flange,
                checked.spec().sides.map(|side| first_wall + side),
                authored.outline.iter().map(|&(_, source)| source),
            )?;
        }
        Ok(pierced.volume_mm3())
    }

    fn pierce(
        &mut self,
        opening: &CertifiedPiercing,
        sweep: SweepFrame,
        mouths: [usize; 2],
        sources: impl DoubleEndedIterator<Item = usize> + ExactSizeIterator,
    ) -> Result<(), Error> {
        let spec = opening.spec();
        let height = opening.thickness_mm();
        let [u, v, axial] = spec.origin;
        let frame = SweepFrame {
            origin: sweep.point([u, v], axial),
            u: sweep.vector(spec.u),
            v: sweep.vector(spec.v),
            direction: sweep.vector(scale(spec.through, height.recip())),
        };
        // The inward wall of a hole uses the opposite winding to a
        // positive extrusion. The resulting end loops face into the cut.
        let edges = spec
            .outline
            .iter()
            .zip(sources)
            .rev()
            .map(|(&edge, source)| (edge.reversed(), source, false));
        let loops = self.sweep_ring(edges, frame, height)?;
        for (mouth, ring) in mouths.into_iter().zip(loops) {
            let wall = &mut self.faces[mouth];
            wall.holes.push(ring);
            wall.area_mm2 -= opening.area_mm2();
        }
        Ok(())
    }

    fn finish(self, bounds_mm: [Point; 2], volume_mm3: f64) -> Result<Shell, Error> {
        let mut uses = vec![[0; 2]; self.edges.len()];
        for face in &self.faces {
            for ring in std::iter::once(&face.boundary).chain(&face.holes) {
                for edge in ring {
                    uses[edge.edge][usize::from(edge.forward)] += 1;
                }
            }
        }
        if uses.iter().any(|&count| count != [1, 1]) {
            return Err(Error::Invalid);
        }
        Ok(Shell {
            vertices: self.vertices,
            edges: self.edges,
            faces: self.faces,
            bounds_mm,
            volume_mm3,
            solid_count: 1,
            shell_count: 1,
            band: None,
            chart_w: vec![],
            chart_edges: 0,
            matched_quarter_rims: None,
        })
    }
}

fn reverse(edge: OrientedEdge) -> OrientedEdge {
    OrientedEdge {
        forward: !edge.forward,
        ..edge
    }
}

fn reverse_loop(ring: &mut [OrientedEdge]) {
    ring.reverse();
    for edge in ring {
        *edge = reverse(*edge);
    }
}

fn quarter_parts(edge: ProfileEdge) -> usize {
    match edge {
        ProfileEdge::Circle { .. } => 4,
        ProfileEdge::Line { .. } => 1,
        ProfileEdge::Arc { .. } => {
            let sweep = edge.signed_sweep_rad().expect("arc sweep").abs();
            (((sweep - ANGULAR) / FRAC_PI_2).ceil() as usize).max(1)
        }
    }
}

fn segments(edge: ProfileEdge) -> impl Iterator<Item = ProfileEdge> {
    let parts = if let ProfileEdge::Circle {
        center,
        radius_mm: r,
        clockwise,
    } = edge
    {
        let points =
            [[r, 0.], [0., r], [-r, 0.], [0., -r]].map(|[x, y]| [center[0] + x, center[1] + y]);
        let step = if clockwise { 3 } else { 1 };
        std::array::from_fn(|i| {
            Some(ProfileEdge::Arc {
                center,
                a: points[(i * step) % 4],
                b: points[((i + 1) * step) % 4],
                clockwise,
            })
        })
    } else {
        [Some(edge), None, None, None]
    };
    parts.into_iter().flatten()
}

fn curve(
    edge: ProfileEdge,
    fraction: [f64; 2],
    frame: SweepFrame,
    level: f64,
) -> Result<Curve, Error> {
    let ProfileEdge::Arc { center, a, .. } = edge else {
        return Ok(Curve::Line);
    };
    let [u, v] = [a[0] - center[0], a[1] - center[1]];
    let radius = u.portable_hypot(v);
    let ellipse = Ellipse::new(
        frame.point(center, level),
        cross(frame.u, frame.v),
        frame.u,
        [radius, radius],
    )
    .map_err(|_| Error::Invalid)?;
    let sweep = edge.signed_sweep_rad().expect("arc sweep");
    // Keep the authored circular support. Refitting rounded subdivision
    // endpoints could shift the radius or phase at large coordinates.
    EllipseArc::new(
        ellipse,
        sweep.mul_add(fraction[0], v.portable_atan2(u)),
        sweep * (fraction[1] - fraction[0]),
    )
    .map(Curve::Ellipse)
    .map_err(|_| Error::Unsupported)
}

fn side_support(edge: ProfileEdge, frame: SweepFrame, reversed: bool) -> Support {
    if let ProfileEdge::Arc {
        center, clockwise, ..
    } = edge
    {
        return Support::Cylinder {
            axis_origin_mm: frame.point(center, 0.),
            axis_direction: frame.direction,
            radius_mm: edge.radius_mm().expect("arc radius"),
            outward: !clockwise,
        };
    }
    let delta = sub(frame.point(edge.end(), 0.), frame.point(edge.start(), 0.));
    let normal = cross(delta, frame.direction);
    let sign = if reversed { -1. } else { 1. };
    Support::Plane {
        normal: scale(normal, sign / dot(normal, normal).sqrt()),
    }
}

fn add(a: Point, b: Point) -> Point {
    std::array::from_fn(|i| a[i] + b[i])
}

fn scale(point: Point, factor: f64) -> Point {
    point.map(|value| value * factor)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::display::Tessellation;
    const TOL: f64 = 1e-7;
    const MIN: f64 = 1e-5;

    fn line(a: [f64; 2], b: [f64; 2]) -> ProfileEdge {
        ProfileEdge::Line { a, b }
    }

    #[test]
    fn circular_trim_outside_the_coordinate_envelope_refuses_the_sweep() {
        let frame = SweepFrame {
            origin: [0., 0., 1e12],
            u: [1., 0., 0.],
            v: [0., 1., 0.],
            direction: [0., 0., 1.],
        };
        let arc = ProfileEdge::Arc {
            center: [0.; 2],
            a: [1., 0.],
            b: [0., 1.],
            clockwise: false,
        };
        let mut builder = Builder::new(TOL);
        assert_eq!(
            builder.sweep_ring(std::iter::once((arc, 0, false)), frame, 1.),
            Err(Error::Invalid)
        );
    }

    #[test]
    fn sweep_ring_keeps_both_caps_when_a_join_consumes_the_side() {
        let arc = ProfileEdge::Arc {
            center: [0.; 2],
            a: [2., 0.],
            b: [-2., 0.],
            clockwise: false,
        };
        for direction in [-1., 1.] {
            let frame = SweepFrame {
                origin: [0.; 3],
                u: [1., 0., 0.],
                v: [0., 1., 0.],
                direction: [0., 0., direction],
            };
            for (split, omit, count) in [(false, false, 1), (true, false, 2), (true, true, 2)] {
                let mut builder = Builder::new(TOL);
                builder.split_long_arcs = split;
                let caps = builder
                    .sweep_ring(std::iter::once((arc, 17, omit)), frame, 1.)
                    .unwrap();
                assert_eq!(caps.map(|cap| cap.len()), [count; 2]);
                assert_eq!(builder.faces.len(), if omit { 0 } else { count });
                assert!(builder.faces.iter().all(|face| face.source_face == 17));
            }
        }
    }

    fn web() -> Vec<SheetWebEdge> {
        let corners = [
            [0., 0.],
            [2., 0.],
            [8., 0.],
            [10., 0.],
            [10., 8.],
            [8., 8.],
            [2., 8.],
            [0., 8.],
        ];
        (0..corners.len())
            .map(|i| SheetWebEdge {
                ends: [corners[i], corners[(i + 1) % corners.len()]],
                source: i,
                flange: match i {
                    1 => Some(0),
                    5 => Some(1),
                    _ => None,
                },
            })
            .collect()
    }

    fn section() -> Vec<(ProfileEdge, usize)> {
        let arc = |a, b, clockwise| ProfileEdge::Arc {
            center: [0., -1.],
            a,
            b,
            clockwise,
        };
        [
            arc([0., 0.], [1., -1.], true),
            line([1., -1.], [1., -3.]),
            line([1., -3.], [2., -3.]),
            line([2., -3.], [2., -1.]),
            arc([2., -1.], [0., 1.], false),
            line([0., 1.], [0., 0.]),
        ]
        .into_iter()
        .enumerate()
        .map(|(i, edge)| (edge, i + 10))
        .collect()
    }

    #[test]
    fn opposite_partial_bends_share_edges_and_keep_analytic_holes() {
        let web = web();
        let section = section();
        let flanges = [[20, 21], [22, 23]].map(|ends| SheetFlange {
            section: &section,
            ends,
        });
        let hole = ProfileEdge::Circle {
            center: [5., 4.],
            radius_mm: 1.,
            clockwise: true,
        };
        for thickness in [-1., 1.] {
            let shell = formed_sheet(
                &web,
                &[vec![(hole, 50)]],
                [30, 31],
                &flanges,
                thickness,
                TOL,
                MIN,
            )
            .unwrap();
            assert_eq!(shell.solid_count(), 1);
            assert_eq!(shell.shell_count(), 1);
            assert_eq!(shell.faces().len(), 26);
            assert_eq!(shell.bounds_mm()[0][1], -2.);
            let volume = 80. - std::f64::consts::PI + 12. * (2. + std::f64::consts::PI * 0.75);
            assert!((shell.volume_mm3() - volume).abs() < TOL);
            let mesh = shell
                .mesh(Tessellation {
                    deflection: 0.01,
                    angle: 0.1,
                    triangle_limit: 2000,
                })
                .unwrap();
            let mut incidence = BTreeMap::<_, [usize; 2]>::new();
            for t in mesh.indices.as_chunks::<3>().0 {
                for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                    incidence.entry([a.min(b), a.max(b)]).or_default()[usize::from(a < b)] += 1;
                }
            }
            assert!(incidence.values().all(|counts| *counts == [1, 1]));
        }
    }

    fn rectangle() -> [(ProfileEdge, usize); 4] {
        let corners = [[0., 0.], [1.5, 0.], [1.5, 0.6], [0., 0.6]];
        std::array::from_fn(|i| (line(corners[i], corners[(i + 1) % 4]), 60 + i))
    }

    fn aperture<'a>(
        origin_mm: Point,
        up: f64,
        outline: &'a [(ProfileEdge, usize)],
    ) -> SheetFlangeOpening<'a> {
        SheetFlangeOpening {
            walls: [1, 3],
            origin_mm,
            u: [1., 0., 0.],
            v: [0., 0., up],
            through_mm: [0., -1., 0.],
            outline,
        }
    }

    fn aperture_fixture(
        section: &[(ProfileEdge, usize)],
        holes: &[&[SheetFlangeOpening<'_>]],
        signed_thickness: f64,
    ) -> Result<Shell, Error> {
        let flanges = [[20, 21], [22, 23]].map(|ends| SheetFlange { section, ends });
        let web_hole = ProfileEdge::Circle {
            center: [5., 4.],
            radius_mm: 1.,
            clockwise: true,
        };
        formed_sheet_with_openings(
            &web(),
            &[vec![(web_hole, 50)]],
            [30, 31],
            &flanges,
            holes,
            signed_thickness,
            TOL,
            MIN,
        )
    }

    #[test]
    fn rectangular_and_round_apertures_share_a_flange_without_changing_the_bend() {
        let section = section();
        let rectangle = rectangle();
        let circle = [(
            ProfileEdge::Circle {
                center: [0.; 2],
                radius_mm: 0.3,
                clockwise: false,
            },
            70,
        )];
        for up in [-1., 1.] {
            let baseline = aperture_fixture(&section, &[], up).unwrap();
            let rect = aperture([3., -1., -2.2 * up], up, &rectangle);
            let mut round = aperture([6.5, -2., -2. * up], up, &circle);
            round.walls.reverse();
            round.through_mm[1] *= -1.;
            let shell = aperture_fixture(&section, &[&[rect, round], &[]], up).unwrap();

            let rect_area = 1.5 * 0.6;
            let round_area = std::f64::consts::PI * 0.3_f64.powi(2);
            assert_eq!(shell.solid_count(), 1);
            assert_eq!(shell.bounds_mm(), baseline.bounds_mm());
            assert!(
                (baseline.volume_mm3() - shell.volume_mm3() - rect_area - round_area).abs() < TOL
            );
            let area_delta =
                2. * (1.5 + 0.6) + std::f64::consts::TAU * 0.3 - 2. * (rect_area + round_area);
            let area = |shell: &Shell| shell.faces().iter().map(|face| face.area_mm2).sum::<f64>();
            assert!((area(&shell) - area(&baseline) - area_delta).abs() < TOL);

            let rect_walls: Vec<_> = shell
                .faces()
                .iter()
                .filter(|face| (60..64).contains(&face.source_face))
                .collect();
            assert_eq!(rect_walls.len(), 4);
            assert!(
                rect_walls
                    .iter()
                    .all(|face| matches!(face.support, Support::Plane { .. }))
            );
            let circular_walls: Vec<_> = shell
                .faces()
                .iter()
                .filter(|face| face.source_face == 70)
                .collect();
            assert_eq!(circular_walls.len(), 4);
            assert!(
                circular_walls
                    .iter()
                    .all(|face| matches!(face.support, Support::Cylinder { outward: false, .. }))
            );
            for source in [11, 13] {
                assert_eq!(
                    shell
                        .faces()
                        .iter()
                        .find(|face| face.source_face == source)
                        .unwrap()
                        .holes
                        .len(),
                    2
                );
            }

            let mesh = shell
                .mesh(Tessellation {
                    deflection: 0.002,
                    angle: 0.1,
                    triangle_limit: 6000,
                })
                .unwrap();
            let proof = crate::manufacturing::certify(&mesh).unwrap();
            assert!((proof.signed_volume_mm3 - shell.volume_mm3()).abs() < 0.1);
        }
    }

    #[test]
    fn flange_aperture_refuses_bends_overlap_root_and_thickened_legs() {
        let section = section();
        let rectangle = rectangle();
        let valid = || aperture([3., -1., -2.2], 1., &rectangle);
        assert_eq!(
            aperture_fixture(&section, &[&[valid()]], 1.),
            Err(Error::Invalid)
        );
        let bend = aperture([3., -1., -1.4], 1., &rectangle);
        assert_eq!(
            aperture_fixture(&section, &[&[bend], &[]], 1.),
            Err(Error::Invalid)
        );
        assert_eq!(
            aperture_fixture(&section, &[&[valid(), valid()], &[]], 1.),
            Err(Error::Invalid)
        );
        let mut root = valid();
        root.walls = [1, 5];
        assert_eq!(
            aperture_fixture(&section, &[&[root], &[]], 1.),
            Err(Error::Unsupported)
        );

        let mut thick = section;
        thick.splice(
            2..4,
            [
                (line([1., -3.], [2.25, -3.]), 12),
                (line([2.25, -3.], [2.25, -1.5]), 13),
                (line([2.25, -1.5], [2., -1.5]), 15),
                (line([2., -1.5], [2., -1.]), 16),
            ],
        );
        assert!(aperture_fixture(&thick, &[], 1.).is_ok());
        let mut hole = valid();
        hole.through_mm[1] = -1.25;
        assert_eq!(
            aperture_fixture(&thick, &[&[hole], &[]], 1.),
            Err(Error::Unsupported)
        );
    }

    #[test]
    fn flange_aperture_cannot_reuse_an_existing_circular_trim() {
        let outline = section();
        let sweep = checked_sweep(&[&outline], 6., TOL, MIN).unwrap();
        let frame = SweepFrame {
            origin: [2., 0., 0.],
            u: [0., -1., 0.],
            v: [0., 0., 1.],
            direction: [1., 0., 0.],
        };
        let ring = [(
            ProfileEdge::Circle {
                center: [0.; 2],
                radius_mm: 0.3,
                clockwise: false,
            },
            70,
        )];
        let openings = [aperture([5., -1., -2.], 1., &ring)];
        let mut builder = Builder::new(TOL);
        builder
            .sweep(&[&outline], &sweep, frame, [20, 21], std::iter::empty())
            .unwrap();
        builder
            .pierce_flange(&sweep, &openings, frame, 0, 1., MIN)
            .unwrap();

        // The opening can certify in isolation but cannot be appended a
        // second time because its circular trim already belongs to the shell.
        let face_count = builder.faces.len();
        assert_eq!(
            builder.pierce_flange(&sweep, &openings, frame, 0, 1., MIN),
            Err(Error::Invalid),
        );
        assert_eq!(builder.faces.len(), face_count);
        assert_eq!(builder.faces[1].holes.len(), 1);
        assert_eq!(builder.faces[3].holes.len(), 1);
    }

    fn web_from_corners(corners: &[[f64; 2]], roots: &[(usize, usize)]) -> Vec<SheetWebEdge> {
        corners
            .iter()
            .zip(corners.iter().cycle().skip(1))
            .enumerate()
            .map(|(i, (&a, &b))| SheetWebEdge {
                ends: [a, b],
                source: i,
                flange: roots
                    .iter()
                    .find_map(|&(edge, flange)| (edge == i).then_some(flange)),
            })
            .collect()
    }

    fn upward_section() -> Vec<(ProfileEdge, usize)> {
        let arc = |a, b, clockwise| ProfileEdge::Arc {
            center: [0., 2.],
            a,
            b,
            clockwise,
        };
        [
            arc([0., 0.], [2., 2.], false),
            line([2., 2.], [2., 4.]),
            line([2., 4.], [1., 4.]),
            line([1., 4.], [1., 2.]),
            arc([1., 2.], [0., 1.], true),
            line([0., 1.], [0., 0.]),
        ]
        .into_iter()
        .enumerate()
        .map(|(i, edge)| (edge, i + 40))
        .collect()
    }

    fn four_flange_web(inset_long_sides: bool, relief_depth: f64) -> Vec<SheetWebEdge> {
        let corners = [[0., 0.], [10., 0.], [10., 8.], [0., 8.]];
        let mut web = Vec::new();
        for (flange, (&a, &b)) in corners
            .iter()
            .zip(corners.iter().cycle().skip(1))
            .enumerate()
        {
            let [dx, dy]: [f64; 2] = std::array::from_fn(|i| b[i] - a[i]);
            let length = dx.portable_hypot(dy);
            let [u, v] = [dx / length, dy / length];
            let place = |[s, q]: [f64; 2]| [a[0] + u * s + v * q, a[1] + v * s - u * q];
            let inset = flange % 2 == 1 || inset_long_sides;
            let contour = if inset {
                vec![
                    [0., 0.],
                    [0.5, 0.],
                    [0.5, -relief_depth],
                    [1.5, -relief_depth],
                    [1.5, 0.],
                    [length - 1.5, 0.],
                    [length - 1.5, -relief_depth],
                    [length - 0.5, -relief_depth],
                    [length - 0.5, 0.],
                    [length, 0.],
                ]
            } else {
                vec![[0., 0.], [length, 0.]]
            };
            let root = if inset { 4 } else { 0 };
            for (i, pair) in contour.windows(2).enumerate() {
                web.push(SheetWebEdge {
                    ends: [place(pair[0]), place(pair[1])],
                    source: web.len(),
                    flange: (i == root).then_some(flange),
                });
            }
        }
        web
    }

    #[test]
    fn four_adjacent_flange_pairs_keep_open_corners_and_opposite_bend_directions() {
        let down = section();
        let up = upward_section();
        let flanges = [&down, &up, &down, &up]
            .into_iter()
            .enumerate()
            .map(|(i, section)| SheetFlange {
                section,
                ends: [2 * i + 60, 2 * i + 61],
            })
            .collect::<Vec<_>>();
        for (inset_long_sides, edges, relieved_area, span) in
            [(false, 20, 1., 30.), (true, 36, 2., 24.)]
        {
            let web = four_flange_web(inset_long_sides, 0.25);
            assert_eq!(web.len(), edges);
            let volume = 80. - relieved_area + span * (2. + 3. * std::f64::consts::FRAC_PI_4);
            for (thickness, z) in [(1., [-3., 4.]), (-1., [-4., 3.])] {
                let shell =
                    formed_sheet(&web, &[], [80, 81], &flanges, thickness, TOL, MIN).unwrap();
                assert_eq!(shell.solid_count(), 1);
                assert_eq!(shell.bounds_mm(), [[-2., -2., z[0]], [12., 10., z[1]]]);
                assert!((shell.volume_mm3() - volume).abs() < TOL);
                assert_eq!(
                    shell
                        .faces()
                        .iter()
                        .filter(|face| matches!(face.support, Support::Cylinder { .. }))
                        .count(),
                    8,
                );
                let mesh = shell
                    .mesh(Tessellation {
                        deflection: 0.001,
                        angle: 0.07,
                        triangle_limit: 6000,
                    })
                    .unwrap();
                let proof = crate::manufacturing::certify(&mesh).unwrap();
                assert!((proof.signed_volume_mm3 - volume).abs() < 0.1);
            }
        }
        // At 0.5 mm depth, neighboring reliefs touch at an inner corner.
        let touching_reliefs = four_flange_web(true, 0.5);
        assert_eq!(
            formed_sheet(&touching_reliefs, &[], [80, 81], &flanges, 1., TOL, MIN),
            Err(Error::Invalid),
        );
    }

    #[test]
    fn adjacent_root_contact_and_relief_eroding_a_neighbor_fail_closed() {
        let section = section();
        let flanges = [[60, 61], [62, 63]].map(|ends| SheetFlange {
            section: &section,
            ends,
        });
        let touching = web_from_corners(
            &[[0., 0.], [10., 0.], [10., 8.], [0., 8.]],
            &[(0, 0), (1, 1)],
        );
        assert_eq!(
            formed_sheet(&touching, &[], [80, 81], &flanges, 1., TOL, MIN),
            Err(Error::Unsupported),
        );
        // A relief flush with the corner overlaps the previous root's
        // outline. It cannot be spliced into a full-root welded shell.
        let eroded = web_from_corners(
            &[
                [0., 0.],
                [10., 0.],
                [9.5, 0.],
                [9.5, 1.5],
                [10., 1.5],
                [10., 8.],
                [0., 8.],
            ],
            &[(0, 0), (4, 1)],
        );
        assert_eq!(
            formed_sheet(&eroded, &[], [80, 81], &flanges, 1., TOL, MIN),
            Err(Error::Invalid),
        );
    }

    fn hem(slanted: bool) -> Vec<(ProfileEdge, usize)> {
        let arc = |a, b, clockwise| ProfileEdge::Arc {
            center: [0., -0.5],
            a,
            b,
            clockwise,
        };
        // A 180° hem with a 1 mm gap folds its 3 mm leg back under the web.
        // The slanted variant lifts the leg's far end into the web slab.
        let tip = if slanted { [-3., 0.5] } else { [-3., -1.] };
        [
            arc([0., 0.], [0., -1.], true),
            line([0., -1.], tip),
            line(tip, [-3., -2.]),
            line([-3., -2.], [0., -2.]),
            arc([0., -2.], [0., 1.], false),
            line([0., 1.], [0., 0.]),
        ]
        .into_iter()
        .enumerate()
        .map(|(i, edge)| (edge, i + 70))
        .collect()
    }

    fn section_area(section: &[(ProfileEdge, usize)]) -> f64 {
        checked_sweep(&[section], 1., TOL, MIN).unwrap().area_mm2()
    }

    fn certified_volume(shell: &Shell) -> f64 {
        let mesh = shell
            .mesh(Tessellation {
                deflection: 0.001,
                angle: 0.07,
                triangle_limit: 20_000,
            })
            .unwrap();
        crate::manufacturing::certify(&mesh)
            .unwrap()
            .signed_volume_mm3
    }

    #[test]
    fn a_hem_folded_under_the_web_is_separated_by_its_gap() {
        let section = hem(false);
        let flanges = [SheetFlange {
            section: &section,
            ends: [60, 61],
        }];
        let web = web_from_corners(&[[0., 0.], [10., 0.], [10., 8.], [0., 8.]], &[(0, 0)]);
        let shell = formed_sheet(&web, &[], [80, 81], &flanges, 1., TOL, MIN).unwrap();
        let volume = 80. + 10. * (std::f64::consts::PI + 3.);
        assert!((section_area(&section) - (std::f64::consts::PI + 3.)).abs() < 1e-9);
        assert!((shell.volume_mm3() - volume).abs() < 1e-9);
        assert_eq!(shell.bounds_mm(), [[0., -1.5, -2.], [10., 8., 1.]]);
        assert!((certified_volume(&shell) - volume).abs() < 0.05);
        // A leg rising back into the web slab is refused, not fused.
        let slanted = hem(true);
        let flanges = [SheetFlange {
            section: &slanted,
            ends: [60, 61],
        }];
        assert_eq!(
            formed_sheet(&web, &[], [80, 81], &flanges, 1., TOL, MIN),
            Err(Error::Unsupported)
        );
    }

    #[test]
    fn flanges_bend_into_a_notch_only_when_they_clear_its_walls_and_each_other() {
        let section = section();
        let flanges = [[60, 61], [62, 63]].map(|ends| SheetFlange {
            section: &section,
            ends,
        });
        // A 6 mm notch: each 2 mm-deep bend leaves 2 mm between the legs.
        // The roots stop 1 mm short of the notch floor.
        let wide = web_from_corners(
            &[
                [0., 0.],
                [12., 0.],
                [12., 10.],
                [9., 10.],
                [9., 4.],
                [9., 3.],
                [3., 3.],
                [3., 4.],
                [3., 10.],
                [0., 10.],
            ],
            &[(3, 0), (7, 1)],
        );
        let shell = formed_sheet(&wide, &[], [80, 81], &flanges, 1., TOL, MIN).unwrap();
        let volume = 78. + 2. * 6. * section_area(&section);
        assert!((shell.volume_mm3() - volume).abs() < 1e-9);
        assert!((certified_volume(&shell) - volume).abs() < 0.05);
        // Legs that meet in a 4 mm notch, or a root running into the notch
        // floor, are refused.
        let narrow = web_from_corners(
            &[
                [0., 0.],
                [12., 0.],
                [12., 10.],
                [8., 10.],
                [8., 4.],
                [8., 3.],
                [4., 3.],
                [4., 4.],
                [4., 10.],
                [0., 10.],
            ],
            &[(3, 0), (7, 1)],
        );
        let floor = web_from_corners(
            &[
                [0., 0.],
                [12., 0.],
                [12., 10.],
                [9., 10.],
                [9., 3.],
                [3., 3.],
                [3., 4.],
                [3., 10.],
                [0., 10.],
            ],
            &[(3, 0), (6, 1)],
        );
        for web in [narrow, floor] {
            assert_eq!(
                formed_sheet(&web, &[], [80, 81], &flanges, 1., TOL, MIN),
                Err(Error::Unsupported)
            );
        }
    }

    #[test]
    fn a_root_past_a_flush_relief_keeps_an_exposed_face() {
        let section = section();
        let flanges = [[60, 61], [62, 63]].map(|ends| SheetFlange {
            section: &section,
            ends,
        });
        // The relief of the flange on x = 10 runs 1 mm along the edge-0 root.
        let web = web_from_corners(
            &[
                [0., 0.],
                [9., 0.],
                [9., 1.5],
                [10., 1.5],
                [10., 8.],
                [0., 8.],
            ],
            &[(0, 0), (3, 1)],
        );
        let overhangs = [[0., 1.], [0., 0.]];
        let shell = formed_sheet_with_overhangs(
            &web,
            &[],
            [80, 81],
            &flanges,
            &[],
            &overhangs,
            1.,
            TOL,
            MIN,
        )
        .unwrap();
        let volume = 78.5 + 16.5 * section_area(&section);
        assert!((shell.volume_mm3() - volume).abs() < 1e-9);
        assert_eq!(shell.bounds_mm()[1][0], 12.);
        assert_eq!(shell.bounds_mm()[0][0], 0.);
        let exposed: Vec<_> = shell
            .faces()
            .iter()
            .filter(|face| face.source_face == 15)
            .collect();
        assert_eq!(exposed.len(), 1);
        assert!((exposed[0].area_mm2 - 1.).abs() < 1e-12);
        assert!(matches!(exposed[0].support, Support::Plane { normal } if normal == [0., 1., 0.]));
        assert!((certified_volume(&shell) - volume).abs() < 0.05);
        // Without the relief the web still covers that root: refused.
        let covered = web_from_corners(
            &[
                [0., 0.],
                [9., 0.],
                [10., 0.],
                [10., 1.5],
                [10., 8.],
                [0., 8.],
            ],
            &[(0, 0), (3, 1)],
        );
        assert_eq!(
            formed_sheet_with_overhangs(
                &covered,
                &[],
                [80, 81],
                &flanges,
                &[],
                &overhangs,
                1.,
                TOL,
                MIN
            ),
            Err(Error::Unsupported)
        );
        assert_eq!(
            formed_sheet_with_overhangs(
                &web,
                &[],
                [80, 81],
                &flanges,
                &[],
                &[[0., -1.], [0., 0.]],
                1.,
                TOL,
                MIN
            ),
            Err(Error::Invalid)
        );
    }

    #[test]
    fn malformed_spans_and_sections_are_refused() {
        let web = web();
        let section = section();
        let flange = || SheetFlange {
            section: &section,
            ends: [20, 21],
        };
        let build = |web: &[SheetWebEdge], flanges: &[SheetFlange<'_>]| {
            formed_sheet(web, &[], [30, 31], flanges, 1., TOL, MIN)
        };
        assert_eq!(build(&web, &[]), Err(Error::Unsupported));
        assert_eq!(build(&web, &[flange()]), Err(Error::Invalid)); // Unresolved hinge.
        let mut absent = web.clone();
        absent[1].flange = None;
        assert_eq!(build(&absent, &[flange(), flange()]), Err(Error::Invalid));
        let mut duplicate = web.clone();
        duplicate[0].flange = Some(0);
        assert_eq!(
            build(&duplicate, &[flange(), flange()]),
            Err(Error::Invalid)
        );
        let mut adjacent = web.clone();
        adjacent[5].flange = None;
        adjacent[3].flange = Some(1);
        assert!(build(&adjacent, &[flange(), flange()]).is_ok());
        let mut diagonal = web.clone();
        diagonal[1].ends[0][1] = 1.;
        diagonal[0].ends[1][1] = 1.;
        assert_eq!(
            build(&diagonal, &[flange(), flange()]),
            Err(Error::Unsupported)
        );

        let mut one = web;
        one[5].flange = None;
        let wrong_root = &section[..5];
        assert_eq!(
            build(
                &one,
                &[SheetFlange {
                    section: wrong_root,
                    ends: [20, 21]
                }]
            ),
            Err(Error::Invalid)
        );
        let mut rotated = section;
        rotated.rotate_left(1);
        assert_eq!(
            build(
                &one,
                &[SheetFlange {
                    section: &rotated,
                    ends: [20, 21]
                }]
            ),
            Err(Error::Unsupported)
        );
        let mut huge = one;
        huge[1].ends[0][1] = 1e10;
        assert_eq!(
            build(
                &huge,
                &[SheetFlange {
                    section: &rotated,
                    ends: [20, 21]
                }]
            ),
            Err(Error::Limit)
        );
    }

    #[test]
    fn builder_rejects_open_or_ambiguous_edge_incidence() {
        let mut builder = Builder::new(TOL);
        let [a, b] = [[0.; 3], [1., 0., 0.]];
        let edge = builder.edge([a, b], Curve::Line).unwrap();
        assert_eq!(builder.edge([b, a], Curve::Line), Ok(reverse(edge)));
        assert_eq!(builder.edge([a, a], Curve::Line), Err(Error::Invalid));
        let circle = Ellipse::new([0.5, 0., 0.], [0., 0., 1.], [1., 0., 0.], [0.5; 2]).unwrap();
        let arc = EllipseArc::new(circle, std::f64::consts::PI, std::f64::consts::PI).unwrap();
        assert_eq!(
            builder.edge([a, b], Curve::Ellipse(arc)),
            Err(Error::Invalid)
        );
        assert_eq!(builder.finish([[0.; 3], [1.; 3]], 1.), Err(Error::Invalid));
    }
}
