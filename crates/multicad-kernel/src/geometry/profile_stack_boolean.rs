//! Exact Booleans of two extrusions along one common direction.
//!
//! A box is the extrusion of a rectangle along each of its axes, and a
//! normal profile sweep is the extrusion of its line-and-arc profile. When
//! the sweep direction is parallel to the other operand's extrusion axis,
//! both are profiles over spans of one common axis. Their union, difference
//! or common part is then, slab by slab between the span ends, the planar
//! Boolean of the two profiles, and its boundary consists of:
//! - planar caps where the material changes between adjacent slabs;
//! - walls on the profile lines (planes) and arcs (cylinders) where the
//!   material differs across the curve.
//!
//! All slabs share one planar arrangement of both profiles, so every cap
//! and wall is bounded by the same exact pieces. Coplanar or co-cylindrical
//! contributions of both operands merge into one face, and contact faces
//! between them vanish. The boundary is assembled from exact lines,
//! circles, planes and cylinders; the retained solid validates its closure
//! and its integrated volume must reproduce the slab areas from Green's
//! theorem. Nothing is tessellated. Tangent or ambiguous contacts, results
//! with several pieces and touching faces are refused in plain words.
use super::authored_sweep::{AuthoredNaming, FaceName};
use super::*;
use multicad_core::references::{AnalyticGeometry, ReferenceFrame, TopologyKind};
use multicad_rust_kernel::brep::{
    self, Boundary, Curve, Edge as BrepEdge, Face as BrepFace, Frame, Loop, Shell, Surface,
};
use std::collections::{BTreeMap, BTreeSet};

mod arrangement;
use arrangement::{Arrangement, P2, Piece, Refusal, Segment, Shape, Support};

const OPERATION: &str = "profile stack Boolean";
type P3 = [f64; 3];

fn refused(message: impl Into<String>) -> GeometryError {
    error(OPERATION, GeometryErrorKind::Unsupported, message)
}

/// The common chart: profile coordinates along `u`, `v`, spans along `n`.
#[derive(Clone, Copy)]
struct Chart {
    frame: PlaneFrame,
    n: P3,
}

impl Chart {
    fn new(frame: PlaneFrame) -> Self {
        Self {
            frame,
            n: frame.normal(),
        }
    }

    fn world(&self, p: P2, z: f64) -> P3 {
        let o = self.frame.point(p);
        std::array::from_fn(|i| o[i] + z * self.n[i])
    }

    fn vector(&self, d: P2) -> P3 {
        std::array::from_fn(|i| d[0] * self.frame.u[i] + d[1] * self.frame.v[i])
    }

    fn local(&self, p: P3) -> P3 {
        let d = std::array::from_fn(|i| p[i] - self.frame.origin[i]);
        [dot(d, self.frame.u), dot(d, self.frame.v), dot(d, self.n)]
    }
}

/// One operand as a profile over a span of the chart's axis.
struct Extrusion {
    /// Closed boundaries with the material on their left.
    rings: Vec<Vec<Shape>>,
    /// The authored identity of each ring curve's side face.
    sides: Vec<Vec<BTreeSet<Origin>>>,
    span: [f64; 2],
    topology: TopologySnapshot,
}

impl Solid {
    /// A plain box as the extrusion of its cross-section along the axis
    /// parallel to `chart.n`.
    fn box_extrusion(&self, chart: &Chart, policy: GeometryPolicy) -> Option<Extrusion> {
        let (prism, frame) = self.rust_prism_frame()?;
        if prism.rounded_end().is_some() || prism.inner_loops().next().is_some() {
            return None;
        }
        let size = prism.dimensions().map(Rational::as_f64);
        let axes = [frame.u, frame.v, frame.normal()];
        let along = (0..3).find(|&i| (dot(axes[i], chart.n).abs() - 1.).abs() <= 1e-12)?;
        let [j, k] = [(along + 1) % 3, (along + 2) % 3];
        let corner = |s: f64, t: f64| -> P2 {
            let p: P3 = std::array::from_fn(|i| {
                frame.origin[i] + s * size[j] * axes[j][i] + t * size[k] * axes[k][i]
            });
            let l = chart.local(p);
            [l[0], l[1]]
        };
        let mut corners = [
            corner(0., 0.),
            corner(1., 0.),
            corner(1., 1.),
            corner(0., 1.),
        ];
        let area: f64 = (0..4)
            .map(|i| {
                let [a, b] = [corners[i], corners[(i + 1) % 4]];
                a[0] * b[1] - a[1] * b[0]
            })
            .sum();
        if area < 0. {
            corners.reverse();
        }
        let z0 = chart.local(frame.origin)[2];
        let z1 = z0 + size[along] * dot(axes[along], chart.n);
        let topology = self.rust_topology()?.ok()?;
        let ring: Vec<Shape> = (0..4)
            .map(|i| Shape::Line {
                a: corners[i],
                b: corners[(i + 1) % 4],
            })
            .collect();
        let sides = ring
            .iter()
            .map(|&shape| {
                let (point, outward) = wall_plane(chart, shape);
                plane_origins(&topology, point, outward, policy)
            })
            .collect();
        Some(Extrusion {
            rings: vec![ring],
            sides: vec![sides],
            span: [z0.min(z1), z0.max(z1)],
            topology,
        })
    }

    /// An unpierced normal sweep of lines and arcs in `chart`, whose own
    /// frame must share the chart's axis.
    fn sweep_extrusion(&self, chart: &Chart, policy: GeometryPolicy) -> Option<Extrusion> {
        let (body, frame) = self.normal_sweep_body()?;
        let (rings, span) = body.prism_profile()?;
        let n = frame.normal();
        let flip = dot(n, chart.n);
        if (flip.abs() - 1.).abs() > 1e-12 {
            return None;
        }
        let mut origins = body.side_origins();
        let tolerance = policy.linear_tolerance_mm;
        let to_chart = |p: P2| -> P2 {
            let l = chart.local(frame.point(p));
            [l[0], l[1]]
        };
        let offset = chart.local(frame.origin)[2];
        let mut shapes = vec![];
        let mut sides = vec![];
        for ring in &rings {
            let mut ring_sides = vec![];
            let mut ring_shapes = vec![];
            for edge in ring {
                let origin = origins.next()?.clone();
                for shape in arrangement::ring_shapes(std::slice::from_ref(edge), tolerance)? {
                    ring_shapes.push(map_shape(shape, to_chart, flip < 0.));
                    ring_sides.push(BTreeSet::from([origin.clone()]));
                }
            }
            // A mirrored chart reverses the orientation of every ring.
            if flip < 0. {
                ring_shapes.reverse();
                ring_sides.reverse();
                for shape in &mut ring_shapes {
                    *shape = reversed(*shape);
                }
            }
            shapes.push(ring_shapes);
            sides.push(ring_sides);
        }
        let z = span.map(|h| offset + flip * h);
        let topology = self.rust_topology()?.ok()?;
        // Planar sides also carry the names their faces publish.
        for (ring, names) in shapes.iter().zip(&mut sides) {
            for (&shape, names) in ring.iter().zip(names) {
                if matches!(shape, Shape::Line { .. }) {
                    let (point, outward) = wall_plane(chart, shape);
                    names.extend(plane_origins(&topology, point, outward, policy));
                }
            }
        }
        Some(Extrusion {
            rings: shapes,
            sides,
            span: [z[0].min(z[1]), z[0].max(z[1])],
            topology,
        })
    }

    /// A placed circular cylinder whose axis is parallel to the chart's.
    fn cylinder_extrusion(&self, chart: &Chart, policy: GeometryPolicy) -> Option<Extrusion> {
        let (cylinder, frame) = self.axial_cylinder_frame()?;
        let axis = frame.normal();
        let flip = dot(axis, chart.n);
        if (flip.abs() - 1.).abs() > 1e-12 {
            return None;
        }
        let radius = cylinder.circle().radius().as_f64();
        let height = cylinder.height().as_f64();
        let base = chart.local(frame.origin);
        // The frame starts at one rim; a reversed sketch sweeps against it.
        // Either reading must reproduce the solid's own bounds.
        let bounds = self.bounds().ok()?;
        let tolerance = policy.linear_tolerance_mm;
        let span = [1., -1.].into_iter().find_map(|direction| {
            let z = [base[2], base[2] + direction * flip * height];
            let ends: [P3; 2] = [0., 1.].map(|t| {
                std::array::from_fn(|k| frame.origin[k] + t * direction * height * axis[k])
            });
            let fits = (0..3).all(|k| {
                let reach = radius * (1. - axis[k] * axis[k]).max(0.).sqrt();
                let [low, high] = [ends[0][k].min(ends[1][k]), ends[0][k].max(ends[1][k])];
                (low - reach - bounds.min[k]).abs() <= 16. * tolerance
                    && (high + reach - bounds.max[k]).abs() <= 16. * tolerance
            });
            fits.then_some([z[0].min(z[1]), z[0].max(z[1])])
        })?;
        let origin = self.cylinder_sketch().map_or_else(
            || Origin {
                source: EntityId::from_uuid_bits(cylinder.feature()),
                role: "cylinder.wall".into(),
            },
            |sketch| sketch.origin.clone(),
        );
        let ring = arrangement::ring_shapes(
            &[multicad_rust_kernel::normal_sweep::Edge::Circle {
                center: [base[0], base[1]],
                radius_mm: radius,
                clockwise: false,
            }],
            tolerance,
        )?;
        let sides = vec![BTreeSet::from([origin]); ring.len()];
        Some(Extrusion {
            rings: vec![ring],
            sides: vec![sides],
            span,
            topology: self.rust_topology()?.ok()?,
        })
    }

    fn profile_extrusion(&self, chart: &Chart, policy: GeometryPolicy) -> Option<Extrusion> {
        self.sweep_extrusion(chart, policy)
            .or_else(|| self.cylinder_extrusion(chart, policy))
            .or_else(|| self.box_extrusion(chart, policy))
    }

    /// Union, difference or common part of a box, normal sweep or circular
    /// cylinder with a normal sweep or cylinder along a parallel axis.
    pub(super) fn profile_stack_boolean(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        let (operands, chart) = self.profile_operands(other, policy)?;
        let feature = self.rust_feature()?;
        Some(Stack::new(operands, operation, chart, policy).and_then(|stack| stack.solid(feature)))
    }

    /// The exact common volume of two parallel profile extrusions, from the
    /// slab sections of their common part, or `None` outside this provider.
    pub(crate) fn rust_profile_stack_overlap_volume(
        &self,
        other: &Self,
        policy: GeometryPolicy,
    ) -> Option<Result<f64, GeometryError>> {
        let (operands, chart) = self.profile_operands(other, policy)?;
        Some(
            Stack::new(operands, BooleanOperation::Intersect, chart, policy)
                .map(|stack| stack.volume()),
        )
    }

    fn profile_operands(
        &self,
        other: &Self,
        policy: GeometryPolicy,
    ) -> Option<([Extrusion; 2], Chart)> {
        if policy.boolean_fuzzy_mm != 0. || policy != self.policy || policy != other.policy {
            return None;
        }
        // Pairs of placed cylinders stay with the cylinder providers.
        if self.axial_cylinder_frame().is_some() && other.axial_cylinder_frame().is_some() {
            return None;
        }
        // The chart of the curved operand: a sweep's own frame, or a
        // cylinder's frame.
        let frame = [self, other].into_iter().find_map(|solid| {
            solid
                .normal_sweep_body()
                .map(|(_, frame)| frame)
                .or_else(|| solid.axial_cylinder_frame().map(|(_, frame)| frame))
        })?;
        let chart = Chart::new(frame);
        let operands = [
            self.profile_extrusion(&chart, policy)?,
            other.profile_extrusion(&chart, policy)?,
        ];
        Some((operands, chart))
    }
}

fn reversed(shape: Shape) -> Shape {
    match shape {
        Shape::Line { a, b } => Shape::Line { a: b, b: a },
        Shape::Arc {
            center,
            radius,
            start,
            sweep,
        } => Shape::Arc {
            center,
            radius,
            start: start + sweep,
            sweep: -sweep,
        },
    }
}

/// A curve of the sweep's own frame in chart coordinates. Charts share the
/// axis; a reversed axis mirrors the profile, reversing angles.
fn map_shape(shape: Shape, to_chart: impl Fn(P2) -> P2, mirrored: bool) -> Shape {
    match shape {
        Shape::Line { a, b } => Shape::Line {
            a: to_chart(a),
            b: to_chart(b),
        },
        Shape::Arc {
            center,
            radius,
            sweep,
            ..
        } => {
            // The start angle is measured again in the chart.
            let c = to_chart(center);
            let p = to_chart(shape.point(0.));
            let d = arrangement::sub(p, c);
            let angle = d[1].portable_atan2(d[0]);
            Shape::Arc {
                center: c,
                radius,
                start: angle,
                sweep: if mirrored { -sweep } else { sweep },
            }
        }
    }
}

/// A point on the wall of a profile line and its outward normal (to the
/// right of the line, away from the material).
fn wall_plane(chart: &Chart, shape: Shape) -> (P3, P3) {
    let [a, b] = shape.ends();
    let d = arrangement::sub(b, a);
    let length = d[0].portable_hypot(d[1]);
    (
        chart.world(a, 0.),
        chart.vector([d[1] / length, -d[0] / length]),
    )
}

/// The origins of every planar face of `topology` on the plane through
/// `point` with outward `normal`.
fn plane_origins(
    topology: &TopologySnapshot,
    point: P3,
    normal: P3,
    policy: GeometryPolicy,
) -> BTreeSet<Origin> {
    topology
        .elements
        .iter()
        .filter(|e| e.key.kind == TopologyKind::Face)
        .filter(|e| {
            matches!(e.geometry, Some(AnalyticGeometry::Plane { frame })
                if dot(frame.normal(), normal) > 1. - 1e-10
                    && dot(normal, std::array::from_fn(|i| point[i] - frame.origin_mm[i])).abs()
                        <= policy.linear_tolerance_mm)
        })
        .flat_map(|e| e.origins.iter().cloned())
        .collect()
}

/// Key of one boundary edge before simplification.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
enum EdgeKey {
    /// Piece at a level.
    Across(usize, usize),
    /// Node over a slab, running up the axis.
    Up(usize, usize),
}

/// A directed use of an edge key.
type Use = (EdgeKey, bool);

/// A boundary face before assembly.
struct Patch {
    surface: Surface,
    same_sense: bool,
    /// Closed loops, the outer loop of a planar cap first.
    loops: Vec<Vec<Use>>,
    origins: BTreeSet<Origin>,
    support: Option<ReferenceFrame>,
}

struct Stack {
    chart: Chart,
    arrangement: Arrangement,
    levels: Vec<f64>,
    /// Slab range `[first, end)` each operand occupies.
    spans: [[usize; 2]; 2],
    operation: BooleanOperation,
    operands: [Extrusion; 2],
    /// The (operand, ring, curve) each segment came from.
    owners: Vec<[usize; 3]>,
    policy: GeometryPolicy,
}

impl Stack {
    fn new(
        operands: [Extrusion; 2],
        operation: BooleanOperation,
        chart: Chart,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        let tolerance = policy.linear_tolerance_mm;
        let mut segments = vec![];
        let mut owners = vec![];
        for (operand, extrusion) in operands.iter().enumerate() {
            for (r, ring) in extrusion.rings.iter().enumerate() {
                for (c, &shape) in ring.iter().enumerate() {
                    segments.push(Segment { shape, operand });
                    owners.push([operand, r, c]);
                }
            }
        }
        let arrangement = Arrangement::new(&segments, tolerance).map_err(|cause| match cause {
            Refusal::Contact => refused(
                "The sketch region touches the body's outline tangentially or within the \
                 modelling tolerance; move it clear or let it cross cleanly",
            ),
            Refusal::Limit => error(
                OPERATION,
                GeometryErrorKind::Limit,
                "The two profiles have too many curves for an exact profile Boolean",
            ),
        })?;
        // Span ends coincide only up to floating-point noise; a real gap
        // or penetration below the minimum feature is refused below, never
        // snapped into contact or a through cut.
        let magnitude = operands
            .iter()
            .flat_map(|o| o.span)
            .map(f64::abs)
            .fold(1., f64::max);
        let numeric = 64. * f64::EPSILON * magnitude;
        let mut levels: Vec<f64> = vec![];
        for z in operands.iter().flat_map(|o| o.span) {
            if !levels.iter().any(|&l| (l - z).abs() <= numeric) {
                levels.push(z);
            }
        }
        levels.sort_by(f64::total_cmp);
        if levels
            .windows(2)
            .any(|w| w[1] - w[0] < policy.minimum_feature_mm)
        {
            return Err(refused(
                "The extrusion ends closer than the minimum feature size to a face of the body; \
                 end it exactly on that face or clearly away from it",
            ));
        }
        let level = |z: f64| {
            levels
                .iter()
                .position(|&l| (l - z).abs() <= numeric)
                .expect("span end is a level")
        };
        let spans = [0, 1].map(|i| operands[i].span.map(level));
        Ok(Self {
            chart,
            arrangement,
            levels,
            spans,
            operation,
            operands,
            owners,
            policy,
        })
    }

    fn slabs(&self) -> usize {
        self.levels.len() - 1
    }

    fn active(&self, operand: usize, slab: usize) -> bool {
        let [first, end] = self.spans[operand];
        (first..end).contains(&slab)
    }

    /// Each operand's material on `side` (0 left, 1 right) of a piece in a
    /// slab; outside every slab nothing is material.
    fn inputs(&self, piece: &Piece, slab: Option<usize>, side: usize) -> [bool; 2] {
        [0, 1].map(|o| slab.is_some_and(|k| self.active(o, k)) && piece.sides[o][side])
    }

    fn result(&self, [a, b]: [bool; 2]) -> bool {
        match self.operation {
            BooleanOperation::Union => a || b,
            BooleanOperation::Subtract => a && !b,
            BooleanOperation::Intersect => a && b,
        }
    }

    fn material(&self, piece: &Piece, slab: Option<usize>, side: usize) -> bool {
        self.result(self.inputs(piece, slab, side))
    }

    /// The slabs below and above a level.
    fn around(&self, level: usize) -> [Option<usize>; 2] {
        [
            level.checked_sub(1),
            (level < self.slabs()).then_some(level),
        ]
    }

    /// Whether the result's material differs across `piece` in `slab`, and
    /// if so whether it lies on the left.
    fn wall(&self, piece: &Piece, slab: usize) -> Option<bool> {
        let [left, right] = [0, 1].map(|s| self.material(piece, Some(slab), s));
        (left != right).then_some(left)
    }

    /// Exact area of the result's section in `slab`, by Green's theorem.
    fn section_area(&self, slab: usize) -> f64 {
        0.5 * self
            .arrangement
            .pieces
            .iter()
            .filter_map(|p| {
                self.wall(p, slab)
                    .map(|left| if left { 1. } else { -1. } * p.shape.twice_area())
            })
            .sum::<f64>()
    }

    /// The faces of the result before edge simplification.
    fn patches(&self) -> Result<Vec<Patch>, GeometryError> {
        let mut patches = self.caps()?;
        patches.extend(self.walls()?);
        Ok(patches)
    }

    /// Planar faces at each level where the material changes between the
    /// slabs below and above.
    fn caps(&self) -> Result<Vec<Patch>, GeometryError> {
        let mut patches = vec![];
        for level in 0..self.levels.len() {
            let around = self.around(level);
            for up in [true, false] {
                // Whether `side` of `piece` lies under this cap face.
                let covered = |piece: &Piece, side: usize| {
                    let [below, above] = around.map(|slab| self.material(piece, slab, side));
                    if up { below && !above } else { above && !below }
                };
                let mut uses = vec![];
                for (index, piece) in self.arrangement.pieces.iter().enumerate() {
                    let [left, right] = [0, 1].map(|s| covered(piece, s));
                    if left != right {
                        // Counter-clockwise about the outward normal.
                        uses.push((EdgeKey::Across(index, level), left == up));
                    }
                }
                if uses.is_empty() {
                    continue;
                }
                let loops = self.trace(&uses)?;
                for face in self.nest(loops, up)? {
                    let z = self.levels[level];
                    let frame = self.chart.frame;
                    let origin = self.chart.world([0., 0.], z);
                    let support = ReferenceFrame {
                        origin_mm: origin,
                        u: frame.u,
                        v: if up { frame.v } else { frame.v.map(|x| -x) },
                    };
                    let mut origins = BTreeSet::new();
                    // Every operand cap that bounds part of this face.
                    for (index, piece) in self.arrangement.pieces.iter().enumerate() {
                        for side in 0..2 {
                            if !covered(piece, side) || !self.on_face(&face, index, side, up) {
                                continue;
                            }
                            let [below, above] = around.map(|slab| self.inputs(piece, slab, side));
                            for o in 0..2 {
                                if below[o] != above[o] {
                                    let outward = if below[o] { 1. } else { -1. };
                                    let normal = self.chart.n.map(|x| x * outward);
                                    origins.extend(plane_origins(
                                        &self.operands[o].topology,
                                        origin,
                                        normal,
                                        self.policy,
                                    ));
                                }
                            }
                        }
                    }
                    patches.push(Patch {
                        surface: Surface::Plane(Frame {
                            origin,
                            x: frame.u,
                            y: frame.v,
                            z: self.chart.n,
                        }),
                        same_sense: up,
                        loops: face,
                        origins,
                        support: Some(support),
                    });
                }
            }
        }
        Ok(patches)
    }

    /// Whether `side` of piece `index` lies in the region of a cap face
    /// whose loops are given: on its boundary, or strictly inside it.
    fn on_face(&self, face: &[Vec<Use>], index: usize, side: usize, up: bool) -> bool {
        let mut on_boundary = false;
        for &(key, forward) in face.iter().flatten() {
            if key_piece(key) == Some(index) {
                // Up caps keep their region on the left of each use, down
                // caps on the right.
                if (forward == up) == (side == 0) {
                    return true;
                }
                on_boundary = true;
            }
        }
        let point = self.arrangement.pieces[index].shape.point(0.5);
        !on_boundary
            && face
                .iter()
                .filter(|ring| arrangement::contains(self.ring_shapes(ring), point))
                .count()
                % 2
                == 1
    }

    fn ring_shapes<'a>(&'a self, ring: &'a [Use]) -> impl Iterator<Item = Shape> + 'a {
        ring.iter().filter_map(|&(key, _)| match key {
            EdgeKey::Across(p, _) => Some(self.arrangement.pieces[p].shape),
            EdgeKey::Up(..) => None,
        })
    }

    /// Walls on each support where the material differs across a piece,
    /// merged across adjacent pieces and slabs with the same outward side.
    fn walls(&self) -> Result<Vec<Patch>, GeometryError> {
        let pieces = &self.arrangement.pieces;
        // Cells (piece, slab) with the side of the material.
        let mut cells: Vec<(usize, usize, bool)> = vec![];
        let mut index = BTreeMap::new();
        for (p, piece) in pieces.iter().enumerate() {
            for slab in 0..self.slabs() {
                if let Some(left) = self.wall(piece, slab) {
                    index.insert((p, slab), cells.len());
                    cells.push((p, slab, left));
                }
            }
        }
        // Whether a cell's outward normal agrees with its support's.
        let sense = |&(p, _, left): &(usize, usize, bool)| -> bool {
            let piece = &pieces[p];
            let along = match (self.arrangement.supports[piece.support], piece.shape) {
                (Support::Line { direction, .. }, Shape::Line { a, b }) => {
                    let d = arrangement::sub(b, a);
                    d[0] * direction[0] + d[1] * direction[1] > 0.
                }
                (_, Shape::Arc { sweep, .. }) => sweep > 0.,
                _ => unreachable!("lines lie on line supports"),
            };
            along == left
        };
        let mut group: Vec<usize> = (0..cells.len()).collect();
        fn find(group: &mut [usize], mut i: usize) -> usize {
            while group[i] != i {
                group[i] = group[group[i]];
                i = group[i];
            }
            i
        }
        let join = |group: &mut Vec<usize>, a: usize, b: usize| {
            let (a, b) = (find(group, a), find(group, b));
            group[a] = b;
        };
        // Pieces of each support meeting at each node.
        let mut at_node: BTreeMap<(usize, usize), Vec<usize>> = BTreeMap::new();
        for (p, piece) in pieces.iter().enumerate() {
            for node in piece.nodes {
                at_node.entry((piece.support, node)).or_default().push(p);
            }
        }
        for (c, cell) in cells.iter().enumerate() {
            let (p, slab, _) = *cell;
            if let Some(&next) = index.get(&(p, slab + 1))
                && sense(&cells[next]) == sense(cell)
            {
                join(&mut group, c, next);
            }
            for node in pieces[p].nodes {
                for &q in &at_node[&(pieces[p].support, node)] {
                    if q != p
                        && let Some(&other) = index.get(&(q, slab))
                        && sense(&cells[other]) == sense(cell)
                    {
                        join(&mut group, c, other);
                    }
                }
            }
        }
        let mut faces: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for c in 0..cells.len() {
            let root = find(&mut group, c);
            faces.entry(root).or_default().push(c);
        }
        faces
            .into_values()
            .map(|members| {
                let first = cells[members[0]];
                let piece = &pieces[first.0];
                let same_sense = sense(&first);
                let mut uses = vec![];
                let mut origins = BTreeSet::new();
                for &c in &members {
                    let (p, slab, left) = cells[c];
                    let [a, b] = pieces[p].nodes;
                    // Counter-clockwise seen from outside the material.
                    let ring = [
                        (EdgeKey::Across(p, slab), true),
                        (EdgeKey::Up(b, slab), true),
                        (EdgeKey::Across(p, slab + 1), false),
                        (EdgeKey::Up(a, slab), false),
                    ];
                    if left {
                        uses.extend(ring);
                    } else {
                        uses.extend(ring.iter().rev().map(|&(k, f)| (k, !f)));
                    }
                    for &(segment, _) in &pieces[p].marks {
                        let [o, r, k] = self.owners[segment];
                        if self.active(o, slab) {
                            origins.extend(self.operands[o].sides[r][k].iter().cloned());
                        }
                    }
                }
                let surface = match self.arrangement.supports[piece.support] {
                    Support::Line { point, direction } => {
                        let x = self.chart.vector(direction);
                        Surface::Plane(Frame {
                            origin: self.chart.world(point, 0.),
                            x,
                            y: self.chart.n,
                            z: cross(x, self.chart.n),
                        })
                    }
                    Support::Circle { center, radius } => Surface::Cylinder {
                        frame: Frame {
                            origin: self.chart.world(center, 0.),
                            x: self.chart.frame.u,
                            y: self.chart.frame.v,
                            z: self.chart.n,
                        },
                        radius,
                    },
                };
                let support = match &surface {
                    Surface::Plane(frame) => Some(ReferenceFrame {
                        origin_mm: frame.origin,
                        u: if same_sense {
                            frame.x
                        } else {
                            frame.x.map(|x| -x)
                        },
                        v: frame.y,
                    }),
                    _ => None,
                };
                Ok(Patch {
                    surface,
                    same_sense,
                    loops: self.trace(&uses)?,
                    origins,
                    support,
                })
            })
            .collect()
    }

    fn start(&self, (key, forward): Use) -> (usize, usize) {
        let [s, e] = self.ends(key);
        if forward { s } else { e }
    }

    fn end(&self, (key, forward): Use) -> (usize, usize) {
        let [s, e] = self.ends(key);
        if forward { e } else { s }
    }

    /// The (node, level) vertices of an edge key in its own direction.
    fn ends(&self, key: EdgeKey) -> [(usize, usize); 2] {
        match key {
            EdgeKey::Across(p, level) => self.arrangement.pieces[p].nodes.map(|n| (n, level)),
            EdgeKey::Up(n, slab) => [(n, slab), (n, slab + 1)],
        }
    }

    /// Closed loops of directed uses; uses of one key in both directions
    /// cancel. A vertex with two ways on is a touching contact, refused.
    fn trace(&self, uses: &[Use]) -> Result<Vec<Vec<Use>>, GeometryError> {
        let mut count: BTreeMap<EdgeKey, i32> = BTreeMap::new();
        for &(key, forward) in uses {
            *count.entry(key).or_default() += if forward { 1 } else { -1 };
        }
        let mut leaving: BTreeMap<(usize, usize), Vec<Use>> = BTreeMap::new();
        for (&key, &n) in &count {
            match n {
                0 => {}
                1 | -1 => {
                    let used = (key, n == 1);
                    leaving.entry(self.start(used)).or_default().push(used);
                }
                _ => return Err(touching()),
            }
        }
        if leaving.values().any(|list| list.len() != 1) {
            return Err(touching());
        }
        let mut loops = vec![];
        while let Some((&first, _)) = leaving.iter().next() {
            let mut ring = vec![];
            let mut at = first;
            while let Some(list) = leaving.remove(&at) {
                let used = list[0];
                ring.push(used);
                at = self.end(used);
            }
            if at != first {
                return Err(invalid(
                    "A boundary loop of the profile Boolean did not close",
                ));
            }
            loops.push(ring);
        }
        Ok(loops)
    }

    /// Group cap loops into faces, outer loop first. For a cap facing up
    /// its outer loops run counter-clockwise in the chart; facing down they
    /// run clockwise.
    fn nest(&self, loops: Vec<Vec<Use>>, up: bool) -> Result<Vec<Vec<Vec<Use>>>, GeometryError> {
        let area = |ring: &Vec<Use>| -> f64 {
            ring.iter()
                .map(|&(key, forward)| match key {
                    EdgeKey::Across(p, _) => {
                        let shape = self.arrangement.pieces[p].shape;
                        let shape = if forward { shape } else { reversed(shape) };
                        shape.twice_area()
                    }
                    EdgeKey::Up(..) => 0.,
                })
                .sum::<f64>()
                * if up { 0.5 } else { -0.5 }
        };
        let (mut outers, holes): (Vec<_>, Vec<_>) = loops
            .into_iter()
            .map(|r| (area(&r), r))
            .partition(|(a, _)| *a > 0.);
        outers.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut faces: Vec<Vec<Vec<Use>>> = outers.iter().map(|(_, r)| vec![r.clone()]).collect();
        for (_, hole) in holes {
            let EdgeKey::Across(p, _) = hole[0].0 else {
                unreachable!("cap loops are horizontal")
            };
            let point = self.arrangement.pieces[p].shape.point(0.5);
            // The smallest outer loop around the hole holds it.
            let owner = outers
                .iter()
                .position(|(_, outer)| arrangement::contains(self.ring_shapes(outer), point))
                .ok_or_else(|| invalid("A cap hole of the profile Boolean has no outer loop"))?;
            faces[owner].push(hole);
        }
        Ok(faces)
    }

    /// The exact result volume: each slab's section area by its height.
    fn volume(&self) -> f64 {
        (0..self.slabs())
            .map(|k| self.section_area(k) * (self.levels[k + 1] - self.levels[k]))
            .sum()
    }

    fn solid(&self, feature: EntityId) -> Result<Option<Solid>, GeometryError> {
        let expected = self.volume();
        if (0..self.slabs()).all(|k| {
            self.arrangement
                .pieces
                .iter()
                .all(|p| self.wall(p, k).is_none())
        }) {
            return Ok(None);
        }
        let patches = self.patches()?;
        let (boundary, names) = self.assemble(patches, feature)?;
        let tolerance = self.policy.linear_tolerance_mm;
        let solid = brep::Solid::new(boundary, tolerance, None).map_err(|cause| match cause {
            brep::Error::Limit => error(
                OPERATION,
                GeometryErrorKind::Limit,
                "The profile Boolean result exceeds the retained topology limit",
            ),
            _ => refused(
                "The profile Boolean result is not one closed solid; it would split into \
                 several pieces or have faces that only touch",
            ),
        })?;
        let volume = solid.volume_mm3();
        if (volume - expected).abs() > solid.volume_error_mm3().max(1e-9 * expected.abs()) {
            return Err(invalid(format!(
                "The profile Boolean boundary integrates to {volume} mm3, not the exact {expected} mm3"
            )));
        }
        Solid::from_authored_brep(solid, feature, names, self.policy).map(Some)
    }

    /// Exact edges, vertices and faces, with collinear and co-circular edge
    /// runs through two-edge vertices merged into single edges.
    fn assemble(
        &self,
        patches: Vec<Patch>,
        feature: EntityId,
    ) -> Result<(Boundary, AuthoredNaming), GeometryError> {
        let faces: Vec<&Vec<Vec<Use>>> = patches.iter().map(|patch| &patch.loops).collect();
        let mut keys: BTreeMap<EdgeKey, usize> = BTreeMap::new();
        let mut edges: Vec<Run> = vec![];
        for ring in faces.iter().copied().flatten() {
            for &(key, _) in ring {
                keys.entry(key).or_insert_with(|| {
                    edges.push(self.run(key));
                    edges.len() - 1
                });
            }
        }
        let mut loops: Vec<Vec<Vec<(usize, bool)>>> = faces
            .iter()
            .map(|face| {
                face.iter()
                    .map(|ring| {
                        ring.iter()
                            .map(|&(key, forward)| {
                                let edge = keys[&key];
                                (edge, forward == edges[edge].along)
                            })
                            .collect()
                    })
                    .collect()
            })
            .collect();
        merge_runs(&mut edges, &mut loops);
        // Exact vertices and edges.
        let mut vertex_of: BTreeMap<(usize, usize), usize> = BTreeMap::new();
        let mut vertices = vec![];
        let mut ids: BTreeMap<usize, usize> = BTreeMap::new();
        let mut brep_edges = vec![];
        let tolerance = self.policy.linear_tolerance_mm;
        for ring in loops.iter().flatten() {
            for &(e, _) in ring {
                if ids.contains_key(&e) {
                    continue;
                }
                let run = &edges[e];
                let ends = run.ends.map(|v| {
                    *vertex_of.entry(v).or_insert_with(|| {
                        vertices.push(
                            self.chart
                                .world(self.arrangement.nodes[v.0], self.levels[v.1]),
                        );
                        vertices.len() - 1
                    })
                });
                let points = ends.map(|v| vertices[v]);
                let edge = BrepEdge::new(self.curve(run.curve), ends, points, 16. * tolerance)
                    .map_err(|_| invalid("A profile Boolean edge does not lie on its curve"))?;
                ids.insert(e, brep_edges.len());
                brep_edges.push(edge);
            }
        }
        let mut brep_faces = vec![];
        let mut names = vec![];
        for (index, (patch, face)) in patches.into_iter().zip(loops).enumerate() {
            let mut origins: Vec<Origin> = patch.origins.into_iter().collect();
            if origins.is_empty() {
                origins.push(Origin {
                    source: feature,
                    role: format!("profile-stack/face/{index}"),
                });
            }
            brep_faces.push(BrepFace {
                surface: patch.surface,
                same_sense: patch.same_sense,
                loops: face
                    .into_iter()
                    .map(|ring| Loop {
                        uses: ring.into_iter().map(|(e, f)| (ids[&e], f)).collect(),
                        vertex: None,
                    })
                    .collect(),
            });
            names.push(FaceName {
                origins,
                support: patch.support,
                anchor: patch.support,
            });
        }
        Ok((
            Boundary {
                vertices,
                edges: brep_edges,
                shells: vec![Shell {
                    faces: (0..brep_faces.len()).collect(),
                    reversed: false,
                }],
                faces: brep_faces,
            },
            AuthoredNaming { faces: names },
        ))
    }

    /// The edge of a key along its curve's own direction.
    fn run(&self, key: EdgeKey) -> Run {
        match key {
            EdgeKey::Up(n, slab) => Run {
                curve: CurveKey::Axis(n),
                ends: [(n, slab), (n, slab + 1)],
                along: true,
            },
            EdgeKey::Across(p, level) => {
                let piece = &self.arrangement.pieces[p];
                let along = match (self.arrangement.supports[piece.support], piece.shape) {
                    (Support::Line { direction, .. }, Shape::Line { a, b }) => {
                        let d = arrangement::sub(b, a);
                        d[0] * direction[0] + d[1] * direction[1] > 0.
                    }
                    (_, Shape::Arc { sweep, .. }) => sweep > 0.,
                    _ => unreachable!("lines lie on line supports"),
                };
                let [a, b] = piece.nodes.map(|n| (n, level));
                Run {
                    curve: CurveKey::Level(piece.support, level),
                    ends: if along { [a, b] } else { [b, a] },
                    along,
                }
            }
        }
    }

    fn curve(&self, key: CurveKey) -> Curve {
        match key {
            CurveKey::Axis(n) => Curve::Line {
                origin: self.chart.world(self.arrangement.nodes[n], 0.),
                direction: self.chart.n,
            },
            CurveKey::Level(support, level) => {
                let z = self.levels[level];
                match self.arrangement.supports[support] {
                    Support::Line { point, direction } => Curve::Line {
                        origin: self.chart.world(point, z),
                        direction: self.chart.vector(direction),
                    },
                    Support::Circle { center, radius } => Curve::Circle {
                        frame: Frame {
                            origin: self.chart.world(center, z),
                            x: self.chart.frame.u,
                            y: self.chart.frame.v,
                            z: self.chart.n,
                        },
                        radius,
                    },
                }
            }
        }
    }
}

fn key_piece(key: EdgeKey) -> Option<usize> {
    match key {
        EdgeKey::Across(piece, _) => Some(piece),
        EdgeKey::Up(..) => None,
    }
}

fn touching() -> GeometryError {
    refused(
        "The profile Boolean result would have faces or edges that only touch along a line or \
         at a point; this exact Boolean does not support that contact",
    )
}

fn invalid(message: impl Into<String>) -> GeometryError {
    error(OPERATION, GeometryErrorKind::InvalidResult, message)
}

/// Which exact curve an edge runs along.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum CurveKey {
    /// A vertical line through a node.
    Axis(usize),
    /// A support's line or circle at a level.
    Level(usize, usize),
}

/// An edge between (node, level) vertices, in its curve's direction.
/// `along` records whether the key it came from runs the same way.
#[derive(Clone, Debug)]
struct Run {
    curve: CurveKey,
    ends: [(usize, usize); 2],
    along: bool,
}

/// Merge consecutive edges on one curve through vertices that bound
/// exactly those two edges, rewriting every loop that uses them.
fn merge_runs(edges: &mut Vec<Run>, loops: &mut [Vec<Vec<(usize, bool)>>]) {
    let mut alive = vec![true; edges.len()];
    loop {
        let mut incident: BTreeMap<(usize, usize), BTreeSet<usize>> = BTreeMap::new();
        for (e, run) in edges.iter().enumerate() {
            if alive[e] {
                for v in run.ends {
                    incident.entry(v).or_default().insert(e);
                }
            }
        }
        let candidate = incident.iter().find_map(|(&v, list)| {
            let mut list = list.iter().copied();
            let (first, second) = (list.next()?, list.next()?);
            if list.next().is_some() || first == second {
                return None;
            }
            let (a, b) = (&edges[first], &edges[second]);
            if a.curve != b.curve || a.ends[0] == a.ends[1] || b.ends[0] == b.ends[1] {
                return None;
            }
            // `before` ends at v, `after` starts there.
            if a.ends[1] == v && b.ends[0] == v {
                Some((first, second))
            } else if b.ends[1] == v && a.ends[0] == v {
                Some((second, first))
            } else {
                None
            }
        });
        let Some((before, after)) = candidate else {
            break;
        };
        let merged = edges.len();
        edges.push(Run {
            curve: edges[before].curve,
            ends: [edges[before].ends[0], edges[after].ends[1]],
            along: true,
        });
        alive.push(true);
        alive[before] = false;
        alive[after] = false;
        for ring in loops.iter_mut().flatten() {
            let n = ring.len();
            let Some(i) = (0..n).find(|&i| {
                let (x, y) = (ring[i], ring[(i + 1) % n]);
                (x == (before, true) && y == (after, true))
                    || (x == (after, false) && y == (before, false))
            }) else {
                continue;
            };
            let forward = ring[i].1;
            let j = (i + 1) % n;
            ring[i] = (merged, forward);
            ring.remove(j);
        }
    }
}
