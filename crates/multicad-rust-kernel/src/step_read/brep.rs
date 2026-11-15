//! Part 21 reader for the retained boundary of general solids.
//!
//! The reader keeps what the file states: support surfaces, edge curves,
//! oriented loops, shells and voids. It reconstructs nothing and heals
//! nothing; validation and integration belong to [`crate::brep::Solid`].
use super::part21::{self, Entity, Value};
use super::{ReadError, Result, numbers, reference, refs, sense};
use crate::brep::{
    Boundary, Curve, Edge, Face, Frame, Loop, NurbsCurve, NurbsSurface, Point, Shell, Solid,
    Surface,
};
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;

/// The largest number of solids read from one file.
pub const MAX_SOLIDS: usize = 1024;

/// Product structure that places or repeats shapes. It belongs to the
/// assembly layer, which calls [`Source::boundary`] per solid.
const ASSEMBLY: [&str; 5] = [
    "NEXT_ASSEMBLY_USAGE_OCCURRENCE",
    "CONTEXT_DEPENDENT_SHAPE_REPRESENTATION",
    "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION",
    "ITEM_DEFINED_TRANSFORMATION",
    "MAPPED_ITEM",
];

/// Geometry that is not a solid. A file holding any keeps its meaning and
/// is refused rather than imported without it.
const NOT_SOLID: [&str; 8] = [
    "SHELL_BASED_SURFACE_MODEL",
    "FACE_BASED_SURFACE_MODEL",
    "OPEN_SHELL",
    "GEOMETRIC_SET",
    "GEOMETRIC_CURVE_SET",
    "FACETED_BREP",
    "SHELL_BASED_WIREFRAME_MODEL",
    "EDGE_BASED_WIREFRAME_MODEL",
];

fn check(ok: bool, error: ReadError) -> Result<()> {
    if ok { Ok(()) } else { Err(error) }
}

fn number(value: &Value<'_>) -> Result<f64> {
    match value {
        Value::Number(x) if x.is_finite() => Ok(*x),
        _ => Err(ReadError::Invalid),
    }
}

fn count(value: &Value<'_>) -> Result<usize> {
    let x = number(value)?;
    check(x >= 0. && x.fract() == 0. && x <= 1e6, ReadError::Invalid)?;
    Ok(x as usize)
}

fn list<'a, 'i>(value: &'a Value<'i>) -> Result<&'a [Value<'i>]> {
    match value {
        Value::List(items) => Ok(items),
        _ => Err(ReadError::Invalid),
    }
}

fn reals(value: &Value<'_>) -> Result<Vec<f64>> {
    list(value)?.iter().map(number).collect()
}

fn counts(value: &Value<'_>) -> Result<Vec<usize>> {
    list(value)?.iter().map(count).collect()
}

fn length(value: &Value<'_>) -> Result<f64> {
    let x = number(value)?;
    check(x.abs() <= 1e12, ReadError::Limit)?;
    Ok(x)
}

fn radius(value: &Value<'_>, zero: bool) -> Result<f64> {
    let radius = length(value)?;
    check(radius > 0. || (zero && radius == 0.), ReadError::Invalid)?;
    Ok(radius)
}

fn normalized(v: Point) -> Result<Point> {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    check(n.is_finite() && n > 0., ReadError::Invalid)?;
    Ok(v.map(|x| x / n))
}

fn cross(a: Point, b: Point) -> Point {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// A parsed exchange file whose solids can be read one at a time.
pub struct Source<'i> {
    entities: BTreeMap<u32, Vec<Entity<'i>>>,
    tolerance: f64,
    angle: f64,
    solids: Vec<u32>,
}

impl<'i> Source<'i> {
    /// Parse the file and check its units. Lengths must be millimetres;
    /// plane angles may be radians or a stated multiple of them.
    pub fn parse(input: &'i str, tolerance: f64) -> Result<Self> {
        let entities = part21::parse(input)?;
        let mut solids = Vec::new();
        let mut angle = None;
        for (&id, components) in &entities {
            for entity in components {
                check(!NOT_SOLID.contains(&entity.name), ReadError::Unsupported)?;
                if matches!(entity.name, "MANIFOLD_SOLID_BREP" | "BREP_WITH_VOIDS") {
                    solids.push(id);
                }
                if entity.name == "GLOBAL_UNIT_ASSIGNED_CONTEXT" {
                    let [units] = entity.args.as_slice() else {
                        return Err(ReadError::Invalid);
                    };
                    let factor = units_of(&entities, &refs(units)?)?;
                    check(
                        angle.is_none_or(|known| known == factor),
                        ReadError::Unsupported,
                    )?;
                    angle = Some(factor);
                }
            }
        }
        check(!solids.is_empty(), ReadError::Unsupported)?;
        check(solids.len() <= MAX_SOLIDS, ReadError::Limit)?;
        let angle = angle.ok_or(ReadError::Unsupported)?;
        Ok(Self {
            entities,
            tolerance,
            angle,
            solids,
        })
    }

    /// True when the file carries product structure that places shapes.
    pub fn is_assembly(&self) -> bool {
        self.entities
            .values()
            .flatten()
            .any(|entity| ASSEMBLY.contains(&entity.name))
    }

    /// Entity ids of the solids, in file order.
    pub fn solids(&self) -> &[u32] {
        &self.solids
    }

    /// Every face of the file belongs to one of its solids.
    pub fn owns_every_face(&self, boundaries: &[Boundary]) -> bool {
        let faces = self
            .entities
            .values()
            .flatten()
            .filter(|e| e.name == "ADVANCED_FACE")
            .count();
        faces
            == boundaries
                .iter()
                .map(|boundary| boundary.faces.len())
                .sum::<usize>()
    }

    /// The retained boundary of one solid.
    pub fn boundary(&self, solid: u32) -> Result<Boundary> {
        let mut reader = Reader {
            source: self,
            boundary: Boundary::default(),
            vertices: BTreeMap::new(),
            edges: BTreeMap::new(),
            faces: BTreeSet::new(),
        };
        let entity = self.simple(solid)?;
        let shells = match (entity.name, entity.args.as_slice()) {
            ("MANIFOLD_SOLID_BREP", [_, outer]) => vec![(reference(outer)?, false)],
            ("BREP_WITH_VOIDS", [_, outer, voids]) => {
                let mut shells = vec![(reference(outer)?, false)];
                for void in refs(voids)? {
                    let entity = self.simple(void)?;
                    shells.push(match (entity.name, entity.args.as_slice()) {
                        ("ORIENTED_CLOSED_SHELL", [_, _, shell, orientation]) => {
                            (reference(shell)?, !sense(orientation)?)
                        }
                        _ => (void, false),
                    });
                }
                shells
            }
            _ => return Err(ReadError::Invalid),
        };
        for (shell, reversed) in shells {
            let [_, faces] = self.args(shell, "CLOSED_SHELL")? else {
                return Err(ReadError::Invalid);
            };
            let mut indices = Vec::new();
            for face in refs(faces)? {
                indices.push(reader.face(face)?);
            }
            reader.boundary.shells.push(Shell {
                faces: indices,
                reversed,
            });
        }
        Ok(reader.boundary)
    }

    fn components(&self, id: u32) -> Result<&[Entity<'i>]> {
        self.entities
            .get(&id)
            .map(Vec::as_slice)
            .ok_or(ReadError::Invalid)
    }

    fn simple(&self, id: u32) -> Result<&Entity<'i>> {
        match self.components(id)? {
            [entity] => Ok(entity),
            _ => Err(ReadError::Unsupported),
        }
    }

    fn args(&self, id: u32, name: &str) -> Result<&[Value<'i>]> {
        let entity = self.simple(id)?;
        check(entity.name == name, ReadError::Unsupported)?;
        Ok(&entity.args)
    }

    fn point(&self, id: u32) -> Result<Point> {
        let [_, values] = self.args(id, "CARTESIAN_POINT")? else {
            return Err(ReadError::Invalid);
        };
        numbers(values)
    }

    fn direction(&self, id: u32) -> Result<Point> {
        let [_, values] = self.args(id, "DIRECTION")? else {
            return Err(ReadError::Invalid);
        };
        normalized(numbers(values)?)
    }

    /// `VECTOR`: a direction scaled by its magnitude.
    fn vector(&self, id: u32) -> Result<Point> {
        let [_, direction, magnitude] = self.args(id, "VECTOR")? else {
            return Err(ReadError::Invalid);
        };
        let magnitude = length(magnitude)?;
        check(magnitude > 0., ReadError::Invalid)?;
        Ok(self
            .direction(reference(direction)?)?
            .map(|x| x * magnitude))
    }

    /// `AXIS2_PLACEMENT_3D` with the standard defaults for omitted axes.
    fn frame(&self, id: u32) -> Result<Frame> {
        let [_, origin, axis, reference_direction] = self.args(id, "AXIS2_PLACEMENT_3D")? else {
            return Err(ReadError::Invalid);
        };
        let optional = |value: &Value<'i>, default: Point| match value {
            Value::Omitted => Ok(default),
            value => self.direction(reference(value)?),
        };
        let z = optional(axis, [0., 0., 1.])?;
        let default = if z[0].abs() < 0.9 {
            [1., 0., 0.]
        } else {
            [0., 1., 0.]
        };
        let hint = optional(reference_direction, default)?;
        let y = normalized(cross(z, hint))?;
        Ok(Frame {
            origin: self.point(reference(origin)?)?,
            x: cross(y, z),
            y,
            z,
        })
    }

    fn points(&self, value: &Value<'i>) -> Result<Vec<Point>> {
        refs(value)?.into_iter().map(|id| self.point(id)).collect()
    }

    fn component(&self, id: u32, name: &str) -> Result<Option<&[Value<'i>]>> {
        Ok(self
            .components(id)?
            .iter()
            .find(|e| e.name == name)
            .map(|e| e.args.as_slice()))
    }

    fn spline_curve(&self, id: u32) -> Result<Curve> {
        let parts = self.components(id)?;
        let (base, knots, weights) = match parts {
            [entity] => match (entity.name, entity.args.as_slice()) {
                ("B_SPLINE_CURVE_WITH_KNOTS", [_, base @ .., multiplicities, knots, _])
                    if base.len() == 5 =>
                {
                    (base, Some((multiplicities, knots)), None)
                }
                ("BEZIER_CURVE", [_, base @ ..]) if base.len() == 5 => (base, None, None),
                _ => return Err(ReadError::Unsupported),
            },
            _ => {
                let base = self
                    .component(id, "B_SPLINE_CURVE")?
                    .ok_or(ReadError::Unsupported)?;
                let knots = match self.component(id, "B_SPLINE_CURVE_WITH_KNOTS")? {
                    Some([multiplicities, knots, _]) => Some((multiplicities, knots)),
                    Some(_) => return Err(ReadError::Invalid),
                    None => {
                        check(
                            self.component(id, "BEZIER_CURVE")?.is_some(),
                            ReadError::Unsupported,
                        )?;
                        None
                    }
                };
                let weights = match self.component(id, "RATIONAL_B_SPLINE_CURVE")? {
                    Some([weights]) => Some(reals(weights)?),
                    Some(_) => return Err(ReadError::Invalid),
                    None => None,
                };
                (base, knots, weights)
            }
        };
        let [degree, points, ..] = base else {
            return Err(ReadError::Invalid);
        };
        let degree = count(degree)?;
        let (knots, multiplicities) = match knots {
            Some((multiplicities, knots)) => (reals(knots)?, counts(multiplicities)?),
            None => (vec![0., 1.], vec![degree + 1; 2]),
        };
        let curve = NurbsCurve::new(degree, self.points(points)?, weights, knots, multiplicities)?;
        Ok(Curve::nurbs(curve, self.tolerance))
    }

    fn curve(&self, id: u32, depth: usize) -> Result<Curve> {
        let parts = self.components(id)?;
        let [entity] = parts else {
            return self.spline_curve(id);
        };
        match (entity.name, entity.args.as_slice()) {
            ("LINE", [_, origin, direction]) => Ok(Curve::Line {
                origin: self.point(reference(origin)?)?,
                direction: self.vector(reference(direction)?)?,
            }),
            ("CIRCLE", [_, frame, r]) => Ok(Curve::Circle {
                frame: self.frame(reference(frame)?)?,
                radius: radius(r, false)?,
            }),
            ("ELLIPSE", [_, frame, a, b]) => {
                let semiaxes = [length(a)?, length(b)?];
                check(semiaxes[0] > 0. && semiaxes[1] > 0., ReadError::Invalid)?;
                Ok(Curve::Ellipse {
                    frame: self.frame(reference(frame)?)?,
                    semiaxes,
                })
            }
            ("SURFACE_CURVE" | "SEAM_CURVE" | "INTERSECTION_CURVE", [_, curve, ..]) => {
                check(depth < 8, ReadError::Invalid)?;
                self.curve(reference(curve)?, depth + 1)
            }
            ("B_SPLINE_CURVE_WITH_KNOTS" | "BEZIER_CURVE", _) => self.spline_curve(id),
            _ => Err(ReadError::Unsupported),
        }
    }

    fn spline_surface(&self, id: u32) -> Result<Surface> {
        let parts = self.components(id)?;
        let (base, knots, weights) = match parts {
            [entity] => match (entity.name, entity.args.as_slice()) {
                ("B_SPLINE_SURFACE_WITH_KNOTS", [_, base @ .., mu, mv, ku, kv, _])
                    if base.len() == 7 =>
                {
                    (base, Some([mu, mv, ku, kv]), None)
                }
                ("BEZIER_SURFACE", [_, base @ ..]) if base.len() == 7 => (base, None, None),
                _ => return Err(ReadError::Unsupported),
            },
            _ => {
                let base = self
                    .component(id, "B_SPLINE_SURFACE")?
                    .ok_or(ReadError::Unsupported)?;
                let knots = match self.component(id, "B_SPLINE_SURFACE_WITH_KNOTS")? {
                    Some([mu, mv, ku, kv, _]) => Some([mu, mv, ku, kv]),
                    Some(_) => return Err(ReadError::Invalid),
                    None => {
                        check(
                            self.component(id, "BEZIER_SURFACE")?.is_some(),
                            ReadError::Unsupported,
                        )?;
                        None
                    }
                };
                let weights = match self.component(id, "RATIONAL_B_SPLINE_SURFACE")? {
                    Some([weights]) => {
                        let mut flat = Vec::new();
                        for row in list(weights)? {
                            flat.extend(reals(row)?);
                        }
                        Some(flat)
                    }
                    Some(_) => return Err(ReadError::Invalid),
                    None => None,
                };
                (base, knots, weights)
            }
        };
        let [du, dv, rows, ..] = base else {
            return Err(ReadError::Invalid);
        };
        let degrees = [count(du)?, count(dv)?];
        let rows = list(rows)?;
        let mut points = Vec::new();
        for row in rows {
            let row = self.points(row)?;
            check(points.len() % row.len().max(1) == 0, ReadError::Invalid)?;
            points.extend(row);
        }
        let counts_uv = [rows.len(), points.len() / rows.len().max(1)];
        let (knots, multiplicities) = match knots {
            Some([mu, mv, ku, kv]) => ([reals(ku)?, reals(kv)?], [counts(mu)?, counts(mv)?]),
            None => (
                [vec![0., 1.], vec![0., 1.]],
                degrees.map(|d| vec![d + 1; 2]),
            ),
        };
        let surface =
            NurbsSurface::new(degrees, counts_uv, points, weights, knots, multiplicities)?;
        Ok(Surface::nurbs(surface, self.tolerance))
    }

    fn surface(&self, id: u32) -> Result<Surface> {
        let parts = self.components(id)?;
        let [entity] = parts else {
            return self.spline_surface(id);
        };
        match (entity.name, entity.args.as_slice()) {
            ("PLANE", [_, frame]) => Ok(Surface::Plane(self.frame(reference(frame)?)?)),
            ("CYLINDRICAL_SURFACE", [_, frame, r]) => Ok(Surface::Cylinder {
                frame: self.frame(reference(frame)?)?,
                radius: radius(r, false)?,
            }),
            ("CONICAL_SURFACE", [_, frame, r, angle]) => {
                let half_angle = number(angle)? * self.angle;
                check(half_angle > 0. && half_angle < PI / 2., ReadError::Invalid)?;
                Ok(Surface::Cone {
                    frame: self.frame(reference(frame)?)?,
                    radius: radius(r, true)?,
                    half_angle,
                })
            }
            ("SPHERICAL_SURFACE", [_, frame, r]) => Ok(Surface::Sphere {
                frame: self.frame(reference(frame)?)?,
                radius: radius(r, false)?,
            }),
            ("TOROIDAL_SURFACE", [_, frame, major, minor]) => {
                let (major, minor) = (radius(major, false)?, radius(minor, false)?);
                // A spindle torus intersects itself along its axis.
                check(major > minor, ReadError::Unsupported)?;
                Ok(Surface::Torus {
                    frame: self.frame(reference(frame)?)?,
                    major,
                    minor,
                })
            }
            ("SURFACE_OF_LINEAR_EXTRUSION", [_, curve, direction]) => Ok(Surface::Extrusion {
                curve: self.curve(reference(curve)?, 0)?,
                direction: self.vector(reference(direction)?)?,
            }),
            ("SURFACE_OF_REVOLUTION", [_, curve, axis]) => {
                let [_, origin, direction] = self.args(reference(axis)?, "AXIS1_PLACEMENT")? else {
                    return Err(ReadError::Invalid);
                };
                Ok(Surface::Revolution {
                    curve: self.curve(reference(curve)?, 0)?,
                    origin: self.point(reference(origin)?)?,
                    axis: self.direction(reference(direction)?)?,
                })
            }
            ("B_SPLINE_SURFACE_WITH_KNOTS" | "BEZIER_SURFACE", _) => self.spline_surface(id),
            _ => Err(ReadError::Unsupported),
        }
    }
}

/// Check the length unit and return the plane-angle unit in radians.
fn units_of(entities: &BTreeMap<u32, Vec<Entity<'_>>>, units: &[u32]) -> Result<f64> {
    let mut millimetres = 0;
    let mut angle = None;
    for id in units {
        let unit = entities.get(id).ok_or(ReadError::Invalid)?;
        let has = |name: &str| unit.iter().find(|e| e.name == name);
        let si = |prefix: Value<'_>, name: &str| {
            has("SI_UNIT").is_some_and(|e| e.args.as_slice() == [prefix, Value::Enum(name)])
        };
        if has("LENGTH_UNIT").is_some() {
            check(si(Value::Enum("MILLI"), "METRE"), ReadError::Unsupported)?;
            millimetres += 1;
        }
        if has("PLANE_ANGLE_UNIT").is_some() {
            let factor = if si(Value::Omitted, "RADIAN") {
                1.
            } else {
                let conversion = has("CONVERSION_BASED_UNIT").ok_or(ReadError::Unsupported)?;
                let [_, measure] = conversion.args.as_slice() else {
                    return Err(ReadError::Invalid);
                };
                let measure = entities
                    .get(&reference(measure)?)
                    .ok_or(ReadError::Invalid)?;
                let [Entity { args, .. }] = measure.as_slice() else {
                    return Err(ReadError::Unsupported);
                };
                let [Value::Typed(Entity { args: value, .. }), _] = args.as_slice() else {
                    return Err(ReadError::Invalid);
                };
                let [value] = value.as_slice() else {
                    return Err(ReadError::Invalid);
                };
                let factor = number(value)?;
                check(factor > 0., ReadError::Invalid)?;
                if (factor - PI / 180.).abs() <= 1e-9 {
                    PI / 180.
                } else {
                    factor
                }
            };
            check(angle.is_none(), ReadError::Unsupported)?;
            angle = Some(factor);
        }
    }
    check(millimetres == 1, ReadError::Unsupported)?;
    angle.ok_or(ReadError::Unsupported)
}

struct Reader<'s, 'i> {
    source: &'s Source<'i>,
    boundary: Boundary,
    vertices: BTreeMap<u32, usize>,
    /// Edge index and whether the file's edge runs against its curve.
    edges: BTreeMap<u32, (usize, bool)>,
    faces: BTreeSet<u32>,
}

impl Reader<'_, '_> {
    fn vertex(&mut self, id: u32) -> Result<usize> {
        if let Some(&index) = self.vertices.get(&id) {
            return Ok(index);
        }
        let [_, point] = self.source.args(id, "VERTEX_POINT")? else {
            return Err(ReadError::Invalid);
        };
        self.boundary
            .vertices
            .push(self.source.point(reference(point)?)?);
        self.vertices.insert(id, self.boundary.vertices.len() - 1);
        Ok(self.boundary.vertices.len() - 1)
    }

    fn edge(&mut self, id: u32) -> Result<(usize, bool)> {
        if let Some(&known) = self.edges.get(&id) {
            return Ok(known);
        }
        let [_, start, end, curve, same_sense] = self.source.args(id, "EDGE_CURVE")? else {
            return Err(ReadError::Invalid);
        };
        let against = !sense(same_sense)?;
        let mut vertices = [
            self.vertex(reference(start)?)?,
            self.vertex(reference(end)?)?,
        ];
        if against {
            vertices.swap(0, 1);
        }
        let points = vertices.map(|vertex| self.boundary.vertices[vertex]);
        let curve = self.source.curve(reference(curve)?, 0)?;
        let edge = Edge::new(curve, vertices, points, self.source.tolerance)?;
        self.boundary.edges.push(edge);
        self.edges
            .insert(id, (self.boundary.edges.len() - 1, against));
        Ok((self.boundary.edges.len() - 1, against))
    }

    fn trim(&mut self, id: u32) -> Result<Loop> {
        let entity = self.source.simple(id)?;
        let (bound, orientation) = match (entity.name, entity.args.as_slice()) {
            ("FACE_BOUND" | "FACE_OUTER_BOUND", [_, bound, orientation]) => {
                (reference(bound)?, sense(orientation)?)
            }
            _ => return Err(ReadError::Invalid),
        };
        let entity = self.source.simple(bound)?;
        match (entity.name, entity.args.as_slice()) {
            ("VERTEX_LOOP", [_, vertex]) => Ok(Loop {
                uses: vec![],
                vertex: Some(self.vertex(reference(vertex)?)?),
            }),
            ("EDGE_LOOP", [_, edges]) => {
                let mut uses = Vec::new();
                for oriented in refs(edges)? {
                    let [_, _, _, edge, direction] = self.source.args(oriented, "ORIENTED_EDGE")?
                    else {
                        return Err(ReadError::Invalid);
                    };
                    let (index, against) = self.edge(reference(edge)?)?;
                    uses.push((index, sense(direction)? != against));
                }
                check(!uses.is_empty(), ReadError::Invalid)?;
                if !orientation {
                    uses.reverse();
                    for used in &mut uses {
                        used.1 = !used.1;
                    }
                }
                Ok(Loop { uses, vertex: None })
            }
            _ => Err(ReadError::Unsupported),
        }
    }

    fn face(&mut self, id: u32) -> Result<usize> {
        check(self.faces.insert(id), ReadError::Invalid)?;
        let [_, bounds, surface, same_sense] = self.source.args(id, "ADVANCED_FACE")? else {
            return Err(ReadError::Invalid);
        };
        let mut loops = Vec::new();
        for bound in refs(bounds)? {
            loops.push(self.trim(bound)?);
        }
        self.boundary.faces.push(Face {
            surface: self.source.surface(reference(surface)?)?,
            same_sense: sense(same_sense)?,
            loops,
        });
        Ok(self.boundary.faces.len() - 1)
    }
}

/// Read every solid of a flat file. Assemblies and files holding anything
/// but solids are refused as unsupported.
pub fn read_boundaries(input: &str, tolerance: f64) -> Result<Vec<Boundary>> {
    let source = Source::parse(input, tolerance)?;
    check(!source.is_assembly(), ReadError::Unsupported)?;
    let boundaries = source
        .solids()
        .iter()
        .map(|&solid| source.boundary(solid))
        .collect::<Result<Vec<_>>>()?;
    check(source.owns_every_face(&boundaries), ReadError::Unsupported)?;
    Ok(boundaries)
}

/// Read, validate and integrate every solid of a flat file.
pub fn read_solids(
    input: &str,
    tolerance: f64,
    interrupt: crate::brep::Interrupt<'_>,
) -> std::result::Result<Vec<Solid>, crate::brep::Error> {
    let boundaries = read_boundaries(input, tolerance).map_err(|error| match error {
        ReadError::Invalid => crate::brep::Error::Invalid,
        ReadError::Unsupported => crate::brep::Error::Unsupported,
        ReadError::Limit => crate::brep::Error::Limit,
    })?;
    boundaries
        .into_iter()
        .map(|boundary| Solid::new(boundary, tolerance, interrupt))
        .collect()
}

/// Volume, bounds and face count of each solid from its retained boundary.
pub(super) fn readbacks(input: &str, tolerance: f64) -> Result<Vec<super::Readback>> {
    let solids = read_solids(input, tolerance, None)?;
    Ok(solids
        .iter()
        .map(|solid| super::Readback {
            bounds_mm: solid.bounds_mm(),
            volume_mm3: solid.volume_mm3(),
            face_count: solid.faces().len(),
            edges_mm: vec![],
            circles: vec![],
            semicircles: vec![],
            arcs: vec![],
            ellipses: vec![],
            splines: vec![],
            solid_count: 1,
        })
        .collect())
}

impl From<crate::brep::Error> for ReadError {
    fn from(error: crate::brep::Error) -> Self {
        match error {
            crate::brep::Error::Unsupported => Self::Unsupported,
            crate::brep::Error::Limit => Self::Limit,
            crate::brep::Error::Invalid | crate::brep::Error::Cancelled => Self::Invalid,
        }
    }
}
