use super::*;
use crate::{Rational, step::write_box};

const TOL: f64 = 1e-7;

#[derive(Clone)]
enum Ring {
    Polygon(Vec<usize>),
    /// A full circle about +z through one vertex.
    Circle {
        vertex: usize,
        center: Point,
        forward: bool,
    },
}

#[derive(Clone)]
struct Patch {
    normal: Point,
    rings: Vec<Ring>,
    /// Transmit the opposite support normal with reversed face and bound flags.
    reversed: bool,
    vertex_loop: bool,
}

#[derive(Clone, Copy, Default)]
struct Tweak {
    /// Move the line of this edge off its vertices.
    displaced_line: Option<[usize; 2]>,
    /// Use this edge forward in every loop.
    one_way: Option<[usize; 2]>,
    /// Transmit this edge as a longer straight spline.
    spline: Option<[usize; 2]>,
}

struct Text(Vec<String>);

impl Text {
    fn add(&mut self, record: String) -> u32 {
        self.0.push(record);
        self.0.len() as u32
    }
    fn triple(&mut self, kind: &str, [x, y, z]: Point) -> u32 {
        self.add(format!("{kind}('',({x:.12},{y:.12},{z:.12}))"))
    }
    fn axis(&mut self, origin: Point, normal: Point, reference: Point) -> u32 {
        let origin = self.triple("CARTESIAN_POINT", origin);
        let normal = self.triple("DIRECTION", normal);
        let reference = self.triple("DIRECTION", reference);
        self.add(format!(
            "AXIS2_PLACEMENT_3D('',#{origin},#{normal},#{reference})"
        ))
    }
}

fn flag(value: bool) -> &'static str {
    if value { ".T." } else { ".F." }
}

fn list(ids: &[u32]) -> String {
    let ids = ids.iter().map(|id| format!("#{id}")).collect::<Vec<_>>();
    ids.join(",")
}

fn step(points: &[Point], patches: &[Patch], tweak: Tweak) -> String {
    let mut text = Text(Vec::new());
    // Only the vertices of the listed rings are transmitted.
    let mut vertices = vec![0; points.len()];
    for ring in patches.iter().flat_map(|patch| &patch.rings) {
        let used = match ring {
            Ring::Polygon(ring) => ring.as_slice(),
            Ring::Circle { vertex, .. } => std::slice::from_ref(vertex),
        };
        for &vertex in used {
            if vertices[vertex] == 0 {
                let point = text.triple("CARTESIAN_POINT", points[vertex]);
                vertices[vertex] = text.add(format!("VERTEX_POINT('',#{point})"));
            }
        }
    }
    let mut edges = BTreeMap::<[usize; 2], u32>::new();
    let mut faces = Vec::new();
    for patch in patches {
        let mut bounds = Vec::new();
        for (index, ring) in patch.rings.iter().enumerate() {
            let mut uses = Vec::new();
            match ring {
                Ring::Polygon(ring) => {
                    for (&a, &b) in ring.iter().zip(ring.iter().cycle().skip(1)) {
                        let key = [a.min(b), a.max(b)];
                        let edge = *edges.entry(key).or_insert_with(|| {
                            let [low, high] = key.map(|vertex| points[vertex]);
                            let along = unit(sub(high, low)).unwrap();
                            let curve = if tweak.spline == Some(key) {
                                let ends = [-1., 1.].map(|side: f64| {
                                    let from = if side < 0. { low } else { high };
                                    let pole = std::array::from_fn(|i| from[i] + side * along[i]);
                                    text.triple("CARTESIAN_POINT", pole)
                                });
                                text.add(format!(
                                    "B_SPLINE_CURVE_WITH_KNOTS('',1,({}),.UNSPECIFIED.,.F.,.F.,(2,2),(0.,1.),.PIECEWISE_BEZIER_KNOTS.)",
                                    list(&ends)
                                ))
                            } else {
                                let shift = f64::from(u8::from(tweak.displaced_line == Some(key)));
                                let origin = text.triple("CARTESIAN_POINT", low.map(|x| x + shift));
                                let direction = text.triple("DIRECTION", along);
                                let vector = text.add(format!("VECTOR('',#{direction},1.)"));
                                text.add(format!("LINE('',#{origin},#{vector})"))
                            };
                            let [low, high] = key.map(|vertex| vertices[vertex]);
                            text.add(format!("EDGE_CURVE('',#{low},#{high},#{curve},.T.)"))
                        });
                        uses.push((edge, a < b || tweak.one_way == Some(key)));
                    }
                }
                Ring::Circle {
                    vertex,
                    center,
                    forward,
                } => {
                    let edge = *edges.entry([*vertex; 2]).or_insert_with(|| {
                        let radial = sub(points[*vertex], *center);
                        let axis = text.axis(*center, [0., 0., 1.], unit(radial).unwrap());
                        let circle = text.add(format!("CIRCLE('',#{axis},{:.12})", length(radial)));
                        let vertex = vertices[*vertex];
                        text.add(format!("EDGE_CURVE('',#{vertex},#{vertex},#{circle},.T.)"))
                    });
                    uses.push((edge, *forward));
                }
            }
            if patch.reversed {
                uses.reverse();
            }
            let oriented = uses
                .into_iter()
                .map(|(edge, forward)| {
                    let forward = flag(forward != patch.reversed);
                    text.add(format!("ORIENTED_EDGE('',*,*,#{edge},{forward})"))
                })
                .collect::<Vec<_>>();
            let wire = text.add(format!("EDGE_LOOP('',({}))", list(&oriented)));
            let kind = if index == 0 {
                "FACE_OUTER_BOUND"
            } else {
                "FACE_BOUND"
            };
            bounds.push(text.add(format!("{kind}('',#{wire},{})", flag(!patch.reversed))));
        }
        let (origin, reference) = match &patch.rings[0] {
            Ring::Polygon(ring) => (points[ring[0]], sub(points[ring[1]], points[ring[0]])),
            Ring::Circle { vertex, center, .. } => (*center, sub(points[*vertex], *center)),
        };
        if patch.vertex_loop {
            let pole = text.triple("CARTESIAN_POINT", origin);
            let pole = text.add(format!("VERTEX_POINT('',#{pole})"));
            let wire = text.add(format!("VERTEX_LOOP('',#{pole})"));
            bounds.push(text.add(format!("FACE_BOUND('',#{wire},.T.)")));
        }
        let side = if patch.reversed { -1. } else { 1. };
        let axis = text.axis(
            origin,
            patch.normal.map(|x| x * side),
            unit(reference).unwrap(),
        );
        let plane = text.add(format!("PLANE('',#{axis})"));
        faces.push(text.add(format!(
            "ADVANCED_FACE('',({}),#{plane},{})",
            list(&bounds),
            flag(!patch.reversed)
        )));
    }
    let shell = text.add(format!("CLOSED_SHELL('',({}))", list(&faces)));
    let solid = text.add(format!("MANIFOLD_SOLID_BREP('',#{shell})"));
    let length = text.add("(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.))".into());
    let context = text.add(format!(
        "(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT((#{length})) REPRESENTATION_CONTEXT('',''))"
    ));
    let representation = text.add(format!(
        "ADVANCED_BREP_SHAPE_REPRESENTATION('',(#{solid}),#{context})"
    ));
    let shape = text.add("PRODUCT_DEFINITION_SHAPE('','',$)".into());
    text.add(format!(
        "SHAPE_DEFINITION_REPRESENTATION(#{shape},#{representation})"
    ));
    let unit = [Rational::new(1, 1).unwrap(); 3];
    let header = write_box(unit, [0.; 3], [1., 0., 0.], [0., 1., 0.]).unwrap();
    let mut file = header[..header.find("DATA;").unwrap() + 6].to_owned();
    for (index, record) in text.0.iter().enumerate() {
        file += &format!("#{}={record};\n", index + 1);
    }
    file + "ENDSEC;\nEND-ISO-10303-21;\n"
}

fn points() -> Vec<Point> {
    vec![
        [0., 0., 0.],
        [4., 0., 0.],
        [4., 3., 0.],
        [0., 3., 0.],
        [0., 0., 2.],
        [4., 0., 2.],
        [4., 3., 2.],
        [0., 3., 2.],
        [1., 0., 2.],
        [3., 0., 2.],
        [1., 3., 2.],
        [3., 3., 2.],
        [2.5, 1.5, 2.],
    ]
}

fn patch(normal: Point, ring: &[usize]) -> Patch {
    Patch {
        normal,
        rings: vec![Ring::Polygon(ring.to_vec())],
        reversed: false,
        vertex_loop: false,
    }
}

/// A 4 x 3 x 2 box. The partitioned variant divides the top into three
/// strips and a marker disc, which also splits the long top edges.
fn patches(partitioned: bool) -> Vec<Patch> {
    let z = [0., 0., 1.];
    let center = [2., 1.5, 2.];
    let circle = |forward| Ring::Circle {
        vertex: 12,
        center,
        forward,
    };
    let mut faces = vec![
        patch([0., 0., -1.], &[0, 3, 2, 1]),
        patch([1., 0., 0.], &[1, 2, 6, 5]),
        patch([-1., 0., 0.], &[3, 0, 4, 7]),
    ];
    if !partitioned {
        faces.push(patch([0., -1., 0.], &[0, 1, 5, 4]));
        faces.push(patch([0., 1., 0.], &[2, 3, 7, 6]));
        faces.push(patch(z, &[4, 5, 6, 7]));
        return faces;
    }
    faces.push(patch([0., -1., 0.], &[0, 1, 5, 9, 8, 4]));
    faces.push(patch([0., 1., 0.], &[2, 3, 7, 10, 11, 6]));
    faces.push(Patch {
        reversed: true,
        ..patch(z, &[4, 8, 10, 7])
    });
    let mut middle = patch(z, &[8, 9, 11, 10]);
    middle.rings.push(circle(false));
    faces.push(middle);
    faces.push(patch(z, &[9, 5, 6, 11]));
    faces.push(Patch {
        rings: vec![circle(true)],
        ..patch(z, &[])
    });
    faces
}

fn count(entities: &BTreeMap<u32, Vec<Entity<'_>>>, kind: &str) -> usize {
    let named = entities.values().filter(|parts| parts[0].name == kind);
    named.count()
}

/// Merge the single solid and report whether the arena changed, with its
/// face, edge and vertex counts.
fn merged(source: &str) -> (bool, [usize; 3]) {
    let mut entities = parse_entities(source, TOL).unwrap();
    let solids = Shell::new(&entities, TOL, true)
        .represented_solids(true)
        .unwrap();
    let changed = unify(&mut entities, &solids, TOL);
    let counts = ["ADVANCED_FACE", "EDGE_CURVE", "VERTEX_POINT"].map(|kind| count(&entities, kind));
    (changed, counts)
}

#[test]
fn a_partitioned_top_imports_as_the_same_box() {
    let plain = step(&points(), &patches(false), Tweak::default());
    let split = step(&points(), &patches(true), Tweak::default());
    // The marker circle joins two faces of one plane.
    assert_eq!(merged(&plain), (false, [6, 12, 8]));
    assert_eq!(merged(&split), (true, [6, 12, 8]));
    let (expected, body) = read_import_solid(&plain, TOL).unwrap();
    let (proof, merged) = read_import_solid(&split, TOL).unwrap();
    assert_eq!(expected.volume_mm3, 24.);
    assert_eq!(proof.volume_mm3, 24.);
    assert_eq!(proof.bounds_mm, [[0.; 3], [4., 3., 2.]]);
    // The exact rational prism of the same extents is the mass oracle.
    let [length, width, height] = [4, 3, 2].map(|extent| Rational::new(extent, 1).unwrap());
    let exact = crate::prism::Prism::rectangular(1, [length, width, height]).unwrap();
    assert_eq!(exact.volume().as_f64(), proof.volume_mm3);
    assert_eq!(exact.area().as_f64() * 2., proof.volume_mm3);
    assert_eq!([proof.face_count, proof.edges_mm.len()], [6, 12]);
    assert!(proof.circles.is_empty());
    assert_eq!(
        std::mem::discriminant(&merged.profile),
        std::mem::discriminant(&body.profile)
    );
    assert_eq!(read_analytic_solid(&split, TOL).unwrap(), proof);
    assert_eq!(read_planar_solid(&split, TOL).unwrap(), proof);
    assert_eq!(read_source_compound(&split, TOL).unwrap().len(), 1);
}

#[test]
fn irregular_partitions_keep_their_transmitted_topology_and_verdict() {
    type Edit<'a> = &'a dyn Fn(&mut Vec<Point>, &mut Vec<Patch>, &mut Tweak);
    let build = |edit: Edit<'_>| {
        let (mut points, mut faces, mut tweak) = (points(), patches(true), Tweak::default());
        edit(&mut points, &mut faces, &mut tweak);
        step(&points, &faces, tweak)
    };
    // A bound without edges stays on its rewritten face, which has no plane.
    let pole = build(&|_, faces, _| faces[3].vertex_loop = true);
    assert_eq!(merged(&pole), (true, [6, 12, 9]));
    assert!(read_import_solid(&pole, TOL).is_err());
    // An edge that fails the edge reader is never merged.
    let displaced = build(&|_, _, tweak| tweak.displaced_line = Some([8, 9]));
    assert_eq!(merged(&displaced), (true, [6, 14, 10]));
    assert_eq!(
        read_import_solid(&displaced, TOL).err(),
        Some(ReadError::Invalid)
    );
    // A bent chain keeps its vertex.
    let bent = build(&|points, _, _| points[8][1] = 0.25);
    assert_eq!(merged(&bent), (true, [6, 14, 10]));
    assert!(read_import_solid(&bent, TOL).is_err());
    // Two coplanar faces that use their common edge in one direction are
    // not a partition.
    let one_way = build(&|_, _, tweak| tweak.one_way = Some([8, 10]));
    assert_eq!(merged(&one_way), (false, [9, 19, 13]));
    assert_eq!(
        read_import_solid(&one_way, TOL).err(),
        Some(ReadError::Invalid)
    );
}

#[test]
fn regular_families_are_screened_independently_of_face_order() {
    let split = step(&points(), &patches(true), Tweak::default());
    let entities = parse_entities(&split, TOL).unwrap();
    let shell = Shell::new(&entities, TOL, true);
    let faces = shell.shell_faces().unwrap();
    assert_eq!(shell.outside_regular_families(&faces), Ok(()));
    // An unreadable face is left to the member reader.
    assert_eq!(shell.outside_regular_families(&[1]), Ok(()));
    // Four faces on one edge in balanced directions are a pinched contact.
    let pinched = [faces[5], faces[6], faces[5], faces[6]];
    assert_eq!(
        shell.outside_regular_families(&pinched),
        Err(ReadError::Unsupported)
    );
    let mut writer = crate::step::Writer::new();
    let sphere = crate::sphere::Sphere::new(2.).unwrap();
    writer
        .append_sphere(sphere, [0.; 3], [1., 0., 0.], [0., 1., 0.])
        .unwrap();
    let source = writer.finish().unwrap();
    let entities = parse_entities(&source, TOL).unwrap();
    let shell = Shell::new(&entities, TOL, true);
    let faces = shell.shell_faces().unwrap();
    assert_eq!(
        shell.outside_regular_families(&faces),
        Err(ReadError::Unsupported)
    );
    assert_eq!(
        Shell::new(&entities, TOL, false).outside_regular_families(&faces),
        Ok(())
    );
}

#[test]
fn an_edge_may_trim_a_longer_straight_spline() {
    let plain = step(&points(), &patches(false), Tweak::default());
    let tweak = Tweak {
        spline: Some([4, 5]),
        ..Tweak::default()
    };
    let trimmed = step(&points(), &patches(false), tweak);
    assert_ne!(trimmed, plain);
    assert_eq!(
        read_planar_solid(&trimmed, TOL).unwrap(),
        read_planar_solid(&plain, TOL).unwrap()
    );
    // A vertex beyond the transmitted segment is not on the curve.
    let mut far = points();
    far[5][0] = 5.5;
    far[1][0] = 5.5;
    far[2][0] = 5.5;
    far[6][0] = 5.5;
    let beyond = step(&far, &patches(false), tweak);
    assert!(read_planar_solid(&beyond, TOL).is_ok());
    let short = beyond.replacen("(6.5", "(5.25", 1);
    assert_ne!(short, beyond);
    assert_eq!(read_planar_solid(&short, TOL), Err(ReadError::Invalid));
}
