//! Independent expectations for straight sketch extrusions read directly
//! from a saved document: the sketch plane, the saved (solved) sketch point
//! positions, the profile region and the one-sided extrusion distance. No
//! geometry kernel is involved.
#![allow(dead_code)]
use multicad_core::document::{Document, EntityData, EntityId};
use multicad_precision::PortableMath as _;
use serde_json::Value;

/// A straight prism: a planar polygon swept along `sweep`.
#[derive(Clone, Debug)]
pub struct Prism {
    /// The profile loop in world millimeters at the sketch plane.
    pub loop_mm: Vec<[f64; 3]>,
    /// The profile loop in sketch (u, v) coordinates.
    pub loop_uv: Vec<[f64; 2]>,
    /// Per edge: the (u, v) center of a circular-arc edge (the minor arc
    /// between its ends), or none for a line.
    pub arcs: Vec<Option<[f64; 2]>>,
    pub origin: [f64; 3],
    pub u: [f64; 3],
    pub v: [f64; 3],
    pub n: [f64; 3],
    /// The world sweep vector: sketch normal times the signed distance.
    pub sweep: [f64; 3],
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn number(value: &Value) -> f64 {
    value
        .as_f64()
        .unwrap_or_else(|| panic!("expected a number, found {value}"))
}

fn literal_mm(value: &Value) -> f64 {
    assert!(
        value["unit"] == "mm" || value["unit"] == "deg",
        "expected a literal dimension, found {value}"
    );
    number(&value["value"])
}

/// Read a one-sided, untapered, normal extrusion of a polygon sketch region.
pub fn extruded_prism(document: &Document, feature: EntityId) -> Prism {
    let EntityData::Feature { definition, .. } = &document.get(feature).unwrap().data else {
        panic!("{feature} is not a feature");
    };
    assert_eq!(definition.type_id, "model.extrude");
    let data = serde_json::to_value(definition).unwrap()["data"].clone();
    let parameters = &data["parameters"];
    assert_eq!(literal_mm(&parameters["direction"]["azimuth"]), 0.);
    assert_eq!(literal_mm(&parameters["direction"]["elevation"]), 90.);
    assert_eq!(parameters["direction"]["reverse"], false);
    assert_eq!(literal_mm(&parameters["start_offset"]), 0.);
    assert_eq!(literal_mm(&parameters["taper"]), 0.);
    assert_eq!(parameters["extent"]["kind"], "one_sided");
    assert_eq!(parameters["extent"]["end"]["kind"], "distance");
    let distance: multicad_core::constraint_sketch::Dimension =
        serde_json::from_value(parameters["extent"]["end"]["distance"].clone()).unwrap();
    let evaluated = multicad_core::parameters::evaluate(document).unwrap();
    let distance = distance
        .value_in(
            multicad_core::parameters::Unit::Millimeter,
            Some(&evaluated),
        )
        .unwrap();
    let profile = &data["profile"];
    let region = sketch_region(document, profile);
    assert!(region.holes.is_empty());
    let Boundary::Loop { vertices, arcs } = region.outer.clone() else {
        panic!("a circular extrusion profile");
    };
    let loop_mm = vertices.iter().map(|p| region.world(*p)).collect();
    Prism {
        loop_mm,
        loop_uv: vertices,
        arcs,
        origin: region.origin,
        u: region.u,
        v: region.v,
        n: region.n,
        sweep: region.n.map(|x| x * distance),
    }
}

/// One closed boundary of a saved sketch region, in sketch (u, v)
/// coordinates: a loop of lines and minor circular arcs, or a full circle.
#[derive(Clone, Debug)]
pub enum Boundary {
    Loop {
        vertices: Vec<[f64; 2]>,
        /// Per edge (vertex i to i + 1): the arc center, or none for a line.
        arcs: Vec<Option<[f64; 2]>>,
    },
    Circle {
        /// The saved sketch curve identity.
        id: String,
        center: [f64; 2],
        radius: f64,
    },
}

/// A saved sketch region on its resolved plane.
#[derive(Clone, Debug)]
pub struct Region {
    pub origin: [f64; 3],
    pub u: [f64; 3],
    pub v: [f64; 3],
    pub n: [f64; 3],
    pub outer: Boundary,
    pub holes: Vec<Boundary>,
}

impl Region {
    pub fn world(&self, p: [f64; 2]) -> [f64; 3] {
        std::array::from_fn(|i| self.origin[i] + p[0] * self.u[i] + p[1] * self.v[i])
    }
}

/// Read a profile `{sketch, region}` reference from the saved document.
pub fn sketch_region(document: &Document, profile: &Value) -> Region {
    let sketch = EntityId::parse(profile["sketch"].as_str().unwrap()).unwrap();
    let EntityData::Sketch {
        plane, definition, ..
    } = &document.get(sketch).unwrap().data
    else {
        panic!("{sketch} is not a sketch");
    };
    let sketch_data = serde_json::to_value(definition).unwrap()["data"].clone();
    let local_u: [f64; 3] = std::array::from_fn(|i| number(&sketch_data["u"][i]));
    // An unattached plane takes the sketch's world U directly; an offset
    // plane chain stores each distance as origin Z and the sketch U in the
    // unattached base plane's frame.
    let EntityData::Plane {
        origin_mm,
        normal,
        anchor,
    } = &document.get(*plane).unwrap().data
    else {
        panic!("{plane} is not a plane");
    };
    let (origin_mm, u, v, n) = if anchor.is_none() {
        let n = *normal;
        (*origin_mm, local_u, cross(n, local_u), n)
    } else {
        let mut current = *plane;
        let mut offset = 0.;
        let base = loop {
            let EntityData::Plane {
                origin_mm,
                normal,
                anchor,
            } = &document.get(current).unwrap().data
            else {
                panic!("{current} is not a plane");
            };
            match anchor {
                Some(link) if link.role == "offset_plane" => {
                    assert!(link.geometry.is_none());
                    assert_eq!(*normal, [0., 0., 1.]);
                    assert_eq!([origin_mm[0], origin_mm[1]], [0., 0.]);
                    offset += origin_mm[2];
                    current = link.target;
                }
                None => {
                    break multicad_core::sketch::inspect_plane(document, current, 0, None)
                        .unwrap();
                }
                Some(link) => panic!("attached plane {current} ({}) needs evaluation", link.role),
            }
        };
        (
            std::array::from_fn(|i| base.origin[i] + base.n[i] * offset),
            base.vector(local_u),
            base.vector(cross([0., 0., 1.], local_u)),
            base.n,
        )
    };
    assert!((0..3).all(|i| (cross(u, v)[i] - n[i]).abs() < 1e-12));
    let points = sketch_data["sketch"]["points"].as_array().unwrap();
    let curves = sketch_data["sketch"]["curves"].as_array().unwrap();
    let position = |id: &str| -> [f64; 2] {
        let point = points.iter().find(|point| point["id"] == id).unwrap();
        [
            number(&point["position_mm"][0]),
            number(&point["position_mm"][1]),
        ]
    };
    let chain = |ids: &Vec<Value>| -> (Vec<[f64; 2]>, Vec<Option<[f64; 2]>>) {
        let mut segments: Vec<(String, String, Option<[f64; 2]>)> = ids
            .iter()
            .map(|id| {
                let curve = curves.iter().find(|curve| curve["id"] == *id).unwrap();
                assert_eq!(curve["construction"], false);
                let center = match curve["kind"]["kind"].as_str().unwrap() {
                    "line" => None,
                    "arc" => Some(position(curve["kind"]["center"].as_str().unwrap())),
                    other => panic!("unsupported profile curve {other}"),
                };
                (
                    curve["kind"]["start"].as_str().unwrap().to_owned(),
                    curve["kind"]["end"].as_str().unwrap().to_owned(),
                    center,
                )
            })
            .collect();
        // Chain the region's curves into one closed loop; edge i runs from
        // vertex i to vertex i + 1.
        let (first, mut next, center) = segments.remove(0);
        let mut order = vec![first.clone()];
        let mut arcs = vec![center];
        while next != first {
            order.push(next.clone());
            let index = segments
                .iter()
                .position(|(a, b, _)| *a == next || *b == next)
                .expect("an open profile");
            let (a, b, center) = segments.remove(index);
            arcs.push(center);
            next = if a == next { b } else { a };
        }
        assert!(segments.is_empty());
        (order.iter().map(|id| position(id)).collect(), arcs)
    };
    let boundary = |curves_json: &Value| -> Boundary {
        let ids = curves_json.as_array().unwrap();
        if let [id] = ids.as_slice() {
            let curve = curves.iter().find(|curve| curve["id"] == *id).unwrap();
            if curve["kind"]["kind"] == "circle" {
                return Boundary::Circle {
                    id: curve["id"].as_str().unwrap().to_owned(),
                    center: position(curve["kind"]["center"].as_str().unwrap()),
                    radius: number(&curve["kind"]["radius_mm"]),
                };
            }
        }
        let (vertices, arcs) = chain(ids);
        Boundary::Loop { vertices, arcs }
    };
    Region {
        origin: origin_mm,
        u,
        v,
        n,
        outer: boundary(&profile["region"]["outer"]["curves"]),
        holes: profile["region"]["holes"]
            .as_array()
            .unwrap()
            .iter()
            .map(|hole| boundary(&hole["curves"]))
            .collect(),
    }
}

impl Boundary {
    fn shoelace(vertices: &[[f64; 2]]) -> f64 {
        (0..vertices.len())
            .map(|i| {
                let [a, b] = [vertices[i], vertices[(i + 1) % vertices.len()]];
                a[0] * b[1] - b[0] * a[1]
            })
            .sum::<f64>()
            / 2.
    }
    /// One arc edge of a loop: its radius, minor sweep angle, and whether it
    /// bulges out of the chord polygon (center on the interior side).
    fn arc(&self, i: usize) -> Option<(f64, f64, bool)> {
        let Self::Loop { vertices, arcs } = self else {
            return None;
        };
        let c = arcs[i]?;
        let [a, b] = [vertices[i], vertices[(i + 1) % vertices.len()]];
        let [ra, rb] = [a, b].map(|p| (p[0] - c[0]).portable_hypot(p[1] - c[1]));
        assert!((ra - rb).abs() < 1e-9, "arc ends off its circle");
        let dot = (a[0] - c[0]) * (b[0] - c[0]) + (a[1] - c[1]) * (b[1] - c[1]);
        let theta = (dot / (ra * rb)).clamp(-1., 1.).portable_acos();
        let side = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
        Some((ra, theta, (side > 0.) == (Self::shoelace(vertices) > 0.)))
    }
    fn edges(&self) -> usize {
        match self {
            Self::Loop { vertices, .. } => vertices.len(),
            Self::Circle { .. } => 0,
        }
    }
    /// Enclosed area: the chord polygon adjusted by each circular segment
    /// r² (θ - sin θ) / 2.
    pub fn area_mm2(&self) -> f64 {
        match self {
            Self::Circle { radius, .. } => std::f64::consts::PI * radius * radius,
            Self::Loop { vertices, .. } => {
                let mut area = Self::shoelace(vertices).abs();
                for i in 0..self.edges() {
                    if let Some((r, theta, outward)) = self.arc(i) {
                        let segment = r * r * (theta - theta.portable_sin()) / 2.;
                        area += if outward { segment } else { -segment };
                    }
                }
                area
            }
        }
    }
    pub fn perimeter_mm(&self) -> f64 {
        match self {
            Self::Circle { radius, .. } => 2. * std::f64::consts::PI * radius,
            Self::Loop { vertices, .. } => (0..vertices.len())
                .map(|i| match self.arc(i) {
                    Some((r, theta, _)) => r * theta,
                    None => {
                        let [a, b] = [vertices[i], vertices[(i + 1) % vertices.len()]];
                        (b[0] - a[0]).portable_hypot(b[1] - a[1])
                    }
                })
                .sum(),
        }
    }
    /// (u, v) bounds, including the axis extremes reached inside each arc.
    pub fn bounds_uv(&self) -> [[f64; 2]; 2] {
        let mut b = [[f64::INFINITY; 2], [f64::NEG_INFINITY; 2]];
        let mut add = |p: [f64; 2]| {
            for k in 0..2 {
                b[0][k] = b[0][k].min(p[k]);
                b[1][k] = b[1][k].max(p[k]);
            }
        };
        match self {
            Self::Circle { center, radius, .. } => {
                add([center[0] - radius, center[1] - radius]);
                add([center[0] + radius, center[1] + radius]);
            }
            Self::Loop { vertices, arcs } => {
                for (i, p) in vertices.iter().enumerate() {
                    add(*p);
                    let Some(c) = arcs[i] else { continue };
                    let [a, q] = [*p, vertices[(i + 1) % vertices.len()]];
                    let r = (a[0] - c[0]).portable_hypot(a[1] - c[1]);
                    // An axis extreme of the circle lies on the minor arc when
                    // it is on the arc's side of the chord.
                    for d in [[1., 0.], [-1., 0.], [0., 1.], [0., -1.]] {
                        let e = [c[0] + r * d[0], c[1] + r * d[1]];
                        let chord = |p: [f64; 2]| {
                            (q[0] - a[0]) * (p[1] - a[1]) - (q[1] - a[1]) * (p[0] - a[0])
                        };
                        let mid = [
                            c[0] - chord_offset(a, q, c)[0],
                            c[1] - chord_offset(a, q, c)[1],
                        ];
                        if chord(e) * chord(mid) > 0. && chord(e).abs() > 0. {
                            add(e);
                        }
                    }
                }
            }
        }
        b
    }
}

/// The vector from the arc midpoint to the center, used to tell the minor
/// arc's side of its chord: the midpoint is the center moved away from the
/// chord midpoint direction.
fn chord_offset(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> [f64; 2] {
    let m = [(a[0] + b[0]) / 2. - c[0], (a[1] + b[1]) / 2. - c[1]];
    let length = m[0].portable_hypot(m[1]);
    let r = (a[0] - c[0]).portable_hypot(a[1] - c[1]);
    [-m[0] / length * r, -m[1] / length * r]
}

impl Region {
    /// Net area: outer less holes.
    pub fn area_mm2(&self) -> f64 {
        self.outer.area_mm2() - self.holes.iter().map(Boundary::area_mm2).sum::<f64>()
    }
    pub fn perimeter_mm(&self) -> f64 {
        self.outer.perimeter_mm() + self.holes.iter().map(Boundary::perimeter_mm).sum::<f64>()
    }
    /// World bounds of the region swept by `sweep`; u and v must be world axes.
    pub fn swept_bounds(&self, sweep: [f64; 3]) -> [[f64; 3]; 2] {
        for d in [self.u, self.v] {
            assert_eq!(
                d.iter().filter(|x| **x != 0.).count(),
                1,
                "axis-aligned sketch"
            );
        }
        let [lo, hi] = self.outer.bounds_uv();
        let mut b = [[f64::INFINITY; 3], [f64::NEG_INFINITY; 3]];
        for p in [lo, [lo[0], hi[1]], [hi[0], lo[1]], hi] {
            let w = self.world(p);
            for q in [w, std::array::from_fn(|i| w[i] + sweep[i])] {
                for i in 0..3 {
                    b[0][i] = b[0][i].min(q[i]);
                    b[1][i] = b[1][i].max(q[i]);
                }
            }
        }
        b
    }
}

impl Prism {
    fn outer(&self) -> Boundary {
        Boundary::Loop {
            vertices: self.loop_uv.clone(),
            arcs: self.arcs.clone(),
        }
    }
    fn region(&self) -> Region {
        Region {
            origin: self.origin,
            u: self.u,
            v: self.v,
            n: self.n,
            outer: self.outer(),
            holes: vec![],
        }
    }
    pub fn area_mm2(&self) -> f64 {
        self.outer().area_mm2()
    }
    pub fn length_mm(&self) -> f64 {
        self.sweep.iter().map(|x| x * x).sum::<f64>().sqrt()
    }
    pub fn volume_mm3(&self) -> f64 {
        self.area_mm2() * self.length_mm()
    }
    pub fn world(&self, p: [f64; 2]) -> [f64; 3] {
        self.region().world(p)
    }
    pub fn bounds(&self) -> [[f64; 3]; 2] {
        self.region().swept_bounds(self.sweep)
    }
}
