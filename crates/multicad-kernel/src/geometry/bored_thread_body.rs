//! One exact internal thread band buried in an existing coaxial through bore
//! of certified orthogonal stock.
//!
//! The retained stock keeps every original face, cut, edge and vertex. Only
//! the selected bore wall is split: its two analytic collars stay at the
//! original radius, the thread tool's planar caps become steps, and the
//! exact helical bands come from the retained thread tool itself.
use super::helical_body::stock_frame;
use super::*;
use crate::HelicalThread;
use crate::rust_prism_topology::position;
use multicad_core::references::{
    AnalyticGeometry, AnchoredFrame, Origin, ReferenceFrame, TopologyElement, TopologyKey,
    TopologyKind,
};
use multicad_rust_kernel::{
    display::{BandWall, ThreadedThroughBore, threaded_bored_orthogonal},
    helical_thread::Thread,
    orthogonal_pockets::OrthogonalThroughBore,
};
use std::collections::BTreeSet;
use std::f64::consts::{PI, TAU};

type Point = [f64; 3];
const OPERATION: &str = "Rust bored internal thread";

pub(crate) struct BoredThreadBody {
    pub(super) feature: EntityId,
    /// The original orthogonal stock solid, sharing its retained body.
    stock: Solid,
    /// The exact plain helical tool, placed in the same world frame.
    tool: Solid,
    wall: Wall,
    bore: usize,
    /// Tool origin and U/V axes in stock coordinates.
    local: (Point, [Point; 2]),
    /// Original wall origins of each bore when the stock was a bored prism.
    inherited: Rc<[Origin]>,
}

/// The retained band law: exact helical material, or a display-only smooth
/// envelope at the internal crest radius.
#[derive(Clone, Copy)]
enum Wall {
    Thread(Thread),
    /// Crest radius and tool length.
    Cosmetic(f64, f64),
}
/// The authored band before placement.
#[derive(Clone, Copy)]
enum Band {
    Thread(HelicalThread, Option<f64>),
    Cosmetic {
        crest_radius_mm: f64,
        major_radius_mm: f64,
        length_mm: f64,
    },
}
impl Band {
    fn major(self) -> f64 {
        match self {
            Self::Thread(spec, _) => spec.major_radius_mm,
            Self::Cosmetic {
                major_radius_mm, ..
            } => major_radius_mm,
        }
    }
    fn length(self) -> f64 {
        match self {
            Self::Thread(spec, _) => spec.length_mm,
            Self::Cosmetic { length_mm, .. } => length_mm,
        }
    }
}

fn unsupported(message: &str) -> GeometryError {
    error(OPERATION, GeometryErrorKind::Unsupported, message)
}
fn unit(v: Point) -> Point {
    let length = dot(v, v).sqrt();
    v.map(|x| x / length)
}

impl Solid {
    /// Cut an internal thread whose complete axial interval lies inside an
    /// existing coaxial through bore. The bore must lie inside the thread
    /// root, and the protected collar around the crest must be original
    /// material untouched by any other cut.
    pub fn rust_cut_internal_thread_in_bore(
        &self,
        spec: HelicalThread,
        root_radius_mm: Option<f64>,
        minimum_wall_mm: f64,
        feature: EntityId,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        self.cut_band_in_bore(
            Band::Thread(spec, root_radius_mm),
            minimum_wall_mm,
            feature,
            frame,
            policy,
        )
    }

    /// The explicitly cosmetic receiving envelope of the same buried thread:
    /// a smooth coaxial band at its internal crest radius. It never earns
    /// modeled-thread manufacturing or STEP exchange.
    #[allow(clippy::too_many_arguments)]
    pub fn rust_cut_cosmetic_internal_thread_in_bore(
        &self,
        crest_radius_mm: f64,
        major_radius_mm: f64,
        length_mm: f64,
        minimum_wall_mm: f64,
        feature: EntityId,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        for value in [crest_radius_mm, major_radius_mm, length_mm] {
            policy.length(value)?;
        }
        if crest_radius_mm >= major_radius_mm {
            return Err(error(
                OPERATION,
                GeometryErrorKind::Input,
                "The cosmetic crest must lie inside the thread major radius",
            ));
        }
        self.cut_band_in_bore(
            Band::Cosmetic {
                crest_radius_mm,
                major_radius_mm,
                length_mm,
            },
            minimum_wall_mm,
            feature,
            frame,
            policy,
        )
    }

    fn cut_band_in_bore(
        &self,
        band: Band,
        minimum_wall_mm: f64,
        feature: EntityId,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        if self.orthogonal_body().is_none()
            && let Some(stock) = self.bored_prism_orthogonal_stock()
        {
            // The prism's bores are redrilled in their original loop order.
            let (prism, _) = self.rust_prism_frame().expect("bored prism");
            let inherited = prism
                .inner_loops()
                .map(|(source, _)| Origin {
                    source: EntityId::from_uuid_bits(source),
                    role: "bore.wall".into(),
                })
                .collect();
            return stock?.cut_internal_thread_in_bore(
                band,
                minimum_wall_mm,
                feature,
                frame,
                policy,
                inherited,
            );
        }
        self.cut_internal_thread_in_bore(
            band,
            minimum_wall_mm,
            feature,
            frame,
            policy,
            Rc::from([]),
        )
    }

    fn cut_internal_thread_in_bore(
        &self,
        band: Band,
        minimum_wall_mm: f64,
        feature: EntityId,
        frame: PlaneFrame,
        policy: GeometryPolicy,
        inherited: Rc<[Origin]>,
    ) -> Result<Self, GeometryError> {
        policy.validate()?;
        frame.checked(policy)?;
        let Some((body, stock_frame)) = self.orthogonal_body() else {
            return Err(unsupported(
                "A buried receiving thread requires certified orthogonal stock",
            ));
        };
        let (pockets, bores, _, slots) = body.round_cuts();
        if self.policy != policy
            || policy.boolean_fuzzy_mm != 0.
            || body.is_imported()
            || slots != 0
            || body.feature == feature
            || !(inherited.is_empty() || inherited.len() == bores.len())
        {
            return Err(unsupported(
                "A buried receiving thread requires distinct native drilled stock under the same exact policy",
            ));
        }
        if !minimum_wall_mm.is_finite() || minimum_wall_mm < 0. {
            return Err(error(
                OPERATION,
                GeometryErrorKind::Input,
                "Minimum receiving wall must be finite and non-negative",
            ));
        }
        let tolerance = policy.linear_tolerance_mm;
        let chart = stock_frame::Chart::new(stock_frame);
        let origin = chart.point(frame.origin);
        let normal = unit(chart.vector(frame.normal()));
        // The single bore whose analytic axis carries the requested tool axis.
        let mut matches = bores.iter().enumerate().filter_map(|(index, bore)| {
            let [c0, c1] = bore.centers_mm;
            let depth = bore.depth_mm();
            let d: Point = std::array::from_fn(|i| (c1[i] - c0[i]) / depth);
            let offset: Point = std::array::from_fn(|i| origin[i] - c0[i]);
            let s = dot(offset, d);
            let radial: Point = std::array::from_fn(|i| offset[i] - s * d[i]);
            ((dot(normal, d).abs() - 1.).abs() <= 1e-10 && dot(radial, radial).sqrt() <= tolerance)
                .then_some((index, *bore, d, s))
        });
        let Some((index, bore, d, s)) = matches.next() else {
            return Err(unsupported(
                "A buried receiving thread requires one existing coaxial analytic through bore",
            ));
        };
        if matches.next().is_some() {
            return Err(unsupported(
                "The receiving thread axis is shared by more than one through bore",
            ));
        }
        let [c0, _] = bore.centers_mm;
        let depth = bore.depth_mm();
        let along = dot(normal, d).signum();
        let end = s + along * band.length();
        let [low, high] = [s.min(end), s.max(end)];
        if low < policy.minimum_feature_mm || depth - high < policy.minimum_feature_mm {
            return Err(unsupported(
                "A buried receiving thread must leave a finite original collar at both bore mouths",
            ));
        }
        // Snap onto the exact bore axis and keep the requested clock.
        let u_local = chart.vector(frame.u);
        let u_local = unit(std::array::from_fn(|i| u_local[i] - dot(u_local, d) * d[i]));
        let n_local = d.map(|x| along * x);
        let v_local = cross(n_local, u_local);
        let origin_local: Point = std::array::from_fn(|i| c0[i] + s * d[i]);
        let placed = PlaneFrame {
            origin: position(stock_frame, origin_local),
            u: {
                let [x, y, z] = u_local;
                std::array::from_fn(|i| {
                    x * stock_frame.u[i] + y * stock_frame.v[i] + z * stock_frame.normal()[i]
                })
            },
            v: {
                let [x, y, z] = v_local;
                std::array::from_fn(|i| {
                    x * stock_frame.u[i] + y * stock_frame.v[i] + z * stock_frame.normal()[i]
                })
            },
        };
        stock_frame::certify_tool_domain(frame, placed, band.major(), band.length(), policy)?;
        let (tool, wall, inner) = match band {
            Band::Thread(spec, root) => {
                let tool = Self::rust_helical_thread(spec, root, feature, placed, policy)?;
                let (thread, _) = tool
                    .rust_helical_thread_data()
                    .expect("plain helical thread tool");
                (tool, Wall::Thread(thread), thread.root_radius_mm())
            }
            Band::Cosmetic {
                crest_radius_mm,
                length_mm,
                ..
            } => (
                Self::cosmetic_band_tool(crest_radius_mm, length_mm, feature, placed, policy)?,
                Wall::Cosmetic(crest_radius_mm, length_mm),
                crest_radius_mm,
            ),
        };
        if inner - bore.radius_mm < policy.minimum_feature_mm {
            return Err(unsupported(
                "The existing bore must lie inside the receiving thread root",
            ));
        }
        let protected = band.major() + minimum_wall_mm.max(policy.minimum_feature_mm);
        policy.length(protected)?;
        // The protected collar prism around the thread interval.
        let axis = d
            .iter()
            .position(|x| x.abs() == 1.)
            .ok_or_else(|| unsupported("The receiving bore axis is not an exact stock axis"))?;
        let mut zone = [[0.; 3]; 2];
        for i in 0..3 {
            let [a, b] = if i == axis {
                [c0[i] + d[i] * low, c0[i] + d[i] * high]
            } else {
                [origin_local[i] - protected, origin_local[i] + protected]
            };
            zone[0][i] = a.min(b);
            zone[1][i] = a.max(b);
        }
        let overlap =
            |min: Point, max: Point| (0..3).all(|i| min[i] < zone[1][i] && zone[0][i] < max[i]);
        let covered: f64 = body
            .geometry
            .cells()
            .iter()
            .map(|cell| {
                (0..3)
                    .map(|i| (cell.max[i].min(zone[1][i]) - cell.min[i].max(zone[0][i])).max(0.))
                    .product::<f64>()
            })
            .sum();
        let volume: f64 = (0..3).map(|i| zone[1][i] - zone[0][i]).product();
        let cuts_clear = bores
            .iter()
            .enumerate()
            .filter(|&(i, _)| i != index)
            .all(|(_, other)| {
                let [a, b] = other.centers_mm;
                !overlap(
                    std::array::from_fn(|i| a[i].min(b[i]) - other.radius_mm),
                    std::array::from_fn(|i| a[i].max(b[i]) + other.radius_mm),
                )
            })
            && pockets.iter().all(|pocket| {
                let normal = body.boundary.faces[pocket.face].normal;
                let floor = pocket.floor_center_mm(normal);
                let [a, b] = [pocket.center_mm, floor];
                !overlap(
                    std::array::from_fn(|i| a[i].min(b[i]) - pocket.radius_mm),
                    std::array::from_fn(|i| a[i].max(b[i]) + pocket.radius_mm),
                )
            });
        if covered < volume * (1. - 1e-12) || !cuts_clear {
            return Err(error(
                OPERATION,
                GeometryErrorKind::Input,
                "The receiving stock has insufficient original wall around the thread crest",
            ));
        }
        let result = BoredThreadBody {
            feature,
            stock: self.share()?,
            tool,
            wall,
            bore: index,
            local: (origin_local, [u_local, v_local]),
            inherited,
        };
        Ok(Self {
            body: SolidBody::BoredThread(Rc::new(result)),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        })
    }

    fn cosmetic_band_tool(
        radius: f64,
        length: f64,
        feature: EntityId,
        placed: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        use multicad_rust_kernel::{Circle, Rational, cylinder::AxialCylinder};
        let exact = |value: f64| {
            Rational::parse_decimal(&value.to_string())
                .map_err(|_| unsupported("Cosmetic band dimensions are not exact decimals"))
        };
        let zero = Rational::new(0, 1).expect("zero is exact");
        let circle = Circle::new([zero; 2], exact(radius)?)
            .map_err(|_| unsupported("Invalid cosmetic band circle"))?;
        let cylinder = AxialCylinder::new(feature.uuid_bits(), circle, zero, exact(length)?)
            .map_err(|_| unsupported("Invalid cosmetic band cylinder"))?;
        Self::from_axial_cylinder_sweep(
            cylinder,
            placed,
            [0.; 2],
            false,
            Origin {
                source: feature,
                role: "thread.cosmetic.envelope".into(),
            },
            policy,
        )
    }

    pub(crate) fn bored_thread_body(&self) -> Option<&BoredThreadBody> {
        match &self.body {
            SolidBody::BoredThread(body) => Some(body),
            _ => None,
        }
    }
}

impl BoredThreadBody {
    fn stock_parts(&self) -> (&orthogonal_body::OrthogonalBody, PlaneFrame) {
        self.stock
            .orthogonal_body()
            .expect("retained orthogonal stock")
    }
    fn bore(&self) -> OrthogonalThroughBore {
        self.stock_parts().0.round_cuts().1[self.bore]
    }
    /// Axis, bore-local tool start and the sorted collar interval.
    fn axial(&self) -> (Point, f64, f64, [f64; 2]) {
        let bore = self.bore();
        let [c0, c1] = bore.centers_mm;
        let depth = bore.depth_mm();
        let d: Point = std::array::from_fn(|i| (c1[i] - c0[i]) / depth);
        let (origin, [u, v]) = self.local;
        let s = dot(std::array::from_fn(|i| origin[i] - c0[i]), d);
        let along = dot(cross(u, v), d).signum();
        let end = s + along * self.length();
        (d, s, along, [s.min(end), s.max(end)])
    }
    fn local_point(&self, s: f64) -> Point {
        let (d, ..) = self.axial();
        let c0 = self.bore().centers_mm[0];
        std::array::from_fn(|i| c0[i] + s * d[i])
    }
    fn world(&self, local: Point) -> Point {
        position(self.stock_parts().1, local)
    }
    fn stock_face_count(&self) -> usize {
        let (body, _) = self.stock_parts();
        body.faces(PlaneFrame::default()).len()
    }
    /// Bore walls are the last stock faces when there are no slots.
    fn wall_face(&self) -> usize {
        let (body, _) = self.stock_parts();
        self.stock_face_count() - body.round_cuts().1.len() + self.bore
    }

    /// Stock faces keep their original loops and edge ordinals. Steps are
    /// bounded by helical rims, which have no exact planar line/arc loop.
    pub(crate) fn planar_loops(&self, face: u32) -> Vec<Vec<(u32, bool)>> {
        if (face as usize) < self.stock_face_count() {
            self.stock_parts().0.loops(face)
        } else {
            vec![]
        }
    }

    fn length(&self) -> f64 {
        match self.wall {
            Wall::Thread(thread) => thread.spec().length_mm,
            Wall::Cosmetic(_, length) => length,
        }
    }
    pub(crate) fn is_cosmetic(&self) -> bool {
        matches!(self.wall, Wall::Cosmetic(..))
    }

    pub(super) fn transformed(&self, transform: RigidTransform) -> Result<Self, GeometryError> {
        Ok(Self {
            feature: self.feature,
            stock: self.stock.transformed(transform)?,
            tool: self.tool.transformed(transform)?,
            wall: self.wall,
            bore: self.bore,
            local: self.local,
            inherited: Rc::clone(&self.inherited),
        })
    }
    pub(super) fn volume_mm3(&self) -> Result<f64, GeometryError> {
        let radius = self.bore().radius_mm;
        Ok(self.stock.volume_mm3()?
            - (self.tool.volume_mm3()? - PI * radius * radius * self.length()))
    }
    pub(super) fn bounds(&self) -> Result<Bounds, GeometryError> {
        self.stock.bounds()
    }

    /// Stock faces, the collar split, steps at tool caps 0 and 1, then the
    /// tool's helical bands in their original order.
    pub(super) fn faces(&self) -> Result<Vec<Face>, GeometryError> {
        let mut faces = self.stock.faces()?;
        let tool = self.tool.faces()?;
        let radius = self.bore().radius_mm;
        let (_, s, along, [low, high]) = self.axial();
        let depth = self.bore().depth_mm();
        let wall = self.wall_face();
        let collar = |a: f64, b: f64| {
            (
                self.world(self.local_point(a.midpoint(b))),
                TAU * radius * (b - a),
            )
        };
        (faces[wall].center_mm, faces[wall].area_mm2) = collar(0., low);
        let (center_mm, area_mm2) = collar(high, depth);
        faces.push(Face {
            index: faces.len() as u32,
            plane: None,
            center_mm,
            area_mm2,
        });
        let disk = PI * radius * radius;
        for (cap, face) in tool[..2].iter().enumerate() {
            let z = if cap == 0 {
                s
            } else {
                s + along * self.length()
            };
            let axis_point = self.world(self.local_point(z));
            let area = face.area_mm2 - disk;
            let center_mm = std::array::from_fn(|i| {
                (face.center_mm[i] * face.area_mm2 - axis_point[i] * disk) / area
            });
            faces.push(Face {
                index: faces.len() as u32,
                // Material lies outside the tool: each step faces into it.
                plane: face
                    .plane
                    .map(|(origin, normal)| (origin, normal.map(|x| -x))),
                center_mm,
                area_mm2: area,
            });
        }
        for face in &tool[2..] {
            faces.push(Face {
                index: faces.len() as u32,
                ..face.clone()
            });
        }
        Ok(faces)
    }

    fn origin(&self, role: &str) -> Origin {
        Origin {
            source: self.feature,
            role: format!("thread.bore/{role}"),
        }
    }
    fn circle(&self, s: f64) -> AnalyticGeometry {
        let (_, frame) = self.stock_parts();
        let (d, ..) = self.axial();
        let (_, [u, _]) = self.local;
        let to_world = |x: Point| -> Point {
            std::array::from_fn(|i| {
                x[0] * frame.u[i] + x[1] * frame.v[i] + x[2] * frame.normal()[i]
            })
        };
        let u = to_world(u);
        let n = to_world(d);
        AnalyticGeometry::CircularArc {
            frame: ReferenceFrame {
                origin_mm: self.world(self.local_point(s)),
                u,
                v: cross(n, u),
            },
            radius_mm: self.bore().radius_mm,
            start_rad: 0.,
            sweep_rad: TAU,
        }
    }

    /// Each redrilled wall, and every element it bounds, also keeps the
    /// original prism wall origin of the same bore.
    fn inherit_bore_walls(&self, stock: &mut TopologySnapshot) {
        let (body, _) = self.stock_parts();
        let (_, bores, origins, _) = body.round_cuts();
        let first = self.stock_face_count() - bores.len();
        for (i, original) in self.inherited.iter().enumerate() {
            let Some(wall) = stock
                .elements
                .iter()
                .find(|e| e.key.kind == TopologyKind::Face && e.key.ordinal as usize == first + i)
            else {
                continue;
            };
            let marker = wall.origins.clone();
            debug_assert!(marker.contains(&origins[i]));
            for element in &mut stock.elements {
                if !element.origins.is_disjoint(&marker) {
                    element.origins.insert(original.clone());
                }
                if !element.boundary.is_disjoint(&marker) {
                    element.boundary.insert(original.clone());
                }
            }
        }
    }

    pub(crate) fn topology(
        &self,
        owner: EntityId,
        revision: u64,
    ) -> Result<TopologySnapshot, GeometryError> {
        let mut stock = self.stock.topology_snapshot(owner, revision)?;
        self.inherit_bore_walls(&mut stock);
        let tool = self.tool.topology_snapshot(owner, revision)?;
        let bore = self.bore();
        let wall = self.wall_face() as u32;
        let face_count = self.stock_face_count() as u32;
        let count = |kind| stock.elements.iter().filter(|e| e.key.kind == kind).count() as u32;
        let (edge_base, vertex_base) = (count(TopologyKind::Edge), count(TopologyKind::Vertex));
        let (_, s, _, [low, _]) = self.axial();
        let [mouth_rim, exit_rim] = ["collar.mouth/rim", "collar.exit/rim"].map(|r| self.origin(r));
        let [mouth_collar, exit_collar] = ["collar.mouth", "collar.exit"].map(|r| self.origin(r));
        // Rims at the low and high collar ends; cap 0 lies at the tool start.
        let cap_rims = if s == low {
            [mouth_rim.clone(), exit_rim.clone()]
        } else {
            [exit_rim.clone(), mouth_rim.clone()]
        };
        let mouth_faces = bore.faces.map(|face| {
            stock
                .elements
                .iter()
                .find(|e| e.key.kind == TopologyKind::Face && e.key.ordinal == face as u32)
                .map(|e| e.origins.clone())
                .unwrap_or_default()
        });
        let mut elements = Vec::with_capacity(stock.elements.len() + tool.elements.len() + 4);
        let mut exit = None;
        for mut element in stock.elements {
            if element.key.kind == TopologyKind::Face && element.key.ordinal == wall {
                let mut other = element.clone();
                element.origins.insert(mouth_collar.clone());
                element.boundary = mouth_faces[0]
                    .iter()
                    .cloned()
                    .chain([mouth_rim.clone()])
                    .collect();
                other.key.ordinal = face_count;
                other.origins.insert(exit_collar.clone());
                other.boundary = mouth_faces[1]
                    .iter()
                    .cloned()
                    .chain([exit_rim.clone()])
                    .collect();
                exit = Some(other);
            }
            elements.push(element);
        }
        elements.push(exit.ok_or_else(|| {
            error(
                OPERATION,
                GeometryErrorKind::InvalidResult,
                "Retained stock topology lacks the receiving bore wall",
            )
        })?);
        for mut element in tool.elements {
            match element.key.kind {
                TopologyKind::Face if element.key.ordinal < 2 => {
                    let cap = element.key.ordinal as usize;
                    let own = self.origin(["step.start", "step.end"][cap]);
                    element.key.ordinal = face_count + 1 + cap as u32;
                    element.origins.insert(own.clone());
                    element.boundary.insert(cap_rims[cap].clone());
                    if let Some(AnalyticGeometry::Plane { frame }) = element.geometry {
                        let flipped = ReferenceFrame {
                            v: frame.v.map(|x| -x),
                            ..frame
                        };
                        element.geometry = Some(AnalyticGeometry::Plane { frame: flipped });
                        element.frames = vec![AnchoredFrame {
                            origin: own,
                            frame: flipped,
                        }];
                    }
                }
                TopologyKind::Face => element.key.ordinal += face_count + 1,
                TopologyKind::Edge => element.key.ordinal += edge_base + 2,
                TopologyKind::Vertex => element.key.ordinal += vertex_base,
            }
            elements.push(element);
        }
        let (_, _, _, [low, high]) = self.axial();
        for (ordinal, (rim, z)) in [(mouth_rim, low), (exit_rim, high)].into_iter().enumerate() {
            elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Edge,
                    ordinal: edge_base + ordinal as u32,
                },
                origins: BTreeSet::from([rim]),
                boundary: BTreeSet::new(),
                frames: vec![],
                geometry: Some(self.circle(z)),
            });
        }
        let snapshot = TopologySnapshot {
            owner,
            revision,
            elements,
        };
        snapshot
            .validate()
            .map_err(|cause| error(OPERATION, GeometryErrorKind::InvalidResult, cause))?;
        Ok(snapshot)
    }

    pub(super) fn mesh(
        &self,
        options: multicad_rust_kernel::display::Tessellation,
    ) -> Result<Mesh, GeometryError> {
        let (body, frame) = self.stock_parts();
        let (pockets, bores, ..) = body.round_cuts();
        let (origin, axes) = self.local;
        let local = threaded_bored_orthogonal(
            &body.boundary,
            pockets,
            bores,
            ThreadedThroughBore {
                bore: self.bore,
                wall: match &self.wall {
                    Wall::Thread(thread) => BandWall::Thread(thread),
                    &Wall::Cosmetic(radius_mm, length_mm) => BandWall::Cylinder {
                        radius_mm,
                        length_mm,
                    },
                },
                origin,
                axes,
            },
            options,
        );
        crate::rust_display::placed_mesh(local, frame, false)
    }
}

#[cfg(test)]
mod tests;
