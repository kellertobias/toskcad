//! Exact cylindrical stock minus one through tapered internal pipe thread.
//! The cavity is the unchanged certified radial boundary in the stock chart:
//! its caps lie on the stock's end planes, its axis is the stock axis.
use super::*;
use crate::geometry::internal_thread_treatment_body::interval::Interval as I;
use crate::rust_prism_topology::position;
use multicad_core::references::{
    AnalyticGeometry, Origin, TopologyElement, TopologyKey, TopologyKind,
};
use multicad_rust_kernel::{
    axial_pipe_thread::authored::{
        AuthoredSpec, AuthoredThread,
        hardware::mesh::{MeshOptions as JoinOptions, join_sleeve},
    },
    radial_tensor::{RadialBoundary, Surface},
    step_read::{ImportedPipeSleeve, SleeveStock},
};
use std::collections::BTreeSet;
use std::f64::consts::PI;
mod clearance;
#[cfg(test)]
mod tests;

const OPERATION: &str = "Rust tapered pipe sleeve";

fn failure(kind: GeometryErrorKind, message: impl Into<String>) -> GeometryError {
    error(OPERATION, kind, message)
}
fn pi() -> I {
    I::new(PI.next_down(), PI.next_up())
}

/// The authored law of the cavity and its complete Euclidean bound to the
/// actual retained cavity boundary, both in the stock chart.
#[derive(Clone, Copy, Debug)]
pub(crate) struct CavityLaw {
    pub(crate) transport: AuthoredThread,
    pub(crate) error_mm: f64,
    /// Sleeve face of every law role; the extra role is the stock wall.
    pub(crate) faces: super::law_print::LawFaces,
}

pub(crate) struct PipeSleeveBody {
    pub(super) feature: EntityId,
    /// Certified cavity material at the identity occurrence of the chart.
    cavity: Solid,
    radius_mm: f64,
    height_mm: f64,
    /// Cavity face indices of the bottom (z=0) and top (z=height) caps.
    caps: [usize; 2],
    /// Present only for an explicit authored factory, never for imports.
    law: Option<CavityLaw>,
}

impl Solid {
    /// Cylindrical stock `radius_mm` over the thread's own length minus one
    /// internal tapered pipe thread, mouth at `z=0` of `frame`.
    pub fn rust_tapered_pipe_sleeve_with_interrupt(
        spec: AuthoredSpec,
        radius_mm: f64,
        feature: EntityId,
        frame: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, GeometryError> {
        if spec.radial_taper >= 0. || spec.crest_radius_limit_mm.is_some() {
            return Err(failure(
                GeometryErrorKind::Unsupported,
                "A tapered pipe sleeve receives an internal thread whose taper narrows from its mouth, without a crest limit",
            ));
        }
        let cavity = Self::rust_authored_axial_pipe_thread_with_interrupt(
            spec,
            feature,
            PlaneFrame::default(),
            policy,
            interrupt,
        )?;
        let (source, occurrence) =
            cavity
                .rust_authored_axial_pipe_thread_source()
                .ok_or_else(|| {
                    failure(
                        GeometryErrorKind::Unsupported,
                        "The internal pipe thread has no retained authored cavity",
                    )
                })?;
        if occurrence != PlaneFrame::default() || source.construction_frame != PlaneFrame::default()
        {
            return Err(failure(
                GeometryErrorKind::Unsupported,
                "The retained cavity is not in the stock chart",
            ));
        }
        let (body, at) = cavity
            .radial_boundary_body()
            .expect("authored radial cavity");
        let boundary = body.placed_boundary(at)?;
        let faces = law_faces(source, boundary.input().faces.len())?;
        let law = CavityLaw {
            transport: source.transport,
            error_mm: (source.local_source_error_mm + boundary.boundary_error_mm()).next_up(),
            faces,
        };
        let caps = source.cap_faces;
        Self::pipe_sleeve(
            cavity,
            caps,
            radius_mm,
            spec.length_mm,
            Some(law),
            feature,
            frame,
            policy,
        )
    }

    /// Independent STEP readback: the cavity re-earns its radial material
    /// certificate in the file's own world chart before it is moved into
    /// the stock chart. No authored law is inferred from its poles.
    pub(super) fn import_pipe_sleeve(
        imported: ImportedPipeSleeve,
        feature: EntityId,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        let world = Self::qualify_imported_radial(imported.cavity, policy)?;
        let SleeveStock {
            origin,
            normal,
            reference,
            radius_mm,
            height_mm,
        } = imported.stock;
        let frame = PlaneFrame {
            origin,
            u: reference,
            v: cross(normal, reference),
        };
        frame.checked(policy)?;
        // Rows of the frame rotation map world offsets into the chart.
        let rows = [frame.u, frame.v, frame.normal()];
        let column = |k: usize| rows.map(|r| r[k]);
        let inverse_origin =
            std::array::from_fn(|i| -(0..3).map(|k| rows[i][k] * origin[k]).sum::<f64>());
        let local = world
            .transformed_with_interrupt(
                [inverse_origin, column(0), column(1)],
                policy.linear_tolerance_mm,
                policy.minimum_feature_mm,
                None,
            )
            .map_err(|cause| {
                failure(
                    GeometryErrorKind::Unsupported,
                    format!("Imported cavity does not transfer into the stock chart: {cause:?}"),
                )
            })?;
        let mut caps = [None; 2];
        for (index, face) in local.input().faces.iter().enumerate() {
            if let Surface::Plane(axis) = local.input().surfaces[face.surface] {
                let end = usize::from(axis.origin[2] > height_mm / 2.);
                if caps[end].replace(index).is_some() {
                    return Err(failure(
                        GeometryErrorKind::Unsupported,
                        "The imported cavity has more than one cap per stock end",
                    ));
                }
            }
        }
        let [Some(bottom), Some(top)] = caps else {
            return Err(failure(
                GeometryErrorKind::Unsupported,
                "The imported cavity does not reach both stock ends",
            ));
        };
        let cavity = Self::from_rust_radial_boundary(local, feature, policy)?;
        Self::pipe_sleeve(
            cavity,
            [bottom, top],
            radius_mm,
            height_mm,
            None,
            feature,
            frame,
            policy,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn pipe_sleeve(
        cavity: Solid,
        caps: [usize; 2],
        radius_mm: f64,
        height_mm: f64,
        law: Option<CavityLaw>,
        feature: EntityId,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        policy.length(radius_mm)?;
        policy.length(height_mm)?;
        frame.checked(policy)?;
        frame.checked_metric_extent([radius_mm, radius_mm, height_mm], policy)?;
        let (body, at) = cavity.radial_boundary_body().ok_or_else(|| {
            failure(
                GeometryErrorKind::Unsupported,
                "A pipe sleeve cavity must be retained radial material",
            )
        })?;
        let boundary = body.placed_boundary(at)?;
        let tolerance = policy.linear_tolerance_mm;
        // The cavity caps are the stock's end planes and its axis the stock
        // axis; otherwise the stock caps would not be annuli of this cavity.
        let axes = boundary.cap_axes();
        let input = boundary.input();
        for (end, (&face, z)) in caps.iter().zip([0., height_mm]).enumerate() {
            let Surface::Plane(axis) = input.surfaces[input.faces[face].surface] else {
                return Err(failure(
                    GeometryErrorKind::Unsupported,
                    "A pipe sleeve cap is not a retained plane",
                ));
            };
            let matches = |a: multicad_rust_kernel::radial_tensor::Axis| {
                a.origin[0].portable_hypot(a.origin[1]) <= tolerance
                    && (a.origin[2] - z).abs() <= tolerance
                    && a.normal[0].portable_hypot(a.normal[1]) <= 1e-12
            };
            if !matches(axis) || !matches(axes[end]) {
                return Err(failure(
                    GeometryErrorKind::Unsupported,
                    "The cavity does not run through the stock between its end planes",
                ));
            }
        }
        // The complete retained cavity is strictly inside the stock wall:
        // the authored radius plus the actual-boundary bound, or for an
        // import the support of the actual boundary in 32 directions.
        let reach = if let Some(law) = law {
            let spec = law.transport.spec();
            let outer = law
                .transport
                .radial_bounds([0., spec.pitch_mm], [0., spec.length_mm])
                .map_err(|cause| {
                    failure(
                        GeometryErrorKind::Unsupported,
                        format!("Cavity radius is not resolved: {cause:?}"),
                    )
                })?[1];
            (outer + law.error_mm).next_up()
        } else {
            let mut support = 0f64;
            for k in 0..32 {
                let angle = std::f64::consts::TAU * k as f64 / 32.;
                let [_, hi] = boundary
                    .support_bounds(
                        [angle.portable_cos(), angle.portable_sin(), 0.],
                        tolerance,
                        None,
                    )
                    .map_err(|cause| {
                        failure(
                            GeometryErrorKind::Unsupported,
                            format!("Cavity support is not resolved: {cause:?}"),
                        )
                    })?;
                support = support.max(hi);
            }
            // A circumscribed 32-gon of the supports: r <= h / cos(pi/32).
            (support / (PI / 32.).portable_cos() * (1. + 1e-12)).next_up()
        };
        let wall = radius_mm - reach;
        if wall.is_nan() || wall <= policy.minimum_feature_mm {
            return Err(failure(
                GeometryErrorKind::Input,
                "The receiving stock does not clear the complete internal thread",
            ));
        }
        let body = PipeSleeveBody {
            feature,
            cavity,
            radius_mm,
            height_mm,
            caps,
            law,
        };
        let bounds = body.bounds(frame);
        policy.point(bounds.min)?;
        policy.point(bounds.max)?;
        body.volume_interval()?;
        Ok(Self {
            body: SolidBody::PipeSleeve(Rc::new(body), frame),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        })
    }

    pub(crate) fn pipe_sleeve_body(&self) -> Option<(&PipeSleeveBody, PlaneFrame)> {
        match &self.body {
            SolidBody::PipeSleeve(body, frame) => Some((body, *frame)),
            _ => None,
        }
    }
}

/// Face of every law role of a retained authored pipe source.
pub(super) fn law_faces(
    source: &super::radial_boundary_body::RetainedAxialPipeThreadSource,
    extra: usize,
) -> Result<super::law_print::LawFaces, GeometryError> {
    let bad = || {
        failure(
            GeometryErrorKind::InvalidResult,
            "The retained pipe source has no face for a law role",
        )
    };
    let band_count = source.band_faces.len();
    if band_count > 5 || band_count != source.transport.band_count() {
        return Err(bad());
    }
    let mut bands = [u32::MAX; 5];
    for (b, &f) in source.band_faces.iter().enumerate() {
        bands[b] = f as u32;
    }
    let root = *source.root_faces.first().ok_or_else(bad)? as u32;
    Ok(super::law_print::LawFaces {
        band_count: band_count as u32,
        bands,
        root,
        // No crest limit: a stock cylinder role never occurs.
        clip: u32::MAX,
        caps: source.cap_faces.map(|f| f as u32),
        extra: extra as u32,
    })
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn rotate(frame: PlaneFrame, d: [f64; 3]) -> [f64; 3] {
    let n = frame.normal();
    std::array::from_fn(|i| frame.u[i] * d[0] + frame.v[i] * d[1] + n[i] * d[2])
}

impl PipeSleeveBody {
    fn boundary(&self) -> Result<Rc<RadialBoundary>, GeometryError> {
        let (body, at) = self.cavity.radial_boundary_body().expect("sealed cavity");
        body.placed_boundary(at)
    }
    pub(crate) fn law(&self) -> Option<CavityLaw> {
        self.law
    }
    pub(crate) fn radius_mm(&self) -> f64 {
        self.radius_mm
    }
    pub(crate) fn height_mm(&self) -> f64 {
        self.height_mm
    }
    fn volume_interval(&self) -> Result<I, GeometryError> {
        let (body, at) = self.cavity.radial_boundary_body().expect("sealed cavity");
        let [lo, hi] = body.volume_properties(at, None)?.volume_mm3;
        let stock = pi() * I::point(self.radius_mm).square() * I::point(self.height_mm);
        let v = stock - I::new(lo, hi);
        if !v.valid() || v.lo <= 0. {
            return Err(failure(
                GeometryErrorKind::Unsupported,
                "The pipe sleeve volume is not resolved",
            ));
        }
        Ok(v)
    }
    pub(super) fn volume_mm3(&self) -> Result<f64, GeometryError> {
        let v = self.volume_interval()?;
        Ok(v.lo + (v.hi - v.lo) * 0.5)
    }
    pub(super) fn bounds(&self, frame: PlaneFrame) -> Bounds {
        let r = self.radius_mm;
        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];
        for corner in 0..8 {
            let local = [
                if corner & 1 == 0 { -r } else { r },
                if corner & 2 == 0 { -r } else { r },
                if corner & 4 == 0 { 0. } else { self.height_mm },
            ];
            let p = position(frame, local);
            for k in 0..3 {
                min[k] = min[k].min(p[k]);
                max[k] = max[k].max(p[k]);
            }
        }
        // Rounded placement of a box corner moves by at most a few ulps.
        let slack = 64. * f64::EPSILON * (r + self.height_mm + 1.);
        Bounds {
            min: min.map(|v| v - slack - v.abs() * f64::EPSILON),
            max: max.map(|v| v + slack + v.abs() * f64::EPSILON),
        }
    }
    /// Cavity faces keep their indices (caps become annuli), the cylinder
    /// wall is the last face.
    pub(super) fn faces(&self, frame: PlaneFrame) -> Result<Vec<Face>, GeometryError> {
        let (body, at) = self.cavity.radial_boundary_body().expect("sealed cavity");
        let properties = body.face_properties(at, None)?;
        let disc = pi() * I::point(self.radius_mm).square();
        let n = frame.normal();
        let mut faces = Vec::with_capacity(properties.len() + 1);
        for face in properties {
            let area = I::new(face.area_mm2[0], face.area_mm2[1]);
            let centre = face.centroid_mm.map(|[a, b]| I::new(a, b));
            let (area, centre, plane) =
                if let Some(end) = self.caps.iter().position(|&c| c == face.face_index) {
                    let z = [0., self.height_mm][end];
                    let annulus = disc - area;
                    let c: [I; 3] = std::array::from_fn(|k| {
                        let full = if k == 2 { disc * I::point(z) } else { I::ZERO };
                        (full - area * centre[k]) / annulus
                    });
                    let sign = if end == 0 { -1. } else { 1. };
                    (
                        annulus,
                        c,
                        Some((position(frame, [0., 0., z]), n.map(|v| v * sign))),
                    )
                } else {
                    (area, centre, None)
                };
            if !area.valid() || area.lo <= 0. || centre.iter().any(|v| !v.valid()) {
                return Err(failure(
                    GeometryErrorKind::Unsupported,
                    "A pipe sleeve face measure is not resolved",
                ));
            }
            faces.push(Face {
                index: face.face_index as u32,
                plane,
                center_mm: position(frame, centre.map(|v| v.lo + (v.hi - v.lo) * 0.5)),
                area_mm2: area.lo + (area.hi - area.lo) * 0.5,
            });
        }
        let wall = I::point(2.) * pi() * I::point(self.radius_mm) * I::point(self.height_mm);
        faces.push(Face {
            index: faces.len() as u32,
            plane: None,
            center_mm: position(frame, [0., 0., self.height_mm / 2.]),
            area_mm2: wall.lo + (wall.hi - wall.lo) * 0.5,
        });
        Ok(faces)
    }
    pub(super) fn mesh(
        &self,
        frame: PlaneFrame,
        options: multicad_rust_kernel::display::Tessellation,
    ) -> Result<Mesh, GeometryError> {
        if let Some(law) = self.law {
            return self.law_mesh(law, frame, options);
        }
        let boundary = self.boundary()?;
        let mesh_failure = |what: &str, cause: String| {
            failure(
                if cause.contains("Limit") {
                    GeometryErrorKind::Limit
                } else if cause.contains("Cancelled") {
                    GeometryErrorKind::Cancelled
                } else {
                    GeometryErrorKind::Unsupported
                },
                format!("{what}: {cause}"),
            )
        };
        let trace = boundary
            .tessellate_with_cap_trace(options, Default::default(), None)
            .map_err(|e| mesh_failure("Cavity trace", format!("{e:?}")))?;
        let joined = join_sleeve(
            &trace,
            self.radius_mm,
            JoinOptions {
                tessellation: options,
                maximum_work: 512_000_000,
                maximum_bytes: 64 * 1024 * 1024,
                additional_resident_bytes: boundary.retained_budget_bytes(),
            },
            None,
        )
        .map_err(|e| mesh_failure("Stock join", format!("{e:?}")))?;
        crate::rust_display::placed_mesh(Ok(joined.into_mesh()), frame, false)
    }
    /// Dedicated authored-law bore joined to the stock cylinder. Its facet
    /// bound plus the cavity's own law-to-boundary bound meets the request.
    fn law_mesh(
        &self,
        law: CavityLaw,
        frame: PlaneFrame,
        options: multicad_rust_kernel::display::Tessellation,
    ) -> Result<Mesh, GeometryError> {
        use multicad_rust_kernel::axial_pipe_thread::authored::{
            LawClosure, LawRole, hardware::mesh::join_sleeve_law,
        };
        let bore = law
            .transport
            .tessellate_law(options, LawClosure::Bore, None)
            .map_err(super::law_print::failure)?;
        let roles = [
            bore.role_index(LawRole::StartCap),
            bore.role_index(LawRole::EndCap),
            law.faces.extra_role(),
        ];
        let (mesh, _, distance) = join_sleeve_law(
            bore,
            [0., self.height_mm],
            roles,
            self.radius_mm,
            JoinOptions {
                tessellation: options,
                maximum_work: 512_000_000,
                maximum_bytes: 64 * 1024 * 1024,
                additional_resident_bytes: 0,
            },
            None,
        )
        .map_err(super::law_print::failure)?;
        super::law_print::admit(distance, law.error_mm, options.deflection)?;
        let mut mesh = mesh;
        law.faces.remap(&mut mesh.face_indices)?;
        crate::rust_display::placed_mesh(Ok(mesh), frame, false)
    }
    fn stock(&self, frame: PlaneFrame) -> SleeveStock {
        SleeveStock {
            origin: frame.origin,
            normal: frame.normal(),
            reference: frame.u,
            radius_mm: self.radius_mm,
            height_mm: self.height_mm,
        }
    }
    fn world_cavity(&self, frame: PlaneFrame) -> Result<Rc<RadialBoundary>, GeometryError> {
        let (body, _) = self.cavity.radial_boundary_body().expect("sealed cavity");
        body.placed_boundary(frame)
    }
    pub(super) fn supports_step(&self, frame: PlaneFrame, tolerance: f64) -> bool {
        self.world_cavity(frame).is_ok_and(|cavity| {
            multicad_rust_kernel::step::supports_pipe_sleeve(&cavity, self.stock(frame), tolerance)
        })
    }
    pub(super) fn append_step(
        &self,
        writer: &mut multicad_rust_kernel::step::Writer,
        frame: PlaneFrame,
        tolerance: f64,
    ) -> Result<(), multicad_rust_kernel::step::StepError> {
        let cavity = self.world_cavity(frame).map_err(|cause| {
            if cause.kind == GeometryErrorKind::Limit {
                multicad_rust_kernel::step::StepError::Limit
            } else {
                multicad_rust_kernel::step::StepError::Unsupported
            }
        })?;
        writer.append_pipe_sleeve(&cavity, self.stock(frame), tolerance)
    }
    /// Cavity topology at this occurrence with inward thread roles, plus the
    /// stock wall, its two rim circles and its seam-free boundary sets.
    pub(crate) fn topology(
        &self,
        frame: PlaneFrame,
        owner: EntityId,
        revision: u64,
    ) -> Result<TopologySnapshot, GeometryError> {
        let (body, _) = self.cavity.radial_boundary_body().expect("sealed cavity");
        let mut result = body.topology(frame, owner, revision)?;
        let inward = |origin: &Origin| Origin {
            source: origin.source,
            role: if origin.role.starts_with("thread/") || origin.role == "thread.root" {
                format!("receiving.inward/{}", origin.role)
            } else {
                origin.role.clone()
            },
        };
        for element in &mut result.elements {
            element.origins = element.origins.iter().map(inward).collect();
            element.boundary = element.boundary.iter().map(inward).collect();
        }
        let origin = |role: &str| Origin {
            source: self.feature,
            role: role.into(),
        };
        let wall = origin("stock.wall");
        let rims = [origin("stock.rim/start"), origin("stock.rim/end")];
        let count = |kind: TopologyKind| {
            result
                .elements
                .iter()
                .filter(|e| e.key.kind == kind)
                .map(|e| e.key.ordinal + 1)
                .max()
                .unwrap_or(0)
        };
        let first_face = count(TopologyKind::Face);
        let first_edge = count(TopologyKind::Edge);
        let mut caps: [Option<Origin>; 2] = [const { None }; 2];
        for element in &mut result.elements {
            if element.key.kind == TopologyKind::Face
                && let Some(end) = self
                    .caps
                    .iter()
                    .position(|&c| c as u32 == element.key.ordinal)
            {
                element.boundary.insert(rims[end].clone());
                caps[end] = element.origins.iter().next().cloned();
            }
        }
        let [Some(start), Some(end)] = caps else {
            return Err(failure(
                GeometryErrorKind::InvalidResult,
                "The pipe sleeve caps have no cavity topology",
            ));
        };
        let mut added = vec![TopologyElement {
            key: TopologyKey {
                kind: TopologyKind::Face,
                ordinal: first_face,
            },
            origins: BTreeSet::from([wall.clone()]),
            boundary: BTreeSet::from(rims.clone()),
            frames: vec![],
            geometry: None,
        }];
        for (level, cap) in [start, end].into_iter().enumerate() {
            let centre = position(frame, [0., 0., [0., self.height_mm][level]]);
            added.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Edge,
                    ordinal: first_edge + level as u32,
                },
                origins: BTreeSet::from([rims[level].clone()]),
                boundary: BTreeSet::from([cap, wall.clone()]),
                frames: vec![],
                geometry: Some(AnalyticGeometry::CircularArc {
                    frame: ReferenceFrame {
                        origin_mm: centre,
                        u: frame.u,
                        v: rotate(frame, [0., 1., 0.]),
                    },
                    radius_mm: self.radius_mm,
                    start_rad: 0.,
                    sweep_rad: std::f64::consts::TAU,
                }),
            });
        }
        result.elements.extend(added);
        result
            .validate()
            .map_err(|e| failure(GeometryErrorKind::InvalidResult, e))?;
        Ok(result)
    }
    pub(super) fn display_edges(
        &self,
        frame: PlaneFrame,
        deflection: f64,
    ) -> Result<Vec<DisplayEdge>, GeometryError> {
        let mut topology = self.topology(frame, self.feature, 0)?;
        topology.elements.retain(|e| {
            e.key.kind == TopologyKind::Edge
                && matches!(e.geometry, Some(AnalyticGeometry::CircularArc { .. }))
        });
        crate::rust_display::edges(None, &topology, deflection)
    }
}
