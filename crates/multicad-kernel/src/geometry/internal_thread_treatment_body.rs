//! Actual cylindrical stock minus an explicitly authored shared-boundary MAX
//! cavity. No display mesh, legacy thread approximation or recipe is material.
pub(super) mod interval;
mod properties;
use super::arithmetic::{dot_error, sum_error};
use super::*;
use crate::HelicalThread;
use multicad_core::references::TopologySnapshot;
use multicad_rust_kernel::{
    authored_flat_thread::{AuthoredThread, SourceOptions},
    authored_internal_treatment::{self as source, Cavity, Enclosure, Material, QueryOptions},
    display::Tessellation,
};
pub use properties::InternalThreadFaceProperties;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InternalThreadTreatmentCut {
    pub profile: HelicalThread,
    pub intent: source::Intent,
    pub minimum_wall_mm: f64,
    pub frame: PlaneFrame,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InternalThreadTreatmentData {
    /// Source and candidates use the stock-relative chart. The original
    /// cylinder origin remains exact data in material.stock_chart_origin().
    pub material: Material,
    /// Active qualified placement of that chart, with independent axial-WCS
    /// or explicitly framed stock authority; never inferred from mesh points.
    pub frame: PlaneFrame,
}
pub(crate) struct InternalThreadTreatmentBody {
    material: Material,
    pub(super) feature: EntityId,
    requested_wall: f64,
    stock_topology: TopologySnapshot,
}
fn unsupported(message: &str) -> GeometryError {
    error(
        "Rust internal thread treatment",
        GeometryErrorKind::Unsupported,
        message,
    )
}
fn failure(cause: source::Error) -> GeometryError {
    error(
        "Rust internal thread treatment",
        match cause {
            source::Error::Invalid => GeometryErrorKind::Input,
            source::Error::Limit => GeometryErrorKind::Limit,
            source::Error::Unresolved => GeometryErrorKind::Unsupported,
            source::Error::Cancelled => GeometryErrorKind::Cancelled,
        },
        format!("Actual internal treatment source proof failed: {cause:?}"),
    )
}
fn checkpoint(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), GeometryError> {
    if interrupt.is_some_and(|f| f()) {
        Err(failure(source::Error::Cancelled))
    } else {
        Ok(())
    }
}
fn up(a: f64, b: f64) -> f64 {
    (a + b).next_up()
}
fn product(a: f64, b: f64) -> f64 {
    (a * b).next_up()
}
// Literal stored columns, including the rounded cross product used by the
// public frame, define the actual affine placement. Never normalize them.
fn metric(frame: PlaneFrame) -> Result<([f64; 2], f64), GeometryError> {
    let axes = [frame.u, frame.v, frame.normal()];
    let mut gram = 0.;
    for (i, a) in axes.iter().enumerate() {
        for (j, b) in axes.iter().enumerate() {
            let (value, error) = dot_error(*a, *b);
            let (value, residual) = sum_error(value, if i == j { -1. } else { 0. });
            gram = up(gram, up(value.abs(), up(error, residual)));
        }
    }
    let pairs = [(1, 2), (2, 0), (0, 1)];
    let (mut normal, mut errors) = ([0.; 3], [0.; 3]);
    for (i, (j, k)) in pairs.into_iter().enumerate() {
        (normal[i], errors[i]) =
            dot_error([axes[0][j], -axes[0][k], 0.], [axes[1][k], axes[1][j], 0.]);
    }
    let (det, mut error) = dot_error(normal, axes[2]);
    for i in 0..3 {
        error = up(error, product(errors[i], axes[2][i].abs()));
    }
    let determinant = [(det - error).next_down(), (det + error).next_up()];
    if !gram.is_finite() || gram > 1e-6 || determinant[0] <= 0. || !determinant[1].is_finite() {
        return Err(unsupported("Actual stock chart metric cannot be certified"));
    }
    Ok((determinant, (1. - gram).next_down().sqrt().next_down()))
}
fn radius_bounds(material: &Material) -> [f64; 2] {
    let (n, d) = material.stock().circle().radius().ratio();
    let integer = |v: i64| {
        let value = v as f64;
        if v.unsigned_abs() <= 1 << 53 {
            [value; 2]
        } else {
            [value.next_down(), value.next_up()]
        }
    };
    let [n, d] = [integer(n), integer(d)];
    [(n[0] / d[1]).next_down(), (n[1] / d[0]).next_up()]
}
fn radius_upper(material: &Material) -> f64 {
    radius_bounds(material)[1]
}
// Reduced bounded Rational has no subnormal/overflow range. A dyadic ratio
// is exactly binary64 iff its odd numerator fits the53-bit significand.
fn exact_binary(value: Rational, actual: f64) -> bool {
    let (n, d) = value.ratio();
    let magnitude = n.unsigned_abs();
    (d as u64).is_power_of_two()
        && (magnitude == 0 || 64 - magnitude.leading_zeros() - magnitude.trailing_zeros() <= 53)
        && value.as_f64() == actual
}

impl Solid {
    /// Restricted full-height receiving cut. An arbitrary retained legacy
    /// Thread is not admitted as authored shared-boundary material.
    pub fn rust_cut_authored_internal_thread_treatments_with_interrupt(
        &self,
        cut: InternalThreadTreatmentCut,
        feature: EntityId,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, GeometryError> {
        checkpoint(interrupt)?;
        policy.validate()?;
        if self.policy != policy || policy.boolean_fuzzy_mm != 0. {
            return Err(unsupported(
                "Internal treatment stock requires its original policy without Boolean fuzz",
            ));
        }
        let Some((stock, frame)) = self.axial_cylinder_frame() else {
            return Err(unsupported(
                "Internal treatment requires one actual unmachined cylinder",
            ));
        };
        if feature.is_nil() || feature.uuid_bits() == stock.feature() {
            return Err(input(
                "Internal treatment requires a distinct non-nil feature",
            ));
        }
        // No rounded center/clock/reverse reconciliation. Existing stored
        // stock axes and origin are the actual common source chart.
        if frame != cut.frame {
            return Err(unsupported(
                "Internal treatment requires the exact stored cylinder frame",
            ));
        }
        if self.axial_cylinder().is_some() {
            let [x, y] = stock.circle().center();
            if ![x, y, stock.bottom()]
                .into_iter()
                .zip(frame.origin)
                .all(|(value, actual)| exact_binary(value, actual))
                || frame.u != [1., 0., 0.]
                || frame.v != [0., 1., 0.]
            {
                return Err(unsupported(
                    "Axial stock requires exactly binary64 source coordinates and a cardinal WCS chart",
                ));
            }
        }
        policy.length(cut.minimum_wall_mm)?;
        if let Some(depth) = cut.intent.entry_depth_mm {
            policy.length(depth)?;
        }
        if let Some(tail) = cut.intent.terminal {
            policy.length(tail.width_mm)?;
            policy.length(tail.radial_depth_mm)?;
        }
        let legacy = helical_body::checked_thread(cut.profile, None, policy)?;
        let thread = AuthoredThread::new_with_interrupt(
            legacy.spec(),
            SourceOptions {
                position_limit_mm: policy.linear_tolerance_mm / 16.,
                normal_angle_limit_rad: policy.linear_tolerance_mm
                    / (16. * cut.profile.major_radius_mm.max(1.)),
            },
            interrupt,
        )
        // The identical Spec already passed checked_thread. Root isNone,
        // and both caller-derived SourceOptions are positive finite under the
        // admitted policy/radius bounds. Authored construction can now only
        // cancel or conservatively fail its independent transfer proof.
        .map_err(|cause| {
            if cause == multicad_rust_kernel::authored_flat_thread::Error::Cancelled {
                failure(source::Error::Cancelled)
            } else {
                unsupported("Authored flat-thread transfer cannot be certified")
            }
        })?;
        let cavity = Cavity::new(thread, cut.intent).map_err(failure)?;
        let (_, scale) = metric(frame)?;
        let local_wall = (cut.minimum_wall_mm / scale).next_up();
        let material =
            Material::new_with_interrupt(*stock, cavity, local_wall, interrupt).map_err(failure)?;
        let stock_topology =
            self.topology_snapshot(EntityId::from_uuid_bits(stock.feature()), 0)?;
        let body = InternalThreadTreatmentBody {
            material,
            feature,
            requested_wall: cut.minimum_wall_mm,
            stock_topology,
        };
        body.validate_frame(frame, policy)?;
        let bounds = body.bounds(frame);
        policy.point(bounds.min)?;
        policy.point(bounds.max)?;
        checkpoint(interrupt)?;
        Ok(Self {
            body: SolidBody::InternalThreadTreatment(Rc::new(body), frame),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        })
    }
    pub(crate) fn internal_thread_treatment_body(
        &self,
    ) -> Option<(&InternalThreadTreatmentBody, PlaneFrame)> {
        match &self.body {
            SolidBody::InternalThreadTreatment(body, frame) => Some((body, *frame)),
            _ => None,
        }
    }
    pub fn rust_internal_thread_treatment_data(&self) -> Option<InternalThreadTreatmentData> {
        self.internal_thread_treatment_body()
            .map(|(body, frame)| InternalThreadTreatmentData {
                material: body.material,
                frame,
            })
    }
    pub fn rust_internal_thread_treatment_volume_bounds(
        &self,
    ) -> Result<Option<Enclosure>, GeometryError> {
        self.internal_thread_treatment_body()
            .map(|(body, frame)| body.volume(frame))
            .transpose()
    }
    /// Immutable ancestry receipt for the original stock, not a topology
    /// snapshot of the changed finite cap wires or the new receiving walls.
    pub fn rust_internal_thread_treatment_original_stock_topology(
        &self,
    ) -> Option<&TopologySnapshot> {
        self.internal_thread_treatment_body()
            .map(|(body, _)| &body.stock_topology)
    }
}
impl InternalThreadTreatmentBody {
    pub(super) fn validate_frame(
        &self,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<(), GeometryError> {
        frame.checked_metric_extent(
            [
                radius_upper(&self.material),
                radius_upper(&self.material),
                self.material.cavity().full_length_mm(),
            ],
            policy,
        )?;
        let (_, scale) = metric(frame)?;
        let gap = (radius_bounds(&self.material)[0] - self.material.cavity().maximum_radius_mm())
            .next_down();
        if (gap * scale).next_down() < self.requested_wall {
            return Err(unsupported(
                "Placed stock no longer proves the requested material wall",
            ));
        }
        Ok(())
    }
    pub(super) fn bounds(&self, frame: PlaneFrame) -> Bounds {
        let radius = radius_upper(&self.material);
        let length = self.material.cavity().full_length_mm();
        let n = frame.normal();
        let reach: [f64; 3] =
            std::array::from_fn(|i| product(radius, up(frame.u[i].abs(), frame.v[i].abs())));
        Bounds {
            min: std::array::from_fn(|i| {
                ((frame.origin[i] + (length * n[i]).next_down().min(0.)).next_down() - reach[i])
                    .next_down()
            }),
            max: std::array::from_fn(|i| {
                up(
                    up(frame.origin[i], (length * n[i]).next_up().max(0.)),
                    reach[i],
                )
            }),
        }
    }
    pub(super) fn volume(&self, frame: PlaneFrame) -> Result<Enclosure, GeometryError> {
        let (det, _) = metric(frame)?;
        let local = self
            .material
            .volume_bounds(QueryOptions::default())
            .map_err(failure)?;
        let value = Enclosure {
            lower: (local.lower * det[0]).next_down(),
            upper: (local.upper * det[1]).next_up(),
        };
        if value.upper - value.lower > QueryOptions::default().tolerance_mm3 {
            return Err(unsupported(
                "Placed material volume reaches the numerical precision floor",
            ));
        }
        Ok(value)
    }
    pub(super) fn mesh(
        &self,
        frame: PlaneFrame,
        options: Tessellation,
    ) -> Result<Mesh, GeometryError> {
        entry_chamfer_body::placement::qualified_source_mesh(frame, options, |quality| {
            let candidate = self
                .material
                .tessellate_candidate(quality, None)
                .map_err(failure)?;
            Ok(Mesh {
                vertices: candidate.mesh.vertices,
                indices: candidate.mesh.indices,
                face_indices: candidate.mesh.face_indices,
            })
        })
    }
    pub(crate) fn unavailable_boundary(&self) -> GeometryError {
        unsupported(
            "Exact finite internal-treatment boundary properties are not yet qualified; mesh is not topology",
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use multicad_rust_kernel::{Circle, Rational};
    fn fixture_stock(height: Rational, frame: PlaneFrame) -> Solid {
        let zero = Rational::new(0, 1).unwrap();
        let cylinder = AxialCylinder::new(
            7,
            Circle::new([zero; 2], Rational::new(5, 2).unwrap()).unwrap(),
            zero,
            height,
        )
        .unwrap();
        Solid::build_axial_cylinder(cylinder, frame, None, GeometryPolicy::default()).unwrap()
    }
    fn cut(frame: PlaneFrame, left_hand: bool, starts: u32) -> InternalThreadTreatmentCut {
        InternalThreadTreatmentCut {
            profile: HelicalThread {
                major_radius_mm: 1.5,
                pitch_mm: 0.5,
                length_mm: 4.,
                depth_mm: 0.25,
                crest_width_mm: 0.125,
                lower_flank_deg: 30.,
                upper_flank_deg: 30.,
                radial_taper: 0.,
                starts,
                left_hand,
            },
            intent: source::Intent {
                entry_depth_mm: Some(0.5),
                terminal: Some(source::Terminal {
                    width_mm: 0.5,
                    radial_depth_mm: 0.125,
                }),
            },
            minimum_wall_mm: 0.5,
            frame,
        }
    }
    fn options() -> MeshOptions {
        MeshOptions {
            linear_deflection_mm: 0.2,
            angular_deflection_rad: 0.4,
            maximum_triangles: 1_000_000,
        }
    }
    #[test]
    fn actual_stock_source_mesh_and_placed_interval_are_owned_without_fabricated_properties() {
        let policy = GeometryPolicy::default();
        for (left, starts, posed) in [
            (false, 1, false),
            (true, 2, false),
            (false, 1, true),
            (true, 2, true),
        ] {
            let frame = if posed {
                PlaneFrame {
                    origin: [10., 20., 30.],
                    u: [0.8, 0.6, 0.],
                    v: [-0.36, 0.48, 0.8],
                }
            } else {
                PlaneFrame::default()
            };
            let stock = fixture_stock(Rational::new(4, 1).unwrap(), frame);
            let original = stock.volume_mm3().unwrap();
            let request = cut(frame, left, starts);
            let result = stock
                .rust_cut_authored_internal_thread_treatments_with_interrupt(
                    request,
                    EntityId::from_uuid_bits(8),
                    policy,
                    None,
                )
                .unwrap();
            let actual = result.rust_internal_thread_treatment_data().unwrap();
            assert_eq!(actual.frame, frame);
            assert_eq!(
                actual.material.stock(),
                stock.axial_cylinder_frame().unwrap().0
            );
            assert_eq!(
                actual.material.cavity().source_thread().spec().left_hand,
                left
            );
            assert_eq!(
                actual.material.cavity().source_thread().spec().starts,
                starts
            );
            assert_eq!(actual.material.cavity().intent(), request.intent);
            let interval = result
                .rust_internal_thread_treatment_volume_bounds()
                .unwrap()
                .unwrap();
            let scalar = result.volume_mm3().unwrap();
            assert!(interval.lower <= scalar && scalar <= interval.upper && scalar < original);
            let receipt = result.manufacturing_mesh_with_receipt(options()).unwrap();
            assert!(receipt.edge_pairs > 0 && receipt.vertex_links > 0);
            assert_eq!(
                result.mesh(options()).unwrap().indices,
                receipt.mesh.indices
            );
            assert_eq!(result.solid_count().unwrap(), 1);
            assert!(
                result
                    .rust_internal_thread_treatment_original_stock_topology()
                    .unwrap()
                    .elements
                    .iter()
                    .any(|e| e
                        .origins
                        .iter()
                        .any(|o| o.source == EntityId::from_uuid_bits(7)))
            );
            assert_eq!(
                result.faces().unwrap_err().kind,
                GeometryErrorKind::Unsupported
            );
            assert_eq!(
                result
                    .topology_snapshot(EntityId::new(), 0)
                    .unwrap_err()
                    .kind,
                GeometryErrorKind::Unsupported
            );
            assert!(!result.supports_step_export());
            assert_eq!(result.share().unwrap().volume_mm3().unwrap(), scalar);
            assert_eq!(stock.volume_mm3().unwrap(), original);
        }
    }
    #[test]
    fn exact_height_pose_wall_policy_and_cancel_fail_atomically() {
        let policy = GeometryPolicy::default();
        let frame = PlaneFrame::default();
        let stock = fixture_stock(Rational::new(4, 1).unwrap(), frame);
        let original = stock.volume_mm3().unwrap();
        for bad in [
            InternalThreadTreatmentCut {
                minimum_wall_mm: 1.,
                ..cut(frame, false, 1)
            },
            InternalThreadTreatmentCut {
                frame: PlaneFrame {
                    origin: [0., 0., 1e-9],
                    ..frame
                },
                ..cut(frame, false, 1)
            },
            InternalThreadTreatmentCut {
                profile: HelicalThread {
                    length_mm: 3.5,
                    ..cut(frame, false, 1).profile
                },
                ..cut(frame, false, 1)
            },
        ] {
            assert!(
                stock
                    .rust_cut_authored_internal_thread_treatments_with_interrupt(
                        bad,
                        EntityId::new(),
                        policy,
                        None
                    )
                    .is_err()
            );
        }
        assert_eq!(
            stock
                .rust_cut_authored_internal_thread_treatments_with_interrupt(
                    cut(frame, false, 1),
                    EntityId::new(),
                    policy,
                    Some(&|| true)
                )
                .err()
                .unwrap()
                .kind,
            GeometryErrorKind::Cancelled
        );
        let nonbinary = fixture_stock(Rational::new(33, 10).unwrap(), frame);
        let request = InternalThreadTreatmentCut {
            profile: HelicalThread {
                length_mm: 3.3,
                ..cut(frame, false, 1).profile
            },
            ..cut(frame, false, 1)
        };
        assert_eq!(
            nonbinary
                .rust_cut_authored_internal_thread_treatments_with_interrupt(
                    request,
                    EntityId::new(),
                    policy,
                    None
                )
                .err()
                .unwrap()
                .kind,
            GeometryErrorKind::Unsupported
        );
        assert!(stock.rust_internal_thread_treatment_data().is_none());
        assert!(
            stock
                .rust_internal_thread_treatment_volume_bounds()
                .unwrap()
                .is_none()
        );
        assert!(
            stock
                .rust_internal_thread_treatment_original_stock_topology()
                .is_none()
        );
        assert_eq!(stock.volume_mm3().unwrap(), original);
    }
    #[test]
    fn actual_occurrence_transform_preserves_source_and_export_refuses_without_touching_destination()
     {
        let policy = GeometryPolicy::default();
        let stock = fixture_stock(Rational::new(4, 1).unwrap(), PlaneFrame::default());
        let original = stock
            .rust_cut_authored_internal_thread_treatments_with_interrupt(
                cut(PlaneFrame::default(), false, 1),
                EntityId::new(),
                policy,
                None,
            )
            .unwrap();
        let transform = RigidTransform {
            rotation: [[0.8, -0.36, 0.48], [0.6, 0.48, -0.64], [0., 0.8, 0.6]],
            translation_mm: [1000., -2000., 3000.],
        };
        let moved = original.transformed(transform).unwrap();
        assert_eq!(
            moved
                .rust_internal_thread_treatment_data()
                .unwrap()
                .material,
            original
                .rust_internal_thread_treatment_data()
                .unwrap()
                .material
        );
        assert_eq!(
            moved.rust_internal_thread_treatment_original_stock_topology(),
            original.rust_internal_thread_treatment_original_stock_topology()
        );
        let interval = moved
            .rust_internal_thread_treatment_volume_bounds()
            .unwrap()
            .unwrap();
        assert!(
            interval.lower <= moved.volume_mm3().unwrap()
                && moved.volume_mm3().unwrap() <= interval.upper
        );
        let mesh = moved.manufacturing_mesh_with_receipt(options()).unwrap();
        assert!(mesh.edge_pairs > 0);
        let bounds = moved.bounds().unwrap();
        assert!(
            mesh.mesh
                .vertices
                .as_chunks::<3>()
                .0
                .iter()
                .all(|p| (0..3).all(|i| bounds.min[i] <= p[i] && p[i] <= bounds.max[i]))
        );
        let destination = std::env::temp_dir().join(format!(
            "multicad-internal-treatment-{}.step",
            EntityId::new()
        ));
        std::fs::write(&destination, b"preserved previous file").unwrap();
        let refused = moved.write_step(&destination).unwrap_err();
        assert_eq!(refused.kind, GeometryErrorKind::Unsupported);
        assert_eq!(
            std::fs::read(&destination).unwrap(),
            b"preserved previous file"
        );
        std::fs::remove_file(destination).unwrap();
    }
    #[test]
    fn malformed_requests_source_cancellation_limits_and_large_ratio_stock_keep_real_guards() {
        let policy = GeometryPolicy::default();
        let frame = PlaneFrame::default();
        let stock = fixture_stock(Rational::new(4, 1).unwrap(), frame);
        let request = cut(frame, false, 1);
        for feature in [EntityId::nil(), EntityId::from_uuid_bits(7)] {
            assert_eq!(
                stock
                    .rust_cut_authored_internal_thread_treatments_with_interrupt(
                        request, feature, policy, None
                    )
                    .err()
                    .unwrap()
                    .kind,
                GeometryErrorKind::Input
            );
        }
        let box_stock = Solid::box_solid([5., 5., 4.], policy).unwrap();
        assert_eq!(
            box_stock
                .rust_cut_authored_internal_thread_treatments_with_interrupt(
                    request,
                    EntityId::new(),
                    policy,
                    None
                )
                .err()
                .unwrap()
                .kind,
            GeometryErrorKind::Unsupported
        );
        assert!(
            stock
                .rust_cut_authored_internal_thread_treatments_with_interrupt(
                    request,
                    EntityId::new(),
                    GeometryPolicy {
                        minimum_feature_mm: 0.01,
                        ..policy
                    },
                    None
                )
                .is_err()
        );
        for bad in [
            InternalThreadTreatmentCut {
                minimum_wall_mm: 0.,
                ..request
            },
            InternalThreadTreatmentCut {
                profile: HelicalThread {
                    radial_taper: 0.01,
                    ..request.profile
                },
                ..request
            },
            InternalThreadTreatmentCut {
                intent: source::Intent {
                    entry_depth_mm: None,
                    terminal: None,
                },
                ..request
            },
            InternalThreadTreatmentCut {
                profile: HelicalThread {
                    major_radius_mm: 1e6,
                    ..request.profile
                },
                ..request
            },
        ] {
            assert!(
                stock
                    .rust_cut_authored_internal_thread_treatments_with_interrupt(
                        bad,
                        EntityId::new(),
                        policy,
                        None
                    )
                    .is_err()
            );
        }
        let calls = std::cell::Cell::new(0);
        let count = || {
            calls.set(calls.get() + 1);
            false
        };
        let body = stock
            .rust_cut_authored_internal_thread_treatments_with_interrupt(
                request,
                EntityId::new(),
                policy,
                Some(&count),
            )
            .unwrap();
        let total = calls.get();
        for stop in 0..total {
            calls.set(0);
            let cancel = || {
                calls.set(calls.get() + 1);
                calls.get() > stop
            };
            assert_eq!(
                stock
                    .rust_cut_authored_internal_thread_treatments_with_interrupt(
                        request,
                        EntityId::new(),
                        policy,
                        Some(&cancel)
                    )
                    .err()
                    .unwrap()
                    .kind,
                GeometryErrorKind::Cancelled
            );
        }
        assert_eq!(
            body.manufacturing_mesh_with_receipt(MeshOptions {
                maximum_triangles: 32,
                ..options()
            })
            .err()
            .unwrap()
            .kind,
            GeometryErrorKind::Limit
        );
        assert!(
            metric(PlaneFrame {
                u: [0.; 3],
                ..frame
            })
            .is_err()
        );
        let zero = Rational::new(0, 1).unwrap();
        let radius = Rational::new(900_000_000_000_000_001, 300_000_000_000_000_000).unwrap();
        let actual = AxialCylinder::new(
            7,
            Circle::new([zero; 2], radius).unwrap(),
            zero,
            Rational::new(4, 1).unwrap(),
        )
        .unwrap();
        let source_stock = Solid::build_axial_cylinder(actual, frame, None, policy).unwrap();
        let result = source_stock
            .rust_cut_authored_internal_thread_treatments_with_interrupt(
                request,
                EntityId::new(),
                policy,
                None,
            )
            .unwrap();
        assert_eq!(
            result
                .rust_internal_thread_treatment_data()
                .unwrap()
                .material
                .stock()
                .circle()
                .radius(),
            radius
        );
        assert!(result.manufacturing_mesh_with_receipt(options()).is_ok());
        assert_eq!(
            stock.volume_mm3().unwrap(),
            fixture_stock(Rational::new(4, 1).unwrap(), frame)
                .volume_mm3()
                .unwrap()
        );
    }
    #[test]
    fn axial_rational_coordinates_and_framed_literal_origins_have_distinct_authority() {
        let policy = GeometryPolicy::default();
        let zero = Rational::new(0, 1).unwrap();
        for (center, bottom, expected) in [
            ([Rational::new(33, 10).unwrap(), zero], zero, false),
            ([zero; 2], Rational::new(33, 10).unwrap(), false),
            (
                [
                    Rational::new(200_000_000_000_000_001, 1_i64 << 58).unwrap(),
                    zero,
                ],
                zero,
                false,
            ),
            (
                [Rational::new(13, 4).unwrap(), Rational::new(-7, 2).unwrap()],
                Rational::new(9, 4).unwrap(),
                true,
            ),
        ] {
            let actual = AxialCylinder::new(
                7,
                Circle::new(center, Rational::new(5, 2).unwrap()).unwrap(),
                bottom,
                Rational::new(4, 1).unwrap(),
            )
            .unwrap();
            let frame = PlaneFrame {
                origin: [center[0].as_f64(), center[1].as_f64(), bottom.as_f64()],
                ..PlaneFrame::default()
            };
            let stock = Solid::from_axial_cylinder(actual, frame, policy).unwrap();
            let result = stock.rust_cut_authored_internal_thread_treatments_with_interrupt(
                cut(frame, false, 1),
                EntityId::new(),
                policy,
                None,
            );
            assert_eq!(result.is_ok(), expected);
            if let Ok(body) = result {
                assert!(body.manufacturing_mesh_with_receipt(options()).is_ok());
            }
            assert_eq!(stock.axial_cylinder().unwrap(), &actual);
        }
        let actual = AxialCylinder::new(
            7,
            Circle::new([zero; 2], Rational::new(5, 2).unwrap()).unwrap(),
            zero,
            Rational::new(4, 1).unwrap(),
        )
        .unwrap();
        let frame = PlaneFrame {
            origin: [3.3, -2.7, 0.1],
            u: [0., 1., 0.],
            v: [0., 0., 1.],
        };
        let stock = Solid::from_axial_cylinder_sweep(
            actual,
            frame,
            [0.; 2],
            false,
            multicad_core::references::Origin {
                source: EntityId::new(),
                role: "stock.circle".into(),
            },
            policy,
        )
        .unwrap();
        assert!(stock.axial_cylinder().is_none());
        let result = stock
            .rust_cut_authored_internal_thread_treatments_with_interrupt(
                cut(frame, true, 2),
                EntityId::new(),
                policy,
                None,
            )
            .unwrap();
        assert_eq!(
            result.rust_internal_thread_treatment_data().unwrap().frame,
            frame
        );
        assert!(result.manufacturing_mesh_with_receipt(options()).is_ok());
    }
    #[test]
    fn whole_stock_extent_wall_and_volume_precision_are_reproved_for_literal_placement() {
        let policy = GeometryPolicy::default();
        let frame = PlaneFrame::default();
        let stock_radius = |radius| {
            let zero = Rational::new(0, 1).unwrap();
            let source = AxialCylinder::new(
                7,
                Circle::new([zero; 2], Rational::new(radius, 1).unwrap()).unwrap(),
                zero,
                Rational::new(4, 1).unwrap(),
            )
            .unwrap();
            Solid::from_axial_cylinder(source, frame, policy).unwrap()
        };
        let wide = stock_radius(1_000_000)
            .rust_cut_authored_internal_thread_treatments_with_interrupt(
                cut(frame, false, 1),
                EntityId::new(),
                policy,
                None,
            )
            .unwrap();
        let nearly_proper = RigidTransform {
            rotation: [[1. + 4e-11, 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            ..RigidTransform::default()
        };
        assert_eq!(
            wide.transformed(nearly_proper).err().unwrap().kind,
            GeometryErrorKind::Unsupported
        );
        let stock = fixture_stock(Rational::new(4, 1).unwrap(), frame);
        let tight = stock
            .rust_cut_authored_internal_thread_treatments_with_interrupt(
                InternalThreadTreatmentCut {
                    minimum_wall_mm: 0.75 - 1e-12,
                    ..cut(frame, false, 1)
                },
                EntityId::new(),
                policy,
                None,
            )
            .unwrap();
        let contraction = RigidTransform {
            rotation: [[1. - 4e-11, 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            ..RigidTransform::default()
        };
        assert_eq!(
            tight.transformed(contraction).err().unwrap().message,
            "Placed stock no longer proves the requested material wall"
        );
        let mut placed_floor = false;
        for radius in [300, 500, 700, 900, 1100, 1300, 1500, 1800, 2100] {
            let body = stock_radius(radius)
                .rust_cut_authored_internal_thread_treatments_with_interrupt(
                    cut(frame, false, 1),
                    EntityId::new(),
                    policy,
                    None,
                )
                .unwrap();
            if let Err(error) = body.rust_internal_thread_treatment_volume_bounds() {
                placed_floor |=
                    error.message == "Placed material volume reaches the numerical precision floor";
            }
        }
        assert!(
            placed_floor,
            "actual placed source must refuse unrepresentable volume width"
        );
        let clock = PlaneFrame {
            u: [0.8, 0.6, 0.],
            v: [-0.6, 0.8, 0.],
            ..frame
        };
        let actual = *stock.axial_cylinder().unwrap();
        let wcs_stock = Solid::from_axial_cylinder(actual, clock, policy).unwrap();
        assert_eq!(
            wcs_stock
                .rust_cut_authored_internal_thread_treatments_with_interrupt(
                    cut(clock, false, 1),
                    EntityId::new(),
                    policy,
                    None
                )
                .err()
                .unwrap()
                .kind,
            GeometryErrorKind::Unsupported
        );
    }
}
