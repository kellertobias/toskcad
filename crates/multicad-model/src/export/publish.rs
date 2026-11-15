use super::*;
use multicad_core::{
    commands::Actor,
    exchange::{WrittenFile, filename_token},
    project_file::ContentHash,
    references::TopologySnapshot,
};
use multicad_kernel::{Solid, rust_kernel::step_read};
use serde::{Deserialize, Serialize};
use std::{
    collections::BTreeSet,
    fs::{self, OpenOptions},
    io::{Read, Write},
    path::Path,
};
const MANIFEST: &str = "multicad-export.json";
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum StepVerificationMethod {
    DecodedBoundary,
    DecodedBoundaryEnclosure,
    SourceKnownBoundary,
    MixedMemberBoundary,
    DecodedSupportHealing,
}
/// Mass of independently decoded supports after the reported bounded sewing.
/// Raw-file occupied mass and healed-cap derivative bounds are not asserted.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StepHealedMemberReceipt {
    pub member_index: u32,
    pub healed_volume_bounds_mm3: [f64; 2],
    pub source_volume_bounds_mm3: [f64; 2],
    pub source_to_raw_position_error_mm: f64,
    pub sewing_correspondence_bound_mm: f64,
    pub source_to_healed_position_error_mm: f64,
    pub source_support_derivative_error_mm: [f64; 2],
    pub original_coedge_derivative_error_mm: f64,
    pub verified_tolerance_mm: f64,
    pub matched_entities: usize,
    pub literal_binding_work: usize,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StepRetainedMemberReceipt {
    pub member_index: u32,
    pub decoded_volume_bounds_mm3: [f64; 2],
    pub source_volume_bounds_mm3: [f64; 2],
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StepSourceCorrespondenceReceipt {
    pub member_index: u32,
    pub maximum_position_error_mm: f64,
    /// Absolute error in the original whole-row parameter; not normal angle.
    pub row_derivative_error_mm: f64,
    pub ruling_derivative_error_mm: f64,
}
/// A source-known member has a checked actual finite boundary, but its
/// volume bounds come from the authored source rather than independent file mass.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StepSourceKnownMemberReceipt {
    pub member_index: u32,
    pub source_volume_bounds_mm3: [f64; 2],
    pub verified_tolerance_mm: f64,
    pub maximum_polynomial_trim_error_mm: f64,
    pub analytic_closure_tolerance_mm: f64,
}
/// Distinguishes decoded file mass from certified authored-source mass.
/// Missing decoded mass is explicit; source support tolerance is not a
/// substitute for an independently integrated file-volume error bound.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StepVerificationReceipt {
    pub method: StepVerificationMethod,
    pub decoded_volume_mm3: Option<f64>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub retained_members: Vec<StepRetainedMemberReceipt>,
    /// Separate existing decoded scalar contracts; never treated as exact
    /// zero-width contributions to a retained whole-file enclosure.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub decoded_scalar_members: Vec<(u32, f64)>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub source_correspondence: Vec<StepSourceCorrespondenceReceipt>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub source_known_members: Vec<StepSourceKnownMemberReceipt>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub healed_members: Vec<StepHealedMemberReceipt>,
    pub source_volume_bounds_mm3: Option<[f64; 2]>,
    /// Actual full-checker parameter, not a measured combined boundary error.
    pub verified_tolerance_mm: Option<f64>,
    /// Measured polynomial trim-to-support binding only.
    pub maximum_polynomial_trim_error_mm: Option<f64>,
    /// Separate analytic-closure check budget; not added to source mass.
    pub analytic_closure_tolerance_mm: Option<f64>,
}
impl StepVerificationReceipt {
    pub fn decoded(volume: f64) -> Self {
        Self {
            method: StepVerificationMethod::DecodedBoundary,
            decoded_volume_mm3: Some(volume),
            retained_members: vec![],
            decoded_scalar_members: vec![],
            source_correspondence: vec![],
            source_known_members: vec![],
            healed_members: vec![],
            source_volume_bounds_mm3: None,
            verified_tolerance_mm: None,
            maximum_polynomial_trim_error_mm: None,
            analytic_closure_tolerance_mm: None,
        }
    }
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Manifest {
    schema: u32,
    provider: String,
    set: EntityId,
    source_hash: ContentHash,
    plan_hash: ContentHash,
    actor: Actor,
    files: Vec<WrittenFile>,
}
fn io(e: impl ToString) -> ExportError {
    failure(ErrorKind::Io, None, e)
}
fn verify(
    exporter: &Exporter,
    built: &Built,
    artifact: &ExportArtifact,
    path: &Path,
) -> Result<(), ExportError> {
    if let Some(evidence) = &artifact.schedule {
        let expected = built
            .generated
            .get(&artifact.output.filename)
            .ok_or_else(|| io("Current schedule output is missing"))?;
        let written = fs::read(path).map_err(io)?;
        if written != *expected || ContentHash::of(&written) != evidence.csv_hash {
            return Err(io(
                "Written schedule CSV differs from the generated schedule",
            ));
        }
        return verify_csv(&written, evidence)
            .map_err(|e| failure(ErrorKind::Invalid, Some(artifact.producer), e));
    }
    if let Some(evidence) = &artifact.drawing {
        let expected = built
            .generated
            .get(&artifact.output.filename)
            .ok_or_else(|| io("Current drawing output is missing"))?;
        let written = fs::read(path).map_err(io)?;
        if written != *expected {
            return Err(io("Written drawing bytes differ from the generated output"));
        }
        let text = String::from_utf8_lossy(&written);
        let structure_ok = if artifact.output.filename.ends_with(".pdf") {
            written.starts_with(b"%PDF-1.4")
                && text.ends_with("%%EOF\n")
                && text.matches("/Type /Page ").count() == evidence.pages.len()
        } else {
            text.starts_with("0\nSECTION\n2\nHEADER")
                && text.ends_with("0\nEOF\n")
                && evidence.pages.len() == 1
        };
        if !structure_ok {
            return Err(failure(
                ErrorKind::Invalid,
                Some(artifact.producer),
                "Drawing readback structure or page count differs from the plan",
            ));
        }
        return Ok(());
    }
    if let Some(evidence) = &artifact.mesh {
        let expected = built
            .generated
            .get(&artifact.output.filename)
            .ok_or_else(|| io("Current exact mesh output is missing"))?;
        return verify_mesh(&fs::read(path).map_err(io)?, expected, artifact, evidence);
    }
    if artifact.planar.is_some() {
        let expected = built
            .generated
            .get(&artifact.output.filename)
            .ok_or_else(|| io("Current exact planar output is missing"))?;
        let (_, _) = hash_file(path)?;
        return dxf::verify(&fs::read(path).map_err(io)?, expected).map_err(io);
    }
    let tolerance = exporter.compiler.policy.linear_tolerance_mm;
    let actual = if built.backend == GeometryBackend::Rust {
        let model = built
            .model
            .as_ref()
            .ok_or_else(|| io("Missing source model"))?;
        let topology = rust_step_topology(model, artifact)?;
        let solid = artifact_solid(model, artifact)?;
        if let Some(proof) = verify_partitioned_step(
            path,
            solid,
            (
                artifact.solid_count,
                artifact.volume_mm3,
                artifact.bounds_mm,
            ),
            artifact.producer,
            tolerance,
        ) {
            return proof.map(|_| ());
        }
        if let Some(proof) = verify_retained_step(
            path,
            solid,
            (
                artifact.solid_count,
                artifact.volume_mm3,
                artifact.bounds_mm,
            ),
            artifact.producer,
            tolerance,
        ) {
            return proof.map(|_| ());
        }
        if let Some(proof) = verify_source_known_step(
            path,
            solid,
            (
                artifact.solid_count,
                artifact.volume_mm3,
                artifact.bounds_mm,
            ),
            artifact.producer,
            tolerance,
        ) {
            return proof.map(|_| ());
        }
        read_rust_step(path, solid, topology, artifact.producer, tolerance)?
    } else {
        let read = Solid::read_step(path, exporter.compiler.policy).map_err(io)?;
        let bounds = read.bounds().map_err(io)?;
        (
            read.solid_count().map_err(io)?,
            read.volume_mm3().map_err(io)?,
            [bounds.min, bounds.max],
        )
    };
    verify_step_dimensions(
        artifact.producer,
        (
            artifact.solid_count,
            artifact.volume_mm3,
            artifact.bounds_mm,
        ),
        actual,
        tolerance,
    )
}

type StepMetrics = multicad_kernel::StepCompoundMetrics;

fn verify_step_dimensions(
    producer: EntityId,
    (expected_count, expected_volume, [a, b]): StepMetrics,
    (count, volume, bounds): StepMetrics,
    tolerance: f64,
) -> Result<(), ExportError> {
    let [x, y, z] = std::array::from_fn::<_, 3, _>(|i| b[i] - a[i]);
    let volume_tolerance = (expected_volume * 1e-8).max(tolerance * 2. * (x * y + y * z + z * x));
    if count != expected_count
        || (volume - expected_volume).abs() > volume_tolerance
        || bounds
            .into_iter()
            .flatten()
            .zip(a.into_iter().chain(b))
            .any(|(x, y)| (x - y).abs() > tolerance * 2.)
    {
        return Err(failure(
            ErrorKind::Invalid,
            Some(producer),
            "STEP readback differs from the current solid's count, volume or world bounds",
        ));
    }
    Ok(())
}

/// Decoded-mass compatibility helper. Source-known-only verification requires
/// the discriminated receipt API and never returns expected mass here.
pub fn verify_rust_step_shape(
    path: &Path,
    solid: &Solid,
    shape: &crate::ShapeSnapshot,
    tolerance: f64,
) -> Result<f64, ExportError> {
    if solid.rust_thread_entry_chamfer_data().is_some()
        || solid.has_retained_radial_source()
        || solid.has_partitioned_step_source()
        || solid.has_literal_npt_source()
    {
        return Err(failure(
            ErrorKind::Unsupported,
            Some(shape.producer),
            "STEP mass needs a discriminated verification receipt; use the verification receipt API",
        ));
    }
    let actual = read_rust_step(path, solid, &shape.topology, shape.producer, tolerance)?;
    verify_step_dimensions(
        shape.producer,
        (
            shape.solid_count,
            shape.volume_mm3,
            [shape.bounds.min, shape.bounds.max],
        ),
        actual,
        tolerance,
    )?;
    Ok(actual.1)
}

/// Authored export verification with an explicit mass provenance contract.
pub fn verify_rust_step_shape_receipt(
    path: &Path,
    solid: &Solid,
    shape: &crate::ShapeSnapshot,
    tolerance: f64,
) -> Result<StepVerificationReceipt, ExportError> {
    if !tolerance.is_finite() || tolerance <= 0. {
        return Err(failure(
            ErrorKind::Invalid,
            Some(shape.producer),
            "STEP verification needs a finite positive tolerance",
        ));
    }
    if let Some(proof) = verify_partitioned_step(
        path,
        solid,
        (
            shape.solid_count,
            shape.volume_mm3,
            [shape.bounds.min, shape.bounds.max],
        ),
        shape.producer,
        tolerance,
    ) {
        return proof;
    }
    if let Some(proof) = verify_retained_step(
        path,
        solid,
        (
            shape.solid_count,
            shape.volume_mm3,
            [shape.bounds.min, shape.bounds.max],
        ),
        shape.producer,
        tolerance,
    ) {
        return proof;
    }
    if let Some(proof) = verify_source_known_step(
        path,
        solid,
        (
            shape.solid_count,
            shape.volume_mm3,
            [shape.bounds.min, shape.bounds.max],
        ),
        shape.producer,
        tolerance,
    ) {
        return proof;
    }
    verify_rust_step_shape(path, solid, shape, tolerance).map(StepVerificationReceipt::decoded)
}
fn verify_partitioned_step(
    path: &Path,
    solid: &Solid,
    expected: StepMetrics,
    producer: EntityId,
    tolerance: f64,
) -> Option<Result<StepVerificationReceipt, ExportError>> {
    if !solid.has_partitioned_step_source() {
        return None;
    }
    Some((|| {
        if !tolerance.is_finite()
            || tolerance <= 0.
            || !expected.1.is_finite()
            || expected.1 <= 0.
            || expected.2.into_iter().flatten().any(|v| !v.is_finite())
        {
            return Err(failure(
                ErrorKind::Invalid,
                Some(producer),
                "Insert source snapshot needs finite bounds, positive mass and verification tolerance",
            ));
        }
        let bounds = solid.bounds().map_err(io)?;
        verify_step_dimensions(
            producer,
            expected,
            (
                solid.solid_count().map_err(io)?,
                solid.volume_mm3().map_err(io)?,
                [bounds.min, bounds.max],
            ),
            tolerance,
        )?;
        let input = read_step_source(path, producer)?;
        let proof = solid
            .verify_partitioned_step_source(&input, tolerance, None)
            .expect("actual retained insert source")
            .map_err(|e| {
                failure(
                    if e.kind == multicad_kernel::GeometryErrorKind::Limit {
                        ErrorKind::Budget
                    } else {
                        ErrorKind::Invalid
                    },
                    Some(producer),
                    format!("Actual insert file proof failed:{e}"),
                )
            })?;
        Ok(StepVerificationReceipt {
            method: if !proof.healed_members.is_empty() {
                if proof.source_known_members.is_empty()
                    && proof.members.retained_members.is_empty()
                    && proof.members.decoded_scalar_members.is_empty()
                {
                    StepVerificationMethod::DecodedSupportHealing
                } else {
                    StepVerificationMethod::MixedMemberBoundary
                }
            } else if proof.source_known_members.is_empty() {
                StepVerificationMethod::DecodedBoundaryEnclosure
            } else {
                StepVerificationMethod::MixedMemberBoundary
            },
            decoded_volume_mm3: None,
            retained_members: proof
                .members
                .retained_members
                .into_iter()
                .map(|member| StepRetainedMemberReceipt {
                    member_index: member.member_index,
                    decoded_volume_bounds_mm3: member.decoded_volume_bounds_mm3,
                    source_volume_bounds_mm3: member.source_volume_bounds_mm3,
                })
                .collect(),
            decoded_scalar_members: proof.members.decoded_scalar_members,
            source_correspondence: proof
                .insert_correspondence
                .into_iter()
                .map(|(member_index, p)| StepSourceCorrespondenceReceipt {
                    member_index,
                    maximum_position_error_mm: p.maximum_position_error_mm,
                    row_derivative_error_mm: p.maximum_row_derivative_error_mm,
                    ruling_derivative_error_mm: p.maximum_ruling_derivative_error_mm,
                })
                .collect(),
            source_known_members: proof
                .source_known_members
                .into_iter()
                .map(|(member_index, p)| StepSourceKnownMemberReceipt {
                    member_index,
                    source_volume_bounds_mm3: p.source_volume_bounds_mm3,
                    verified_tolerance_mm: p.verified_tolerance_mm,
                    maximum_polynomial_trim_error_mm: p.maximum_polynomial_trim_error_mm,
                    analytic_closure_tolerance_mm: p.analytic_closure_tolerance_mm,
                })
                .collect(),
            source_volume_bounds_mm3: None,
            verified_tolerance_mm: Some(tolerance),
            healed_members: proof
                .healed_members
                .into_iter()
                .map(|(member_index, p)| StepHealedMemberReceipt {
                    member_index,
                    healed_volume_bounds_mm3: p.healed_volume_bounds_mm3,
                    source_volume_bounds_mm3: p.source_volume_bounds_mm3,
                    source_to_raw_position_error_mm: p.source_to_raw_position_error_mm,
                    sewing_correspondence_bound_mm: p.sewing_correspondence_bound_mm,
                    source_to_healed_position_error_mm: p.source_to_healed_position_error_mm,
                    source_support_derivative_error_mm: p.source_support_derivative_error_mm,
                    original_coedge_derivative_error_mm: p.original_coedge_derivative_error_mm,
                    verified_tolerance_mm: p.verified_tolerance_mm,
                    matched_entities: p.matched_entities,
                    literal_binding_work: p.literal_binding_work,
                })
                .collect(),
            maximum_polynomial_trim_error_mm: None,
            analytic_closure_tolerance_mm: None,
        })
    })())
}
fn verify_retained_step(
    path: &Path,
    solid: &Solid,
    expected: StepMetrics,
    producer: EntityId,
    tolerance: f64,
) -> Option<Result<StepVerificationReceipt, ExportError>> {
    if !solid.has_retained_radial_source() && !solid.has_literal_npt_source() {
        return None;
    }
    Some((|| {
        if !tolerance.is_finite() || tolerance <= 0. {
            return Err(failure(
                ErrorKind::Invalid,
                Some(producer),
                "STEP verification needs a finite positive tolerance",
            ));
        }
        if !expected.1.is_finite()
            || expected.1 <= 0.
            || expected
                .2
                .into_iter()
                .flatten()
                .any(|value| !value.is_finite())
        {
            return Err(failure(
                ErrorKind::Invalid,
                Some(producer),
                "Retained source snapshot needs finite bounds and positive finite mass",
            ));
        }
        // Check source snapshot freshness separately from file mass. The
        // source scalar is a display value, not a decoded exact-volume proof.
        let bounds = solid.bounds().map_err(io)?;
        verify_step_dimensions(
            producer,
            expected,
            (
                solid.solid_count().map_err(io)?,
                solid.volume_mm3().map_err(io)?,
                [bounds.min, bounds.max],
            ),
            tolerance,
        )?;
        let input = read_step_source(path, producer)?;
        let proof = if solid.has_literal_npt_source() {
            solid.verify_literal_npt_step_compound(&input, tolerance)
        } else {
            solid.verify_retained_step_source(&input, tolerance)
        }
        .expect("typed retained source")
        .map_err(|cause| {
            failure(
                if cause == step_read::ReadError::Limit {
                    ErrorKind::Budget
                } else {
                    ErrorKind::Invalid
                },
                Some(producer),
                format!("Retained STEP boundary proof failed: {cause:?}"),
            )
        })?;
        if proof.solid_count != expected.0 {
            return Err(failure(
                ErrorKind::Invalid,
                Some(producer),
                "Retained STEP member count differs from the current source",
            ));
        }
        Ok(StepVerificationReceipt {
            method: StepVerificationMethod::DecodedBoundaryEnclosure,
            decoded_volume_mm3: None,
            retained_members: proof
                .retained_members
                .into_iter()
                .map(|member| StepRetainedMemberReceipt {
                    member_index: member.member_index,
                    decoded_volume_bounds_mm3: member.decoded_volume_bounds_mm3,
                    source_volume_bounds_mm3: member.source_volume_bounds_mm3,
                })
                .collect(),
            decoded_scalar_members: proof.decoded_scalar_members,
            source_correspondence: vec![],
            source_known_members: vec![],
            healed_members: vec![],
            source_volume_bounds_mm3: None,
            verified_tolerance_mm: Some(tolerance),
            maximum_polynomial_trim_error_mm: None,
            analytic_closure_tolerance_mm: None,
        })
    })())
}
fn verify_source_known_step(
    path: &Path,
    solid: &Solid,
    expected: StepMetrics,
    producer: EntityId,
    tolerance: f64,
) -> Option<Result<StepVerificationReceipt, ExportError>> {
    solid.rust_thread_entry_chamfer_data()?;
    Some((|| {
        if !tolerance.is_finite() || tolerance <= 0. {
            return Err(failure(
                ErrorKind::Invalid,
                Some(producer),
                "STEP verification needs a finite positive tolerance",
            ));
        }
        let source = read_step_source(path, producer)?;
        let proof = solid
            .verify_entry_chamfer_step_with_tolerance(&source, tolerance)
            .expect("typed entry source")
            .map_err(|cause| {
                failure(
                    if cause.kind == multicad_kernel::GeometryErrorKind::Limit {
                        ErrorKind::Budget
                    } else {
                        ErrorKind::Invalid
                    },
                    Some(producer),
                    format!("STEP authored boundary proof failed: {cause}"),
                )
            })?;
        if proof.verified_tolerance_mm > tolerance
            || !proof.maximum_polynomial_trim_error_mm.is_finite()
            || proof.maximum_polynomial_trim_error_mm > tolerance
        {
            return Err(failure(
                ErrorKind::Invalid,
                Some(producer),
                "STEP boundary binding exceeds the requested tolerance",
            ));
        }
        let (count, volume, bounds) = expected;
        let current = solid.bounds().map_err(io)?;
        if count != 1
            || !volume.is_finite()
            || volume < proof.source_volume_bounds_mm3[0]
            || volume > proof.source_volume_bounds_mm3[1]
            || bounds
                .into_iter()
                .flatten()
                .zip(current.min.into_iter().chain(current.max))
                .any(|(a, b)| !a.is_finite() || (a - b).abs() > tolerance * 2.)
        {
            return Err(failure(
                ErrorKind::Invalid,
                Some(producer),
                "Current source snapshot differs from the certified authored STEP source",
            ));
        }
        Ok(StepVerificationReceipt {
            method: StepVerificationMethod::SourceKnownBoundary,
            decoded_volume_mm3: None,
            retained_members: vec![],
            decoded_scalar_members: vec![],
            source_correspondence: vec![],
            source_known_members: vec![],
            healed_members: vec![],
            source_volume_bounds_mm3: Some(proof.source_volume_bounds_mm3),
            verified_tolerance_mm: Some(proof.verified_tolerance_mm),
            maximum_polynomial_trim_error_mm: Some(proof.maximum_polynomial_trim_error_mm),
            analytic_closure_tolerance_mm: Some(proof.analytic_closure_tolerance_mm),
        })
    })())
}

fn read_step_source(path: &Path, producer: EntityId) -> Result<String, ExportError> {
    let mut bytes = Vec::new();
    fs::File::open(path)
        .map_err(io)?
        .take(step_read::MAX_PLANAR_STEP_BYTES as u64 + 1)
        .read_to_end(&mut bytes)
        .map_err(io)?;
    if bytes.len() > step_read::MAX_PLANAR_STEP_BYTES {
        return Err(failure(
            ErrorKind::Budget,
            Some(producer),
            "STEP readback exceeds the file budget",
        ));
    }
    String::from_utf8(bytes).map_err(|e| failure(ErrorKind::Invalid, Some(producer), e.to_string()))
}

fn rust_step_topology<'a>(
    model: &'a crate::Model,
    artifact: &ExportArtifact,
) -> Result<&'a TopologySnapshot, ExportError> {
    (|| {
        let OutputSubject::Part(part) = &artifact.output.subject else {
            return None;
        };
        if let Some((&root, members)) = part.occurrence_path.split_first() {
            let leaf = model.leaf_snapshot(&multicad_core::model::OccurrencePath {
                root,
                members: members.to_vec(),
            })?;
            (root == artifact.solid
                && leaf.component == part.part
                && leaf.producer == artifact.producer)
                .then_some(&leaf.topology)
        } else {
            model
                .snapshot()
                .shapes
                .get(&artifact.solid)
                .map(|shape| &shape.topology)
        }
    })()
    .ok_or_else(|| {
        failure(
            ErrorKind::Invalid,
            Some(artifact.producer),
            "Rust STEP verification requires the exact source topology",
        )
    })
}

fn read_rust_step(
    path: &Path,
    solid: &Solid,
    topology: &TopologySnapshot,
    producer: EntityId,
    tolerance: f64,
) -> Result<StepMetrics, ExportError> {
    let invalid = |message: String| failure(ErrorKind::Invalid, Some(producer), message);
    let mut bytes = Vec::new();
    fs::File::open(path)
        .map_err(io)?
        .take(step_read::MAX_PLANAR_STEP_BYTES as u64 + 1)
        .read_to_end(&mut bytes)
        .map_err(io)?;
    if bytes.len() > step_read::MAX_PLANAR_STEP_BYTES {
        return Err(failure(
            ErrorKind::Budget,
            Some(producer),
            "STEP readback exceeds the file budget",
        ));
    }
    let source = std::str::from_utf8(&bytes).map_err(|e| invalid(e.to_string()))?;
    let readback_error = |e| {
        let kind = match e {
            step_read::ReadError::Limit => ErrorKind::Budget,
            _ => ErrorKind::Invalid,
        };
        failure(
            kind,
            Some(producer),
            format!("Rust STEP readback rejected the geometry: {e:?}"),
        )
    };
    if let Some(proof) = solid.verify_authored_step_compound(source, tolerance) {
        return proof.map_err(readback_error);
    }
    if let Some(proof) = solid
        .verify_helical_step(source)
        .or_else(|| solid.verify_normal_relief_step(source))
        .or_else(|| solid.verify_herringbone_step(source))
        .or_else(|| solid.verify_bevel_step(source))
        .or_else(|| solid.verify_worm_step(source))
    {
        let volume = proof.map_err(readback_error)?;
        // The verified boundary gives the source bounds within its uniform
        // position error. Volume comes from the independently decoded boundary.
        let bounds = solid.bounds().map_err(io)?;
        return Ok((1, volume, [bounds.min, bounds.max]));
    }
    let decoded = step_read::read_analytic_compound_with_helical_cuts(source, tolerance)
        .map_err(readback_error)?;
    if !solid.step_readback_matches(topology, &decoded, tolerance * 2.) {
        return Err(invalid(
            "STEP readback differs from the source member topology".into(),
        ));
    }
    Ok(decoded.into_iter().fold(
        (0, 0., [[f64::INFINITY; 3], [f64::NEG_INFINITY; 3]]),
        |(count, volume, [min, max]), (member, _)| {
            let [lo, hi] = member.bounds_mm;
            (
                count + member.solid_count,
                volume + member.volume_mm3,
                [
                    std::array::from_fn(|i| min[i].min(lo[i])),
                    std::array::from_fn(|i| max[i].max(hi[i])),
                ],
            )
        },
    ))
}

fn verify_mesh(
    written: &[u8],
    expected: &[u8],
    artifact: &ExportArtifact,
    evidence: &MeshEvidence,
) -> Result<(), ExportError> {
    if written != expected {
        return Err(io("Written mesh bytes differ from the generated output"));
    }
    let name = &artifact.output.filename;
    // STL and GLB store single-precision coordinates; STL keeps no objects.
    let (read, single, flat) = if name.ends_with(".obj") {
        (super::mesh::read_obj(written), false, false)
    } else if name.ends_with(".stl") {
        (super::mesh::read_stl(written), true, true)
    } else if name.ends_with(".glb") {
        (super::mesh::read_glb(written), true, false)
    } else {
        (super::mesh::read_three_mf(written), false, false)
    };
    let read = read.map_err(io)?;
    let bounds = read
        .bounds_mm()
        .ok_or_else(|| io("Mesh readback is empty"))?;
    let [a, b] = artifact.bounds_mm;
    let extent = a.iter().chain(&b).fold(1f64, |m, v| m.max(v.abs()));
    let precision = if single { extent * 1e-6 } else { 0. };
    let slack = evidence.chord_deviation_mm * 2. + 1e-6 + precision;
    let relative = if single { 1e-5 } else { 1e-6 };
    let volume_ok = (read.volume_mm3() - evidence.mesh_volume_mm3).abs()
        <= evidence.mesh_volume_mm3.abs() * relative + 1e-6;
    // The tessellation must also enclose the exact solid's volume, up to the
    // chord deviation over its surface; a mis-oriented or open mesh fails.
    let [x, y, z] = std::array::from_fn::<_, 3, _>(|i| b[i] - a[i]);
    let chord_allowance = 2. * (x * y + y * z + z * x) * evidence.chord_deviation_mm * 1.5
        + artifact.volume_mm3 * 1e-6;
    let exact_ok = (read.volume_mm3() - artifact.volume_mm3).abs() <= chord_allowance;
    if (!flat && (read.objects != evidence.objects || read.items != evidence.build_items))
        || read.triangles != evidence.triangles
        || !volume_ok
        || !exact_ok
        || bounds[0].iter().zip(a).any(|(x, y)| (x - y).abs() > slack)
        || bounds[1].iter().zip(b).any(|(x, y)| (x - y).abs() > slack)
    {
        return Err(failure(
            ErrorKind::Invalid,
            Some(artifact.producer),
            "Mesh readback differs from the generated objects, placements, volume or exact world bounds",
        ));
    }
    Ok(())
}
fn hash_file(path: &Path) -> Result<(ContentHash, u64), ExportError> {
    let metadata = fs::symlink_metadata(path).map_err(io)?;
    if !metadata.is_file() || metadata.file_type().is_symlink() || metadata.len() > MAX_STEP_BYTES {
        return Err(io("Output is not a regular bounded file"));
    }
    let bytes = fs::read(path).map_err(io)?;
    Ok((ContentHash::of(&bytes), bytes.len() as u64))
}
fn reuse(
    exporter: &Exporter,
    built: &Built,
    directory: &Path,
    actor: &Actor,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<RunReceipt, ExportError> {
    let metadata = fs::symlink_metadata(directory).map_err(io)?;
    if !metadata.is_dir() || metadata.file_type().is_symlink() {
        return Err(io("Existing export destination is not an owned directory"));
    }
    let manifest_path = directory.join(MANIFEST);
    let metadata = fs::symlink_metadata(&manifest_path).map_err(|_| {
        io("Existing export destination has no ownership manifest; choose another destination")
    })?;
    if !metadata.is_file() || metadata.file_type().is_symlink() || metadata.len() > 2 * 1024 * 1024
    {
        return Err(io("Invalid export manifest file"));
    }
    let mut manifest: Manifest =
        serde_json::from_slice(&fs::read(manifest_path).map_err(io)?).map_err(io)?;
    let preview = &built.preview;
    if manifest.schema != 1
        || manifest.provider != preview.provider
        || manifest.set != preview.key.set
        || manifest.source_hash != preview.key.source.source_hash
        || manifest.plan_hash != preview.key.plan_hash
        || manifest.files.len() != preview.artifacts.len()
    {
        return Err(io(
            "Existing export bundle belongs to different source or provider; nothing was overwritten",
        ));
    }
    let expected: BTreeSet<_> = preview
        .artifacts
        .iter()
        .map(|a| a.output.filename.as_str())
        .chain(std::iter::once(MANIFEST))
        .collect();
    let entries: Vec<_> = fs::read_dir(directory)
        .map_err(io)?
        .take(expected.len() + 1)
        .collect::<Result<_, _>>()
        .map_err(io)?;
    if entries.len() != expected.len() {
        return Err(io(
            "Export bundle contains additional or missing files; nothing was overwritten",
        ));
    }
    for entry in entries {
        let name = entry.file_name();
        let name = name
            .to_str()
            .ok_or_else(|| io("Unexpected filename in export bundle"))?;
        if !expected.contains(name) || !entry.file_type().map_err(io)?.is_file() {
            return Err(io(
                "Export bundle contains unowned files or symlinks; nothing was overwritten",
            ));
        }
    }
    for (stored, current) in manifest.files.iter_mut().zip(&preview.artifacts) {
        checkpoint(interrupt)?;
        if stored.artifact.output != current.output
            || stored.artifact.producer != current.producer
            || stored.artifact.solid != current.solid
        {
            return Err(io("Export bundle output bindings differ from current plan"));
        }
        let path = directory.join(&current.output.filename);
        let (hash, bytes) = hash_file(&path)?;
        if hash != stored.content_hash || bytes != stored.bytes {
            return Err(io("Export bundle was modified; nothing was overwritten"));
        }
        verify(exporter, built, current, &path)?;
        stored.artifact = current.clone();
    }
    Ok(RunReceipt {
        actor: actor.clone(),
        key: preview.key.clone(),
        directory: directory
            .to_str()
            .ok_or_else(|| io("Export path must be UTF-8"))?
            .into(),
        files: manifest.files,
        reused: true,
    })
}
pub(super) fn run(
    exporter: &Exporter,
    built: &Built,
    root: &Path,
    actor: &Actor,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<RunReceipt, ExportError> {
    checkpoint(interrupt)?;
    if !root.is_absolute() || root.to_str().is_none() || !fs::metadata(root).map_err(io)?.is_dir() {
        return Err(failure(
            ErrorKind::Invalid,
            None,
            "Choose an existing absolute UTF-8 export directory",
        ));
    }
    let preview = &built.preview;
    let prefix = filename_token(&preview.selection.name);
    if prefix.len() > 160 {
        return Err(failure(
            ErrorKind::Invalid,
            Some(preview.key.set),
            "Export-set name is too long for its bundle directory",
        ));
    }
    let name = format!(
        "{}_{}_{}_{}",
        prefix,
        preview.key.set,
        &preview.key.source.source_hash.to_string()[..12],
        &ContentHash::of(preview.provider.as_bytes()).to_string()[..8]
    );
    let directory = root.join(&name);
    let lock_path = root.join(format!(".multicad-{}.lock", preview.key.set));
    if fs::symlink_metadata(&lock_path).is_ok_and(|m| m.file_type().is_symlink()) {
        return Err(io("Export lock cannot be a symlink"));
    }
    // A persistent OS lock releases on crash; deleting lock files would allow
    // different inode locks to race across application instances.
    let lock = OpenOptions::new()
        .create(true)
        .truncate(false)
        .read(true)
        .write(true)
        .open(&lock_path)
        .map_err(io)?;
    lock.try_lock().map_err(|e| {
        io(format!(
            "Export destination is busy or cannot be locked: {e}"
        ))
    })?;
    if fs::symlink_metadata(&directory).is_ok() {
        return reuse(exporter, built, &directory, actor, interrupt);
    }
    let staging = tempfile::Builder::new()
        .prefix(".multicad-export-")
        .tempdir_in(root)
        .map_err(io)?;
    let model = built.model.as_ref().ok_or_else(|| {
        failure(
            ErrorKind::Invalid,
            None,
            "Export has no owned current model",
        )
    })?;
    let mut files = Vec::new();
    for artifact in &preview.artifacts {
        checkpoint(interrupt)?;
        let path = staging.path().join(&artifact.output.filename);
        if let Some(bytes) = built.generated.get(&artifact.output.filename) {
            let mut file = OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&path)
                .map_err(io)?;
            file.write_all(bytes).map_err(io)?;
        } else {
            artifact_solid(model, artifact)?
                .write_step(&path)
                .map_err(io)?;
        }
        checkpoint(interrupt)?;
        verify(exporter, built, artifact, &path)?;
        let (content_hash, bytes) = hash_file(&path)?;
        OpenOptions::new()
            .read(true)
            .write(true)
            .open(&path)
            .map_err(io)?
            .sync_all()
            .map_err(io)?;
        files.push(WrittenFile {
            artifact: artifact.clone(),
            content_hash,
            bytes,
        });
    }
    let manifest = Manifest {
        schema: 1,
        provider: preview.provider.clone(),
        set: preview.key.set,
        source_hash: preview.key.source.source_hash,
        plan_hash: preview.key.plan_hash,
        actor: actor.clone(),
        files: files.clone(),
    };
    let mut manifest_file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(staging.path().join(MANIFEST))
        .map_err(io)?;
    manifest_file
        .write_all(&serde_json::to_vec_pretty(&manifest).map_err(io)?)
        .map_err(io)?;
    manifest_file.sync_all().map_err(io)?;
    drop(manifest_file);
    checkpoint(interrupt)?;
    fs::rename(staging.path(), &directory).map_err(io)?;
    let _ = staging.keep();
    // Publication succeeded. A later cancellation cannot relabel a completed run.
    Ok(RunReceipt {
        actor: actor.clone(),
        key: preview.key.clone(),
        directory: directory
            .to_str()
            .ok_or_else(|| io("Export path must be UTF-8"))?
            .into(),
        files,
        reused: false,
    })
}
/// Read a written schedule CSV back: its kind and revision comments, a
/// header and one record per row, each with every column.
fn verify_csv(written: &[u8], evidence: &ScheduleEvidence) -> Result<(), String> {
    let text = std::str::from_utf8(written).map_err(|_| "Schedule CSV is not UTF-8")?;
    let (comments, records) = crate::schedules::read_csv(text)?;
    let kind = format!("MultiCAD {} schedule", evidence.spec.kind.token());
    let revision = format!("revision {}", evidence.revision);
    if comments.first() != Some(&kind)
        || comments.get(1) != Some(&revision)
        || records.len() != evidence.rows + 1
        || records.iter().any(|r| r.len() != evidence.columns)
    {
        return Err("Schedule CSV readback differs from the planned kind, revision or size".into());
    }
    Ok(())
}
