//! Per-member original-file verification for assemblies containing inserts.
use super::step::{fail, read_error};
use super::*;
use multicad_rust_kernel::{step_read, step_read::owned_insert as read};
const BYTES: usize = 64 * 1024 * 1024;
#[derive(Clone, Debug, PartialEq)]
pub struct PartitionedStepVerification {
    pub members: RetainedStepVerification,
    pub insert_correspondence: Vec<(u32, WoodInsertStepVerification)>,
    /// Authored boundary binding with SOURCE-integral mass, never decoded mass.
    pub source_known_members: Vec<(u32, EntryChamferBoundaryProof)>,
}
fn leaves<'a>(solid: &'a Solid, output: &mut Vec<&'a Solid>) {
    if let Some(compound) = solid.rust_compound() {
        for member in compound.members() {
            leaves(member, output);
        }
    } else {
        output.push(solid);
    }
}
impl Solid {
    pub fn has_wood_insert_step_source(&self) -> bool {
        if self.wood_insert_body().is_some() {
            return true;
        }
        self.rust_compound()
            .is_some_and(|c| c.members().any(Self::has_wood_insert_step_source))
    }
    /// Source families requiring explicit per-member mass provenance.
    pub fn has_partitioned_step_source(&self) -> bool {
        self.has_wood_insert_step_source()
            || (self.rust_compound().is_some() && contains_entry(self))
    }
    pub fn verify_wood_insert_step_compound(
        &self,
        text: &str,
        tolerance: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Option<Result<PartitionedStepVerification, GeometryError>> {
        if !self.has_wood_insert_step_source() {
            return None;
        }
        self.verify_partitioned_step_source(text, tolerance, interrupt)
    }
    /// Keep decoded scalar contracts distinct from independently integrated
    /// insert/radial intervals. Every represented occurrence must be matched.
    pub fn verify_partitioned_step_source(
        &self,
        text: &str,
        tolerance: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Option<Result<PartitionedStepVerification, GeometryError>> {
        if !self.has_partitioned_step_source() {
            return None;
        }
        Some((|| {
            if !tolerance.is_finite() || tolerance <= 0. {
                return Err(fail(GeometryErrorKind::Input, "Invalid STEP tolerance"));
            }
            let file =
                read::CompoundDocument::parse(text, read::DocumentOptions::default(), interrupt)
                    .map_err(read_error)?;
            let mut sources = Vec::new();
            leaves(self, &mut sources);
            let mut result = PartitionedStepVerification {
                members: RetainedStepVerification {
                    solid_count: 0,
                    bounds_mm: [[f64::INFINITY; 3], [f64::NEG_INFINITY; 3]],
                    retained_members: vec![],
                    decoded_scalar_members: vec![],
                },
                insert_correspondence: vec![],
                source_known_members: vec![],
            };
            let mut cursor = 0usize;
            for leaf in sources {
                if interrupt.is_some_and(|f| f()) {
                    return Err(fail(
                        GeometryErrorKind::Cancelled,
                        "STEP member verification superseded",
                    ));
                }
                let count = leaf.solid_count()? as usize;
                let end = cursor
                    .checked_add(count)
                    .ok_or_else(|| fail(GeometryErrorKind::Limit, "Member count overflow"))?;
                let actual = file.represented_members().get(cursor..end).ok_or_else(|| {
                    fail(
                        GeometryErrorKind::InvalidResult,
                        "Missing represented STEP member",
                    )
                })?;
                if let Some((body, frame)) = leaf.wood_insert_body() {
                    if count != 1
                        || (actual[0].mapping_link().is_none() && frame != PlaneFrame::default())
                    {
                        return Err(fail(
                            GeometryErrorKind::Unsupported,
                            "Placed insert verification awaits mapped emission",
                        ));
                    }
                    let policy = GeometryPolicy {
                        linear_tolerance_mm: leaf.policy.linear_tolerance_mm.min(tolerance),
                        ..leaf.policy
                    };
                    let proof = body.with_plan(policy, interrupt, |topology| {
                        let bytes = BYTES
                            .checked_sub(topology.retained_bytes().saturating_add(8192))
                            .filter(|&n| n > 0)
                            .ok_or_else(|| {
                                fail(
                                    GeometryErrorKind::Limit,
                                    "Insert member proof storage limit",
                                )
                            })?;
                        let packet = file
                            .member_packet(
                                actual[0].item_entity(),
                                read::Options {
                                    retained_budget_bytes: bytes,
                                    ..Default::default()
                                },
                                interrupt,
                            )
                            .map_err(read_error)?;
                        let parsed = packet.parsed();
                        let proof =
                            body.verify_parsed(topology, parsed, bytes, 8192, policy, interrupt)?;
                        if let Some(link) = actual[0].mapping_link() {
                            let mapping = parsed
                                .physical_mapping(
                                    link.operator_entity(),
                                    read::Options {
                                        retained_budget_bytes: bytes,
                                        ..Default::default()
                                    },
                                    interrupt,
                                )
                                .map_err(read_error)?;
                            super::step_pose::transport(
                                &body.source,
                                frame,
                                topology,
                                &mapping,
                                proof,
                                interrupt,
                            )
                        } else {
                            Ok(proof)
                        }
                    })?;
                    result
                        .members
                        .retained_members
                        .push(RetainedStepMemberMass {
                            member_index: cursor as u32,
                            decoded_volume_bounds_mm3: proof.decoded_volume_bounds_mm3,
                            source_volume_bounds_mm3: proof.source_volume_bounds_mm3,
                        });
                    include(&mut result.members.bounds_mm, proof.decoded_bounds_mm);
                    result.insert_correspondence.push((cursor as u32, proof));
                } else {
                    let mut decoded = Vec::with_capacity(count);
                    for (offset, member) in actual.iter().enumerate() {
                        let view = file
                            .local_member_document(
                                member.item_entity(),
                                read::Options::default(),
                                interrupt,
                            )
                            .map_err(read_error)?;
                        if let Some(proof) =
                            leaf.verify_retained_step_source(view.as_str(), tolerance)
                        {
                            if count != 1 {
                                return Err(fail(
                                    GeometryErrorKind::Unsupported,
                                    "Multi-member retained leaf not supported",
                                ));
                            }
                            let proof = proof.map_err(|e| {
                                fail(
                                    GeometryErrorKind::InvalidResult,
                                    format!("Retained member verification failed:{e:?}"),
                                )
                            })?;
                            for mut mass in proof.retained_members {
                                mass.member_index += cursor as u32;
                                result.members.retained_members.push(mass);
                            }
                            for (i, v) in proof.decoded_scalar_members {
                                result
                                    .members
                                    .decoded_scalar_members
                                    .push((i + cursor as u32, v));
                            }
                            include(&mut result.members.bounds_mm, proof.bounds_mm);
                            continue;
                        }
                        if let Some(proof) =
                            leaf.verify_entry_chamfer_step_with_tolerance(view.as_str(), tolerance)
                        {
                            let proof = proof?;
                            include(
                                &mut result.members.bounds_mm,
                                source_enclosure(leaf, proof.verified_tolerance_mm)?,
                            );
                            result
                                .source_known_members
                                .push(((cursor + offset) as u32, proof));
                            continue;
                        }
                        // Reuse the ordinary part's independent finite-file
                        // verifier, with the requested policy bound on a shared
                        // immutable source. No source geometry is regenerated.
                        let mut scoped = leaf.share()?;
                        scoped.policy.linear_tolerance_mm =
                            scoped.policy.linear_tolerance_mm.min(tolerance);
                        if let Some(volume) = scoped
                            .verify_helical_step(view.as_str())
                            .or_else(|| scoped.verify_normal_relief_step(view.as_str()))
                            .or_else(|| scoped.verify_herringbone_step(view.as_str()))
                            .or_else(|| scoped.verify_bevel_step(view.as_str()))
                            .or_else(|| scoped.verify_worm_step(view.as_str()))
                        {
                            let volume = volume.map_err(|e| {
                                fail(
                                    GeometryErrorKind::InvalidResult,
                                    format!("Authored member verification failed:{e:?}"),
                                )
                            })?;
                            include(
                                &mut result.members.bounds_mm,
                                source_enclosure(leaf, scoped.policy.linear_tolerance_mm)?,
                            );
                            result
                                .members
                                .decoded_scalar_members
                                .push(((cursor + offset) as u32, volume));
                            continue;
                        }
                        let members = step_read::read_analytic_compound_with_helical_cuts(
                            view.as_str(),
                            tolerance,
                        )
                        .map_err(|e| {
                            fail(
                                GeometryErrorKind::InvalidResult,
                                format!("Actual member decode failed:{e:?}"),
                            )
                        })?;
                        if members.len() != 1 {
                            return Err(fail(
                                GeometryErrorKind::InvalidResult,
                                "Member view must decode exactly one solid",
                            ));
                        }
                        for (readback, imported) in members {
                            include(&mut result.members.bounds_mm, readback.bounds_mm);
                            result
                                .members
                                .decoded_scalar_members
                                .push(((cursor + offset) as u32, readback.volume_mm3));
                            decoded.push((readback, imported));
                        }
                    }
                    if !decoded.is_empty() {
                        let topology = leaf.rust_topology().ok_or_else(|| {
                            fail(GeometryErrorKind::Unsupported, "Missing source topology")
                        })??;
                        if !leaf.step_readback_matches(&topology, &decoded, tolerance * 2.) {
                            return Err(fail(
                                GeometryErrorKind::InvalidResult,
                                "Actual member differs from source geometry",
                            ));
                        }
                    }
                }
                cursor = end;
            }
            if cursor != file.represented_members().len() {
                return Err(fail(
                    GeometryErrorKind::InvalidResult,
                    "Unexpected represented STEP member",
                ));
            }
            if interrupt.is_some_and(|f| f()) {
                return Err(fail(
                    GeometryErrorKind::Cancelled,
                    "STEP verification superseded",
                ));
            }
            result.members.solid_count = cursor as u32;
            Ok(result)
        })())
    }
}
fn include(bounds: &mut [[f64; 3]; 2], other: [[f64; 3]; 2]) {
    for k in 0..3 {
        bounds[0][k] = bounds[0][k].min(other[0][k]);
        bounds[1][k] = bounds[1][k].max(other[1][k]);
    }
}

// Finite boundary verifiers prove source correspondence. Expand the source
// envelope by that actual position budget; do not call it an exact file box.
fn source_enclosure(source: &Solid, tolerance: f64) -> Result<[[f64; 3]; 2], GeometryError> {
    let b = source.bounds()?;
    Ok([
        b.min.map(|x| (x - tolerance).next_down()),
        b.max.map(|x| (x + tolerance).next_up()),
    ])
}

fn contains_entry(solid: &Solid) -> bool {
    solid.rust_thread_entry_chamfer_data().is_some()
        || solid
            .rust_compound()
            .is_some_and(|c| c.members().any(contains_entry))
}
