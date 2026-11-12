//! Source-derived numerical seating for automatically generated hardware.
//!
//! Canonical drivers express the nominal engineering datum. The actual finite
//! head source supplies a bounded, outward-safe seating datum. The coupled
//! axial/roll adjustment preserves the helical phase and is applied before
//! publishing occurrence geometry, references, or display.
use super::*;
use multicad_kernel::{GeometryPolicy, rust_kernel::conical_screw_seat};

pub(crate) struct Inputs<'a> {
    pub document: &'a Document,
    pub parameters: &'a EvaluatedParameters,
    pub topologies: &'a BTreeMap<EntityId, TopologySnapshot>,
    pub source: &'a Solid,
    pub policy: GeometryPolicy,
    pub instance: EntityId,
    pub revision: u64,
}

pub(crate) fn refine(input: Inputs<'_>, nominal: RigidPose) -> Result<RigidPose, ModelError> {
    let Some((screw, source_frame)) = input.source.rust_profiled_screw_geometry() else {
        return Ok(nominal);
    };
    // Flat/crowned heads keep the existing planar seating path.
    if !screw
        .profile()
        .segments()
        .any(|s| s.role == multicad_kernel::rust_kernel::head_profile::Role::Cone)
    {
        return Ok(nominal);
    }
    let fail = |message: &str| error(Some(input.instance), message);
    let Some(EntityData::Instance {
        placement_source: Some(driver),
        ..
    }) = input.document.get(input.instance).map(|e| &e.data)
    else {
        return Ok(nominal);
    };
    let placement = model::read_screw_hardware_placement(input.document, *driver)
        .map_err(|e| error(Some(input.instance), e))?;
    let Some(binding) = placement.connection else {
        return Ok(nominal);
    };
    let request = model::read_screw_connection_request(input.document, binding.connection)
        .map_err(|e| error(Some(input.instance), e))?;
    let row = request
        .hardware
        .iter()
        .find(|row| row.role_id == binding.role)
        .ok_or_else(|| fail("The generated seating role is missing"))?;
    if row.role != Role::Screw {
        return Ok(nominal);
    }
    if placement.owner != input.instance
        || row.instance != input.instance
        || placement.axis != request.axis.id
        || !matches!(
            row.source,
            model::ScrewConnectionHardwareSource::Generated { .. }
        )
    {
        return Err(fail(
            "Physical seating requires the current generated screw owner and axis",
        ));
    }
    let [receiver] = request.receivers.as_slice() else {
        return Err(fail(
            "Physical conical seating requires one current receiver",
        ));
    };
    let resolved = receiver
        .source
        .resolve(input.parameters)
        .map_err(|e| error(Some(input.instance), e))?;
    let Some((true, diameter, depth)) = resolved.seat else {
        return Err(fail(
            "A retained countersunk head requires its actual conical receiving seat",
        ));
    };
    let seat = conical_screw_seat::Seat {
        bore_radius_mm: resolved.diameter_mm / 2.,
        mouth_radius_mm: diameter / 2.,
        depth_mm: depth,
    };
    let options = conical_screw_seat::Options {
        tolerance_mm: input.policy.linear_tolerance_mm,
        minimum_land_mm: input.policy.minimum_feature_mm,
    };
    let certified =
        conical_screw_seat::certify(screw.profile(), seat, options).map_err(|cause| {
            fail(&format!(
                "Actual conical seating cannot be certified: {cause:?}"
            ))
        })?;
    let axis = request
        .axis
        .source
        .resolve(
            input.document,
            input.parameters,
            &request.axis.face,
            input.revision,
            input.topologies,
        )
        .map_err(|e| error(Some(input.instance), e))?
        .frame(ReferenceFrame {
            origin_mm: [0.; 3],
            u: [1., 0., 0.],
            v: [0., 1., 0.],
        });
    let normal = axis.normal();
    let actual_normal: [f64; 3] =
        std::array::from_fn(|i| dot(nominal.rotation[i], source_frame.normal()));
    let tilt = cross(actual_normal, normal);
    let extent = screw.thread().spec().length_mm
        + screw.spec().shank_length_mm
        + screw.profile().top()[1].abs()
        + screw.profile().support([1., 0., 0.]);
    if dot(actual_normal, normal) <= 0.
        || dot(tilt, tilt).sqrt() * extent > input.policy.linear_tolerance_mm / 16.
    {
        return Err(fail(
            "The actual conical head axis differs from the receiving axis",
        ));
    }
    let local_neck = std::array::from_fn(|i| {
        (-screw.spec().shank_length_mm).mul_add(source_frame.normal()[i], source_frame.origin[i])
    });
    let actual_neck = nominal.point(local_neck);
    let (current, radial) = station(axis, actual_neck);
    let tolerance = input.policy.linear_tolerance_mm;
    let [lo, hi] = certified.ideal_neck_z_mm;
    if radial > tolerance / 8. || current < lo - tolerance || current > hi + tolerance {
        return Err(fail(
            "Actual head mounting differs from the canonical nominal seating datum",
        ));
    }
    let delta = certified.safe_neck_z_mm - current;
    if !delta.is_finite() || delta.abs() > tolerance {
        return Err(fail(
            "Physical seating correction exceeds the numerical placement budget",
        ));
    }
    // theta=(z-q)/rate: translating by delta and rotating by delta/rate
    // leaves the global thread phase unchanged.
    let roll = delta / screw.thread().helical_rate();
    let (sin, sin_half) = (roll.sin(), (roll / 2.).sin());
    let one_minus_cos = 2. * sin_half * sin_half;
    let rotate = |v: [f64; 3]| {
        let first = cross(normal, v);
        let second = cross(normal, first);
        std::array::from_fn(|i| one_minus_cos.mul_add(second[i], sin.mul_add(first[i], v[i])))
    };
    let columns: [[f64; 3]; 3] =
        std::array::from_fn(|i| rotate(std::array::from_fn(|j| nominal.rotation[j][i])));
    let rotation: [[f64; 3]; 3] = std::array::from_fn(|i| std::array::from_fn(|j| columns[j][i]));
    let target: [f64; 3] = std::array::from_fn(|i| {
        certified
            .safe_neck_z_mm
            .mul_add(normal[i], axis.origin_mm[i])
    });
    let translation_mm = std::array::from_fn(|i| {
        (-rotation[i][2]).mul_add(
            local_neck[2],
            (-rotation[i][1]).mul_add(
                local_neck[1],
                (-rotation[i][0]).mul_add(local_neck[0], target[i]),
            ),
        )
    });
    let result = RigidPose {
        rotation,
        translation_mm,
    };
    result
        .validate()
        .map_err(|e| error(Some(input.instance), e))?;
    // Validate the exact frame the kernel will publish, including its rounded
    // axes and origin. Distance to the proposal alone is not clearance proof.
    let placed = input
        .source
        .transformed(multicad_kernel::RigidTransform {
            rotation: result.rotation,
            translation_mm: result.translation_mm,
        })
        .map_err(|e| error(Some(input.instance), e))?;
    let (actual_source, actual_frame) = placed
        .rust_profiled_screw_geometry()
        .ok_or_else(|| fail("The placed screw lost its retained head source"))?;
    conical_screw_seat::certify_placement(
        actual_source.profile(),
        seat,
        options,
        conical_screw_seat::Placement {
            origin_mm: actual_frame.origin,
            axes: [actual_frame.u, actual_frame.v, actual_frame.normal()],
            axial_offset_mm: -actual_source.spec().shank_length_mm,
        },
        conical_screw_seat::Axis {
            origin_mm: axis.origin_mm,
            direction: normal,
        },
    )
    .map_err(|cause| {
        fail(&format!(
            "Rounded generated head placement cannot be certified: {cause:?}"
        ))
    })?;
    Ok(result)
}

fn cross([x, y, z]: [f64; 3], [u, v, w]: [f64; 3]) -> [f64; 3] {
    [y * w - z * v, z * u - x * w, x * v - y * u]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_screw_connection_fixture as fixture;
    use multicad_core::{
        commands::{Actor, Channel, Session},
        model::*,
        parameters,
    };
    use multicad_kernel::{
        HelicalThread, PlaneFrame, RigidTransform,
        rust_kernel::{head_profile::Head, profiled_screw::Spec},
    };

    #[test]
    fn generated_safe_seating_physically_moves_the_head_and_preserves_world_thread_phase() {
        for left in [false, true] {
            let (document, mut request, stock_id) = fixture::fixture();
            request.receivers[0].source.seat = Some(ConnectionReceiverSeat::Countersink {
                diameter: length(7.),
                depth: length(1.8),
            });
            for row in &mut request.hardware {
                let ScrewConnectionHardwareSource::Generated { part } = &mut row.source else {
                    unreachable!()
                };
                match &mut part.parameters {
                    MechanicalPartParameters::MachineScrew {
                        head,
                        drive,
                        thread,
                        ..
                    } => {
                        *head = ScrewHead::Countersunk {
                            profile: Box::new(CountersunkScrewHeadParameters {
                                diameter: length(6.),
                                height: length(1.65),
                                included_angle: angle(90.),
                                edge_chamfer: length(0.1),
                                neck_radius: length(0.1),
                            }),
                        };
                        *drive = ScrewDrive::HexSocket {
                            across_flats: length(2.5),
                            depth: length(0.8),
                        };
                        thread.length = length(6.);
                        thread.hand = if left {
                            ThreadHand::Left
                        } else {
                            ThreadHand::Right
                        };
                    }
                    MechanicalPartParameters::HexNut { thread, .. } => {
                        thread.hand = if left {
                            ThreadHand::Left
                        } else {
                            ThreadHand::Right
                        };
                    }
                    _ => unreachable!(),
                }
            }
            let session = Session::new(document).unwrap();
            let prepared =
                prepare_screw_connection(&session, 0, Actor::local(Channel::Ui), &request).unwrap();
            let document = prepared.prepared.document();
            let parameters = parameters::evaluate(document).unwrap();
            let stock = crate::Compiler {
                omit_display: true,
                ..Default::default()
            }
            .compile_sources_with_backend(
                document,
                0,
                &[stock_id],
                crate::GeometryBackend::Rust,
                None,
            )
            .unwrap();
            let topologies = stock
                .snapshot()
                .shapes
                .iter()
                .map(|(&id, s)| (id, s.topology.clone()))
                .collect();
            let row = &request.hardware[0];
            assert_eq!(row.role, ScrewConnectionHardwareRole::Screw);
            let source = Solid::rust_profiled_machine_screw(
                HelicalThread {
                    major_radius_mm: 1.5,
                    pitch_mm: 0.5,
                    length_mm: 6.,
                    depth_mm: 0.3,
                    crest_width_mm: 0.0625,
                    lower_flank_deg: 30.,
                    upper_flank_deg: 30.,
                    radial_taper: 0.,
                    starts: 1,
                    left_hand: left,
                },
                None,
                Spec {
                    head: Head::Countersunk {
                        radius_mm: 3.,
                        height_mm: 1.65,
                        included_angle_deg: 90.,
                        edge_chamfer_mm: 0.1,
                        neck_radius_mm: 0.1,
                    },
                    shank_length_mm: 2.,
                    socket_across_flats_mm: 2.5,
                    socket_depth_mm: 0.8,
                    minimum_wall_mm: 0.1,
                },
                EntityId::new(),
                PlaneFrame::default(),
                GeometryPolicy::default(),
            )
            .unwrap();
            let nominal = resolve_instance_pose_with_snapshots(
                document,
                row.instance,
                &parameters,
                0,
                &topologies,
            )
            .unwrap();
            let run = |source: &Solid| {
                refine(
                    Inputs {
                        document,
                        parameters: &parameters,
                        topologies: &topologies,
                        source,
                        policy: GeometryPolicy::default(),
                        instance: row.instance,
                        revision: 0,
                    },
                    nominal,
                )
            };
            let pose = run(&source).unwrap();
            assert_ne!(pose, nominal);
            let (screw, frame) = source.rust_profiled_screw_geometry().unwrap();
            let axis = request
                .axis
                .source
                .resolve(document, &parameters, &request.axis.face, 0, &topologies)
                .unwrap()
                .frame(ReferenceFrame {
                    origin_mm: [0.; 3],
                    u: [1., 0., 0.],
                    v: [0., 1., 0.],
                });
            let old = nominal.frame(ReferenceFrame {
                origin_mm: frame.origin,
                u: frame.u,
                v: frame.v,
            });
            let new = pose.frame(ReferenceFrame {
                origin_mm: frame.origin,
                u: frame.u,
                v: frame.v,
            });
            let phi = |f: ReferenceFrame| dot(f.u, axis.v).atan2(dot(f.u, axis.u));
            let delta = station(axis, new.origin_mm).0 - station(axis, old.origin_mm).0;
            assert!(delta < 0. && delta.abs() < GeometryPolicy::default().linear_tolerance_mm);
            let phase = phi(new) - phi(old) - delta / screw.thread().helical_rate();
            assert!(phase.abs() < 1e-12, "phase changed by {phase}");
            source
                .transformed(RigidTransform {
                    rotation: pose.rotation,
                    translation_mm: pose.translation_mm,
                })
                .unwrap();
            let shifted = source
                .transformed(RigidTransform {
                    translation_mm: [0., 0., 1.],
                    ..Default::default()
                })
                .unwrap();
            assert!(run(&shifted).is_err());
            assert_eq!(source.rust_profiled_screw_geometry().unwrap().1, frame);
        }
    }
}
