//! Preserve the original pipe law and canonical neck/mouth placement.
use super::*;
use multicad_core::model::ThreadFamily;
use multicad_kernel::PipeHardwareKind;

pub(super) fn build(
    feature: EntityId,
    resolved: &ResolvedMechanicalPart,
    policy: GeometryPolicy,
    factory: Option<&mut multicad_kernel::PipeThreadFactory>,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Option<Result<Solid, ModelError>> {
    let thread = match resolved {
        ResolvedMechanicalPart::MachineScrew { thread, .. }
        | ResolvedMechanicalPart::HexNut { thread, .. } => thread,
        _ => return None,
    };
    let family = match (&thread.npt, &thread.pipe) {
        (Some(_), None) => ThreadFamily::Npt,
        (None, Some(_)) => ThreadFamily::Bspt,
        (None, None) => return None,
        (Some(_), Some(_)) => {
            return Some(Err(error(
                Some(feature),
                "Pipe hardware has conflicting physical thread sources",
            )));
        }
    };
    Some((|| {
        if thread.representation != ThreadRepresentation::Modeled {
            return Err(error(
                Some(feature),
                "Pipe hardware requires its modeled original thread source",
            ));
        }
        let (kind, frame) = match resolved {
            ResolvedMechanicalPart::MachineScrew {
                underhead_length_mm,
                head:
                    ResolvedScrewHead::ExternalHex {
                        across_flats_mm,
                        height_mm,
                    },
                drive: ResolvedScrewDrive::ExternalHex,
                ..
            } => {
                if *underhead_length_mm != thread.length_mm {
                    return Err(error(
                        Some(feature),
                        "Pipe hex screws currently require a fully threaded underhead length; a plain shank needs its own physical union",
                    ));
                }
                // Native retains the original neck frame and exact coupled L
                // in its proper tip-to-neck construction. Do not round a new
                // tip datum here or reflect the original thread hand.
                let frame = PlaneFrame {
                    origin: thread.frame.origin,
                    u: thread.frame.u,
                    v: thread.frame.v,
                };
                (
                    PipeHardwareKind::HexHead {
                        across_flats_mm: *across_flats_mm,
                        height_mm: *height_mm,
                    },
                    frame,
                )
            }
            ResolvedMechanicalPart::HexNut {
                across_flats_mm, ..
            } => (
                PipeHardwareKind::HexNut {
                    across_flats_mm: *across_flats_mm,
                },
                PlaneFrame {
                    origin: thread.frame.origin,
                    u: thread.frame.u,
                    v: thread.frame.v,
                },
            ),
            _ => {
                return Err(error(
                    Some(feature),
                    "Pipe screws currently require a physical external hex head and drive",
                ));
            }
        };
        let spec = crate::threads::pipe::authored_spec(feature, family, thread)?;
        match factory {
            Some(factory) => factory.build(
                spec,
                kind,
                thread.minimum_wall_mm,
                feature,
                frame,
                policy,
                interrupt,
            ),
            None => Solid::rust_authored_axial_pipe_hardware_with_interrupt(
                spec,
                kind,
                thread.minimum_wall_mm,
                feature,
                frame,
                policy,
                interrupt,
            ),
        }
        .map_err(|cause| error(Some(feature), cause))
    })())
}
