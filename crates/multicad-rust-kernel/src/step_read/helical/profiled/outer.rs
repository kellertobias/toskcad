//! Shared independently checked outer supports and original thread junction.
use super::*;
use crate::step::helical::HeadSource;
pub(super) struct Verified {
    pub frame: Frame,
    pub proof: Proof,
    pub top: Rim,
    pub mouth: u32,
    pub faces: Vec<u32>,
    pub segment_count: usize,
}
pub(super) fn verify(
    proof: Proof,
    shell: &mut Shell<'_, '_>,
    all_faces: &[u32],
    expected: HeadSource<'_>,
    rings: &[Vec<u32>; 2],
    corners: &[[[u32; 2]; 2]],
    placement: [Point; 3],
) -> Result<Verified> {
    let hex = expected
        .existing()
        .is_some_and(|s| s.recess().flat_hex_geometry().is_some());
    let phi = expected.phillips_geometry().is_some();
    closed_incidence(shell, all_faces)?;
    // The generic tensor proof has its own approximation budget. Reconstructing
    // an analytic head may use only its intrinsic source arithmetic reserve.
    let proof = Proof {
        error: if phi {
            Some(expected.profile().rounding_error_mm())
        } else {
            proof.error.map(|_| expected.profile().rounding_error_mm())
        },
        ..proof
    };
    let frame = Frame::new(placement)?;
    let thread = expected.thread();
    let count = 4 * thread.spec().starts as usize;
    let tail = &all_faces[count..];
    let join = junction(
        proof,
        shell,
        tail,
        thread,
        rings,
        corners,
        (frame, expected.visible_thread_interval_mm()[0]),
    )?;
    let faces = &tail[join.shoulders + 1..];
    let segments: Vec<_> = expected.profile().segments().collect();
    require(faces.len() == expected.face_count() - count - 2)?;
    if let Some(error) = proof.error {
        let analytic = if phi {
            &faces[1..segments.len() + 1]
        } else if hex {
            faces
        } else {
            &faces[..segments.len() + 1]
        };
        head_edge_residuals(shell, analytic, frame, error)?;
    }
    let shift = expected.plain_length_mm();
    let (surface, forward, bounds) = face(shell, faces[0])?;
    require(forward && bounds.len() == 2)?;
    let previous = single_bound_rim(shell, bounds[0], true, false, frame)?;
    if phi {
        let axis = hex_head::plane(shell, surface)?;
        let actual = shell.placement(reference(axis)?)?;
        proof.axis(shell, axis, previous.center, frame.n.map(|v| -v), actual.u)?;
    } else {
        plane_at(proof, shell, surface, previous, frame, false)?;
    }
    let top = previous;
    let socket_mouth = bounds[1];
    verify_walls(
        proof,
        shell,
        &faces[1..=segments.len()],
        MeridianWalls {
            segments: &segments,
            frame,
            shift,
            top: previous,
            join: &join,
        },
    )?;
    Ok(Verified {
        frame,
        proof,
        top,
        mouth: socket_mouth,
        faces: faces.to_vec(),
        segment_count: segments.len(),
    })
}

/// Shared actual outer wall/finite-curve proof. The caller supplies a retained
/// meridian, never a dummy recess or reconstructed authoring head recipe.
pub(super) struct MeridianWalls<'a> {
    pub segments: &'a [Segment],
    pub frame: Frame,
    pub shift: f64,
    pub top: Rim,
    pub join: &'a Junction,
}
pub(super) fn verify_walls(
    proof: Proof,
    shell: &mut Shell<'_, '_>,
    faces: &[u32],
    walls: MeridianWalls<'_>,
) -> Result<()> {
    let MeridianWalls {
        segments,
        frame,
        shift,
        top: mut previous,
        join,
    } = walls;
    require(faces.len() == segments.len())?;
    for (i, source) in segments.iter().enumerate() {
        require(
            previous.radius == source.start[0]
                && proof.equal(
                    previous.center,
                    frame.world([0.; 2], source.start[1] - shift),
                ),
        )?;
        let (surface, forward, bounds) = face(shell, faces[i])?;
        let last = i + 1 == segments.len();
        let end = if source.role == Role::Bearing {
            require(forward && bounds.len() == 2)?;
            require(
                ring(shell, bounds[0], true)?
                    == [(previous.edge.ok_or(ReadError::Invalid)?, previous.positive)],
            )?;
            plane_at(proof, shell, surface, previous, frame, true)?;
            if last {
                require(ring(shell, bounds[1], false)? == orient(&join.positive, false))?;
                join.end
            } else {
                single_bound_rim(shell, bounds[1], false, false, frame)?
            }
        } else {
            require(bounds.len() == 1)?;
            let (end, seam) = periodic_end(shell, bounds[0], previous, last, join, frame)?;
            match source.curve {
                HeadCurve::Line => line_surface(
                    proof,
                    shell,
                    (surface, forward, seam),
                    [previous, end],
                    source,
                    frame,
                )?,
                HeadCurve::Arc {
                    center,
                    radius_mm,
                    sweep,
                } => {
                    let arc = ArcSource {
                        center: [center[0], center[1] - shift],
                        radius: radius_mm,
                        angle: (source.start[1] - center[1])
                            .portable_atan2(source.start[0] - center[0]),
                        sweep,
                        role: source.role,
                    };
                    curved_surface(
                        proof,
                        shell,
                        (surface, forward, seam),
                        [previous, end],
                        arc,
                        frame,
                    )?;
                }
            }
            end
        };
        require(
            end.radius == source.end[0]
                && proof.equal(end.center, frame.world([0.; 2], source.end[1] - shift)),
        )?;
        previous = end;
    }
    Ok(())
}
