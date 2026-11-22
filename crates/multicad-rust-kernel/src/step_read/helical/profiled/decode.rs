//! Retain independently decoded finite carriers; never infer a head recipe.
use super::*;
use crate::profiled_screw::SourceSpec;
use crate::{head_profile::TransportError, tapered_thread::interval::I};
mod header;

/// Retained meridian/header extraction only. Literal tensors and all finite
/// trims stay in the caller's immutable file arena for independent matching.
pub(in super::super::super) struct PhillipsOuterHeader {
    pub profile: Profile,
    pub head_top_face: u32,
    pub head_top_surface: u32,
    pub outer_bound: u32,
    pub removed_hole: u32,
    pub profile_faces: Vec<u32>,
}
pub(in super::super::super) fn read_phillips_outer_header(
    shell: &mut Shell<'_, '_>,
    faces: &[u32],
    thread: &Thread,
    placement: [Point; 3],
    visible_start: f64,
    cavity: &crate::step_read::ParsedPhillipsBoundary,
) -> Result<PhillipsOuterHeader> {
    let layout = header::layout(shell, faces, thread, placement, visible_start)?;
    let cavity_faces: BTreeSet<_> = cavity
        .faces()
        .iter()
        .filter_map(|f| f.material_entity)
        .collect();
    let count = layout
        .head
        .len()
        .checked_sub(cavity_faces.len() + 1)
        .ok_or(ReadError::Unsupported)?;
    require((1..=5).contains(&count))?;
    require(
        layout.head[0] == cavity.mouth().material_face
            && layout.head[count + 1..]
                .iter()
                .copied()
                .collect::<BTreeSet<_>>()
                == cavity_faces,
    )?;
    let decoded = header::meridian(shell, &layout, count)?;
    // This helper's sole caller constructs cavity from this same Shell before
    // extracting this header. The verified head-face identity above therefore
    // means meridian reads exactly the same raw surface and bounds[1] that own
    // Mouth. Shell caches do not mutate those immutable entity fields; repeating
    // their equality as a fallible cross-arena check would be redundant.
    Ok(PhillipsOuterHeader {
        profile: decoded.profile,
        head_top_face: layout.head[0],
        head_top_surface: decoded.top_surface,
        outer_bound: decoded.outer_bound,
        removed_hole: decoded.mouth,
        profile_faces: layout.head[1..=count].to_vec(),
    })
}

// Decimal STEP numbers are rounded to f64 once by the parser. One adjacent
// representable value on each side encloses that numeric conversion; this is
// independent of the caller's geometry tolerance. Directed intervals then
// enclose every subtraction, cross product, determinant and inverse projection.
#[derive(Clone, Copy)]
pub(super) struct ChartEnclosure {
    origin: [I; 3],
    rows: [[I; 3]; 3],
    determinant: I,
}
pub(super) fn icross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
pub(super) fn idot(a: [I; 3], b: [I; 3]) -> I {
    a.into_iter()
        .zip(b)
        .fold(I::point(0.), |sum, (x, y)| sum + x * y)
}
impl ChartEnclosure {
    fn new(frame: Frame, conversion: bool) -> Result<Self> {
        let enclose = |v: f64| {
            if conversion {
                I::new(v.next_down(), v.next_up())
            } else {
                I::point(v)
            }
        };
        let [u, v] = [frame.u.map(enclose), frame.v.map(enclose)];
        let n = icross(u, v);
        Self::from_axes(frame.origin.map(enclose), [u, v, n])
    }
    pub(super) fn from_axes(origin: [I; 3], [u, v, n]: [[I; 3]; 3]) -> Result<Self> {
        let rows = [icross(v, n), icross(n, u), icross(u, v)];
        let determinant = idot(u, rows[0]);
        require(determinant.valid() && determinant.lo > 0.)?;
        Ok(Self {
            origin,
            rows,
            determinant,
        })
    }
    pub(super) fn vector(&self, vector: [I; 3]) -> Result<[I; 3]> {
        let local = self.rows.map(|row| idot(row, vector) / self.determinant);
        require(local.into_iter().all(I::valid))?;
        Ok(local)
    }
    pub(super) fn local(&self, point: Point, conversion: bool) -> Result<[I; 3]> {
        let point = point.map(|v| {
            if conversion {
                I::new(v.next_down(), v.next_up())
            } else {
                I::point(v)
            }
        });
        self.point(point)
    }
    pub(super) fn point(&self, point: [I; 3]) -> Result<[I; 3]> {
        self.vector(std::array::from_fn(|i| point[i] - self.origin[i]))
    }
}
pub(super) fn references(value: &Value<'_>, out: &mut Vec<u32>) {
    match value {
        Value::Reference(id) => out.push(*id),
        Value::List(values) => {
            for value in values {
                references(value, out);
            }
        }
        _ => {}
    }
}
fn transport(shell: &Shell<'_, '_>, head: &[u32], frame: Frame) -> Result<TransportError> {
    let literal = ChartEnclosure::new(frame, false)?;
    let numeric = ChartEnclosure::new(frame, true)?;
    let mut pending = head.to_vec();
    let mut visited = BTreeSet::new();
    let mut error = TransportError::ZERO;
    while let Some(id) = pending.pop() {
        if !visited.insert(id) {
            continue;
        }
        for entity in shell.components(id)? {
            if entity.name == "CARTESIAN_POINT"
                && matches!(entity.args.as_slice(), [Value::String(_),Value::List(values)] if values.len()==3)
            {
                let point = shell.point(id)?;
                let actual = frame.local(point);
                let distance = |bounds: [I; 3]| {
                    actual.into_iter().zip(bounds).fold(0_f64, |sum, (a, i)| {
                        (sum + (a - i.lo).abs().max((a - i.hi).abs())).next_up()
                    })
                };
                error.arithmetic_mm = error
                    .arithmetic_mm
                    .max(distance(literal.local(point, false)?));
                error.numeric_conversion_mm = error
                    .numeric_conversion_mm
                    .max(distance(numeric.local(point, true)?));
            }
            for value in &entity.args {
                references(value, &mut pending);
            }
        }
    }
    require(error.total_mm().is_finite())?;
    Ok(error)
}

pub(super) fn volume_scale_error(frame: Frame) -> Result<f64> {
    let matrix = ChartEnclosure::new(frame, true)?;
    let error = (matrix.determinant - I::point(1.)).magnitude();
    require(error.is_finite())?;
    Ok(error)
}

pub(in super::super::super) fn read(
    shell: &mut Shell<'_, '_>,
    faces: &[u32],
    thread: &Thread,
    placement: [Point; 3],
    visible_start: f64,
) -> Result<Screw> {
    let layout = header::layout(shell, faces, thread, placement, visible_start)?;
    let frame = layout.frame;
    let count = layout.count;
    let head = &layout.head;
    let torx = (15..=19).contains(&head.len());
    require(torx || (9..=13).contains(&head.len()))?;
    let segment_count = head.len() - if torx { 14 } else { 8 };
    let decoded = header::meridian(shell, &layout, segment_count)?;
    let profile = decoded.profile;
    let top = decoded.top;
    let mouth = decoded.mouth;
    let reserve = profile.rounding_error_mm();
    if torx {
        let recess = super::torx::read(shell, faces, head, segment_count, frame, count, top)?;
        let error = recess.source_error_mm();
        let radius = (recess.maximum_radius_bound_mm() + error).next_up();
        let depth = (recess.maximum_depth_mm() + error).next_up();
        let wall = certified_wall(thread, &profile, radius, depth, reserve + error)?;
        return Screw::from_source_with_union(
            *thread,
            profile,
            SourceSpec {
                head: None,
                shank_length_mm: 0.,
                minimum_wall_mm: wall,
            },
            recess,
        )
        .map_err(|_| ReadError::Invalid);
    }
    let edges = ring(shell, mouth, false)?;
    require(edges.len() == 6)?;
    let mut corners = [[0.; 2]; 6];
    for (i, (edge, forward)) in edges.into_iter().enumerate() {
        let vertices = shell.edge(edge)?;
        let point = frame.local(shell.vertex(vertices[usize::from(!forward)])?);
        require((point[2] - top.z).abs() <= reserve)?;
        corners[i] = [point[0], point[1]];
    }
    let across = corners
        .into_iter()
        .map(|p| p[0].portable_hypot(p[1]))
        .fold(0., f64::max)
        * 3_f64.sqrt();
    let drive = crate::hex_drive::HexDrive::from_retained_corners(corners, across, reserve)
        .map_err(|_| ReadError::Invalid)?;
    let (surface, forward, bounds) = face(shell, head[segment_count + 1])?;
    require(forward && bounds.len() == 1)?;
    let [Value::String(_), axis] = shell.args(surface, "PLANE")? else {
        return Err(ReadError::Invalid);
    };
    let floor = shell.placement(reference(axis)?)?;
    let depth = frame.local(floor.origin)[2] - top.z;
    let wall = certified_wall(thread, &profile, drive.maximum_radius_mm(), depth, reserve)?;
    Screw::from_source_with_union(
        *thread,
        profile,
        SourceSpec {
            head: None,
            shank_length_mm: 0.,
            minimum_wall_mm: wall,
        },
        crate::screw_recess::Recess::flat_hex(drive, depth).map_err(|_| ReadError::Invalid)?,
    )
    .map_err(|_| ReadError::Invalid)
}

/// Certify actual surrounding radial material over the complete cavity depth.
/// The reserve starts the proof; it is not a dimensional feature wall. Callers
/// expand the recess by its own spatial error before using this head contract.
fn certified_wall(
    thread: &Thread,
    profile: &Profile,
    radius: f64,
    depth: f64,
    reserve: f64,
) -> Result<f64> {
    let back = (I::point(profile.top()[1]) + I::point(depth)).hi;
    let admits = |wall| Screw::admits_envelope(thread, profile, 0., radius, depth, wall, back);
    let mut wall = 4. * reserve;
    require(admits(wall))?;
    let mut upper = profile.support([1., 0., 0.]);
    for _ in 0..52 {
        let candidate = wall + (upper - wall) * 0.5;
        if admits(candidate) {
            wall = candidate;
        } else {
            upper = candidate;
        }
    }
    Ok(wall)
}

/// Independently match the retained meridian against its actual finite outer
/// supports. The supplied Profile was extracted once from this same arena; no
/// authored Head or recess recipe is reconstructed or supplied here.
pub(in crate::step_read) fn match_phillips_outer(
    shell: &mut Shell<'_, '_>,
    faces: &[u32],
    thread: &Thread,
    (placement, visible): ([Point; 3], f64),
    profile: &Profile,
    prefix: &super::super::Prefix,
) -> Result<()> {
    let frame = Frame::new(placement)?;
    let proof = Proof {
        error: Some(profile.rounding_error_mm()),
        extent: thread.spec().length_mm.max(thread.spec().major_radius_mm),
    };
    let count = 4 * thread.spec().starts as usize;
    let join = junction(
        proof,
        shell,
        &faces[count..],
        thread,
        &prefix.rings,
        &prefix.corners,
        (frame, visible),
    )?;
    let head = &faces[count + join.shoulders + 1..];
    let segments: Vec<_> = profile.segments().collect();
    require(head.len() == segments.len() + 29)?;
    head_edge_residuals(
        shell,
        &head[1..=segments.len()],
        frame,
        profile.rounding_error_mm(),
    )?;
    let (surface, forward, bounds) = face(shell, head[0])?;
    require(forward && bounds.len() == 2)?;
    let top = single_bound_rim(shell, bounds[0], true, false, frame)?;
    let axis = hex_head::plane(shell, surface)?;
    let actual = shell.placement(reference(axis)?)?;
    proof.axis(shell, axis, top.center, frame.n.map(|v| -v), actual.u)?;
    outer::verify_walls(
        proof,
        shell,
        &head[1..=segments.len()],
        outer::MeridianWalls {
            segments: &segments,
            frame,
            shift: 0.,
            top,
            join: &join,
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn physical_wall_certificate_has_real_radial_thickness_and_rejects_breakthrough() {
        let profile = Profile::new(
            crate::head_profile::Head::Crowned {
                radius_mm: 2.8,
                height_mm: 2.4,
                crown_radius_mm: 5.,
                top_land_radius_mm: 2.05,
                neck_radius_mm: 0.125,
            },
            1.5,
            1.,
        )
        .unwrap();
        let thread = Thread::new(crate::helical_thread::Spec {
            major_radius_mm: 1.5,
            pitch_mm: 0.5,
            length_mm: 7.5,
            depth_mm: 0.2,
            crest_width_mm: 0.0625,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            starts: 1,
            left_hand: false,
            root_radius_mm: None,
        })
        .unwrap();
        let reserve = profile.rounding_error_mm();
        let wall = certified_wall(&thread, &profile, 1.4, 1.9, reserve).unwrap();
        // The mouth land is the limiting radius throughout this cavity:
        // 2.05-1.4=0.65mm, independent of the precision reserve.
        assert!(wall > 0.64 && wall <= 0.65);
        assert!(profile.admits_socket(1.4, 1.9, wall));
        assert!(!profile.admits_socket(1.4, 1.9, wall + 1e-5));
        assert_eq!(
            certified_wall(&thread, &profile, 2.1, 1.9, reserve),
            Err(ReadError::Invalid)
        );
    }
    #[test]
    fn transport_encloses_data_and_arithmetic_but_refuses_unresolved_inverse() {
        let frame =
            Frame::new([[100000., -100000., 100000.], [0.6, 0.8, 0.], [0., 0., 1.]]).unwrap();
        let point = [100001.2, -99998.4, 99999.5];
        let literal = ChartEnclosure::new(frame, false)
            .unwrap()
            .local(point, false)
            .unwrap();
        let numeric = ChartEnclosure::new(frame, true)
            .unwrap()
            .local(point, true)
            .unwrap();
        for ((nominal, literal), numeric) in
            frame.local(point).into_iter().zip(literal).zip(numeric)
        {
            assert!(literal.lo <= nominal && nominal <= literal.hi);
            assert!(numeric.lo <= literal.lo && numeric.hi >= literal.hi);
        }
        assert!(
            ChartEnclosure::new(frame, false)
                .unwrap()
                .local([f64::INFINITY, 0., 0.], false)
                .is_err()
        );
        let ill = Frame::new([[0.; 3], [1., 1., 0.], [1., 1. + 1e-15, 0.]]).unwrap();
        assert!(ChartEnclosure::new(ill, true).is_err());
        assert!(Frame::new([[0.; 3], [0.; 3], [0.; 3]]).is_err());
    }
}
