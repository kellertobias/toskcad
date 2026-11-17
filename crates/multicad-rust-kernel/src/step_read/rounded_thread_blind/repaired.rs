//! A repaired interpretation, bound to the real member and retained records.
//! Canonical stock rationals are derived from decoded binary geometry. They
//! are not the original REAL values, family, fit, or an exact raw mass proof.
use super::*;
use crate::{Rational, rounded_thread::BlindBore};
use std::{collections::BTreeSet, sync::Arc};
const RECORD_BYTES: usize = 64 * 1024 * 1024;
type SourceRecords = Arc<[(u32, Vec<Box<str>>)]>;
#[derive(Clone, Debug, PartialEq)]
pub struct RepairedRoundedBlindCandidate {
    bore: BlindBore,
    frame: [Point; 3],
    source_member: u32,
    records: SourceRecords,
    c0: f64,
}
impl RepairedRoundedBlindCandidate {
    pub fn canonical_bore(&self) -> &BlindBore {
        &self.bore
    }
    pub fn frame(&self) -> [Point; 3] {
        self.frame
    }
    pub fn source_member_entity(&self) -> u32 {
        self.source_member
    }
    /// Original parsed component spellings and IDs of this member's complete
    /// referenced geometry graph. Representation membership/units were checked
    /// by the containing flat-document reader before this value was returned.
    pub fn source_records(&self) -> impl Iterator<Item = (u32, &[Box<str>])> {
        self.records
            .iter()
            .map(|(id, parts)| (*id, parts.as_slice()))
    }
    /// Earned complete-boundary positional repair allowance, matching each
    /// actual full-multiplicity spline cell to its normalized canonical cell.
    /// Positive raw knot spans may reparameterize that image. No original-chart
    /// derivative or independent exact-REAL mass claim is attached.
    pub fn maximum_source_c0_error_mm(&self) -> f64 {
        self.c0
    }
}
fn binary(value: f64) -> Result<Rational> {
    if !value.is_finite() || value <= 0. {
        return Err(ReadError::Invalid);
    }
    let bits = value.to_bits();
    let exp = ((bits >> 52) & 2047) as i32;
    if exp == 0 {
        return Err(ReadError::Limit);
    }
    let mut n = (bits & ((1u64 << 52) - 1)) | (1u64 << 52);
    let mut power = exp - 1023 - 52;
    let zeros = n.trailing_zeros();
    n >>= zeros;
    power += zeros as i32;
    let (n, d) = if power >= 0 {
        if power >= 63 || n > (i64::MAX as u64 >> power) {
            return Err(ReadError::Limit);
        }
        ((n << power) as i64, 1)
    } else {
        if -power >= 63 {
            return Err(ReadError::Limit);
        }
        (n as i64, 1i64 << (-power))
    };
    Rational::new(n, d).map_err(|_| ReadError::Limit)
}
fn raw_records(shell: &Shell<'_, '_>, member: u32) -> Result<SourceRecords> {
    // Before any graph/record copies, bound them by the actual parsed arena.
    // Components borrow nonoverlapping source slices. The fixed per-record
    // margin covers DFS Vec growth, visited BTree nodes and output metadata;
    // parser arena remains a separately bounded pre-existing phase.
    let mut total = 8192usize;
    total = total
        .checked_add(
            shell
                .entities
                .len()
                .checked_mul(384)
                .ok_or(ReadError::Limit)?,
        )
        .ok_or(ReadError::Limit)?;
    for parts in shell.entities.values() {
        for e in parts {
            total = total
                .checked_add(32 + e.source.ok_or(ReadError::Unsupported)?.len())
                .ok_or(ReadError::Limit)?;
        }
    }
    if total > RECORD_BYTES {
        return Err(ReadError::Limit);
    }
    let mut seen = BTreeSet::from([member]);
    let mut pending = vec![member];
    while let Some(id) = pending.pop() {
        let parts = shell.entities.get(&id).ok_or(ReadError::Invalid)?;
        fn refs_in(v: &Value<'_>, seen: &mut BTreeSet<u32>, pending: &mut Vec<u32>) {
            match v {
                Value::Reference(id) => {
                    if seen.insert(*id) {
                        pending.push(*id);
                    }
                }
                Value::List(values) => {
                    for v in values {
                        refs_in(v, seen, pending);
                    }
                }
                _ => {}
            }
        }
        for e in parts {
            for v in &e.args {
                refs_in(v, &mut seen, &mut pending);
            }
        }
    }
    let mut records = Vec::with_capacity(seen.len());
    for id in seen {
        let parts = shell.entities.get(&id).ok_or(ReadError::Invalid)?;
        records.push((
            id,
            parts
                .iter()
                .map(|e| e.source.map(Box::<str>::from).ok_or(ReadError::Unsupported))
                .collect::<Result<Vec<_>>>()?,
        ));
    }
    Ok(records.into())
}
pub fn read_repaired_rounded_blind_bore(
    input: &str,
    tolerance: f64,
) -> Result<(Readback, RepairedRoundedBlindCandidate)> {
    let mut members = compound::read_members_by(input, tolerance, true, true, |shell, id, _| {
        read_member(shell, id)
    })?;
    members.pop().ok_or(ReadError::Invalid)
}
fn read_member(
    shell: &mut Shell<'_, '_>,
    id: u32,
) -> Result<(Readback, RepairedRoundedBlindCandidate)> {
    raw_axis_class(shell)?;
    let faces = shell.member_faces(id)?;
    if faces.len() != 9 {
        return Err(ReadError::Unsupported);
    }
    let mut splines = Vec::with_capacity(5);
    let mut planes = Vec::with_capacity(3);
    let mut cylinders = Vec::with_capacity(1);
    let mut starts = Vec::with_capacity(5);
    for face_id in faces {
        let (surface, forward, _) = helical::face(shell, face_id)?;
        match shell.simple(surface)?.name {
            "B_SPLINE_SURFACE_WITH_KNOTS" => {
                let args = shell.args(surface, "B_SPLINE_SURFACE_WITH_KNOTS")?;
                let Value::List(rows) = args.get(3).ok_or(ReadError::Invalid)? else {
                    return Err(ReadError::Invalid);
                };
                let ids = refs(rows.first().ok_or(ReadError::Invalid)?)?;
                starts.push(shell.point(*ids.first().ok_or(ReadError::Invalid)?)?);
                splines.push(face_id);
            }
            "PLANE" => {
                helical::require(forward)?;
                let [Value::String(_), axis] = shell.args(surface, "PLANE")? else {
                    return Err(ReadError::Invalid);
                };
                planes.push((face_id, shell.placement(reference(axis)?)?));
            }
            "CYLINDRICAL_SURFACE" => {
                let [Value::String(_), _, Value::Number(r)] =
                    shell.args(surface, "CYLINDRICAL_SURFACE")?
                else {
                    return Err(ReadError::Invalid);
                };
                cylinders.push(*r);
            }
            _ => return Err(ReadError::Unsupported),
        }
    }
    helical::require(splines.len() == 5 && planes.len() == 3 && cylinders.len() == 1)?;
    let tolerance = shell.tolerance.min(1e-7);
    let error = tolerance / 16.;
    let entries = planes
        .iter()
        .enumerate()
        .filter(|(_, (_, p))| {
            starts
                .iter()
                .all(|&s| dot(sub(s, p.origin), p.normal).abs() <= error)
        })
        .map(|(i, _)| i)
        .collect::<Vec<_>>();
    let [entry] = entries[..] else {
        return Err(ReadError::Invalid);
    };
    let floors = planes
        .iter()
        .enumerate()
        .filter(|(i, (_, p))| *i != entry && dot(p.normal, planes[entry].1.normal) > 1. - 1e-10)
        .map(|(i, _)| i)
        .collect::<Vec<_>>();
    let [floor] = floors[..] else {
        return Err(ReadError::Invalid);
    };
    let top = (0..3)
        .find(|&i| i != entry && i != floor)
        .ok_or(ReadError::Invalid)?;
    splines.extend([planes[entry].0, planes[floor].0]);
    let thread = rounded_thread::import::candidate_with_tip_orientation(shell, &splines, true)?;
    let frame = [thread.origin_mm, thread.u, thread.v];
    let height = dot(
        sub(planes[top].1.origin, frame[0]),
        cross(frame[1], frame[2]),
    );
    let bore = BlindBore::new(thread.thread, binary(cylinders[0])?, binary(height)?)
        .map_err(|_| ReadError::Invalid)?;
    crate::step::rounded_thread_blind::precision(&bore, frame[0], tolerance)
        .map_err(|_| ReadError::Limit)?;
    let volume = rounded_thread::verify_repaired_blind_member(shell, id, &bore, frame, error)?;
    let mut bounds = [[0.; 3]; 2];
    let normal = cross(frame[1], frame[2]);
    for (end, sign) in [-1., 1.].into_iter().enumerate() {
        for axis in 0..3 {
            let support = bore
                .support_bounds(
                    [
                        sign * frame[1][axis],
                        sign * frame[2][axis],
                        sign * normal[axis],
                    ],
                    tolerance / 8.,
                )
                .map_err(|_| ReadError::Limit)?[1];
            let x = frame[0][axis] + sign * (support + tolerance);
            bounds[end][axis] = if end == 0 { x.next_down() } else { x.next_up() };
        }
    }
    let records = raw_records(shell, id)?;
    Ok((
        Readback {
            bounds_mm: bounds,
            volume_mm3: volume,
            face_count: 9,
            solid_count: 1,
            edges_mm: vec![],
            circles: vec![],
            semicircles: vec![],
            arcs: vec![],
            ellipses: vec![],
            splines: vec![],
        },
        RepairedRoundedBlindCandidate {
            bore,
            frame,
            source_member: id,
            records,
            c0: tolerance,
        },
    ))
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::rounded_thread::{Spec, Thread};
    fn source(left: bool) -> BlindBore {
        let p = 1.27;
        let a = 27.5_f64.to_radians();
        let r = p * a.cos() / (12. * (1. - a.sin()));
        BlindBore::new(
            Thread::new(Spec {
                major_radius_mm: 3.175,
                pitch_mm: p,
                length_mm: 3.81,
                depth_mm: p / (3. * a.tan()),
                crest_radius_mm: r,
                root_radius_mm: r,
                flank_angle_deg: 27.5,
                starts: 1,
                left_hand: left,
            })
            .unwrap(),
            Rational::new(127, 20).unwrap(),
            Rational::new(127, 25).unwrap(),
        )
        .unwrap()
    }
    #[test]
    fn actual_member_repair_binds_full_graph_and_retains_original_records_after_owner_drop() {
        for left in [false, true] {
            for posed in [false, true] {
                let source = source(left);
                let [origin, u, v] = if posed {
                    [[17., -8., 13.], [0.8, 0.6, 0.], [-0.36, 0.48, 0.8]]
                } else {
                    [[0.; 3], [1., 0., 0.], [0., 1., 0.]]
                };
                let text =
                    crate::step::write_rounded_thread_blind_bore(&source, origin, u, v, 1e-7)
                        .unwrap();
                let (readback, candidate) = read_repaired_rounded_blind_bore(&text, 1e-7).unwrap();
                assert_eq!(readback.face_count, 9);
                assert_eq!(candidate.canonical_bore().cutter().spec().left_hand, left);
                assert!(
                    (candidate.canonical_bore().cutter().spec().length_mm - 3.81).abs() < 1e-12
                );
                assert!(
                    (candidate.canonical_bore().stock_height_mm().as_f64() - 5.08).abs() < 1e-12
                );
                assert!(
                    candidate
                        .source_records()
                        .any(|(id, _)| id == candidate.source_member_entity())
                );
                assert_eq!(candidate.maximum_source_c0_error_mm(), 1e-7);
                let cylinder = candidate
                    .source_records()
                    .flat_map(|(_, p)| p)
                    .find(|p| p.starts_with("CYLINDRICAL_SURFACE"))
                    .unwrap()
                    .to_string();
                let (prefix, _) = cylinder.rsplit_once(',').unwrap();
                let changed = format!("{prefix},6.36)");
                let altered = text.replace(&cylinder, &changed);
                assert_ne!(altered, text);
                assert!(read_repaired_rounded_blind_bore(&altered, 1e-7).is_err());
                let (_, retry) = read_repaired_rounded_blind_bore(&text, 1e-7).unwrap();
                assert_eq!(retry, candidate);
                let bytes = candidate
                    .source_records()
                    .flat_map(|(_, p)| p)
                    .map(|p| p.len())
                    .sum::<usize>();
                drop(text);
                assert_eq!(
                    candidate
                        .source_records()
                        .flat_map(|(_, p)| p)
                        .map(|p| p.len())
                        .sum::<usize>(),
                    bytes
                );
                assert!(
                    candidate
                        .canonical_bore()
                        .support_bounds([1., 0., 0.], 1e-9)
                        .unwrap()[1]
                        > 6.35 - 1e-12
                );
            }
        }
    }
    #[test]
    fn decoded_binary_stock_rational_is_explicitly_canonical_and_range_checked() {
        for value in [1., 6.35, 5.08, 1e7] {
            assert_eq!(binary(value).unwrap().as_f64(), value);
        }
        for value in [0., -1., f64::NAN, f64::INFINITY, f64::MIN_POSITIVE, 1e100] {
            assert!(binary(value).is_err());
        }
    }
}
