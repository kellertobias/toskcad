//! Atomic authored export with independently read, support-derived healed
//! material/mass. Qualification never asserts the raw trim occupied mass.
use super::*;
mod replay;
use crate::{
    authored_hanger_bolt::{BoundaryOptions, Hanger},
    step_read::owned_hanger as read,
};
pub use replay::QualifiedAppendReceipt;
#[derive(Clone, Copy)]
pub struct QualifiedExportOptions {
    pub maximum_work: usize,
    pub maximum_bytes: usize,
    pub maximum_healed_volume_width_mm3: f64,
}
impl Default for QualifiedExportOptions {
    fn default() -> Self {
        Self {
            maximum_work: 640_000_000,
            maximum_bytes: 64 * 1024 * 1024,
            maximum_healed_volume_width_mm3: 1e-6,
        }
    }
}
pub struct QualifiedHangerStep {
    source: Hanger,
    candidate: StepCandidate,
    work: usize,
    bytes: usize,
    volume: [f64; 2],
    sewing: f64,
    source_position: f64,
    support_derivatives: [f64; 2],
    original_derivative: f64,
}
impl QualifiedHangerStep {
    pub fn source(&self) -> &Hanger {
        &self.source
    }
    pub fn as_str(&self) -> &str {
        self.candidate.as_str()
    }
    pub fn into_string(self) -> String {
        self.candidate.into_string()
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    /// Conservative phase-wide accounted upper bound; not process RSS.
    pub fn peak_accounted_bytes(&self) -> usize {
        self.bytes
    }
    pub fn healed_volume_bounds_mm3(&self) -> [f64; 2] {
        self.volume
    }
    pub fn sewing_correspondence_bound_mm(&self) -> f64 {
        self.sewing
    }
    pub fn source_to_raw_position_error_mm(&self) -> f64 {
        self.source_position
    }
    pub fn source_support_derivative_error_mm(&self) -> [f64; 2] {
        self.support_derivatives
    }
    /// Original-source to transmitted raw coedge D1; not healed-cap D1.
    pub fn original_coedge_derivative_error_mm(&self) -> f64 {
        self.original_derivative
    }
    /// Composed C0 envelope, not a healed-cap derivative certificate.
    pub fn source_to_healed_position_error_mm(&self) -> f64 {
        (I::point(self.source_position) + I::point(self.sewing)).hi
    }
}
struct Ledger<'a> {
    used: usize,
    left: usize,
    bytes: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Ledger<'_> {
    fn charge(&mut self, n: usize) -> Result<(), Error> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        self.left = self.left.checked_sub(n).ok_or(Error::Limit)?;
        self.used = self.used.checked_add(n).ok_or(Error::Limit)?;
        Ok(())
    }
    fn allowance(&self, maximum: usize, prefix: usize) -> Result<usize, Error> {
        if self.left == 0 {
            return Err(Error::Limit);
        }
        Ok(maximum.min(prefix.saturating_add(self.left)))
    }
    fn bytes(&self, extra: usize) -> Result<usize, Error> {
        self.bytes
            .checked_sub(extra)
            .filter(|v| *v > 0)
            .ok_or(Error::Limit)
    }
    fn curve(&self, prefix: usize, extra: usize) -> Result<CurveOptions, Error> {
        Ok(CurveOptions {
            work: self.allowance(16_000_000, prefix)?,
            bytes: self.bytes(extra)?,
            ..CurveOptions::default()
        })
    }
    fn domain(&self, prefix: usize, extra: usize) -> Result<DomainOptions, Error> {
        Ok(DomainOptions {
            maximum_work: self.allowance(16_000_000, prefix)?,
            maximum_bytes: self.bytes(extra)?,
        })
    }
}
fn source_error(e: crate::authored_hanger_bolt::Error) -> Error {
    match e {
        crate::authored_hanger_bolt::Error::Invalid => Error::Invalid,
        crate::authored_hanger_bolt::Error::Limit => Error::Limit,
        crate::authored_hanger_bolt::Error::Cancelled => Error::Cancelled,
        _ => Error::Unresolved,
    }
}
fn read_error(e: read::Error) -> Error {
    match e {
        read::Error::Invalid => Error::Invalid,
        read::Error::Limit => Error::Limit,
        read::Error::Cancelled => Error::Cancelled,
        read::Error::Unsupported => Error::Unresolved,
    }
}
fn produce(
    source: &Hanger,
    l: &mut Ledger<'_>,
) -> Result<(StepCandidate, f64, [f64; 2], f64), Error> {
    let i = l.interrupt;
    let boundary = source
        .boundary(
            BoundaryOptions {
                maximum_work: l.allowance(16_000_000, 0)?,
                maximum_bytes: l.bytes.min(1 << 20),
            },
            i,
        )
        .map_err(source_error)?;
    l.charge(boundary.work_used())?;
    let original = boundary.retained_bytes();
    let field = EncodedField::prepare(&boundary, l.curve(0, original)?, i)?;
    l.charge(field.work_used())?;
    let domains = FiniteDomains::prepare(&field, l.domain(field.work_used(), original)?, i)?;
    l.charge(domains.work_used() - field.work_used())?;
    let closed = ClosedField::prepare(&domains, l.domain(domains.work_used(), original)?, i)?;
    l.charge(closed.work_used() - domains.work_used())?;
    let surfaces = CompactSurfaces::prepare(
        &field,
        l.domain(
            field.work_used(),
            original + closed.retained_bytes() - field.retained_bytes(),
        )?,
        i,
    )?;
    l.charge(surfaces.work_used() - field.work_used())?;
    let compact_prefix = closed.work_used() + surfaces.work_used() - field.work_used();
    let compact =
        CompactBoundary::prepare(&closed, &surfaces, l.domain(compact_prefix, original)?, i)?;
    l.charge(compact.work_used() - compact_prefix)?;
    let extra = original + compact.retained_bytes() + surfaces.retained_bytes();
    let vertices = closed.literal_vertices(l.curve(closed.work_used(), extra)?, i)?;
    l.charge(vertices.work_used() - closed.work_used())?;
    let events = vertices.cap_events(l.domain(vertices.work_used(), extra)?, i)?;
    l.charge(events.work_used() - vertices.work_used())?;
    let topo_prefix = compact.work_used() + events.work_used() - closed.work_used();
    let topology = compact.trim_caps(&events, l.domain(topo_prefix, original)?, i)?;
    l.charge(topology.work_used() - topo_prefix)?;
    // Caller-owned collections remain live beside the helper ancestors. Their
    // full capacity is removed from every subsequent helper's byte allowance.
    let slots = events.segments().len()
        * (size_of::<BoundAffineCap<'_, '_, '_, '_, '_, '_>>()
            + size_of::<BoundRootCap<'_, '_, '_, '_, '_, '_>>()
            + size_of::<&dyn CapLiteral>())
        + compact.edges().len()
            * (size_of::<CompactHelix<'_, '_, '_, '_, '_, '_>>()
                + size_of::<&CompactHelix<'_, '_, '_, '_, '_, '_>>())
        + compact.faces().len()
            * (size_of::<PlainWall<'_, '_, '_, '_, '_>>()
                + size_of::<&PlainWall<'_, '_, '_, '_, '_>>());
    let mut held = original + slots + topology.retained_bytes();
    let mut affine = Vec::with_capacity(events.segments().len());
    let mut roots = Vec::with_capacity(events.segments().len());
    for segment in 0..events.segments().len() {
        l.charge(1)?;
        let options = l.curve(
            events.work_used(),
            held + compact.retained_bytes() + surfaces.retained_bytes(),
        )?;
        if let Some(cap) = events.bound_affine_cap(segment, options, i)? {
            l.charge(cap.work_used() - events.work_used())?;
            held += cap.retained_bytes();
            affine.push(cap)
        } else {
            let cap = events
                .bound_root_cap(segment, options, i)?
                .ok_or(Error::Unresolved)?;
            l.charge(cap.work_used() - events.work_used())?;
            held += cap.retained_bytes();
            roots.push(cap)
        }
    }
    let mut helices = Vec::with_capacity(compact.edges().len());
    let helix_prefix = compact.work_used() + vertices.work_used() - closed.work_used();
    for edge in 0..compact.edges().len() {
        l.charge(1)?;
        let extra = held + events.retained_bytes() - vertices.retained_bytes();
        if let Some(helix) = compact.helix(edge, &vertices, l.domain(helix_prefix, extra)?, i)? {
            l.charge(helix.work_used() - helix_prefix)?;
            held += helix.retained_bytes();
            helices.push(helix)
        }
    }
    let mut plains = Vec::with_capacity(compact.faces().len());
    for face in 0..compact.faces().len() {
        l.charge(1)?;
        let extra = held + events.retained_bytes() - closed.retained_bytes();
        if let Some(plain) = compact.plain_wall(face, l.domain(compact.work_used(), extra)?, i)? {
            l.charge(plain.work_used() - compact.work_used())?;
            held += plain.retained_bytes();
            plains.push(plain)
        }
    }
    let caps = affine
        .iter()
        .map(|c| c as &dyn CapLiteral)
        .chain(roots.iter().map(|c| c as &dyn CapLiteral))
        .collect::<Vec<_>>();
    let hr = helices.iter().collect::<Vec<_>>();
    let pr = plains.iter().collect::<Vec<_>>();
    let mut source_position = 0_f64;
    let mut support_derivatives = [0_f64; 2];
    let mut original_derivative = 0_f64;
    for surface in surfaces.surfaces() {
        l.charge(1)?;
        let bound = surface.source_bounds();
        source_position = source_position.max(bound.position_error_mm);
        let turns = match field.meridians()[surface.meridian()].component() {
            crate::authored_hanger_bolt::Component::Metric => {
                I::point(source.metric().full_length_mm())
                    / (I::point(source.metric().spec().pitch_mm)
                        * I::point(f64::from(source.metric().spec().starts)))
            }
            crate::authored_hanger_bolt::Component::Wood => {
                I::point(source.wood().full_length_mm())
                    / (I::point(source.wood().spec().pitch_mm)
                        * I::point(f64::from(source.wood().spec().starts)))
            }
        };
        original_derivative = original_derivative
            .max((I::point(bound.derivative_error_mm[1]) * (turns + I::point(2.))).hi);
        for (to, value) in support_derivatives
            .iter_mut()
            .zip(bound.derivative_error_mm)
        {
            *to = to.max(value);
        }
    }
    // An original whole helix coedge spans no more than this many turns;
    // compact V stays measured in turns. Cap receipts carry their own original
    // parameter conversion, including inherited root-parameter scales.
    for cap in &caps {
        l.charge(1)?;
        let b = cap.source_bounds();
        source_position = source_position.max(b[0]);
        original_derivative = original_derivative.max(b[2]);
    }
    let file = encode_candidate(
        &topology,
        &caps,
        &hr,
        &pr,
        FileOptions {
            maximum_encoding_work: l.allowance(128_000_000, 0)?,
            maximum_stage_work: 16_000_000,
            maximum_bytes: l.bytes(original + slots)?,
        },
        i,
    )?;
    l.charge(file.encoding_work_used())?;
    Ok((
        file,
        source_position,
        support_derivatives,
        original_derivative,
    ))
}
/// The complete operation is bounded, cancellable and atomic. Source encoding
/// arenas are dropped before readback. Returned bytes passed literal incidence,
/// vertex links, whole finite lifts, exact polar factorization, true cap-plane
/// closure and healed support-derived volume checks. No header recipe is read.
pub fn export_qualified_hanger(
    source: &Hanger,
    options: QualifiedExportOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<QualifiedHangerStep, Error> {
    if options.maximum_work == 0
        || options.maximum_work > 1_000_000_000
        || options.maximum_bytes == 0
        || options.maximum_bytes > 64 * 1024 * 1024
        || !options.maximum_healed_volume_width_mm3.is_finite()
        || options.maximum_healed_volume_width_mm3 <= 0.
    {
        return Err(Error::Invalid);
    }
    let mut l = Ledger {
        used: 0,
        left: options.maximum_work,
        bytes: options.maximum_bytes,
        interrupt,
    };
    l.charge(1)?;
    let (candidate, source_position, support_derivatives, original_derivative) =
        produce(source, &mut l)?;
    let extra =
        candidate.retained_bytes() - candidate.as_str().len() + size_of::<QualifiedHangerStep>();
    let bytes = l.bytes(extra)?;
    let decoded = read::decode(
        candidate.as_str(),
        read::ReadOptions {
            maximum_work: l.allowance(128_000_000, 0)?,
            retained_budget_bytes: bytes,
        },
        interrupt,
    )
    .map_err(read_error)?;
    l.charge(decoded.work_used())?;
    let grid = read::ExactPoleGrid::prepare(
        &decoded,
        read::GridOptions {
            maximum_preparation_work: l.allowance(224_000_000, 0)?,
            maximum_bytes: bytes,
        },
        interrupt,
    )
    .map_err(read_error)?;
    l.charge(grid.preparation_work_used())?;
    let polar = read::PolarSupports::prepare(
        &grid,
        read::PolarOptions {
            maximum_work: l.allowance(64_000_000, 0)?,
            maximum_bytes: bytes,
        },
        interrupt,
    )
    .map_err(read_error)?;
    l.charge(polar.preparation_work_used())?;
    let lifts = decoded
        .surface_lifts(
            read::LiftOptions {
                maximum_work: l.allowance(160_000_000, decoded.work_used())?,
                maximum_bytes: l
                    .bytes(extra + polar.retained_bytes() - decoded.retained_budget_bytes())?,
                ..read::LiftOptions::default()
            },
            interrupt,
        )
        .map_err(read_error)?;
    l.charge(lifts.work_used() - decoded.work_used())?;
    let healed = read::HealedMaterial::prepare(
        &polar,
        &lifts,
        read::HealingOptions {
            maximum_work: l.allowance(64_000_000, 0)?,
            maximum_bytes: bytes,
            ..read::HealingOptions::default()
        },
        interrupt,
    )
    .map_err(read_error)?;
    l.charge(healed.admission_work_used())?;
    let volume = healed
        .healed_volume(
            read::VolumeOptions {
                maximum_work: l.allowance(16_000_000, 0)?,
                maximum_bytes: bytes,
                maximum_width_mm3: options.maximum_healed_volume_width_mm3,
            },
            interrupt,
        )
        .map_err(read_error)?;
    l.charge(volume.work_used())?;
    let volume = volume.bounds_mm3();
    let sewing = healed.sewing_error_bound_mm();
    Ok(QualifiedHangerStep {
        source: *source,
        candidate,
        work: l.used,
        bytes: options.maximum_bytes,
        volume,
        sewing,
        source_position,
        support_derivatives,
        original_derivative,
    })
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn actual_authored_hangers_export_only_after_independent_healed_qualification() {
        for (left, nb) in [(false, false), (true, true)] {
            let source =
                Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
            let started = std::time::Instant::now();
            let file =
                export_qualified_hanger(&source, QualifiedExportOptions::default(), None).unwrap();
            eprintln!(
                "atomic healed hanger export left={left} seconds={} bytes={} work={} accounted_bound={} healed_volume={:?} sewing={}",
                started.elapsed().as_secs_f64(),
                file.as_str().len(),
                file.work_used(),
                file.peak_accounted_bytes(),
                file.healed_volume_bounds_mm3(),
                file.sewing_correspondence_bound_mm()
            );
            assert!(file.as_str().starts_with("ISO-10303-21;"));
            assert!(file.healed_volume_bounds_mm3()[0] > 0.);
            assert!(file.sewing_correspondence_bound_mm() <= 1e-7);
            assert!(file.peak_accounted_bytes() <= QualifiedExportOptions::default().maximum_bytes);
            assert!(matches!(
                export_qualified_hanger(
                    &source,
                    QualifiedExportOptions {
                        maximum_work: file.work_used() - 1,
                        ..QualifiedExportOptions::default()
                    },
                    None
                ),
                Err(Error::Limit)
            ));
            assert!(matches!(
                export_qualified_hanger(
                    &source,
                    QualifiedExportOptions {
                        maximum_work: 1,
                        ..QualifiedExportOptions::default()
                    },
                    None
                ),
                Err(Error::Limit)
            ));
            assert!(matches!(
                export_qualified_hanger(
                    &source,
                    QualifiedExportOptions {
                        maximum_bytes: 1,
                        ..QualifiedExportOptions::default()
                    },
                    None
                ),
                Err(Error::Limit)
            ));
            assert!(matches!(
                export_qualified_hanger(&source, QualifiedExportOptions::default(), Some(&|| true)),
                Err(Error::Cancelled)
            ));
            let calls = std::cell::Cell::new(0);
            assert!(matches!(
                export_qualified_hanger(
                    &source,
                    QualifiedExportOptions::default(),
                    Some(&|| {
                        let n = calls.get() + 1;
                        calls.set(n);
                        n >= 20000
                    })
                ),
                Err(Error::Cancelled)
            ));
            assert!(calls.get() >= 20000);
            assert!(matches!(
                export_qualified_hanger(
                    &source,
                    QualifiedExportOptions {
                        maximum_work: 0,
                        ..QualifiedExportOptions::default()
                    },
                    None
                ),
                Err(Error::Invalid)
            ));
            eprintln!(
                "qualified source residual left={left} raw={} healed={} supportd1={:?} originald1={}",
                file.source_to_raw_position_error_mm(),
                file.source_to_healed_position_error_mm(),
                file.source_support_derivative_error_mm(),
                file.original_coedge_derivative_error_mm()
            );
            assert!(file.source_to_healed_position_error_mm() <= 1e-7);
            assert!(
                file.source_support_derivative_error_mm()
                    .iter()
                    .all(|v| *v <= 1e-6)
            );
            assert!(file.original_coedge_derivative_error_mm() <= 1e-6);
            assert_eq!(file.source(), &source);
            let mut single = crate::step::Writer::new();
            let appended = single
                .append_qualified_hanger(&file, QualifiedExportOptions::default(), None)
                .unwrap();
            assert_eq!(
                appended.source_to_raw_position_error_mm(),
                file.source_to_raw_position_error_mm()
            );
            assert_eq!(
                appended.source_to_healed_position_error_mm(),
                file.source_to_healed_position_error_mm()
            );
            assert_eq!(
                appended.source_support_derivative_error_mm(),
                file.source_support_derivative_error_mm()
            );
            assert_eq!(
                appended.original_coedge_derivative_error_mm(),
                file.original_coedge_derivative_error_mm()
            );
            assert!(appended.work_used() > file.work_used());
            assert!(appended.append_work_used() > 0);
            assert!(
                appended.peak_accounted_bytes() <= QualifiedExportOptions::default().maximum_bytes
            );
            assert_eq!(
                appended.healed_volume_bounds_mm3(),
                file.healed_volume_bounds_mm3()
            );
            assert_eq!(
                appended.sewing_correspondence_bound_mm(),
                file.sewing_correspondence_bound_mm()
            );
            let baseline = single.finish().unwrap();
            read::decode(&baseline, read::ReadOptions::default(), None).unwrap();
            let mut writer = crate::step::Writer::new();
            writer
                .append_qualified_hanger(&file, QualifiedExportOptions::default(), None)
                .unwrap();
            let calls = std::cell::Cell::new(0);
            assert!(matches!(
                writer.append_qualified_hanger(
                    &file,
                    QualifiedExportOptions::default(),
                    Some(&|| {
                        let n = calls.get() + 1;
                        calls.set(n);
                        n >= 40000
                    })
                ),
                Err(Error::Cancelled)
            ));
            assert!(calls.get() >= 40000);
            assert_eq!(
                writer.finish().unwrap(),
                baseline,
                "cancel after partial literal append rolls back all geometry/IDs"
            );
            let mut writer = crate::step::Writer::new();
            assert!(matches!(
                writer.append_qualified_hanger(
                    &file,
                    QualifiedExportOptions {
                        maximum_work: appended.work_used() - 1,
                        ..QualifiedExportOptions::default()
                    },
                    None
                ),
                Err(Error::Limit)
            ));
            assert!(matches!(
                writer.append_qualified_hanger(
                    &file,
                    QualifiedExportOptions {
                        maximum_bytes: 1,
                        ..QualifiedExportOptions::default()
                    },
                    None
                ),
                Err(Error::Limit)
            ));
            writer
                .append_qualified_hanger(&file, QualifiedExportOptions::default(), None)
                .unwrap();
            writer
                .append_qualified_hanger(&file, QualifiedExportOptions::default(), None)
                .unwrap();
            let compound = writer.finish().unwrap();
            assert_eq!(compound.matches("MANIFOLD_SOLID_BREP(").count(), 2);
            let mut ids = std::collections::BTreeSet::new();
            let data = compound.split_once("\nDATA;\n").unwrap().1;
            for line in data.lines().take_while(|l| l.starts_with('#')) {
                assert!(
                    ids.insert(
                        line.split_once('=').unwrap().0[1..]
                            .parse::<usize>()
                            .unwrap()
                    )
                );
            }
            let mut quoted = false;
            let bytes = data.as_bytes();
            let mut cursor = 0;
            while cursor < bytes.len() {
                if bytes[cursor] == b'\'' {
                    if quoted && bytes.get(cursor + 1) == Some(&b'\'') {
                        cursor += 2;
                        continue;
                    }
                    quoted = !quoted;
                    cursor += 1;
                    continue;
                }
                if bytes[cursor] == b'#' && !quoted {
                    let start = cursor + 1;
                    cursor = start;
                    while bytes.get(cursor).is_some_and(u8::is_ascii_digit) {
                        cursor += 1
                    }
                    assert!(ids.contains(&data[start..cursor].parse::<usize>().unwrap()));
                } else {
                    cursor += 1
                }
            }
            eprintln!(
                "qualified replay left={left} appendwork={} peak={} compoundbytes={}",
                appended.append_work_used(),
                appended.peak_accounted_bytes(),
                compound.len()
            );
            let text = file.into_string();
            assert!(text.ends_with("END-ISO-10303-21;\n"));
        }
    }
}
