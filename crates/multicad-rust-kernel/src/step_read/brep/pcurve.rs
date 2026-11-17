//! Whole-curve agreement of every stated two-dimensional surface trim.
//! A chosen 3D master does not erase supplementary PCURVE geometry.
use super::*;
use crate::tapered_thread::interval::I;
mod affine;
mod analytic;
mod precision;
mod sweep;
mod uv;

const MAX_WORK: usize = 512_000_000;
const MAX_PAIRS: usize = 100_000;
const MAX_WORKSPACE: usize = 64 * 1024 * 1024;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PcurveValidationReport {
    pub edges_checked: usize,
    pub pcurves_checked: usize,
    pub certified_maximum_deviation_mm: f64,
    /// Fresh geometry/checking actions plus one visit per input byte. The
    /// existing bounded Part21 parser is a separate delegated parse scope.
    pub work_used: usize,
}
#[derive(Clone, Debug, PartialEq)]
pub struct PcurveDisagreement {
    pub surface_entity: u32,
    pub edge_entity: u32,
    pub pcurve_entity: u32,
    pub surface_kind: &'static str,
    /// A proved discrepancy witness in millimetres, never the requested limit.
    pub deviation_mm: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub struct PcurveValidationError {
    pub kind: ReadError,
    pub disagreements: Vec<PcurveDisagreement>,
    pub cancelled: bool,
}
impl From<ReadError> for PcurveValidationError {
    fn from(kind: ReadError) -> Self {
        Self {
            kind,
            disagreements: Vec::new(),
            cancelled: false,
        }
    }
}

pub(super) enum Proof {
    /// Upper bound on the complete finite geometric correspondence.
    Certified(f64),
    /// Lower witness for a genuine geometric discrepancy. A failed upper
    /// estimate alone can NEVER construct this variant.
    Disagreement(f64),
    Unsupported,
}
pub(super) struct Budget<'a> {
    remaining: usize,
    initial: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
    cancelled: std::cell::Cell<bool>,
}
impl<'a> Budget<'a> {
    pub(super) fn new(work: usize, interrupt: Option<&'a dyn Fn() -> bool>) -> Self {
        Self {
            remaining: work,
            initial: work,
            interrupt,
            cancelled: std::cell::Cell::new(false),
        }
    }
    pub(super) fn poll(&self) -> Result<()> {
        if self.interrupt.is_some_and(|f| f()) {
            self.cancelled.set(true);
            Err(ReadError::Limit)
        } else {
            Ok(())
        }
    }
    pub(super) fn charge(&mut self, count: usize) -> Result<()> {
        self.poll()?;
        self.remaining = self.remaining.checked_sub(count).ok_or(ReadError::Limit)?;
        Ok(())
    }
    fn used(&self) -> usize {
        self.initial - self.remaining
    }
}

fn norm(p: [I; 3]) -> I {
    p.into_iter()
        .fold(I::point(0.), |s, x| s + x.square())
        .sqrt()
}
fn vector(p: Point) -> [I; 3] {
    p.map(I::point)
}
fn difference(a: [I; 3], b: [I; 3]) -> [I; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}

/// Proof of the supplementary-curve checks on THIS immutable text at THIS
/// original tolerance. It is not a whole-solid/topology/manufacturing receipt.
pub struct PcurveValidatedInput<'i> {
    input: &'i str,
    tolerance: f64,
    report: PcurveValidationReport,
}
impl std::fmt::Debug for PcurveValidatedInput<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PcurveValidatedInput")
            .field("input_bytes", &self.input.len())
            .field("tolerance", &self.tolerance)
            .field("report", &self.report)
            .finish()
    }
}
impl<'i> PcurveValidatedInput<'i> {
    pub fn input(&self) -> &'i str {
        self.input
    }
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }
    pub fn report(&self) -> &PcurveValidationReport {
        &self.report
    }
}
pub fn validate_pcurves(
    input: &str,
    tolerance: f64,
) -> std::result::Result<PcurveValidationReport, PcurveValidationError> {
    validate_pcurves_with_interrupt(input, tolerance, None)
}
pub fn validate_pcurves_with_interrupt(
    input: &str,
    tolerance: f64,
    interrupt: Option<&dyn Fn() -> bool>,
) -> std::result::Result<PcurveValidationReport, PcurveValidationError> {
    validate_pcurve_input(input, tolerance, interrupt).map(|p| p.report)
}
pub fn validate_pcurve_input<'i>(
    input: &'i str,
    tolerance: f64,
    interrupt: Option<&dyn Fn() -> bool>,
) -> std::result::Result<PcurveValidatedInput<'i>, PcurveValidationError> {
    let mut budget = Budget::new(MAX_WORK, interrupt);
    let result: std::result::Result<PcurveValidationReport, PcurveValidationError> = (|| {
        budget.poll()?;
        check(tolerance.is_finite() && tolerance > 0., ReadError::Invalid)?;
        check(input.len() <= part21::MAX_BYTES, ReadError::Limit)?;
        budget.charge(input.len())?;
        // Necessary lexical cue only: false positives merely cause parsing.
        // Every spelling accepted as a PCURVE contains these ASCII bytes.
        // Chunk overlap retains keywords crossing the polling boundary.
        let bytes = input.as_bytes();
        let mut cue = false;
        for start in (0..bytes.len()).step_by(4096) {
            budget.poll()?;
            let end = (start + 4096 + 6).min(bytes.len());
            if bytes[start..end].windows(7).any(|s| {
                s[..6].eq_ignore_ascii_case(b"PCURVE")
                    && !s[6].is_ascii_alphanumeric()
                    && s[6] != b'_'
            }) {
                cue = true;
                break;
            }
        }
        if !cue {
            return Ok(PcurveValidationReport::default());
        }
        let source = Source::parse(input, tolerance)?;
        budget.poll()?;
        let report = validate_source(&source, tolerance, &mut budget)?;
        drop(source);
        budget.poll()?;
        Ok(report)
    })();
    match result {
        Ok(mut report) => {
            if let Err(kind) = budget.poll() {
                return Err(PcurveValidationError {
                    kind,
                    disagreements: Vec::new(),
                    cancelled: true,
                });
            }
            report.work_used = budget.used();
            Ok(PcurveValidatedInput {
                input,
                tolerance,
                report,
            })
        }
        Err(mut error) => {
            error.cancelled |= budget.cancelled.get();
            Err(error)
        }
    }
}
#[derive(Clone, Copy)]
struct Binding {
    edge: u32,
    curve: u32,
    pcurve: u32,
}
fn add_refs(value: &Value<'_>, pending: &mut Vec<u32>) {
    match value {
        Value::Reference(id) => pending.push(*id),
        Value::List(items) => items.iter().for_each(|v| add_refs(v, pending)),
        Value::Typed(entity) => entity.args.iter().for_each(|v| add_refs(v, pending)),
        _ => {}
    }
}
fn geometry_preflight(source: &Source<'_>, id: u32, budget: &mut Budget<'_>) -> Result<()> {
    let mut seen = BTreeSet::new();
    let mut pending = vec![id];
    while let Some(id) = pending.pop() {
        budget.charge(256)?;
        if !seen.insert(id) {
            continue;
        }
        check(
            seen.len() <= crate::brep::nurbs::MAX_CONTROL_POINTS + 128,
            ReadError::Limit,
        )?;
        for part in source.components(id)? {
            let before = pending.len();
            if matches!(
                part.name,
                "SURFACE_CURVE" | "SEAM_CURVE" | "INTERSECTION_CURVE"
            ) {
                let [_, curve, ..] = part.args.as_slice() else {
                    return Err(ReadError::Invalid);
                };
                pending.push(reference(curve)?);
            } else {
                for arg in &part.args {
                    add_refs(arg, &mut pending);
                }
            }
            let added = pending.len() - before;
            budget.charge(added.saturating_mul(32))?;
            check(
                added <= crate::brep::nurbs::MAX_CONTROL_POINTS + 128
                    && pending.capacity() * 4 <= MAX_WORKSPACE / 8,
                ReadError::Limit,
            )?;
        }
    }
    Ok(())
}
fn unwrap(
    source: &Source<'_>,
    id: u32,
    depth: usize,
    pcs: &mut Vec<u32>,
    budget: &mut Budget<'_>,
) -> Result<u32> {
    budget.charge(256)?;
    check(depth < 8, ReadError::Invalid)?;
    let parts = source.components(id)?;
    let [part] = parts else {
        return Ok(id);
    };
    if !matches!(
        part.name,
        "SURFACE_CURVE" | "SEAM_CURVE" | "INTERSECTION_CURVE"
    ) {
        return Ok(id);
    }
    let [_, curve, associated, Value::Enum(master)] = part.args.as_slice() else {
        return Err(ReadError::Invalid);
    };
    check(
        matches!(*master, "CURVE_3D" | "PCURVE_S1" | "PCURVE_S2"),
        ReadError::Unsupported,
    )?;
    let associated = refs(associated)?;
    check((1..=2).contains(&associated.len()), ReadError::Invalid)?;
    for pc in associated {
        if source.components(pc)?.iter().any(|p| p.name == "PCURVE") {
            pcs.push(pc);
        } else {
            // Associated geometry may name the support surface directly. It
            // states no two-dimensional trim, and is not fabricated into one.
            geometry_preflight(source, pc, budget)?;
            source.surface(pc)?;
        }
    }
    unwrap(source, reference(curve)?, depth + 1, pcs, budget)
}
fn bindings(source: &Source<'_>, budget: &mut Budget<'_>) -> Result<Vec<Binding>> {
    let mut result = Vec::new();
    for (&edge, parts) in &source.entities {
        budget.charge(32)?;
        if !parts
            .iter()
            .any(|e| e.name.eq_ignore_ascii_case("EDGE_CURVE"))
        {
            continue;
        }
        let [_, _, _, curve, _] = source.args(edge, "EDGE_CURVE")? else {
            return Err(ReadError::Invalid);
        };
        let mut pcs = Vec::new();
        let curve = unwrap(source, reference(curve)?, 0, &mut pcs, budget)?;
        for pcurve in pcs {
            check(result.len() < MAX_PAIRS, ReadError::Limit)?;
            result.try_reserve(1).map_err(|_| ReadError::Limit)?;
            check(
                result.capacity() * std::mem::size_of::<Binding>() <= MAX_WORKSPACE / 8,
                ReadError::Limit,
            )?;
            result.push(Binding {
                edge,
                curve,
                pcurve,
            });
        }
    }
    Ok(result)
}
fn edge(source: &Source<'_>, binding: Binding, budget: &mut Budget<'_>) -> Result<Edge> {
    budget.charge(1024)?;
    let [_, a, b, _, same] = source.args(binding.edge, "EDGE_CURVE")? else {
        return Err(ReadError::Invalid);
    };
    let mut vertices = [reference(a)?, reference(b)?];
    if !sense(same)? {
        vertices.reverse();
    }
    let points = vertices.map(|v| {
        let [_, p] = source.args(v, "VERTEX_POINT")? else {
            return Err(ReadError::Invalid);
        };
        source.point(reference(p)?)
    });
    let [a, b] = points;
    geometry_preflight(source, binding.curve, budget)?;
    let curve = source.curve(binding.curve, 0)?;
    // This tolerance is ONLY for locating the existing 3D finite range. All
    // supplementary geometry comparisons still use the ORIGINAL caller limit.
    let incidence = source.declared_uncertainty().map_or(source.tolerance, |d| {
        source.tolerance.max(d.min(MAX_DECLARED_TOLERANCE_MM))
    });
    Edge::new(curve, vertices.map(|v| v as usize), [a?, b?], incidence).map_err(|e| match e {
        crate::brep::Error::Invalid => ReadError::Unsupported,
        other => other.into(),
    })
}
fn angular_axes(surface: &Surface) -> [bool; 2] {
    let angular = |c: &Curve| matches!(c, Curve::Circle { .. } | Curve::Ellipse { .. });
    match surface {
        Surface::Cylinder { .. } | Surface::Cone { .. } => [true, false],
        Surface::Sphere { .. } | Surface::Torus { .. } => [true, true],
        Surface::Extrusion { curve, .. } => [angular(curve), false],
        Surface::Revolution { curve, .. } => [true, angular(curve)],
        _ => [false, false],
    }
}
fn zero_axis(uv: &Curve, axis: usize) -> bool {
    match uv {
        Curve::Line { origin, direction } => origin[axis] == 0. && direction[axis] == 0.,
        Curve::Nurbs { curve, .. } => curve.points().iter().all(|p| p[axis] == 0.),
        _ => false,
    }
}
fn si_radians(source: &Source<'_>, budget: &mut Budget<'_>) -> Result<bool> {
    let mut found = false;
    for parts in source.entities.values() {
        budget.charge(16)?;
        for entity in parts {
            if entity.name != "GLOBAL_UNIT_ASSIGNED_CONTEXT" {
                continue;
            }
            let [units] = entity.args.as_slice() else {
                return Err(ReadError::Invalid);
            };
            for id in refs(units)? {
                if source.component(id, "PLANE_ANGLE_UNIT")?.is_none() {
                    continue;
                }
                found = true;
                if !matches!(
                    source.component(id, "SI_UNIT")?,
                    Some([Value::Omitted, Value::Enum("RADIAN")])
                ) {
                    return Ok(false);
                }
            }
        }
    }
    Ok(found)
}
fn surface_kind(surface: &Surface) -> &'static str {
    match surface {
        Surface::Plane(_) => "PLANE",
        Surface::Cylinder { .. } => "CYLINDRICAL_SURFACE",
        Surface::Cone { .. } => "CONICAL_SURFACE",
        Surface::Sphere { .. } => "SPHERICAL_SURFACE",
        Surface::Torus { .. } => "TOROIDAL_SURFACE",
        Surface::Nurbs { .. } => "B_SPLINE_SURFACE",
        Surface::Extrusion { .. } => "SURFACE_OF_LINEAR_EXTRUSION",
        Surface::Revolution { .. } => "SURFACE_OF_REVOLUTION",
    }
}
fn validate_source(
    source: &Source<'_>,
    tol: f64,
    budget: &mut Budget<'_>,
) -> std::result::Result<PcurveValidationReport, PcurveValidationError> {
    let pairs = bindings(source, budget)?;
    let radians = si_radians(source, budget)?;
    let mut edges = BTreeSet::new();
    let mut disagreements = Vec::new();
    let mut report = PcurveValidationReport::default();
    let mut unqualified = false;
    for binding in pairs {
        budget.charge(256)?;
        edges.insert(binding.edge);
        let attempt = (|| -> Result<(Proof, u32, &'static str, f64)> {
            let edge = edge(source, binding, budget)?;
            let definition = uv::definition(source, binding.pcurve, budget)?;
            geometry_preflight(source, definition.surface, budget)?;
            let surface = source.surface(definition.surface)?;
            if !radians
                && (matches!(surface, Surface::Cone { .. })
                    || angular_axes(&surface)
                        .into_iter()
                        .enumerate()
                        .any(|(i, yes)| yes && !zero_axis(&definition.geometry, i)))
            {
                return Err(ReadError::Unsupported);
            }
            if let Some(proof) = precision::plane_nurbs_certificate(
                source,
                binding.curve,
                &definition,
                &surface,
                &edge,
                tol,
                budget,
            )? {
                return Ok((proof, definition.surface, surface_kind(&surface), 0.));
            }
            if let Some(proof) = precision::plane_analytic_certificate(
                source,
                binding.curve,
                &definition,
                &surface,
                &edge,
                tol,
                budget,
            )? {
                return Ok((proof, definition.surface, surface_kind(&surface), 0.));
            }
            let reserve =
                precision::reserve(source, binding.curve, &definition, &surface, &edge, budget)?;
            if !reserve.is_finite() || reserve < 0. {
                return Err(ReadError::Unsupported);
            }
            let mut proof = affine::certify(&surface, &definition.geometry, &edge, tol, budget)?;
            if matches!(proof, Proof::Unsupported) {
                proof = analytic::certify(&surface, &definition.geometry, &edge, tol, budget)?;
            }
            if matches!(proof, Proof::Unsupported) {
                proof = sweep::certify(&surface, &definition.geometry, &edge, tol, budget)?;
            }
            if matches!(proof, Proof::Certified(_)) {
                precision::sense(source, binding.curve, &definition, &surface, &edge, budget)?;
            }
            Ok((proof, definition.surface, surface_kind(&surface), reserve))
        })();
        match attempt {
            Ok((Proof::Certified(upper), _, _, reserve))
                if upper.is_finite()
                    && upper >= 0.
                    && (I::point(upper) + I::point(reserve)).hi <= tol =>
            {
                report.pcurves_checked += 1;
                report.certified_maximum_deviation_mm = report
                    .certified_maximum_deviation_mm
                    .max((I::point(upper) + I::point(reserve)).hi);
            }
            Ok((Proof::Disagreement(lower), surface, kind, reserve))
                if lower.is_finite() && (I::point(lower) - I::point(reserve)).lo > tol =>
            {
                disagreements.try_reserve(1).map_err(|_| ReadError::Limit)?;
                check(
                    disagreements.capacity() * std::mem::size_of::<PcurveDisagreement>()
                        <= MAX_WORKSPACE / 8,
                    ReadError::Limit,
                )?;
                disagreements.push(PcurveDisagreement {
                    surface_entity: surface,
                    edge_entity: binding.edge,
                    pcurve_entity: binding.pcurve,
                    surface_kind: kind,
                    deviation_mm: (I::point(lower) - I::point(reserve)).lo,
                });
            }
            Err(ReadError::Limit) => return Err(ReadError::Limit.into()),
            Err(ReadError::Invalid) => return Err(ReadError::Invalid.into()),
            _ => unqualified = true,
        }
    }
    // Unused construction PCURVEs are not finite trims of an owned EDGE_CURVE.
    // Every association of each actual edge was examined above, irrespective
    // of master choice; no global construction-entity ownership is invented.
    budget.charge(1)?;
    if !disagreements.is_empty() {
        return Err(PcurveValidationError {
            kind: ReadError::Invalid,
            disagreements,
            cancelled: false,
        });
    }
    if unqualified {
        return Err(ReadError::Unsupported.into());
    }
    report.edges_checked = edges.len();
    budget.charge(1)?;
    Ok(report)
}

#[cfg(test)]
mod tests;
