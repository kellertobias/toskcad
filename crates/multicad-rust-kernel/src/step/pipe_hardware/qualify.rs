//! Independent actual-file occupied admission plus complete original-operator
//! correspondence. A writer trace is comparison data, never Source admission.
use super::*;
use crate::{
    axial_pipe_thread::{
        Error,
        authored::hardware::shell::{self, occupied::OccupiedHead},
    },
    radial_tensor::{self as radial, Surface},
    step_read::{literal_member_binding as binding, pipe_hardware as reader},
};
use reader::{ParameterMapping as Mapping, Target};
type Result<T> = std::result::Result<T, Error>;
// Metered binding/operator work uses the checked remainder of one explicit
// aggregate. Unchanged delegated reader costs remain separately reserved.
const MAX_WORK: usize = 272_000_000;
#[derive(Clone, Copy, Debug)]
pub struct QualificationOptions {
    pub maximum_work: usize,
    pub maximum_bytes: usize,
    pub additional_resident_bytes: usize,
}
impl Default for QualificationOptions {
    fn default() -> Self {
        Self {
            maximum_work: MAX_WORK,
            maximum_bytes: 64 * 1024 * 1024,
            additional_resident_bytes: 0,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct QualificationReceipt {
    pub proposal_work_used: usize,
    pub independent_read_work_reserved: usize,
    pub correspondence_work_used: usize,
    pub binding_work_used: usize,
    pub operator_check_work_used: usize,
    pub aggregate_work_reserved: usize,
    pub matched_original_records: usize,
    pub matched_curves: usize,
    pub matched_surfaces: usize,
    /// Exact operator identity under earned parameter maps, NOT zero raw-to-
    /// sewn error and NOT an authored-law/standard/Core-preset certificate.
    pub raw_operator_position_delta_mm: f64,
    /// Spline/tensor, canonical circle/cone and finite line derivatives under
    /// their earned maps only. Plane support is implicit; no Plane UV or Plane
    /// parameter-derivative identity is certified by cached Axis equality.
    pub mapped_operator_derivative_delta_mm: f64,
    pub source_join_collar_mm: f64,
    pub decoded_join_collar_mm: f64,
    pub proposal_peak_bytes_reserved: usize,
    pub read_peak_bytes_reserved: usize,
    pub correspondence_peak_bytes_reserved: usize,
}
pub struct QualifiedStep {
    actual: reader::ReadHead,
    receipt: QualificationReceipt,
}
impl QualifiedStep {
    pub fn as_str(&self) -> &str {
        self.actual.original_input()
    }
    pub fn material(&self) -> &reader::ReadHead {
        &self.actual
    }
    pub fn receipt(&self) -> QualificationReceipt {
        self.receipt
    }
    pub fn into_material(self) -> reader::ReadHead {
        self.actual
    }
}
fn binding_error(e: binding::Error) -> Error {
    match e {
        binding::Error::Invalid => Error::Invalid,
        binding::Error::Different => Error::Unsupported,
        binding::Error::Cancelled => Error::Cancelled,
        _ => Error::Limit,
    }
}
pub fn qualify_actual_file(
    source: &OccupiedHead,
    actual: &str,
    tolerance: f64,
    minimum: f64,
    options: QualificationOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<QualifiedStep> {
    if options.maximum_work == 0
        || options.maximum_work > MAX_WORK
        || options.maximum_bytes == 0
        || options.maximum_bytes > 64 * 1024 * 1024
    {
        return Err(Error::Invalid);
    }
    let outside = options
        .additional_resident_bytes
        .checked_add(actual.len())
        .ok_or(Error::Limit)?;
    let expected = propose(
        source,
        Options {
            maximum_work: options.maximum_work.min(32_000_000),
            maximum_bytes: options.maximum_bytes,
            additional_resident_bytes: outside,
        },
        interrupt,
    )?;
    let proposal = expected.receipt();
    let allowance = options
        .maximum_work
        .checked_sub(proposal.work_used)
        .ok_or(Error::Limit)?
        .min(224_000_000);
    let reader_resident = source
        .receipt()
        .retained_result_bytes_upper_bound
        .checked_add(expected.retained_bytes()?)
        .and_then(|n| n.checked_add(options.additional_resident_bytes))
        .ok_or(Error::Limit)?;
    let decoded = reader::read(
        actual,
        tolerance,
        minimum,
        reader::Options {
            maximum_work: allowance,
            maximum_bytes: options.maximum_bytes,
            additional_resident_bytes: reader_resident,
        },
        interrupt,
    )?;
    let live = reader_resident
        .checked_add(decoded.retained_bytes_upper_bound())
        .and_then(|n| n.checked_add(actual.len()))
        .ok_or(Error::Limit)?;
    let scratch = options
        .maximum_bytes
        .checked_sub(live)
        .ok_or(Error::Limit)?;
    let correspondence = options
        .maximum_work
        .checked_sub(proposal.work_used)
        .and_then(|n| n.checked_sub(decoded.receipt().aggregate_work_reserved))
        .ok_or(Error::Limit)?;
    if correspondence == 0 {
        return Err(Error::Limit);
    }
    let bound = binding::bind_retained(
        &expected,
        expected.root,
        &decoded.records,
        decoded.records.member,
        binding::Options {
            maximum_work: correspondence.min(512_000_000),
            maximum_bytes: scratch,
            maximum_entities: 65_536,
        },
        interrupt,
    )
    .map_err(|e| {
        #[cfg(test)]
        eprintln!("BIND_REFUSAL {e:?} live={live} scratch={scratch}");
        binding_error(e)
    })?;
    let mut work = Work {
        remaining: correspondence
            .checked_sub(bound.work_used())
            .ok_or(Error::Limit)?,
        interrupt,
    };
    compare(
        source.shell().input(),
        decoded.material().shell().input(),
        decoded.parameter_maps(),
        &mut work,
    )?;
    let read = decoded.receipt();
    let peak = read
        .parser_peak_bytes_reserved
        .max(read.decode_peak_bytes_reserved)
        .max(read.occupied_peak_bytes_reserved);
    let correspondence_peak = live
        .checked_add(bound.query_peak_bytes())
        .and_then(|n| n.checked_add(32 * 1024))
        .ok_or(Error::Limit)?;
    if correspondence_peak > options.maximum_bytes {
        return Err(Error::Limit);
    }
    let receipt = QualificationReceipt {
        proposal_work_used: proposal.work_used,
        independent_read_work_reserved: read.aggregate_work_reserved,
        correspondence_work_used: correspondence - work.remaining,
        binding_work_used: bound.work_used(),
        operator_check_work_used: correspondence - bound.work_used() - work.remaining,
        aggregate_work_reserved: proposal
            .work_used
            .checked_add(read.aggregate_work_reserved)
            .and_then(|n| n.checked_add(correspondence - work.remaining))
            .ok_or(Error::Limit)?,
        matched_original_records: bound.matched_entities(),
        matched_curves: source.shell().input().curves.len(),
        matched_surfaces: source.shell().input().surfaces.len(),
        raw_operator_position_delta_mm: 0.,
        mapped_operator_derivative_delta_mm: 0.,
        source_join_collar_mm: source.join_sewing_mm(),
        decoded_join_collar_mm: decoded.material().join_sewing_mm(),
        proposal_peak_bytes_reserved: proposal.phase_peak_bytes_reserved,
        read_peak_bytes_reserved: peak,
        correspondence_peak_bytes_reserved: correspondence_peak,
    };
    if receipt.aggregate_work_reserved > options.maximum_work {
        return Err(Error::Limit);
    }
    Ok(QualifiedStep {
        actual: decoded,
        receipt,
    })
}
struct Work<'a> {
    remaining: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Work<'_> {
    fn step(&mut self) -> Result<()> {
        self.charge(1)
    }
    fn charge(&mut self, n: usize) -> Result<()> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        self.remaining = self.remaining.checked_sub(n).ok_or(Error::Limit)?;
        Ok(())
    }
}
struct Pairs<const N: usize> {
    forward: [usize; N],
    reverse: [usize; N],
}
impl<const N: usize> Pairs<N> {
    fn new() -> Self {
        Self {
            forward: [usize::MAX; N],
            reverse: [usize::MAX; N],
        }
    }
    fn pair(&mut self, a: usize, b: usize) -> Result<()> {
        let f = self.forward.get_mut(a).ok_or(Error::Limit)?;
        let r = self.reverse.get_mut(b).ok_or(Error::Limit)?;
        if (*f != usize::MAX && *f != b) || (*r != usize::MAX && *r != a) {
            return Err(Error::Unsupported);
        }
        *f = b;
        *r = a;
        Ok(())
    }
    fn covers(&self, a: usize, b: usize) -> bool {
        a <= N
            && b <= N
            && self.forward[..a].iter().all(|x| *x != usize::MAX)
            && self.reverse[..b].iter().all(|x| *x != usize::MAX)
    }
}
fn map<'a>(
    maps: &'a [reader::ParameterMap],
    target: Target,
    work: &mut Work<'_>,
) -> Result<&'a Mapping> {
    let mut found = None;
    for m in maps {
        work.step()?;
        if m.target() == target && found.replace(m.mapping()).is_some() {
            return Err(Error::Unsupported);
        }
    }
    found.ok_or(Error::Unsupported)
}
fn compare(
    a: &shell::Input,
    b: &shell::Input,
    maps: &[reader::ParameterMap],
    work: &mut Work<'_>,
) -> Result<()> {
    if a.faces.len() != b.faces.len()
        || a.edges.len() != b.edges.len()
        || a.vertices.len() != b.vertices.len()
        || a.curves.len() != b.curves.len()
        || a.surfaces.len() != b.surfaces.len()
    {
        return Err(Error::Unsupported);
    }
    let (mut vs, mut es, mut cs, mut ss) = (
        Pairs::<96>::new(),
        Pairs::<128>::new(),
        Pairs::<128>::new(),
        Pairs::<48>::new(),
    );
    for (x, y) in a.faces.iter().zip(&b.faces) {
        work.step()?;
        if x.forward != y.forward || x.inner.len() != y.inner.len() {
            return Err(Error::Unsupported);
        }
        ss.pair(x.surface, y.surface)?;
        if let Surface::Tensor(t) = &a.surfaces[x.surface] {
            work.charge(t.poles().len())?;
        }
        if a.surfaces[x.surface] != b.surfaces[y.surface] {
            return Err(Error::Unsupported);
        }
        match &a.surfaces[x.surface] {
            Surface::Tensor(t) => {
                if map(maps, Target::Surface(y.surface), work)?
                    != &(Mapping::Tensor { spans: t.spans() })
                {
                    return Err(Error::Unsupported);
                }
            }
            Surface::Cone { semi_angle_rad, .. } => {
                let Mapping::Cone {
                    axial_to_slant_scale,
                } = map(maps, Target::Surface(y.surface), work)?
                else {
                    return Err(Error::Unsupported);
                };
                let c = crate::tapered_thread::interval::I::point(*semi_angle_rad)
                    .cos()
                    .positive()?;
                let scale = (crate::tapered_thread::interval::I::point(1.) / c).checked()?;
                if *axial_to_slant_scale != [scale.lo, scale.hi] {
                    return Err(Error::Unsupported);
                }
            }
            Surface::Plane(_) => {}
        }
        for (ar, br) in std::iter::once(&x.outer)
            .chain(&x.inner)
            .zip(std::iter::once(&y.outer).chain(&y.inner))
        {
            if ar.len() != br.len() {
                return Err(Error::Unsupported);
            }
            for (&(e, d), &(f, g)) in ar.iter().zip(br) {
                work.step()?;
                if d != g {
                    return Err(Error::Unsupported);
                }
                es.pair(e, f)?;
                let (x, y) = (a.edges[e], b.edges[f]);
                if x.forward != y.forward {
                    return Err(Error::Unsupported);
                }
                cs.pair(x.curve, y.curve)?;
                for k in 0..2 {
                    vs.pair(x.vertices[k], y.vertices[k])?;
                    if a.vertices[x.vertices[k]] != b.vertices[y.vertices[k]] {
                        return Err(Error::Unsupported);
                    }
                }
            }
        }
    }
    if !vs.covers(a.vertices.len(), b.vertices.len())
        || !es.covers(a.edges.len(), b.edges.len())
        || !cs.covers(a.curves.len(), b.curves.len())
        || !ss.covers(a.surfaces.len(), b.surfaces.len())
    {
        return Err(Error::Unsupported);
    }
    for (i, c) in a.curves.iter().enumerate() {
        work.step()?;
        let j = cs.forward[i];
        match (c, &b.curves[j]) {
            (shell::Curve::Original(x), shell::Curve::Original(y)) => {
                if let radial::Curve::Spline { poles, .. } = x {
                    work.charge(poles.len())?;
                }
                if x != y {
                    return Err(Error::Unsupported);
                }
                match x {
                    radial::Curve::Spline { spans, .. } => {
                        if map(maps, Target::Curve(j), work)?
                            != &(Mapping::Spline { spans: *spans })
                        {
                            return Err(Error::Unsupported);
                        }
                    }
                    radial::Curve::Circle { .. } => {
                        if map(maps, Target::Curve(j), work)? != &Mapping::AngularIdentity {
                            return Err(Error::Unsupported);
                        }
                    }
                }
            }
            (shell::Curve::Line(x), shell::Curve::Line(y)) => {
                for k in 0..2 {
                    if vs.forward[x[k]] != y[k] {
                        return Err(Error::Unsupported);
                    }
                }
                let Mapping::Line {
                    delta,
                    length_bounds_mm,
                } = map(maps, Target::Curve(j), work)?
                else {
                    return Err(Error::Unsupported);
                };
                let [p, q] = x.map(|v| a.vertices[v]);
                for k in 0..3 {
                    if !numbers::exact_difference(q[k], p[k], delta[k])? {
                        return Err(Error::Unsupported);
                    }
                }
                let n = delta
                    .map(crate::tapered_thread::interval::I::point)
                    .into_iter()
                    .fold(crate::tapered_thread::interval::I::point(0.), |s, d| {
                        s + d.square()
                    })
                    .sqrt()
                    .positive()?;
                if *length_bounds_mm != [n.lo, n.hi] {
                    return Err(Error::Unsupported);
                }
            }
            _ => return Err(Error::Unsupported),
        }
    }
    // Full coefficient equality plus actual integer knots establishes the
    // identity on EVERY cell, not samples: STEP(n*u)=source(u). Differentiating
    // includes n in each spline/tensor axis. Canonical analytic frames and the
    // retained cone sec(angle)/line-length maps give their full pullbacks.
    Ok(())
}
