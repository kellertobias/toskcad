//! Source-bound sign domains of the actual repaired thread rim and head circle.
//! This is not a sewn material, occupied proof, manufacturing mesh or file mass.
use super::*;
mod degree11;
pub use degree11::{
    Degree11ContactOptions, Degree11DomainEnd, Degree11RimContact, Degree11RimContactId,
    Degree11RimKnot, Degree11SignDomain, Degree11ThreadContacts,
};
use std::mem::size_of;
#[cfg(test)]
mod tests;
use num_bigint::BigInt as Z;
use num_traits::{One, Signed, Zero};
pub(super) mod exact;
pub(super) mod roots;
use exact::{Arithmetic, Poly};
pub use roots::DyadicParameter;

const ROOTS: usize = 512;
// Every defining candidate net has height>=37 and width=9U+1. The existing
// aggregate262144-control guard thus bounds sum(U)<=787; no new source cap.
const CELLS: usize = MAX_POINTS / (9 * 37);
// Accounted exact-solver scratch, not numeric-bitlength/RSS estimation.
// Arithmetic preflights every operand/result at <=8192bits. With pinned
// num-bigint0.4.8, BigDigits::clone copies normalized LEN, not prior capacity;
// normalize/shrink never grows it. Add/multiply/shift allocate <=two operand
// extents. Division first normalizes to <=129u64/257u32 limbs; Burnikel's
// initial block level is <=128/256, and its concat/shift/recursive blocks and
// Vec growth stay <=four levels (4096bytes on either target). Thus8KiB per
// retained BigInt covers even a prior larger operation allocation plus header.
// No uncharged larger integer enters this sealed degree18 solver.
// Depth60 DFS retains <=61 sibling polynomials, plus original/current/split/
// evaluation grids: <=66*19 integers (<10MiB at8KiB each). The square-free
// fallback DROPS that frontier before its <=19-row Sturm chain (<400 live
// coefficients including remainder copies). <=9 recursive arithmetic levels
// use <1MiB extra; frontier/Vec headers, outputs, conversion and operation
// temporaries fit the remaining>5MiB. These are accounted allocation bounds,
// with unchanged bit/depth/work/cancellation guards, not allocator RSS.
// Reserve before construction; actual compact receipt capacities are separate.
const EXACT_SCRATCH: usize = 16 * 1024 * 1024;
/// Bounded construction and contact query, sharing the original source work
/// and byte limits. The event quota may only lower the hard512-contact cap.
#[derive(Clone, Copy)]
pub struct ThreadContactOptions {
    /// C0/D1/control limits; work/bytes cover the ENTIRE contact request,
    /// including construction and the separately described exact scratch.
    pub source: DerivedThreadOptions,
    pub maximum_contacts: usize,
}
impl Default for ThreadContactOptions {
    fn default() -> Self {
        Self {
            source: Default::default(),
            maximum_contacts: ROOTS,
        }
    }
}
/// Exact source ownership, not a coordinate-proximity equivalence. IDs are
/// local to the owning ContactThreadCarrier, never persistent part references.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RimKnot {
    Vertex(usize),
    SplineKnot { curve: usize, knot: usize },
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RimContactId {
    Knot(RimKnot),
    Interior {
        curve: usize,
        span: usize,
        ordinal: usize,
    },
}
/// A distinct root of the original exact binary64 coefficient polynomial.
/// Its isolator does not replace the algebraic endpoint by a sampled point.
#[derive(Clone, Copy, Debug)]
pub struct RimContact {
    id: RimContactId,
    curve: usize,
    span: usize,
    bracket: [DyadicParameter; 2],
}
impl RimContact {
    pub fn source_identity(&self) -> RimContactId {
        self.id
    }
    pub fn source_curve(&self) -> usize {
        self.curve
    }
    pub fn source_span(&self) -> usize {
        self.span
    }
    /// Isolates the exact algebraic endpoint in the indicated normalized
    /// Bezier span. The whole isolator remains part of the adjoining domains;
    /// its endpoints are NOT replacement trim parameters.
    pub fn isolating_parameter(&self) -> [DyadicParameter; 2] {
        self.bracket
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RimDomainEnd {
    Knot(RimKnot),
    Contact(usize),
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RimDomainRole {
    /// Gamma lies inside the circle: head-bottom annulus candidate, intended+Z after occupied union proof.
    HeadBottom,
    /// Gamma lies outside the circle: thread-top annulus candidate, intended-Z after occupied union proof.
    ThreadTop,
}
#[derive(Clone, Copy, Debug)]
pub struct RimSignDomain {
    curve: usize,
    span: usize,
    ends: [RimDomainEnd; 2],
    role: RimDomainRole,
}
impl RimSignDomain {
    pub fn source_curve(&self) -> usize {
        self.curve
    }
    pub fn source_span(&self) -> usize {
        self.span
    }
    pub fn ends(&self) -> [RimDomainEnd; 2] {
        self.ends
    }
    pub fn role(&self) -> RimDomainRole {
        self.role
    }
}
struct Report {
    contacts: Vec<RimContact>,
    domains: Vec<RimSignDomain>,
}
impl Report {
    fn bytes(&self) -> usize {
        self.contacts.capacity() * size_of::<RimContact>()
            + self.domains.capacity() * size_of::<RimSignDomain>()
    }
}
/// Owns the same measured degree9 thread candidate and exact sign-domain proof.
/// The strict fitted factory remains separate. This receipt alone cannot
/// certify occupied material, finite face filling, global sewing or STEP mass.
pub struct ContactThreadCarrier<'a, 'input> {
    carrier: DerivedThreadCarrier<'a, 'input>,
    report: Report,
    head_radius: f64,
    work: usize,
    retained: usize,
}
impl<'a, 'input> ContactThreadCarrier<'a, 'input> {
    pub fn body(&self) -> &'a Body<'input> {
        self.carrier.body()
    }
    pub fn thread_carrier(&self) -> &DerivedThreadCarrier<'a, 'input> {
        &self.carrier
    }
    /// Exact retained circle radius in the common local chart, not a new
    /// Euclidean-world radius after its literal four-column transport.
    pub fn head_circle_radius_mm(&self) -> f64 {
        self.head_radius
    }
    pub fn contacts(&self) -> &[RimContact] {
        &self.report.contacts
    }
    pub fn sign_domains(&self) -> &[RimSignDomain] {
        &self.report.domains
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.retained
    }
}
fn descriptor_reserve(maximum_contacts: usize) -> usize {
    size_of::<ContactThreadCarrier<'_, '_>>() - size_of::<DerivedThreadCarrier<'_, '_>>()
        + 2 * maximum_contacts * size_of::<RimContact>()
        + 2 * (CELLS + maximum_contacts) * size_of::<RimSignDomain>()
}
fn knot(rim: usize, knot: usize, spans: usize, rims: usize) -> RimKnot {
    if knot == 0 {
        RimKnot::Vertex(2 * rim)
    } else if knot == spans {
        RimKnot::Vertex(2 * ((rim + 1) % rims))
    } else {
        RimKnot::SplineKnot {
            curve: 3 * rim,
            knot,
        }
    }
}
// Parameter ownership only: this does not assert that a supplied parameter
// is a geometric root. The sole production caller supplies the earned root
// of this SAME immutable rim cell; no coordinate comparison is involved.
fn root_identity(
    rim: usize,
    span: usize,
    spans: usize,
    rims: usize,
    ordinal: usize,
    bracket: [DyadicParameter; 2],
) -> (RimContactId, bool) {
    let at_start = bracket == [DyadicParameter::ZERO; 2];
    let at_end = bracket == [DyadicParameter::ONE; 2];
    let id = if at_start {
        RimContactId::Knot(knot(rim, span, spans, rims))
    } else if at_end {
        RimContactId::Knot(knot(rim, span + 1, spans, rims))
    } else {
        RimContactId::Interior {
            curve: 3 * rim,
            span,
            ordinal,
        }
    };
    (id, at_start)
}
fn intern(
    result: &mut Report,
    contact: RimContact,
    maximum: usize,
    work: &mut Work<'_>,
) -> Result<usize> {
    // One candidate/Body owns the entire builder. Repeated endpoint IDs thus
    // refer to the exact shared source vertex/knot earned by graph::build.
    work.step(1 + result.contacts.len())?;
    if let Some(i) = result.contacts.iter().position(|r| r.id == contact.id) {
        return Ok(i);
    }
    if result.contacts.len() == maximum {
        return Err(ReadError::Limit);
    }
    let index = result.contacts.len();
    result.contacts.push(contact);
    Ok(index)
}
fn report(
    nets: &[Net],
    radius: f64,
    maximum_contacts: usize,
    work: &mut Work<'_>,
) -> Result<Report> {
    work.step(1)?; // Poll/charge before the first descriptor allocation.
    let cells = nets.iter().map(|n| n.spans[0]).sum::<usize>();
    // Source MAX_POINTS and each net's minimum height prove cells<=CELLS.
    let mut result = Report {
        contacts: Vec::with_capacity(maximum_contacts),
        domains: Vec::with_capacity(cells + maximum_contacts),
    };
    let mut arithmetic = Arithmetic { work };
    let mut budget = roots::Budget::new();
    for (rim, net) in nets.iter().enumerate() {
        let h = 9 * net.spans[1] + 1;
        for span in 0..net.spans[0] {
            let poles = std::array::from_fn::<_, 10, _>(|i| net.poles[(9 * span + i) * h]);
            let polynomial = arithmetic.radius_polynomial(poles, radius)?;
            let certificate = roots::certify(polynomial, &mut arithmetic, &mut budget)?;
            let mut ends = Vec::with_capacity(certificate.roots.len() + 2);
            ends.push(RimDomainEnd::Knot(knot(
                rim,
                span,
                net.spans[0],
                nets.len(),
            )));
            for (ordinal, root) in certificate.roots.iter().enumerate() {
                let (id, at_start) =
                    root_identity(rim, span, net.spans[0], nets.len(), ordinal, root.bracket);
                let index = intern(
                    &mut result,
                    RimContact {
                        id,
                        curve: 3 * rim,
                        span,
                        bracket: root.bracket,
                    },
                    maximum_contacts,
                    arithmetic.work,
                )?;
                if at_start {
                    ends[0] = RimDomainEnd::Contact(index);
                }
                ends.push(RimDomainEnd::Contact(index));
            }
            let final_knot = knot(rim, span + 1, net.spans[0], nets.len());
            let final_end = certificate
                .roots
                .last()
                .filter(|r| r.bracket == [DyadicParameter::ONE; 2])
                .map(|_| *ends.last().unwrap())
                .unwrap_or(RimDomainEnd::Knot(final_knot));
            ends.push(final_end);
            for domain in certificate.domains {
                arithmetic.work.step(1)?;
                result.domains.push(RimSignDomain {
                    curve: 3 * rim,
                    span,
                    ends: domain.ends.map(|i| ends[i]),
                    // Strict signed domains only; an identically zero cell
                    // refuses in the finite-event solver before binding.
                    role: if domain.sign < 0 {
                        RimDomainRole::HeadBottom
                    } else {
                        RimDomainRole::ThreadTop
                    },
                });
            }
        }
    }
    Ok(result)
}
fn publish<'a, 'input>(
    carrier: DerivedThreadCarrier<'a, 'input>,
    report: Report,
    head_radius: f64,
    maximum_bytes: usize,
    original_work: usize,
    mut work: Work<'_>,
) -> std::result::Result<ContactThreadCarrier<'a, 'input>, PhillipsReadError> {
    let retained = carrier
        .retained_budget_bytes()
        .saturating_add(
            size_of::<ContactThreadCarrier<'_, '_>>() - size_of::<DerivedThreadCarrier<'_, '_>>(),
        )
        .saturating_add(report.bytes());
    if retained.saturating_add(EXACT_SCRATCH) > maximum_bytes {
        return Err(PhillipsReadError::Limit);
    }
    work.step(1).map_err(|e| bound_error(e, &work))?;
    Ok(ContactThreadCarrier {
        carrier,
        report,
        head_radius,
        work: original_work - work.remaining,
        retained,
    })
}
impl<'h, 'q, 'r, 'a, 'input> ClosedRepairedHead<'h, 'q, 'r, 'a, 'input> {
    /// Explicitly distinct measured degree9 candidate with nominal-contact
    /// domains. U and V quarters are exact restrictions of the original nets;
    /// every C0/D1 transfer remains measured in the ORIGINAL cell parameter.
    /// All construction/contact work shares the existing16M maximum. The64MiB
    /// accounted live budget reserves16MiB for bounded exact scratch, not RSS.
    pub fn contact_thread_carrier(
        &self,
        options: ThreadContactOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<ContactThreadCarrier<'a, 'input>, PhillipsReadError> {
        let mut initial = Work {
            remaining: options.source.maximum_work,
            interrupt,
            cancelled: Cell::new(false),
        };
        initial.step(1).map_err(|e| bound_error(e, &initial))?;
        if !valid_options(options.source) || options.maximum_contacts > ROOTS {
            return Err(PhillipsReadError::Invalid);
        }
        let maximum_bytes = options
            .source
            .maximum_bytes
            .checked_sub(EXACT_SCRATCH + descriptor_reserve(options.maximum_contacts))
            .ok_or(PhillipsReadError::Limit)?;
        let component = DerivedThreadOptions {
            maximum_work: initial.remaining,
            maximum_bytes,
            ..options.source
        };
        let (carrier, report) =
            self.with_candidate(component, interrupt, 4, |nets, radius, _, work| {
                report(nets, radius, options.maximum_contacts, work)
            })?;
        initial.remaining -= carrier.work_used();
        publish(
            carrier,
            report,
            self.virtual_bottom_station()[0],
            options.source.maximum_bytes,
            options.source.maximum_work,
            initial,
        )
    }
}
