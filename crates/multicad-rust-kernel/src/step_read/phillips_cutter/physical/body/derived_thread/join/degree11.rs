//! Original degree11 rim/head-circle sign domains. This is not material or STEP
//! admission. Algebraic roots remain owned by exact decoded coefficient cells.
use super::super::degree11::CoherentDegree11ThreadCarrier;
use super::*;
const TOTAL_WORK: usize = 64_000_000;
const CELLS11: usize = MAX_POINTS / (11 * 12);
// Same8192-bit arithmetic preflights/depth60 as the degree18 solver. Its DFS
// retains<=66*23 integers (<13MiB at the previously audited8KiB perBigInt).
// The frontier is dropped before Sturm's<=23-row chain:<=529 coefficients,
// even with current/remainder copies and recursive arithmetic<8MiB. Reserve
//20MiB for either path plus fixed headers/temporaries; allocator RSS is separate.
const SCRATCH22: usize = 20 * 1024 * 1024;
#[derive(Clone, Copy)]
pub struct Degree11ContactOptions {
    /// NEW contact-stage work only; never renews a whole-import quota.
    pub maximum_work: usize,
    /// Aggregate head + measured-carrier + new-contact stages. The caller's
    /// other import-stage charges remain separate, explicitly accounted work.
    pub maximum_total_work: usize,
    pub maximum_bytes: usize,
    pub maximum_contacts: usize,
}
impl Default for Degree11ContactOptions {
    fn default() -> Self {
        Self {
            maximum_work: WORK,
            maximum_total_work: TOTAL_WORK,
            maximum_bytes: BYTES,
            maximum_contacts: ROOTS,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Degree11RimKnot {
    /// Original VERTEX_POINT owner of the measured derived rim endpoint. This
    /// names topology provenance, not equality to the raw coordinate after
    /// the carrier's explicitly bounded end/rail correction.
    OriginalVertex(u32),
    SurfaceKnot {
        surface_entity: u32,
        original_u_knot: usize,
    },
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Degree11RimContactId {
    Knot(Degree11RimKnot),
    Interior {
        surface_entity: u32,
        original_u_cell: usize,
        ordinal: usize,
    },
}
#[derive(Clone, Copy, Debug)]
pub struct Degree11RimContact {
    pub source_identity: Degree11RimContactId,
    pub face_entity: u32,
    pub surface_entity: u32,
    pub original_u_cell: usize,
    pub isolating_parameter: [DyadicParameter; 2],
}
#[derive(Clone, Copy, Debug)]
pub enum Degree11DomainEnd {
    Knot(Degree11RimKnot),
    Contact(usize),
}
#[derive(Clone, Copy, Debug)]
pub struct Degree11SignDomain {
    pub face_entity: u32,
    pub surface_entity: u32,
    pub original_u_cell: usize,
    pub ends: [Degree11DomainEnd; 2],
    pub role: RimDomainRole,
}
pub struct Degree11ThreadContacts<'c, 'a, 'input> {
    carrier: &'c CoherentDegree11ThreadCarrier<'a, 'input>,
    radius: f64,
    contacts: Vec<Degree11RimContact>,
    domains: Vec<Degree11SignDomain>,
    new_work: usize,
    retained: usize,
}
impl<'c, 'a, 'input> Degree11ThreadContacts<'c, 'a, 'input> {
    pub fn carrier(&self) -> &'c CoherentDegree11ThreadCarrier<'a, 'input> {
        self.carrier
    }
    pub fn head_circle_radius_mm(&self) -> f64 {
        self.radius
    }
    pub fn contacts(&self) -> &[Degree11RimContact] {
        &self.contacts
    }
    pub fn sign_domains(&self) -> &[Degree11SignDomain] {
        &self.domains
    }
    pub fn new_work_used(&self) -> usize {
        self.new_work
    }
    pub fn total_work_used(&self) -> usize {
        self.carrier.total_work_used() + self.new_work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.retained
    }
}
fn knot11(
    head: &ClosedRepairedHead<'_, '_, '_, '_, '_>,
    carrier: &CoherentDegree11ThreadCarrier<'_, '_>,
    rim: usize,
    k: usize,
) -> Degree11RimKnot {
    let faces = carrier.faces();
    if k == 0 {
        Degree11RimKnot::OriginalVertex(
            head.common_mouth().tensor_transfers()[rim].entry_vertex_entity,
        )
    } else if k == faces[rim].source().spans()[0] {
        Degree11RimKnot::OriginalVertex(
            head.common_mouth().tensor_transfers()[(rim + 1) % faces.len()].entry_vertex_entity,
        )
    } else {
        Degree11RimKnot::SurfaceKnot {
            surface_entity: faces[rim].transfer().surface_entity,
            original_u_knot: k,
        }
    }
}
impl<'h, 'q, 'r, 'a, 'input> ClosedRepairedHead<'h, 'q, 'r, 'a, 'input> {
    pub fn degree11_thread_contacts<'c>(
        &self,
        carrier: &'c CoherentDegree11ThreadCarrier<'a, 'input>,
        options: Degree11ContactOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<Degree11ThreadContacts<'c, 'a, 'input>, PhillipsReadError> {
        let mut work = Work {
            remaining: options.maximum_work,
            interrupt,
            cancelled: Cell::new(false),
        };
        let result = (|| {
            work.step(1)?;
            if options.maximum_work > WORK
                || options.maximum_total_work > TOTAL_WORK
                || options.maximum_bytes > BYTES
                || options.maximum_contacts > ROOTS
                || !std::ptr::eq(self.body(), carrier.body())
                || self.chart() != carrier.chart()
                || carrier.cap_stations_mm()[0] != self.virtual_bottom_station()[1]
            {
                return Err(ReadError::Invalid);
            }
            let known = carrier.total_work_used();
            let allowance = options
                .maximum_total_work
                .checked_sub(known)
                .ok_or(ReadError::Limit)?
                .min(options.maximum_work);
            work.remaining = allowance.checked_sub(1).ok_or(ReadError::Limit)?;
            let reserve = std::mem::size_of::<Degree11ThreadContacts<'_, '_, '_>>()
                + 2 * options.maximum_contacts * std::mem::size_of::<Degree11RimContact>()
                + 2 * (CELLS11 + options.maximum_contacts)
                    * std::mem::size_of::<Degree11SignDomain>();
            if carrier
                .retained_budget_bytes()
                .saturating_add(reserve)
                .saturating_add(SCRATCH22)
                > options.maximum_bytes
            {
                return Err(ReadError::Limit);
            }
            let mut contacts: Vec<Degree11RimContact> =
                Vec::with_capacity(options.maximum_contacts);
            let mut domains = Vec::with_capacity(CELLS11 + options.maximum_contacts);
            let radius = self.virtual_bottom_station()[0];
            let mut arithmetic = Arithmetic { work: &mut work };
            let mut root_budget = roots::Budget::new();
            for (rim, face) in carrier.faces().iter().enumerate() {
                let transfer = face.transfer();
                let height = face.source().dimensions()[1];
                for span in 0..face.source().spans()[0] {
                    let poles =
                        std::array::from_fn(|u| face.source().poles()[(11 * span + u) * height]);
                    let polynomial = arithmetic.radius_polynomial11(poles, radius)?;
                    let certificate =
                        roots::certify22(polynomial, &mut arithmetic, &mut root_budget)?;
                    let mut ends = Vec::with_capacity(certificate.roots.len() + 2);
                    ends.push(Degree11DomainEnd::Knot(knot11(self, carrier, rim, span)));
                    for (ordinal, root) in certificate.roots.iter().enumerate() {
                        arithmetic.work.step(1 + contacts.len())?;
                        let at_start = root.bracket == [DyadicParameter::ZERO; 2];
                        let at_end = root.bracket == [DyadicParameter::ONE; 2];
                        let id = if at_start {
                            Degree11RimContactId::Knot(knot11(self, carrier, rim, span))
                        } else if at_end {
                            Degree11RimContactId::Knot(knot11(self, carrier, rim, span + 1))
                        } else {
                            Degree11RimContactId::Interior {
                                surface_entity: transfer.surface_entity,
                                original_u_cell: span,
                                ordinal,
                            }
                        };
                        let index = if let Some(index) =
                            contacts.iter().position(|c| c.source_identity == id)
                        {
                            index
                        } else {
                            if contacts.len() == options.maximum_contacts {
                                return Err(ReadError::Limit);
                            }
                            let index = contacts.len();
                            contacts.push(Degree11RimContact {
                                source_identity: id,
                                face_entity: transfer.face_entity,
                                surface_entity: transfer.surface_entity,
                                original_u_cell: span,
                                isolating_parameter: root.bracket,
                            });
                            index
                        };
                        if at_start {
                            ends[0] = Degree11DomainEnd::Contact(index);
                        }
                        ends.push(Degree11DomainEnd::Contact(index));
                    }
                    let final_end = certificate
                        .roots
                        .last()
                        .filter(|r| r.bracket == [DyadicParameter::ONE; 2])
                        .map(|_| *ends.last().unwrap())
                        .unwrap_or(Degree11DomainEnd::Knot(knot11(
                            self,
                            carrier,
                            rim,
                            span + 1,
                        )));
                    ends.push(final_end);
                    for d in certificate.domains {
                        arithmetic.work.step(1)?;
                        domains.push(Degree11SignDomain {
                            face_entity: transfer.face_entity,
                            surface_entity: transfer.surface_entity,
                            original_u_cell: span,
                            ends: d.ends.map(|i| ends[i]),
                            role: if d.sign < 0 {
                                RimDomainRole::HeadBottom
                            } else {
                                RimDomainRole::ThreadTop
                            },
                        });
                    }
                }
            }
            let retained = carrier.retained_budget_bytes()
                + std::mem::size_of::<Degree11ThreadContacts<'_, '_, '_>>()
                + contacts.capacity() * std::mem::size_of::<Degree11RimContact>()
                + domains.capacity() * std::mem::size_of::<Degree11SignDomain>();
            work.step(1)?;
            Ok(Degree11ThreadContacts {
                carrier,
                radius,
                contacts,
                domains,
                new_work: allowance - work.remaining,
                retained,
            })
        })();
        result.map_err(|e| bound_error(e, &work))
    }
}
#[cfg(test)]
mod tests;
