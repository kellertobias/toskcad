//! Independent fixed-frame physical head packets over actual clipped D11 nets.
//! The auxiliary standalone file is DERIVED verification data, never a raw
//! physical cap. Source-law matching and Native capability remain separate.
use super::head::{HeadInput, OccupiedClippedHead};
use super::*;
use crate::step_read::owned_insert::DocumentOptions;
use crate::step_read::radial_boundary::pipe_hardware::budget as cold;
use std::ops::Range;
mod authored;
mod auxiliary;
pub use authored::{
    AuthoredHeadExport, AuthoredHeadReceipt, OccurrenceFactor, export_authored_hex_head,
    qualify_authored_hex_head_file,
};
mod graph;
pub mod measure;
pub mod mesh;
mod owned;
pub mod query;
pub use owned::{AxialSewingTransfer, OwnedActualHead, OwnedHeadOptions, OwnedHeadReceipt};
mod source;
#[cfg(test)]
pub(in crate::step_read) mod source_tests;
pub mod transport;
pub use source::{HexHeadIntent, QualifiedFile, SourceOptions, SourceReceipt, qualify_actual_file};
#[cfg(test)]
mod tests;

const MAX_WORK: usize = 272_000_000;
const DECODER_WORK: usize = 64_000_000;
const MAX_BYTES: usize = 64 * 1024 * 1024;
#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub maximum_work: usize,
    pub maximum_bytes: usize,
    pub additional_resident_bytes: usize,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            maximum_work: MAX_WORK,
            maximum_bytes: MAX_BYTES,
            additional_resident_bytes: 0,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Receipt {
    pub physical_cold_work_reserved: usize,
    pub fresh_decoding_work: usize,
    pub derived_cold_work_reserved: usize,
    pub aggregate_work: usize,
    pub accounted_peak_bytes: usize,
}
#[derive(Clone, Debug)]
pub struct ComponentBinding {
    entity: u32,
    component: usize,
    actual: Range<usize>,
    derived: Range<usize>,
}
impl ComponentBinding {
    pub fn entity(&self) -> u32 {
        self.entity
    }
    pub fn component(&self) -> usize {
        self.component
    }
    pub fn actual_range(&self) -> Range<usize> {
        self.actual.clone()
    }
    pub fn derived_range(&self) -> Range<usize> {
        self.derived.clone()
    }
}
/// Actual physical record identities, distinct from the derived cap entity.
#[derive(Clone, Copy, Debug)]
pub struct PhysicalRoles {
    pub solid_entity: u32,
    pub annulus_face_entity: u32,
    pub thread_face_entities: [u32; 9],
    /// Orders correspond exactly to HeadInput faces/edges/vertices.
    pub head_face_entities: [u32; 8],
    pub head_edge_entities: [u32; 18],
    pub head_vertex_entities: [u32; 12],
}
#[derive(Clone, Copy, Debug)]
pub struct CapChart {
    pub physical_support_entity: u32,
    pub derived_support_entity: u32,
    pub physical_height_mm: f64,
    pub derived_height_bounds_mm: [f64; 2],
    /// XY is the identical actual polynomial chart. Only this directed Z
    /// displacement separates physical annulus and auxiliary verification cap.
    pub z_displacement_bounds_mm: [f64; 2],
}
pub struct DecodedFile<'input> {
    actual: &'input str,
    derived: String,
    head: HeadInput,
    bindings: Vec<ComponentBinding>,
    chart: CapChart,
    roles: PhysicalRoles,
    receipt: Receipt,
    retained_beyond_auxiliary_input: usize,
    maximum_work: usize,
    maximum_bytes: usize,
    derived_parser_work_reserved: usize,
}
pub struct OccupiedFile<'packet, 'parsed, 'document, 'auxiliary, 'actual> {
    file: &'packet DecodedFile<'actual>,
    material: OccupiedClippedHead<'parsed, 'document, 'auxiliary>,
    receipt: Receipt,
}
impl OccupiedFile<'_, '_, '_, '_, '_> {
    pub fn material(&self) -> &OccupiedClippedHead<'_, '_, '_> {
        &self.material
    }
    pub fn actual_file(&self) -> &str {
        self.file.actual
    }
    pub fn receipt(&self) -> Receipt {
        self.receipt
    }
    pub fn cap_chart(&self) -> CapChart {
        self.file.chart
    }
}
fn geometry_error(e: crate::axial_pipe_thread::Error) -> Error {
    match e {
        crate::axial_pipe_thread::Error::Invalid => Error::Invalid,
        crate::axial_pipe_thread::Error::Unsupported
        | crate::axial_pipe_thread::Error::Unresolved => Error::Unsupported,
        crate::axial_pipe_thread::Error::Limit => Error::Limit,
        crate::axial_pipe_thread::Error::Cancelled => Error::Cancelled,
    }
}
fn neutral_error(e: crate::step_read::owned_insert::Error) -> Error {
    match e {
        crate::step_read::owned_insert::Error::Invalid => Error::Invalid,
        crate::step_read::owned_insert::Error::Unsupported => Error::Unsupported,
        crate::step_read::owned_insert::Error::Limit => Error::Limit,
        crate::step_read::owned_insert::Error::Cancelled => Error::Cancelled,
    }
}
impl<'input> DecodedFile<'input> {
    pub fn parse(
        input: &'input str,
        options: Options,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self> {
        Self::parse_with_capture(input, options, interrupt, |_, _| Ok(((), 0)))
            .map(|(packet, ())| packet)
    }
    fn parse_with_capture<T>(
        input: &'input str,
        options: Options,
        interrupt: Option<&dyn Fn() -> bool>,
        capture: impl FnOnce(&CompoundDocument<'input>, &mut Budget<'_>) -> Result<(T, usize)>,
    ) -> Result<(Self, T)> {
        if options.maximum_work == 0
            || options.maximum_work > MAX_WORK
            || options.maximum_bytes == 0
            || options.maximum_bytes > MAX_BYTES
        {
            return Err(Error::Invalid);
        }
        let bytes = options
            .maximum_bytes
            .checked_sub(options.additional_resident_bytes)
            .filter(|n| *n > 0)
            .ok_or(Error::Limit)?;
        let original_cold = cold::reserve(
            input,
            options.maximum_work.min(64_000_000),
            bytes,
            interrupt,
        )
        .map_err(geometry_error)?;
        let remaining = options
            .maximum_work
            .checked_sub(original_cold.aggregate_work_reserved)
            .ok_or(Error::Limit)?;
        let document = CompoundDocument::parse(
            input,
            DocumentOptions {
                maximum_work: original_cold.delegated_parser_work_reserved,
                retained_budget_bytes: bytes,
            },
            interrupt,
        )
        .map_err(neutral_error)?;
        // This decoder owns a distinct64M ledger; standalone Parsed retains
        // its existing16M ceiling. Only the checked storage/step primitives are
        // shared, not the standalone constructor's scope or geometry proof.
        let fresh = remaining.min(DECODER_WORK);
        let resident = document
            .retained_budget_bytes()
            .checked_add(options.additional_resident_bytes)
            .and_then(|n| n.checked_add(SCRATCH))
            .ok_or(Error::Limit)?;
        if resident > options.maximum_bytes {
            return Err(Error::Limit);
        }
        let mut budget = Budget {
            used: 0,
            remaining: fresh,
            bytes: resident,
            peak: resident,
            limit: options.maximum_bytes,
            interrupt,
        };
        let plan = graph::decode(&document, &mut budget)?;
        let result = auxiliary::derive(&document, plan, &mut budget)?;
        let (captured, captured_bytes) = capture(&document, &mut budget)?;
        let decode_work = budget.used;
        let physical_peak = budget
            .peak
            .max(original_cold.parser_peak_bytes_reserved + options.additional_resident_bytes);
        // All original-arena consumers end here. Never retain two cold ASTs.
        drop(document);
        let used = original_cold
            .aggregate_work_reserved
            .checked_add(decode_work)
            .ok_or(Error::Limit)?;
        let left = options
            .maximum_work
            .checked_sub(used)
            .filter(|n| *n > 0)
            .ok_or(Error::Limit)?;
        let derived_live = input
            .len()
            .checked_add(result.retained_without_derived)
            .and_then(|n| n.checked_add(captured_bytes))
            .and_then(|n| n.checked_add(result.text.capacity() - result.text.len()))
            .and_then(|n| n.checked_add(size_of::<Self>()))
            .and_then(|n| n.checked_add(options.additional_resident_bytes))
            .ok_or(Error::Limit)?;
        let derived_bytes = options
            .maximum_bytes
            .checked_sub(derived_live)
            .filter(|n| *n > 0)
            .ok_or(Error::Limit)?;
        let derived_cold =
            cold::reserve(&result.text, left.min(64_000_000), derived_bytes, interrupt)
                .map_err(geometry_error)?;
        let aggregate = used
            .checked_add(derived_cold.aggregate_work_reserved)
            .ok_or(Error::Limit)?;
        if aggregate > options.maximum_work {
            return Err(Error::Limit);
        }
        let retained = input
            .len()
            .checked_add(result.retained_without_derived)
            .and_then(|n| n.checked_add(captured_bytes))
            .and_then(|n| n.checked_add(result.text.capacity() - result.text.len()))
            .and_then(|n| n.checked_add(size_of::<Self>()))
            .and_then(|n| n.checked_add(options.additional_resident_bytes))
            .ok_or(Error::Limit)?;
        let peak = physical_peak.max(derived_live + derived_cold.parser_peak_bytes_reserved);
        if peak > options.maximum_bytes {
            return Err(Error::Limit);
        }
        budget.step(0)?;
        Ok((
            Self {
                actual: input,
                derived: result.text,
                head: result.head,
                bindings: result.bindings,
                chart: result.chart,
                roles: result.roles,
                receipt: Receipt {
                    physical_cold_work_reserved: original_cold.aggregate_work_reserved,
                    fresh_decoding_work: decode_work,
                    derived_cold_work_reserved: derived_cold.aggregate_work_reserved,
                    aggregate_work: aggregate,
                    accounted_peak_bytes: peak,
                },
                retained_beyond_auxiliary_input: retained,
                maximum_work: options.maximum_work,
                maximum_bytes: options.maximum_bytes,
                derived_parser_work_reserved: derived_cold.delegated_parser_work_reserved,
            },
            captured,
        ))
    }
    pub fn actual_file(&self) -> &'input str {
        self.actual
    }
    /// Explicitly DERIVED standalone verification graph, not a physical file.
    pub fn derived_auxiliary_text(&self) -> &str {
        &self.derived
    }
    pub fn component_bindings(&self) -> &[ComponentBinding] {
        &self.bindings
    }
    pub fn cap_chart(&self) -> CapChart {
        self.chart
    }
    pub fn physical_roles(&self) -> PhysicalRoles {
        self.roles
    }
    pub fn receipt(&self) -> Receipt {
        self.receipt
    }
    pub fn certify<'packet, 'parsed, 'document, 'auxiliary>(
        &'packet self,
        parsed: &'parsed Parsed<'document, 'auxiliary>,
        tolerance_mm: f64,
        minimum_wall_mm: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<OccupiedFile<'packet, 'parsed, 'document, 'auxiliary, 'input>> {
        let literal = parsed.document.literal_input();
        if literal.as_ptr() != self.derived.as_ptr() || literal.len() != self.derived.len() {
            return Err(Error::Invalid);
        }
        // Pointer identity only seals packet association. Physical/raw source
        // correspondence was earned by full graph decoding and exact component
        // copies, including the separately checked support/chart rehoming.
        let clone_work = 128
            + self.head.vertices.len() * 3
            + self.head.edges.len() * 2
            + self
                .head
                .faces
                .iter()
                .map(|f| f.outer.len() * 2 + f.inner.len() * 2 + 9)
                .sum::<usize>();
        if interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        let inherited = self
            .receipt
            .aggregate_work
            .checked_add(clone_work)
            .ok_or(Error::Limit)?;
        let head_max = self
            .maximum_work
            .checked_sub(inherited)
            .and_then(|n| n.checked_add(parsed.document.work_used()))
            .ok_or(Error::Limit)?
            .min(96_000_000);
        let head_bytes = head::storage(&self.head)?;
        let clone_peak = parsed
            .document
            .retained_budget_bytes()
            .checked_add(self.retained_beyond_auxiliary_input)
            .and_then(|n| n.checked_add(parsed.query_budget_bytes()))
            .and_then(|n| n.checked_add(head_bytes))
            .ok_or(Error::Limit)?;
        if clone_peak > self.maximum_bytes {
            return Err(Error::Limit);
        }
        let material = head::certify(
            parsed,
            self.head.clone(),
            tolerance_mm,
            minimum_wall_mm,
            head::HeadOptions {
                maximum_work: head_max,
                maximum_bytes: self.maximum_bytes,
                additional_resident_bytes: self.retained_beyond_auxiliary_input,
            },
            interrupt,
        )?;
        let aggregate = inherited
            .checked_add(material.receipt().aggregate_work)
            .and_then(|n| n.checked_sub(parsed.document.work_used()))
            .ok_or(Error::Limit)?;
        if aggregate > self.maximum_work {
            return Err(Error::Limit);
        }
        let peak = self
            .receipt
            .accounted_peak_bytes
            .max(material.receipt().accounted_peak_bytes);
        if interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        Ok(OccupiedFile {
            file: self,
            material,
            receipt: Receipt {
                aggregate_work: aggregate,
                accounted_peak_bytes: peak,
                ..self.receipt
            },
        })
    }
}
