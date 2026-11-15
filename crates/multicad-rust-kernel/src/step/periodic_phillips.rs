//! Export-only proposal followed by independent actual-file qualification.
//! The admitted exchange interpretation is the distinctly repaired decoded
//! material, not exact occupied material of the unmodified literal surfaces.
use crate::{step::PhillipsStepError, step_read};
use step_read::{PhillipsOwnedPeriodicMaterial as Material, PhillipsReadError};

#[derive(Clone, Copy)]
pub struct Options {
    pub maximum_proposal_work: usize,
    pub maximum_output_bytes: usize,
    pub reader: step_read::PhillipsPeriodicReadOptions,
    pub correspondence: step_read::PhillipsPeriodicMaterialMatchOptions,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            maximum_proposal_work: 16_000_000,
            maximum_output_bytes: 64 * 1024 * 1024,
            reader: Default::default(),
            correspondence: Default::default(),
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Invalid,
    Limit,
    Cancelled,
    Unsupported,
    Reader(step_read::PhillipsPeriodicReadError),
}
fn source_error(e: PhillipsReadError) -> Error {
    match e {
        PhillipsReadError::Invalid => Error::Invalid,
        PhillipsReadError::Limit => Error::Limit,
        PhillipsReadError::Cancelled => Error::Cancelled,
        PhillipsReadError::Unsupported => Error::Unsupported,
    }
}
fn checkpoint(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), Error> {
    if interrupt.is_some_and(|f| f()) {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Verification {
    /// Complete boundary correspondence between the expected owned source and
    /// the independently proved decoded repair, including cap annuli/world map.
    pub maximum_world_distance_mm: f64,
    /// Parameter derivatives, not unit-normal angles.
    pub maximum_world_chart_derivative_mm: f64,
    pub independent_volume_discrepancy_mm3: f64,
    pub decoded_volume_bounds_mm3: [f64; 2],
    pub reserved_reader_work: usize,
    pub decoded_material_proof_work: usize,
    pub correspondence_fixed_work: usize,
    pub correspondence_cavity_cells: usize,
}
pub struct QualifiedMaterial {
    material: Material,
    verification: Verification,
}
impl QualifiedMaterial {
    pub fn material(&self) -> &Material {
        &self.material
    }
    pub fn verification(&self) -> Verification {
        self.verification
    }
    pub fn into_material(self) -> Material {
        self.material
    }
}
pub struct QualifiedStep {
    text: String,
    qualification: QualifiedMaterial,
    proposal_work: usize,
}
impl QualifiedStep {
    pub fn as_str(&self) -> &str {
        &self.text
    }
    pub fn qualification(&self) -> &QualifiedMaterial {
        &self.qualification
    }
    pub fn proposal_work_used(&self) -> usize {
        self.proposal_work
    }
}
/// Qualifies these ACTUAL bytes independently; a valid but stale/different file
/// is refused by whole-source correspondence, not a writer marker or mass alone.
pub fn qualify_actual_file(
    source: &Material,
    text: &str,
    occurrence: [[f64; 3]; 4],
    options: Options,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<QualifiedMaterial, Error> {
    checkpoint(interrupt)?;
    if options.maximum_output_bytes == 0 || options.maximum_output_bytes > 64 * 1024 * 1024 {
        return Err(Error::Invalid);
    }
    if text.len() > options.maximum_output_bytes {
        return Err(Error::Limit);
    }
    let read = step_read::read_periodic_phillips_material(text, options.reader, interrupt)
        .map_err(Error::Reader)?;
    let matched = source
        .match_material(
            read.material(),
            occurrence,
            options.correspondence,
            interrupt,
        )
        .map_err(source_error)?;
    let verification = Verification {
        maximum_world_distance_mm: matched.maximum_world_distance_mm(),
        maximum_world_chart_derivative_mm: matched.maximum_world_chart_derivative_mm(),
        independent_volume_discrepancy_mm3: matched.independent_volume_discrepancy_mm3(),
        decoded_volume_bounds_mm3: read.material().world_volume_mm3(),
        reserved_reader_work: read.reserved_delegated_work(),
        decoded_material_proof_work: read.material_proof_work_used(),
        correspondence_fixed_work: matched.fixed_work_used(),
        correspondence_cavity_cells: matched.cavity_cells_used(),
    };
    checkpoint(interrupt)?;
    Ok(QualifiedMaterial {
        material: read.into_material(),
        verification,
    })
}
/// Normal authoring geometry remains the input owned source. This serializes
/// only an untrusted export proposal, then earns its distinct decoded material.
/// The byte cap controls publication; original writer limits control encoding.
/// No aggregate process-memory or scalar total-work claim combines the separate
/// proposal/reader/fixed-correspondence/adaptive-cell resource dimensions.
pub fn export(
    source: &Material,
    occurrence: [[f64; 3]; 4],
    options: Options,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<QualifiedStep, Error> {
    checkpoint(interrupt)?;
    let proposal = source
        .mesh_approximation(options.maximum_proposal_work, interrupt)
        .map_err(source_error)?;
    let chart = source.chart();
    // Ordinary floating arithmetic produces ONLY a candidate placement. The
    // complete actual symbolic chart+occurrence is independently compared below.
    let linear = |p: [f64; 3]| {
        std::array::from_fn(|i| (0..3).map(|j| occurrence[j + 1][i] * p[j]).sum::<f64>())
    };
    let offset = linear(std::array::from_fn(|i| {
        chart.columns[0][i] - chart.columns[3][i] * chart.outer_top_z
    }));
    let origin = std::array::from_fn(|i| occurrence[0][i] + offset[i]);
    let u = linear(chart.columns[1]);
    let v = linear(chart.columns[2]);
    let text = crate::step::write_phillips_screw_with_interrupt(
        proposal.proposal(),
        origin,
        u,
        v,
        1e-7,
        interrupt,
    )
    .map_err(|e| match e {
        PhillipsStepError::InvalidGeometry | PhillipsStepError::InvalidPlacement => Error::Invalid,
        PhillipsStepError::Limit => Error::Limit,
        PhillipsStepError::Cancelled => Error::Cancelled,
        PhillipsStepError::Unsupported => Error::Unsupported,
    })?;
    let qualification = qualify_actual_file(source, &text, occurrence, options, interrupt)?;
    checkpoint(interrupt)?;
    Ok(QualifiedStep {
        text,
        qualification,
        proposal_work: proposal.work_used(),
    })
}

#[cfg(test)]
mod tests;
