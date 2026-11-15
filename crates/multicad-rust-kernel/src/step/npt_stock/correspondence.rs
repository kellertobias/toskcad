//! Whole original-source to actual retained decimal geometry correspondence.
//! This gate is local; Native-pose versus STEP-operator transport is separate.
use super::*;
use crate::step_read::npt_stock::OwnedMaterial;
mod topology;
const CORRESPONDENCE_SCRATCH: usize = 16 * 1024;
#[derive(Clone, Copy, Debug)]
pub struct CorrespondenceOptions {
    pub position_error_mm: f64,
    pub original_derivative_error_mm: f64,
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for CorrespondenceOptions {
    fn default() -> Self {
        Self {
            position_error_mm: 1e-7,
            original_derivative_error_mm: 1e-6,
            maximum_work: 96_000_000,
            maximum_bytes: 64 * 1024 * 1024,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct CorrespondenceReceipt {
    matching: [usize; 8],
    edges: [usize; 20],
    vertices: [usize; 12],
    position: f64,
    derivative: [f64; 2],
    edge_position: f64,
    edge_derivative: f64,
    query_work: usize,
    aggregate_work: usize,
    bytes: usize,
}
impl CorrespondenceReceipt {
    /// For each original source patch, the actual retained tensor index.
    pub fn tensor_mapping(&self) -> [usize; 8] {
        self.matching
    }
    pub fn edge_mapping(&self) -> [usize; 20] {
        self.edges
    }
    pub fn vertex_mapping(&self) -> [usize; 12] {
        self.vertices
    }
    pub fn edge_position_error_mm(&self) -> f64 {
        self.edge_position
    }
    pub fn edge_derivative_error_mm(&self) -> f64 {
        self.edge_derivative
    }
    pub fn position_error_mm(&self) -> f64 {
        self.position
    }
    pub fn original_derivative_error_mm(&self) -> [f64; 2] {
        self.derivative
    }
    pub fn query_work_used(&self) -> usize {
        self.query_work
    }
    pub fn aggregate_work_used(&self) -> usize {
        self.aggregate_work
    }
    pub fn accounted_peak_bytes(&self) -> usize {
        self.bytes
    }
}
#[derive(Clone, Copy)]
struct Match {
    position: f64,
    derivative: [f64; 2],
}
fn scalar_delta(source: [f64; 3], target: crate::step_read::npt_stock::Point<'_>) -> f64 {
    source
        .into_iter()
        .zip(target.coordinates())
        .fold(0_f64, |sum, (source, target)| {
            let [lo, hi] = target.bounds();
            let delta = I::new(lo, hi) - I::point(source);
            (sum + delta.magnitude()).next_up()
        })
}
fn compare(
    plan: &SurfacePlan<'_>,
    source: &SurfacePatch,
    target: &crate::step_read::npt_stock::Tensor,
    owned: &OwnedMaterial,
    options: CorrespondenceOptions,
    work: &mut Work<'_>,
) -> Result<Option<Match>, Error> {
    work.charge(64)?;
    if source.spans() != target.spans() {
        return Ok(None);
    }
    // A corner is only a cheap refusal proposal. Admission compares EVERY
    // original pole, with actual decimal leaves enclosed by their own caches.
    let first = owned
        .point(target.point_indices()[0])
        .expect("retained tensor point");
    if scalar_delta(source.poles_mm()[0], first) > options.position_error_mm {
        return Ok(None);
    }
    let mut delta = 0_f64;
    for (source, index) in source.poles_mm().iter().zip(target.point_indices()) {
        work.charge(64)?;
        let target = owned.point(*index).expect("retained tensor point");
        delta = delta.max(scalar_delta(*source, target));
    }
    // Uniform exact dyadic knots and degree11 are independently parsed. A
    // tensor Bernstein convex hull bounds C0 by max_i ||deltaP_i||1. Within
    // ANY original knot cell, its global-u derivative coefficients are
    // 11*spans_u*(deltaP[i+1]-deltaP[i]), hence bounded by22*spans_u*delta.
    // This is a whole-domain L1 bound; no sampled/source-recipe reconstruction.
    let derivative = source
        .spans()
        .map(|n| (I::point(22. * n as f64) * I::point(delta)).hi);
    let mut result = Match {
        position: 0.,
        derivative: [0.; 2],
    };
    for certificate in source.cells() {
        work.charge(32)?;
        result.position = result
            .position
            .max((I::point(certificate.position_error_mm()) + I::point(delta)).hi);
        for (axis, error) in certificate
            .original_derivative_error_mm()
            .into_iter()
            .enumerate()
        {
            result.derivative[axis] =
                result.derivative[axis].max((I::point(error) + I::point(derivative[axis])).hi);
        }
    }
    // Include the independently certified ORIGINAL incident edge laws in
    // every proposal, before owner matching chooses it. Their defining poles
    // were copied into this source net and exact reader boundary binding earns
    // the same actual decimal polynomial; global delta dominates edge delta.
    for (side, edge) in source.source_patch().coedges().into_iter().enumerate() {
        work.charge(64)?;
        let curve = plan
            .curves()
            .curves()
            .iter()
            .find(|curve| curve.source_curve() == edge.curve)
            .expect("original incident edge");
        let axis = usize::from(side % 2 == 0);
        for span in curve.spans() {
            work.charge(32)?;
            result.position = result
                .position
                .max((I::point(span.position_error_mm()) + I::point(delta)).hi);
            result.derivative[axis] = result.derivative[axis].max(
                (I::point(span.original_derivative_error_mm()) + I::point(derivative[axis])).hi,
            );
        }
    }
    if result.position > options.position_error_mm
        || result
            .derivative
            .iter()
            .any(|e| *e > options.original_derivative_error_mm)
    {
        return Ok(None);
    }
    Ok(Some(result))
}
fn assign(
    matches: &[[Option<Match>; 8]; 8],
    row: usize,
    used: u16,
    matching: &mut [usize; 8],
    plan: &SurfacePlan<'_>,
    owned: &OwnedMaterial,
    work: &mut Work<'_>,
) -> Result<Option<topology::Topology>, Error> {
    work.charge(1)?;
    if row == 8 {
        return topology::qualify(plan, owned, *matching, work);
    }
    for (column, value) in matches[row].iter().enumerate() {
        work.charge(1)?;
        if value.is_none() || used & (1 << column) != 0 {
            continue;
        }
        matching[row] = column;
        if let Some(topology) = assign(
            matches,
            row + 1,
            used | (1 << column),
            matching,
            plan,
            owned,
            work,
        )? {
            return Ok(Some(topology));
        }
    }
    Ok(None)
}
fn edges(
    plan: &SurfacePlan<'_>,
    owned: &OwnedMaterial,
    topology: &topology::Topology,
    options: CorrespondenceOptions,
    work: &mut Work<'_>,
) -> Result<(f64, f64), Error> {
    let (mut position, mut derivative) = (0_f64, 0_f64);
    for (index, target) in topology.edges.into_iter().enumerate() {
        let source = &plan.curves().curves()[index];
        let target = &owned.curves()[owned.edges()[target].curve_index()];
        // Exact rectangle-side bindings and matched tensor dimensions earn
        // this span equality, before comparing original curve parameter jets.
        assert_eq!(
            source.spans().len(),
            target.spans(),
            "earned boundary span identity"
        );
        let mut delta = 0_f64;
        for (span_index, span) in source.spans().iter().enumerate() {
            for (i, pole) in span.poles_mm().iter().enumerate() {
                work.charge(64)?;
                delta = delta.max(scalar_delta(
                    *pole,
                    owned
                        .point(target.point_indices()[11 * span_index + i])
                        .expect("retained curve pole"),
                ));
            }
        }
        let d = (I::point(22. * target.spans() as f64) * I::point(delta)).hi;
        for span in source.spans() {
            work.charge(32)?;
            position = position.max((I::point(span.position_error_mm()) + I::point(delta)).hi);
            derivative =
                derivative.max((I::point(span.original_derivative_error_mm()) + I::point(d)).hi);
        }
    }
    // Every accepted incident-patch proposal already included these original
    // curve bounds with a pole delta at least as large. Owner matching and
    // exact literal boundary binding preserve the same parameter/spans.
    assert!(
        position <= options.position_error_mm && derivative <= options.original_derivative_error_mm,
        "earned incident-edge tolerance dominance"
    );
    Ok((position, derivative))
}
impl SurfacePlan<'_> {
    /// The original source construction and independent owned-file construction
    /// are inherited ONCE; query work is also returned separately for a larger
    /// encode/readback pipeline that already paid those prerequisites.
    ///
    /// Edge/vertex owner bijection pairs each continuous original cap rim with
    /// the actual rim and preserves its Start/End station and orientation.
    /// Material projection independently earns the actual planar cap disk.
    /// Pairing both disks by radial interpolation from their common XY axis
    /// anchor multiplies XY rim error by t in[0,1], while the plane-height error
    /// stays constant: the same L1 C0 bound covers every cap point. Both exact
    /// horizontal planes use the same affine XY tangents, so their derivative
    /// error is zero; the reported C1 bounds cover original lateral parameters
    /// and separately reported original edge parameters. No world-pose proof.
    pub fn correspond_literal(
        &self,
        owned: &OwnedMaterial,
        options: CorrespondenceOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<CorrespondenceReceipt, Error> {
        if !options.position_error_mm.is_finite()
            || options.position_error_mm <= 0.
            || !options.original_derivative_error_mm.is_finite()
            || options.original_derivative_error_mm <= 0.
            || options.maximum_work == 0
            || options.maximum_work > 96_000_000
            || options.maximum_bytes == 0
            || options.maximum_bytes > 64 * 1024 * 1024
        {
            return Err(Error::Invalid);
        }
        let inherited = self.work_used() + owned.work_used();
        let remaining = options
            .maximum_work
            .checked_sub(inherited)
            .filter(|n| *n > 0)
            .ok_or(Error::Limit)?;
        let mut work = Work::new(remaining.min(16_000_000), interrupt)?;
        work.charge(1)?;
        let bytes =
            self.retained_budget_bytes() + owned.retained_budget_bytes() + CORRESPONDENCE_SCRATCH;
        if bytes > options.maximum_bytes {
            return Err(Error::Limit);
        }
        assert_eq!(self.patches().len(), 8, "sealed source atlas");
        assert_eq!(owned.tensors().len(), 8, "sealed literal atlas");
        let mut matches = [[None; 8]; 8];
        for (row, source) in self.patches().iter().enumerate() {
            for (column, target) in owned.tensors().iter().enumerate() {
                matches[row][column] = compare(self, source, target, owned, options, &mut work)?;
            }
        }
        let mut matching = [usize::MAX; 8];
        let topology = assign(&matches, 0, 0, &mut matching, self, owned, &mut work)?
            .ok_or(Error::Unresolved)?;
        let (edge_position, edge_derivative) = edges(self, owned, &topology, options, &mut work)?;
        let mut position = 0_f64;
        let mut derivative = [0_f64; 2];
        for (row, column) in matching.into_iter().enumerate() {
            work.charge(16)?;
            let value = matches[row][column].expect("earned matching");
            position = position.max(value.position);
            for (result, value) in derivative.iter_mut().zip(value.derivative) {
                *result = result.max(value);
            }
        }
        work.charge(1)?;
        Ok(CorrespondenceReceipt {
            matching,
            edges: topology.edges,
            vertices: topology.vertices,
            position: position.max(edge_position),
            derivative,
            edge_position,
            edge_derivative,
            query_work: work.used,
            aggregate_work: inherited + work.used,
            bytes,
        })
    }
}
#[cfg(test)]
mod tests;
