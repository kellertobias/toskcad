//! Closed finite embedding of the actual authored expression boundary.
//! This certificate does not certify a literal foreign STEP shell or its mass.
use super::*;
#[cfg(test)]
mod tests;
const MAX_WORK: usize = 16_000_000;
const MAX_BYTES: usize = 64 * 1024 * 1024;
// Fixed four-corner source: <=128 nodes,256edges and88chart cells. The original
// topology, refined incidence, per-vertex links and simultaneous growing Vecs
// fit below2MiB even allowing fourfold payload during growth/reallocation.
const SCRATCH: usize = 2 * 1024 * 1024;
const FIXED_WORK: usize = 524288;
#[derive(Clone, Copy)]
pub struct ClosedBoundaryOptions {
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for ClosedBoundaryOptions {
    fn default() -> Self {
        Self {
            maximum_work: MAX_WORK,
            maximum_bytes: MAX_BYTES,
        }
    }
}
/// Proof-local shared semantic curve/node labels, not persisted topology names.
/// Sides are ordered u=0,u=1,v=0,v=1. Only labels from this certificate compare.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SourceChartSide {
    Curve(usize),
    Collapsed(usize),
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CertifiedPatchDomain {
    pub source_patch: PatchKey,
    pub source_face: FaceKey,
    pub winding: i8,
    pub normalized_domain: [[f64; 2]; 2],
    pub sides: [SourceChartSide; 4],
    /// At the last envelope/contact circle z_y=0, while its geometric boundary
    /// remains a genuine shared circle. It is not a collapsed source side.
    pub degenerate_contact_parameter: bool,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CertifiedTrimDomain {
    pub envelope: EnvelopeRole,
    pub wall: WallKind,
    pub events: [EventRole; 2],
    pub parameters: [Bounds; 2],
    pub y_mm: [Bounds; 2],
    pub depth_mm: [Bounds; 2],
    pub parameter_width: Bounds,
    pub y_decrease_mm: Bounds,
    pub depth_increase_mm: Bounds,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ClosedBoundaryReserves {
    pub front_intercept_mm: Bounds,
    pub front_slope_below_one: Bounds,
    pub adjacent_corner_diagonal_mm: Bounds,
    pub finite_far_mm: Bounds,
    pub finite_start_mm: Bounds,
    pub finite_end_mm: Bounds,
    pub depth_transversality: [Bounds; 2],
    pub cone_footpoints_mm: [Bounds; 2],
    pub rounded_pole_mm2: Option<Bounds>,
    pub floor_tip_depth_mm: Bounds,
}
/// Geometric outward anchor on a regular interior point of the true floor.
/// Derivatives/gradient are scaled by the positive authored b-width; this dot
/// is dimensionless, not a volume or a tolerance-derived residual.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OutwardOrientationWitness {
    pub patch_index: usize,
    pub normalized_point: [f64; 2],
    pub positive_gradient_dot: Bounds,
}
/// A closed oriented embedded sphere bounding the CURRENT Cutter expression.
/// Global proof: each horizontal section is a positive star disk. Corner fronts
/// x=f(|y|,z) have 0<=f_y<1 and f-y*f_y>0; the diagonal reserve keeps their four
/// rotated removed regions disjoint. Along every ray the front residual is
/// strictly increasing, so min(envelopeR,frontRoot) is a single continuous
/// positive radius until the unique true cone apex. Its side is injective in
/// (theta,z), with one circle quotient at the apex; the mouth polar disk closes
/// it. Finite extrusion caps/closure are proved inactive over the whole box.
///
/// Actual charts cover that boundary: trim y/depth order partitions the walls
/// and envelope, quarter angles partition their sectors; cap trace/gap sectors
/// partition the mouth. Wall charts are graphs in(y,z), envelope charts in
/// (theta,z), cap charts are polar disks and the floor is a cone fan. Shared
/// expressions, NOT approximate point matching, glue their source boundaries.
/// Legitimate quotient endpoints are retained explicitly: wall-mouth vertex,
/// cap center and cone apex. The last envelope's contact circle has a vanishing
/// depth derivative but strict whole-function monotonicity, not a false
/// everywhere-positive Jacobian. Torus z'(0)=0 is only its boundary, and its
/// radial derivative there is nonzero. No sampled Jacobian proves this result.
///
/// This is separate from raw-to-derived transfer, outer material containment,
/// manufacturing approximation, literal embedding, flux and imported-file mass.
pub struct ClosedBoundaryCertificate<'a> {
    source: &'a Cutter,
    reserves: ClosedBoundaryReserves,
    trims: [CertifiedTrimDomain; 4],
    domains: [CertifiedPatchDomain; 88],
    topology: TopologyCertificate,
    charts: TopologyCertificate,
    work: usize,
    outward: OutwardOrientationWitness,
}
impl<'a> ClosedBoundaryCertificate<'a> {
    pub(crate) fn scratch_budget_bytes() -> usize {
        SCRATCH
    }
    pub fn source(&self) -> &'a Cutter {
        self.source
    }
    pub fn reserves(&self) -> ClosedBoundaryReserves {
        self.reserves
    }
    pub fn trim_domains(&self) -> &[CertifiedTrimDomain; 4] {
        &self.trims
    }
    pub fn patch_domains(&self) -> &[CertifiedPatchDomain; 88] {
        &self.domains
    }
    pub fn topology(&self) -> TopologyCertificate {
        self.topology
    }
    pub fn chart_topology(&self) -> TopologyCertificate {
        self.charts
    }
    pub fn outward_orientation(&self) -> OutwardOrientationWitness {
        self.outward
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    /// Includes this fixed certificate and its borrowed fixed-size source once.
    pub fn retained_budget_bytes(&self) -> usize {
        std::mem::size_of::<Self>() + std::mem::size_of::<Cutter>()
    }
}
struct Work<'a> {
    remaining: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Work<'_> {
    fn checkpoint(&self) -> Result<(), Error> {
        if self.interrupt.is_some_and(|f| f()) {
            Err(Error::Cancelled)
        } else {
            Ok(())
        }
    }
    fn charge(&mut self, n: usize) -> Result<(), Error> {
        self.checkpoint()?;
        self.remaining = self.remaining.checked_sub(n).ok_or(Error::Limit)?;
        Ok(())
    }
}
fn b([lo, hi]: Bounds) -> I {
    I::new(lo, hi)
}
fn positive(values: impl IntoIterator<Item = I>) -> Result<(), Error> {
    if values.into_iter().all(|v| v.valid() && v.lo > 0.) {
        Ok(())
    } else {
        Err(Error::Unresolved)
    }
}
fn reserves(source: &Cutter) -> Result<ClosedBoundaryReserves, Error> {
    let k = &source.corners[0];
    let e = &source.envelope;
    let beta = k.sin_beta / k.cos_beta;
    let outer_slope = k.cos_theta / k.sin_theta / k.cos_beta;
    let [rb, _offset, outer] = volume::radial_intercepts(source);
    let inner_slope = match k.inner {
        SourceInner::Rounded { .. } => outer_slope,
        SourceInner::ShallowV { sine, cosine } => cosine / sine / k.cos_beta,
    };
    let front = I::new(rb.lo.min(outer.lo), rb.hi.min(outer.hi));
    let slope = I::point(1.)
        - I::new(
            inner_slope.lo.max(outer_slope.lo),
            inner_slope.hi.max(outer_slope.hi),
        );
    let diagonal = e.entry / I::point(2.).sqrt();
    let f = if diagonal.hi < k.join_y.lo {
        match k.inner {
            SourceInner::ShallowV { sine, cosine } => diagonal * cosine / sine,
            SourceInner::Rounded { radius_mm } => {
                let r = I::point(radius_mm);
                diagonal.square() / (r + (r.square() - diagonal.square()).sqrt())
            }
        }
    } else if diagonal.lo > k.join_y.hi {
        k.join_s + (diagonal - k.join_y) * k.cos_theta / k.sin_theta
    } else {
        return Err(Error::Unresolved);
    };
    let separation = rb + f / k.cos_beta - diagonal;
    let smax = (k.cos_beta * (I::new(-e.entry.hi, e.entry.hi) - k.b_half)
        + k.sin_beta * (I::new(0., e.tip.hi) - k.b_depth))
        .hi
        .max(0.);
    let far = k.far_s - I::point(smax);
    let start = -I::point(source.domain.corner_start_depth_mm) - I::point(smax) * k.sin_beta;
    let end = I::point(source.domain.corner_end_depth_mm) - e.tip;
    let slopes = [e.sin_a / e.cos_a, e.cos_c / e.sin_c];
    let intercepts = [
        I::point(source.spec.mouth_diameter_mm) / I::point(2.),
        e.tip * slopes[1],
    ];
    let trans = slopes.map(|s| s - beta);
    let base = k.b_half + k.b_depth * beta;
    let feet = std::array::from_fn::<_, 2, _>(|i| base - beta / slopes[i] * intercepts[i]);
    let pole = match k.inner {
        SourceInner::Rounded { radius_mm } => {
            Some(I::point(radius_mm).square() * k.sin_theta.square())
        }
        SourceInner::ShallowV { .. } => None,
    };
    // Exact authored contact identity: tip-b_depth=(b/2)*tan(floor).
    // Retain that coupling rather than subtracting two near-equal stations.
    let floor_tip = k.b_half * e.sin_c / e.cos_c;
    positive(
        [front, slope, separation, far, start, end, floor_tip]
            .into_iter()
            .chain(trans)
            .chain(feet)
            .chain(pole),
    )?;
    // f-y*f_y is rb for a V and >=min(rb,outer) for a rounded front,
    // including its exact tangent join; no C0 tube is being used as area proof.
    Ok(ClosedBoundaryReserves {
        front_intercept_mm: bounds(front),
        front_slope_below_one: bounds(slope),
        adjacent_corner_diagonal_mm: bounds(separation),
        finite_far_mm: bounds(far),
        finite_start_mm: bounds(start),
        finite_end_mm: bounds(end),
        depth_transversality: trans.map(bounds),
        cone_footpoints_mm: feet.map(bounds),
        rounded_pole_mm2: pole.map(bounds),
        floor_tip_depth_mm: bounds(floor_tip),
    })
}
fn outward_witness(atlas: &PatchAtlas<'_>) -> Result<OutwardOrientationWitness, Error> {
    let index = atlas
        .patches()
        .iter()
        .position(|p| p.source_identity() == PatchKey::FloorFan { quarter: 0 })
        .ok_or(Error::Unresolved)?;
    let p = atlas.patch_bounds(index, [[0.5; 2]; 2])?;
    let source = atlas.topology().events().source();
    let scale = I::point(source.spec.b_width_mm);
    let du = p.du_mm.map(|v| b(v) / scale);
    let dv = p.dv_mm.map(|v| b(v) / scale);
    let x = b(p.point_mm[0]) / scale;
    let y = b(p.point_mm[1]) / scale;
    let rho = (x.square() + y.square()).sqrt();
    // True floor implicit function F=rho+cot(c)*(z-tip); inside is F<0.
    // rho*grad(F) avoids any normalization/division by the radial coordinate.
    // This regular reference lies strictly above contact, where corners cannot
    // cut the cone by the earned depth/front ordering.
    let gradient = [x, y, (source.envelope.cos_c / source.envelope.sin_c) * rho];
    let cross: [I; 3] = std::array::from_fn(|i| {
        du[(i + 1) % 3] * dv[(i + 2) % 3] - du[(i + 2) % 3] * dv[(i + 1) % 3]
    });
    let winding = I::point(f64::from(atlas.patches()[index].winding()));
    let dot = (0..3).fold(I::point(0.), |sum, i| {
        sum + winding * cross[i] * gradient[i]
    });
    positive([dot])?;
    Ok(OutwardOrientationWitness {
        patch_index: index,
        normalized_point: [0.5; 2],
        positive_gradient_dot: bounds(dot),
    })
}
impl Cutter {
    pub fn closed_boundary_certificate(
        &self,
        options: ClosedBoundaryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<ClosedBoundaryCertificate<'_>, Error> {
        let mut work = Work {
            remaining: options.maximum_work,
            interrupt,
        };
        work.checkpoint()?;
        if options.maximum_work == 0
            || options.maximum_work > MAX_WORK
            || options.maximum_bytes == 0
            || options.maximum_bytes > MAX_BYTES
        {
            return Err(Error::Invalid);
        }
        if std::mem::size_of::<ClosedBoundaryCertificate<'_>>()
            + std::mem::size_of::<Cutter>()
            + SCRATCH
            > options.maximum_bytes
        {
            return Err(Error::Limit);
        }
        work.charge(FIXED_WORK)?;
        let reserves = reserves(self)?;
        let atlas = self.surface_patches_with_cancel(&|| work.checkpoint().is_err())?;
        let events = atlas.topology().events();
        let mut trims = [CertifiedTrimDomain {
            envelope: EnvelopeRole::Floor,
            wall: WallKind::Inner,
            events: [EventRole::FloorApex; 2],
            parameters: [[0.; 2]; 2],
            y_mm: [[0.; 2]; 2],
            depth_mm: [[0.; 2]; 2],
            parameter_width: [0.; 2],
            y_decrease_mm: [0.; 2],
            depth_increase_mm: [0.; 2],
        }; 4];
        for (i, t) in events.trims().iter().enumerate() {
            work.charge(1024)?;
            let ends = t.end_events();
            let p = ends.map(|e| events.event(e).point_mm_bounds().map(b));
            let parameters = t.parameter_bounds();
            let width = if t.envelope_role() == EnvelopeRole::EntryRoundover {
                b(parameters[1]) - b(parameters[0])
            } else {
                b(parameters[0]) - b(parameters[1])
            };
            let dy = p[0][1] - p[1][1];
            let dz = p[1][2] - p[0][2];
            positive([width, dy, dz])?;
            trims[i] = CertifiedTrimDomain {
                envelope: t.envelope_role(),
                wall: t.wall_kind(),
                events: ends,
                parameters,
                y_mm: [bounds(p[0][1]), bounds(p[1][1])],
                depth_mm: [bounds(p[0][2]), bounds(p[1][2])],
                parameter_width: bounds(width),
                y_decrease_mm: bounds(dy),
                depth_increase_mm: bounds(dz),
            };
        }
        let chart = atlas.chart_incidence(&|| work.checkpoint().is_err())?;
        // Sealed atlas construction enumerates all4trims*4quarters*5charts+8 caps.
        // Its typed face-owner table and incidence assembly have already validated
        // complete owners, opposite uses, links and shared quotient source identities.
        let domains = std::array::from_fn(|i| {
            let p = &atlas.patches()[i];
            let key = p.source_identity();
            CertifiedPatchDomain {
                source_patch: key,
                source_face: atlas.topology().faces()[p.face_index()].source_identity(),
                winding: p.winding(),
                normalized_domain: [[0., 1.]; 2],
                sides: chart.sides[i],
                degenerate_contact_parameter: matches!(
                    key,
                    PatchKey::Envelope {
                        envelope: EnvelopeRole::Floor,
                        wall: WallKind::Inner,
                        ..
                    }
                ),
            }
        });
        let topology = atlas.topology().certificate();
        work.charge(4096)?;
        let outward = outward_witness(&atlas)?;
        // This actual outward reference plus the connected, consistently
        // oriented source chart complex fixes the global orientation.
        work.charge(1)?;
        Ok(ClosedBoundaryCertificate {
            source: self,
            reserves,
            trims,
            domains,
            topology,
            charts: chart.topology,
            work: options.maximum_work - work.remaining,
            outward,
        })
    }
}
