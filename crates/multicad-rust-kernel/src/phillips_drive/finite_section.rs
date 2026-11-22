//! Sealed wires of the authored horizontal section. Evaluations enclose shared
//! source expressions; these are not literal floating support intersections or
//! a manufacturing/foreign-file incidence certificate.
use super::*;
use crate::tapered_thread::interval::J;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VertexRole {
    Axis,
    LowerEnvelope,
    LowerJoin,
    InnerTip,
    UpperJoin,
    UpperEnvelope,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EdgeRole {
    Envelope,
    OuterLower,
    InnerLower,
    InnerUpper,
    OuterUpper,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SectionCurve {
    /// Centered on the authored cutter axis in section coordinates.
    CircleArc {
        radius_mm_bounds: Bounds,
        start_angle_rad_bounds: Bounds,
        sweep_rad_bounds: Bounds,
    },
    /// In the exact quarter-rotated local frame. The horizontal semiaxis of
    /// the true inclined cylinder is R0/cos(beta); the vertical one is R0.
    EllipseArc {
        center_mm_bounds: Point2Bounds,
        semiaxes_mm_bounds: [Bounds; 2],
        quarter_turn: u8,
        start_angle_rad_bounds: Bounds,
        sweep_rad_bounds: Bounds,
    },
    /// Exact source equation n.x*x+n.y*y=offset, evaluated with bounds.
    Line {
        normal_bounds: [Bounds; 2],
        offset_mm_bounds: Bounds,
    },
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SectionVertex {
    point: Point2Bounds,
    quarter: u8,
    role: VertexRole,
}
impl SectionVertex {
    pub fn point_mm_bounds(&self) -> Point2Bounds {
        self.point
    }
    pub fn source_identity(&self) -> (u8, VertexRole) {
        (self.quarter, self.role)
    }
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SectionEdge {
    vertices: [usize; 2],
    quarter: u8,
    role: EdgeRole,
    curve: SectionCurve,
}
impl SectionEdge {
    pub fn vertices(&self) -> [usize; 2] {
        self.vertices
    }
    pub fn source_identity(&self) -> (u8, EdgeRole) {
        (self.quarter, self.role)
    }
    pub fn curve(&self) -> SectionCurve {
        self.curve
    }
}
/// Shared IDs denote the same authored endpoint expression. Bounds do not
/// assert that arbitrary independently rounded support data intersect exactly.
#[derive(Clone, Debug)]
pub struct SectionBoundary<'a> {
    source: &'a Cutter,
    depth: f64,
    envelope_role: EnvelopeRole,
    vertices: Vec<SectionVertex>,
    edges: Vec<SectionEdge>,
}
impl<'a> SectionBoundary<'a> {
    pub(crate) fn retained_heap_bytes(&self) -> usize {
        self.vertices.capacity() * std::mem::size_of::<SectionVertex>()
            + self.edges.capacity() * std::mem::size_of::<SectionEdge>()
    }

    pub fn source(&self) -> &'a Cutter {
        self.source
    }
    pub fn depth_mm(&self) -> f64 {
        self.depth
    }
    pub fn envelope_role(&self) -> EnvelopeRole {
        self.envelope_role
    }
    pub fn vertices(&self) -> &[SectionVertex] {
        &self.vertices
    }
    pub fn edges(&self) -> &[SectionEdge] {
        &self.edges
    }
}
fn c(v: I) -> J {
    J {
        v,
        ..J::constant(0.)
    }
}
fn pi() -> I {
    I::new(PI.next_down(), PI.next_up())
}
fn checkpoint(cancel: &dyn Fn() -> bool) -> Result<(), Error> {
    if cancel() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn line(slope: I, offset: I, lower: bool, quarter: u8) -> SectionCurve {
    let normal = rotate([I::point(1.), if lower { slope } else { -slope }], quarter);
    SectionCurve::Line {
        normal_bounds: point2(normal),
        offset_mm_bounds: bounds(offset),
    }
}
// Certify the numerical chart of a source endpoint independently of its
// carrier. A positive unresolved/tangent family must not become a fake edge.
fn intersection_angle(x: I, y: I) -> Result<I, Error> {
    if x.lo <= 0. || y.lo <= 0. || !x.valid() || !y.valid() {
        return Err(Error::Unresolved);
    }
    let phi = (y / x).atan();
    if phi.lo <= 0. || phi.hi >= (pi() / I::point(4.)).lo {
        return Err(Error::Unresolved);
    }
    Ok(phi)
}
pub(super) fn construct<'a>(
    source: &'a Cutter,
    depth: f64,
    cancel: &dyn Fn() -> bool,
) -> Result<SectionBoundary<'a>, Error> {
    checkpoint(cancel)?;
    if !depth.is_finite() || depth < 0. || depth > source.envelope.tip.hi {
        return Err(Error::Invalid);
    }
    let eligibility = volume::Sections::new(source)?;
    let e = &source.envelope;
    let k = &source.corners[0];
    let branch = if depth == 0. || depth <= e.join.lo {
        0
    } else if depth >= e.join.hi && depth <= e.flank.lo {
        1
    } else if depth >= e.flank.hi && depth < e.tip.lo {
        2
    } else {
        return Err(Error::Unresolved);
    };
    let radius = if depth == 0. {
        e.entry
    } else {
        eligibility
            .radius(J::d(I::point(depth)), branch)
            .ok_or(Error::Unresolved)?
            .v
    };
    let envelope_role = match branch {
        0 => EnvelopeRole::EntryRoundover,
        1 => EnvelopeRole::OuterFlank,
        _ => EnvelopeRole::Floor,
    };
    // The mouth expression is admitted by Cutter construction; radius()
    // checks finite strictly positive values before returning any other depth.
    let rb = k.b_half + (k.b_depth - I::point(depth)) * k.sin_beta / k.cos_beta;
    checkpoint(cancel)?;
    let mut result = SectionBoundary {
        source,
        depth,
        envelope_role,
        vertices: Vec::new(),
        edges: Vec::new(),
    };
    if radius.hi < rb.lo {
        for q in 0..4 {
            checkpoint(cancel)?;
            result.vertices.push(SectionVertex {
                point: point2(rotate([radius, I::point(0.)], q)),
                quarter: q,
                role: VertexRole::Axis,
            });
            result.edges.push(SectionEdge {
                vertices: [q as usize, (q as usize + 1) % 4],
                quarter: q,
                role: EdgeRole::Envelope,
                curve: SectionCurve::CircleArc {
                    radius_mm_bounds: bounds(radius),
                    start_angle_rad_bounds: bounds(I::point(q as f64) * pi() / I::point(2.)),
                    sweep_rad_bounds: bounds(pi() / I::point(2.)),
                },
            });
        }
        return Ok(result);
    }
    if radius.lo <= rb.hi {
        return Err(Error::Unresolved);
    }
    let outer_slope = k.cos_theta / k.sin_theta / k.cos_beta;
    let outer_offset = rb + (k.join_s - k.join_y * k.cos_theta / k.sin_theta) / k.cos_beta;
    let join_x = rb + k.join_s / k.cos_beta;
    let circle_sq = radius.square();
    let join_sq = join_x.square() + k.join_y.square();
    let outer = if circle_sq.lo > join_sq.hi {
        true
    } else if circle_sq.hi < join_sq.lo {
        false
    } else {
        return Err(Error::Unresolved);
    };
    let [x, y] = if outer {
        volume::line_intersection(c(radius), c(outer_offset), c(outer_slope))
    } else {
        match k.inner {
            SourceInner::ShallowV { sine, cosine } => {
                volume::line_intersection(c(radius), c(rb), c(cosine / sine / k.cos_beta))
            }
            SourceInner::Rounded { radius_mm } => {
                volume::rounded_intersection(c(radius), c(rb), radius_mm, k)
            }
        }
    }
    .ok_or(Error::Unresolved)?
    .map(|v| v.v);
    let phi = intersection_angle(x, y)?;
    let mut points = vec![([x, -y], VertexRole::LowerEnvelope)];
    if outer {
        points.push(([join_x, -k.join_y], VertexRole::LowerJoin))
    }
    points.push(([rb, I::point(0.)], VertexRole::InnerTip));
    if outer {
        points.push(([join_x, k.join_y], VertexRole::UpperJoin))
    }
    points.push(([x, y], VertexRole::UpperEnvelope));
    let inner_y = if outer { k.join_y } else { y };
    checkpoint(cancel)?;
    for q in 0..4 {
        checkpoint(cancel)?;
        let base = result.vertices.len();
        for (p, role) in &points {
            result.vertices.push(SectionVertex {
                point: point2(rotate(*p, q)),
                quarter: q,
                role: *role,
            })
        }
        let count = points.len();
        for i in 0..count - 1 {
            let lower = i < (count - 1) / 2;
            let is_outer = outer && (i == 0 || i == count - 2);
            let role = match (is_outer, lower) {
                (true, true) => EdgeRole::OuterLower,
                (true, false) => EdgeRole::OuterUpper,
                (false, true) => EdgeRole::InnerLower,
                (false, false) => EdgeRole::InnerUpper,
            };
            let curve = if is_outer {
                line(outer_slope, outer_offset, lower, q)
            } else {
                match k.inner {
                    SourceInner::ShallowV { sine, cosine } => {
                        line(cosine / sine / k.cos_beta, rb, lower, q)
                    }
                    SourceInner::Rounded { radius_mm } => {
                        let alpha = (inner_y
                            / (I::point(radius_mm).square() - inner_y.square()).sqrt())
                        .atan();
                        if !alpha.valid() || alpha.lo <= 0. || alpha.hi >= (pi() / I::point(2.)).lo
                        {
                            return Err(Error::Unresolved);
                        }
                        SectionCurve::EllipseArc {
                            center_mm_bounds: point2(rotate(
                                [rb + I::point(radius_mm) / k.cos_beta, I::point(0.)],
                                q,
                            )),
                            semiaxes_mm_bounds: [
                                bounds(I::point(radius_mm) / k.cos_beta),
                                bounds(I::point(radius_mm)),
                            ],
                            quarter_turn: q,
                            start_angle_rad_bounds: bounds(if lower { pi() + alpha } else { pi() }),
                            sweep_rad_bounds: bounds(-alpha),
                        }
                    }
                }
            };
            result.edges.push(SectionEdge {
                vertices: [base + i, base + i + 1],
                quarter: q,
                role,
                curve,
            });
        }
        result.edges.push(SectionEdge {
            vertices: [base + count - 1, ((q as usize + 1) % 4) * count],
            quarter: q,
            role: EdgeRole::Envelope,
            curve: SectionCurve::CircleArc {
                radius_mm_bounds: bounds(radius),
                start_angle_rad_bounds: bounds(phi + I::point(q as f64) * pi() / I::point(2.)),
                sweep_rad_bounds: bounds(pi() / I::point(2.) - I::point(2.) * phi),
            },
        });
    }
    Ok(result)
}
#[cfg(test)]
mod tests;
