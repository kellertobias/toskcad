//! Cap UV traces derived from the SAME chosen exact phase/weight/lead owners.
//! Algebraic endpoint enclosures are not rounded literal BRep trim admission.
use super::*;
use std::cmp::Ordering;
#[derive(Debug)]
pub struct CapEndpoint {
    pub phase: EncodedScalar,
    /// Exact canonical PHYSICAL unit turn in[0,1), within this source cap.
    pub physical_turn: EncodedScalar,
    /// Exact parameter of the ORIGINAL Boundary cap curve (linear phase).
    pub original_parameter: EncodedScalar,
    /// Parameter of the encoded meridian; rounded root endpoints may require
    /// a genuine algebraic enclosure. Not an exact f64 trim parameter claim.
    pub meridian_parameter: [f64; 2],
}
pub struct CapTrace {
    pub source_key: MeridianRowKey,
    pub source_curve: CurveKey,
    pub cap_face: FaceKey,
    pub end: End,
    pub plane_z: EncodedScalar,
    pub endpoints: [CapEndpoint; 2],
    /// Exact homogeneous UV control triples, degree2 or3. A common degree
    /// elevation factor removes the otherwise introduced halves/thirds.
    controls: Vec<[EncodedScalar; 3]>,
    pub work_used: usize,
}
impl CapTrace {
    pub fn controls(&self) -> &[[EncodedScalar; 3]] {
        &self.controls
    }
    pub fn degree(&self) -> usize {
        self.controls.len() - 1
    }
}
impl EncodedField<'_> {
    /// None means no POSITIVE-LENGTH cap segment in this source cell. Point
    /// contacts remain owned by the original Boundary vertex/event registry.
    pub fn cap_trace(
        &self,
        index: usize,
        end: End,
        options: DomainOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Option<CapTrace>, Error> {
        let mut w = super::super::domain::work(options, self.bytes, interrupt)?;
        let mut memory = Memory::new(self.bytes, options.maximum_bytes)?;
        memory.add(SCRATCH)?;
        let cell = self.cells.get(index).ok_or(Error::Invalid)?;
        let row = &self.meridians[cell.meridian];
        let source = self.source.source_boundary().source();
        let (pitch, starts, left, length) = if row.component == Component::Metric {
            let s = source.metric().spec();
            (s.pitch_mm, s.starts, s.left_hand, s.length_mm)
        } else {
            let s = source.wood().spec();
            (s.pitch_mm, s.starts, s.left_hand, s.length_mm)
        };
        let pitch = binary(pitch, &mut w)?;
        let lead = pitch.mul(&N::integer(i64::from(starts)), &mut w)?;
        let station = if end == End::Start {
            N::integer(0)
        } else {
            binary(length, &mut w)?
        };
        let shifted = station.sub(
            &pitch.mul(&N::integer(i64::from(cell.source_key.total_period)), &mut w)?,
            &mut w,
        )?;
        let q = row
            .phase
            .iter()
            .map(|v| v.number(&mut w))
            .collect::<Result<Vec<_>, _>>()?;
        let weights = row
            .weights
            .iter()
            .map(|v| v.number(&mut w))
            .collect::<Result<Vec<_>, _>>()?;
        let a = q[0].div(&weights[0], &mut w)?;
        let b = q[row.degree()].div(&weights[row.degree()], &mut w)?;
        let lower = shifted.sub(&lead, &mut w)?;
        let lo = if lower.cmp(&a, &mut w)? == Ordering::Greater {
            lower
        } else {
            a.clone()
        };
        let hi = if shifted.cmp(&b, &mut w)? == Ordering::Less {
            shifted.clone()
        } else {
            b.clone()
        };
        if lo.cmp(&hi, &mut w)? != Ordering::Less {
            w.charge(1)?;
            return Ok(None);
        }
        w.charge(self.source.source_boundary().curves().len())?;
        let source_curve = self
            .source
            .source_boundary()
            .curves()
            .iter()
            .find_map(|c| match c.key {
                CurveKey::Cap {
                    component,
                    end: owner_end,
                    band,
                    ..
                } if component == row.component
                    && owner_end == end
                    && band == row.band
                    && c.faces.contains(&cell.source_key.face) =>
                {
                    Some(c.key)
                }
                _ => None,
            })
            .ok_or(Error::Invalid)?;

        let parameter = |phase: &N, w: &mut Work<'_>| -> Result<I, Error> {
            if phase.cmp(&a, w)? == Ordering::Equal {
                return Ok(p(0.));
            }
            if phase.cmp(&b, w)? == Ordering::Equal {
                return Ok(p(1.));
            }
            let value = if row.component == Component::Metric && row.band == ThreadBandRole::Root {
                let d = phase
                    .div(&pitch, w)?
                    .sub(&N::integer(1).div(&N::integer(2), w)?, w)?
                    .mul(&N::integer(8), w)?
                    .enclosure(w)?;
                w.charge(64)?;
                let x = p(3.) * d / (p(2.) + (p(4.) - p(3.) * d.square()).sqrt());
                (p(1.) + x) / p(2.)
            } else {
                phase.sub(&a, w)?.div(&b.sub(&a, w)?, w)?.enclosure(w)?
            };
            let clipped = I::new(value.lo.max(0.), value.hi.min(1.));
            // Earned domain: exact q is in[a,b]. The linear ratio lies
            // in[0,1]; root d lies in[-1,1], its radicand>0.99 and
            // denominator>2.99. N::enclosure has already propagated every
            // conversion/precision failure. This finite outward result
            // contains the actual u∈[0,1], so intersection is nonempty.
            Ok(clipped)
        };
        let p0 = parameter(&lo, &mut w)?;
        let p1 = parameter(&hi, &mut w)?;
        if p0.hi >= p1.lo {
            return Err(Error::Unresolved);
        }
        let (original_a, original_b) = original_domain(self, row, &pitch, &a, &b, &mut w)?;
        let mut endpoint = |phase: &N, parameter: I| -> Result<CapEndpoint, Error> {
            let v = shifted.sub(phase, &mut w)?.div(&lead, &mut w)?;
            let mut physical = if left {
                N::integer(1).sub(&v, &mut w)?
            } else {
                v
            };
            if physical.cmp(&N::integer(1), &mut w)? == Ordering::Equal {
                physical = N::integer(0);
            }
            let unwrapped = if row.piece == BandPiece::CrestLowPhase {
                phase.add(&pitch, &mut w)?
            } else {
                phase.clone()
            };
            let original = unwrapped
                .sub(&original_a, &mut w)?
                .div(&original_b.sub(&original_a, &mut w)?, &mut w)?;
            Ok(CapEndpoint {
                phase: EncodedScalar::retain(phase, &mut memory, &mut w)?,
                physical_turn: EncodedScalar::retain(&physical, &mut memory, &mut w)?,
                original_parameter: EncodedScalar::retain(&original, &mut memory, &mut w)?,
                meridian_parameter: [parameter.lo, parameter.hi],
            })
        };
        let endpoints = [endpoint(&lo, p0)?, endpoint(&hi, p1)?];
        let mut plane = station.clone();
        if row.component == Component::Metric {
            plane = plane
                .sub(&binary(source.metric().full_length_mm(), &mut w)?, &mut w)?
                .sub(&binary(source.spec().plain_length_mm, &mut w)?, &mut w)?;
        }
        let plane_z = EncodedScalar::retain(&plane, &mut memory, &mut w)?;
        let mut den = Vec::with_capacity(q.len());
        let mut num = Vec::with_capacity(q.len());
        for (q, wgt) in q.iter().zip(&weights) {
            den.push(lead.mul(wgt, &mut w)?);
            num.push(shifted.mul(wgt, &mut w)?.sub(q, &mut w)?);
        }
        let degree = row.degree() + 1;
        let mut controls = memory.allocate(degree + 1, &mut w)?;
        for i in 0..=degree {
            let mut u = N::integer(0);
            let mut v = N::integer(0);
            let mut weight = N::integer(0);
            if i > 0 {
                u = den[i - 1].mul(&N::integer(i as i64), &mut w)?;
                v = num[i - 1].mul(&N::integer(i as i64), &mut w)?;
                weight = u.clone();
            }
            if i < degree {
                let scale = N::integer((degree - i) as i64);
                v = v.add(&num[i].mul(&scale, &mut w)?, &mut w)?;
                weight = weight.add(&den[i].mul(&scale, &mut w)?, &mut w)?;
            }
            controls.push([
                EncodedScalar::retain(&u, &mut memory, &mut w)?,
                EncodedScalar::retain(&v, &mut memory, &mut w)?,
                EncodedScalar::retain(&weight, &mut memory, &mut w)?,
            ]);
        }
        let cap_face = match (row.component, end) {
            (Component::Metric, End::Start) => FaceKey::MetricTip,
            (Component::Metric, End::End) => FaceKey::MetricShoulder,
            (Component::Wood, End::Start) => FaceKey::WoodShoulder,
            (Component::Wood, End::End) => FaceKey::WoodTip,
        };
        w.charge(1)?;
        Ok(Some(CapTrace {
            source_key: cell.source_key,
            source_curve,
            cap_face,
            end,
            plane_z,
            endpoints,
            controls,
            work_used: w.used,
        }))
    }
}
fn original_domain(
    field: &EncodedField<'_>,
    row: &EncodedMeridian,
    pitch: &N,
    a: &N,
    b: &N,
    w: &mut Work<'_>,
) -> Result<(N, N), Error> {
    if row.band != ThreadBandRole::Crest {
        return Ok((a.clone(), b.clone()));
    }
    if row.component == Component::Wood {
        let bands = field.source.source_boundary().source().wood().bands();
        Ok((
            binary(bands[4].phase_mm[0], w)?,
            pitch.add(&binary(bands[0].phase_mm[1], w)?, w)?,
        ))
    } else {
        let c = binary(
            field
                .source
                .source_boundary()
                .source()
                .metric()
                .spec()
                .crest_width_mm,
            w,
        )?
        .div(&N::integer(2), w)?;
        Ok((pitch.sub(&c, w)?, pitch.add(&c, w)?))
    }
}
#[cfg(test)]
pub(super) mod tests;
