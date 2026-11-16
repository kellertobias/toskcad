//! Source -> actual finite encoded-field volume enclosure. No file is decoded
//! here. Exact affine axial stations and the common angular graph are required.
use super::*;
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FieldVolume {
    pub source_mm3: [f64; 2],
    pub encoded_mm3: [f64; 2],
    pub difference_bound_mm3: f64,
    pub work_used: usize,
}
fn source_error(e: crate::authored_hanger_bolt::Error) -> Error {
    use crate::authored_hanger_bolt::Error as E;
    match e {
        E::Invalid => Error::Invalid,
        E::Unsupported | E::Unresolved => Error::Unresolved,
        E::Limit => Error::Limit,
        E::Cancelled => Error::Cancelled,
    }
}
impl ClosedField<'_, '_, '_> {
    /// Actual filled finite field; includes geometry preparation in the ledger.
    /// A later serializer/decoder must earn its OWN transfer/measurement proof.
    pub fn volume_bounds(
        &self,
        error_mm3: f64,
        options: DomainOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<FieldVolume, Error> {
        let mut w = super::super::super::domain::work(options, self.bytes, interrupt)?;
        w.charge(self.work_used())?;
        if !error_mm3.is_finite() || error_mm3 <= 0. {
            return Err(Error::Invalid);
        }
        if error_mm3 / 4. == 0. {
            return Err(Error::Unresolved);
        }
        let field = self.domains.field();
        let source = field.source.source_boundary().source();
        // Eleven fixed analytic original-source cells. This precharge covers
        // their bounded interval arithmetic; original cancellation is retained.
        w.charge(11 * 512)?;
        let original = source
            .volume(
                crate::authored_hanger_bolt::QueryOptions {
                    volume_error_mm3: error_mm3 / 4.,
                    maximum_cells: 11,
                },
                interrupt,
            )
            .map_err(source_error)?;
        let ep = p(field.angular.position_error());
        let tau = I::new(
            std::f64::consts::TAU.next_down(),
            std::f64::consts::TAU.next_up(),
        );
        let pi = tau / p(2.);
        // The earned metric-root normal proof implies ep<1/(2sqrt3)<1/2:
        // Du_error>=root_r_du*ep and |Pv|>=tau*Rmax, while the normal lower
        // bound<=tau*Rmin*p/8 and root_r_du>=sqrt3*p/8. Thus the elementary
        // angle inequality |arg(B)-trueAngle|<=asin(ep)<=2ep applies, without
        // assuming a libm inverse-trig rounding bound.
        let phase_angle = p(2.) * ep / tau;
        let mut reserve = p(0.);
        for component in [Component::Metric, Component::Wood] {
            let mut radial = 0_f64;
            for m in field.meridians() {
                w.charge(1)?;
                if m.component == component {
                    radial = radial.max(m.bounds.position_error_mm);
                }
            }
            let (radius, length, pitch, n, slope) = if component == Component::Metric {
                let m = source.metric().spec();
                (
                    m.major_radius_mm,
                    m.length_mm,
                    m.pitch_mm,
                    m.starts,
                    p(3.).sqrt(),
                )
            } else {
                let m = source.wood().spec();
                let mut slope = 0_f64;
                for b in source.wood().bands() {
                    w.charge(128)?;
                    let s = ((p(b.radius_mm[1]) - p(b.radius_mm[0]))
                        / (p(b.phase_mm[1]) - p(b.phase_mm[0])))
                    .magnitude();
                    slope = slope.max(s);
                }
                (
                    m.major_radius_mm,
                    m.length_mm,
                    m.pitch_mm,
                    m.starts,
                    p(slope),
                )
            };
            // Whole-phase global Lipschitz continuity includes every actual
            // knot and wrap. Lead (NOT pitch) transfers the angle error to q.
            // Stored C0 bound conservatively includes radial-proxy error too.
            let delta = p(radial)
                + (p(radius) + p(radial)) * ep
                + slope * p(pitch) * p(f64::from(n)) * phase_angle;
            reserve = reserve + pi * p(length) * (p(2.) * p(radius) * delta + delta.square());
        }
        let plain = source.spec();
        let delta = p(plain.plain_radius_mm) * ep;
        reserve = reserve
            + pi * p(plain.plain_length_mm)
                * (p(2.) * p(plain.plain_radius_mm) * delta + delta.square());
        let actual = hull(original.material_mm3) + I::new(-reserve.hi, reserve.hi);
        if !reserve.valid() || !actual.valid() || actual.lo <= 0. || actual.width() > error_mm3 {
            return Err(Error::Unresolved);
        }
        w.charge(1)?;
        Ok(FieldVolume {
            source_mm3: original.material_mm3,
            encoded_mm3: [actual.lo, actual.hi],
            difference_bound_mm3: reserve.hi,
            work_used: w.used,
        })
    }
}
#[cfg(test)]
mod tests;
