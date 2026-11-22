//! Intentional original-thread bite below an actual two-radius blind seat.
use super::*;
use crate::axial_profile_stack::BlindCounterboreSource;

/// Actual receiving dimensions, coupled to the factory-owned material proof.
pub struct SeatedContact {
    pub contact: Contact,
    pub bore_radius_mm: f64,
    pub bore_depth_mm: f64,
    pub seat_radius_mm: f64,
    pub seat_depth_mm: f64,
}

/// Prove contact with the actual rectangular material minus its two cylinders.
/// The source chart has bore center XY=0 and depth increasing along +Z.
/// No receiving geometry is inferred from bounds or labels.
pub fn certify_seated(
    stock: BlindCounterboreSource<'_>,
    stock_frame: Frame,
    screw: &Screw,
    screw_frame: Frame,
    req: Requirements,
    interrupt: &dyn Fn() -> bool,
) -> R<SeatedContact> {
    checkpoint(interrupt)?;
    if [
        req.minimum_wall_mm,
        req.minimum_floor_mm,
        req.minimum_engagement_mm,
        req.tip_clearance_mm,
        req.linear_tolerance_mm,
        req.minimum_feature_mm,
    ]
    .into_iter()
    .any(|v| !v.is_finite() || !(0. ..=1e7).contains(&v))
        || req.minimum_wall_mm == 0.
        || req.minimum_floor_mm == 0.
        || req.minimum_engagement_mm == 0.
        || req.minimum_feature_mm == 0.
        || !(1e-7..=1e-3).contains(&req.linear_tolerance_mm)
    {
        return Err(Error::Invalid);
    }
    if screw.thread().spec().root_radius_mm.is_some() {
        return Err(Error::Unsupported);
    }
    // This sealed wrapper is emitted only by the actual two-cylinder factory.
    // Its descriptor denotes the complete material, not a discarded tool.
    let s = stock.spec();
    let metric = scale(stock_frame)?;
    let sigma = metric.minimum_scale;
    let screw_metric = scale(screw_frame)?;
    let (origin, columns) = relative(stock_frame, screw_frame, interrupt)?;
    if columns[2][2] != 1 {
        return Err(Error::Unsupported);
    }
    let spec = screw.spec();
    let t = screw.thread();
    if &origin[2] - q(spec.shank_length_mm)? != q(s.seat_depth_mm)? {
        return Err(Error::Unsupported);
    }
    let offsets = [enclosure(&origin[0])?, enclosure(&origin[1])?];
    let eccentricity = offsets
        .into_iter()
        .fold(ip(0.), |a, b| a + ip(b.magnitude()));
    let exact_offset = origin[0].abs() + origin[1].abs();
    let core = (q(t.spec().major_radius_mm)? - q(t.spec().depth_mm)?).max(q(t.root_radius_mm())?);
    if core + &exact_offset > q(s.bore_radius_mm)?
        || spec.shank_length_mm > 0.
            && q(spec.shank_radius_mm)? + &exact_offset > q(s.bore_radius_mm)?
        || s.bore_radius_mm >= t.spec().major_radius_mm
    {
        return Err(Error::Unsupported);
    }
    let (head_inner, head_outer) = head_radii(screw)?;
    // The head ends exactly at the shoulder. Every head point at positive
    // depth is therefore in the actual seat void; above Z=0 there is no stock.
    if q(head_outer)? + &exact_offset > q(s.seat_radius_mm)? {
        return Err(Error::Unsupported);
    }
    let radius = thread_radius(screw)?;
    let protected = (ip(radius) + eccentricity).hi.max(s.seat_radius_mm);
    let wall = (ip(req.minimum_wall_mm) / sigma).hi;
    for side in [-s.min_mm[0], -s.min_mm[1], s.max_mm[0], s.max_mm[1]] {
        if q(protected)? + q(wall)? > q(side)? {
            return Err(Error::Unsupported);
        }
    }
    if (q(s.height_mm)? - q(s.bore_depth_mm)?) * q(sigma.lo)? < q(req.minimum_floor_mm)? {
        return Err(Error::Unsupported);
    }
    let gap = bounded(q(s.bore_depth_mm)? - &origin[2] - q(t.spec().length_mm)?)?;
    let n = vector(screw_frame.axes[2])?;
    if gap < Q::zero()
        || &gap * q(sigma.lo)? < q(req.tip_clearance_mm)?
        || bounded(q(t.spec().length_mm)?.pow(2) * dot(&n, &n)?)?
            < q(req.minimum_engagement_mm)?.pow(2)
    {
        return Err(Error::Unsupported);
    }
    let length = screw_frame.axes[2]
        .map(ip)
        .into_iter()
        .fold(ip(0.), |a, b| a + b.square())
        .sqrt();
    let engagement = finite(ip(t.spec().length_mm) * length)?;
    let gap = enclosure(&gap)?;
    // Source dimensions/gap are bounded by 1e7 and scale() proves both
    // positive scales near one, so these products are finite. The exact
    // gap*sigma.lo comparison above and outward gap enclosure also prove
    // required_tip <= upper; clamping the lower bound preserves order.
    let clearance = I::new(
        (gap * sigma).lo.max(req.tip_clearance_mm),
        (gap * metric.maximum_scale).hi,
    );
    // The disk inside the actual head's retained hex/circular boundary, minus
    // the actual pilot and screw shaft, supplies a complete flat bearing ring.
    let land = [
        (ip(radius.max(spec.shank_radius_mm).max(s.bore_radius_mm)) + eccentricity).hi,
        head_inner.min((ip(s.seat_radius_mm) - eccentricity).lo),
    ];
    if (ip(land[1]) - ip(land[0])).lo * sigma.lo < req.minimum_feature_mm {
        return Err(Error::Unsupported);
    }
    let z = t.spec().length_mm / 2.;
    let theta =
        ((z - t.spec().crest_width_mm / 4.) / t.helical_rate()).rem_euclid(std::f64::consts::TAU);
    let r = s.bore_radius_mm.midpoint(t.spec().major_radius_mm);
    let center = [r * theta.cos(), r * theta.sin(), z];
    let mut half = ((t.spec().major_radius_mm - s.bore_radius_mm) / 32.)
        .min(t.spec().crest_width_mm / 32.)
        .min(t.spec().length_mm / 32.)
        .min(1e-3);
    for _ in 0..12 {
        checkpoint(interrupt)?;
        let witness = center.map(|c| [(c - half).next_down(), (c + half).next_up()]);
        if t.classify_box(witness).map_err(thread_error)? == BoxClassification::Inside {
            let mut local = [ip(0.); 3];
            for i in 0..3 {
                local[i] = enclosure(&origin[i])?;
                for j in 0..3 {
                    local[i] = local[i]
                        + I::new(witness[j][0], witness[j][1]) * ip(f64::from(columns[j][i]));
                }
            }
            let radial_squared = local[0].square() + local[1].square();
            let v = ip(2. * half) * ip(2. * half) * ip(2. * half) * screw_metric.determinant;
            if radial_squared.lo > ip(s.bore_radius_mm).square().hi
                && local[2].lo > s.seat_depth_mm
                && local[2].hi < s.bore_depth_mm
                && (0..2).all(|i| local[i].lo > s.min_mm[i] && local[i].hi < s.max_mm[i])
                && v.valid()
                && v.lo > 0.
            {
                checkpoint(interrupt)?;
                return Ok(SeatedContact {
                    contact: Contact {
                        engagement_mm: [
                            engagement.lo.max(req.minimum_engagement_mm),
                            engagement.hi,
                        ],
                        tip_clearance_mm: [clearance.lo, clearance.hi],
                        thread_witness_mm: witness,
                        witness_volume_lower_mm3: v.lo,
                        source_chart_bearing_radial_land_mm: land,
                    },
                    bore_radius_mm: s.bore_radius_mm,
                    bore_depth_mm: s.bore_depth_mm,
                    seat_radius_mm: s.seat_radius_mm,
                    seat_depth_mm: s.seat_depth_mm,
                });
            }
        }
        half /= 2.;
    }
    Err(Error::Unsupported)
}

#[cfg(test)]
mod tests;
