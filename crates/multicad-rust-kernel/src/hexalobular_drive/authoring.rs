//! Explicit authored tangent construction, never inferred from foreign arcs.
use super::*;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TangentSpec {
    pub major_diameter_mm: f64,
    pub minor_diameter_mm: f64,
    pub outer_lobe_radius_mm: f64,
}

fn midpoint(p: [I; 2]) -> [f64; 2] {
    p.map(I::midpoint)
}
fn distance(p: [I; 2], value: [f64; 2]) -> I {
    p.into_iter().zip(value).fold(I::point(0.), |sum, (x, v)| {
        sum + I::point((x - I::point(v)).magnitude())
    })
}
fn rotate([x, y]: [I; 2], index: usize, root: I) -> [I; 2] {
    // Exact sixfold directions in Q(sqrt3); no rounded-angle rotation recipe.
    let half = I::point(0.5);
    let directions = [
        (I::point(1.), I::point(0.)),
        (half, root * half),
        (-half, root * half),
        (-I::point(1.), I::point(0.)),
        (-half, -root * half),
        (half, -root * half),
    ];
    let (c, s) = directions[index % 6];
    [c * x - s * y, s * x + c * y]
}

impl Drive {
    /// The inputs define the exact tangent family in Q(sqrt3). The evaluated
    /// descriptors carry an enclosure; they do not replace this source with
    /// independent rounded circles. For c=1/2-re/A, b=B/(2A) and r=re/A,
    /// `(b+ri)^2+c^2-2*c*(b+ri)*cos30=(r+ri)^2` determines ri. Each contact
    /// lies on the exact center line in ratio r/(r+ri), hence on both supports.
    pub fn from_authored_tangent(spec: TangentSpec, cone: Cone) -> Result<Self, Error> {
        let [a, b, re] = [
            spec.major_diameter_mm,
            spec.minor_diameter_mm,
            spec.outer_lobe_radius_mm,
        ];
        if ![a, b, re].into_iter().all(|x| x.is_finite() && x > 0.) || a <= b || re >= a / 2. {
            return Err(Error::Invalid);
        }
        if a > 1e7 {
            return Err(Error::Limit);
        }
        let root = I::point(3.).sqrt();
        let cos30 = root * I::point(0.5);
        let major = I::point(a);
        let b = I::point(b) / major * I::point(0.5);
        let r = I::point(re) / major;
        let c = I::point(0.5) - r;
        let denominator = I::point(2.) * (b - c * cos30 - r);
        let ri =
            (r.square() - b.square() - c.square() + I::point(2.) * b * c * cos30) / denominator;
        let d = b + ri;
        if !ri.valid() || ri.lo <= 0. || !d.valid() || d.lo <= 0. {
            return Err(Error::Unresolved);
        }
        let contact = [
            c + (d * cos30 - c) * r / (r + ri),
            d * I::point(0.5) * r / (r + ri),
        ];
        let contacts: [[I; 2]; 12] = std::array::from_fn(|i| {
            let p = if i.is_multiple_of(2) {
                rotate(contact, i / 2, root)
            } else {
                rotate([contact[0], -contact[1]], i / 2 + 1, root)
            };
            p.map(|x| x * major)
        });
        let mut error = 0_f64;
        let arcs = std::array::from_fn(|i| {
            let convex = i.is_multiple_of(2);
            let center = rotate(
                if convex {
                    [c, I::point(0.)]
                } else {
                    [d * cos30, d * I::point(0.5)]
                },
                i / 2,
                root,
            )
            .map(|x| x * major);
            let radius = if convex { I::point(re) } else { ri * major };
            let start = contacts[(i + 11) % 12];
            let end = contacts[i];
            let [center_mm, start_mm, end_mm] = [center, start, end].map(midpoint);
            let radius_mm = radius.midpoint();
            let a = (start_mm[1] - center_mm[1]).atan2(start_mm[0] - center_mm[0]);
            let b = (end_mm[1] - center_mm[1]).atan2(end_mm[0] - center_mm[0]);
            let sweep = if convex {
                (b - a).rem_euclid(TAU)
            } else {
                -(a - b).rem_euclid(TAU)
            };
            let center_error = distance(center, center_mm);
            let radius_error = I::point((radius - I::point(radius_mm)).magnitude());
            let mut phase_error = 0_f64;
            for (point, angle) in [(start, I::point(a)), (end, I::point(a) + I::point(sweep))] {
                let literal = [
                    I::point(radius_mm) * angle.cos(),
                    I::point(radius_mm) * angle.sin(),
                ];
                let exact: [I; 2] = std::array::from_fn(|axis| point[axis] - center[axis]);
                let delta = exact
                    .into_iter()
                    .zip(literal)
                    .fold(I::point(0.), |sum, (a, b)| {
                        sum + I::point((a - b).magnitude())
                    });
                // For vector displacement delta<r/2, angular error is at
                // most2*delta/r. Arc branch stability is checked below.
                let angle_error = (I::point(2.) * delta / radius).hi;
                phase_error = phase_error.max(angle_error);
            }
            let bound = center_error + radius_error + radius * I::point(phase_error);
            error = error.max(bound.hi);
            Arc {
                center_mm,
                radius_mm,
                start_mm,
                end_mm,
                start_angle_rad: a,
                sweep_rad: sweep,
            }
        });
        if !error.is_finite()
            || arcs.iter().any(|arc| {
                arc.radius_mm <= 32. * error
                    || arc.sweep_rad.abs() + 4. * error / arc.radius_mm >= PI
            })
        {
            return Err(Error::Unresolved);
        }
        let mut result = Self::from_retained(
            arcs,
            cone,
            TransportError {
                numeric_conversion_mm: 0.,
                arithmetic_mm: error,
            },
        )?;
        result.boundary_source = BoundarySource::AuthoredTangent(spec);
        Ok(result)
    }
}
