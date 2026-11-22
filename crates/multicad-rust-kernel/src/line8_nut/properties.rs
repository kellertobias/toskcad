//! Source area and first moments of all seventeen physical faces. Crown
//! intersections are integrated on their actual finite disk/roof domains.
use super::*;
use crate::authored_internal_treatment::{FaceKey, FacePropertyOptions, ThreadBandRole};
#[derive(Clone, Copy, Debug)]
pub struct FaceMeasure {
    face: u32,
    area: [f64; 2],
    moments: [[f64; 2]; 3],
    centroid: [[f64; 2]; 3],
}
impl FaceMeasure {
    pub fn face(&self) -> u32 {
        self.face
    }
    pub fn area_mm2(&self) -> [f64; 2] {
        self.area
    }
    pub fn first_moments_mm3(&self) -> [[f64; 2]; 3] {
        self.moments
    }
    pub fn centroid_mm(&self) -> [[f64; 2]; 3] {
        self.centroid
    }
}
fn pair(v: I) -> [f64; 2] {
    [v.lo, v.hi]
}
fn integral_error(e: crate::authored_internal_treatment::Error) -> Error {
    match e {
        crate::authored_internal_treatment::Error::Invalid => Error::Invalid,
        crate::authored_internal_treatment::Error::Limit => Error::Limit,
        crate::authored_internal_treatment::Error::Cancelled => Error::Cancelled,
        crate::authored_internal_treatment::Error::Unresolved => Error::Unresolved,
    }
}
impl Nut {
    pub fn face_measure(
        &self,
        face: u32,
        options: FacePropertyOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<FaceMeasure, Error> {
        checkpoint(interrupt)?;
        if [
            options.area_tolerance_mm2,
            options.first_moment_tolerance_mm3,
            options.centroid_tolerance_mm,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.)
            || options.maximum_cells == 0
            || options.maximum_cells > 1_000_000
        {
            return Err(Error::Invalid);
        }
        let pi = I::new(
            std::f64::consts::PI.next_down(),
            std::f64::consts::PI.next_up(),
        );
        let length = point(self.stock.height_mm());
        let left = point(-self.bore_from_start_mm);
        let right = length - point(self.bore_from_start_mm);
        let mid = (left + right) / point(2.);
        let top = point(self.height[3]);
        let join = point(self.relief.start_mm);
        let r = point(self.relief.radius_mm);
        let a = point(self.width[3]);
        let d = (r * r - a * a).sqrt();
        let alpha = (a / d).atan();
        let j0 = pi * r * r / point(4.) - (a * d + r * r * alpha) / point(2.);
        let j1 = d * d * d / point(3.);
        let k = (top - point(self.height[2])) / (point(self.width[2]) - a);
        let metric = (point(1.) + k * k).sqrt();
        let original = |side: usize| {
            let a = self.stock.base()[0][side];
            let b = self.stock.base()[0][(side + 1) % 8];
            let dy = point(b[0]) - point(a[0]);
            let dz = point(b[1]) - point(a[1]);
            let area = length * (dy * dy + dz * dz).sqrt();
            [
                area,
                area * mid,
                area * (point(a[0]) + point(b[0])) / point(2.),
                area * (point(a[1]) + point(b[1])) / point(2.),
            ]
        };
        let source = |key: FaceKey| -> Result<[I; 4], Error> {
            let request = FacePropertyOptions {
                area_tolerance_mm2: options.area_tolerance_mm2 / 8.,
                first_moment_tolerance_mm3: options.first_moment_tolerance_mm3 / 8.,
                centroid_tolerance_mm: options.centroid_tolerance_mm / 8.,
                ..options
            };
            let p = self
                .thread
                .local_prefix_face_properties(self.relief.start_mm, key, request, interrupt)
                .map_err(integral_error)?;
            let m = p.first_moment_mm3.map(|v| I::new(v.lower, v.upper));
            Ok([I::new(p.area_mm2.lower, p.area_mm2.upper), m[0], m[1], m[2]])
        };
        let m = self.phase_rotation.map(|r| r.map(point));
        let cap = |key| -> Result<[I; 4], Error> {
            let p = source(key)?;
            Ok([
                p[0] * self.determinant,
                (m[0][0] * p[1] + m[0][1] * p[2]) * self.determinant,
                (m[1][0] * p[1] + m[1][1] * p[2]) * self.determinant,
                p[3] * self.determinant,
            ])
        };
        let values = match face {
            0 => {
                let hole = cap(FaceKey::StartCap)?;
                let area = point(2.) * point(self.width[0]) * length - hole[0];
                [
                    area,
                    point(2.) * point(self.width[0]) * length * mid - hole[1],
                    -hole[2],
                    point(0.),
                ]
            }
            1 | 2 | 6 | 7 => original(face as usize),
            3 | 5 => {
                let mut full = original(face as usize);
                let removed = point(2.) * metric * j0;
                full[0] = full[0] - removed;
                let y = point(2.) * metric * j1;
                full[2] = full[2] - if face == 3 { y } else { -y };
                full[3] = full[3] - point(2.) * metric * ((top + k * a) * j0 - k * j1);
                full
            }
            4 | 10 => {
                let half_strip = a * d + r * r * alpha;
                let removed_x = r * r * a - a * a * a / point(3.);
                let (area, x) = if face == 4 {
                    (
                        point(2.) * a * (-left) - half_strip,
                        -a * left * left + removed_x,
                    )
                } else {
                    (
                        point(2.) * a * right - half_strip,
                        a * right * right - removed_x,
                    )
                };
                [area, x, point(0.), area * top]
            }
            8 | 9 => {
                let mut area = point(0.);
                let mut z = point(0.);
                for i in 0..3 {
                    let w = point(self.width[i]);
                    let dw = point(self.width[i + 1]) - w;
                    let h = point(self.height[i]);
                    let dh = point(self.height[i + 1]) - h;
                    area = area + (w + point(self.width[i + 1])) * dh;
                    z = z + point(2.)
                        * dh
                        * (w * h + (w * dh + dw * h) / point(2.) + dw * dh / point(3.));
                }
                [
                    area,
                    area * if face == 8 { left } else { right },
                    point(0.),
                    z,
                ]
            }
            11..=14 => {
                let role = [
                    ThreadBandRole::Lower,
                    ThreadBandRole::Crest,
                    ThreadBandRole::Upper,
                    ThreadBandRole::Root,
                ][(face - 11) as usize];
                let p = source(FaceKey::Thread {
                    start: 0,
                    band: role,
                })?;
                let lo = (self.minimum_thread / point(self.thread.minimum_radius_mm()))
                    .lo
                    .min(self.determinant.lo);
                let hi = (self.maximum_thread / point(self.thread.maximum_radius_mm()))
                    .hi
                    .max(self.determinant.hi);
                let factor = I::new(lo, hi);
                let error = point((lo - 1.).abs().max((hi - 1.).abs()));
                let xy_error = error * point(self.maximum_thread.hi) * p[0];
                let z_error = error * join * p[0];
                let spread = |v: I, e: I| v + I::new(-e.hi, e.hi);
                [
                    p[0] * factor,
                    spread(m[0][0] * p[1] + m[0][1] * p[2], xy_error),
                    spread(m[1][0] * p[1] + m[1][1] * p[2], xy_error),
                    spread(p[3], z_error),
                ]
            }
            15 => {
                let hole = cap(FaceKey::EndCap)?;
                let area = pi * r * r - hole[0];
                [area, -hole[1], -hole[2], area * join]
            }
            16 => {
                let span = pi - point(2.) * alpha;
                let g1 = point(2.) * (point(2.) * d - a * span);
                let g2 = point(2.) * ((r * r / point(2.) + a * a) * span - point(3.) * a * d);
                let area = r * (point(2.) * pi * (top - join) - k * g1);
                let z = r / point(2.)
                    * (point(2.) * pi * (top * top - join * join) - point(2.) * top * k * g1
                        + k * k * g2);
                [area, point(0.), point(0.), z]
            }
            _ => return Err(Error::Invalid),
        };
        for v in values {
            valid(v)?;
        }
        if values[0].lo <= 0.
            || values[0].width() > options.area_tolerance_mm2
            || values[1..]
                .iter()
                .any(|v| v.width() > options.first_moment_tolerance_mm3)
        {
            return Err(Error::Unresolved);
        }
        let centroid = values[1..]
            .iter()
            .map(|v| *v / values[0])
            .collect::<Vec<_>>();
        if centroid.iter().any(|v| {
            !v.lo.is_finite() || !v.hi.is_finite() || v.width() > options.centroid_tolerance_mm
        }) {
            return Err(Error::Unresolved);
        }
        checkpoint(interrupt)?;
        Ok(FaceMeasure {
            face,
            area: pair(values[0]),
            moments: [values[1], values[2], values[3]].map(pair),
            centroid: [centroid[0], centroid[1], centroid[2]].map(pair),
        })
    }
}
