//! Actual raw analytic SUPPORT maps in the repaired source charts. This leaf
//! does not identify finite raw face fillings, occupancy, or raw material mass.
//! t is the repaired line/normalized-ray coordinate; u is a complete angular
//! revolution with theta=2*pi*u. Raw curved supports use the SAME normalized
//! ray, not a uniform-angle t. Both partial derivative bounds include this map.
use super::*;
use crate::head_profile::{Curve, Segment};
use supports::{Kind, RawSupport, dot, norm, phase, vector};
type V = [I; 3];
#[derive(Clone, Copy, Debug)]
pub struct RawHeadFaceCorrespondence {
    pub material_face_entity: u32,
    pub support_entity: u32,
    pub position_mm: f64,
    pub normalized_meridian_derivative_mm: f64,
    pub normalized_revolution_derivative_mm: f64,
}
/// Source-bound support comparison only. Original finite trim matching is a
/// separate prerequisite. No manufacturing/provider capability is conferred.
pub struct RawHeadCorrespondence<'s, 'h, 'q, 'r, 'a, 'input> {
    head: &'s ClosedRepairedHead<'h, 'q, 'r, 'a, 'input>,
    faces: [Option<RawHeadFaceCorrespondence>; 6],
    work: usize,
}
impl<'s, 'h, 'q, 'r, 'a, 'input> RawHeadCorrespondence<'s, 'h, 'q, 'r, 'a, 'input> {
    pub fn head(&self) -> &'s ClosedRepairedHead<'h, 'q, 'r, 'a, 'input> {
        self.head
    }
    pub fn faces(&self) -> impl Iterator<Item = RawHeadFaceCorrespondence> + '_ {
        self.faces.iter().flatten().copied()
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.head
            .retained_budget_bytes()
            .saturating_add(std::mem::size_of::<Self>())
    }
}
fn sub(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i] - b[i])
}
fn mul(a: V, b: I) -> V {
    a.map(|v| v * b)
}
fn point(c: CommonMouthChart, z: f64) -> V {
    std::array::from_fn(|i| {
        I::point(c.columns[0][i])
            + I::point(c.columns[3][i]) * (I::point(z) - I::point(c.outer_top_z))
    })
}
fn sum(values: &[f64]) -> Result<f64> {
    let v = values.iter().fold(I::point(0.), |a, &b| a + I::point(b));
    require(v.valid() && v.hi.is_finite() && v.lo >= 0.)?;
    Ok(v.hi)
}
fn revolution(v: f64) -> f64 {
    // Adjacent encloses real pi, unlike treating cached binary TAU as exact.
    (support::adjacent(std::f64::consts::PI) * I::point(2.) * I::point(v)).hi
}
fn basis(raw: &RawSupport, c: CommonMouthChart, sign: f64, axial: f64) -> [V; 3] {
    let e = phase(vector(&raw.frame, c.columns[1]));
    [
        std::array::from_fn(|i| raw.frame.axes[0][i] * e[0] + raw.frame.axes[1][i] * e[1]),
        std::array::from_fn(|i| {
            I::point(sign) * (-raw.frame.axes[0][i] * e[1] + raw.frame.axes[1][i] * e[0])
        }),
        mul(raw.frame.axes[2], I::point(axial)),
    ]
}
fn plane(raw: &RawSupport, c: CommonMouthChart, s: Segment) -> Result<[f64; 3]> {
    // Orthogonal projection onto the actual raw plane, retaining all in-plane
    // source coordinates. Residual is n*(n dot(P-O)); n is a true unit vector.
    let n = raw.frame.axes[2];
    let ends = [s.start, s.end].map(|p| dot(n, sub(point(c, p[1]), raw.frame.origin)));
    let radial = [c.columns[1], c.columns[2]].map(|a| dot(n, a.map(I::point)));
    let radius = s.start[0].abs().max(s.end[0].abs());
    let tilt = sum(&radial.map(I::magnitude))?;
    let dr = (I::point(s.end[0]) - I::point(s.start[0])).magnitude();
    Ok([
        sum(&[
            ends[0].magnitude().max(ends[1].magnitude()),
            (I::point(radius) * I::point(tilt)).hi,
        ])?,
        sum(&[
            (ends[1] - ends[0]).magnitude(),
            (I::point(dr) * I::point(tilt)).hi,
        ])?,
        revolution((I::point(radius) * I::point(tilt)).hi),
    ])
}
fn tube(
    raw: &RawSupport,
    c: CommonMouthChart,
    s: Segment,
    r: I,
    k: I,
    sign: f64,
) -> Result<[f64; 3]> {
    let [x, y, n] = basis(raw, c, sign, 1.);
    let a = s.start;
    let b = s.end;
    let center = [a, b].map(|p| sub(point(c, p[1]), raw.frame.origin));
    let z = center.map(|p| dot(n, p));
    let radii = z.map(|z| r + k * z);
    require(radii.into_iter().all(|r| r.valid() && r.lo > 0.))?;
    let d = std::array::from_fn::<_, 2, _>(|j| sub(center[j], mul(n, z[j])));
    let coeff = std::array::from_fn::<_, 2, _>(|j| {
        let rho = I::point([a, b][j][0]);
        [
            sub(mul(c.columns[1].map(I::point), rho), mul(x, radii[j])),
            sub(mul(c.columns[2].map(I::point), rho), mul(y, radii[j])),
        ]
    });
    let ddt = sub(
        mul(c.columns[3].map(I::point), I::point(b[1]) - I::point(a[1])),
        mul(
            n,
            dot(
                n,
                mul(c.columns[3].map(I::point), I::point(b[1]) - I::point(a[1])),
            ),
        ),
    );
    let dr = I::point(b[0]) - I::point(a[0]);
    let dz = dot(
        n,
        mul(c.columns[3].map(I::point), I::point(b[1]) - I::point(a[1])),
    );
    let radial_d = [x, y];
    let derivatives = std::array::from_fn::<_, 2, _>(|j| {
        sub(
            mul(c.columns[j + 1].map(I::point), dr),
            mul(radial_d[j], k * dz),
        )
    });
    let radial = sum(&[
        norm(coeff[0][0])?.max(norm(coeff[1][0])?),
        norm(coeff[0][1])?.max(norm(coeff[1][1])?),
    ])?;
    Ok([
        sum(&[norm(d[0])?.max(norm(d[1])?), radial])?,
        sum(&[norm(ddt)?, norm(derivatives[0])?, norm(derivatives[1])?])?,
        revolution(radial),
    ])
}
fn curved(
    raw: &RawSupport,
    c: CommonMouthChart,
    center: [f64; 2],
    radius: f64,
    envelope: [f64; 3],
    radii: [I; 2],
    orientation: [f64; 2],
) -> Result<[f64; 3]> {
    let [x, y, n] = basis(raw, c, orientation[0], orientation[1]);
    let [major, minor] = radii;
    let [speed, repair, repair_dt] = envelope;
    let u = c.columns[1].map(I::point);
    let v = c.columns[2].map(I::point);
    let z = c.columns[3].map(I::point);
    let d = sub(point(c, center[1]), raw.frame.origin);
    let a = sub(mul(u, I::point(center[0])), mul(x, major));
    let b = sub(mul(v, I::point(center[0])), mul(y, major));
    let ar = sub(mul(u, I::point(radius)), mul(x, minor));
    let br = sub(mul(v, I::point(radius)), mul(y, minor));
    let az = sub(mul(z, I::point(radius)), mul(n, minor));
    let major_error = sum(&[norm(a)?, norm(b)?])?;
    let radial_error = sum(&[norm(ar)?, norm(br)?])?;
    let minor_error = sum(&[radial_error, norm(az)?])?;
    // ||affine map|| <= Frobenius; no source frame is silently normalized.
    let stretch = norm([I::point(norm(u)?), I::point(norm(v)?), I::point(norm(z)?)])?;
    let correction = (I::point(stretch) * I::point(repair)).hi;
    let correction_dt = (I::point(stretch) * I::point(repair_dt)).hi;
    let correction_du = (I::point(repair) * (I::point(norm(u)?) + I::point(norm(v)?))).hi;
    Ok([
        sum(&[norm(d)?, major_error, minor_error, correction])?,
        sum(&[(I::point(speed) * I::point(minor_error)).hi, correction_dt])?,
        revolution(sum(&[major_error, radial_error, correction_du])?),
    ])
}
fn compare(
    raw: &RawSupport,
    c: CommonMouthChart,
    s: Segment,
    envelope: Option<[f64; 3]>,
) -> Result<[f64; 3]> {
    let mut best = None;
    for sign in [-1., 1.] {
        for axial in [-1., 1.] {
            let b = match (raw.kind, s.curve) {
                (Kind::Plane, Curve::Line) => plane(raw, c, s)?,
                (Kind::Cylinder(r), Curve::Line) => tube(raw, c, s, r, I::point(0.), sign)?,
                (Kind::Cone(r, k), Curve::Line) => tube(raw, c, s, r, k, sign)?,
                (
                    Kind::Sphere(r),
                    Curve::Arc {
                        center, radius_mm, ..
                    },
                ) => curved(
                    raw,
                    c,
                    center,
                    radius_mm,
                    envelope.ok_or(ReadError::Unsupported)?,
                    [I::point(0.), r],
                    [sign, axial],
                )?,
                (
                    Kind::Torus(r, t),
                    Curve::Arc {
                        center, radius_mm, ..
                    },
                ) => curved(
                    raw,
                    c,
                    center,
                    radius_mm,
                    envelope.ok_or(ReadError::Unsupported)?,
                    [r, t],
                    [sign, axial],
                )?,
                _ => return Err(ReadError::Unsupported),
            };
            require(b.into_iter().all(|v| v.is_finite() && v >= 0.))?;
            // Select ONE actual parameter map by max bound. Never combine
            // separate component minima from inconsistent frame choices.
            let score = b.into_iter().fold(0_f64, f64::max);
            if best.is_none_or(|(_, old)| score < old) {
                best = Some((b, score));
            }
        }
    }
    Ok(best.unwrap().0)
}
impl<'h, 'q, 'r, 'a, 'input> ClosedRepairedHead<'h, 'q, 'r, 'a, 'input> {
    pub fn raw_support_correspondence<'s>(
        &'s self,
        options: ClosedHeadOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<RawHeadCorrespondence<'s, 'h, 'q, 'r, 'a, 'input>, PhillipsReadError>
    {
        let mut w = Work {
            remaining: options.maximum_work,
            interrupt,
            cancelled: Cell::new(false),
        };
        let result = (|| {
            w.step(1)?;
            require(
                options.maximum_distance_mm.is_finite()
                    && options.maximum_distance_mm > 0.
                    && options.maximum_distance_mm <= 1e-7
                    && options.maximum_work <= WORK
                    && options.maximum_bytes <= BYTES,
            )?;
            w.step(self.work_used())?;
            let retained = self
                .retained_budget_bytes()
                .saturating_add(std::mem::size_of::<
                    RawHeadCorrespondence<'_, '_, '_, '_, '_, '_>,
                >());
            if retained.saturating_add(SCRATCH) > options.maximum_bytes {
                return Err(ReadError::Limit);
            }
            let body = self.body();
            let shell = Shell::new(&body.entities, body.tolerance, true);
            let mut faces = [None; 6];
            let top = body.profile.top();
            for (index, slot) in faces
                .iter_mut()
                .enumerate()
                .take(self.source.segment_count() + 1)
            {
                // <=4 fixed sign/axis candidates, each <4096 elementary
                // scalar/vector actions; parser + arc envelope <4096 extra.
                w.step(20_480)?;
                let (face, s, envelope) = if index == 0 {
                    (
                        body.cap.material_face,
                        Segment {
                            start: [0., top[1]],
                            end: top,
                            curve: Curve::Line,
                            role: Role::Bearing,
                        },
                        None,
                    )
                } else {
                    (
                        body.profile_faces[index - 1],
                        body.profile.segments().nth(index - 1).unwrap(),
                        self.source.raw_support_arc_envelope(index - 1),
                    )
                };
                let (surface, _, _) = helical::face(&shell, face)?;
                let raw = supports::read(&shell, surface)?;
                let b = compare(&raw, self.chart(), s, envelope)?;
                require(b.into_iter().all(|v| v <= options.maximum_distance_mm))?;
                *slot = Some(RawHeadFaceCorrespondence {
                    material_face_entity: face,
                    support_entity: surface,
                    position_mm: b[0],
                    normalized_meridian_derivative_mm: b[1],
                    normalized_revolution_derivative_mm: b[2],
                });
            }
            w.step(1)?;
            Ok(faces)
        })();
        let faces = result.map_err(|e| read_error(e, &w))?;
        Ok(RawHeadCorrespondence {
            head: self,
            faces,
            work: options.maximum_work - w.remaining,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::step_read::phillips_cutter::physical::body::{
        common_mouth::tests::derived_thread_probe::with_head, derived_thread::actual_m3,
    };
    #[test]
    fn independent_cylinder_affine_and_radius_controls_include_normalized_revolution_factor() {
        let identity = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
        let mut c = CommonMouthChart {
            columns: [[0.; 3], identity[0], identity[1], identity[2]],
            outer_top_z: 0.,
        };
        c.columns[1][0] += 1e-9;
        c.columns[3][0] = 1e-9;
        let mut raw = RawSupport {
            kind: Kind::Cylinder(I::point(2.)),
            frame: support::Frame {
                origin: [I::point(0.); 3],
                axes: identity.map(|a| a.map(I::point)),
            },
        };
        let segment = Segment {
            start: [2., 0.],
            end: [2., 1.],
            curve: Curve::Line,
            role: Role::Shank,
        };
        let b = compare(&raw, c, segment, None).unwrap();
        assert!(b[0] >= 2e-9 && b[0] < 1e-7);
        assert!(b[1] >= 1e-9 && b[1] < 1e-7);
        assert!(b[2] >= 2. * std::f64::consts::PI * 2e-9 && b[2] < 1e-7);
        raw.kind = Kind::Cylinder(I::point(2.01));
        let wrong = compare(&raw, c, segment, None).unwrap();
        assert!(wrong[0] > 0.009 && wrong[2] > 0.009 * 2. * std::f64::consts::PI);
        assert!(matches!(
            compare(
                &raw,
                c,
                Segment {
                    curve: Curve::Arc {
                        center: [0., 0.],
                        radius_mm: 2.,
                        sweep: 1.
                    },
                    ..segment
                },
                Some([1., 0., 0.])
            ),
            Err(ReadError::Unsupported)
        ));
    }
    #[test]
    fn actual_raw_head_support_c1_covers_both_hands_fits_and_resource_frontiers() {
        for fitted in [false, true] {
            for left in [false, true] {
                let text = actual_m3(fitted, left);
                with_head(&text, &mut |head| {
                    let proof = head
                        .raw_support_correspondence(Default::default(), None)
                        .unwrap();
                    assert!(std::ptr::eq(proof.head(), head));
                    assert_eq!(proof.faces().count(), head.source().segment_count() + 1);
                    for face in proof.faces() {
                        assert!(face.position_mm > 0. && face.position_mm <= 1e-7);
                        assert!(
                            face.normalized_meridian_derivative_mm >= 0.
                                && face.normalized_meridian_derivative_mm <= 1e-7
                        );
                        assert!(
                            face.normalized_revolution_derivative_mm > 0.
                                && face.normalized_revolution_derivative_mm <= 1e-7
                        );
                        eprintln!("RAW_HEAD_C1 fitted={fitted} left={left} face={face:?}");
                    }
                    let used = proof.work_used();
                    assert!(matches!(
                        head.raw_support_correspondence(
                            ClosedHeadOptions {
                                maximum_work: used - 1,
                                ..Default::default()
                            },
                            None
                        ),
                        Err(PhillipsReadError::Limit)
                    ));
                    assert_eq!(
                        head.raw_support_correspondence(
                            ClosedHeadOptions {
                                maximum_work: used,
                                ..Default::default()
                            },
                            None
                        )
                        .unwrap()
                        .work_used(),
                        used
                    );
                    assert!(matches!(
                        head.raw_support_correspondence(
                            ClosedHeadOptions {
                                maximum_bytes: proof.retained_budget_bytes() + SCRATCH - 1,
                                ..Default::default()
                            },
                            None
                        ),
                        Err(PhillipsReadError::Limit)
                    ));
                    assert!(
                        head.raw_support_correspondence(
                            ClosedHeadOptions {
                                maximum_bytes: proof.retained_budget_bytes() + SCRATCH,
                                ..Default::default()
                            },
                            None
                        )
                        .is_ok()
                    );
                    assert!(matches!(
                        head.raw_support_correspondence(Default::default(), Some(&|| true)),
                        Err(PhillipsReadError::Cancelled)
                    ));
                    assert!(matches!(
                        head.raw_support_correspondence(
                            ClosedHeadOptions {
                                maximum_distance_mm: f64::NAN,
                                ..Default::default()
                            },
                            None
                        ),
                        Err(PhillipsReadError::Invalid)
                    ));
                    assert!(
                        head.raw_support_correspondence(Default::default(), None)
                            .is_ok()
                    );
                });
            }
        }
    }
}
