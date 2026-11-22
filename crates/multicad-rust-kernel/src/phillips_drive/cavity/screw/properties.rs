//! Measures of the retained thread, meridian supports and finite Phillips faces.
use super::*;
use crate::head_profile::Curve;
#[derive(Clone, Copy, Debug)]
pub struct ScrewFaceProperties {
    pub face: ScrewFace,
    pub area_mm2: Bounds,
    pub first_moment_mm3: [Bounds; 3],
    pub centroid_mm: [Bounds; 3],
}
type Values = [I; 4];
fn bounds(x: I) -> Bounds {
    [x.lo, x.hi]
}
fn raw(area: Bounds, moment: [Bounds; 3]) -> Values {
    [
        I::new(area[0], area[1]),
        I::new(moment[0][0], moment[0][1]),
        I::new(moment[1][0], moment[1][1]),
        I::new(moment[2][0], moment[2][1]),
    ]
}
fn checkpoint(f: &dyn Fn() -> bool) -> Result<(), Error> {
    if f() { Err(Error::Cancelled) } else { Ok(()) }
}
pub(super) fn thread_error(e: crate::helical_thread::FacePropertyError) -> Error {
    use crate::helical_thread::FacePropertyError as E;
    match e {
        E::Invalid => Error::Invalid,
        E::Unresolved => Error::Unresolved,
        E::Limit => Error::Limit,
        E::Cancelled => Error::Cancelled,
    }
}
fn shifted(mut v: Values, z: I) -> Values {
    v[3] = v[3] + z * v[0];
    v
}
fn receipt(
    face: ScrewFace,
    v: Values,
    o: FacePropertyOptions,
    transported: Option<[I; 3]>,
) -> Result<ScrewFaceProperties, Error> {
    if v.iter().any(|v| !v.valid())
        || v[0].lo <= 0.
        || v[0].width() > o.area_tolerance_mm2
        || v[1..]
            .iter()
            .any(|v| v.width() > o.first_moment_tolerance_mm3)
    {
        return Err(Error::Unresolved);
    }
    let mut centroid = [v[1], v[2], v[3]].map(|m| m / v[0]);
    if let Some(known) = transported {
        for (value, bound) in centroid.iter_mut().zip(known) {
            if !value.valid() || !bound.valid() {
                return Err(Error::Unresolved);
            }
            *value = I::new(value.lo.max(bound.lo), value.hi.min(bound.hi));
        }
    }
    if centroid
        .iter()
        .any(|v| !v.valid() || v.width() > o.centroid_tolerance_mm)
    {
        return Err(Error::Unresolved);
    }
    Ok(ScrewFaceProperties {
        face,
        area_mm2: bounds(v[0]),
        first_moment_mm3: [v[1], v[2], v[3]].map(bounds),
        centroid_mm: centroid.map(bounds),
    })
}
fn finish(
    raw: Vec<(ScrewFace, Values, Option<[I; 3]>)>,
    options: FacePropertyOptions,
    interrupt: &dyn Fn() -> bool,
) -> Result<Vec<ScrewFaceProperties>, Error> {
    let result = raw
        .into_iter()
        .map(|(face, value, centroid)| receipt(face, value, options, centroid))
        .collect::<Result<Vec<_>, _>>()?;
    checkpoint(interrupt)?;
    Ok(result)
}

fn meridian(segment: crate::head_profile::Segment) -> Result<Values, Error> {
    let p = I::new(
        std::f64::consts::PI.next_down(),
        std::f64::consts::PI.next_up(),
    );
    let [r, z] = segment.start.map(I::point);
    let (area, mz) = match segment.curve {
        Curve::Line => {
            let dr = I::point(segment.end[0]) - r;
            let dz = I::point(segment.end[1]) - z;
            let factor = I::point(2.) * p * (dr.square() + dz.square()).sqrt();
            (
                factor * (r + dr / I::point(2.)),
                factor * (r * z + (r * dz + dr * z) / I::point(2.) + dr * dz / I::point(3.)),
            )
        }
        Curve::Arc {
            center,
            radius_mm,
            sweep,
        } => {
            let [cr, cz] = center.map(I::point);
            let radius = I::point(radius_mm);
            let w = I::point(sweep);
            let x = r - cr;
            let y = z - cz;
            let length = (x.square() + y.square()).sqrt();
            if length.lo <= 0. || w.lo == 0. {
                return Err(Error::Unresolved);
            }
            // Reuse the source-angle enclosure, including actual atan2
            // evaluation allowance. The literal support radius is never fitted.
            let angle = crate::head_profile::source_volume::angle(segment, center);
            let a = angle.cos();
            let b = angle.sin();
            let mean_cos = (a * w.sin() + b * (w.cos() - I::point(1.))) / w;
            let mean_sin = (a * (I::point(1.) - w.cos()) + b * w.sin()) / w;
            let twice = I::point(2.) * w;
            let mean_product = a * b * twice.sin() / twice
                + (a.square() - b.square()) * (I::point(1.) - twice.cos()) / (I::point(4.) * w);
            let factor = I::point(2.) * p * radius * I::point(sweep.abs());
            (
                factor * (cr + radius * mean_cos),
                factor
                    * (cr * cz
                        + cr * radius * mean_sin
                        + cz * radius * mean_cos
                        + radius.square() * mean_product),
            )
        }
    };
    Ok([area, I::point(0.), I::point(0.), mz])
}
fn segment_bounds(segment: crate::head_profile::Segment, shift: I) -> [I; 2] {
    match segment.curve {
        Curve::Line => [
            I::new(
                segment.start[0].min(segment.end[0]),
                segment.start[0].max(segment.end[0]),
            ),
            I::new(segment.start[1], segment.end[1]) - shift,
        ],
        Curve::Arc {
            center,
            radius_mm,
            sweep,
        } => {
            let x = I::point(segment.start[0]) - I::point(center[0]);
            let y = I::point(segment.start[1]) - I::point(center[1]);
            let len = (x.square() + y.square()).sqrt();
            if !len.valid() || len.lo <= 0. {
                [
                    I::point(center[0]) + I::new(-radius_mm, radius_mm),
                    I::point(center[1]) + I::new(-radius_mm, radius_mm) - shift,
                ]
            } else {
                let angle = crate::head_profile::source_volume::angle(segment, center)
                    + I::new(sweep.min(0.), sweep.max(0.));
                [
                    I::point(center[0]) + I::point(radius_mm) * angle.cos(),
                    I::point(center[1]) + I::point(radius_mm) * angle.sin() - shift,
                ]
            }
        }
    }
}

impl PhillipsScrew {
    /// Directed source-coordinate box from the original finite supports.
    pub fn source_coordinate_bounds_mm(&self) -> [Bounds; 3] {
        let mut domain = self
            .thread
            .source_coordinate_bounds_mm()
            .map(|v| I::new(v[0], v[1]));
        let shift = I::point(self.spec.shank_length_mm);
        for s in self.profile.segments() {
            let values = segment_bounds(s, shift);
            let r = values[0].magnitude();
            domain[0] = I::new(domain[0].lo.min(-r), domain[0].hi.max(r));
            domain[1] = domain[0];
            domain[2] = I::new(
                domain[2].lo.min(values[1].lo),
                domain[2].hi.max(values[1].hi),
            );
        }
        let error = I::point(self.profile.rounding_error_mm());
        domain.map(|v| bounds(v + I::new(-error.hi, error.hi)))
    }
    /// Local source-support measures. The temporary cutter Mouth cap is
    /// subtracted from the true head land and never published as a screw face.
    /// The cell budget is shared conservatively between the two integrations.
    pub fn face_properties_with_cancel(
        &self,
        options: FacePropertyOptions,
        interrupt: &dyn Fn() -> bool,
    ) -> Result<Vec<ScrewFaceProperties>, Error> {
        checkpoint(interrupt)?;
        if options.maximum_cells == 0
            || options.maximum_cells > 65536
            || [
                options.area_tolerance_mm2,
                options.first_moment_tolerance_mm3,
                options.centroid_tolerance_mm,
            ]
            .into_iter()
            .any(|v| !v.is_finite() || v <= 0.)
        {
            return Err(Error::Invalid);
        }
        if options.maximum_cells < 2 {
            return Err(Error::Limit);
        }
        let head_extent = self
            .profile
            .support([1., 0., 0.])
            .max(self.thread.spec().length_mm)
            .max(self.spec.shank_length_mm.abs())
            .max(1.);
        let internal = FacePropertyOptions {
            area_tolerance_mm2: options
                .area_tolerance_mm2
                .min(options.first_moment_tolerance_mm3 / (head_extent + 1.))
                .min(options.centroid_tolerance_mm / (head_extent + 1.))
                / 8.,
            first_moment_tolerance_mm3: options
                .first_moment_tolerance_mm3
                .min(options.centroid_tolerance_mm)
                / 8.,
            centroid_tolerance_mm: options.centroid_tolerance_mm / 8.,
            maximum_cells: options.maximum_cells / 2,
        };
        let thread = self
            .thread
            .face_properties_with_cancel(
                self.visible_thread_interval_mm(),
                crate::helical_thread::FacePropertyOptions {
                    area_tolerance_mm2: internal.area_tolerance_mm2,
                    first_moment_tolerance_mm3: internal.first_moment_tolerance_mm3,
                    centroid_tolerance_mm: internal.centroid_tolerance_mm,
                    maximum_cells: internal.maximum_cells,
                },
                interrupt,
            )
            .map_err(thread_error)?;
        let cavity = self
            .cavity
            .face_properties_with_cancel(internal, interrupt)?;
        let pi = I::new(
            std::f64::consts::PI.next_down(),
            std::f64::consts::PI.next_up(),
        );
        let start = I::point(self.visible_start_mm);
        let radius = I::point(self.profile.end()[0]);
        let cap = raw(thread[0].area_mm2, thread[0].first_moment_mm3);
        let disk = pi * radius.square();
        let mut result = Vec::with_capacity(self.faces().len());
        result.push((
            ScrewFace::ThreadShoulder,
            [disk - cap[0], -cap[1], -cap[2], start * disk - cap[3]],
            None,
        ));
        for f in thread.into_iter().skip(1) {
            let role = if f.face == 1 {
                ScrewFace::ThreadTip
            } else {
                ScrewFace::ThreadWall(f.face - 2)
            };
            result.push((
                role,
                raw(f.area_mm2, f.first_moment_mm3),
                Some(f.centroid_mm.map(|b| I::new(b[0], b[1]))),
            ));
        }
        let mouth = cavity
            .iter()
            .find(|f| f.face == FaceKey::Mouth)
            .expect("sealed cavity Mouth");
        let mouth = raw(mouth.area_mm2, mouth.first_moment_mm3);
        let top = I::point(self.head_top_z());
        let disk = pi * I::point(self.profile.top()[0]).square();
        result.push((
            ScrewFace::HeadTop,
            [
                disk - mouth[0],
                -mouth[1],
                -mouth[2],
                top * (disk - mouth[0]),
            ],
            None,
        ));
        for (i, segment) in self.profile.segments().enumerate() {
            checkpoint(interrupt)?;
            let values = meridian(segment)?;
            let shift = I::point(-self.spec.shank_length_mm);
            let mut centroid = [values[1], values[2], values[3]].map(|v| v / values[0]);
            centroid[2] = centroid[2] + shift;
            result.push((
                ScrewFace::ProfileWall(i as u32),
                shifted(values, shift),
                Some(centroid),
            ));
        }
        for f in cavity {
            if f.face != FaceKey::Mouth {
                result.push((
                    ScrewFace::Cavity(f.face),
                    shifted(raw(f.area_mm2, f.first_moment_mm3), top),
                    // Exact translation preserves each earned centroid's
                    // correlation with its SAME face area/first moments.
                    Some(std::array::from_fn(|i| {
                        I::new(f.centroid_mm[i][0], f.centroid_mm[i][1])
                            + if i == 2 { top } else { I::point(0.) }
                    })),
                ));
            }
        }
        finish(result, options, interrupt)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn numerical_receipts_and_actual_arc_support_precision_refuse_unresolved_measures() {
        let options = FacePropertyOptions::default();
        let face = ScrewFace::HeadTop;
        let mut v = [I::point(1.), I::point(0.), I::point(0.), I::point(0.)];
        v[0] = I::new(1., 2.);
        assert!(receipt(face, v, options, None).is_err());
        assert!(finish(vec![(face, v, None)], options, &|| false).is_err());
        v[0] = I::point(1e-10);
        v[1] = I::new(0., 1e-8);
        assert!(receipt(face, v, options, None).is_err());
        let segment = crate::head_profile::Segment {
            start: [1., 0.],
            end: [1., 0.],
            curve: Curve::Arc {
                center: [0., 0.],
                radius_mm: 1.,
                sweep: 0.,
            },
            role: crate::head_profile::Role::Crown,
        };
        assert!(meridian(segment).is_err());
        let degenerate = crate::head_profile::Segment {
            start: [0., 0.],
            ..segment
        };
        let full = segment_bounds(degenerate, I::point(0.));
        assert!(full[0].lo <= -1. && full[0].hi >= 1. && full[1].lo <= -1. && full[1].hi >= 1.);
        use crate::helical_thread::FacePropertyError as E;
        for error in [E::Invalid, E::Limit, E::Cancelled, E::Unresolved] {
            assert!(matches!(
                thread_error(error),
                Error::Invalid | Error::Limit | Error::Cancelled | Error::Unresolved
            ));
        }
    }
    #[test]
    fn nominal_h0_long_plain_cavity_centroids_preserve_exact_translation() {
        use crate::helical_thread;
        let options = FacePropertyOptions {
            area_tolerance_mm2: 1e-7 * 2. / 16.,
            first_moment_tolerance_mm3: 1e-7 * 4. / 16.,
            centroid_tolerance_mm: 1e-7 / 16.,
            maximum_cells: 65536,
        };
        for left in [false, true] {
            let mut previous: Option<(f64, Vec<ScrewFaceProperties>)> = None;
            for plain in [1., 5.2, 16.] {
                let fixture = super::super::tests::fixture(0, plain, left);
                let thread = Thread::new(helical_thread::Spec {
                    major_radius_mm: 1.,
                    pitch_mm: 0.4,
                    length_mm: 0.8,
                    depth_mm: 0.21650635094610965,
                    crest_width_mm: 0.05,
                    ..fixture.thread().spec()
                })
                .unwrap();
                let source = PhillipsScrew::new(
                    thread,
                    *fixture.profile(),
                    fixture.spec(),
                    fixture.cavity().clone(),
                    100_000,
                    &|| false,
                )
                .unwrap();
                let before = source.thread().spec();
                let faces = source
                    .face_properties_with_cancel(options, &|| false)
                    .unwrap();
                let local = source
                    .cavity()
                    .face_properties_with_cancel(
                        FacePropertyOptions {
                            area_tolerance_mm2: 1e-10,
                            first_moment_tolerance_mm3: 1e-10,
                            centroid_tolerance_mm: 1e-10,
                            maximum_cells: 65536,
                        },
                        &|| false,
                    )
                    .unwrap();
                let top = source.head_top_z();
                for face in &faces {
                    let ScrewFace::Cavity(key) = face.face else {
                        continue;
                    };
                    let oracle = local.iter().find(|f| f.face == key).unwrap();
                    for axis in 0..3 {
                        let expected =
                            I::new(oracle.centroid_mm[axis][0], oracle.centroid_mm[axis][1])
                                + I::point(if axis == 2 { top } else { 0. });
                        let actual = I::new(face.centroid_mm[axis][0], face.centroid_mm[axis][1]);
                        assert!(actual.lo <= expected.hi && expected.lo <= actual.hi);
                        assert!(actual.width() <= options.centroid_tolerance_mm);
                    }
                    if let Some((old_top, old_faces)) = &previous {
                        let old = old_faces.iter().find(|f| f.face == face.face).unwrap();
                        let delta = I::new(face.centroid_mm[2][0], face.centroid_mm[2][1])
                            - I::new(old.centroid_mm[2][0], old.centroid_mm[2][1]);
                        let exact = I::point(top) - I::point(*old_top);
                        assert!(delta.lo <= exact.lo && exact.hi <= delta.hi);
                    }
                }
                assert_eq!(source.thread().spec(), before);
                previous = Some((top, faces));
            }
        }
    }
    #[test]
    fn transported_centroid_intersection_rejects_inconsistent_or_invalid_bounds() {
        let options = FacePropertyOptions::default();
        let values = [I::point(1.), I::point(0.), I::point(0.), I::point(0.)];
        assert!(receipt(ScrewFace::HeadTop, values, options, Some([I::point(1.); 3])).is_err());
        assert!(
            receipt(
                ScrewFace::HeadTop,
                values,
                options,
                Some([I::new(f64::NAN, 1.); 3])
            )
            .is_err()
        );
    }
    #[test]
    fn complete_retained_screw_faces_have_directed_source_measures_without_a_floor_or_mouth() {
        for index in 0..3 {
            for plain in [0., 1.] {
                let s = super::super::tests::fixture(index, plain, false);
                let before = s.thread().spec();
                assert!(
                    s.source_volume_bounds()
                        .unwrap()
                        .thread_phase_domain_error_mm3()
                        > 0.
                );
                let bbox = s.source_coordinate_bounds_mm();
                for segment in s.profile().segments() {
                    for t in [0., 0.5, 1.] {
                        let [r, z] = segment.point(t);
                        assert!(
                            bbox[0][1] >= r && bbox[2][0] <= z - plain && z - plain <= bbox[2][1]
                        );
                    }
                }
                let options = FacePropertyOptions {
                    maximum_cells: 65536,
                    ..FacePropertyOptions::default()
                };
                let p = s.face_properties_with_cancel(options, &|| false).unwrap();
                assert_eq!(p.iter().map(|p| p.face).collect::<Vec<_>>(), s.faces());
                assert!(
                    !p.iter()
                        .any(|p| p.face == ScrewFace::Cavity(FaceKey::Mouth))
                );
                assert!(p.iter().all(|p| p.area_mm2[0] > 0.
                    && p.area_mm2[1] - p.area_mm2[0] <= options.area_tolerance_mm2));
                assert_eq!(s.thread().spec(), before);
                let surface = p.iter().map(|p| p.area_mm2[0]).sum::<f64>();
                assert!(surface > 20.);
                let top = p.iter().find(|p| p.face == ScrewFace::HeadTop).unwrap();
                assert!(top.centroid_mm[0][0] <= 0. && top.centroid_mm[0][1] >= 0.);
                assert!(
                    top.centroid_mm[2][0] <= s.head_top_z()
                        && s.head_top_z() <= top.centroid_mm[2][1]
                );
                for (i, seg) in s.profile().segments().enumerate() {
                    let p = p
                        .iter()
                        .find(|p| p.face == ScrewFace::ProfileWall(i as u32))
                        .unwrap();
                    assert!((seg.area_mm2() - p.area_mm2[0]).abs() < 1e-7);
                    assert!(p.centroid_mm[0][0] <= 0. && p.centroid_mm[0][1] >= 0.);
                }
            }
        }
    }
    #[test]
    fn source_queries_cancel_refuse_budget_and_precision_then_retry_without_mutation() {
        use std::cell::Cell;
        let s = super::super::tests::fixture(1, 1., false);
        let options = FacePropertyOptions {
            maximum_cells: 65536,
            ..FacePropertyOptions::default()
        };
        assert!(matches!(
            s.face_properties_with_cancel(options, &|| true),
            Err(Error::Cancelled)
        ));
        let count = Cell::new(0);
        assert!(matches!(
            s.face_properties_with_cancel(options, &|| {
                count.set(count.get() + 1);
                count.get() > 100
            }),
            Err(Error::Cancelled)
        ));
        assert!(matches!(
            s.face_properties_with_cancel(
                FacePropertyOptions {
                    maximum_cells: 1,
                    ..options
                },
                &|| false
            ),
            Err(Error::Limit)
        ));
        assert!(matches!(
            s.face_properties_with_cancel(
                FacePropertyOptions {
                    maximum_cells: 0,
                    ..options
                },
                &|| false
            ),
            Err(Error::Invalid)
        ));
        assert!(
            s.face_properties_with_cancel(
                FacePropertyOptions {
                    area_tolerance_mm2: 1e-30,
                    ..options
                },
                &|| false
            )
            .is_err()
        );
        let actual = s.face_properties_with_cancel(options, &|| false).unwrap();
        let face = actual[0];
        assert!(matches!(
            finish(
                vec![(face.face, raw(face.area_mm2, face.first_moment_mm3), None)],
                options,
                &|| true
            ),
            Err(Error::Cancelled)
        ));
    }
}
