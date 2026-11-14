//! Strict outward source envelopes. Overlapping envelopes prove nothing;
//! these bounds never replace stock material or positive overlap predicates.
use super::*;

fn world_box(local: [I; 3], frame: Frame) -> R<[I; 3]> {
    scale(frame)?;
    let bounds: [I; 3] = std::array::from_fn(|i| {
        (0..3).fold(ip(frame.origin_mm[i]), |sum, j| {
            sum + ip(frame.axes[j][i]) * local[j]
        })
    });
    Ok([finite(bounds[0])?, finite(bounds[1])?, finite(bounds[2])?])
}
fn screw_box(screw: &Screw, frame: Frame) -> R<[I; 3]> {
    let (_, head) = head_radii(screw)?;
    // Include original line coefficients and the complete retained root circle
    // bound. No phase mean, recipe or sampled extrema is a material authority.
    let radius = thread_radius(screw)?
        .max(screw.thread().source_coordinate_bounds_mm()[0][1])
        .max(head)
        .max(screw.spec().shank_radius_mm);
    let r = positive(ip(radius))?.hi;
    let top = (-ip(screw.spec().shank_length_mm) - ip(screw.spec().head.height_mm()))
        .lo
        .min(screw.head_top_z());
    world_box(
        [
            I::new(-r, r),
            I::new(-r, r),
            I::new(top, screw.thread().spec().length_mm),
        ],
        frame,
    )
}
fn separated(a: [I; 3], b: [I; 3]) -> bool {
    (0..3).any(|i| a[i].hi < b[i].lo || b[i].hi < a[i].lo)
}

/// The outer envelope also contains any declared subtracted pockets/slots.
/// A false result is unresolved, not a positive collision classification.
pub fn separated_from_stock(
    material: &OrthogonalSolid,
    material_frame: Frame,
    screw: &Screw,
    screw_frame: Frame,
    interrupt: &dyn Fn() -> bool,
) -> R<bool> {
    checkpoint(interrupt)?;
    let a = screw_box(screw, screw_frame)?;
    let bounds = material.bounds();
    let b = world_box(
        std::array::from_fn(|i| I::new(bounds.min[i], bounds.max[i])),
        material_frame,
    )?;
    let result = separated(a, b);
    checkpoint(interrupt)?;
    Ok(result)
}
/// Only the exact retained rectangular outer source is admitted here. Its
/// circular holes remove material and cannot invalidate this OUTER envelope.
/// Other curved profiles remain unknown rather than using nominal bounds.
pub fn separated_from_rectangular_prism(
    prism: &crate::prism::Prism,
    material_frame: Frame,
    screw: &Screw,
    screw_frame: Frame,
    interrupt: &dyn Fn() -> bool,
) -> R<Option<bool>> {
    checkpoint(interrupt)?;
    if prism.rounded_end().is_some() {
        return Ok(None);
    }
    let a = screw_box(screw, screw_frame)?;
    let mut local = [ip(0.); 3];
    for (i, dimension) in prism.dimensions().into_iter().enumerate() {
        let (n, d) = dimension.ratio();
        local[i] = I::new(0., positive(enclosure(&Q::new(n.into(), d.into()))?)?.hi);
    }
    let b = world_box(local, material_frame)?;
    let result = separated(a, b);
    checkpoint(interrupt)?;
    Ok(Some(result))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        Rational,
        orthogonal::{self, Box3, Operation},
        prism::Prism,
    };
    #[test]
    fn actual_hex_vertices_and_raw_affine_endpoints_are_enclosed_outward() {
        let (_, _, base, _, _) = super::super::tests::fixture(false, 1);
        for af in [10., 10.3] {
            let mut spec = base.spec();
            spec.head = ScrewHead::ExternalHex {
                across_flats_mm: af,
                height_mm: 4.,
            };
            let screw = Screw::new(*base.thread(), spec).unwrap();
            for origin in [12., 1000.3, 1000000.3] {
                let frame = Frame {
                    origin_mm: [origin, -2000.7, 20.],
                    axes: [[0.8, 0.6, 0.], [-0.6, 0.8, 0.], [0., 0., 1.]],
                };
                let bounds = screw_box(&screw, frame).unwrap();
                for [x, y] in crate::headed_screw::hex_vertices(af) {
                    for z in [screw.head_top_z(), screw.head_seat_z()] {
                        for (i, bound) in bounds.iter().enumerate() {
                            let actual = q(frame.origin_mm[i]).unwrap()
                                + q(frame.axes[0][i]).unwrap() * q(x).unwrap()
                                + q(frame.axes[1][i]).unwrap() * q(y).unwrap()
                                + q(frame.axes[2][i]).unwrap() * q(z).unwrap();
                            assert!(
                                q(bound.lo).unwrap() <= actual && q(bound.hi).unwrap() >= actual
                            );
                        }
                    }
                }
            }
        }
    }
    #[test]
    fn exact_stock_and_rational_prism_separation_refuse_contact_and_ulp_ambiguity() {
        let (_, _, screw, frame, _) = super::super::tests::fixture(false, 1);
        let material = orthogonal::boolean(
            &[Box3 {
                min: [0.; 3],
                max: [2.; 3],
            }],
            &[],
            Operation::Union,
            1e-7,
        )
        .unwrap()
        .unwrap();
        let near = Frame {
            origin_mm: [9., 9., 14.],
            axes: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        };
        assert!(!separated_from_stock(&material, near, &screw, frame, &|| false).unwrap());
        let far = Frame {
            origin_mm: [30.; 3],
            ..near
        };
        assert!(separated_from_stock(&material, far, &screw, frame, &|| false).unwrap());
        let prism = Prism::rectangular(7, [Rational::new(21, 10).unwrap(); 3]).unwrap();
        assert_eq!(
            separated_from_rectangular_prism(&prism, far, &screw, frame, &|| false).unwrap(),
            Some(true)
        );
        assert_eq!(
            separated_from_rectangular_prism(&prism, near, &screw, frame, &|| false).unwrap(),
            Some(false)
        );
        let bounds = screw_box(&screw, frame).unwrap();
        for x in [
            bounds[0].hi.next_down(),
            bounds[0].hi,
            bounds[0].hi.next_up(),
        ] {
            let touch = Frame {
                origin_mm: [x, 10., 20.],
                ..near
            };
            assert!(
                !separated_from_stock(&material, touch, &screw, frame, &|| false).unwrap(),
                "an interval-wide contact is unresolved"
            );
        }
        let clear = Frame {
            origin_mm: [bounds[0].hi + 1e-8, 10., 20.],
            ..near
        };
        assert!(separated_from_stock(&material, clear, &screw, frame, &|| false).unwrap());
        let rounded = Prism::extrude(
            7,
            Rational::new(10, 1).unwrap(),
            Rational::new(2, 1).unwrap(),
            Rational::new(2, 1).unwrap(),
        )
        .unwrap();
        assert_eq!(
            separated_from_rectangular_prism(&rounded, far, &screw, frame, &|| false).unwrap(),
            None
        );
    }
    #[test]
    fn envelope_publication_polls_cancellation_and_never_normalizes_bad_frames() {
        use std::cell::Cell;
        let (material, _, screw, frame, _) = super::super::tests::fixture(false, 1);
        let stock = Frame {
            origin_mm: [0.; 3],
            axes: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        };
        let calls = Cell::new(0);
        separated_from_stock(&material, stock, &screw, frame, &|| {
            calls.set(calls.get() + 1);
            false
        })
        .unwrap();
        for stop in [1, calls.get()] {
            let n = Cell::new(0);
            assert_eq!(
                separated_from_stock(&material, stock, &screw, frame, &|| {
                    n.set(n.get() + 1);
                    n.get() >= stop
                })
                .unwrap_err(),
                Error::Cancelled
            );
        }
        let prism = Prism::rectangular(7, [Rational::new(2, 1).unwrap(); 3]).unwrap();
        let n = Cell::new(0);
        assert_eq!(
            separated_from_rectangular_prism(&prism, stock, &screw, frame, &|| {
                n.set(n.get() + 1);
                n.get() >= 2
            })
            .unwrap_err(),
            Error::Cancelled
        );
        let invalid = Frame {
            axes: [[2., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            ..stock
        };
        assert!(separated_from_stock(&material, invalid, &screw, frame, &|| false).is_err());
    }
}
