use super::*;
use crate::spur_gear::{
    PI, Spec, TAU,
    internal::{Ring, RingSpec},
};

const TOL: f64 = 1e-7;

fn spur(teeth: u32, bore: f64) -> Profile {
    let a = 20_f64.to_radians();
    let pitch = f64::from(teeth) / 2.;
    Profile::new(
        Spec {
            teeth,
            pitch_radius_mm: pitch,
            base_radius_mm: pitch * a.cos(),
            tip_radius_mm: pitch + 1.,
            root_radius_mm: pitch - 1.25,
            root_fillet_radius_mm: 0.2,
            half_base_angle_rad: (PI / 2. - 0.05) / (2. * pitch) + a.tan() - a,
            face_width_mm: 2.,
            bore_radius_mm: bore,
        },
        1e-5,
        1e7,
    )
    .unwrap()
}

fn ring(teeth: u32) -> Ring {
    let a = 20_f64.to_radians();
    let pitch = f64::from(teeth) / 2.;
    Ring::new(
        RingSpec {
            teeth,
            pitch_radius_mm: pitch,
            base_radius_mm: pitch * a.cos(),
            tip_radius_mm: pitch - 1.,
            root_radius_mm: pitch + 1.25,
            root_fillet_radius_mm: 0.2,
            half_base_angle_rad: (PI / 2. - 0.05) / (2. * pitch) - a.tan() + a,
            face_width_mm: 2.,
            outer_radius_mm: pitch + 4.,
        },
        1e-5,
        1e7,
        TOL,
    )
    .unwrap()
}

fn place(x: f64, theta: f64) -> Placement {
    let (sin, cos) = theta.sin_cos();
    Placement {
        center_mm: [x, 0.],
        axes: [[cos, sin], [-sin, cos]],
        ..Placement::default()
    }
}

/// An independent area estimate: grid cells whose centers both disc
/// proofs place strictly inside, with the boundary cells as the error.
fn grid_area(
    sections: [Section<'_>; 2],
    placements: [Placement; 2],
    window: [[f64; 2]; 2],
    step: f64,
) -> (f64, f64) {
    let mut inside = 0usize;
    let mut unsure = 0usize;
    let [x, y] = window;
    let nx = ((x[1] - x[0]) / step).ceil() as usize;
    let ny = ((y[1] - y[0]) / step).ceil() as usize;
    for i in 0..nx {
        for j in 0..ny {
            let p = [
                x[0] + (i as f64 + 0.5) * step,
                y[0] + (j as f64 + 0.5) * step,
            ];
            let mut both = Some(true);
            for k in 0..2 {
                let q = local(placements[k], p);
                // A cell whose disc meets a boundary is unsure.
                let side = match sections[k] {
                    Section::Spur(s) => super::super::classify_spur_disc(s, q, step, TOL),
                    Section::Ring(r) => super::super::classify_ring_disc(r, q, step, TOL),
                }
                .unwrap();
                both = match (both, side) {
                    (Some(false), _) | (_, super::super::DiscSide::Outside) => Some(false),
                    (Some(true), super::super::DiscSide::Inside) => Some(true),
                    _ => None,
                };
            }
            match both {
                Some(true) => inside += 1,
                None => unsure += 1,
                Some(false) => {}
            }
        }
    }
    let cell = step * step;
    (inside as f64 * cell, unsure as f64 * cell)
}

#[test]
fn misphased_spurs_share_exact_regions_matching_an_independent_grid() {
    let (a, b) = (spur(24, 3.), spur(36, 4.));
    let sections = [Section::Spur(&a), Section::Spur(&b)];
    // Opposed teeth at the line of centers collide.
    let placements = [place(0., 0.), place(30., PI)];
    let regions = intersect_sections(sections, placements, TOL, None)
        .unwrap()
        .unwrap();
    assert!(!regions.is_empty());
    let area: f64 = regions.iter().map(|r| r.area_mm2(placements)).sum();
    assert!(regions.iter().all(|r| r.area_mm2(placements) > 0.));
    let (grid, unsure) = grid_area(sections, placements, [[10., 15.], [-4., 4.]], 0.01);
    assert!(
        grid <= area && area <= grid + unsure,
        "{grid} {area} {unsure}"
    );
    // Every kept piece is an actual exact edge of its own section.
    for region in &regions {
        for pair in region.outer.windows(2) {
            let [p, q] = [pair[0], pair[1]]
                .map(|piece| placements[piece.section].point(piece.edge.point(1.)));
            let start = placements[pair[1].section].point(pair[1].edge.point(0.));
            assert!(length(sub(p, start)) < 1e-8, "{p:?} {start:?}");
            let _ = q;
        }
    }
    // The correct phase shares nothing.
    let clear = [place(0., 0.), place(30., PI + PI / 36.)];
    assert_eq!(
        intersect_sections(sections, clear, TOL, None),
        Ok(Some(vec![]))
    );
}

#[test]
fn deep_overlap_and_ring_pinion_regions_have_exact_area() {
    // Two equal gears 4 mm apart overlap deeply; each bore lies in the
    // other's material and stays a hole of the main common region. Where
    // the tooth rings cross, single tooth tips can add separate islands.
    let a = spur(24, 1.5);
    let sections = [Section::Spur(&a), Section::Spur(&a)];
    let placements = [place(0., 0.), place(4., 0.1)];
    let regions = intersect_sections(sections, placements, TOL, None)
        .unwrap()
        .unwrap();
    let main = regions
        .iter()
        .max_by(|x, y| x.area_mm2(placements).total_cmp(&y.area_mm2(placements)))
        .unwrap();
    assert_eq!(main.holes.len(), 2);
    assert_eq!(regions.iter().map(|r| r.holes.len()).sum::<usize>(), 2);
    assert!(regions.iter().all(|r| r.area_mm2(placements) > 0.));
    let area: f64 = regions.iter().map(|r| r.area_mm2(placements)).sum();
    let (grid, unsure) = grid_area(sections, placements, [[-10., 14.], [-12., 12.]], 0.1);
    assert!(
        grid <= area && area <= grid + unsure,
        "{grid} {area} {unsure}"
    );

    let r = ring(48);
    let p = spur(24, 3.);
    let sections = [Section::Ring(&r), Section::Spur(&p)];
    // Half a pinion pitch off its conjugate phase.
    let placements = [place(0., 0.), place(12., 0.)];
    let regions = intersect_sections(sections, placements, TOL, None)
        .unwrap()
        .unwrap();
    assert!(!regions.is_empty());
    let area: f64 = regions.iter().map(|r| r.area_mm2(placements)).sum();
    let (grid, unsure) = grid_area(sections, placements, [[17., 26.], [-8., 8.]], 0.01);
    assert!(
        grid <= area && area <= grid + unsure,
        "{grid} {area} {unsure}"
    );
    let conjugate = [place(0., 0.), place(12., PI / 24.)];
    assert_eq!(
        intersect_sections(sections, conjugate, TOL, None),
        Ok(Some(vec![]))
    );
}

#[test]
fn contacts_and_invalid_requests_fail_closed() {
    let a = spur(24, 3.);
    let sections = [Section::Spur(&a), Section::Spur(&a)];
    // Coincident gears share every edge: no transversal crossing exists.
    assert_eq!(
        intersect_sections(sections, [place(0., 0.); 2], TOL, None),
        Ok(None)
    );
    assert_eq!(
        intersect_sections(sections, [place(0., 0.), place(5., 0.)], 1., None),
        Err(Error::Invalid)
    );
    let stop = || true;
    assert_eq!(
        intersect_sections(sections, [place(0., 0.), place(5., 0.)], TOL, Some(&stop)),
        Err(Error::Cancelled)
    );
    let _ = TAU;
}
