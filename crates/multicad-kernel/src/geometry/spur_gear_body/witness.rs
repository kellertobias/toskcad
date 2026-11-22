//! Actual tooth-flank proximity witnesses for placed Rust gear pairs. Each
//! witness is a pair of points evaluated on the exact flank curves of the
//! two placed solids, so its distance bounds the true flank gap from above.
use super::*;
use crate::rust_kernel::spur_gear::{
    Edge as ToothCurve, ToothEdge, ToothRole,
    clearance::{self, Placement},
};

const OPERATION: &str = "gear flank witness";

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GearFlankWitness {
    pub distance_mm: f64,
    /// World points on the first and second solid.
    pub points_mm: [[f64; 3]; 2],
    /// Evaluation-local wall face ordinals, when both bodies number them.
    pub faces: Option<[u32; 2]>,
}

fn flanks(outline: &[ToothEdge]) -> Vec<(usize, ToothCurve)> {
    outline
        .iter()
        .enumerate()
        .filter(|(_, item)| matches!(item.role, ToothRole::LeftFlank | ToothRole::RightFlank))
        .map(|(index, item)| (index, item.edge))
        .collect()
}

fn failure(cause: clearance::Error) -> GeometryError {
    let kind = match cause {
        clearance::Error::Invalid => GeometryErrorKind::Unsupported,
        clearance::Error::Limit => GeometryErrorKind::Limit,
        clearance::Error::Cancelled => GeometryErrorKind::Cancelled,
    };
    error(OPERATION, kind, format!("Rust flank witness: {cause:?}"))
}

fn distance(points: [[f64; 3]; 2]) -> f64 {
    (0..3)
        .map(|i| (points[0][i] - points[1][i]).powi(2))
        .sum::<f64>()
        .sqrt()
}

/// Local height of `frame` on the plane at `height` along `base`'s normal.
fn local_height(base: PlaneFrame, frame: PlaneFrame, height: f64) -> f64 {
    let normal = base.normal();
    let offset = (0..3)
        .map(|i| (frame.origin[i] - base.origin[i]) * normal[i])
        .sum::<f64>();
    let vertical = (0..3).map(|i| frame.normal()[i] * normal[i]).sum::<f64>();
    (height - offset) / vertical
}

impl Solid {
    /// Search the tooth flanks of an aligned Rust spur, rack, helical or
    /// chevron pair for actual boundary points within `reach_mm`, in the
    /// transverse section at the middle of each common tooth band: one band
    /// for straight and single-hand gears, lower then upper for chevrons.
    /// An entry is `None` when that section has no flank points in reach.
    /// This is an existence witness only; clearance is a separate proof.
    pub fn rust_gear_flank_witnesses(
        &self,
        other: &Self,
        policy: GeometryPolicy,
        reach_mm: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Vec<Option<GearFlankWitness>>, GeometryError> {
        policy.validate()?;
        let unsupported = |message| error(OPERATION, GeometryErrorKind::Unsupported, message);
        if self.policy != policy || other.policy != policy {
            return Err(unsupported("Rust flank witnesses require a shared policy"));
        }
        if !(reach_mm.is_finite() && reach_mm > 0.) {
            return Err(input("Flank witness reach must be positive and finite"));
        }
        if let (Some(a), Some(b)) = (self.spur_gear_body(), other.spur_gear_body()) {
            return spur_witness(a, b, policy, reach_mm, interrupt).map(|w| vec![w]);
        }
        if let Some(pair) = self.rack_pair(other) {
            return rack_witness(&pair, policy, reach_mm, interrupt).map(|w| vec![w]);
        }
        if let (Some(a), Some(b)) = (self.herringbone_body(), other.herringbone_body()) {
            return twisted_witnesses(a, b, policy, reach_mm, interrupt);
        }
        Err(unsupported(
            "Rust flank witnesses require two spur, two twisted, or a spur and rack solid",
        ))
    }
}

fn spur_witness(
    (a, af): (&SpurGearBody, PlaneFrame),
    (b, bf): (&SpurGearBody, PlaneFrame),
    policy: GeometryPolicy,
    reach_mm: f64,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Option<GearFlankWitness>, GeometryError> {
    let profiles = [&a.profile, &b.profile];
    let (placements, shared) = coaxial_chart(
        (OPERATION, "Rust spur flank witness"),
        policy,
        [af, bf],
        profiles.map(|profile| Band {
            height_mm: profile.height_mm(),
            extent_mm: profile.axial_bounds(),
            support: Box::new(|direction| profile.support(direction)),
        }),
    )?;
    let edges = profiles.map(|profile| flanks(profile.outline()));
    let Some(found) =
        clearance::nearest_boundary([&edges[0], &edges[1]], placements, reach_mm, interrupt)
            .map_err(failure)?
    else {
        return Ok(None);
    };
    let middle = (shared[0] + shared[1]) / 2.;
    let points = [(af, found.local_mm[0]), (bf, found.local_mm[1])]
        .map(|(frame, [x, y])| position(frame, [x, y, local_height(af, frame, middle)]));
    Ok(within(
        points,
        reach_mm,
        Some(found.edges.map(|e| e as u32 + 2)),
    ))
}

fn rack_witness(
    pair: &RackPair<'_>,
    policy: GeometryPolicy,
    reach_mm: f64,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Option<GearFlankWitness>, GeometryError> {
    let (placements, shared) = pair.chart(policy)?;
    let pinion = flanks(pair.pinion.0.profile.outline());
    let rack = pair
        .rack
        .0
        .iter()
        .enumerate()
        .filter(|(index, _)| (pair.roles)(*index).contains("/flank/"))
        .map(|(index, &edge)| (index, ToothCurve::Planar(edge)))
        .collect::<Vec<_>>();
    if rack.is_empty() {
        return Err(error(
            OPERATION,
            GeometryErrorKind::Unsupported,
            "Rust rack flank witness requires named rack flanks",
        ));
    }
    let Some(found) =
        clearance::nearest_boundary([&pinion, &rack], placements, reach_mm, interrupt)
            .map_err(failure)?
    else {
        return Ok(None);
    };
    let base = pair.pinion.1;
    let middle = (shared[0] + shared[1]) / 2.;
    let mut points = [(base, found.local_mm[0]), (pair.rack.2, found.local_mm[1])]
        .map(|(frame, [x, y])| position(frame, [x, y, local_height(base, frame, middle)]));
    if pair.swapped {
        points.reverse();
    }
    Ok(within(points, reach_mm, None))
}

fn twisted_witnesses(
    (a, af): (&super::super::herringbone_body::HerringboneBody, PlaneFrame),
    (b, bf): (&super::super::herringbone_body::HerringboneBody, PlaneFrame),
    policy: GeometryPolicy,
    reach_mm: f64,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Vec<Option<GearFlankWitness>>, GeometryError> {
    let gears = [&a.geometry, &b.geometry];
    let height = gears[0].profile().height_mm();
    if gears[1].profile().height_mm() != height || gears[0].gap_mm() != gears[1].gap_mm() {
        return Err(error(
            OPERATION,
            GeometryErrorKind::Unsupported,
            "Rust twisted flank witnesses require equal face widths and gaps",
        ));
    }
    let base =
        super::super::herringbone_body::aligned_chart(OPERATION, policy, [(a, af), (b, bf)])?;
    let band = gears[0].band_length_mm();
    // A single-hand sweep has one band; its lower-half wall holds the
    // quarter-height section. Chevron band middles stay inside each band.
    let sections: &[(f64, u32)] = if gears.iter().any(|gear| gear.is_helical()) {
        &[(height * 0.25, 0)]
    } else {
        &[(band * 0.5, 0), (height - band * 0.5, 1)]
    };
    let edges = gears.map(|gear| flanks(gear.profile().outline()));
    let counts = gears.map(|gear| gear.profile().outline().len() as u32);
    let mut result = Vec::with_capacity(sections.len());
    for &(z, wall) in sections {
        let placed = std::array::from_fn(|i| {
            let [u, v] = base[i].axes;
            let (sin, cos) = gears[i].phase_at(z).portable_sin_cos();
            Placement {
                axes: [
                    [u[0] * cos + v[0] * sin, u[1] * cos + v[1] * sin],
                    [v[0] * cos - u[0] * sin, v[1] * cos - u[1] * sin],
                ],
                ..base[i]
            }
        });
        let found =
            clearance::nearest_boundary([&edges[0], &edges[1]], placed, reach_mm, interrupt)
                .map_err(failure)?;
        result.push(found.and_then(|found| {
            let points = [
                (af, gears[0], found.local_mm[0]),
                (bf, gears[1], found.local_mm[1]),
            ]
            .map(|(frame, gear, local)| position(frame, gear.section_point(local, z)));
            let faces = std::array::from_fn(|i| 2 + wall * counts[i] + found.edges[i] as u32);
            within(points, reach_mm, Some(faces))
        }));
    }
    Ok(result)
}

/// Report only actual 3D points whose distance is within reach.
fn within(
    points: [[f64; 3]; 2],
    reach_mm: f64,
    faces: Option<[u32; 2]>,
) -> Option<GearFlankWitness> {
    let distance_mm = distance(points);
    (distance_mm <= reach_mm).then_some(GearFlankWitness {
        distance_mm,
        points_mm: points,
        faces,
    })
}

/// The exact section boundary of a spur: tooth outline, then any bore.
fn section(profile: &multicad_rust_kernel::spur_gear::Profile) -> Vec<(usize, ToothCurve)> {
    profile
        .outline()
        .iter()
        .map(|item| item.edge)
        .chain(profile.bore().map(ToothCurve::Planar))
        .enumerate()
        .collect()
}

impl Solid {
    /// Certified minimum distance between a clear Rust spur pair or a spur
    /// and a rack sweep sharing an axial direction. Both solids project
    /// into their exact sections, so the section gap bounds the spatial gap
    /// from below; the witness is a pair of actual points in the common band
    /// and bounds it from above. `None` leaves other pairs to other measures.
    pub(crate) fn rust_gear_section_distance(
        &self,
        other: &Self,
        tolerance_mm: f64,
    ) -> Option<Result<crate::MinimumDistance, GeometryError>> {
        const DISTANCE: &str = "gear minimum distance";
        let policy = self.policy;
        let measured = || -> Result<crate::MinimumDistance, GeometryError> {
            let unsupported =
                |message: &str| error(DISTANCE, GeometryErrorKind::Unsupported, message);
            let (placements, shared, edges, frames, swapped, clear) =
                if let (Some(a), Some(b)) = (self.spur_gear_body(), other.spur_gear_body()) {
                    let profiles = [&a.0.profile, &b.0.profile];
                    let (placements, shared) = coaxial_chart(
                        (DISTANCE, "Rust gear distance"),
                        policy,
                        [a.1, b.1],
                        profiles.map(|profile| Band {
                            height_mm: profile.height_mm(),
                            extent_mm: profile.axial_bounds(),
                            support: Box::new(|direction| profile.support(direction)),
                        }),
                    )?;
                    let clear = self.certify_rust_spur_pair(other, policy, None)?;
                    (
                        placements,
                        shared,
                        profiles.map(section),
                        [a.1, b.1],
                        false,
                        clear,
                    )
                } else {
                    let pair = self
                        .rack_pair(other)
                        .ok_or_else(|| unsupported("Rust gear distance requires a spur pair"))?;
                    let (placements, shared) = pair.chart(policy)?;
                    let rack = pair
                        .rack
                        .0
                        .iter()
                        .map(|&edge| ToothCurve::Planar(edge))
                        .enumerate()
                        .collect();
                    let clear = self.certify_rust_rack_pair(other, policy, None)?;
                    (
                        placements,
                        shared,
                        [section(&pair.pinion.0.profile), rack],
                        [pair.pinion.1, pair.rack.2],
                        pair.swapped,
                        clear,
                    )
                };
            if clear != clearance::ResultKind::Clear {
                return Err(unsupported(
                    "Rust gear distance requires a certified clear pair",
                ));
            }
            let found = clearance::certify_boundary_distance(
                [&edges[0], &edges[1]],
                placements,
                tolerance_mm.clamp(1e-7, 1e-3),
                None,
            )
            .map_err(failure)?
            .ok_or_else(|| unsupported("Rust gear distance bounds did not meet"))?;
            let middle = (shared[0] + shared[1]) / 2.;
            let base = frames[0];
            let mut points = [
                (frames[0], found.local_mm[0]),
                (frames[1], found.local_mm[1]),
            ]
            .map(|(frame, [x, y])| position(frame, [x, y, local_height(base, frame, middle)]));
            if swapped {
                points.reverse();
            }
            Ok(crate::MinimumDistance {
                distance_mm: distance(points),
                witnesses: vec![crate::DistanceWitness {
                    a_mm: points[0],
                    b_mm: points[1],
                }],
                witnesses_complete: false,
                total_solutions: None,
                inner_solution: false,
            })
        };
        let gear_pair = (self.spur_gear_body().is_some() && other.spur_gear_body().is_some())
            || self.rack_pair(other).is_some();
        gear_pair.then(measured)
    }
}
