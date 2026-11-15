use super::Error;
use super::*;
use crate::{
    authored_hanger_bolt::{BoundaryOptions, Hanger},
    step::hanger::*,
};
#[test]
fn actual_trimmed_candidate_is_independently_decoded_from_literal_step_entities() {
    for (left, nb) in [(false, false), (true, true)] {
        let source =
            Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
        let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
        let field = EncodedField::prepare(&boundary, CurveOptions::default(), None).unwrap();
        let domains = FiniteDomains::prepare(&field, DomainOptions::default(), None).unwrap();
        let closed = ClosedField::prepare(&domains, DomainOptions::default(), None).unwrap();
        let surfaces = CompactSurfaces::prepare(&field, DomainOptions::default(), None).unwrap();
        let compact =
            CompactBoundary::prepare(&closed, &surfaces, DomainOptions::default(), None).unwrap();
        let vertices = closed
            .literal_vertices(CurveOptions::default(), None)
            .unwrap();
        let events = vertices.cap_events(DomainOptions::default(), None).unwrap();
        let topology = compact
            .trim_caps(&events, DomainOptions::default(), None)
            .unwrap();
        let mut affine = Vec::new();
        let mut roots = Vec::new();
        for segment in 0..events.segments().len() {
            if let Some(cap) = events
                .bound_affine_cap(segment, CurveOptions::default(), None)
                .unwrap()
            {
                affine.push(cap)
            } else {
                roots.push(
                    events
                        .bound_root_cap(segment, CurveOptions::default(), None)
                        .unwrap()
                        .unwrap(),
                )
            }
        }
        let caps = affine
            .iter()
            .map(|c| c as &dyn CapLiteral)
            .chain(roots.iter().map(|c| c as &dyn CapLiteral))
            .collect::<Vec<_>>();
        let helices = (0..compact.edges().len())
            .filter_map(|e| {
                compact
                    .helix(e, &vertices, DomainOptions::default(), None)
                    .unwrap()
            })
            .collect::<Vec<_>>();
        let plains = (0..compact.faces().len())
            .filter_map(|f| {
                compact
                    .plain_wall(f, DomainOptions::default(), None)
                    .unwrap()
            })
            .collect::<Vec<_>>();
        let file = encode_candidate(
            &topology,
            &caps,
            &helices.iter().collect::<Vec<_>>(),
            &plains.iter().collect::<Vec<_>>(),
            FileOptions::default(),
            None,
        )
        .unwrap();
        let started = std::time::Instant::now();
        let decoded = decode(file.as_str(), ReadOptions::default(), None).unwrap();
        assert_eq!(decoded.literal_input(), file.as_str());
        assert_eq!(decoded.faces().len(), topology.faces().len());
        assert_eq!(decoded.edges().len(), topology.edges().len());
        assert_eq!(
            decoded.vertices().len(),
            topology
                .edges()
                .iter()
                .flat_map(|e| e.vertices)
                .collect::<BTreeSet<_>>()
                .len()
        );
        assert_eq!(decoded.planes().len(), 4);
        assert_eq!(decoded.tensors().len(), 14);
        assert!(!decoded.points().is_empty());
        assert!(!decoded.curves().is_empty());
        assert_ne!(decoded.solid_entity(), decoded.shell_entity());
        assert!(decoded.work_used() <= ReadOptions::default().maximum_work);
        assert!(decoded.retained_budget_bytes() <= ReadOptions::default().retained_budget_bytes);
        for point in decoded.points().values() {
            assert!(point.entity() > 0);
            for coordinate in point.coordinates() {
                let offset =
                    coordinate.spelling().as_ptr() as usize - file.as_str().as_ptr() as usize;
                assert_eq!(
                    file.as_str()
                        .get(offset..offset + coordinate.spelling().len()),
                    Some(coordinate.spelling())
                );
                assert!(coordinate.bounds()[0] <= coordinate.bounds()[1]);
            }
        }
        for curve in decoded.curves().values() {
            assert!(curve.entity() > 0);
            assert!((1..=32).contains(&curve.degree()));
            assert!(!curve.point_entities().is_empty());
            assert!(
                curve.weights().is_empty() || curve.weights().len() == curve.point_entities().len()
            );
            assert!(curve.knots().len() >= 2);
        }
        for tensor in decoded.tensors() {
            assert!(tensor.entity() > 0);
            assert_eq!(tensor.degrees()[1], 10);
            assert_eq!(
                tensor.point_entities().len(),
                tensor.dimensions().into_iter().product()
            );
            assert_eq!(tensor.weights().len(), tensor.point_entities().len());
            assert!(tensor.knots_v().len() >= 2);
        }
        for plane in decoded.planes() {
            assert!(plane.entity() > 0);
            assert!(decoded.points().contains_key(&plane.origin_entity()));
            assert!(plane.axis_entity() > 0);
            assert_ne!(plane.direction_entities()[0], plane.direction_entities()[1]);
        }
        for face in decoded.faces() {
            assert!(face.entity() > 0);
            assert!(face.surface_entity() > 0);
            let _ = face.surface();
            let _ = face.same_sense();
            assert!(!face.loops().is_empty());
        }
        for edge in decoded.edges().values() {
            assert!(edge.entity() > 0);
            assert!(decoded.curves().contains_key(&edge.curve_entity()));
            assert_ne!(edge.vertices()[0], edge.vertices()[1]);
            assert_ne!(edge.pcurves()[0].surface, edge.pcurves()[1].surface);
            let _ = edge.same_sense();
        }
        let grid = ExactPoleGrid::prepare(&decoded, GridOptions::default(), None).unwrap();
        assert!(std::ptr::eq(grid.source(), &decoded));
        assert!(grid.denominator().signum() == 1);
        use num_bigint::{BigInt, Sign};
        use num_rational::BigRational;
        let integer = |v: &PoleInteger| {
            BigInt::from_bytes_le(
                match v.signum() {
                    -1 => Sign::Minus,
                    0 => Sign::NoSign,
                    _ => Sign::Plus,
                },
                v.numerator_magnitude_le(),
            )
        };
        let denominator = integer(grid.denominator());
        assert_eq!(
            denominator,
            (BigInt::from(1) << grid.denominator_powers()[0])
                * BigInt::from(5).pow(grid.denominator_powers()[1])
        );
        let mut powers = BTreeMap::new();
        for (&id, point) in grid.points() {
            for (numerator, original) in point.iter().zip(decoded.points()[&id].coordinates()) {
                let text = original.spelling();
                let (mantissa, exp) = text.split_once('E').unwrap_or((text, "0"));
                let fraction = mantissa.split_once('.').map_or(0, |(_, s)| s.len());
                let exp = exp.parse::<i32>().unwrap() - fraction as i32;
                let digits = mantissa.replace('.', "");
                let n = BigInt::parse_bytes(digits.as_bytes(), 10).unwrap();
                let power = powers
                    .entry(exp.unsigned_abs())
                    .or_insert_with(|| BigInt::from(10).pow(exp.unsigned_abs()));
                let expected = if exp >= 0 {
                    BigRational::from_integer(n * &*power)
                } else {
                    BigRational::new(n, power.clone())
                };
                assert_eq!(
                    BigRational::new(integer(numerator), denominator.clone()),
                    expected
                );
            }
        }
        assert!(matches!(
            ExactPoleGrid::prepare(
                &decoded,
                GridOptions {
                    maximum_preparation_work: grid.preparation_work_used() - 1,
                    ..GridOptions::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        assert!(matches!(
            ExactPoleGrid::prepare(
                &decoded,
                GridOptions {
                    maximum_bytes: 1,
                    ..GridOptions::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        assert!(matches!(
            ExactPoleGrid::prepare(&decoded, GridOptions::default(), Some(&|| true)),
            Err(Error::Cancelled)
        ));
        eprintln!(
            "actual physical grid left={left} powers={:?} work={} aggregate={} bytes={} points={}",
            grid.denominator_powers(),
            grid.preparation_work_used(),
            grid.aggregate_known_work_used(),
            grid.retained_bytes(),
            grid.points().len()
        );
        let polar = PolarSupports::prepare(&grid, PolarOptions::default(), None).unwrap();
        assert!(matches!(
            PolarSupports::prepare(
                &grid,
                PolarOptions {
                    maximum_work: polar.preparation_work_used() - 1,
                    ..PolarOptions::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        assert!(matches!(
            PolarSupports::prepare(
                &grid,
                PolarOptions {
                    maximum_bytes: 1,
                    ..PolarOptions::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        assert!(matches!(
            PolarSupports::prepare(&grid, PolarOptions::default(), Some(&|| true)),
            Err(Error::Cancelled)
        ));
        assert!(std::ptr::eq(polar.grid(), &grid));
        assert_eq!(polar.strips().len(), 10);
        for strip in polar.strips() {
            assert!(strip.tensor_index() < decoded.tensors().len());
            assert!(strip.radius_controls_mm().iter().all(|v| v[0] > 0.));
            assert_eq!(strip.phase_controls_mm()[0], [0., 0.]);
            assert!(strip.phase_start_mm()[0] >= 0.);
            assert!(strip.phase_width_mm()[0] > 0.);
            assert!(matches!(strip.lower_plane_index(), 0 | 2));
            assert_eq!(strip.left_hand(), left);
            assert!(strip.weight_bounds().iter().all(|v| v[0] > 0.));
        }
        assert_eq!(polar.plane_entities().len(), 4);
        assert_eq!(
            polar.angular_piece_count(),
            polar.angular_coefficients().len()
        );
        eprintln!(
            "actual polar supports left={left} work={} aggregate={} bytes={} plain={:?} leads={:?} origins={:?} planes={:?}",
            polar.preparation_work_used(),
            polar.aggregate_known_work_used(),
            polar.retained_bytes(),
            polar.plain_radius_mm().unwrap(),
            polar.thread_leads_mm().unwrap(),
            polar.phase_origins_mm().unwrap(),
            polar.plane_heights_mm().unwrap()
        );
        let lifts = decoded.surface_lifts(LiftOptions::default(), None).unwrap();
        let occupied =
            HealedMaterial::prepare(&polar, &lifts, HealingOptions::default(), None).unwrap();
        assert!(matches!(
            HealedMaterial::prepare(
                &polar,
                &lifts,
                HealingOptions {
                    maximum_work: occupied.admission_work_used() - 1,
                    ..HealingOptions::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        assert!(matches!(
            HealedMaterial::prepare(
                &polar,
                &lifts,
                HealingOptions {
                    maximum_bytes: 1,
                    ..HealingOptions::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        assert!(matches!(
            HealedMaterial::prepare(&polar, &lifts, HealingOptions::default(), Some(&|| true)),
            Err(Error::Cancelled)
        ));
        assert!(matches!(
            HealedMaterial::prepare(
                &polar,
                &lifts,
                HealingOptions {
                    maximum_sewing_error_mm: occupied.sewing_error_bound_mm() * 0.5,
                    ..HealingOptions::default()
                },
                None
            ),
            Err(Error::Unsupported)
        ));
        assert!(std::ptr::eq(occupied.polar_supports(), &polar));
        assert!(std::ptr::eq(occupied.finite_lifts(), &lifts));
        assert!(
            occupied
                .component_radius_bounds_mm()
                .iter()
                .all(|r| r[0] > 0.)
        );
        assert!(occupied.angular_norm_bounds()[0] > 0.);
        assert!(occupied.healed_volume_bounds_mm3()[0] > 0.);
        let volume = occupied
            .healed_volume(VolumeOptions::default(), None)
            .unwrap();
        let reference = source
            .volume(
                crate::authored_hanger_bolt::QueryOptions {
                    volume_error_mm3: 1e-6,
                    maximum_cells: 65536,
                },
                None,
            )
            .unwrap()
            .material_mm3;
        assert!(
            volume.bounds_mm3()[0] <= reference[1] + 1e-6
                && volume.bounds_mm3()[1] >= reference[0] - 1e-6
        );
        eprintln!(
            "actual decoded healed volume left={left} bounds={:?} work={}",
            volume.bounds_mm3(),
            volume.work_used()
        );
        assert!(volume.work_used() > 0);
        assert!(matches!(
            occupied.healed_volume(
                VolumeOptions {
                    maximum_work: volume.work_used() - 1,
                    ..VolumeOptions::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        assert!(matches!(
            occupied.healed_volume(
                VolumeOptions {
                    maximum_bytes: 1,
                    ..VolumeOptions::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        assert!(matches!(
            occupied.healed_volume(VolumeOptions::default(), Some(&|| true)),
            Err(Error::Cancelled)
        ));
        assert!(matches!(
            occupied.healed_volume(
                VolumeOptions {
                    maximum_width_mm3: 1e-14,
                    ..VolumeOptions::default()
                },
                None
            ),
            Err(Error::Unsupported)
        ));
        let mut resolved = 0;
        let mut tested = 0;
        let mut most_work = (0, [0.; 3]);
        for z in [-12.713, -7.132, -0.713, 0.713, 5.132, 11.713] {
            for r in [0.7, 1.13, 1.37, 1.63, 1.83, 2.03] {
                for angle in [0., 0.173_f64, 1.371, 2.713, 4.113, 5.839] {
                    let point = [r * angle.cos(), r * angle.sin(), z];
                    let oracle = source.classify_box(point.map(|v| [v, v]), None).unwrap();
                    let query = occupied
                        .classify_point(point, QueryOptions::default(), None)
                        .unwrap();
                    tested += 1;
                    if query.work_used() > most_work.0 {
                        most_work = (query.work_used(), point);
                    }
                    if query.location() != PointLocation::Unresolved {
                        resolved += 1;
                        let expected = match oracle {
                            crate::authored_hanger_bolt::Classification::Inside => {
                                PointLocation::Inside
                            }
                            crate::authored_hanger_bolt::Classification::Outside => {
                                PointLocation::Outside
                            }
                            crate::authored_hanger_bolt::Classification::UnresolvedBoundary => {
                                continue;
                            }
                        };
                        assert_eq!(
                            query.location(),
                            expected,
                            "decoded point {point:?} left={left}"
                        );
                    }
                    assert!(query.work_used() > 0);
                }
            }
        }
        eprintln!("healed independent point query left={left} resolved={resolved}/{tested}");
        assert!(resolved * 10 > tested * 9);
        assert!(matches!(
            occupied.classify_point(
                most_work.1,
                QueryOptions {
                    maximum_work: most_work.0 - 1,
                    ..QueryOptions::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        assert!(matches!(
            occupied.classify_point(
                most_work.1,
                QueryOptions {
                    maximum_bytes: 1,
                    ..QueryOptions::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        assert!(matches!(
            occupied.classify_point(most_work.1, QueryOptions::default(), Some(&|| true)),
            Err(Error::Cancelled)
        ));
        assert!(matches!(
            occupied.classify_point([f64::NAN, 0., 0.], QueryOptions::default(), None),
            Err(Error::Invalid)
        ));
        assert_eq!(
            occupied
                .classify_point([0., 0., 1000.], QueryOptions::default(), None)
                .unwrap()
                .location(),
            PointLocation::Outside
        );
        assert_eq!(
            occupied
                .classify_point([0., 0., 0.], QueryOptions::default(), None)
                .unwrap()
                .location(),
            PointLocation::Unresolved
        );
        eprintln!("healed point maximum work left={left} {}", most_work.0);
        eprintln!(
            "healed actual field left={left} work={} aggregate={} bytes={} volume={:?} sewing={} jacobian={}",
            occupied.admission_work_used(),
            occupied.aggregate_known_work_used(),
            occupied.retained_budget_bytes(),
            occupied.healed_volume_bounds_mm3(),
            occupied.sewing_error_bound_mm(),
            occupied.minimum_angular_jacobian()
        );
        assert!(std::ptr::eq(lifts.source(), &decoded));
        assert_eq!(lifts.edges().len(), decoded.edges().len());
        for lift in lifts.edges() {
            assert!(decoded.edges().contains_key(&lift.edge_entity()));
            assert!(lift.field_bounds()[0] <= 1e-7 && lift.field_bounds()[1] <= 1e-6);
            assert!(lift.vertex_errors_mm().iter().all(|&e| e <= 1e-7));
            assert!(lift.parameter_bounds()[0][1] < lift.parameter_bounds()[1][0]);
            assert!(lift.work_used() <= LiftOptions::default().maximum_edge_work);
        }
        assert!(matches!(
            decoded.surface_lifts(
                LiftOptions {
                    maximum_work: lifts.work_used() - 1,
                    ..LiftOptions::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        assert!(matches!(
            decoded.surface_lifts(
                LiftOptions {
                    maximum_edge_work: 1,
                    ..LiftOptions::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        assert!(matches!(
            decoded.surface_lifts(
                LiftOptions {
                    maximum_bytes: 1,
                    ..LiftOptions::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        assert!(matches!(
            decoded.surface_lifts(LiftOptions::default(), Some(&|| true)),
            Err(Error::Cancelled)
        ));
        let calls = std::cell::Cell::new(0_usize);
        let count = || {
            calls.set(calls.get() + 1);
            false
        };
        decoded
            .surface_lifts(LiftOptions::default(), Some(&count))
            .unwrap();
        for checkpoint in 1..=calls.get().min(64) {
            let current = std::cell::Cell::new(0_usize);
            let cancel = || {
                current.set(current.get() + 1);
                current.get() == checkpoint
            };
            assert!(matches!(
                decoded.surface_lifts(LiftOptions::default(), Some(&cancel)),
                Err(Error::Cancelled)
            ));
        }
        if !left {
            let candidate = decoded
                .edges()
                .values()
                .find(|e| decoded.curves()[&e.curve_entity()].degree() > 10)
                .unwrap();
            let c = &decoded.curves()[&candidate.curve_entity()];
            let pole = *c
                .point_entities()
                .iter()
                .skip(1)
                .take(c.point_entities().len() - 2)
                .find(|&&id| {
                    decoded
                        .tensors()
                        .iter()
                        .all(|s| !s.point_entities().contains(&id))
                })
                .unwrap();
            let coordinates = decoded.points()[&pole].coordinates();
            let begin = file
                .as_str()
                .find(&format!("#{pole}=CARTESIAN_POINT("))
                .unwrap();
            let end = begin + file.as_str()[begin..].find(';').unwrap() + 1;
            let mut damaged = file.as_str().to_owned();
            damaged.replace_range(
                begin..end,
                &format!(
                    "#{pole}=CARTESIAN_POINT('',(999.,{},{}));",
                    coordinates[1].spelling(),
                    coordinates[2].spelling()
                ),
            );
            let corrupted = decode(&damaged, ReadOptions::default(), None).unwrap();
            assert_eq!(corrupted.faces().len(), decoded.faces().len());
            assert!(
                corrupted
                    .surface_lifts(LiftOptions::default(), None)
                    .is_err()
            );
            // Weld two nonadjacent topology vertices in the actual file.
            // Edge ownership and face-loop directions remain otherwise intact;
            // the resulting point has two local cycles and is not a manifold.
            let first = *decoded.vertices().keys().next().unwrap();
            let second = *decoded
                .vertices()
                .keys()
                .find(|&&v| {
                    v != first
                        && decoded
                            .edges()
                            .values()
                            .all(|e| !(e.vertices().contains(&first) && e.vertices().contains(&v)))
                })
                .unwrap();
            let mut pinched = file.as_str().to_owned();
            let lines = decoded
                .edges()
                .values()
                .filter(|e| e.vertices().contains(&second))
                .map(|e| e.entity())
                .collect::<Vec<_>>();
            for id in lines {
                let prefix = format!("#{id}=EDGE_CURVE(");
                let begin = pinched.find(&prefix).unwrap();
                let end = begin + pinched[begin..].find(';').unwrap() + 1;
                let old = pinched[begin..end].to_owned();
                let changed = old.replace(&format!("#{second},"), &format!("#{first},"));
                assert_ne!(old, changed);
                pinched.replace_range(begin..end, &changed);
            }
            assert!(matches!(
                decode(&pinched, ReadOptions::default(), None),
                Err(Error::Invalid)
            ));
            for dimension in [2, 3] {
                let changed = file.as_str().replacen(
                    &format!("GEOMETRIC_REPRESENTATION_CONTEXT({dimension})"),
                    &format!("GEOMETRIC_REPRESENTATION_CONTEXT({dimension}.00000000000000000001)"),
                    1,
                );
                assert!(matches!(
                    decode(&changed, ReadOptions::default(), None),
                    Err(Error::Invalid)
                ));
            }
            let centimetres =
                file.as_str()
                    .replacen("SI_UNIT(.MILLI.,.METRE.)", "SI_UNIT(.CENTI.,.METRE.)", 1);
            assert!(matches!(
                decode(&centimetres, ReadOptions::default(), None),
                Err(Error::Unsupported)
            ));
        }
        eprintln!(
            "actual finite lifts left={left} work={} bytes={} c0={} d1={}",
            lifts.work_used(),
            lifts.retained_bytes(),
            lifts
                .edges()
                .iter()
                .map(|e| e.field_bounds()[0])
                .fold(0., f64::max),
            lifts
                .edges()
                .iter()
                .map(|e| e.field_bounds()[1])
                .fold(0., f64::max)
        );
        eprintln!(
            "decoded left={left} work={} retained={} seconds={}",
            decoded.work_used(),
            decoded.retained_budget_bytes(),
            started.elapsed().as_secs_f64()
        );
        assert!(matches!(
            decode(
                file.as_str(),
                ReadOptions {
                    maximum_work: decoded.work_used() - 1,
                    ..ReadOptions::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        assert!(matches!(
            decode(file.as_str(), ReadOptions::default(), Some(&|| true)),
            Err(Error::Cancelled)
        ));
    }
}
