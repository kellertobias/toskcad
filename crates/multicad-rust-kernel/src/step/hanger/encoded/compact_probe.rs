//! Exact manufacturing-chart tiling oracle; not a file/material receipt.
use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
use std::collections::BTreeMap;
#[test]
fn multi_turn_meridian_strand_charts_exactly_reproduce_every_original_cell_pole() {
    for (left, nb, metric_starts, wood_starts) in [
        (false, false, 1, 1),
        (true, true, 1, 1),
        (false, true, 2, 3),
    ] {
        let mut spec = crate::authored_hanger_bolt::tests::fixture(left, nb);
        spec.metric = crate::authored_metric_thread::AuthoredMetricThread::new(
            crate::authored_metric_thread::Spec {
                starts: metric_starts,
                ..spec.metric.spec()
            },
        )
        .unwrap();
        spec.wood = crate::authored_flat_thread::AuthoredThread::new(
            crate::helical_thread::Spec {
                starts: wood_starts,
                ..spec.wood.spec()
            },
            crate::authored_flat_thread::SourceOptions {
                position_limit_mm: 1e-7,
                normal_angle_limit_rad: 1e-6,
            },
        )
        .unwrap();
        let source = Hanger::new(spec, None).unwrap();
        let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
        let field = EncodedField::prepare(&boundary, CurveOptions::default(), None).unwrap();
        let mut groups: BTreeMap<(usize, i32), Vec<usize>> = BTreeMap::new();
        for (index, cell) in field.cells.iter().enumerate() {
            let starts = if field.meridians[cell.meridian].component == Component::Metric {
                metric_starts
            } else {
                wood_starts
            } as i32;
            groups
                .entry((
                    cell.meridian,
                    cell.source_key.total_period.rem_euclid(starts),
                ))
                .or_default()
                .push(index);
        }
        let mut seen = vec![false; field.cells.len()];
        let mut poles = 0;
        let mut oracle_work = 0;
        for ((meridian, strand), cells) in &groups {
            let row = &field.meridians[*meridian];
            let (starts, pitch) = if row.component == Component::Metric {
                (metric_starts, spec.metric.spec().pitch_mm)
            } else {
                (wood_starts, spec.wood.spec().pitch_mm)
            };
            let starts = starts as i32;
            let reference = *cells
                .iter()
                .min_by_key(|&&i| field.cells[i].source_key.total_period)
                .unwrap();
            let reference_period = field.cells[reference].source_key.total_period;
            let start = reference_period.div_euclid(starts);
            let end = cells
                .iter()
                .map(|&i| field.cells[i].source_key.total_period.div_euclid(starts))
                .max()
                .unwrap()
                + 1;
            assert_eq!(cells.len(), (end - start) as usize);
            let mut prep = Work::new(16_000_000, None).unwrap();
            let lead = binary(pitch, &mut prep)
                .unwrap()
                .mul(&N::integer(starts.into()), &mut prep)
                .unwrap();
            let mut template = Vec::new();
            for u in 0..=row.degree() {
                for piece in 0..field.angular.coefficients().len() {
                    for v in 0..=10 {
                        template.push(
                            field
                                .homogeneous_control(reference, u, piece, v, &mut prep)
                                .unwrap()
                                .1,
                        );
                    }
                }
            }
            oracle_work += prep.used;
            for &index in cells {
                assert!(!seen[index]);
                seen[index] = true;
                let cell = field.cells[index];
                assert_eq!(cell.source_key.total_period.rem_euclid(starts), *strand);
                let turn = cell.source_key.total_period.div_euclid(starts);
                let mut w = Work::new(16_000_000, None).unwrap();
                let shift = lead
                    .mul(&N::integer((turn - start).into()), &mut w)
                    .unwrap();
                let mut control = 0;
                for u in 0..=row.degree() {
                    for piece in 0..field.angular.coefficients().len() {
                        for v in 0..=10 {
                            let h = field
                                .homogeneous_control(index, u, piece, v, &mut w)
                                .unwrap()
                                .1;
                            let base = &template[control];
                            for a in [0, 1, 3] {
                                assert_eq!(h[a].0, base[a].0)
                            }
                            assert_eq!(
                                h[2].0,
                                base[2]
                                    .add(&base[3].mul(&shift, &mut w).unwrap(), &mut w)
                                    .unwrap()
                                    .0
                            );
                            control += 1;
                            poles += 1;
                        }
                    }
                }
                oracle_work += w.used;
            }
        }
        assert!(seen.into_iter().all(|x| x));
        eprintln!(
            "exact compact charts left={left} nb={nb} starts={metric_starts}/{wood_starts} original_cells={} merged_faces={} checked_poles={poles} oracle_standalone_work={oracle_work}",
            field.cells.len(),
            groups.len()
        );
        assert!(groups.len() < field.cells.len());
    }
}
