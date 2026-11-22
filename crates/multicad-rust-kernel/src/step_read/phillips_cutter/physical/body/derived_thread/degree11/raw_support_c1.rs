//! Actual Core-authored, saved-file support correspondence, not file admission.
use super::*;
use crate::step_read::PhillipsRawCavitySupportOptions;
use crate::step_read::phillips_cutter::physical::body::common_mouth::tests::derived_thread_probe::with_head;
#[test]
fn actual_m3_flat_and_fitted_raw_cavity_supports_both_hands_have_whole_chart_c1() {
    for fitted in [false, true] {
        for left in [false, true] {
            let text = tests::actual_m3(fitted, left);
            with_head(&text, &mut |head| {
                let representative = head.common_mouth().closed_cavity().representative();
                let receipt = representative
                    .raw_support_correspondence(PhillipsRawCavitySupportOptions::default(), None)
                    .unwrap();
                let max_position = receipt
                    .faces()
                    .iter()
                    .map(|f| f.maximum_position_mm)
                    .fold(0f64, f64::max);
                let max_derivative = receipt
                    .faces()
                    .iter()
                    .flat_map(|f| f.maximum_chart_derivative_mm)
                    .fold(0f64, f64::max);
                eprintln!(
                    "actual Core M3 raw cavity fitted={fitted} left={left} C0={max_position:.17e} C1={max_derivative:.17e} fresh={} inherited={} total={} peak={}",
                    receipt.new_work_used(),
                    representative.work_used(),
                    receipt.work_used(),
                    receipt.accounted_peak_bytes()
                );
                assert!(max_position < 1e-7 && max_derivative < 1e-7);
                assert_eq!(
                    receipt
                        .faces()
                        .iter()
                        .map(|f| f.source_patches)
                        .sum::<usize>(),
                    88
                );
            });
        }
    }
}
