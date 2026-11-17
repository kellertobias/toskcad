//! Actual must-earn source windows; no coupled surface agreement claim.
use super::*;
use crate::step_read::brep::pcurve;

#[test]
fn installed_sma_three_general_partial_windows_use_actual_original_endpoints() {
    let path = "/Applications/KiCad/KiCad.app/Contents/SharedSupport/3dmodels/Connector_Coaxial.3dshapes/SMA_Samtec_SMA-J-P-H-ST-EM1_EdgeMount.step";
    let Ok(text) = std::fs::read_to_string(path) else {
        eprintln!("unavailable: external KiCad SMA fixture");
        return;
    };
    let source = Source::parse(&text, 1e-7).unwrap();
    let mut discovery = Budget::new(512_000_000, None);
    let bindings = pcurve::bindings(&source, &mut discovery).unwrap();
    for (edge_id, pc, target_id, uv_id) in [
        (2907, 2922, 2909, 2924),
        (2962, 2983, 2966, 2985),
        (8320, 8330, 8322, 8332),
    ] {
        // This is the actual public binding, not a constructed Edge/range.
        let binding = bindings
            .iter()
            .find(|binding| binding.edge == edge_id && binding.pcurve == pc)
            .copied()
            .unwrap();
        let mut budget = Budget::new(512_000_000, None);
        let actual_target = super::super::unwrapped(&source, binding.curve, &mut budget).unwrap();
        assert_eq!(actual_target, target_id);
        let edge = pcurve::edge(&source, binding, &mut budget).unwrap();
        let definition = pcurve::uv::definition(&source, pc, &mut budget).unwrap();
        assert_eq!(definition.curve_id, uv_id);
        let Curve::Nurbs { curve, .. } = &edge.curve else {
            panic!("actual SMA target is a spline")
        };
        let Curve::Nurbs { curve: uv, .. } = &definition.geometry else {
            panic!("actual SMA UV is a spline")
        };
        let packet = super::super::bezier::original_curve(
            &source,
            target_id,
            curve,
            Coordinates::Raw(3),
            &mut budget,
        )
        .unwrap();
        assert!(packet.matches(&source, target_id));
        let prefix = budget.used();
        match prepare(&source, target_id, &definition, &edge, &packet, &mut budget) {
            Ok(window) => {
                assert!(window.matches(&source, target_id, &definition, &edge));
                assert!(window.gamma() > 0. && window.speed() >= window.gamma());
                assert!(window.raw_uv_width().lo > window.closest_delta_sum());
                eprintln!(
                    "SMA_PARTIAL_WINDOW edge={edge_id} pc={pc} target={target_id} uv={uv_id} status=earned cached_edge={:?} cached_uv={:?} raw_uv={:?} gamma={:.17e} speed={:.17e} acceleration={:.17e} closest={:.17e} closest_sum={:.17e} cached_delta={:.17e} numeric={:.17e} prefix_work={prefix} window_work={} aggregate_work={} packet_retained={} packet_peak={}",
                    edge.range,
                    uv.range(),
                    window.raw_uv_domain().map(|x| [x.lo, x.hi]),
                    window.gamma(),
                    window.speed(),
                    window.acceleration(),
                    window.closest_delta(),
                    window.closest_delta_sum(),
                    window.cached_window_delta(),
                    window.numeric_reserve(),
                    budget.used() - prefix,
                    budget.used(),
                    packet.retained_bytes(),
                    packet.peak_bytes(),
                );
            }
            Err(error) => {
                // Preserve diagnostics, but every pinned case must earn its
                // source-window receipt. A refusal is a regression here.
                eprintln!(
                    "SMA_PARTIAL_WINDOW edge={edge_id} pc={pc} target={target_id} uv={uv_id} status=refused error={error:?} cached_edge={:?} cached_uv={:?} prefix_work={prefix} window_work={} aggregate_work={}",
                    edge.range,
                    uv.range(),
                    budget.used() - prefix,
                    budget.used(),
                );
                let mut diagnostic = Budget::new(512_000_000, None);
                eprintln!(
                    "SMA_PARTIAL_WINDOW_PREREQUISITES pc={pc} derivative_bounds={:?} target_domain={:?} uv_domain={:?}",
                    bounds(&packet, &mut diagnostic),
                    knot_domain(&source, target_id, curve, true, &mut diagnostic)
                        .map(|x| x.map(|x| [x.lo, x.hi])),
                    knot_domain(&source, uv_id, uv, false, &mut diagnostic)
                        .map(|x| x.map(|x| [x.lo, x.hi])),
                );
                panic!("actual SMA PC{pc} no longer earns its source window: {error:?}");
            }
        }
    }
}
