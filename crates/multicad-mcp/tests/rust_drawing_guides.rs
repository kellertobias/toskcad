//! External control generates Rust drawing views of a spur gear through
//! guides, edits the gear, regenerates the views and reports explicit
//! unresolved findings. Stale and cancelled requests change nothing.
use multicad_core::{
    commands::{Actor, Channel},
    document::{Document, EntityId},
    drawings::{self, VectorCurve, VectorRole, ViewEvaluation},
};
use multicad_mcp::{
    Headless, Workspace,
    tools::{Call, dispatch},
};
use multicad_precision::PortableMath as _;
use serde_json::{Value, json};
use std::sync::{Arc, atomic::AtomicBool};

fn call(cancelled: bool) -> Call {
    Call {
        actor: Actor {
            id: "mcp:drawing-guides".into(),
            channel: Channel::Mcp,
        },
        cancel: Arc::new(AtomicBool::new(cancelled)),
    }
}

fn mm(v: f64) -> Value {
    json!({ "value": v, "unit": "mm" })
}

/// Evaluate the host's current drawings with the default Rust backend. The
/// compile and the drawing evaluation share one local revision stamp.
fn views(host: &Headless) -> std::collections::BTreeMap<EntityId, ViewEvaluation> {
    let document = host.document();
    let revision = 0;
    let model = multicad_model::Compiler::default()
        .compile(&document, revision)
        .unwrap();
    model
        .drawings(
            &document,
            revision,
            &multicad_model::DrawingProviders::default(),
        )
        .unwrap()
        .views
}

fn bore_radii(evaluation: &ViewEvaluation) -> Vec<f64> {
    let ViewEvaluation::Ready { vectors, .. } = evaluation else {
        panic!("{evaluation:?}");
    };
    vectors
        .iter()
        .filter_map(|v| match v.curve {
            VectorCurve::Ellipse {
                major_mm,
                minor_mm,
                start_rad,
                end_rad,
                ..
            } if (major_mm - minor_mm).abs() < 1e-9
                && ((end_rad - start_rad).abs() - std::f64::consts::TAU).abs() < 1e-9 =>
            {
                Some(major_mm)
            }
            _ => None,
        })
        .collect()
}

/// Guide sections carry a lateral window, so their rims are the
/// tolerance-bounded sampled section. Count section points on a bore ring
/// around the gear axis, which projects to (60, 0) in the section frame.
fn sampled_bore(evaluation: &ViewEvaluation, radius: f64) -> usize {
    let ViewEvaluation::Ready { vectors, .. } = evaluation else {
        panic!("{evaluation:?}");
    };
    vectors
        .iter()
        .filter(|v| v.role == VectorRole::Section)
        .flat_map(|v| &v.preview_mm)
        .filter(|p| ((p[0] - 60.).portable_hypot(p[1]) - radius).abs() < 1e-6)
        .count()
}

fn id(value: &Value) -> EntityId {
    serde_json::from_value(value.clone()).unwrap()
}

#[test]
fn mcp_guides_draw_and_regenerate_rust_gear_sections() {
    let host = Headless::new(
        Workspace::new(Document::new("Gear guides").unwrap()).unwrap(),
        None,
        false,
    );
    let run = |name: &str, arguments: Value| dispatch(&host, &call(false), name, arguments);
    // A spur gear whose axis runs along +X, with an 8 mm tooth band.
    let gear = run(
        "set_gear",
        json!({ "expected_revision": 0, "name": "Spur", "family": "spur",
        "face_width": mm(8.), "bore_diameter": mm(6.),
        "axis_u": [0., 1., 0.], "axis_v": [0., 0., 1.] }),
    )
    .unwrap();
    let part = gear["part"].as_str().unwrap().to_owned();
    let feature = gear["gear"].clone();

    // A vertical section line along Y cuts the band across the gear axis.
    let section = |revision: u64, guide: Option<&str>| {
        let mut args = json!({ "expected_revision": revision,
            "view": { "kind": "section", "start_mm": [3., -60.], "end_mm": [3., 60.] },
            "sources": [part], "scale_denominator": 1, "generate_pages": true });
        match guide {
            Some(guide) => args["guide"] = json!(guide),
            None => args["name"] = json!("Across the axis"),
        }
        run("set_drawing_guide", args)
    };
    let created = section(1, None).unwrap();
    let generation = &created["generation"];
    let [view] = generation["created"].as_array().unwrap().as_slice() else {
        panic!("{created}");
    };
    let view = id(view);
    assert!(
        generation["unresolved"].as_array().unwrap().is_empty(),
        "{created}"
    );
    let guide = created["guide"].as_str().unwrap().to_owned();
    let evaluated = views(&host);
    assert!(sampled_bore(&evaluated[&view], 3.) > 16);

    // An elevation along the axis is the exact axial gear view.
    let elevation = run(
        "set_drawing_guide",
        json!({ "expected_revision": 2, "name": "Along the axis",
        "view": { "kind": "elevation", "direction_deg": 0 },
        "sources": [part], "scale_denominator": 1, "generate_pages": true }),
    )
    .unwrap();
    let axial = id(&elevation["generation"]["created"][0]);
    assert!(
        elevation["generation"]["unresolved"]
            .as_array()
            .unwrap()
            .is_empty()
    );
    assert_eq!(bore_radii(&views(&host)[&axial]), [3.]);

    // A section containing the gear axis draws through the exact-curve scene
    // engine (CAD-346): the cut walls of the bore are two Section lines 6 mm
    // apart along the 8 mm tooth band.
    let along = run(
        "set_drawing_guide",
        json!({ "expected_revision": 3, "name": "Through the axis",
        "view": { "kind": "section", "start_mm": [-40., 0.], "end_mm": [40., 0.] },
        "sources": [part], "scale_denominator": 1, "generate_pages": true }),
    )
    .unwrap();
    assert!(
        along["generation"]["unresolved"]
            .as_array()
            .unwrap()
            .is_empty(),
        "{along}"
    );
    let through = id(&along["generation"]["created"][0]);
    let ViewEvaluation::Ready { vectors, .. } = &views(&host)[&through] else {
        panic!("{along}");
    };
    let walls: Vec<[[f64; 2]; 2]> = vectors
        .iter()
        .filter(|v| v.role == VectorRole::Section)
        .filter_map(|v| match v.curve {
            VectorCurve::Line { a, b }
                if ((b[0] - a[0]).portable_hypot(b[1] - a[1]) - 8.).abs() < 1e-9 =>
            {
                Some([a, b])
            }
            _ => None,
        })
        .collect();
    let bore = walls.iter().enumerate().any(|(i, p)| {
        walls[i + 1..].iter().any(|q| {
            let d = [p[1][0] - p[0][0], p[1][1] - p[0][1]];
            let offset = [q[0][0] - p[0][0], q[0][1] - p[0][1]];
            ((offset[0] * d[1] - offset[1] * d[0]).abs() / 8. - 6.).abs() < 1e-9
        })
    });
    assert!(bore, "{walls:?}");

    // Editing the gear changes the regenerated section.
    run(
        "set_gear",
        json!({ "expected_revision": 4, "gear": feature, "bore_diameter": mm(8.) }),
    )
    .unwrap();
    let updated = section(5, Some(&guide)).unwrap();
    assert_eq!(
        updated["generation"]["updated"]
            .as_array()
            .unwrap()
            .as_slice(),
        [json!(view)]
    );
    assert!(
        updated["generation"]["created"]
            .as_array()
            .unwrap()
            .is_empty()
    );
    let evaluated = views(&host);
    assert_eq!(sampled_bore(&evaluated[&view], 3.), 0);
    assert!(sampled_bore(&evaluated[&view], 4.) > 16);
    assert_eq!(bore_radii(&evaluated[&axial]), [4.]);

    // Stale and cancelled requests leave the document untouched.
    let current = updated["revision"].as_u64().unwrap();
    let before = host.document();
    let stale = section(current - 1, Some(&guide)).unwrap_err();
    assert!(format!("{stale:?}").contains("stale_revision"), "{stale:?}");
    let cancelled = dispatch(
        &host,
        &call(true),
        "set_drawing_guide",
        json!({ "expected_revision": current, "guide": guide,
        "view": { "kind": "section", "start_mm": [3., -60.], "end_mm": [3., 60.] },
        "sources": [part], "scale_denominator": 1, "generate_pages": true }),
    );
    assert!(cancelled.is_err());
    assert_eq!(host.document(), before);
    let book = updated["generation"]["book"].clone();
    assert!(drawings::read(&before, id(&book)).is_ok());
}
