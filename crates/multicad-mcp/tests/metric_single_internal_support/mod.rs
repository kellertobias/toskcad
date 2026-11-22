//! Actual normal MCP positive gates for the two single receiving treatments.
use multicad_core::{
    commands::{Actor, Channel, Session},
    document::Document,
    model::{self, PrimitiveParameters, PrimitiveRequest, length},
};
use multicad_mcp::{
    Headless, Workspace,
    tools::{Call, dispatch},
};
use multicad_model::GeometryBackend;
use serde_json::{Value, json};
use std::sync::{Arc, atomic::AtomicBool};
fn mm(v: f64) -> Value {
    json!({"value":v,"unit":"mm"})
}
pub fn require_positive(entry: bool) {
    let mut failures = Vec::new();
    for hand in ["right", "left"] {
        for height in [4., 3.3000000000000003] {
            let actor = Actor::local(Channel::Mcp);
            let mut session =
                Session::new(Document::new("MCP single receiving treatment").unwrap()).unwrap();
            let stock = PrimitiveRequest::new(
                "Actual cylinder stock",
                PrimitiveParameters::Cylinder {
                    radius: length(3.5),
                    height: length(height),
                    frame: Default::default(),
                },
            );
            let p = model::prepare_primitive(&session, 0, actor.clone(), &stock).unwrap();
            session.commit(p.prepared).unwrap();
            let host = Headless::with_geometry_backend(
                Workspace::new(session.document().clone()).unwrap(),
                None,
                false,
                GeometryBackend::Rust,
            );
            let ends = if entry {
                json!({"start":{"chamfer45":{"radial_depth":mm(0.4)}},"end":"flat_cut"})
            } else {
                json!({"start":"flat_cut","end":{"cylindrical_relief":{"length":mm(0.3),"radial_depth":mm(0.05)}}})
            };
            let args = json!({"expected_revision":0,"name":"M3 single receiving treatment","mode":"cut_internal","target_part":stock.body,"family":"iso_metric","fit":"metric_6h_internal","major_diameter":mm(3.),"pitch":mm(0.5),"length":mm(height),"minimum_wall":mm(0.5),"hand":hand,"end_treatments":ends});
            let before = host.document();
            match dispatch(
                &host,
                &Call {
                    actor,
                    cancel: Arc::new(AtomicBool::new(false)),
                },
                "set_thread",
                args,
            ) {
                Ok(created) => assert!(!created["thread"].is_null()),
                Err(e) => {
                    assert_eq!(host.document(), before);
                    failures.push(format!(
                        "entry={entry} hand={hand} height={height:?}: {e:?}"
                    ));
                }
            }
        }
    }
    assert!(
        failures.is_empty(),
        "Normal MCP single receiving modes: {failures:#?}"
    );
}
