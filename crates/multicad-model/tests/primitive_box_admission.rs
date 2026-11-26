//! The UI keeps dragged box sizes inside exact admission (CAD-401); the typed
//! API still refuses genuinely invalid or inadmissible boxes explicitly.
use multicad_core::{
    commands::{Actor, Channel, Session},
    document::Document,
    model::{self, PrimitiveFrame, PrimitiveParameters, PrimitiveRequest, length},
};
use multicad_model::{Compiler, GeometryBackend};

fn request(size: [f64; 3]) -> PrimitiveRequest {
    PrimitiveRequest::new(
        "Typed box",
        PrimitiveParameters::Box {
            size: size.map(length),
            frame: PrimitiveFrame::default(),
        },
    )
}

/// The native result of committing this box through the typed API.
fn build(size: [f64; 3]) -> Result<(), String> {
    let mut session = Session::new(Document::new("Typed box admission").unwrap()).unwrap();
    let request = request(size);
    let prepared = model::prepare_primitive(
        &session,
        session.revision(),
        Actor::local(Channel::Mcp),
        &request,
    )
    .map_err(|error| error.to_string())?;
    session.commit(prepared.prepared).unwrap();
    let model = Compiler::default()
        .compile_with_backend(
            session.document(),
            session.revision(),
            GeometryBackend::Rust,
            None,
        )
        .map_err(|error| error.to_string())?;
    model
        .solid(request.feature)
        .map(|_| ())
        .ok_or_else(|| "the box produced no solid".into())
}

#[test]
fn snapped_and_free_resolution_drag_sizes_build_natively() {
    build([90., 40., 20.]).unwrap();
    build([91.393, 41.678, 20.]).unwrap();
}

#[test]
fn typed_api_still_refuses_unadmissible_and_invalid_boxes_explicitly() {
    // The full binary64 values a free drag used to send (CAD-401): the exact
    // volume exceeds the bounded rational domain and stays refused.
    let error = build([
        91.393_323_123_456_78,
        41.678_413_987_654_32,
        20.123_456_789_012_34,
    ])
    .unwrap_err();
    assert!(error.contains("OutOfRange"), "{error}");
    for size in [[0., 10., 10.], [-5., 10., 10.], [10., f64::NAN, 10.]] {
        let error = build(size).unwrap_err();
        assert!(!error.is_empty(), "{size:?}");
    }
}
