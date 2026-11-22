//! Deliberately altered literal arrays, never forged source/valid-body receipts.
use super::*;
#[test]
fn actual_literal_cap_mutations_refuse_angular_transfer_and_normal_certification() {
    let source = super::super::super::tests::canonical(false, 1.);
    let atlas = source.stock_clip_plan(Default::default(), None).unwrap();
    let plan = atlas
        .cell_plan(128, 128, 512_000_000, 64 * 1024 * 1024, None)
        .unwrap();
    let quality = Tessellation {
        deflection: 0.1,
        angle: 0.3,
        triangle_limit: 1_000_000,
    };
    let candidate = plan
        .mesh_candidate(quality, 512_000_000, 64 * 1024 * 1024, None)
        .unwrap();
    let edge = &plan.caps[0];
    let ids = edge.ids;
    let mut mesh = candidate.mesh().clone();
    assert!(cap_edge(&plan, edge, &mesh, quality).is_ok());
    let mut nonfinite = mesh.clone();
    nonfinite.vertices[ids[0] as usize * 3] = f64::NAN;
    assert!(matches!(
        cap_edge(&plan, edge, &nonfinite, quality),
        Err(Error::Unresolved)
    ));
    // Reverse the actual directed edge by swapping its real coordinate leaves.
    for i in 0..3 {
        mesh.vertices
            .swap(ids[0] as usize * 3 + i, ids[1] as usize * 3 + i);
    }
    assert!(matches!(
        cap_edge(&plan, edge, &mesh, quality),
        Err(Error::Unresolved)
    ));
    // Keep its literal angular orientation but put it in the opposite halfplane.
    mesh = candidate.mesh().clone();
    for id in ids {
        for i in 0..2 {
            mesh.vertices[id as usize * 3 + i] *= -1.;
        }
    }
    assert!(matches!(
        cap_edge(&plan, edge, &mesh, quality),
        Err(Error::Unresolved)
    ));
    // Same halfplane/orientation, genuine geometric displacement outside C0.
    mesh = candidate.mesh().clone();
    for id in ids {
        for i in 0..2 {
            mesh.vertices[id as usize * 3 + i] *= 2.;
        }
    }
    assert!(matches!(
        cap_edge(&plan, edge, &mesh, quality),
        Err(Error::Unresolved)
    ));
    // A changed literal cap center breaks its required plane normal.
    mesh = candidate.mesh().clone();
    let start = edge.station == StockClipStation::Start;
    let center = plan.vertex_count() + usize::from(!start);
    let triangle = [center as u32, ids[0], ids[1]];
    let allowed = 2. * (quality.angle / 2.).portable_sin();
    assert!(cap_normal(&mesh, triangle, start, allowed).is_ok());
    mesh.vertices[center * 3 + 2] += 100.;
    assert!(matches!(
        cap_normal(&mesh, triangle, start, allowed),
        Err(Error::Unresolved)
    ));
    // Actual degenerate raw facet, independently rejected by the unchanged
    // indexed manufacturing certifier as well as this numerical normal stage.
    mesh = candidate.mesh().clone();
    let triangle: [u32; 3] = mesh.indices[..3].try_into().unwrap();
    for i in 0..3 {
        mesh.vertices[triangle[1] as usize * 3 + i] = mesh.vertices[triangle[0] as usize * 3 + i];
    }
    assert!(matches!(facet(&mesh, triangle), Err(Error::Unresolved)));
    assert!(manufacturing::certify(&mesh).is_err());
    mesh = candidate.mesh().clone();
    mesh.vertices[triangle[0] as usize * 3] = f64::NAN;
    assert!(matches!(facet(&mesh, triangle), Err(Error::Unresolved)));
    assert!(manufacturing::certify(&mesh).is_err());
    assert_eq!(*plan.source().source(), source);
    assert_eq!(
        plan.mesh_candidate(quality, 512_000_000, 64 * 1024 * 1024, None)
            .unwrap()
            .mesh(),
        candidate.mesh()
    );
}
#[test]
fn actual_raw_radius_and_parameter_winding_fail_closed_on_changed_geometry() {
    let source = super::super::super::tests::canonical(false, 1.);
    let atlas = source.stock_clip_plan(Default::default(), None).unwrap();
    let plan = atlas
        .cell_plan(128, 128, 512_000_000, 64 * 1024 * 1024, None)
        .unwrap();
    let candidate = plan
        .mesh_candidate(
            Tessellation {
                deflection: 0.1,
                angle: 0.3,
                triangle_limit: 1_000_000,
            },
            512_000_000,
            64 * 1024 * 1024,
            None,
        )
        .unwrap();
    assert!(candidate.maximum_literal_radius_mm() <= atlas.stock_radius_limit_mm());
    let original: [f64; 3] = candidate.mesh().vertices[..3].try_into().unwrap();
    assert!(stock_radius(original, atlas.stock_radius_limit_mm()).is_ok());
    assert!(matches!(
        stock_radius(
            [f64::NAN, original[1], original[2]],
            atlas.stock_radius_limit_mm()
        ),
        Err(Error::Unresolved)
    ));
    let outside = [original[0] * 2., original[1] * 2., original[2]];
    assert!(matches!(
        stock_radius(outside, atlas.stock_radius_limit_mm()),
        Err(Error::Unresolved)
    ));
    let cell = &plan.cells[0];
    let points = [cell.parameters[0], cell.parameters[1], cell.parameters[2]];
    assert!(positive_area(points).is_ok());
    assert!(matches!(
        positive_area([points[0], points[2], points[1]]),
        Err(Error::Unresolved)
    ));
    assert_eq!(*plan.source().source(), source);
}
