use multicad_core::{
    constraint_sketch::Dimension,
    document::EntityId,
    model::{
        self, MechanicalPartParameters, ScrewDrive, ScrewHead, ScrewLengthInput,
        ThreadRepresentation,
    },
    project_file,
    references::{AnalyticGeometry, TopologyKind},
};
use multicad_kernel::{GeometryErrorKind, MeshOptions, PlaneFrame};
use multicad_model::{Compiler, GeometryBackend};
use std::f64::consts::PI;

fn mm(dimension: &Dimension) -> f64 {
    let Dimension::Literal { value, .. } = dimension else {
        panic!("Expected a literal fastener dimension, found {dimension:?}");
    };
    *value
}

/// Independent cosmetic-envelope expectation: volume and world bounds.
///
/// A nut is its hexagonal prism less a bore at the ISO basic minor diameter
/// D1 = D - 5√3/8 P, extruded from the frame origin along u × v. A
/// countersunk screw is a flat top land, a cone to the shank, a full-length
/// Ø major shank and a flat hexagon socket, less a 45° ring at the head's
/// top edge; its head reaches back from the frame origin.
fn expected(parameters: &MechanicalPartParameters) -> (f64, [[f64; 3]; 2]) {
    let (volume, thread, half_u, half_v, axial) = match parameters {
        MechanicalPartParameters::HexNut {
            across_flats,
            thread,
        } => {
            let flats = mm(across_flats);
            let thickness = mm(&thread.length);
            let minor = mm(&thread.major_diameter) - 5. * 3_f64.sqrt() / 8. * mm(&thread.pitch);
            let volume = thickness * (3_f64.sqrt() / 2. * flats * flats - PI * minor * minor / 4.);
            (
                volume,
                thread,
                flats / 2.,
                flats / 3_f64.sqrt(),
                [0., thickness],
            )
        }
        MechanicalPartParameters::MachineScrew {
            underhead_length: None,
            length_input: Some(ScrewLengthInput::Overall { value }),
            shank_diameter,
            head: ScrewHead::Countersunk { profile },
            drive:
                ScrewDrive::HexSocket {
                    across_flats,
                    depth,
                },
            thread,
        } => {
            let overall = mm(value);
            let r = mm(shank_diameter) / 2.;
            assert_eq!(r, mm(&thread.major_diameter) / 2.);
            let radius = mm(&profile.diameter) / 2.;
            let height = mm(&profile.height);
            let chamfer = mm(&profile.edge_chamfer);
            assert_eq!(mm(&profile.neck_radius), 0.);
            let cone = (radius - r) / (mm(&profile.included_angle).to_radians() / 2.).tan();
            let land = height - cone;
            assert!(land > chamfer);
            let flats = mm(across_flats);
            let head = PI * radius * radius * land
                + PI * cone / 3. * (radius * radius + radius * r + r * r)
                - PI * chamfer * chamfer * (radius - chamfer / 3.)
                - 3_f64.sqrt() / 2. * flats * flats * mm(depth);
            let volume = head + PI * r * r * (overall - height);
            (volume, thread, radius, radius, [-height, overall - height])
        }
        other => panic!("Unexpected fastener {other:?}"),
    };
    let frame = &thread.frame;
    let origin = frame.origin.each_ref().map(mm);
    let [u, v] = [frame.u, frame.v];
    let n = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    let mut bounds = [[f64::INFINITY; 3], [f64::NEG_INFINITY; 3]];
    for a in [-half_u, half_u] {
        for b in [-half_v, half_v] {
            for c in axial {
                for i in 0..3 {
                    let x = origin[i] + a * u[i] + b * v[i] + c * n[i];
                    bounds[0][i] = bounds[0][i].min(x);
                    bounds[1][i] = bounds[1][i].max(x);
                }
            }
        }
    }
    (volume, bounds)
}

#[test]
fn the_toolcase_cosmetic_screws_and_nuts_match_analytic_envelopes() {
    let mut bytes = include_bytes!("../../../examples/toolcase.tcd").as_slice();
    let opened = project_file::read(
        &mut bytes,
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap();
    let document = opened.project.document;
    let mut roots = vec![];
    let mut counts = [0, 0];
    for entity in document.entities() {
        let Ok(source) = model::read_mechanical_part(&document, entity.id) else {
            continue;
        };
        let (thread, slot) = match source.parameters {
            MechanicalPartParameters::MachineScrew { thread, .. } => (thread, 0),
            MechanicalPartParameters::HexNut { thread, .. } => (thread, 1),
            _ => continue,
        };
        assert_eq!(
            thread.representation,
            ThreadRepresentation::CosmeticEnvelope
        );
        roots.push((entity.id, slot));
        counts[slot] += 1;
    }
    assert_eq!(counts, [43, 25]);
    let sources: Vec<EntityId> = roots.iter().map(|&(feature, _)| feature).collect();
    let compile = |backend| {
        Compiler::default()
            .compile_sources_with_backend(&document, 0, &sources, backend, None)
            .unwrap()
    };
    let rust = compile(GeometryBackend::Rust);
    for (feature, slot) in roots {
        let actual = &rust.snapshot().shapes[&feature];
        let (volume, [min, max]) = expected(
            &model::read_mechanical_part(&document, feature)
                .unwrap()
                .parameters,
        );
        let label = &document.get(feature).unwrap().name;
        assert_eq!(actual.solid_count, 1, "{label}");
        assert!(
            (actual.volume_mm3 - volume).abs() < 1e-7,
            "{label}: {} vs {volume}",
            actual.volume_mm3,
        );
        for axis in 0..3 {
            assert!(
                (actual.bounds.min[axis] - min[axis]).abs() < 1e-6,
                "{label}"
            );
            assert!(
                (actual.bounds.max[axis] - max[axis]).abs() < 1e-6,
                "{label}"
            );
        }
        assert!(!actual.mesh.indices.is_empty());
        let role = if slot == 0 {
            "screw.head/top.land"
        } else {
            "cap.start"
        };
        let top = actual
            .topology
            .elements
            .iter()
            .find(|element| {
                element.key.kind == TopologyKind::Face
                    && element
                        .origins
                        .iter()
                        .any(|origin| origin.source == feature && origin.role == role)
            })
            .unwrap();
        let Some(AnalyticGeometry::Plane { frame }) = top.geometry else {
            panic!("Expected planar ring on {label}");
        };
        let solid = rust.solid(feature).unwrap();
        assert!(solid.uses_rust_geometry());
        assert!(!solid.supports_step_export());
        let error = solid
            .manufacturing_mesh(MeshOptions::default())
            .unwrap_err();
        assert_eq!(error.kind, GeometryErrorKind::Unsupported);
        assert!(error.message.contains("modeled thread geometry"));
        let profile = solid
            .planar_face_contours(
                top.key.ordinal,
                PlaneFrame {
                    origin: frame.origin_mm,
                    u: frame.u,
                    v: frame.v,
                },
            )
            .unwrap();
        assert_eq!(profile.loops.len(), 2);
        assert!(profile.loops[0].outer);
        assert!(!profile.loops[1].outer);
    }
}
