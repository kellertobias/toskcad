use super::tests::{argument_change, file, record_change};
use super::*;
use multicad_precision::PortableMath as _;
fn refused(text: &str, label: &str) {
    let p = parse(text, Options::default(), None)
        .unwrap_or_else(|e| panic!("{label} structural {e:?}"));
    assert!(
        p.support_bindings(Options::default(), None).is_err(),
        "{label}"
    );
}
#[test]
fn actual_dual_support_owner_heights_knots_and_uv_changes_are_rejected() {
    let text = file(false);
    let p = parse(&text, Options::default(), None).unwrap();
    let origin = p.point(p.planes()[1].origin_entity()).unwrap();
    let c = origin.coordinates();
    refused(
        &argument_change(
            &text,
            origin.entity(),
            1,
            &format!("({},{},0.)", c[0].spelling(), c[1].spelling()),
        ),
        "coincident cap heights",
    );
    let f0 = &p.faces()[0];
    let f1 = &p.faces()[1];
    refused(
        &argument_change(&text, f1.entity(), 2, &format!("#{}", f0.surface_entity())),
        "duplicate support owner",
    );
    let edge = &p.edges()[0];
    let pc = edge.pcurves()[0];
    refused(
        &argument_change(
            &text,
            edge.pcurves()[1].entity,
            1,
            &format!("#{}", pc.support_entity),
        ),
        "ambiguous dual support",
    );
    let cap = p
        .faces()
        .iter()
        .find(|f| matches!(f.surface(), SurfaceIndex::Plane(_)))
        .unwrap();
    let e = p.edge(cap.loops()[0].coedges()[0].edge).unwrap();
    let pc = e
        .pcurves()
        .iter()
        .find(|pc| pc.support_entity == cap.surface_entity())
        .unwrap();
    let c = p.curve(pc.curve_entity).unwrap();
    let mut poles = c.pole_entities().to_vec();
    poles.extend_from_slice(&c.pole_entities()[1..]);
    let ids = poles
        .iter()
        .map(|id| format!("#{id}"))
        .collect::<Vec<_>>()
        .join(",");
    let mut changed = argument_change(&text, c.entity(), 2, &format!("({ids})"));
    changed = argument_change(&changed, c.entity(), 6, "(12,11,12)");
    changed = argument_change(&changed, c.entity(), 7, "(0.,0.5,1.)");
    refused(&changed, "different finite cap spline basis");
    let tensor_face = p
        .faces()
        .iter()
        .find(|f| matches!(f.surface(), SurfaceIndex::Tensor(_)))
        .unwrap();
    let e = p.edge(tensor_face.loops()[0].coedges()[0].edge).unwrap();
    let pc = e
        .pcurves()
        .iter()
        .find(|c| c.support_entity == tensor_face.surface_entity())
        .unwrap();
    let uv = p.curve(pc.curve_entity).unwrap();
    let a = p.point(uv.pole_entities()[0]).unwrap();
    refused(
        &argument_change(&text, a.entity(), 1, "(0.5,0.)"),
        "partial UV border",
    );
    let ruling = p
        .edges()
        .iter()
        .find(|e| p.curve(e.curve_entity()).unwrap().degree() == 1)
        .unwrap();
    let c = p.curve(ruling.curve_entity()).unwrap();
    let reversed = c
        .pole_entities()
        .iter()
        .rev()
        .map(|id| format!("#{id}"))
        .collect::<Vec<_>>()
        .join(",");
    // Preserve EDGE_CURVE endpoint aliases too: the finite ruling reverses but its
    // parameter-space definition does not, so only incidence qualification refuses.
    let mut changed = argument_change(&text, c.entity(), 2, &format!("({reversed})"));
    changed = argument_change(
        &changed,
        ruling.entity(),
        1,
        &format!("#{}", ruling.vertex_entities()[1]),
    );
    changed = argument_change(
        &changed,
        ruling.entity(),
        2,
        &format!("#{}", ruling.vertex_entities()[0]),
    );
    refused(&changed, "reversed ruling without PCurve reversal");
    let changed = record_change(&text, tensor_face.loops()[0].entity(), |line| {
        let a = line.find("(#").unwrap();
        let prefix = &line[..a + 1];
        let entries = &line[a + 1..line.len() - 3];
        let mut ids = entries.split(',').collect::<Vec<_>>();
        ids.swap(0, 1);
        format!("{prefix}{}));", ids.join(","))
    });
    refused(&changed, "disconnected ordered wire");
}
#[test]
fn actual_literal_interior_height_and_extreme_numeric_enclosures_are_fail_closed() {
    let text = file(false);
    let p = parse(&text, Options::default(), None).unwrap();
    let pole = p
        .tensors()
        .iter()
        .flat_map(|t| t.rows())
        .flat_map(|r| r[1..r.len() - 1].iter())
        .map(|id| p.point(*id).unwrap())
        .find(|point| {
            let z = point.coordinates()[2].value();
            z > 0. && z < 4.
        })
        .unwrap();
    let xyz = pole.coordinates();
    let changed = argument_change(
        &text,
        pole.entity(),
        1,
        &format!("({},{},-1.)", xyz[0].spelling(), xyz[1].spelling()),
    );
    let pp = parse(&changed, Options::default(), None).unwrap();
    let b = pp.support_bindings(Options::default(), None).unwrap();
    assert!(matches!(
        b.material_bounds(Options::default(), None),
        Err(Error::Unsupported)
    ));
    // An unrelated retained 3D numeric leaf still participates in the actual
    // outward-cache preflight; finite parse is not finite interval authority.
    let changed=text.replacen("ENDSEC;\nEND-ISO-10303-21;","#4294967294=CARTESIAN_POINT('',(1.7976931348623157E308,0.,1.));\nENDSEC;\nEND-ISO-10303-21;",1);
    let pp = parse(&changed, Options::default(), None).unwrap();
    let b = pp.support_bindings(Options::default(), None).unwrap();
    assert!(matches!(
        b.material_bounds(Options::default(), None),
        Err(Error::Unsupported)
    ));
}
#[test]
fn actual_reference_dimensions_contexts_and_unrepresented_sources_are_rejected() {
    let text = file(false);
    let p = parse(&text, Options::default(), None).unwrap();
    let plane = &p.planes()[0];
    let edge = &p.edges()[0];
    let pc = edge.pcurves()[0];
    let vertex = &p.vertices()[0];
    let point2 = p
        .points()
        .iter()
        .find(|x| x.coordinates().len() == 2)
        .unwrap();
    let curve3 = p.curve(edge.curve_entity()).unwrap();
    let curve2 = p.curve(pc.curve_entity).unwrap();
    let context_record = p.entities[&pc.context]
        .iter()
        .find(|e| e.name == "GEOMETRIC_REPRESENTATION_CONTEXT")
        .unwrap();
    let mut cases = vec![
        argument_change(&text, vertex.entity, 0, ".T."),
        argument_change(&text, vertex.entity, 1, &format!("#{}", point2.entity())),
        argument_change(
            &text,
            plane.axis_entity(),
            1,
            &format!("#{}", point2.entity()),
        ),
        argument_change(&text, plane.direction_entities()[0], 0, ".T."),
        argument_change(&text, plane.direction_entities()[0], 1, "(0.,1.)"),
        argument_change(&text, plane.direction_entities()[0], 1, "(0.,'x',1.)"),
        argument_change(&text, pc.entity, 2, "#4294967295"),
        argument_change(&text, pc.representation, 2, "#4294967295"),
        argument_change(&text, pc.representation, 2, &format!("#{}", vertex.entity)),
        argument_change(&text, pc.entity, 0, ".T."),
    ];
    let poles = curve3
        .pole_entities()
        .iter()
        .enumerate()
        .map(|(i, id)| format!("#{}", if i == 1 { point2.entity() } else { *id }))
        .collect::<Vec<_>>()
        .join(",");
    cases.push(argument_change(
        &text,
        curve3.entity(),
        2,
        &format!("({poles})"),
    ));
    let t = &p.tensors()[0];
    // One dimensionally invalid tensor pole, leaving curve metadata intact.
    let mut rows = t.rows().map(|row| row.to_vec());
    rows[0][0] = point2.entity();
    let net = rows.map(|row| {
        format!(
            "({})",
            row.iter()
                .map(|id| format!("#{id}"))
                .collect::<Vec<_>>()
                .join(",")
        )
    });
    cases.push(argument_change(
        &text,
        t.entity(),
        3,
        &format!("({},{})", net[0], net[1]),
    ));
    let poles = vec![format!("#{}", point2.entity()); curve3.pole_entities().len()].join(",");
    cases.push(argument_change(
        &text,
        curve3.entity(),
        2,
        &format!("({poles})"),
    ));
    cases.push(argument_change(
        &text,
        pc.representation,
        1,
        &format!("(#{})", curve3.entity()),
    ));
    let source = context_record.source.unwrap();
    cases.push(text.replacen(
        source,
        &source.replace(
            "GEOMETRIC_REPRESENTATION_CONTEXT(2)",
            "GEOMETRIC_REPRESENTATION_CONTEXT()",
        ),
        1,
    ));
    cases.push(text.replacen(
        source,
        &source.replace(
            "GEOMETRIC_REPRESENTATION_CONTEXT(2)",
            "GEOMETRIC_REPRESENTATION_CONTEXT(4)",
        ),
        1,
    ));
    cases.push(text.replacen(
        source,
        &source.replace(
            "GEOMETRIC_REPRESENTATION_CONTEXT(2)",
            "GEOMETRIC_REPRESENTATION_CONTEXT(3)",
        ),
        1,
    ));
    let id = p
        .entities
        .iter()
        .find_map(|(&id, e)| {
            matches!(e.as_slice(),[e]if e.name=="MANIFOLD_SOLID_BREP").then_some(id)
        })
        .unwrap();
    let raw = text
        .lines()
        .find(|line| line.starts_with(&format!("#{id}=")))
        .unwrap();
    cases.push(text.replacen(
        "ENDSEC;\nEND-ISO-10303-21;",
        &format!(
            "{}\nENDSEC;\nEND-ISO-10303-21;",
            raw.replacen(&format!("#{id}="), "#4294967294=", 1)
        ),
        1,
    ));
    cases.push(record_change(&text, vertex.entity, |_| {
        format!(
            "#{}=(VERTEX_POINT('',#{})REPRESENTATION_ITEM(''));",
            vertex.entity, vertex.point_entity
        )
    }));
    cases.push(record_change(&text, vertex.entity, |s| {
        s.replace("VERTEX_POINT", "NOT_A_VERTEX")
    }));
    // The 2D curve is genuinely referenced by a finite 3D EDGE_CURVE.
    cases.push(argument_change(
        &text,
        edge.surface_curve_entity(),
        1,
        &format!("#{}", curve2.entity()),
    ));
    for (i, changed) in cases.iter().enumerate() {
        assert!(
            parse(changed, Options::default(), None).is_err(),
            "actual raw reference mutation {i}"
        );
    }
}
#[test]
fn actual_repeated_reversed_and_unused_finite_owners_cannot_earn_incidence() {
    let text = file(false);
    let p = parse(&text, Options::default(), None).unwrap();
    let cap = p
        .faces()
        .iter()
        .find(|f| matches!(f.surface(), SurfaceIndex::Plane(_)))
        .unwrap();
    let wire = &cap.loops()[0];
    let ids = refs(&p.entities[&wire.entity()][0].args[1]).unwrap();
    let doubled = ids
        .iter()
        .chain(&ids)
        .map(|id| format!("#{id}"))
        .collect::<Vec<_>>()
        .join(",");
    refused(
        &argument_change(&text, wire.entity(), 1, &format!("({doubled})")),
        "third literal coedge use",
    );
    let mut changed = text.clone();
    for (id, coedge) in ids.iter().zip(wire.coedges()) {
        changed = argument_change(&changed, *id, 4, if coedge.forward { ".F." } else { ".T." });
    }
    let reverse = ids
        .iter()
        .rev()
        .map(|id| format!("#{id}"))
        .collect::<Vec<_>>()
        .join(",");
    changed = argument_change(&changed, wire.entity(), 1, &format!("({reverse})"));
    refused(&changed, "same-direction shared cap uses");
    let tensor = &p.tensors()[0];
    let line = text
        .lines()
        .find(|line| line.starts_with(&format!("#{}=", tensor.entity())))
        .unwrap();
    let changed = text.replacen(
        "ENDSEC;\nEND-ISO-10303-21;",
        &format!(
            "{}\nENDSEC;\nEND-ISO-10303-21;",
            line.replacen(&format!("#{}=", tensor.entity()), "#4294967294=", 1)
        ),
        1,
    );
    refused(&changed, "unrepresented tensor support");
    let shell = p
        .entities
        .iter()
        .find_map(|(&id, e)| matches!(e.as_slice(),[e]if e.name=="CLOSED_SHELL").then_some(id))
        .unwrap();
    let selected = p
        .faces()
        .iter()
        .filter(|f| f.entity() != cap.entity())
        .map(|f| format!("#{}", f.entity()))
        .collect::<Vec<_>>()
        .join(",");
    refused(
        &argument_change(&text, shell, 1, &format!("({selected})")),
        "unrepresented cap support",
    );
}
#[test]
fn actual_literal_double_angular_cover_is_not_a_single_occupied_body() {
    let text = file(false);
    let p = parse(&text, Options::default(), None).unwrap();
    let mut points = std::collections::BTreeSet::<u32>::new();
    for face in p
        .faces()
        .iter()
        .filter(|f| matches!(f.surface(), SurfaceIndex::Plane(_)))
    {
        for wire in face.loops() {
            for co in wire.coedges() {
                let e = p.edge(co.edge).unwrap();
                let pc = e
                    .pcurves()
                    .iter()
                    .find(|pc| pc.support_entity == face.surface_entity())
                    .unwrap();
                points.extend(
                    p.curve(pc.curve_entity)
                        .unwrap()
                        .pole_entities()
                        .iter()
                        .copied(),
                );
            }
        }
    }
    let mut replacements = std::collections::BTreeMap::new();
    for point in p.points() {
        let c = point.coordinates();
        if c.len() != 3 && !points.contains(&point.entity()) {
            continue;
        }
        let [x, y] = [c[0].value(), c[1].value()];
        let r = x.portable_hypot(y);
        if r == 0. {
            continue;
        }
        let xx = (x * x - y * y) / r;
        let yy = 2. * x * y / r;
        let xyz = if c.len() == 3 {
            format!("({xx:.16E},{yy:.16E},{})", c[2].spelling())
        } else {
            format!("({xx:.16E},{yy:.16E})")
        };
        replacements.insert(
            point.entity(),
            format!("#{}=CARTESIAN_POINT('',{xyz});", point.entity()),
        );
    }
    let changed = text
        .lines()
        .map(|line| {
            let id = line
                .strip_prefix('#')
                .and_then(|s| s.split_once('='))
                .and_then(|(s, _)| s.parse::<u32>().ok());
            id.and_then(|id| replacements.get(&id))
                .map_or(line, |s| s.as_str())
        })
        .collect::<Vec<_>>()
        .join("\n");
    let pp = parse(&changed, Options::default(), None).unwrap();
    let b = pp.support_bindings(Options::default(), None).unwrap();
    assert!(matches!(
        b.material_bounds(Options::default(), None),
        Err(Error::Unsupported)
    ));
}
#[test]
fn actual_tensor_parameter_degrees_directions_and_cap_normals_are_not_guessed() {
    let text = file(false);
    let p = parse(&text, Options::default(), None).unwrap();
    let side = p
        .faces()
        .iter()
        .find(|f| matches!(f.surface(), SurfaceIndex::Tensor(_)))
        .unwrap();
    let row_edge = side.loops()[0]
        .coedges()
        .iter()
        .map(|c| p.edge(c.edge).unwrap())
        .find(|e| p.curve(e.curve_entity()).unwrap().degree() == 11)
        .unwrap();
    let row_pc = row_edge
        .pcurves()
        .iter()
        .find(|c| c.support_entity == side.surface_entity())
        .unwrap();
    let uv = p.curve(row_pc.curve_entity).unwrap();
    let a = p.point(uv.pole_entities()[0]).unwrap().coordinates();
    let pole = p.point(uv.pole_entities()[1]).unwrap();
    refused(
        &argument_change(
            &text,
            pole.entity(),
            1,
            &format!("({},{})", a[0].spelling(), a[1].spelling()),
        ),
        "zero-length side UV row",
    );
    let mut ids = vec![uv.pole_entities()[0]; 11];
    ids.push(uv.pole_entities()[1]);
    let mut changed = argument_change(&text, uv.entity(), 1, "11");
    changed = argument_change(
        &changed,
        uv.entity(),
        2,
        &format!(
            "({})",
            ids.iter()
                .map(|id| format!("#{id}"))
                .collect::<Vec<_>>()
                .join(",")
        ),
    );
    changed = argument_change(&changed, uv.entity(), 6, "(12,12)");
    refused(&changed, "nonlinear tensor PCurve");
    let ruling = p
        .edges()
        .iter()
        .find(|e| p.curve(e.curve_entity()).unwrap().degree() == 1)
        .unwrap();
    let pc = ruling.pcurves()[0];
    let uv = p.curve(pc.curve_entity).unwrap();
    let a = p.point(uv.pole_entities()[0]).unwrap();
    let c = a.coordinates();
    refused(
        &argument_change(
            &text,
            a.entity(),
            1,
            &format!("({},{}.)", c[0].spelling(), 1 - u8::from(c[1].is_one())),
        ),
        "slanted UV ruling",
    );
    let reversed = uv
        .pole_entities()
        .iter()
        .rev()
        .map(|id| format!("#{id}"))
        .collect::<Vec<_>>()
        .join(",");
    refused(
        &argument_change(&text, uv.entity(), 2, &format!("({reversed})")),
        "opposed ruling UV parameter",
    );
    let cap = p
        .faces()
        .iter()
        .find(|f| matches!(f.surface(), SurfaceIndex::Plane(_)))
        .unwrap();
    refused(
        &argument_change(
            &text,
            cap.entity(),
            3,
            if cap.same_sense() { ".F." } else { ".T." },
        ),
        "opposed cap material normal",
    );
    let wire = &side.loops()[0];
    let ids = refs(&p.entities[&wire.entity()][0].args[1]).unwrap();
    let doubled = ids
        .iter()
        .chain(&ids)
        .map(|id| format!("#{id}"))
        .collect::<Vec<_>>()
        .join(",");
    refused(
        &argument_change(&text, wire.entity(), 1, &format!("({doubled})")),
        "repeated tensor domain border",
    );
}
#[test]
fn actual_open_side_and_fully_collapsed_patch_cannot_be_closed_material() {
    let text = file(false);
    let p = parse(&text, Options::default(), None).unwrap();
    let side = p
        .faces()
        .iter()
        .find(|f| matches!(f.surface(), SurfaceIndex::Tensor(_)))
        .unwrap();
    let shell = p
        .entities
        .iter()
        .find_map(|(&id, e)| matches!(e.as_slice(),[e]if e.name=="CLOSED_SHELL").then_some(id))
        .unwrap();
    let selected = p
        .faces()
        .iter()
        .filter(|f| f.entity() != side.entity())
        .map(|f| format!("#{}", f.entity()))
        .collect::<Vec<_>>()
        .join(",");
    let changed = argument_change(&text, shell, 1, &format!("({selected})"));
    let prefix = format!("#{}=", side.surface_entity());
    let changed = changed
        .lines()
        .filter(|line| !line.starts_with(&prefix))
        .collect::<Vec<_>>()
        .join("\n");
    refused(&changed, "open finite side shell");
    let t = p
        .tensors()
        .iter()
        .find(|t| {
            t.rows().iter().all(|r| {
                [r[0], *r.last().unwrap()].iter().all(|id| {
                    let z = p.point(*id).unwrap().coordinates()[2].value();
                    z > 0. && z < 4.
                })
            })
        })
        .unwrap();
    let [a, b] = t.rows();
    let aliases = [(b[0], a[0]), (*b.last().unwrap(), *a.last().unwrap())];
    let changed = references_changed(&text, &aliases);
    refused(&changed, "two collapsed source rows");
}
fn references_changed(text: &str, aliases: &[(u32, u32)]) -> String {
    let bytes = text.as_bytes();
    let mut out = String::with_capacity(text.len());
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'#' {
            let start = i + 1;
            let mut end = start;
            while end < bytes.len() && bytes[end].is_ascii_digit() {
                end += 1;
            }
            let id = text[start..end].parse::<u32>().unwrap();
            let mapped = if bytes.get(end) == Some(&b'=') {
                id
            } else {
                aliases
                    .iter()
                    .find_map(|&(a, b)| (a == id).then_some(b))
                    .unwrap_or(id)
            };
            out.push_str(&format!("#{mapped}"));
            i = end;
        } else {
            out.push(bytes[i] as char);
            i += 1;
        }
    }
    out
}
#[test]
fn actual_coincident_cartesian_points_do_not_merge_distinct_parameter_corner_vertices() {
    let text = file(false);
    let p = parse(&text, Options::default(), None).unwrap();
    let face = p
        .faces()
        .iter()
        .find(|f| {
            matches!(f.surface(), SurfaceIndex::Tensor(_)) && f.loops()[0].coedges().len() == 4
        })
        .unwrap();
    let wire = &face.loops()[0];
    let ids = refs(&p.entities[&wire.entity()][0].args[1]).unwrap();
    let vertex_ids = wire
        .coedges()
        .iter()
        .map(|co| p.edge(co.edge).unwrap().vertex_entities()[usize::from(!co.forward)])
        .collect::<Vec<_>>();
    assert_eq!(
        vertex_ids
            .iter()
            .copied()
            .collect::<std::collections::BTreeSet<_>>()
            .len(),
        4
    );
    let point = p.vertex(vertex_ids[0]).unwrap().point_entity;
    let aliases = vertex_ids
        .iter()
        .map(|id| (p.vertex(*id).unwrap().point_entity, point))
        .collect::<Vec<_>>();
    let mut changed = references_changed(&text, &aliases);
    let order = [0, 2, 1, 3];
    for (i, &at) in order.iter().enumerate() {
        let co = wire.coedges()[at];
        let e = p.edge(co.edge).unwrap();
        let [start, end] = [vertex_ids[i], vertex_ids[(i + 1) % 4]];
        let ends = if co.forward {
            [start, end]
        } else {
            [end, start]
        };
        for (j, &id) in ends.iter().enumerate() {
            if e.vertex_entities()[j] != id {
                changed = argument_change(&changed, co.edge, j + 1, &format!("#{id}"));
            }
        }
    }
    let uses = order.map(|i| format!("#{}", ids[i])).join(",");
    changed = argument_change(&changed, wire.entity(), 1, &format!("({uses})"));
    let shell = p
        .entities
        .iter()
        .find_map(|(&id, e)| matches!(e.as_slice(),[e]if e.name=="CLOSED_SHELL").then_some(id))
        .unwrap();
    let faces = std::iter::once(face)
        .chain(p.faces().iter().filter(|f| f.entity() != face.entity()))
        .map(|f| format!("#{}", f.entity()))
        .collect::<Vec<_>>()
        .join(",");
    changed = argument_change(&changed, shell, 1, &format!("({faces})"));
    refused(
        &changed,
        "different raw corner vertices despite identical actual Cartesian aliases",
    );
}
