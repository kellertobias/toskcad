use super::*;
use crate::{CurvedEdge2d, CurvedProfile};

fn policy() -> GeometryPolicy {
    GeometryPolicy::default()
}

fn origin(feature: EntityId, role: &str) -> Origin {
    Origin {
        source: feature,
        role: role.into(),
    }
}

/// A sweep of `outer` less `holes` in the `y, z` world plane along `+x`
/// over `[-length / 2, length / 2]`.
fn sweep(outer: Vec<CurvedEdge2d>, holes: Vec<Vec<CurvedEdge2d>>, length: f64) -> Solid {
    let feature = EntityId::new();
    let count = outer.len() + holes.iter().map(Vec::len).sum::<usize>();
    let origins: Vec<Origin> = (0..count)
        .map(|i| origin(feature, &format!("edge.{i}")))
        .collect();
    Solid::rust_extrude_normal_profile_named(
        &CurvedProfile {
            frame: PlaneFrame {
                origin: [-length / 2., 0., 0.],
                u: [0., 1., 0.],
                v: [0., 0., 1.],
            },
            outer,
            holes,
        },
        length,
        feature,
        &origins,
        policy(),
    )
    .unwrap()
}

fn polygon(points: &[[f64; 2]]) -> Vec<CurvedEdge2d> {
    (0..points.len())
        .map(|i| CurvedEdge2d::Line {
            a: points[i],
            b: points[(i + 1) % points.len()],
        })
        .collect()
}

fn square(half: f64) -> Vec<CurvedEdge2d> {
    polygon(&[[-half, -half], [half, -half], [half, half], [-half, half]])
}

fn tool(steps: &[(f64, f64, f64)]) -> TransverseBoreTool {
    let feature = EntityId::new();
    TransverseBoreTool {
        origin_mm: [0.; 3],
        axis: [0., 0., 1.],
        steps: steps
            .iter()
            .map(|&(start_mm, end_mm, radius_mm)| BoreStep {
                start_mm,
                end_mm,
                radius_mm,
            })
            .collect(),
        wall: origin(feature, "bore.wall"),
        shoulder: origin(feature, "bore.shoulder"),
    }
}

/// Exact classification of a small box about `p` by the retained boundary.
fn material(solid: &Solid, p: [f64; 3]) -> bool {
    let e = 1e-3;
    let body = solid.imported_brep_body().expect("a retained boundary");
    match body
        .classify_box([p.map(|v| v - e), p.map(|v| v + e)], 0.)
        .unwrap()
    {
        brep::BoxSide::Inside => true,
        brep::BoxSide::Outside => false,
        brep::BoxSide::Unresolved => panic!("unresolved at {p:?}"),
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected} within {tolerance}"
    );
}

#[test]
fn through_and_stepped_blind_bores_in_a_box_keep_closed_form_volumes() {
    let stock = sweep(square(20.), vec![], 100.);
    let full = 40. * 40. * 100.;
    let through = stock
        .rust_transverse_bore(&tool(&[(-25., 25., 5.)]), policy())
        .unwrap();
    close(through.volume_mm3().unwrap(), full - PI * 25. * 40., 1e-6);
    assert_eq!(through.solid_count().unwrap(), 1);
    // A counterbore entering the top face over a deeper pilot: a shoulder
    // at z = 4 and a floor at z = -13.
    let stepped = stock
        .rust_transverse_bore(&tool(&[(4., 20.02, 10.), (-13., 20.02, 3.5)]), policy())
        .unwrap();
    let removed = PI * 100. * 16. + PI * 3.5 * 3.5 * 17.;
    close(stepped.volume_mm3().unwrap(), full - removed, 1e-6);
    // The retained boundary decides material exactly.
    let inside = |p: [f64; 3]| material(&stepped, p);
    assert!(!inside([0., 9., 10.]));
    assert!(inside([0., 11., 10.]));
    assert!(!inside([0., 3., 0.]));
    assert!(inside([0., 4., 0.]));
    assert!(inside([0., 0., -14.]));
}

#[test]
fn a_tool_reaching_the_prism_end_is_refused() {
    let stock = sweep(square(20.), vec![], 20.);
    let error = stock
        .rust_transverse_bore(&tool(&[(-25., 25., 10.)]), policy())
        .err()
        .unwrap();
    assert_eq!(error.kind, GeometryErrorKind::Unsupported, "{error}");
}

/// Removed volume of a perpendicular half cross-bore of radius `r` from the
/// tool axis through the wall of a round bar of radius `big`.
fn half_cross_bore(r: f64, big: f64) -> f64 {
    let n = 4096;
    let h = PI / n as f64;
    (0..n)
        .map(|k| {
            let t = -FRAC_PI_2 + (k as f64 + 0.5) * h;
            let y = r * t.sin();
            2. * r * t.cos() * (big * big - y * y).sqrt() * r * t.cos()
        })
        .sum::<f64>()
        * h
}

#[test]
fn stepped_cross_bore_of_a_round_barrel_matches_its_quartic_volume() {
    // The universal-end barrel: a round bar along x with a head seat and a
    // screw passage meeting at its axis.
    let radius = 9.98;
    let length = 15.98;
    let bar = sweep(
        vec![CurvedEdge2d::Circle {
            center: [0., 0.],
            radius_mm: radius,
        }],
        vec![],
        length,
    );
    let cut = bar
        .rust_transverse_bore(&tool(&[(-12., 0., 7.04), (-0.02, 11.98, 4.1)]), policy())
        .unwrap();
    let expected = PI * radius * radius * length
        - half_cross_bore(7.04, radius)
        - half_cross_bore(4.1, radius);
    close(cut.volume_mm3().unwrap(), expected, 1e-6);
    assert_eq!(cut.solid_count().unwrap(), 1);
    let inside = |p: [f64; 3]| material(&cut, p);
    assert!(!inside([0., 0., -5.]));
    assert!(!inside([0., 6., -5.]));
    assert!(inside([0., 6., 5.]));
    assert!(inside([0., 8., -5.]));
}

#[test]
fn a_plain_cylinder_takes_the_same_cross_bore() {
    let feature = EntityId::new();
    let bar = crate::geometry::construction::cylinder(
        PlaneFrame {
            origin: [-7.99, 0., 0.],
            u: [0., 1., 0.],
            v: [0., 0., 1.],
        },
        9.98,
        15.98,
        feature,
        policy(),
    )
    .unwrap();
    let cut = bar
        .rust_transverse_bore(&tool(&[(-12., 0., 7.04), (-0.02, 11.98, 4.1)]), policy())
        .unwrap();
    let expected =
        PI * 9.98 * 9.98 * 15.98 - half_cross_bore(7.04, 9.98) - half_cross_bore(4.1, 9.98);
    close(cut.volume_mm3().unwrap(), expected, 1e-6);
}

const LINE8_FLANK: [[f64; 2]; 4] = [[9.25, 4.5], [7.75, 8.], [5.75, 10.], [2.75, 12.25]];

/// The item Line 8 40x40 section about its centre, with its 6.8 mm core.
fn line8() -> (Vec<CurvedEdge2d>, Vec<Vec<CurvedEdge2d>>) {
    let (w, r, m, lip) = (40., 4., 4., 4.5);
    let c = w / 2.;
    let mut side = vec![[r, 0.], [c - m, 0.], [c - m, lip]];
    side.extend(LINE8_FLANK.map(|[half, depth]| [c - half, depth]));
    side.extend(
        LINE8_FLANK
            .iter()
            .rev()
            .map(|&[half, depth]| [c + half, depth]),
    );
    side.extend([[c + m, lip], [c + m, 0.], [w - r, 0.]]);
    let rotate = |[x, y]: [f64; 2]| [w - y, x];
    let shift = |[x, y]: [f64; 2]| [x - c, y - c];
    let mut center = [w - r, r];
    let mut edges = Vec::new();
    for _ in 0..4 {
        for pair in side.windows(2) {
            edges.push(CurvedEdge2d::Line {
                a: shift(pair[0]),
                b: shift(pair[1]),
            });
        }
        edges.push(CurvedEdge2d::Arc {
            center: shift(center),
            a: shift(side[side.len() - 1]),
            b: shift(rotate(side[0])),
            clockwise: false,
        });
        for point in &mut side {
            *point = rotate(*point);
        }
        center = rotate(center);
    }
    let core = vec![CurvedEdge2d::Circle {
        center: [0., 0.],
        radius_mm: 3.4,
    }];
    (edges, vec![core])
}

fn line8_area() -> f64 {
    let groove = 8. * 4.5
        + LINE8_FLANK
            .windows(2)
            .map(|p| (p[0][0] + p[1][0]) * (p[1][1] - p[0][1]))
            .sum::<f64>();
    1600. - (4. - PI) * 16. - PI * 3.4 * 3.4 - 4. * groove
}

/// Vertical material chords of the section at `y`, as `z` intervals.
fn section_chords(edges: &[CurvedEdge2d], holes: &[Vec<CurvedEdge2d>], y: f64) -> Vec<[f64; 2]> {
    let mut zs = Vec::new();
    for edge in edges.iter().chain(holes.iter().flatten()) {
        match *edge {
            CurvedEdge2d::Line { a, b } => {
                if (a[0] > y) != (b[0] > y) {
                    zs.push(a[1] + (y - a[0]) * (b[1] - a[1]) / (b[0] - a[0]));
                }
            }
            // Corner arcs lie at |y| >= 16, beyond every sampled chord.
            CurvedEdge2d::Arc { center, .. } => assert!(center[0].abs() >= 16. && y.abs() <= 10.),
            CurvedEdge2d::Circle { center, radius_mm } => {
                let dy = y - center[0];
                if dy.abs() < radius_mm {
                    let root = (radius_mm * radius_mm - dy * dy).sqrt();
                    zs.extend([center[1] - root, center[1] + root]);
                }
            }
            CurvedEdge2d::BSpline(_) => unreachable!(),
        }
    }
    zs.sort_by(f64::total_cmp);
    zs.chunks(2).map(|p| [p[0], p[1]]).collect()
}

fn gauss(f: &dyn Fn(f64) -> f64, a: f64, b: f64) -> f64 {
    // Composite 5-point Gauss-Legendre.
    const NODES: [(f64, f64); 5] = [
        (0., 0.5688888888888889),
        (-0.5384693101056831, 0.47862867049936647),
        (0.5384693101056831, 0.47862867049936647),
        (-0.906179845938664, 0.23692688505618908),
        (0.906179845938664, 0.23692688505618908),
    ];
    let n = 400;
    let h = (b - a) / n as f64;
    (0..n)
        .map(|k| {
            let c = a + h * (k as f64 + 0.5);
            NODES
                .iter()
                .map(|(x, w)| w * f(c + 0.5 * h * x))
                .sum::<f64>()
                * 0.5
                * h
        })
        .sum()
}

#[test]
fn double_entry_access_bores_across_a_line8_section_match_an_independent_oracle() {
    let (outer, holes) = line8();
    let stock = sweep(outer.clone(), holes.clone(), 100.);
    // Two counterbores over pilots from opposite groove faces along z: the
    // pilots cross the 6.8 mm core at 7 mm, and the side-groove flank
    // vertices at |y| = 10 touch the counterbore silhouettes.
    let mut tool = tool(&[
        (-0.02, 16., 10.),
        (-0.02, 33., 3.5),
        (24., 40.02, 10.),
        (7., 40.02, 3.5),
    ]);
    tool.origin_mm = [0., 0., -20.];
    let cut = stock.rust_transverse_bore(&tool, policy()).unwrap();
    assert_eq!(cut.solid_count().unwrap(), 1);
    // Slice along y instead of the builder's z: each material chord loses
    // the tool's x widths in the bands it crosses.
    let length = |chords: &[[f64; 2]], range: [f64; 2]| {
        chords
            .iter()
            .map(|c| (c[1].min(range[1]) - c[0].max(range[0])).max(0.))
            .sum::<f64>()
    };
    let removed_at = |y: f64| {
        let chords = section_chords(&outer, &holes, y);
        let half = |r: f64| ((r * r - y * y).max(0.)).sqrt();
        let all = length(&chords, [-20.02, 20.02]);
        let outer_bands = length(&chords, [-20.02, -4.]) + length(&chords, [4., 20.02]);
        2. * half(3.5) * all + 2. * (half(10.) - half(3.5)) * outer_bands
    };
    let mut cuts = vec![-10., -3.5, -3.4, 3.4, 3.5, 10.];
    for edge in &outer {
        if let CurvedEdge2d::Line { a, .. } = edge
            && a[0].abs() < 10.
        {
            cuts.push(a[0]);
        }
    }
    cuts.sort_by(f64::total_cmp);
    cuts.dedup();
    let removed: f64 = cuts
        .windows(2)
        .map(|w| gauss(&removed_at, w[0], w[1]))
        .sum();
    let expected = line8_area() * 100. - removed;
    let volume = cut.volume_mm3().unwrap();
    close(volume, expected, 1e-6 * expected);
    assert!(removed > 4000., "{removed}");
    let inside = |p: [f64; 3]| material(&cut, p);
    // Counterbore, its silhouette pinch, the lip beyond it and the floor.
    assert!(!inside([0., 9., -12.]));
    assert!(inside([0., 11., -16.]));
    assert!(!inside([0., 9.9, -5.]));
    assert!(inside([0., 3.6, -3.]));
    // The pilot crosses the core: the 0.1 mm ring between them is gone in
    // the tool and kept beside it.
    assert!(!inside([0., 3.45, 0.]));
    assert!(!inside([-11., 0., 0.]));
    assert!(inside([-11., 3.6, 0.]));
    assert!(inside([2., 3.45, 0.5]));
    // The display mesh, STEP exchange and a drawing section read the same
    // retained boundary.
    let mesh = cut.mesh(MeshOptions::default()).unwrap();
    let triangles = mesh.indices.chunks(3).map(|t| {
        let p = |i: u32| {
            let i = i as usize * 3;
            [mesh.vertices[i], mesh.vertices[i + 1], mesh.vertices[i + 2]]
        };
        [p(t[0]), p(t[1]), p(t[2])]
    });
    let mesh_volume: f64 = triangles
        .map(|[a, b, c]| {
            (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0]))
                / 6.
        })
        .sum();
    close(mesh_volume, volume, 1e-2 * volume);
    let path = std::env::temp_dir().join(format!("transverse-bore-{}.step", EntityId::new()));
    cut.write_step(&path).unwrap();
    let read = Solid::read_step(&path, policy()).unwrap();
    std::fs::remove_file(&path).unwrap();
    close(read.volume_mm3().unwrap(), volume, 1e-6 * volume);
    let section = crate::drawing::project_drawing(
        &[crate::drawing::DrawingSolid {
            solid: &cut,
            hatch: None,
        }],
        crate::drawing::DrawingProjection {
            frame: PlaneFrame {
                origin: [0.; 3],
                u: [0., 1., 0.],
                v: [0., 0., 1.],
            },
            section_offset_mm: Some(-1.),
            section_only: true,
            chord_tolerance_mm: 0.01,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(!section.is_empty());
}
