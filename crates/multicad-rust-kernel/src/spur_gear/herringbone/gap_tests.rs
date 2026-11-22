//! Gapped double-helical chevrons: two opposite-hand tooth bands joined by
//! the root core through the central gap.
use super::tests::fixture;
use super::*;
use crate::{
    display::{MeshError, Tessellation},
    manufacturing::certify,
};
use std::collections::{BTreeMap, BTreeSet};

const MINIMUM: f64 = 1e-5;
const MAXIMUM: f64 = 1e7;
const OPTIONS: Tessellation = Tessellation {
    deflection: 0.005,
    angle: 0.15,
    triangle_limit: 400_000,
};

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b} ±{tolerance}");
}

fn gapped(teeth: u32, width: f64, gap: f64, bore: f64) -> Herringbone {
    let (mut spec, twist) = fixture(teeth);
    spec.face_width_mm = width;
    spec.bore_radius_mm = bore;
    Herringbone::new_gapped(spec, twist, gap, MINIMUM, MAXIMUM).unwrap()
}

/// Every face boundary is a set of closed, consistently oriented loops,
/// and the shell's Euler characteristic matches one handle per bore.
fn assert_closed_shell(gear: &Herringbone, topology: &Topology) {
    let mut balance: BTreeMap<(u32, u32), i64> = BTreeMap::new();
    let mut used = BTreeSet::new();
    for edge in &topology.edges {
        let [a, b] = edge.vertices;
        for (face, [from, to]) in edge.faces.into_iter().zip([[a, b], [b, a]]) {
            used.insert(face);
            *balance.entry((face, from)).or_default() += 1;
            *balance.entry((face, to)).or_default() -= 1;
        }
        for (fraction, vertex) in [0., 1.].into_iter().zip(edge.vertices) {
            for (actual, expected) in edge
                .curve
                .point(gear, fraction)
                .into_iter()
                .zip(topology.vertices_mm[vertex as usize])
            {
                close(actual, expected, 1e-11);
            }
        }
    }
    assert!(balance.values().all(|&b| b == 0), "open face loop");
    assert_eq!(used.len(), topology.faces.len());
    // Caps around a bore have one inner loop each; hub caps likewise.
    let bore = gear.profile().bore().is_some();
    let holes = if bore { 2 } else { 0 };
    let chi = topology.vertices_mm.len() as i64 - topology.edges.len() as i64
        + topology.faces.len() as i64
        - holes;
    assert_eq!(chi, if bore { 0 } else { 2 });
}

#[test]
fn gapped_shell_names_the_core_and_every_tooth_shoulder() {
    for (teeth, bore) in [(24, 4.), (36, 0.), (24, 0.)] {
        let gear = gapped(teeth, 10., 3., bore);
        assert!(gear.is_gapped());
        assert_eq!(gear.gap_mm(), 3.);
        assert_eq!(gear.band_length_mm(), 3.5);
        assert_eq!(gear.shoulder_heights_mm(), Some([3.5, 6.5]));
        assert_eq!(gear.teeth(), teeth);
        let topology = gear.topology();
        let count = gear.profile().outline().len();
        let bore_faces = usize::from(bore > 0.);
        assert_eq!(
            topology.faces.len(),
            2 + 2 * count + bore_faces + 1 + 2 * teeth as usize
        );
        assert_eq!(topology.vertices_mm.len(), 4 * count + 2 * bore_faces);
        assert_eq!(
            topology.edges.len(),
            6 * count + 2 * teeth as usize + 1 + 3 * bore_faces
        );
        let core = gear.gap_faces_start() as usize;
        assert_eq!(topology.faces[core], Face::RootCore);
        assert_eq!(topology.faces[core].role(), "gear/root.core");
        assert_eq!(
            topology.faces[core + 1],
            Face::Shoulder {
                upper: false,
                tooth: 0
            }
        );
        assert_eq!(topology.faces[core + 1].role(), "gear/gap.lower");
        assert_eq!(
            topology.faces[core + 1 + teeth as usize].role(),
            "gear/gap.upper"
        );
        assert_closed_shell(&gear, &topology);
        // Under-tooth arcs lie on the root circle at both shoulders.
        for edge in &topology.edges {
            if let Curve::UnderTooth { ring, .. } = edge.curve {
                let z = gear.ring_height_mm(ring);
                for fraction in [0., 0.3, 1.] {
                    let [x, y, actual] = edge.curve.point(&gear, fraction);
                    close(x.hypot(y), spec_root(teeth), 1e-11);
                    assert_eq!(actual, z);
                }
            }
            if let Curve::RootSeam { .. } = edge.curve {
                let [a, b] = [0., 1.].map(|f| edge.curve.point(&gear, f));
                assert_eq!([a[2], b[2]], [3.5, 6.5]);
                close(a[0], b[0], 1e-12);
                close(a[1], b[1], 1e-12);
            }
        }
    }
}

fn spec_root(teeth: u32) -> f64 {
    fixture(teeth).0.root_radius_mm
}

#[test]
fn each_band_keeps_its_hand_and_both_shoulders_share_one_phase() {
    let gear = gapped(24, 10., 3., 4.);
    let twist = gear.twist_rad_per_mm();
    for z in [0., 1., 3.5] {
        assert_eq!(gear.phase_at(z), twist * z);
    }
    for z in [3.5, 4., 5., 6.5] {
        assert_eq!(gear.phase_at(z), twist * 3.5, "gap at {z}");
    }
    for z in [6.5, 8., 10.] {
        assert_eq!(gear.phase_at(z), twist * (10. - z));
    }
    assert_eq!(gear.band_twist_rad_per_mm(Band::Lower), twist);
    assert_eq!(gear.band_twist_rad_per_mm(Band::Upper), -twist);
    assert_eq!(gear.band_interval_mm(Band::Lower), [0., 3.5]);
    assert_eq!(gear.band_interval_mm(Band::Upper), [6.5, 10.]);
    // Continuous chevrons keep their apex station.
    let (spec, twist) = fixture(24);
    let continuous = Herringbone::new(spec, twist, MINIMUM, MAXIMUM).unwrap();
    assert_eq!(continuous.rings(), [Ring::Start, Ring::Apex, Ring::End]);
    assert_eq!(continuous.gap_face_count(), 0);
    assert_eq!(continuous.shoulder_heights_mm(), None);
    assert_eq!(
        gear.rings(),
        [
            Ring::Start,
            Ring::LowerShoulder,
            Ring::UpperShoulder,
            Ring::End
        ]
    );
}

#[test]
fn analytic_volume_bounds_and_certified_mesh_agree() {
    for (teeth, bore, gap) in [(24, 4., 3.), (36, 0., 1.), (24, 0., 8.)] {
        let gear = gapped(teeth, 10., gap, bore);
        let (spec, _) = fixture(teeth);
        let band = (10. - gap) / 2.;
        let area = gear.profile().area_mm2();
        let core = PI * (spec.root_radius_mm.powi(2) - bore * bore);
        close(gear.volume_mm3(), 2. * band * area + gap * core, 1e-9);
        let mesh = gear.mesh(OPTIONS).unwrap();
        let proof = certify(&mesh).unwrap();
        // Chords lie inside the exact boundary, within deflection × area.
        let surface = 2. * PI * spec.tip_radius_mm * 10. * 3.;
        assert!(proof.signed_volume_mm3 <= gear.volume_mm3() + 1e-6);
        assert!(
            gear.volume_mm3() - proof.signed_volume_mm3 < OPTIONS.deflection * surface,
            "{} vs {}",
            proof.signed_volume_mm3,
            gear.volume_mm3()
        );
        let topology = gear.topology();
        let faces: BTreeSet<_> = mesh.face_indices.iter().copied().collect();
        assert_eq!(faces.len(), topology.faces.len());
        let [lo, hi] = gear.bounds();
        assert_eq!([lo[2], hi[2]], [0., 10.]);
        let mut extreme = [[f64::INFINITY; 3], [f64::NEG_INFINITY; 3]];
        for vertex in mesh.vertices.as_chunks::<3>().0 {
            for i in 0..3 {
                extreme[0][i] = extreme[0][i].min(vertex[i]);
                extreme[1][i] = extreme[1][i].max(vertex[i]);
                assert!(vertex[i] >= lo[i] - 1e-9 && vertex[i] <= hi[i] + 1e-9);
            }
        }
        for i in 0..2 {
            close(extreme[0][i], lo[i], 0.02);
            close(extreme[1][i], hi[i], 0.02);
        }
        // Gap triangles lie on the core, a shoulder plane or the bore.
        let core_face = gear.gap_faces_start();
        for (triangle, &face) in mesh
            .indices
            .as_chunks::<3>()
            .0
            .iter()
            .zip(&mesh.face_indices)
        {
            let points = triangle.map(|i| {
                let i = 3 * i as usize;
                [mesh.vertices[i], mesh.vertices[i + 1], mesh.vertices[i + 2]]
            });
            if face == core_face {
                for p in points {
                    close(p[0].hypot(p[1]), spec.root_radius_mm, 1e-9);
                    assert!((band..=10. - band).contains(&p[2]));
                }
            } else if face > core_face && face <= core_face + 2 * teeth {
                let z = if face <= core_face + teeth {
                    band
                } else {
                    10. - band
                };
                for p in points {
                    assert_eq!(p[2], z);
                    assert!(p[0].hypot(p[1]) >= spec.root_radius_mm - 1e-9);
                }
            }
        }
    }
}

#[test]
fn hubs_and_bores_close_through_the_gap() {
    for bore in [0., 4.] {
        let plain = gapped(24, 10., 3., bore);
        for span in [[-2., 2.], [8., 12.], [-2., 12.], [4., 6.]] {
            let gear = plain
                .clone()
                .with_hub(
                    super::super::Hub {
                        radius_mm: 12.,
                        axial_start_mm: span[0],
                        axial_end_mm: span[1],
                    },
                    MINIMUM,
                    MAXIMUM,
                )
                .unwrap();
            let topology = gear.topology();
            let mesh = gear.mesh(OPTIONS).unwrap();
            let proof = certify(&mesh).unwrap();
            assert!(proof.signed_volume_mm3 <= gear.volume_mm3() + 1e-6);
            assert!(gear.volume_mm3() - proof.signed_volume_mm3 < 2.);
            let faces: BTreeSet<_> = mesh.face_indices.iter().copied().collect();
            assert_eq!(faces.len(), topology.faces.len());
            let hubs = gear.hub_extensions();
            for extension in hubs.into_iter().flatten() {
                assert!(extension.wall_face >= gear.gap_faces_start() + gear.gap_face_count());
                assert!(matches!(
                    topology.faces[extension.wall_face as usize],
                    Face::HubWall { .. }
                ));
            }
            if span == [4., 6.] {
                assert_eq!(mesh.vertices, plain.mesh(OPTIONS).unwrap().vertices);
            }
            assert_eq!(
                gear.mesh(Tessellation {
                    triangle_limit: mesh.indices.len() / 3 - 1,
                    ..OPTIONS
                })
                .err(),
                Some(MeshError::BudgetExceeded)
            );
        }
    }
}

#[test]
fn invalid_gap_and_band_combinations_are_refused() {
    let (spec, twist) = fixture(24);
    let spec = Spec {
        face_width_mm: 10.,
        ..spec
    };
    for gap in [0., -1., f64::NAN, f64::INFINITY, 1e-6, 10., 12., 10. - 1e-6] {
        assert_eq!(
            Herringbone::new_gapped(spec, twist, gap, MINIMUM, MAXIMUM),
            Err(Error::InvalidDimensions),
            "{gap}"
        );
    }
    // A vanishing band twist cannot resolve a hand.
    assert_eq!(
        Herringbone::new_gapped(spec, 1e-9, 3., MINIMUM, MAXIMUM),
        Err(Error::InvalidDimensions)
    );
    assert_eq!(
        Herringbone::new_gapped(spec, 0., 3., MINIMUM, MAXIMUM),
        Err(Error::InvalidDimensions)
    );
    assert_eq!(
        Herringbone::new_gapped(Spec { teeth: 0, ..spec }, twist, 3., MINIMUM, MAXIMUM),
        Err(Error::Limit)
    );
}

#[test]
fn support_encloses_dense_band_sections_and_the_core_in_oblique_directions() {
    for teeth in [24, 36] {
        let gear = gapped(teeth, 10., 3., 4.);
        let (spec, twist) = fixture(teeth);
        let root = spec.root_radius_mm;
        for [dx, dy, dz] in [
            [1., 0.2, 0.12],
            [0.1, -1., -0.21],
            [-1., 0.9, 0.],
            [0.5, 0.5, 0.6],
            [0.01, 0., 1.],
            [0., 0.01, -1.],
        ] {
            let support = gear.support([dx, dy, dz]);
            let intervals = 2000;
            let radial = dx.hypot(dy);
            let sampled = (0..=intervals)
                .map(|i| {
                    let z = i as f64 * 10. / intervals as f64;
                    if (3.5..=6.5).contains(&z) && !(z == 3.5 || z == 6.5) {
                        radial * root + dz * z
                    } else {
                        let [x, y] = rotate([dx, dy], -gear.phase_at(z));
                        gear.profile.support([x, y, 0.]) + dz * z
                    }
                })
                .fold(f64::NEG_INFINITY, f64::max);
            let maximum_sampling_gap = 10. / (2. * intervals as f64)
                * (radial * spec.tip_radius_mm * twist.abs() + dz.abs());
            assert!(support + 1e-12 >= sampled, "{support} < {sampled}");
            assert!(support - sampled <= maximum_sampling_gap + 1e-12);
        }
        // Axial supports are exactly the end caps.
        assert_eq!(gear.support([0., 0., 1.]), 10.);
        assert_eq!(gear.support([0., 0., -1.]), 0.);
    }
}
