use super::*;

#[test]
fn drive_ring_requires_distinct_corner_stations_before_allocating_geometry() {
    let thread = crate::helical_thread::Thread::new(crate::helical_thread::Spec {
        major_radius_mm: 1.5,
        pitch_mm: 0.5,
        length_mm: 0.5,
        depth_mm: 0.2,
        crest_width_mm: 0.0625,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: false,
        root_radius_mm: None,
    })
    .unwrap();
    let options = Tessellation {
        deflection: 0.1,
        angle: 1.,
        triangle_limit: 10000,
    };
    let strips = phase_strips(&thread, TAU / 24., options).unwrap();
    let mut builder = Builder::new(&thread, 24, strips, options);
    let (s, c) = 0.37_f64.portable_sin_cos();
    let points = crate::headed_screw::hex_vertices(2.).map(|[x, y]| [c * x - s * y, s * x + c * y]);
    let drive = crate::hex_drive::HexDrive::from_retained_corners(points, 2., 1e-12).unwrap();
    let mut head = HeadMesh {
        builder: &mut builder,
        angles: (0..24).map(|i| TAU * f64::from(i) / 24.).collect(),
    };
    assert_eq!(
        head.drive_ring(0., &drive, 1e-12),
        Err(MeshError::UnresolvedBoundary)
    );
    assert!(head.builder.mesh.vertices.is_empty());
    head.angles = vec![0., TAU / 3., 2. * TAU / 3.];
    assert_eq!(
        head.drive_ring(0., &drive, 2.),
        Err(MeshError::UnresolvedBoundary)
    );
    assert!(head.builder.mesh.vertices.is_empty());
    let angles = points.map(|[x, y]| y.portable_atan2(x).rem_euclid(TAU));
    head.builder
        .insert_angular_stations(&angles, 1e-12)
        .unwrap();
    head.angles = (0..head.builder.steps)
        .map(|i| TAU * head.builder.fraction(i))
        .collect();
    let ring = head.drive_ring(0., &drive, 1e-12).unwrap();
    assert_eq!(ring.len(), 30);
    for p in points {
        assert!(
            head.builder
                .mesh
                .vertices
                .as_chunks::<3>()
                .0
                .contains(&[p[0], p[1], 0.])
        );
    }
}
