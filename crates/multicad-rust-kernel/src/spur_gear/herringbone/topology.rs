use super::*;
use multicad_precision::PortableMath as _;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Band {
    Lower,
    Upper,
}

/// A transverse tooth-outline station. Continuous bands share the apex;
/// a gapped chevron has a lower and an upper shoulder instead.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Ring {
    Start,
    Apex,
    LowerShoulder,
    UpperShoulder,
    End,
}

impl Herringbone {
    pub fn ring_height_mm(&self, ring: Ring) -> f64 {
        let height = self.profile.height_mm();
        match ring {
            Ring::Start => 0.,
            Ring::Apex => height * 0.5,
            Ring::LowerShoulder => self.band_length_mm(),
            Ring::UpperShoulder => height - self.band_length_mm(),
            Ring::End => height,
        }
    }

    /// The tooth-outline stations in axial order.
    pub fn rings(&self) -> &'static [Ring] {
        if self.is_gapped() {
            &[
                Ring::Start,
                Ring::LowerShoulder,
                Ring::UpperShoulder,
                Ring::End,
            ]
        } else {
            &[Ring::Start, Ring::Apex, Ring::End]
        }
    }

    /// Unrotated angles of the root-circle arc under one tooth, from its
    /// right fillet foot clockwise to its left fillet foot.
    pub fn under_tooth_angles(&self, tooth: u32) -> [f64; 2] {
        let outline = self.profile.outline();
        let per_tooth = outline.len() / self.teeth() as usize;
        let first = &outline[tooth as usize * per_tooth];
        let last = &outline[tooth as usize * per_tooth + per_tooth - 1];
        let [right, left] =
            [last.edge.point(0.), first.edge.point(0.)].map(|[x, y]| y.portable_atan2(x));
        // The left foot precedes the right foot by less than one pitch.
        [right, right - (right - left).rem_euclid(TAU)]
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Face {
    Cap {
        end: bool,
    },
    Tooth {
        band: Band,
        edge: ToothEdge,
    },
    Bore,
    HubWall {
        extension: HubExtension,
    },
    HubCap {
        extension: HubExtension,
    },
    /// The root cylinder spanning the central gap.
    RootCore,
    /// The planar end of one tooth at a gap shoulder.
    Shoulder {
        upper: bool,
        tooth: u32,
    },
}

impl Face {
    pub fn role(self) -> String {
        match self {
            Self::Cap { end: false } => "gear/end.start".into(),
            Self::Cap { end: true } => "gear/end.end".into(),
            Self::Bore => "gear/bore".into(),
            Self::RootCore => "gear/root.core".into(),
            Self::Shoulder { upper: false, .. } => "gear/gap.lower".into(),
            Self::Shoulder { upper: true, .. } => "gear/gap.upper".into(),
            Self::HubWall { .. } => "gear/hub.side".into(),
            Self::HubCap { extension } => if extension.end_z_mm < extension.gear_z_mm {
                "gear/hub.end.start"
            } else {
                "gear/hub.end.end"
            }
            .into(),
            Self::Tooth { band, edge } => {
                let band = match band {
                    Band::Lower => "lower",
                    Band::Upper => "upper",
                };
                let role = match edge.role {
                    super::super::ToothRole::LeftFillet => "fillet.left",
                    super::super::ToothRole::LeftRadial => "continuation.left",
                    super::super::ToothRole::LeftFlank => "involute.left",
                    super::super::ToothRole::Tip => "tip",
                    super::super::ToothRole::RightFlank => "involute.right",
                    super::super::ToothRole::RightRadial => "continuation.right",
                    super::super::ToothRole::RightFillet => "fillet.right",
                    super::super::ToothRole::Root => "root",
                };
                format!("gear/{band}/tooth/{}/{role}", edge.tooth)
            }
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Curve {
    Rim {
        ring: Ring,
        edge: Edge,
    },
    Rail {
        band: Band,
        point: Point,
    },
    BoreRim {
        end: bool,
        radius_mm: f64,
    },
    BoreSeam {
        radius_mm: f64,
    },
    HubRim {
        radius_mm: f64,
        z_mm: f64,
    },
    HubSeam {
        radius_mm: f64,
        span_mm: [f64; 2],
    },
    /// A root-circle arc under one tooth at a gap shoulder, running from
    /// `angles[0]` to `angles[1]` before the shoulder phase rotation.
    UnderTooth {
        ring: Ring,
        angles: [f64; 2],
        radius_mm: f64,
    },
    /// The straight seam of the root core through one outline vertex.
    RootSeam {
        point: Point,
    },
}
impl Curve {
    pub fn point(self, gear: &Herringbone, fraction: f64) -> [f64; 3] {
        match self {
            Self::Rim { ring, edge } => {
                gear.section_point(edge.point(fraction), gear.ring_height_mm(ring))
            }
            Self::Rail { band, point } => {
                let [start, end] = gear.band_interval_mm(band);
                gear.section_point(point, start + (end - start) * fraction)
            }
            Self::UnderTooth {
                ring,
                angles: [from, to],
                radius_mm,
            } => {
                let angle = (to - from).mul_add(fraction, from);
                gear.section_point(
                    [
                        radius_mm * angle.portable_cos(),
                        radius_mm * angle.portable_sin(),
                    ],
                    gear.ring_height_mm(ring),
                )
            }
            Self::RootSeam { point } => {
                let [lo, hi] = gear.shoulder_heights_mm().expect("gapped root seam");
                gear.section_point(point, (hi - lo).mul_add(fraction, lo))
            }
            Self::BoreRim { end, radius_mm } => {
                let [x, y] = rotate([radius_mm, 0.], -TAU * fraction);
                [x, y, gear.axial_bounds()[usize::from(end)]]
            }
            Self::BoreSeam { radius_mm } => {
                let [lo, hi] = gear.axial_bounds();
                [radius_mm, 0., (hi - lo).mul_add(fraction, lo)]
            }
            Self::HubRim { radius_mm, z_mm } => {
                let [x, y] = rotate([radius_mm, 0.], TAU * fraction);
                [x, y, z_mm]
            }
            Self::HubSeam {
                radius_mm,
                span_mm: [lo, hi],
            } => [radius_mm, 0., (hi - lo).mul_add(fraction, lo)],
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TopologyEdge {
    pub vertices: [u32; 2],
    /// Face using the directed edge, then face using its reversal. A
    /// periodic cylindrical seam is used twice by the bore face.
    pub faces: [u32; 2],
    pub curve: Curve,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Topology {
    pub vertices_mm: Vec<[f64; 3]>,
    pub edges: Vec<TopologyEdge>,
    pub faces: Vec<Face>,
}

impl Herringbone {
    pub fn topology(&self) -> Topology {
        let outline = self.profile.outline();
        let count = outline.len() as u32;
        let extensions = self.hub_extensions();
        let mut result = Topology {
            vertices_mm: Vec::with_capacity(4 * outline.len() + 6),
            edges: Vec::with_capacity(5 * outline.len() + 3),
            faces: vec![Face::Cap { end: false }, Face::Cap { end: true }],
        };
        let wall = |band: u32, index: u32| 2 + band * count + index;
        let rings = self.rings();
        let gap = self.gap_faces_start();
        let teeth = self.teeth();
        let per_tooth = count / teeth.max(1);
        let is_root = |index: u32| outline[index as usize].role == super::super::ToothRole::Root;
        // Above a lower shoulder edge lies its tooth's shoulder or the core.
        let shoulder = |upper: bool, index: u32| {
            if is_root(index) {
                gap
            } else {
                gap + 1 + u32::from(upper) * teeth + outline[index as usize].tooth
            }
        };
        for (station, &ring) in rings.iter().enumerate() {
            let base = station as u32 * count;
            for (index, item) in (0..count).zip(outline) {
                let curve = Curve::Rim {
                    ring,
                    edge: item.edge,
                };
                let faces = match ring {
                    Ring::Start => [wall(0, index), 0],
                    Ring::Apex => [wall(1, index), wall(0, index)],
                    Ring::LowerShoulder => [shoulder(false, index), wall(0, index)],
                    Ring::UpperShoulder => [wall(1, index), shoulder(true, index)],
                    Ring::End => [1, wall(1, index)],
                };
                result.vertices_mm.push(curve.point(self, 0.));
                result.edges.push(TopologyEdge {
                    vertices: [base + index, base + (index + 1) % count],
                    faces,
                    curve,
                });
            }
        }
        let last = rings.len() as u32 - 1;
        for (layer, band) in [Band::Lower, Band::Upper].into_iter().enumerate() {
            let layer = layer as u32;
            // Gapped bands are stations 0→1 and 2→3.
            let from = if layer == 0 { 0 } else { last - 1 };
            for (index, item) in (0..count).zip(outline) {
                result.faces.push(Face::Tooth { band, edge: *item });
                result.edges.push(TopologyEdge {
                    vertices: [from * count + index, (from + 1) * count + index],
                    faces: [wall(layer, (index + count - 1) % count), wall(layer, index)],
                    curve: Curve::Rail {
                        band,
                        point: item.edge.point(0.),
                    },
                });
            }
        }
        let gap_faces = if self.is_gapped() {
            let mut faces = vec![Face::RootCore];
            for upper in [false, true] {
                faces.extend((0..teeth).map(|tooth| Face::Shoulder { upper, tooth }));
            }
            faces
        } else {
            Vec::new()
        };
        if self.is_gapped() {
            let radius_mm = self.profile.root_radius_mm;
            for (station, ring, upper) in [
                (1, Ring::LowerShoulder, false),
                (2, Ring::UpperShoulder, true),
            ] {
                let base = station * count;
                for tooth in 0..teeth {
                    // From the right fillet foot back to the left foot.
                    let first = tooth * per_tooth;
                    let right = first + per_tooth - 1;
                    let face = gap + 1 + u32::from(upper) * teeth + tooth;
                    result.edges.push(TopologyEdge {
                        vertices: [base + right, base + first],
                        faces: if upper { [gap, face] } else { [face, gap] },
                        curve: Curve::UnderTooth {
                            ring,
                            angles: self.under_tooth_angles(tooth),
                            radius_mm,
                        },
                    });
                }
            }
            // The seam runs through the first root edge's start vertex.
            let seam = per_tooth - 1;
            result.edges.push(TopologyEdge {
                vertices: [count + seam, 2 * count + seam],
                faces: [gap; 2],
                curve: Curve::RootSeam {
                    point: outline[seam as usize].edge.point(0.),
                },
            });
        }
        if let Some(radius_mm) = self.profile.bore_radius_mm() {
            let face = result.faces.len() as u32;
            let start = result.vertices_mm.len() as u32;
            result.faces.push(Face::Bore);
            let [lo, hi] = self.axial_bounds();
            result
                .vertices_mm
                .extend([[radius_mm, 0., lo], [radius_mm, 0., hi]]);
            for (index, end) in [false, true].into_iter().enumerate() {
                let vertex = start + index as u32;
                let cap = extensions[index].map_or(index as u32, |extension| extension.cap_face);
                result.edges.push(TopologyEdge {
                    vertices: [vertex, vertex],
                    faces: if end { [cap, face] } else { [face, cap] },
                    curve: Curve::BoreRim { end, radius_mm },
                });
            }
            result.edges.push(TopologyEdge {
                vertices: [start, start + 1],
                faces: [face; 2],
                curve: Curve::BoreSeam { radius_mm },
            });
        }
        result.faces.extend(gap_faces);
        for extension in extensions.into_iter().flatten() {
            result
                .faces
                .extend([Face::HubWall { extension }, Face::HubCap { extension }]);
            let upper = extension.end_z_mm > extension.gear_z_mm;
            let low = extension.gear_z_mm.min(extension.end_z_mm);
            let high = extension.gear_z_mm.max(extension.end_z_mm);
            let start = result.vertices_mm.len() as u32;
            result.vertices_mm.extend([
                [extension.radius_mm, 0., low],
                [extension.radius_mm, 0., high],
            ]);
            let cap = u32::from(upper);
            for (index, z_mm) in [low, high].into_iter().enumerate() {
                let vertex = start + index as u32;
                let adjacent = if z_mm == extension.gear_z_mm {
                    cap
                } else {
                    extension.cap_face
                };
                result.edges.push(TopologyEdge {
                    vertices: [vertex; 2],
                    faces: if index == 0 {
                        [extension.wall_face, adjacent]
                    } else {
                        [adjacent, extension.wall_face]
                    },
                    curve: Curve::HubRim {
                        radius_mm: extension.radius_mm,
                        z_mm,
                    },
                });
            }
            result.edges.push(TopologyEdge {
                vertices: [start, start + 1],
                faces: [extension.wall_face; 2],
                curve: Curve::HubSeam {
                    radius_mm: extension.radius_mm,
                    span_mm: [low, high],
                },
            });
        }
        result
    }
}
