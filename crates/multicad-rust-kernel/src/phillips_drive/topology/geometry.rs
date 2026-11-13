use super::*;
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Surface {
    Plane {
        normal_bounds: [Bounds; 3],
        offset_mm_bounds: Bounds,
    },
    Cylinder {
        origin_mm_bounds: [Bounds; 3],
        axis_bounds: [Bounds; 3],
        radius_mm: f64,
        outward_sign: i8,
    },
    Torus {
        major_radius_mm_bounds: Bounds,
        minor_radius_mm: f64,
        center_depth_mm: f64,
        outward_sign: i8,
    },
    Cone {
        radius_at_zero_mm_bounds: Bounds,
        radius_decrease_per_depth_bounds: Bounds,
    },
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LineSupport {
    pub origin_mm_bounds: [Bounds; 3],
    pub direction_bounds: [Bounds; 3],
}
fn pi() -> I {
    I::new(PI.next_down(), PI.next_up())
}
impl BoundaryTopology<'_> {
    pub fn face_surface(&self, index: usize) -> Result<Surface, Error> {
        let face = self.faces.get(index).ok_or(Error::Invalid)?;
        let c = self.events.source();
        let k = &c.corners[0];
        let e = &c.envelope;
        match face.key {
            FaceKey::Mouth => Ok(Surface::Plane {
                normal_bounds: [
                    bounds(I::point(0.)),
                    bounds(I::point(0.)),
                    bounds(I::point(-1.)),
                ],
                offset_mm_bounds: bounds(I::point(0.)),
            }),
            FaceKey::Envelope { role, .. } => match role {
                EnvelopeRole::EntryRoundover => Ok(Surface::Torus {
                    major_radius_mm_bounds: bounds(e.entry),
                    minor_radius_mm: c.spec.mouth_radius_mm,
                    center_depth_mm: c.spec.mouth_radius_mm,
                    outward_sign: -1,
                }),
                EnvelopeRole::OuterFlank => Ok(Surface::Cone {
                    radius_at_zero_mm_bounds: bounds(
                        I::point(c.spec.mouth_diameter_mm) / I::point(2.),
                    ),
                    radius_decrease_per_depth_bounds: bounds(e.sin_a / e.cos_a),
                }),
                EnvelopeRole::Floor => Ok(Surface::Cone {
                    radius_at_zero_mm_bounds: bounds(e.tip * e.cos_c / e.sin_c),
                    radius_decrease_per_depth_bounds: bounds(e.cos_c / e.sin_c),
                }),
            },
            FaceKey::Wall {
                kind,
                quarter,
                lower,
            } => match (kind, &k.inner) {
                (WallKind::Inner, SourceInner::Rounded { radius_mm }) => {
                    let base = k.b_half + k.b_depth * k.sin_beta / k.cos_beta;
                    let [x, y] = rotate(
                        [base + I::point(*radius_mm) / k.cos_beta, I::point(0.)],
                        quarter,
                    );
                    let [u, v] = rotate([-k.sin_beta, I::point(0.)], quarter);
                    Ok(Surface::Cylinder {
                        origin_mm_bounds: [bounds(x), bounds(y), bounds(I::point(0.))],
                        axis_bounds: [bounds(u), bounds(v), bounds(k.cos_beta)],
                        radius_mm: *radius_mm,
                        outward_sign: -1,
                    })
                }
                (WallKind::Inner, SourceInner::ShallowV { sine, cosine }) => {
                    Ok(wall_plane(k, *sine, *cosine, I::point(0.), quarter, lower))
                }
                (WallKind::Outer, _) => Ok(wall_plane(
                    k,
                    k.sin_theta,
                    k.cos_theta,
                    k.join_s * k.sin_theta - k.join_y * k.cos_theta,
                    quarter,
                    lower,
                )),
            },
        }
    }
    /// Actual authored direction of straight source edges. Consumers must not
    /// replace this with a normalized difference of rounded endpoint values.
    pub fn curve_line_support(&self, index: usize) -> Result<Option<LineSupport>, Error> {
        let edge = *self.edges.get(index).ok_or(Error::Invalid)?;
        let c = self.events.source();
        let k = &c.corners[0];
        let direction = match edge.key {
            CurveKey::Mouth { edge: index } => match self.events.mouth().edges()[index].curve() {
                SectionCurve::Line { normal_bounds, .. } => {
                    let [x, y] = normal_bounds.map(i);
                    let length = (x.square() + y.square()).sqrt();
                    [-y / length, x / length, I::point(0.)]
                }
                _ => return Ok(None),
            },
            CurveKey::WallJoin { quarter, .. } | CurveKey::Center { quarter } => {
                let [x, y] = rotate([-k.sin_beta, I::point(0.)], quarter);
                [x, y, k.cos_beta]
            }
            CurveKey::FloorGenerator { quarter } => {
                let [x, y] = rotate([-c.envelope.cos_c, I::point(0.)], quarter);
                [x, y, c.envelope.sin_c]
            }
            _ => return Ok(None),
        };
        Ok(Some(LineSupport {
            origin_mm_bounds: self.vertices[edge.ends[0]].point_mm_bounds(),
            direction_bounds: direction.map(bounds),
        }))
    }
    /// At the real cone apex the normal is a limiting angular family, not a
    /// single invented plane normal. t spans this quarter's actual cone sector.
    pub fn floor_apex_normal_limit(&self, quarter: u8, t: f64) -> Result<[Bounds; 3], Error> {
        if quarter > 3 || !t.is_finite() || !(0. ..=1.).contains(&t) {
            return Err(Error::Invalid);
        }
        let angle = (I::point(quarter as f64) + I::point(t)) * pi() / I::point(2.);
        let e = &self.events.source().envelope;
        Ok([
            bounds(e.sin_c * angle.cos()),
            bounds(e.sin_c * angle.sin()),
            bounds(e.cos_c),
        ])
    }
    pub fn curve_point_mm_bounds(&self, index: usize, t: f64) -> Result<[Bounds; 3], Error> {
        if index >= self.edges.len() || !t.is_finite() || !(0. ..=1.).contains(&t) {
            return Err(Error::Invalid);
        }
        let edge = self.edges[index];
        if t == 0. || t == 1. {
            return Ok(self.vertices[edge.ends[t as usize]].point_mm_bounds());
        }
        let u = I::point(t);
        let point = match edge.key {
            CurveKey::Trim {
                envelope,
                wall,
                quarter,
                lower,
            } => {
                let index = self
                    .events
                    .trims()
                    .iter()
                    .position(|t| t.envelope_role() == envelope && t.wall_kind() == wall)
                    .ok_or(Error::Unresolved)?;
                let [x, y, z] = self.events.trim_point_mm_bounds(index, t)?.map(i);
                let [x, y] = rotate([x, if lower { -y } else { y }], quarter);
                [x, y, z]
            }
            CurveKey::Mouth { edge: index } => {
                let curve = self.events.mouth().edges()[index].curve();
                match curve {
                    SectionCurve::CircleArc {
                        radius_mm_bounds,
                        start_angle_rad_bounds,
                        sweep_rad_bounds,
                    } => {
                        let angle = i(start_angle_rad_bounds) + u * i(sweep_rad_bounds);
                        let r = i(radius_mm_bounds);
                        [r * angle.cos(), r * angle.sin(), I::point(0.)]
                    }
                    SectionCurve::EllipseArc {
                        center_mm_bounds,
                        semiaxes_mm_bounds,
                        quarter_turn,
                        start_angle_rad_bounds,
                        sweep_rad_bounds,
                    } => {
                        let angle = i(start_angle_rad_bounds) + u * i(sweep_rad_bounds);
                        let [x, y] = rotate(
                            [
                                i(semiaxes_mm_bounds[0]) * angle.cos(),
                                i(semiaxes_mm_bounds[1]) * angle.sin(),
                            ],
                            quarter_turn,
                        );
                        [
                            i(center_mm_bounds[0]) + x,
                            i(center_mm_bounds[1]) + y,
                            I::point(0.),
                        ]
                    }
                    SectionCurve::Line { .. } => self.line_point(edge, u),
                }
            }
            CurveKey::MeridianCircle { station, quarter } => {
                let p = self.events.event(station.event()).point_mm_bounds().map(i);
                let phi = (p[1] / p[0]).atan();
                let angle = phi
                    + I::point(quarter as f64) * pi() / I::point(2.)
                    + u * (pi() / I::point(2.) - I::point(2.) * phi);
                let c = self.events.source();
                let radius = match station {
                    MeridianStation::Entry => {
                        c.envelope.entry - I::point(c.spec.mouth_radius_mm) * c.envelope.cos_a
                    }
                    MeridianStation::Floor => I::point(c.spec.g_diameter_mm) / I::point(2.),
                };
                [radius * angle.cos(), radius * angle.sin(), p[2]]
            }
            CurveKey::WallJoin { .. }
            | CurveKey::Center { .. }
            | CurveKey::FloorGenerator { .. } => self.line_point(edge, u),
        };
        Ok(point.map(bounds))
    }
    fn line_point(&self, edge: BoundaryEdge, u: I) -> [I; 3] {
        let [a, b] = edge.ends.map(|i| self.vertices[i].point);
        std::array::from_fn(|i| a[i] + u * (b[i] - a[i]))
    }
}

fn wall_plane(k: &Corner, sine: I, cosine: I, extra: I, quarter: u8, lower: bool) -> Surface {
    let [x, y] = rotate(
        [k.cos_beta * sine, if lower { cosine } else { -cosine }],
        quarter,
    );
    let z = k.sin_beta * sine;
    let offset = sine * (k.cos_beta * k.b_half + k.sin_beta * k.b_depth) + extra;
    Surface::Plane {
        normal_bounds: [bounds(x), bounds(y), bounds(z)],
        offset_mm_bounds: bounds(offset),
    }
}
