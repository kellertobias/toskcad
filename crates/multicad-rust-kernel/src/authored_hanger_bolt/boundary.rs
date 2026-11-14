//! Finite source incidence, independent of the manufacturing partition.
//! Computational crest midpoints are not physical edges. Each annular shoulder
//! and the periodic plain wall retains two distinct oriented boundary loops.
use super::*;
use crate::authored_internal_treatment::ThreadBandRole;
mod curves;
mod incidence;
#[cfg(test)]
mod tests;
pub use curves::CurveBounds;
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Component {
    Metric,
    Wood,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum End {
    Start,
    End,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Knot {
    CrestEnd,
    UpperRoot,
    LowerRoot,
    CrestStart,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum VertexKey {
    Thread {
        component: Component,
        end: End,
        period: u32,
        knot: Knot,
    },
    Plain {
        end: End,
    },
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CurveKey {
    Rail {
        component: Component,
        period: u32,
        knot: Knot,
    },
    Cap {
        component: Component,
        end: End,
        period: u32,
        band: ThreadBandRole,
    },
    PlainCircle {
        end: End,
    },
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Coedge {
    pub curve: CurveKey,
    pub forward: bool,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Curve {
    pub key: CurveKey,
    pub vertices: [VertexKey; 2],
    pub faces: [FaceKey; 2],
}
#[derive(Clone, Copy)]
struct Face {
    key: FaceKey,
    spans: [(usize, usize); 2],
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BoundaryOptions {
    pub maximum_bytes: usize,
    pub maximum_work: usize,
}
impl Default for BoundaryOptions {
    fn default() -> Self {
        Self {
            maximum_bytes: 1 << 20,
            maximum_work: 16_000_000,
        }
    }
}
/// A sealed incidence receipt for these exact original authored components.
/// It does not certify a world placement, display mesh or exchange file.
pub struct Boundary {
    source: Hanger,
    vertices: Vec<VertexKey>,
    curves: Vec<Curve>,
    faces: Vec<Face>,
    coedges: Vec<Coedge>,
    bytes: usize,
}
// Maximum16originalstarts: four simultaneous rim proposals and each <=32
// coedges, plus fixed validation/iterator scratch, fit well inside64KiB. The
// caller's construction cap includes this transient reserve and owner headers.
const SCRATCH: usize = 64 << 10;
const KNOTS: [Knot; 4] = [
    Knot::CrestEnd,
    Knot::UpperRoot,
    Knot::LowerRoot,
    Knot::CrestStart,
];
const BANDS: [ThreadBandRole; 4] = [
    ThreadBandRole::Upper,
    ThreadBandRole::Root,
    ThreadBandRole::Lower,
    ThreadBandRole::Crest,
];
struct Work<'a> {
    left: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Work<'_> {
    fn charge(&mut self, n: usize) -> Result<(), Error> {
        checkpoint(self.interrupt)?;
        self.left = self.left.checked_sub(n).ok_or(Error::Limit)?;
        Ok(())
    }
}
impl Hanger {
    pub fn boundary(
        &self,
        options: BoundaryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Boundary, Error> {
        checkpoint(interrupt)?;
        if options.maximum_bytes == 0
            || options.maximum_bytes > 1 << 20
            || options.maximum_work == 0
            || options.maximum_work > 16_000_000
        {
            return Err(Error::Invalid);
        }
        let n = (self.metric().spec().starts + self.wood().spec().starts) as usize;
        let (nv, ne, nf) = (8 * n + 2, 12 * n + 2, self.face_count());
        let bytes = std::mem::size_of::<Boundary>()
            + nv * std::mem::size_of::<VertexKey>()
            + ne * std::mem::size_of::<Curve>()
            + nf * std::mem::size_of::<Face>()
            + 2 * ne * std::mem::size_of::<Coedge>();
        if bytes + SCRATCH > options.maximum_bytes {
            return Err(Error::Limit);
        }
        let mut w = Work {
            left: options.maximum_work,
            interrupt,
        };
        w.charge(nv + ne + nf + 2 * ne)?;
        let mut b = Boundary {
            source: *self,
            vertices: Vec::with_capacity(nv),
            curves: Vec::with_capacity(ne),
            faces: Vec::with_capacity(nf),
            coedges: Vec::with_capacity(2 * ne),
            bytes,
        };
        for component in [Component::Metric, Component::Wood] {
            let (starts, left) = b.settings(component);
            for period in 0..starts {
                for knot in KNOTS {
                    for end in [End::Start, End::End] {
                        b.vertices.push(VertexKey::Thread {
                            component,
                            end,
                            period,
                            knot,
                        });
                    }
                }
                for band in BANDS {
                    let (a, c) = b.endpoints(component, period, band);
                    let face = b.side_face(component, period, band);
                    for end in [End::Start, End::End] {
                        b.curves.push(Curve {
                            key: CurveKey::Cap {
                                component,
                                end,
                                period,
                                band,
                            },
                            vertices: [
                                VertexKey::Thread {
                                    component,
                                    end,
                                    period: a.0,
                                    knot: a.1,
                                },
                                VertexKey::Thread {
                                    component,
                                    end,
                                    period: c.0,
                                    knot: c.1,
                                },
                            ],
                            faces: [face, b.cap_face(component, end)],
                        });
                    }
                    let arc = |end| Coedge {
                        curve: CurveKey::Cap {
                            component,
                            end,
                            period,
                            band,
                        },
                        forward: true,
                    };
                    let rail = |(period, knot)| Coedge {
                        curve: CurveKey::Rail {
                            component,
                            period,
                            knot,
                        },
                        forward: true,
                    };
                    let mut loop_ = vec![
                        arc(End::Start),
                        rail(c),
                        Coedge {
                            forward: false,
                            ..arc(End::End)
                        },
                        Coedge {
                            forward: false,
                            ..rail(a)
                        },
                    ];
                    if !left {
                        reverse(&mut loop_);
                    }
                    b.add_face(face, &loop_, &[]);
                }
                for (i, knot) in KNOTS.into_iter().enumerate() {
                    let before = if i == 0 {
                        (period + starts - 1) % starts
                    } else {
                        period
                    };
                    let previous = BANDS[(i + 3) % 4];
                    let next = BANDS[i];
                    b.curves.push(Curve {
                        key: CurveKey::Rail {
                            component,
                            period,
                            knot,
                        },
                        vertices: [
                            VertexKey::Thread {
                                component,
                                end: End::Start,
                                period,
                                knot,
                            },
                            VertexKey::Thread {
                                component,
                                end: End::End,
                                period,
                                knot,
                            },
                        ],
                        faces: [
                            b.side_face(component, before, previous),
                            b.side_face(component, period, next),
                        ],
                    });
                }
            }
        }
        for end in [End::Start, End::End] {
            let vertex = VertexKey::Plain { end };
            b.vertices.push(vertex);
            b.curves.push(Curve {
                key: CurveKey::PlainCircle { end },
                vertices: [vertex; 2],
                faces: [
                    FaceKey::PlainWall,
                    if end == End::Start {
                        FaceKey::MetricShoulder
                    } else {
                        FaceKey::WoodShoulder
                    },
                ],
            });
        }
        let mut metric_bottom = b.rim(Component::Metric, End::Start);
        let mut metric_top = b.rim(Component::Metric, End::End);
        let mut wood_bottom = b.rim(Component::Wood, End::Start);
        let mut wood_top = b.rim(Component::Wood, End::End);
        // Increasing q runs clockwise for a right-hand source. The actual
        // cap normals are -z except for the wood tip; inner loops reverse.
        if b.settings(Component::Metric).1 {
            reverse(&mut metric_bottom);
        } else {
            reverse(&mut metric_top);
        }
        if b.settings(Component::Wood).1 {
            reverse(&mut wood_bottom);
        } else {
            reverse(&mut wood_top);
        }
        let plain_bottom = Coedge {
            curve: CurveKey::PlainCircle { end: End::Start },
            forward: true,
        };
        let plain_top = Coedge {
            curve: CurveKey::PlainCircle { end: End::End },
            forward: true,
        };
        b.add_face(FaceKey::MetricTip, &metric_bottom, &[]);
        b.add_face(
            FaceKey::MetricShoulder,
            &[Coedge {
                forward: false,
                ..plain_bottom
            }],
            &metric_top,
        );
        b.add_face(FaceKey::WoodShoulder, &wood_bottom, &[plain_top]);
        b.add_face(FaceKey::WoodTip, &wood_top, &[]);
        b.add_face(
            FaceKey::PlainWall,
            &[plain_bottom],
            &[Coedge {
                forward: false,
                ..plain_top
            }],
        );
        b.validate(nv, ne, nf, &mut w)?;
        b.check_memory(options.maximum_bytes)?;
        checkpoint(interrupt)?;
        Ok(b)
    }
}
fn reverse(loop_: &mut [Coedge]) {
    loop_.reverse();
    for e in loop_ {
        e.forward = !e.forward;
    }
}
impl Boundary {
    fn check_memory(&mut self, maximum: usize) -> Result<(), Error> {
        self.bytes = std::mem::size_of::<Boundary>()
            + self.vertices.capacity() * std::mem::size_of::<VertexKey>()
            + self.curves.capacity() * std::mem::size_of::<Curve>()
            + self.faces.capacity() * std::mem::size_of::<Face>()
            + self.coedges.capacity() * std::mem::size_of::<Coedge>();
        if self.bytes + SCRATCH > maximum {
            Err(Error::Limit)
        } else {
            Ok(())
        }
    }
    pub fn source(&self) -> &Hanger {
        &self.source
    }
    pub fn vertices(&self) -> &[VertexKey] {
        &self.vertices
    }
    pub fn curves(&self) -> &[Curve] {
        &self.curves
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
    pub fn face_loops(&self, key: FaceKey) -> Option<[&[Coedge]; 2]> {
        let f = self.faces.iter().find(|f| f.key == key)?;
        Some(f.spans.map(|(a, b)| &self.coedges[a..b]))
    }
    fn settings(&self, c: Component) -> (u32, bool) {
        match c {
            Component::Metric => {
                let s = self.source.metric().spec();
                (s.starts, s.left_hand)
            }
            Component::Wood => {
                let s = self.source.wood().spec();
                (s.starts, s.left_hand)
            }
        }
    }
    fn side_face(&self, c: Component, period: u32, band: ThreadBandRole) -> FaceKey {
        let (n, left) = self.settings(c);
        let carry = u32::from(matches!(
            band,
            ThreadBandRole::Lower | ThreadBandRole::Crest
        ));
        let group = (period + carry) % n;
        let start = if left { group } else { (n - group) % n };
        match c {
            Component::Metric => FaceKey::MetricThread { start, band },
            Component::Wood => FaceKey::WoodThread { start, band },
        }
    }
    fn cap_face(&self, c: Component, e: End) -> FaceKey {
        match (c, e) {
            (Component::Metric, End::Start) => FaceKey::MetricTip,
            (Component::Metric, End::End) => FaceKey::MetricShoulder,
            (Component::Wood, End::Start) => FaceKey::WoodShoulder,
            (Component::Wood, End::End) => FaceKey::WoodTip,
        }
    }
    fn endpoints(&self, c: Component, p: u32, band: ThreadBandRole) -> ((u32, Knot), (u32, Knot)) {
        match band {
            ThreadBandRole::Upper => ((p, Knot::CrestEnd), (p, Knot::UpperRoot)),
            ThreadBandRole::Root => ((p, Knot::UpperRoot), (p, Knot::LowerRoot)),
            ThreadBandRole::Lower => ((p, Knot::LowerRoot), (p, Knot::CrestStart)),
            ThreadBandRole::Crest => (
                (p, Knot::CrestStart),
                ((p + 1) % self.settings(c).0, Knot::CrestEnd),
            ),
        }
    }
    fn rim(&self, c: Component, e: End) -> Vec<Coedge> {
        (0..self.settings(c).0)
            .flat_map(|p| {
                BANDS.map(|band| Coedge {
                    curve: CurveKey::Cap {
                        component: c,
                        end: e,
                        period: p,
                        band,
                    },
                    forward: true,
                })
            })
            .collect()
    }
    fn add_face(&mut self, key: FaceKey, a: &[Coedge], b: &[Coedge]) {
        let start = self.coedges.len();
        self.coedges.extend_from_slice(a);
        let mid = self.coedges.len();
        self.coedges.extend_from_slice(b);
        self.faces.push(Face {
            key,
            spans: [(start, mid), (mid, self.coedges.len())],
        });
    }
}
