//! Unions of coaxial right prisms over exact planar profiles: boxes,
//! cylinders, annuli, slots, bored stock, spur gears with hubs and internal
//! rings. Every answer is evaluated on the profile's lines, circular arcs and
//! involutes; nothing is sampled.
use super::paths::{Lift, ProfilePath, Segment, add, dot2, edge_point, scaled};
use super::profile::{Family, Rings, contains, exclusive, intersect, merge};
use super::*;
use multicad_rust_kernel::spur_gear::Edge;
use std::sync::OnceLock;

/// One axial band `axis[0]..axis[1]` of a profile.
pub(in crate::drawing) struct Band {
    pub rings: Rings,
    pub axis: [f64; 2],
}

struct Layer {
    rings: Rings,
    axis: [f64; 2],
    /// Monotone along the paper direction perpendicular to the sight line.
    sight: Family,
    /// Monotone along the profile-plane part of the sight line.
    cut: OnceLock<Family>,
    /// Monotone in the second profile coordinate, for containment.
    inside: Family,
}

pub(in crate::drawing) struct Prism {
    origin: P3,
    axes: [P3; 3],
    /// The direction toward the viewer in profile coordinates.
    toward: P3,
    /// Length of its profile-plane part, and that part normalized.
    slope: f64,
    along: P2,
    across: P2,
    layers: Vec<Layer>,
    tolerance: f64,
}

const AXIAL: f64 = 1e-9;

impl Prism {
    /// `axes` are the orthonormal profile x, profile y and extrusion
    /// directions in view space; bands are stacked without overlap.
    pub(in crate::drawing) fn new(
        origin: P3,
        axes: [P3; 3],
        mut bands: Vec<Band>,
        tolerance: f64,
    ) -> Result<Self, GeometryError> {
        bands.sort_by(|a, b| a.axis[0].total_cmp(&b.axis[0]));
        if bands.is_empty()
            || bands.iter().any(|band| {
                band.axis[1] - band.axis[0] <= tolerance || band.rings.iter().any(Vec::is_empty)
            })
            || bands
                .windows(2)
                .any(|pair| pair[1].axis[0] < pair[0].axis[1] - tolerance)
        {
            return Err(err(
                GeometryErrorKind::InvalidResult,
                "Invalid prism bands for a drawing scene",
            ));
        }
        let toward = axes.map(|axis| axis[2]);
        let slope = toward[0].hypot(toward[1]);
        let along = if slope > AXIAL {
            [toward[0] / slope, toward[1] / slope]
        } else {
            [1., 0.]
        };
        let across = [-along[1], along[0]];
        let layers = bands
            .into_iter()
            .map(|band| Layer {
                sight: Family::new(&band.rings, across),
                inside: Family::new(&band.rings, [0., 1.]),
                cut: OnceLock::new(),
                rings: band.rings,
                axis: band.axis,
            })
            .collect();
        Ok(Self {
            origin,
            axes,
            toward,
            slope,
            along,
            across,
            layers,
            tolerance,
        })
    }

    fn axial(&self) -> bool {
        self.slope <= AXIAL
    }

    fn edge_on(&self) -> bool {
        self.toward[2].abs() <= AXIAL
    }

    fn place(&self, [x, y]: P2, w: f64) -> P3 {
        std::array::from_fn(|i| {
            self.origin[i] + x * self.axes[0][i] + y * self.axes[1][i] + w * self.axes[2][i]
        })
    }

    fn level(&self, w: f64) -> Lift {
        Lift {
            origin: add(self.origin, scaled(self.axes[2], w)),
            x: self.axes[0],
            y: self.axes[1],
        }
    }

    fn cut<'a>(&self, layer: &'a Layer) -> &'a Family {
        layer
            .cut
            .get_or_init(|| Family::new(&layer.rings, self.along))
    }

    /// Profile coordinates of the sight line through a paper point at z = 0.
    fn local(&self, [x, y]: P2) -> P3 {
        let offset = [x - self.origin[0], y - self.origin[1], -self.origin[2]];
        self.axes.map(|axis| paths::dot3(axis, offset))
    }

    /// Material depths of one layer's unbounded profile along a sight line.
    fn profile_spans(&self, layer: &Layer, base: P3) -> Result<Vec<[f64; 2]>, GeometryError> {
        let start = [base[0], base[1]];
        layer
            .sight
            .spans(&layer.rings, dot2(start, self.across), |point| {
                (dot2(point, self.along) - dot2(start, self.along)) / self.slope
            })
    }

    /// Rings of a layer that continue unchanged into a neighboring layer.
    fn shared(&self, layer: usize, ring: usize, end: usize) -> bool {
        let level = self.layers[layer].axis[end];
        self.layers.iter().enumerate().any(|(other, neighbor)| {
            other != layer
                && (neighbor.axis[1 - end] - level).abs() <= self.tolerance
                && neighbor.rings.contains(&self.layers[layer].rings[ring])
        })
    }

    /// Material spans of the layers ending (`side` 1) or starting at a level.
    fn level_spans(
        &self,
        level: f64,
        side: usize,
        spans: impl Fn(&Layer) -> Result<Vec<[f64; 2]>, GeometryError>,
    ) -> Result<Vec<[f64; 2]>, GeometryError> {
        let mut output = Vec::new();
        for layer in &self.layers {
            if (layer.axis[side] - level).abs() <= self.tolerance {
                output.extend(spans(layer)?);
            }
        }
        Ok(merge(output, self.tolerance))
    }

    fn levels(&self) -> Vec<f64> {
        let mut levels: Vec<f64> = self.layers.iter().flat_map(|layer| layer.axis).collect();
        levels.sort_by(f64::total_cmp);
        levels.dedup_by(|a, b| (*a - *b).abs() <= self.tolerance);
        levels
    }
}

fn path(edge: Edge, range: [f64; 2], lift: Lift, straight: bool) -> Arc<dyn Path> {
    Arc::new(ProfilePath {
        edge,
        range,
        lift,
        straight,
    })
}

fn turns(before: Edge, after: Edge) -> bool {
    let [a, b] = [before.tangent(1.), after.tangent(0.)];
    (a[0] * b[1] - a[1] * b[0]).abs() > 1e-9 || dot2(a, b) < 0.
}

impl Member for Prism {
    fn bounds(&self) -> [P3; 2] {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for layer in &self.layers {
            for k in 0..3 {
                let direction = [self.axes[0][k], self.axes[1][k]];
                let reach = |sign: f64| {
                    layer
                        .rings
                        .iter()
                        .flatten()
                        .map(|edge| edge.support(direction.map(|v| v * sign)))
                        .fold(f64::NEG_INFINITY, f64::max)
                };
                let axial = layer.axis.map(|w| self.origin[k] + w * self.axes[2][k]);
                lo[k] = lo[k].min(axial[0].min(axial[1]) - reach(-1.));
                hi[k] = hi[k].max(axial[0].max(axial[1]) + reach(1.));
            }
        }
        [lo, hi]
    }

    fn candidates(&self, _: f64) -> Result<Vec<Candidate>, GeometryError> {
        let mut output = Vec::new();
        let mut push = |kind, path| output.push(Candidate { kind, path });
        for (index, layer) in self.layers.iter().enumerate() {
            for (ring, edges) in layer.rings.iter().enumerate() {
                for end in 0..2 {
                    if self.shared(index, ring, end) {
                        continue;
                    }
                    let lift = self.level(layer.axis[end]);
                    for &edge in edges {
                        // Split where the rim turns back in the paper, which
                        // is also where the wall's outline generator ends.
                        let mut stops = vec![0.];
                        if self.axial() {
                            if matches!(edge, Edge::Planar(e) if e.start() == e.end()) {
                                stops.push(0.5);
                            }
                        } else {
                            stops.extend(profile::folds(edge, self.across));
                        }
                        stops.push(1.);
                        for pair in stops.windows(2) {
                            push(
                                Kind::Edge,
                                path(edge, [pair[0], pair[1]], lift, self.edge_on()),
                            );
                        }
                    }
                }
                if self.axial() {
                    continue;
                }
                let mut generator = |kind, point: P2| {
                    push(
                        kind,
                        Arc::new(Segment(
                            self.place(point, layer.axis[0]),
                            self.place(point, layer.axis[1]),
                        )),
                    );
                };
                for (position, &edge) in edges.iter().enumerate() {
                    for fraction in profile::folds(edge, self.across) {
                        generator(Kind::Outline, edge_point(edge, fraction));
                    }
                    let next = edges[(position + 1) % edges.len()];
                    let closed = edges.len() == 1 && edge_point(edge, 0.) == edge_point(edge, 1.);
                    if closed && dot2(edge.tangent(0.), self.across).abs() <= 1e-12 {
                        // The outline falls on the seam of a closed edge.
                        generator(Kind::Outline, edge_point(edge, 0.));
                    }
                    if !closed {
                        generator(
                            if turns(edge, next) {
                                Kind::Edge
                            } else {
                                Kind::Smooth
                            },
                            edge_point(edge, 1.),
                        );
                    }
                }
            }
        }
        Ok(output)
    }

    fn sight(&self, paper: P2) -> Result<Vec<[f64; 2]>, GeometryError> {
        let base = self.local(paper);
        let rise = self.toward[2];
        let mut output = Vec::new();
        if self.axial() {
            for layer in &self.layers {
                if contains(&layer.inside, &layer.rings, [base[0], base[1]]) {
                    let ends = layer.axis.map(|w| (w - base[2]) / rise);
                    output.push([ends[0].min(ends[1]), ends[0].max(ends[1])]);
                }
            }
            return Ok(merge(output, self.tolerance));
        }
        if self.edge_on() {
            // The sight line stays at one axial level. On a level shared by
            // two layers it is inside the solid only where both have material.
            for level in self.levels() {
                if (base[2] - level).abs() <= self.tolerance {
                    let spans = |layer: &Layer| self.profile_spans(layer, base);
                    return Ok(intersect(
                        &self.level_spans(level, 1, spans)?,
                        &self.level_spans(level, 0, spans)?,
                    ));
                }
            }
            for layer in &self.layers {
                if layer.axis[0] < base[2] && base[2] < layer.axis[1] {
                    output.extend(self.profile_spans(layer, base)?);
                }
            }
            return Ok(merge(output, self.tolerance));
        }
        for layer in &self.layers {
            let ends = layer.axis.map(|w| (w - base[2]) / rise);
            let slab = [[ends[0].min(ends[1]), ends[0].max(ends[1])]];
            output.extend(intersect(&self.profile_spans(layer, base)?, &slab));
        }
        Ok(merge(output, self.tolerance))
    }

    fn line(&self, start: P3, direction: P3) -> Option<Result<Vec<[f64; 2]>, GeometryError>> {
        let base = self
            .axes
            .map(|axis| paths::dot3(axis, paths::add(start, paths::scaled(self.origin, -1.))));
        let rate = self.axes.map(|axis| paths::dot3(axis, direction));
        let planar = rate[0].hypot(rate[1]);
        let length = paths::dot3(direction, direction).sqrt();
        Some((|| {
            let mut output = Vec::new();
            for layer in &self.layers {
                // The axial slab of the layer along the line.
                let slab = if rate[2].abs() <= AXIAL * length {
                    if base[2] <= layer.axis[0] || base[2] >= layer.axis[1] {
                        continue;
                    }
                    [f64::NEG_INFINITY, f64::INFINITY]
                } else {
                    let ends = layer.axis.map(|w| (w - base[2]) / rate[2]);
                    [ends[0].min(ends[1]), ends[0].max(ends[1])]
                };
                let profile = if planar <= AXIAL * length {
                    if contains(&layer.inside, &layer.rings, [base[0], base[1]]) {
                        vec![[f64::NEG_INFINITY, f64::INFINITY]]
                    } else {
                        continue;
                    }
                } else {
                    let along = [rate[0] / planar, rate[1] / planar];
                    let family = Family::new(&layer.rings, [-along[1], along[0]]);
                    let start = [base[0], base[1]];
                    family.spans(&layer.rings, dot2(start, [-along[1], along[0]]), |point| {
                        (dot2(point, along) - dot2(start, along)) / planar
                    })?
                };
                output.extend(intersect(&profile, &[slab]));
            }
            Ok(merge(output, 0.))
        })())
    }

    fn section(&self, depth: f64, _: f64) -> Result<Vec<Loop>, GeometryError> {
        let offset = depth - self.origin[2];
        let rise = self.toward[2];
        let mut output: Loop = Vec::new();
        if self.axial() {
            // Material immediately behind the plane decides at a shared level.
            // On a cap with material only in front, the plane still meets the
            // solid in that cap.
            let level = offset / rise;
            let step = rise.signum() * 2. * self.tolerance;
            let within = |w: f64| {
                self.layers
                    .iter()
                    .position(|layer| layer.axis[0] < w && w < layer.axis[1])
            };
            let chosen = within(level - step).or_else(|| within(level + step));
            if let Some(layer) = chosen.map(|index| &self.layers[index]) {
                let lift = self.level(level);
                for edges in &layer.rings {
                    let mut ring: Loop = Vec::new();
                    for &edge in edges {
                        let closed = matches!(edge, Edge::Planar(e) if e.start() == e.end());
                        let stops: &[f64] = if closed { &[0., 0.5, 1.] } else { &[0., 1.] };
                        for pair in stops.windows(2) {
                            let range = if rise < 0. {
                                [pair[1], pair[0]]
                            } else {
                                [pair[0], pair[1]]
                            };
                            ring.push(path(edge, range, lift, false));
                        }
                    }
                    if rise < 0. {
                        ring.reverse();
                    }
                    output.extend(ring);
                }
            }
            return Ok(if output.is_empty() {
                vec![]
            } else {
                vec![output]
            });
        }
        let levels = self.levels();
        let segment = |a: P3, b: P3| -> Arc<dyn Path> { Arc::new(Segment(a, b)) };
        if self.edge_on() {
            // The plane contains the axis direction: rectangles between the
            // crossings of one profile line, joined across shared levels.
            let line = offset / self.slope;
            let at = |s: f64| [0, 1].map(|i| line * self.along[i] + s * self.across[i]);
            let spans = |layer: &Layer| {
                self.cut(layer)
                    .spans(&layer.rings, line, |point| dot2(point, self.across))
            };
            for layer in &self.layers {
                for stop in spans(layer)?.into_iter().flatten() {
                    output.push(segment(
                        self.place(at(stop), layer.axis[0]),
                        self.place(at(stop), layer.axis[1]),
                    ));
                }
            }
            for level in levels {
                let below = self.level_spans(level, 1, spans)?;
                let above = self.level_spans(level, 0, spans)?;
                for [a, b] in exclusive(&below, &above, self.tolerance) {
                    output.push(segment(self.place(at(a), level), self.place(at(b), level)));
                }
            }
            return Ok(if output.is_empty() {
                vec![]
            } else {
                vec![output]
            });
        }
        // Oblique plane: the affine image of each profile between the lines
        // where the plane meets the layer's two levels.
        let lift = Lift {
            origin: add(self.origin, scaled(self.axes[2], offset / rise)),
            x: add(self.axes[0], scaled(self.axes[2], -self.toward[0] / rise)),
            y: add(self.axes[1], scaled(self.axes[2], -self.toward[1] / rise)),
        };
        let line = |w: f64| (offset - rise * w) / self.slope;
        for layer in &self.layers {
            let family = self.cut(layer);
            let ends = layer.axis.map(line);
            let [lo, hi] = [ends[0].min(ends[1]), ends[0].max(ends[1])];
            for piece in &family.pieces {
                let [a, b] = piece.level;
                if a.max(b) <= lo || a.min(b) >= hi {
                    continue;
                }
                let stop = |index: usize| {
                    let value = piece.level[index];
                    if value < lo {
                        family.solve(&layer.rings, piece, lo).fraction
                    } else if value > hi {
                        family.solve(&layer.rings, piece, hi).fraction
                    } else {
                        piece.range[index]
                    }
                };
                let range = [stop(0), stop(1)];
                if (range[1] - range[0]).abs() > 1e-14 {
                    output.push(path(
                        layer.rings[piece.ring][piece.edge],
                        range,
                        lift,
                        false,
                    ));
                }
            }
        }
        for level in levels {
            let value = line(level);
            let at = |s: f64| [0, 1].map(|i| value * self.along[i] + s * self.across[i]);
            let spans = |layer: &Layer| {
                self.cut(layer)
                    .spans(&layer.rings, value, |point| dot2(point, self.across))
            };
            let below = self.level_spans(level, 1, spans)?;
            let above = self.level_spans(level, 0, spans)?;
            for [a, b] in exclusive(&below, &above, self.tolerance) {
                output.push(segment(self.place(at(a), level), self.place(at(b), level)));
            }
        }
        Ok(if output.is_empty() {
            vec![]
        } else {
            vec![output]
        })
    }
}
