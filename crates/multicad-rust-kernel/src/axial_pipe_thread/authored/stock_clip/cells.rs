//! Source-owned angular/axial clipping cells for the eight retained charts.
//! Grid cuts are auxiliary; every original meridian/event remains present.
use super::*;
use std::mem::size_of;
mod mesh;
mod trace;
pub use mesh::{PlacedStockClipMeshCandidate, StockClipMeshCandidate, StockClipPlacement};
pub use trace::{StockClipCapTraceEdge, StockClipTracedMesh};
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct Tags {
    rail: Option<(usize, bool)>,
    station: Option<StockClipStation>,
    s: Option<i32>,
    z: Option<usize>,
}
#[derive(Clone, Copy, Debug)]
struct V {
    sz: [I; 2],
    tags: Tags,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Cut {
    S(i32),
    Z(usize),
}
/// Exact auxiliary station in the original source coordinates. Angular k
/// denotes s=q0+k*p/N; axial j denotes z=j*L/M. These are expressions,
/// not rounded chart coordinates or replacement authoring datums.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StockClipGridLine {
    Angular(i32),
    Axial(usize),
}
/// Immutable source ownership of a shared physical point. Original periodic
/// root aliases use the canonical RootLower curve and its matching station.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StockClipPointOwner {
    Vertex(StockClipVertex),
    Curve {
        curve: StockClipCurve,
        grid: StockClipGridLine,
    },
    Grid {
        patch: StockClipPatch,
        angular: i32,
        axial: usize,
    },
}
/// A source graph retaining original curve/event identities and per-cell
/// unwrapped coordinates. Work and byte limits bound this graph; literal
/// point/triangle limits apply separately at `mesh_candidate`.
#[derive(Debug)]
pub struct StockClipCellPlan {
    source: StockClipPlan,
    vertices: Vec<SourceVertex>,
    cells: Vec<SourceCell>,
    caps: Vec<CapEdge>,
    angular: usize,
    axial: usize,
    work: usize,
    bytes: usize,
}
#[derive(Clone, Copy, Debug)]
pub(super) struct SourceVertex {
    pub owner: StockClipPointOwner,
    pub radius: I,
    pub position: [I; 3],
}
#[derive(Clone, Copy, Debug)]
pub(super) struct SourceCell {
    pub column: i32,
    pub row: usize,
    pub ids: [u32; 8],
    pub parameters: [[I; 2]; 8],
    pub count: usize,
    pub patch: StockClipPatch,
}
#[derive(Clone, Copy, Debug)]
pub(super) struct CapEdge {
    pub ids: [u32; 2],
    pub parameters: [[I; 2]; 2],
    pub station: StockClipStation,
    pub patch: StockClipPatch,
}
struct Work<'a> {
    used: usize,
    limit: usize,
    bytes: usize,
    peak: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Work<'_> {
    fn charge(&mut self, n: usize) -> Result<(), Error> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        self.used = self
            .used
            .checked_add(n)
            .filter(|n| *n <= self.limit)
            .ok_or(Error::Limit)?;
        Ok(())
    }
    fn memory(&mut self, n: usize) -> Result<(), Error> {
        if n > self.bytes {
            return Err(Error::Limit);
        }
        self.peak = self.peak.max(n);
        Ok(())
    }
    fn observe(&mut self, n: usize) {
        self.peak = self.peak.max(n);
    }
    fn query(&self) -> Result<StockClipOptions, Error> {
        let remaining = (self.limit - self.used).min(16_000_000);
        if remaining == 0 {
            return Err(Error::Limit);
        }
        Ok(StockClipOptions {
            maximum_work: remaining,
            maximum_bytes: 64 * 1024 * 1024,
        })
    }
}
struct Builder<'a> {
    plan: &'a StockClipPlan,
    n: usize,
    m: usize,
    ds: I,
    dz: I,
    kmin: i32,
    kcount: usize,
    registry: Vec<u32>,
    vertices: Vec<SourceVertex>,
    cells: Vec<SourceCell>,
    caps: Vec<CapEdge>,
    work: Work<'a>,
}
impl StockClipCellPlan {
    pub fn source(&self) -> &StockClipPlan {
        &self.source
    }
    pub fn vertex_owner(&self, index: usize) -> Option<StockClipPointOwner> {
        self.vertices.get(index).map(|v| v.owner)
    }
    pub fn cell_grid(&self, index: usize) -> Option<(i32, usize)> {
        self.cells.get(index).map(|c| (c.column, c.row))
    }
    pub fn cell_patch(&self, index: usize) -> Option<StockClipPatch> {
        self.cells.get(index).map(|c| c.patch)
    }
    pub fn cell_vertices(&self, index: usize) -> Option<&[u32]> {
        self.cells.get(index).map(|c| &c.ids[..c.count])
    }
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }
    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }
    pub fn cap_edge_count(&self) -> usize {
        self.caps.len()
    }
    pub fn subdivisions(&self) -> [usize; 2] {
        [self.angular, self.axial]
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
}
fn b(p: [f64; 2]) -> I {
    I::new(p[0], p[1])
}
// The same retained q0/q3/contacts source law, expressed as s=a-beta*z.
// k=1/32, sqrt3>1.7 and 0<slopes<.6 prove .98<beta<1.02.
// Intersect with that earned positive interval before inversion, rather than
// dividing nearly equal endpoint boxes on a very short source slab.
fn rail_law(plan: &StockClipPlan, slab: StockClipSlab, rail: StockClipRail, shift: bool) -> (I, I) {
    let one = I::point(1.);
    let k = I::point(plan.source.spec.radial_taper);
    let delta = I::point(plan.stock_radius_limit_mm()) - I::point(plan.source.spec.major_radius_mm);
    let (a, raw) = match rail {
        StockClipRail::RootLower => (
            plan.q0 + I::point(if shift { plan.source.spec.pitch_mm } else { 0. }),
            one,
        ),
        StockClipRail::UpperRoot => (plan.q3, one),
        StockClipRail::LowerCrest if slab == StockClipSlab::After => (
            -plan.half + delta * plan.slopes[0],
            one + k * plan.slopes[0],
        ),
        StockClipRail::CrestUpper if slab == StockClipSlab::After => {
            (plan.half - delta * plan.slopes[1], one - k * plan.slopes[1])
        }
        StockClipRail::LowerCrest => (-plan.half, one),
        StockClipRail::CrestUpper => (plan.half, one),
    };
    (a, I::new(raw.lo.max(0.98), raw.hi.min(1.02)))
}
// Lowest occupied dyadic exponent of a finite, nonzero binary64 number.
fn lowest_bit(value: f64) -> i32 {
    let bits = value.to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i32;
    let fraction = bits & ((1u64 << 52) - 1);
    let (significand, quantum) = if exponent == 0 {
        (fraction, -1074)
    } else {
        (fraction | (1u64 << 52), exponent - 1075)
    };
    quantum + significand.trailing_zeros() as i32
}
fn exact_product(a: f64, b: f64) -> Option<[f64; 2]> {
    if !a.is_finite() || !b.is_finite() {
        return None;
    }
    let value = a * b;
    if !value.is_finite() {
        return None;
    }
    if a == 0. || b == 0. {
        return Some([value, 0.]);
    }
    // FMA's residual is exact only if its full dyadic tail is representable.
    // In the sealed source calls p/N and L/M (optionally /32) have lowest
    // bits >=-78; bounded integer indices cannot lower that exponent.
    // The explicit applicability check also protects this numeric primitive
    // from tiny future inputs, instead of equating rounded zero products.
    if lowest_bit(a) + lowest_bit(b) < -1074 {
        return None;
    }
    Some([value, a.mul_add(b, -value)])
}
fn same_products(a: [f64; 2], b: [f64; 2]) -> bool {
    exact_product(a[0], a[1])
        .zip(exact_product(b[0], b[1]))
        .is_some_and(|(a, b)| a == b)
}
fn exact_difference(a: f64, b: f64) -> [f64; 2] {
    let value = a - b;
    let bv = a - value;
    let av = value + bv;
    [value, (a - av) + (bv - b)]
}
fn merge(a: Tags, c: Tags) -> Tags {
    Tags {
        rail: (a.rail == c.rail).then_some(a.rail).flatten(),
        station: (a.station == c.station).then_some(a.station).flatten(),
        s: (a.s == c.s).then_some(a.s).flatten(),
        z: (a.z == c.z).then_some(a.z).flatten(),
    }
}
impl Builder<'_> {
    fn retained(&self) -> usize {
        size_of::<Self>()
            + size_of::<StockClipCellPlan>()
            + self.registry.capacity() * 4
            + self.vertices.capacity() * size_of::<SourceVertex>()
            + self.cells.capacity() * size_of::<SourceCell>()
            + self.caps.capacity() * size_of::<CapEdge>()
            + SCRATCH
    }
    fn reserve(&mut self, which: usize, n: usize) -> Result<(), Error> {
        let (len, cap, unit) = match which {
            0 => (
                self.vertices.len(),
                self.vertices.capacity(),
                size_of::<SourceVertex>(),
            ),
            1 => (
                self.cells.len(),
                self.cells.capacity(),
                size_of::<SourceCell>(),
            ),
            _ => (self.caps.len(), self.caps.capacity(), size_of::<CapEdge>()),
        };
        if len + n <= cap {
            return Ok(());
        }
        let target = (cap * 2).max(len + n).max(16);
        self.work.memory(self.retained() + target * unit)?;
        self.work.charge(target * unit)?;
        match which {
            0 => self.vertices.try_reserve_exact(target - len),
            1 => self.cells.try_reserve_exact(target - len),
            _ => self.caps.try_reserve_exact(target - len),
        }
        .map_err(|_| Error::Limit)?;
        // Pinned rustc58cf RawVec::grow_exact records len+additional.
        // The preflight includes BOTH old and requested new non-ZST buffers;
        // observing all actual capacities cannot exceed that earned bound.
        // Re-audit on toolchain change (not a public Vec API guarantee):
        // https://github.com/rust-lang/rust/blob/58cf8f9c35b355199727627d07ae34a4e2846d67/library/alloc/src/raw_vec/mod.rs#L734-L740
        self.work.observe(self.retained());
        Ok(())
    }
    fn coord(&self, cut: Cut) -> I {
        match cut {
            Cut::S(k) => self.plan.q0 + I::point(k as f64) * self.ds,
            Cut::Z(j) => I::point(j as f64) * self.dz,
        }
    }
    // Exact binary-source identities only. Product witnesses explicitly
    // qualify residual representability; lack of a witness never equals a
    // second failure. The bounded source difference cannot overflow.
    fn same_height(&self, tags: Tags, cut: Cut) -> bool {
        if let (Some(StockClipStation::Contact), Cut::Z(j)) = (tags.station, cut) {
            return exact_product(
                self.plan.source.spec.radial_taper
                    * (self.plan.source.spec.length_mm / self.m as f64),
                j as f64,
            )
            .is_some_and(|p| {
                p == exact_difference(
                    self.plan.stock_radius_limit_mm(),
                    self.plan.source.spec.major_radius_mm,
                )
            });
        }
        let Some((0, shift)) = tags.rail else {
            return false;
        };
        let angular = match cut {
            Cut::S(k) => Some(k),
            Cut::Z(_) => tags.s,
        };
        let Some(k) = angular else { return false };
        let coefficient = (if shift { self.n as i32 } else { 0 }) - k;
        let step = self.plan.source.spec.pitch_mm / self.n as f64;
        match (cut, tags.z, tags.station) {
            (Cut::Z(j), _, _) | (_, Some(j), _) => same_products(
                [step, coefficient as f64],
                [self.plan.source.spec.length_mm / self.m as f64, j as f64],
            ),
            // Original Start/End stations always retain their exact z=0/M
            // grid tags; the preceding branch already handles both.
            (_, _, Some(StockClipStation::Contact)) => exact_product(
                step * self.plan.source.spec.radial_taper,
                coefficient as f64,
            )
            .is_some_and(|p| {
                p == exact_difference(
                    self.plan.stock_radius_limit_mm(),
                    self.plan.source.spec.major_radius_mm,
                )
            }),
            _ => false,
        }
    }
    fn side(&self, v: V, cut: Cut, lower: bool) -> Result<i8, Error> {
        let (p, known) = match cut {
            Cut::S(k) => (v.sz[0], v.tags.s == Some(k)),
            Cut::Z(j) => (v.sz[1], v.tags.z == Some(j)),
        };
        if known || self.same_height(v.tags, cut) {
            return Ok(0);
        }
        let bound = self.coord(cut);
        let sign = if p.hi < bound.lo {
            -1
        } else if p.lo > bound.hi {
            1
        } else {
            return Err(Error::Unresolved);
        };
        Ok(if lower { sign } else { -sign })
    }
    fn intersection(&self, a: V, c: V, cut: Cut, patch: StockClipPatch) -> Result<V, Error> {
        let mut tags = merge(a.tags, c.tags);
        let sz = if let Some((rail, shift)) = tags.rail {
            let (at_zero, rate) = rail_law(self.plan, patch.slab, RAILS[rail], shift);
            match cut {
                Cut::S(_) => {
                    let s = self.coord(cut);
                    [s, (at_zero - s) / rate]
                }
                Cut::Z(_) => {
                    let z = self.coord(cut);
                    [at_zero - rate * z, z]
                }
            }
        } else if let Some(station) = tags.station {
            // A fixed-height ring cannot have opposite strict Z signs.
            // Only an angular cut reaches this admitted crossing.
            [self.coord(cut), self.plan.height(station)]
        } else {
            // Z clipping precedes S clipping. A non-source crossing is
            // therefore an angular cut through an existing axial grid edge;
            // a subsequent parallel angular cut cannot cross its own edge.
            let z = tags.z.expect("earned axial edge before angular clipping");
            [self.coord(cut), self.coord(Cut::Z(z))]
        };
        match cut {
            Cut::S(k) => tags.s = Some(k),
            Cut::Z(j) => tags.z = Some(j),
        };
        Ok(V { sz, tags })
    }
    fn clip(
        &mut self,
        input: &[V],
        cut: Cut,
        lower: bool,
        patch: StockClipPatch,
    ) -> Result<Vec<V>, Error> {
        // At most eight classifications, their bounded dyadic witnesses,
        // and two rail/grid intersections. Charge the complete stage before
        // evaluating it; the whole request's 512M ceiling is unchanged.
        self.work.charge(512)?;
        let mut out = Vec::with_capacity(8);
        if input.is_empty() {
            return Ok(out);
        }
        let mut a = input[input.len() - 1];
        let mut sa = self.side(a, cut, lower)?;
        for &point in input {
            let sc = self.side(point, cut, lower)?;
            let mut c = point;
            if sc == 0 {
                match cut {
                    Cut::S(k) => c.tags.s = Some(k),
                    Cut::Z(j) => c.tags.z = Some(j),
                }
            }
            if sa * sc < 0 {
                out.push(self.intersection(a, c, cut, patch)?)
            }
            if sc >= 0 {
                out.push(c)
            }
            a = c;
            sa = sc;
        }
        Ok(out)
    }
    fn address(&self, v: V, patch: StockClipPatch) -> (usize, StockClipPointOwner) {
        let station = |s| match s {
            StockClipStation::Start => 0,
            StockClipStation::Contact => 1,
            StockClipStation::End => 2,
        };
        if let (Some(s), Some((r, _))) = (v.tags.station, v.tags.rail) {
            return (
                station(s) * 4 + r,
                StockClipPointOwner::Vertex(StockClipVertex {
                    station: s,
                    rail: RAILS[r],
                }),
            );
        }
        let mut base = 12;
        if let (Some((r, shift)), Some(k)) = (v.tags.rail, v.tags.s) {
            let k = k - if shift { self.n as i32 } else { 0 };
            let index = k - self.kmin;
            // Actual source q lies within[q0,q0+p] and 0<=z<=L.
            // Thus canonical k lies in[-N*L/p,N]. kmin includes two
            // extra whole turns and the upper allocation includes another;
            // clipping never publishes the conservative search-only margins.
            return (
                base + ((patch.slab == StockClipSlab::After) as usize * 4 + r) * self.kcount
                    + index as usize,
                StockClipPointOwner::Curve {
                    curve: StockClipCurve::Axial {
                        slab: patch.slab,
                        rail: RAILS[r],
                    },
                    grid: StockClipGridLine::Angular(k),
                },
            );
        }
        base += 8 * self.kcount;
        if let (Some((r, _)), Some(j)) = (v.tags.rail, v.tags.z) {
            return (
                base + ((patch.slab == StockClipSlab::After) as usize * 4 + r) * (self.m + 1) + j,
                StockClipPointOwner::Curve {
                    curve: StockClipCurve::Axial {
                        slab: patch.slab,
                        rail: RAILS[r],
                    },
                    grid: StockClipGridLine::Axial(j),
                },
            );
        }
        base += 8 * (self.m + 1);
        if let (Some(s), Some(k)) = (v.tags.station, v.tags.s) {
            return (
                base + (station(s) * 4 + band_index(patch.band)) * self.n
                    + k.rem_euclid(self.n as i32) as usize,
                StockClipPointOwner::Curve {
                    curve: StockClipCurve::Ring {
                        station: s,
                        band: patch.band,
                    },
                    grid: StockClipGridLine::Angular(k),
                },
            );
        }
        base += 12 * self.n;
        // All vertices originate as source corners or intersections of
        // a source edge/grid line or two grid lines. The preceding cases
        // exhaust corners and source edges; both remaining tags are earned.
        let k = v.tags.s.expect("angular grid crossing");
        let j = v.tags.z.expect("axial grid crossing");
        (
            base + j * self.n + k.rem_euclid(self.n as i32) as usize,
            StockClipPointOwner::Grid {
                patch,
                angular: k,
                axial: j,
            },
        )
    }
    fn vertex(&mut self, v: V, patch: StockClipPatch) -> Result<u32, Error> {
        let (address, owner) = self.address(v, patch);
        let old = self.registry[address];
        if old != u32::MAX {
            return Ok(old);
        }
        // This is an owned SOURCE graph, constrained by the work/byte
        // ledger before every growth. The actual literal vertex limit belongs
        // to mesh_candidate, including its two later cap centers. The 64MiB
        // non-ZST source budget also keeps every published ID far below u32.
        let [s, z] = v.sz;
        let q = s + z;
        let k = I::point(self.plan.source.spec.radial_taper);
        let root = self.plan.root_radius + k * z;
        let raw = match patch.band {
            StockClipBand::Lower => root + (I::point(3.).sqrt() - k) * (q - self.plan.q0),
            StockClipBand::Upper => root + (-I::point(3.).sqrt() - k) * (q - self.plan.q3),
            StockClipBand::Crest if patch.slab == StockClipSlab::After => {
                I::point(self.plan.stock_radius_limit_mm())
            }
            StockClipBand::Crest => I::point(self.plan.source.spec.major_radius_mm) + k * z,
            StockClipBand::Root => root,
        };
        let radius = I::new(
            raw.lo.max(self.plan.radius_range.lo),
            raw.hi.min(self.plan.radius_range.hi),
        );
        let physical_s = if let Some(index) = v.tags.s {
            self.plan.q0 + I::point(index.rem_euclid(self.n as i32) as f64) * self.ds
        } else {
            // RootUpper physical aliases were already registered by the
            // earlier RootLower chart, so no new lifted-root point reaches
            // allocation here. The per-cell unwrapped source lift is kept.
            s
        };
        let angle = -physical_s * (I::point(2.) * pi() / I::point(self.plan.source.spec.pitch_mm));
        let vertex = SourceVertex {
            owner,
            radius,
            position: [radius * angle.cos(), radius * angle.sin(), z],
        };
        self.work.charge(256)?;
        self.reserve(0, 1)?;
        let id = self.vertices.len() as u32;
        self.vertices.push(vertex);
        self.registry[address] = id;
        Ok(id)
    }
    fn patch(&mut self, patch: StockClipPatch) -> Result<(), Error> {
        let mut original = [V {
            sz: [I::point(0.); 2],
            tags: Tags::default(),
        }; 4];
        let (start, end) = if patch.slab == StockClipSlab::Before {
            (StockClipStation::Start, StockClipStation::Contact)
        } else {
            (StockClipStation::Contact, StockClipStation::End)
        };
        let index = band_index(patch.band);
        for (i, (u, v)) in [(0., 0.), (1., 0.), (1., 1.), (0., 1.)]
            .into_iter()
            .enumerate()
        {
            let p = self.plan.point_bounds(
                patch,
                [u, u],
                [v, v],
                self.work.query()?,
                self.work.interrupt,
            )?;
            self.work.charge(p.work_used())?;
            let station = if v == 0. { start } else { end };
            let r = if u == 0. { index } else { (index + 1) % 4 };
            let shift = u == 1. && patch.band == StockClipBand::Root;
            let z = b(p.axial_mm());
            let mut tags = Tags {
                rail: Some((r, shift)),
                station: Some(station),
                ..Tags::default()
            };
            tags.z = match station {
                StockClipStation::Start => Some(0),
                StockClipStation::End => Some(self.m),
                _ => None,
            };
            if station == StockClipStation::Start && r == 0 {
                tags.s = Some(if shift { self.n as i32 } else { 0 })
            }
            original[i] = V {
                sz: [b(p.phase_mm()) - z, z],
                tags,
            };
        }
        for row in 0..self.m {
            let lower = self.clip(&original, Cut::Z(row), true, patch)?;
            let horizontal = self.clip(&lower, Cut::Z(row + 1), false, patch)?;
            if horizontal.len() < 3 {
                continue;
            }
            let lo = horizontal
                .iter()
                .map(|p| ((p.sz[0] - self.plan.q0) / self.ds).lo)
                .fold(f64::INFINITY, f64::min)
                .floor() as i32
                - 1;
            let hi = horizontal
                .iter()
                .map(|p| ((p.sz[0] - self.plan.q0) / self.ds).hi)
                .fold(f64::NEG_INFINITY, f64::max)
                .ceil() as i32
                + 1;
            for col in lo..hi {
                let lower = self.clip(&horizontal, Cut::S(col), true, patch)?;
                let polygon = self.clip(&lower, Cut::S(col + 1), false, patch)?;
                if polygon.len() < 3 {
                    continue;
                }
                // A convex quad clipped by four halfplanes has <=8
                // vertices. After the earlier <3 test, this polygon has
                // positive area: original axial rails have strictly negative
                // ds/dz and the two angular cuts are distinct. Horizontal
                // zero-area slices were removed before angular clipping.
                let mut ids = [0; 8];
                let mut parameters = [[I::point(0.); 2]; 8];
                for (i, &vertex) in polygon.iter().enumerate() {
                    ids[i] = self.vertex(vertex, patch)?;
                    parameters[i] = vertex.sz;
                }
                self.reserve(1, 1)?;
                self.cells.push(SourceCell {
                    column: col,
                    row,
                    ids,
                    parameters,
                    count: polygon.len(),
                    patch,
                });
                for i in 0..polygon.len() {
                    let j = (i + 1) % polygon.len();
                    let station = merge(polygon[i].tags, polygon[j].tags).station;
                    if let Some(station @ (StockClipStation::Start | StockClipStation::End)) =
                        station
                    {
                        self.reserve(2, 1)?;
                        self.caps.push(CapEdge {
                            ids: [ids[i], ids[j]],
                            parameters: [parameters[i], parameters[j]],
                            station,
                            patch,
                        });
                    }
                }
            }
        }
        Ok(())
    }
}
impl StockClipPlan {
    pub fn cell_plan(
        &self,
        angular: usize,
        axial: usize,
        maximum_work: usize,
        maximum_bytes: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<StockClipCellPlan, Error> {
        if !(4..=512).contains(&angular)
            || !angular.is_power_of_two()
            || !(4..=2048).contains(&axial)
            || !axial.is_power_of_two()
            || maximum_work == 0
            || maximum_work > 512_000_000
            || maximum_bytes == 0
            || maximum_bytes > 64 * 1024 * 1024
        {
            return Err(Error::Invalid);
        }
        let n = angular;
        let m = axial;
        let turns = (self.source.spec.length_mm / self.source.spec.pitch_mm).ceil() as i32 + 2;
        let kmin = -(n as i32) * turns;
        let kcount = (n as i32 * 2 - kmin + 1) as usize;
        let count = 12 + 8 * kcount + 8 * (m + 1) + 12 * n + n * (m + 1);
        let mut builder = Builder {
            plan: self,
            n,
            m,
            ds: I::point(self.source.spec.pitch_mm) / I::point(n as f64),
            dz: I::point(self.source.spec.length_mm) / I::point(m as f64),
            kmin,
            kcount,
            registry: Vec::new(),
            vertices: Vec::new(),
            cells: Vec::new(),
            caps: Vec::new(),
            work: Work {
                used: 0,
                limit: maximum_work,
                bytes: maximum_bytes,
                peak: 0,
                interrupt,
            },
        };
        builder.work.memory(
            SCRATCH + size_of::<Builder<'_>>() + size_of::<StockClipCellPlan>() + 4 * count,
        )?;
        builder.work.charge(self.work + count)?;
        builder
            .registry
            .try_reserve_exact(count)
            .map_err(|_| Error::Limit)?;
        builder.registry.resize(count, u32::MAX);
        // Initial non-ZST reserve_exact has the requested logical capacity at
        // the pinned RawVec commit audited above. Keep the actual observation.
        builder.work.observe(builder.retained());
        for patch in self.patches() {
            builder.patch(patch)?
        }
        builder.work.charge(1)?;
        Ok(StockClipCellPlan {
            source: *self,
            vertices: builder.vertices,
            cells: builder.cells,
            caps: builder.caps,
            angular: n,
            axial: m,
            work: builder.work.used,
            bytes: builder.work.peak,
        })
    }
}
#[cfg(test)]
mod tests;
