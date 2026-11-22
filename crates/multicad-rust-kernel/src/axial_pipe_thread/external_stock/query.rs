//! Actual original source coordinates and conservative whole-box material.
use super::super::authored::{StockClipClassification, StockClipPatch, StockClipSlab};
use super::*;
fn pair(v: I) -> [f64; 2] {
    [v.lo, v.hi]
}
#[derive(Clone, Copy, Debug)]
pub struct Point {
    position: [I; 3],
    radius: I,
    work: usize,
    bytes: usize,
}
impl Point {
    pub fn position_mm(&self) -> [[f64; 2]; 3] {
        self.position.map(pair)
    }
    pub fn radius_mm(&self) -> [f64; 2] {
        pair(self.radius)
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Query {
    classification: StockClipClassification,
    work: usize,
    bytes: usize,
}
impl Query {
    pub fn classification(&self) -> StockClipClassification {
        self.classification
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
}
impl ExternalStockCut {
    fn point(
        &self,
        vertex: Option<Vertex>,
        curve: Option<(Curve, [f64; 2])>,
        options: StockClipOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Point, Error> {
        let bytes = size_of::<Self>() + size_of::<Point>() + 16 * 1024;
        let mut work = Work::new(options, bytes, interrupt)?;
        work.charge(self.work)?;
        work.charge(128)?;
        let (original, project, top) = if let Some(vertex) = vertex {
            let (key, project, top) = match vertex {
                Vertex::Thread(v) => (v, false, false),
                Vertex::ShoulderRoot(rail) => {
                    if matches!(rail, StockClipRail::LowerCrest | StockClipRail::CrestUpper) {
                        return Err(Error::Invalid);
                    }
                    (
                        StockClipVertex {
                            station: StockClipStation::End,
                            rail,
                        },
                        true,
                        false,
                    )
                }
                Vertex::TailTop(rail) => (
                    StockClipVertex {
                        station: StockClipStation::End,
                        rail,
                    },
                    true,
                    true,
                ),
            };
            (
                self.thread
                    .vertex_bounds(key, work.child(self.thread.work_used())?, interrupt)?,
                project,
                top,
            )
        } else {
            let (key, t) = curve.ok_or(Error::Invalid)?;
            let (key, project, top) = match key {
                Curve::Thread(c) => {
                    if !self.thread.physical_curves().contains(&c) {
                        return Err(Error::Invalid);
                    }
                    (c, false, false)
                }
                Curve::ShoulderOuter(b) => {
                    if b == StockClipBand::Crest {
                        return Err(Error::Invalid);
                    }
                    (
                        StockClipCurve::Ring {
                            station: StockClipStation::End,
                            band: b,
                        },
                        true,
                        false,
                    )
                }
                Curve::TailTop(b) => (
                    StockClipCurve::Ring {
                        station: StockClipStation::End,
                        band: b,
                    },
                    true,
                    true,
                ),
            };
            (
                self.thread.curve_bounds(
                    key,
                    t,
                    work.child(self.thread.work_used())?,
                    interrupt,
                )?,
                project,
                top,
            )
        };
        work.charge(original.work_used() - self.thread.work_used())?;
        let mut position = original.position_mm().map(|v| I::new(v[0], v[1]));
        let mut radius = I::new(original.radius_mm()[0], original.radius_mm()[1]);
        if project {
            if radius.lo <= 0. {
                return Err(Error::Unresolved);
            }
            let factor = self.radius / radius;
            position[0] = position[0] * factor + self.center[0];
            position[1] = position[1] * factor + self.center[1];
            radius = self.radius;
        } else {
            position[0] = position[0] * self.radial + self.center[0];
            position[1] = position[1] * self.radial + self.center[1];
            radius = radius * self.radial;
        }
        position[2] = if top {
            self.top
        } else {
            position[2] + self.bottom
        };
        if position.iter().any(|v| !v.valid()) {
            return Err(Error::Unresolved);
        }
        work.charge(1)?;
        Ok(Point {
            position,
            radius,
            work: work.used,
            bytes,
        })
    }
    pub fn vertex_bounds(
        &self,
        key: Vertex,
        options: StockClipOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Point, Error> {
        self.point(Some(key), None, options, interrupt)
    }
    pub fn curve_bounds(
        &self,
        key: Curve,
        t: [f64; 2],
        options: StockClipOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Point, Error> {
        self.point(None, Some((key, t)), options, interrupt)
    }
    pub fn classify_box(
        &self,
        bounds: [[f64; 2]; 3],
        options: StockClipOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Query, Error> {
        if bounds
            .iter()
            .any(|v| !v[0].is_finite() || !v[1].is_finite() || v[0] > v[1])
        {
            return Err(Error::Invalid);
        }
        let bytes = size_of::<Self>() + size_of::<Query>() + 16 * 1024;
        let mut work = Work::new(options, bytes, interrupt)?;
        work.charge(self.work)?;
        work.charge(256)?;
        let [world_x, world_y, z] = bounds.map(|v| I::new(v[0], v[1]));
        let x = world_x - self.center[0];
        let y = world_y - self.center[1];
        use StockClipClassification::*;
        let classification = if z.hi < self.bottom.lo
            || z.lo > self.top.hi
            || x.lo > self.radius.hi
            || x.hi < -self.radius.hi
            || y.lo > self.radius.hi
            || y.hi < -self.radius.hi
        {
            Outside
        } else {
            let rho = (x.square() + y.square()).sqrt();
            if !rho.valid() {
                Undetermined
            } else if rho.lo > self.radius.hi {
                Outside
            } else if z.lo > self.join.hi {
                if z.hi < self.top.lo && rho.hi < self.radius.lo {
                    Inside
                } else {
                    Undetermined
                }
            } else {
                let root = self.thread.point_bounds(
                    StockClipPatch {
                        slab: StockClipSlab::Before,
                        band: StockClipBand::Root,
                    },
                    [0., 1.],
                    [0., 0.],
                    work.child(self.thread.work_used())?,
                    interrupt,
                )?;
                work.charge(root.work_used() - self.thread.work_used())?;
                if z.lo > self.bottom.hi
                    && z.hi < self.top.lo
                    && rho.hi < (I::point(root.radius_mm()[0]) * self.radial).lo
                {
                    Inside
                } else if z.hi < self.join.lo {
                    let local = [x / self.radial, y / self.radial, z - self.bottom];
                    let result = self.thread.classify_box(
                        local.map(pair),
                        work.child(self.thread.work_used())?,
                        interrupt,
                    )?;
                    work.charge(result.work_used() - self.thread.work_used())?;
                    result.classification()
                } else {
                    Undetermined
                }
            }
        };
        work.charge(1)?;
        Ok(Query {
            classification,
            work: work.used,
            bytes,
        })
    }
}
