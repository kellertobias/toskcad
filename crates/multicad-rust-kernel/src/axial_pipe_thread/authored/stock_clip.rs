//! Authored NPT stock-clipping atlas, source cell plans and bounded literal
//! mesh candidates. Native/STEP/material admission and whole literal embedding
//! remain separate; no independently decoded file mass is supplied here.
use super::*;
use std::mem::size_of;
mod cells;
mod curve_law;
pub use cells::{
    PlacedStockClipMeshCandidate, StockClipCellPlan, StockClipGridLine, StockClipMeshCandidate,
    StockClipPlacement, StockClipPointOwner,
};
pub use curve_law::{StockClipAffineCurveTaylor, StockClipCurveLaw, StockClipCurveTaylor};
mod jets;
mod measures;
mod polyline;
mod query;
pub use measures::{StockClipFaceMeasure, StockClipMeasureOptions};
pub use polyline::StockClipPolyline;
pub use query::{StockClipClassification, StockClipQuery};
mod sample;
#[cfg(test)]
mod tests;
mod topology;
mod volume;
pub use jets::{StockClipJets, StockClipPatchLaw, StockClipPatchTaylor};
pub use topology::StockClipFaceLoop;

pub use sample::{StockClipPoint, StockClipSample};
pub use volume::StockClipVolume;
const SCRATCH: usize = 16 * 1024;
#[derive(Clone, Copy, Debug)]
pub struct StockClipOptions {
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for StockClipOptions {
    fn default() -> Self {
        Self {
            maximum_work: 16_000_000,
            maximum_bytes: 64 * 1024 * 1024,
        }
    }
}
struct Meter<'a> {
    used: usize,
    remaining: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl<'a> Meter<'a> {
    fn new(
        options: StockClipOptions,
        retained: usize,
        interrupt: Option<&'a dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        if options.maximum_work == 0
            || options.maximum_work > 16_000_000
            || options.maximum_bytes == 0
            || options.maximum_bytes > 64 * 1024 * 1024
        {
            return Err(Error::Invalid);
        }
        if retained.saturating_add(SCRATCH) > options.maximum_bytes {
            return Err(Error::Limit);
        }
        Ok(Self {
            used: 0,
            remaining: options.maximum_work,
            interrupt,
        })
    }
    fn child_options(&self, options: StockClipOptions) -> Result<StockClipOptions, Error> {
        if self.remaining == 0 {
            return Err(Error::Limit);
        }
        Ok(StockClipOptions {
            maximum_work: self.remaining,
            ..options
        })
    }
    fn charge(&mut self, n: usize) -> Result<(), Error> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        self.remaining = self.remaining.checked_sub(n).ok_or(Error::Limit)?;
        self.used += n;
        Ok(())
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StockClipSlab {
    Before,
    After,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StockClipBand {
    Lower,
    Crest,
    Upper,
    Root,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StockClipStation {
    Start,
    Contact,
    End,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StockClipRail {
    RootLower,
    LowerCrest,
    CrestUpper,
    UpperRoot,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StockClipFace {
    LowerFlank,
    NaturalCrest,
    StockCylinder,
    UpperFlank,
    Root,
    StartCap,
    EndCap,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct StockClipPatch {
    pub slab: StockClipSlab,
    pub band: StockClipBand,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct StockClipVertex {
    pub station: StockClipStation,
    pub rail: StockClipRail,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StockClipCurve {
    Axial {
        slab: StockClipSlab,
        rail: StockClipRail,
    },
    Ring {
        station: StockClipStation,
        band: StockClipBand,
    },
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct StockClipCoedge {
    pub curve: StockClipCurve,
    pub forward: bool,
}
const BANDS: [StockClipBand; 4] = [
    StockClipBand::Lower,
    StockClipBand::Crest,
    StockClipBand::Upper,
    StockClipBand::Root,
];
const RAILS: [StockClipRail; 4] = [
    StockClipRail::RootLower,
    StockClipRail::LowerCrest,
    StockClipRail::CrestUpper,
    StockClipRail::UpperRoot,
];
fn band_index(b: StockClipBand) -> usize {
    match b {
        StockClipBand::Lower => 0,
        StockClipBand::Crest => 1,
        StockClipBand::Upper => 2,
        StockClipBand::Root => 3,
    }
}
fn rail_index(r: StockClipRail) -> usize {
    match r {
        StockClipRail::RootLower => 0,
        StockClipRail::LowerCrest => 1,
        StockClipRail::CrestUpper => 2,
        StockClipRail::UpperRoot => 3,
    }
}
impl StockClipPatch {
    pub fn face(self) -> StockClipFace {
        match self.band {
            StockClipBand::Lower => StockClipFace::LowerFlank,
            StockClipBand::Crest => {
                if self.slab == StockClipSlab::Before {
                    StockClipFace::NaturalCrest
                } else {
                    StockClipFace::StockCylinder
                }
            }
            StockClipBand::Upper => StockClipFace::UpperFlank,
            StockClipBand::Root => StockClipFace::Root,
        }
    }
    /// Outward oriented chart boundary. Contact rings of lower/upper/root
    /// are auxiliary seams of the SAME physical face, not extra material faces.
    pub fn coedges(self) -> [StockClipCoedge; 4] {
        let i = band_index(self.band);
        let (bottom, top) = match self.slab {
            StockClipSlab::Before => (StockClipStation::Start, StockClipStation::Contact),
            StockClipSlab::After => (StockClipStation::Contact, StockClipStation::End),
        };
        [
            StockClipCoedge {
                curve: StockClipCurve::Axial {
                    slab: self.slab,
                    rail: RAILS[i],
                },
                forward: true,
            },
            StockClipCoedge {
                curve: StockClipCurve::Ring {
                    station: top,
                    band: self.band,
                },
                forward: true,
            },
            StockClipCoedge {
                curve: StockClipCurve::Axial {
                    slab: self.slab,
                    rail: RAILS[(i + 1) % 4],
                },
                forward: false,
            },
            StockClipCoedge {
                curve: StockClipCurve::Ring {
                    station: bottom,
                    band: self.band,
                },
                forward: false,
            },
        ]
    }
}
impl StockClipCurve {
    pub fn endpoints(self) -> [StockClipVertex; 2] {
        match self {
            Self::Axial { slab, rail } => {
                let (a, b) = match slab {
                    StockClipSlab::Before => (StockClipStation::Start, StockClipStation::Contact),
                    StockClipSlab::After => (StockClipStation::Contact, StockClipStation::End),
                };
                [
                    StockClipVertex { station: a, rail },
                    StockClipVertex { station: b, rail },
                ]
            }
            Self::Ring { station, band } => {
                let i = band_index(band);
                [
                    StockClipVertex {
                        station,
                        rail: RAILS[i],
                    },
                    StockClipVertex {
                        station,
                        rail: RAILS[(i + 1) % 4],
                    },
                ]
            }
        }
    }
}
/// The contact event is the retained expression (stock_radius-major_radius)/k.
/// Bounds enclose that expression; no midpoint becomes an authoring datum.
#[derive(Clone, Copy, Debug)]
pub struct StockClipPlan {
    source: AuthoredThread,
    event: I,
    root_radius: I,
    q0: I,
    q3: I,
    half: I,
    width: I,
    slopes: [I; 2],
    root_width: I,
    radius_range: I,
    work: usize,
}
impl AuthoredThread {
    pub fn stock_clip_plan(
        &self,
        options: StockClipOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<StockClipPlan, Error> {
        StockClipPlan::prepare(*self, options, interrupt)
    }
}
impl StockClipPlan {
    fn prepare(
        source: AuthoredThread,
        options: StockClipOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        let mut work = Meter::new(options, size_of::<Self>(), interrupt)?;
        work.charge(64)?;
        let AuthoredSection::NptFlat {
            depth_mm,
            flat_width_mm,
        } = source.spec.section
        else {
            return Err(Error::Unsupported);
        };
        let Some(limit) = source.spec.crest_radius_limit_mm else {
            return Err(Error::Unsupported);
        };
        if source.spec.radial_taper <= 0. {
            return Err(Error::Unsupported);
        }
        let r = I::point(source.spec.major_radius_mm);
        let k = I::point(source.spec.radial_taper);
        // Both sealed radii are <=1e7 and k is exactly +1/32; the
        // complete event enclosure is finite with magnitude <3.3e8.
        let event = (I::point(limit) - r) / k;
        if event.hi <= 0. || event.lo >= source.spec.length_mm {
            return Err(Error::Unsupported);
        }
        if event.lo <= 0. || (I::point(source.spec.length_mm) - event).lo <= 0. {
            return Err(Error::Unresolved);
        }
        let depth = I::point(depth_mm);
        let half = I::point(flat_width_mm) / I::point(2.);
        let m = I::point(3.).sqrt();
        let slopes = [I::point(1.) / (m - k), I::point(1.) / (m + k)];
        let root_radius = r - depth;
        let root_end = root_radius + k * I::point(source.spec.length_mm);
        if (I::point(limit) - root_end).lo <= source.minimum {
            return Err(Error::Unresolved);
        }
        // Reuse the ORIGINAL shared endpoint expressions. The sealed NPT
        // constructor already proved this IDENTICAL root-width enclosure
        // exceeds its positive minimum, so no second fallible admission exists.
        let q0 = source.endpoint(0, false).expect("NPT has its lower band").0;
        let q3 = source.endpoint(2, true).expect("NPT has its upper band").0;
        let root_width = I::point(source.spec.pitch_mm) - (q3 - q0);
        // For every post-event height, qLo-q0=(stock-root(z))/(m-k),
        // q3-qHi=(stock-root(z))/(m+k), qHi-qLo=w+(R+k*z-stock)*S.
        // These are positive; the original root width closes one full period.
        // With r>0 and z increasing in both slabs, (theta,z) has one radial
        // sheet. The two cap loops and fixed opposite uses close that SOURCE
        // atlas. This is not a literal-carrier incidence/embedding receipt.
        work.charge(64)?;
        Ok(Self {
            source,
            event,
            root_radius,
            q0,
            q3,
            half,
            width: I::point(flat_width_mm),
            slopes,
            root_width,
            radius_range: I::new(root_radius.lo, limit),
            work: work.used,
        })
    }
    pub fn source(&self) -> &AuthoredThread {
        &self.source
    }
    /// The sealed NPT plan retains this exact original point interval.
    pub(super) fn flat_width_mm(&self) -> f64 {
        self.width.lo
    }
    pub fn stock_radius_limit_mm(&self) -> f64 {
        self.source
            .spec
            .crest_radius_limit_mm
            .expect("sealed stock event")
    }
    pub fn contact_height_bounds_mm(&self) -> [f64; 2] {
        bounds(self.event)
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    /// Conservative owned-plan plus fixed inline arithmetic scratch. This is
    /// not a process/RSS bound; no heap graph or discretized source is retained.
    pub fn retained_budget_bytes(&self) -> usize {
        size_of::<Self>() + SCRATCH
    }
    pub fn patches(&self) -> [StockClipPatch; 8] {
        std::array::from_fn(|i| StockClipPatch {
            slab: if i < 4 {
                StockClipSlab::Before
            } else {
                StockClipSlab::After
            },
            band: BANDS[i % 4],
        })
    }
    pub fn vertices(&self) -> [StockClipVertex; 12] {
        std::array::from_fn(|i| StockClipVertex {
            station: [
                StockClipStation::Start,
                StockClipStation::Contact,
                StockClipStation::End,
            ][i / 4],
            rail: RAILS[i % 4],
        })
    }
    pub fn curves(&self) -> [StockClipCurve; 20] {
        std::array::from_fn(|i| {
            if i < 8 {
                StockClipCurve::Axial {
                    slab: if i < 4 {
                        StockClipSlab::Before
                    } else {
                        StockClipSlab::After
                    },
                    rail: RAILS[i % 4],
                }
            } else {
                let i = i - 8;
                StockClipCurve::Ring {
                    station: [
                        StockClipStation::Start,
                        StockClipStation::Contact,
                        StockClipStation::End,
                    ][i / 4],
                    band: BANDS[i % 4],
                }
            }
        })
    }
    pub fn cap_coedges(&self, station: StockClipStation) -> Result<[StockClipCoedge; 4], Error> {
        if station == StockClipStation::Contact {
            return Err(Error::Invalid);
        }
        Ok(std::array::from_fn(|i| StockClipCoedge {
            curve: StockClipCurve::Ring {
                station,
                band: BANDS[if station == StockClipStation::Start {
                    i
                } else {
                    3 - i
                }],
            },
            forward: station == StockClipStation::Start,
        }))
    }
    /// Distinct owners denote a physical source edge. Equal owners at the
    /// contact ring are auxiliary chart seams of one continuous source face.
    pub fn curve_faces(&self, curve: StockClipCurve) -> [StockClipFace; 2] {
        match curve {
            StockClipCurve::Axial { slab, rail } => {
                let i = rail_index(rail);
                [
                    StockClipPatch {
                        slab,
                        band: BANDS[(i + 3) % 4],
                    }
                    .face(),
                    StockClipPatch {
                        slab,
                        band: BANDS[i],
                    }
                    .face(),
                ]
            }
            StockClipCurve::Ring { station, band } => {
                let before = StockClipPatch {
                    slab: StockClipSlab::Before,
                    band,
                }
                .face();
                let after = StockClipPatch {
                    slab: StockClipSlab::After,
                    band,
                }
                .face();
                match station {
                    StockClipStation::Start => [StockClipFace::StartCap, before],
                    StockClipStation::Contact => [before, after],
                    StockClipStation::End => [after, StockClipFace::EndCap],
                }
            }
        }
    }
    fn height(&self, station: StockClipStation) -> I {
        match station {
            StockClipStation::Start => I::point(0.),
            StockClipStation::Contact => self.event,
            StockClipStation::End => I::point(self.source.spec.length_mm),
        }
    }
    fn contacts(&self, z: I, station: Option<StockClipStation>) -> [I; 2] {
        if station == Some(StockClipStation::Contact) {
            [-self.half, self.half]
        } else {
            let shift = I::point(self.stock_radius_limit_mm())
                - I::point(self.source.spec.major_radius_mm)
                - I::point(self.source.spec.radial_taper) * z;
            [
                -self.half + shift * self.slopes[0],
                self.half - shift * self.slopes[1],
            ]
        }
    }
}
