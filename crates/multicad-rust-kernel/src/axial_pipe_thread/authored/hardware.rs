//! Actual authored pipe hardware domains. The thread law and its axial span
//! stay immutable; hex stock is added at the large end or surrounds a bore.
//! Source material receipts here do not admit a mesh or a STEP file.
use super::*;
use std::mem::size_of;
pub mod mesh;
mod moments;
pub mod shell;
use moments::Moments;
#[cfg(test)]
mod tests;
const FIXED: usize = 128 * 1024;
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Kind {
    /// Fully threaded: local small tip is z=0, neck is the original z=L.
    HexHead {
        across_flats_mm: f64,
        height_mm: f64,
    },
    /// A through bore uses its original mouth z=0 and negative taper.
    HexNut { across_flats_mm: f64 },
}
#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub limits: StockClipOptions,
    pub volume_width_mm3: f64,
    pub cap_area_width_mm2: f64,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            limits: StockClipOptions::default(),
            volume_width_mm3: 1e-6,
            cap_area_width_mm2: 1e-7,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PlaneKind {
    ThreadTip,
    HeadShoulder,
    HeadTop,
    NutStart,
    NutEnd,
}
#[derive(Clone, Copy, Debug)]
pub struct Plane {
    pub kind: PlaneKind,
    pub height_bounds_mm: [f64; 2],
    pub outward_sign: i8,
    pub area_mm2: [f64; 2],
}
#[derive(Clone, Copy, Debug)]
pub enum Planes {
    Head([Plane; 3]),
    Nut([Plane; 2]),
}
impl Planes {
    pub fn as_slice(&self) -> &[Plane] {
        match self {
            Self::Head(p) => p,
            Self::Nut(p) => p,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Location {
    Inside,
    Outside,
    Unresolved,
}
#[derive(Clone, Copy, Debug)]
pub struct Query {
    location: Location,
    work: usize,
    inherited: usize,
    bytes: usize,
}
impl Query {
    pub fn location(&self) -> Location {
        self.location
    }
    pub fn query_work_used(&self) -> usize {
        self.work
    }
    pub fn inherited_preparation_work_used(&self) -> usize {
        self.inherited
    }
    pub fn accounted_peak_bytes(&self) -> usize {
        self.bytes
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Material {
    thread: AuthoredThread,
    kind: Kind,
    clip: Option<StockClipPlan>,
    volume: [f64; 2],
    thread_volume: [f64; 2],
    caps: [[f64; 2]; 2],
    hex_area: [f64; 2],
    planes: Planes,
    radius: [f64; 2],
    bounds: [[f64; 2]; 3],
    work: usize,
    bytes: usize,
}
struct Work<'a> {
    left: usize,
    used: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Work<'_> {
    fn charge(&mut self, n: usize) -> Result<(), Error> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        self.left = self.left.checked_sub(n).ok_or(Error::Limit)?;
        self.used = self.used.checked_add(n).ok_or(Error::Limit)?;
        Ok(())
    }
}
fn bounded(v: I, width: f64) -> Result<[f64; 2], Error> {
    if !v.valid() || v.lo <= 0. || v.width() > width {
        return Err(Error::Unresolved);
    }
    Ok(bounds(v))
}
fn clipped_area(plan: &StockClipPlan, z: f64) -> I {
    // The sealed plan earns k>0 and limit-b(L)>minimum>0 in prepare.
    // Since b(z)=r-d+kz increases on[0,L], BOTH caps have b(z)<limit:
    // there is no fully-clipped-root branch in this source family. Derive
    // both inputs from that exact plan rather than accepting unrelated data.
    let thread = plan.source();
    let limit = plan.stock_radius_limit_mm();
    let r = I::point(thread.spec.major_radius_mm);
    let d = I::point(thread.depth());
    let k = I::point(thread.spec.radial_taper);
    let w = I::point(plan.flat_width_mm());
    let m = I::point(3.).sqrt();
    let slopes = I::point(1.) / (m - k) + I::point(1.) / (m + k);
    let root = I::point(thread.spec.pitch_mm) - w - d * slopes;
    let b = r - d + k * I::point(z);
    let natural = root * b.square()
        + w * (b + d).square()
        + d * slopes * (b.square() + b * d + d.square() / I::point(3.));
    let limit = I::point(limit);
    let c = limit - b;
    let clipped = root * b.square()
        + slopes * c * (b.square() + limit * b + limit.square()) / I::point(3.)
        + (w + (d - c) * slopes) * limit.square();
    let phase = if (b + d).hi <= limit.lo {
        natural
    } else if (b + d).lo >= limit.hi {
        clipped
    } else {
        I::new(natural.lo.min(clipped.lo), natural.hi.max(clipped.hi))
    };
    phase * pi() / I::point(thread.spec.pitch_mm)
}
impl Material {
    pub fn new(
        thread: AuthoredThread,
        kind: Kind,
        options: Options,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        if options.limits.maximum_work == 0
            || options.limits.maximum_work > 16_000_000
            || options.limits.maximum_bytes == 0
            || options.limits.maximum_bytes > 64 * 1024 * 1024
            || [options.volume_width_mm3, options.cap_area_width_mm2]
                .into_iter()
                .any(|v| !v.is_finite() || v <= 0.)
        {
            return Err(Error::Invalid);
        }
        let bytes = size_of::<Self>() + FIXED;
        if bytes > options.limits.maximum_bytes {
            return Err(Error::Limit);
        }
        let mut work = Work {
            left: options.limits.maximum_work,
            used: 0,
            interrupt,
        };
        work.charge(4096)?;
        let (across, height, head) = match kind {
            Kind::HexHead {
                across_flats_mm,
                height_mm,
            } => (across_flats_mm, height_mm, true),
            Kind::HexNut { across_flats_mm } => (across_flats_mm, 0., false),
        };
        if !across.is_finite() || across <= 0. || !height.is_finite() || (head && height <= 0.) {
            return Err(Error::Invalid);
        }
        if across > 1e7
            || height > 1e7
            || (I::point(thread.spec.length_mm) + I::point(height)).hi > 1e7
        {
            return Err(Error::Limit);
        }
        if (head && thread.spec.radial_taper < 0.) || (!head && thread.spec.radial_taper > 0.) {
            return Err(Error::Unsupported);
        }
        let radius =
            thread.radial_bounds([0., thread.spec.pitch_mm], [0., thread.spec.length_mm])?;
        if (I::point(across) / I::point(2.) - I::new(radius[0], radius[1])).lo <= thread.minimum {
            return Err(Error::Unresolved);
        }
        let hex = I::point(across).square() * I::point(3.).sqrt() / I::point(2.);
        let mut clip = None;
        if work.left == 0 {
            return Err(Error::Limit);
        }
        let (volume, caps) = match thread.stock_clip_plan(
            StockClipOptions {
                maximum_work: work.left,
                ..options.limits
            },
            interrupt,
        ) {
            Ok(plan) => {
                work.charge(plan.work_used())?;
                let receipt = plan.volume_bounds(
                    options.volume_width_mm3,
                    StockClipOptions {
                        maximum_work: work
                            .left
                            .checked_add(plan.work_used())
                            .ok_or(Error::Limit)?,
                        ..options.limits
                    },
                    interrupt,
                )?;
                work.charge(receipt.work_used() - plan.work_used())?;
                let caps = [
                    clipped_area(&plan, 0.),
                    clipped_area(&plan, thread.spec.length_mm),
                ];
                let volume = I::new(
                    receipt.source_volume_mm3()[0],
                    receipt.source_volume_mm3()[1],
                );
                clip = Some(plan);
                (volume, caps)
            }
            Err(Error::Unsupported) => {
                // The existing clip probe charges64 before reporting that
                // its source/event family is unsupported. Preserve that work.
                work.charge(64)?;
                // Do not erase an active ceiling merely because this source
                // family has no moving-clip provider. The common radius
                // envelope earns an inactive ceiling before using moments.
                if thread.spec.crest_radius_limit_mm.is_some_and(|r| {
                    (I::point(thread.spec.major_radius_mm)
                        + I::point(thread.spec.radial_taper) * I::new(0., thread.spec.length_mm))
                    .hi > r
                }) {
                    return Err(Error::Unsupported);
                }
                let moments = Moments::new(&thread);
                (
                    moments.volume(&thread),
                    [
                        moments.area(&thread, 0.),
                        moments.area(&thread, thread.spec.length_mm),
                    ],
                )
            }
            Err(cause) => return Err(cause),
        };
        work.charge(4096)?;
        let thread_volume = bounded(volume, options.volume_width_mm3)?;
        let caps = [
            bounded(caps[0], options.cap_area_width_mm2)?,
            bounded(caps[1], options.cap_area_width_mm2)?,
        ];
        let material = if head {
            volume + hex * I::point(height)
        } else {
            hex * I::point(thread.spec.length_mm) - volume
        };
        let volume = bounded(material, options.volume_width_mm3)?;
        let hex_area = bounded(hex, options.cap_area_width_mm2)?;
        let exposed = [
            bounded(
                hex - I::new(caps[0][0], caps[0][1]),
                options.cap_area_width_mm2,
            )?,
            bounded(
                hex - I::new(caps[1][0], caps[1][1]),
                options.cap_area_width_mm2,
            )?,
        ];
        let circum = I::point(across) / I::point(3.).sqrt();
        let end = I::point(thread.spec.length_mm) + I::point(height);
        if head && (end - I::point(thread.spec.length_mm)).lo <= thread.minimum {
            return Err(Error::Unresolved);
        }
        let plane = |kind, height: I, outward_sign, area_mm2| Plane {
            kind,
            height_bounds_mm: bounds(height),
            outward_sign,
            area_mm2,
        };
        let planes = if head {
            Planes::Head([
                plane(PlaneKind::ThreadTip, I::point(0.), -1, caps[0]),
                plane(
                    PlaneKind::HeadShoulder,
                    I::point(thread.spec.length_mm),
                    -1,
                    exposed[1],
                ),
                plane(PlaneKind::HeadTop, end, 1, hex_area),
            ])
        } else {
            Planes::Nut([
                plane(PlaneKind::NutStart, I::point(0.), -1, exposed[0]),
                plane(
                    PlaneKind::NutEnd,
                    I::point(thread.spec.length_mm),
                    1,
                    exposed[1],
                ),
            ])
        };
        Ok(Self {
            thread,
            kind,
            clip,
            volume,
            thread_volume,
            caps,
            hex_area,
            planes,
            radius,
            bounds: [
                [-circum.hi, circum.hi],
                [-across / 2., across / 2.],
                [0., end.hi],
            ],
            work: work.used,
            bytes,
        })
    }
    pub fn thread(&self) -> &AuthoredThread {
        &self.thread
    }
    pub fn kind(&self) -> Kind {
        self.kind
    }
    pub fn clipped_thread(&self) -> Option<&StockClipPlan> {
        self.clip.as_ref()
    }
    pub fn volume_mm3(&self) -> [f64; 2] {
        self.volume
    }
    pub fn original_thread_volume_mm3(&self) -> [f64; 2] {
        self.thread_volume
    }
    pub fn original_cap_areas_mm2(&self) -> [[f64; 2]; 2] {
        self.caps
    }
    pub fn hex_area_mm2(&self) -> [f64; 2] {
        self.hex_area
    }
    pub fn radial_enclosure_mm(&self) -> [f64; 2] {
        self.radius
    }
    pub fn bounds_mm(&self) -> [[f64; 2]; 3] {
        self.bounds
    }
    pub fn preparation_work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        size_of::<Self>()
    }
    pub fn accounted_peak_bytes(&self) -> usize {
        self.bytes
    }
    pub fn planes(&self) -> &Planes {
        &self.planes
    }
    /// Whole Cartesian-box enclosure using the ORIGINAL radial law. Near a
    /// threaded flank the full-phase enclosure may be unresolved; that is not
    /// a false outside/inside or a sampled classification. Each new query has
    /// its own caller allowance; inherited immutable preparation is explicit.
    pub fn classify_box(
        &self,
        query: [[f64; 2]; 3],
        limits: StockClipOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Query, Error> {
        if limits.maximum_work == 0
            || limits.maximum_work > 16_000_000
            || limits.maximum_bytes == 0
            || limits.maximum_bytes > 64 * 1024 * 1024
            || query
                .iter()
                .any(|q| !q[0].is_finite() || !q[1].is_finite() || q[0] > q[1])
        {
            return Err(Error::Invalid);
        }
        if limits.maximum_bytes < self.bytes {
            return Err(Error::Limit);
        }
        let mut work = Work {
            left: limits.maximum_work,
            used: 0,
            interrupt,
        };
        work.charge(4096)?;
        let result = |location, work: Work<'_>| Query {
            location,
            work: work.used,
            inherited: self.work,
            bytes: self.bytes,
        };
        if query
            .iter()
            .zip(self.bounds)
            .any(|(q, b)| q[1] < b[0] || q[0] > b[1])
        {
            return Ok(result(Location::Outside, work));
        }
        let z = I::new(query[2][0], query[2][1]);
        let end = match self.kind {
            Kind::HexHead { height_mm, .. } => {
                I::point(self.thread.spec.length_mm) + I::point(height_mm)
            }
            Kind::HexNut { .. } => I::point(self.thread.spec.length_mm),
        };
        if z.lo <= 0. || z.hi >= end.lo {
            return Ok(result(Location::Unresolved, work));
        }
        let x = I::new(query[0][0], query[0][1]);
        let y = I::new(query[1][0], query[1][1]);
        let across = match self.kind {
            Kind::HexHead {
                across_flats_mm, ..
            }
            | Kind::HexNut { across_flats_mm } => I::point(across_flats_mm),
        };
        let absolute = |v: I| {
            I::new(
                if v.lo <= 0. && v.hi >= 0. {
                    0.
                } else {
                    v.lo.abs().min(v.hi.abs())
                },
                v.magnitude(),
            )
        };
        let ay = absolute(y);
        let diagonal = I::point(3.).sqrt() * absolute(x) + ay;
        // Exact regular-hex halfspaces for vertices (+/-AF/sqrt3,0)
        // and (+/-AF/(2sqrt3),+/-AF/2); no rounded trig normal is authority.
        if ay.lo > (across / I::point(2.)).hi || diagonal.lo > across.hi {
            return Ok(result(Location::Outside, work));
        }
        let hex_inside = ay.hi < (across / I::point(2.)).lo && diagonal.hi < across.lo;
        if matches!(self.kind, Kind::HexHead { .. }) && z.lo > self.thread.spec.length_mm {
            return Ok(result(
                if hex_inside {
                    Location::Inside
                } else {
                    Location::Unresolved
                },
                work,
            ));
        }
        // Construction already earned this identical immutable full-phase
        // calculation on[0,L]. The guards above restrict z to a subinterval;
        // it only narrows k*z in radial_bounds, whose other expressions and
        // branches depend on the same source/phase constants. There is no
        // callback or resource operation to lose at this sealed handoff.
        let radial = self
            .thread
            .radial_bounds(
                [0., self.thread.spec.pitch_mm],
                [z.lo, z.hi.min(self.thread.spec.length_mm)],
            )
            .expect("material earned the full-phase source on the containing z domain");
        let r = x.square() + y.square();
        let inner = I::point(radial[0].max(0.)).square();
        let outer = I::point(radial[1]).square();
        let location = match self.kind {
            Kind::HexHead { .. } if r.hi < inner.lo => Location::Inside,
            Kind::HexHead { .. } if z.hi < self.thread.spec.length_mm && r.lo > outer.hi => {
                Location::Outside
            }
            Kind::HexNut { .. } if r.hi < inner.lo => Location::Outside,
            Kind::HexNut { .. } if hex_inside && r.lo > outer.hi => Location::Inside,
            _ => Location::Unresolved,
        };
        Ok(result(location, work))
    }
}
