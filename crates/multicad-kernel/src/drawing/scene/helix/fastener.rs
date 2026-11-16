//! Scene members of the plain helical thread, the internal thread cut in
//! cylindrical stock and the threaded hex nut, read from the kernel's exact
//! read-only descriptors (`rust_helical_thread_data`,
//! `rust_internal_thread_data`, `rust_hex_nut_geometry`).
//!
//! All three reuse the thread member: a plain thread is the helicoid's
//! material between its caps; an internal thread is the stock cylinder less
//! the helicoid's void, with the full stock disc as the floor of a
//! partial-depth cut; a nut is the hexagonal prism less the void. Only the
//! stock's reach and chords and the floor plane need their own answers here.
use super::*;
use multicad_rust_kernel::hex_nut::Nut;

fn circle(radius_mm: f64) -> Vec<Edge> {
    vec![Edge::Planar(normal_sweep::Edge::Circle {
        center: [0.; 2],
        radius_mm,
        clockwise: false,
    })]
}

/// The member of one of the three families, or `None` for other solids.
pub(super) fn build(
    solid: &Solid,
    view: &View,
    tolerance: f64,
) -> Option<Result<Thread, GeometryError>> {
    let (thread, frame, stock, floor, hex) =
        if let Some((thread, frame)) = solid.rust_helical_thread_data() {
            (thread, frame, None, None, None)
        } else if let Some(data) = solid.rust_internal_thread_data() {
            let length = data.thread.spec().length_mm;
            let top = data.stock_length_mm;
            if !(length > 0. && length <= top) {
                return Some(Err(err(
                    GeometryErrorKind::InvalidResult,
                    "Rust internal thread void exceeds its stock",
                )));
            }
            let band = Band {
                rings: vec![circle(data.stock_radius_mm)],
                axis: [0., top],
            };
            (
                data.thread,
                data.frame,
                Some((data.stock_radius_mm, band)),
                (top > length).then_some(top),
                None,
            )
        } else if let Some((nut, frame)) = solid.rust_hex_nut_geometry() {
            let band = Band {
                rings: vec![super::hexagon(nut.across_flats_mm())],
                axis: [0., nut.thread().spec().length_mm],
            };
            // The circumradius: the stock's largest reach.
            let reach = nut.across_flats_mm() / 3_f64.sqrt();
            (*nut.thread(), frame, Some((reach, band)), None, Some(*nut))
        } else {
            return None;
        };
    let Some(profile) = Profile::from_helical(&thread) else {
        return Some(Err(err(
            GeometryErrorKind::InvalidResult,
            "Rust thread profile cannot be read for drawing",
        )));
    };
    let origin = view.point(frame.origin);
    let axes = [frame.u, frame.v, frame.normal()].map(|axis| view.vector(axis));
    let spec = thread.spec();
    Some((|| {
        let stock = match stock {
            Some((radius, band)) => {
                Some((radius, Prism::new(origin, axes, vec![band], tolerance)?))
            }
            None => None,
        };
        Ok(Thread {
            law: Arc::new(Law {
                origin,
                axes,
                profile,
                rate: thread.helical_rate(),
                starts: spec.starts,
            }),
            z: [0., spec.length_mm],
            neck: None,
            stock,
            cone: None,
            base: None,
            groove: None,
            hex,
            floor,
            tolerance,
        })
    })())
}

impl Thread {
    /// Whether `level` is the floor of a partial-depth internal thread.
    pub(super) fn is_floor(&self, level: f64) -> bool {
        self.floor.is_some() && (level - self.z[1]).abs() <= self.tolerance
    }

    /// The stock's reach at `level` along the local angle: the nut's
    /// hexagon, else the stock cylinder's (or a groove's tail at its
    /// shoulder).
    pub(super) fn stock_reach(&self, level: f64, angle: f64) -> f64 {
        match &self.hex {
            Some(nut) => nut.outer_radius(angle),
            None => self.cap_reach(level),
        }
    }

    /// The station interval, along the cap chord direction, of the stock's
    /// cap at `level` cut by the plane `n . p = c`.
    pub(super) fn stock_chord(&self, n: P3, c: f64, level: f64) -> Option<[f64; 2]> {
        let Some(nut) = &self.hex else {
            return self.disc_chord(n, c, level, self.cap_reach(level));
        };
        let nxy = n[0].hypot(n[1]);
        let d = [n[1] / nxy, -n[0] / nxy];
        let offset = (c - n[2] * level) / (nxy * nxy);
        let s = [n[0] * offset, n[1] * offset];
        hex_chord(nut, s, d).map(|[a, b]| {
            let center = s[0] * d[0] + s[1] * d[1];
            [center + a, center + b]
        })
    }

    /// Exact bounds of the hexagonal or floored stock, which contains the
    /// whole member.
    pub(super) fn stock_bounds(&self) -> Option<[P3; 2]> {
        if self.hex.is_none() && self.floor.is_none() {
            return None;
        }
        self.stock.as_ref().map(|(_, prism)| prism.bounds())
    }

    /// A line in the floor plane of a partial-depth internal thread: the
    /// material on both sides there is the stock less the void's cap.
    pub(super) fn floor_line(&self, s: P3, d: P3) -> Option<Result<Vec<[f64; 2]>, GeometryError>> {
        if !(d[2].abs() <= 1e-15 && self.is_floor(s[2])) {
            return None;
        }
        let (_, prism) = self.stock.as_ref()?;
        Some((|| {
            let outer = match prism.line(self.law.view(s), self.view_vector(d)) {
                Some(outer) => outer?,
                None => Vec::new(),
            };
            let void = self.cap_line(s, d, self.z[1])?;
            Ok(profile::exclusive(&outer, &void, 0.))
        })())
    }
}

/// Parameters `t`, relative to the foot `s`, where `s + t d` lies inside the
/// nut's hexagon: the six flats `x cos(k pi/3) + y sin(k pi/3) <= w / 2`.
fn hex_chord(nut: &Nut, s: [f64; 2], d: [f64; 2]) -> Option<[f64; 2]> {
    let half = 0.5 * nut.across_flats_mm();
    let [mut lo, mut hi] = [f64::NEG_INFINITY, f64::INFINITY];
    for k in 0..6 {
        let (sin, cos) = (f64::from(k) * PI / 3.).sin_cos();
        let value = s[0] * cos + s[1] * sin;
        let rate = d[0] * cos + d[1] * sin;
        if rate.abs() <= 1e-15 {
            if value >= half {
                return None;
            }
            continue;
        }
        let t = (half - value) / rate;
        if rate > 0. {
            hi = hi.min(t);
        } else {
            lo = lo.max(t);
        }
    }
    (hi > lo).then_some([lo, hi])
}
