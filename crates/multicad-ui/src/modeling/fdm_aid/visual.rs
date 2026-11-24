//! FDM aid schematics (CAD-395): a side elevation above a plan of the
//! selected aid against its finished Part. The finished Part is drawn as
//! dashed context; the aid's own inputs carry the dimension callouts.
use super::super::diagram::{Diagram, Key, Line, Notation, legible};
use multicad_core::model::AidKind;

pub(in crate::modeling) const LENGTH: Key = "fdm.length";
pub(in crate::modeling) const WIDTH: Key = "fdm.width";
pub(in crate::modeling) const THICKNESS: Key = "fdm.thickness";
pub(in crate::modeling) const OVERLAP: Key = "fdm.overlap";
pub(in crate::modeling) const BASE: Key = "fdm.base";
pub(in crate::modeling) const GAP: Key = "fdm.gap";
pub(in crate::modeling) const ANGLE: Key = "fdm.angle";
pub(in crate::modeling) const ORIGIN: [Key; 3] = ["fdm.origin.x", "fdm.origin.y", "fdm.origin.z"];
pub(in crate::modeling) const END: [Key; 3] = ["fdm.end.x", "fdm.end.y", "fdm.end.z"];

/// Keys of each kind's named position section, in display order.
pub(in crate::modeling) fn position_keys(kind: AidKind) -> Vec<Key> {
    match kind {
        AidKind::Tab => vec![ORIGIN[0], ORIGIN[1], ORIGIN[2], ANGLE],
        AidKind::Bridge => ORIGIN.into_iter().chain(END).collect(),
        AidKind::Support => vec![ORIGIN[0], ORIGIN[1], ANGLE],
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub(in crate::modeling) struct Values {
    pub length: Option<f64>,
    pub width: Option<f64>,
    pub thickness: Option<f64>,
    pub overlap: Option<f64>,
    pub gap: Option<f64>,
    /// Bridge anchor distance, when both anchors resolve.
    pub span: Option<f64>,
    /// A bridge joins a second Part.
    pub second_part: bool,
}

pub(in crate::modeling) fn label(kind: AidKind) -> &'static str {
    match kind {
        AidKind::Tab => "Breakaway tab",
        AidKind::Bridge => "Thin bridge",
        AidKind::Support => "Removable support",
    }
}

pub(in crate::modeling) fn diagram(kind: AidKind, v: &Values) -> Diagram {
    match kind {
        AidKind::Tab => tab(v),
        AidKind::Bridge => bridge(v),
        AidKind::Support => support(v),
    }
}

fn tab(v: &Values) -> Diagram {
    let r = legible(v.length, 8., 0.01, 1e6);
    let o = legible(v.overlap, 2., 0.05 * r, 2. * r);
    let w = legible(v.width, 6., 0.15 * r, 4. * r);
    let t = legible(v.thickness, 0.4, 0.07 * r, r);
    let part = 1.6 * r;
    let h = 0.9 * r;
    let plan = -1.2 * r - w;
    Diagram::new()
        .line([-part, 0.], [1.3 * r, 0.], Line::Center)
        .outline_rect(-part, 0., 0., h, Line::Hidden)
        .note([-0.8 * r, 0.6 * h], "Finished Part")
        .rect(-o, 0., r, t)
        .dim(LENGTH, [0., 0.], [r, 0.], -14., Notation::Length, "Reach")
        .dim(OVERLAP, [-o, t], [0., t], 14., Notation::Length, "Overlap")
        .dim(
            THICKNESS,
            [r, 0.],
            [r, t],
            -14.,
            Notation::Length,
            "Thickness",
        )
        .leader(
            ORIGIN[0],
            [0., 0.],
            [-1.1 * r, -0.45 * r],
            Notation::Length,
            "Contact",
        )
        // Plan: the pad leaves the Part edge in the direction angle.
        .line([0., plan - w], [0., plan + w], Line::Hidden)
        .rect(-o, plan - w / 2., r, plan + w / 2.)
        .line([0., plan], [1.25 * r, plan], Line::Center)
        .dim(
            WIDTH,
            [r, plan - w / 2.],
            [r, plan + w / 2.],
            -14.,
            Notation::Length,
            "Width",
        )
        .leader(
            ANGLE,
            [0.7 * r, plan],
            [0.9 * r, plan - 1.1 * w],
            Notation::Angle,
            "Direction",
        )
        .note([0.2 * r, plan + 1.15 * w], "Plan · +Z build")
}

fn bridge(v: &Values) -> Diagram {
    let s = legible(v.span, 20., 0.01, 1e6);
    let t = legible(v.thickness, 0.4, 0.04 * s, 0.3 * s);
    let w = legible(v.width, 6., 0.05 * s, s);
    let h = 0.45 * s;
    let a = 0.4 * h;
    let plan = -0.35 * s - w;
    let mut d = Diagram::new()
        .line([-0.5 * s, 0.], [1.5 * s, 0.], Line::Center)
        .outline_rect(-0.5 * s, 0., 0., h, Line::Hidden)
        .outline_rect(s, 0., 1.5 * s, h, Line::Hidden)
        .rect(0., a, s, a + t)
        .dim(
            THICKNESS,
            [0.5 * s, a],
            [0.5 * s, a + t],
            14.,
            Notation::Length,
            "Thickness",
        )
        .leader(
            ORIGIN[0],
            [0., a + t / 2.],
            [-0.3 * s, 1.25 * h],
            Notation::Length,
            "Start",
        )
        .leader(
            END[0],
            [s, a + t / 2.],
            [1.3 * s, 1.25 * h],
            Notation::Length,
            "End",
        )
        .note(
            [0.5 * s, -0.12 * s],
            match v.span {
                Some(span) => format!("Span {span:.3} mm between anchors"),
                None => "Span between anchors".into(),
            },
        )
        .rect(0., plan - w / 2., s, plan + w / 2.)
        .dim(
            WIDTH,
            [0.5 * s, plan - w / 2.],
            [0.5 * s, plan + w / 2.],
            -14.,
            Notation::Length,
            "Width",
        );
    d = d.note(
        [1.25 * s, 0.5 * h],
        if v.second_part {
            "Second Part"
        } else {
            "Same Part"
        },
    );
    d
}

fn support(v: &Values) -> Diagram {
    let l = legible(v.length, 8., 0.01, 1e6);
    let w = legible(v.width, 6., 0.15 * l, 3. * l);
    let g = legible(v.gap, 0.2, 0.05 * l, 0.5 * l);
    // The underside height is derived natively; the drawing assumes one.
    let under = 1.2 * l;
    let plan = -0.5 * l - w;
    Diagram::new()
        .line([-0.5 * l, 0.], [1.5 * l, 0.], Line::Center)
        .outline_rect(-0.4 * l, under, 1.4 * l, under + 0.4 * l, Line::Hidden)
        .note([0.5 * l, under + 0.55 * l], "Finished Part underside")
        .rect(0., 0., l, under - g)
        .dim(
            GAP,
            [l, under - g],
            [l, under],
            -14.,
            Notation::Length,
            "Top gap",
        )
        .dim(LENGTH, [0., 0.], [l, 0.], -14., Notation::Length, "Length")
        .leader(
            BASE,
            [0.15 * l, 0.],
            [-0.35 * l, 0.45 * l],
            Notation::Length,
            "Base",
        )
        .note([0.5 * l, 0.55 * under], "Height from the underside")
        .rect(0., plan - w, l, plan)
        .dim(
            WIDTH,
            [l, plan - w],
            [l, plan],
            -14.,
            Notation::Length,
            "Width",
        )
        .leader(
            ORIGIN[0],
            [0., plan - w],
            [-0.3 * l, plan - 1.4 * w],
            Notation::Length,
            "Corner",
        )
        .leader(
            ANGLE,
            [0.6 * l, plan - w / 2.],
            [0.9 * l, plan + 0.6 * w],
            Notation::Angle,
            "Direction",
        )
}
