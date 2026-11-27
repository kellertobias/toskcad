//! Profile stock type picker, linked section diagram and plain names
//! (CAD-434).
//!
//! Every picture here is drawn from the canonical section the stock
//! resolves to (`ProfileSection`), never from a separate icon, so the type
//! images, the dimensioned diagram and the member geometry cannot disagree.
//! The diagram is display-only: each callout names one stock field and
//! clicking it focuses that field (the CAD-395 linked-diagram contract).
use crate::{
    egui,
    modeling::diagram::{self, Diagram, DiagramLink, Key, Line, Notation},
};
use multicad_core::{
    materials::{
        ProfileShape, StockKind, corner_radius_problem, format_mm, inner_corner_radius,
        profile_description,
    },
    model::{ProfileSection, resolve_sketch_section, sample_section_loop},
};

pub(crate) const WIDTH: Key = "stock.width";
pub(crate) const HEIGHT: Key = "stock.height";
pub(crate) const DIAMETER: Key = "stock.diameter";
pub(crate) const WALL: Key = "stock.wall";
pub(crate) const RADIUS: Key = "stock.corner_radius";

/// Shape families in picker order, with their accessible label and the
/// plain-language types each covers.
pub(crate) const FAMILIES: [(ProfileShape, &str, &str); 5] = [
    (
        ProfileShape::Rectangular,
        "Rectangular",
        "Square or rectangular tube; flat or square bar without a wall",
    ),
    (
        ProfileShape::Round,
        "Round",
        "Round tube; round bar without a wall",
    ),
    (ProfileShape::Angle, "Angle", "L section with two legs"),
    (
        ProfileShape::ItemLine8_40x40,
        "Slot profile 40×40 (item Line 8)",
        "Aluminium slot profile 40×40 with four 8 mm slots",
    ),
    (
        ProfileShape::Sketch,
        "Based on sketch",
        "Draw a custom cross-section in a sketch",
    ),
];

pub(crate) fn family_label(shape: ProfileShape) -> &'static str {
    FAMILIES
        .iter()
        .find(|(s, ..)| *s == shape)
        .map_or("Profile", |(_, label, _)| label)
}

/// The section a profile specification resolves to, without validating the
/// rest of the stock, so a draft with a typo still has a picture.
pub(crate) fn section(spec: &StockKind) -> Option<ProfileSection> {
    let StockKind::Profile {
        designation,
        length_mm,
        width_mm,
        height_mm,
        wall_mm,
        shape,
        section,
        corner_radius_mm,
    } = spec
    else {
        return None;
    };
    let positive = |v: f64| v.is_finite() && v > 0.;
    if !positive(*width_mm) || !positive(*height_mm) {
        return None;
    }
    let custom = match (shape, section) {
        (ProfileShape::Sketch, Some(section)) => {
            Some(resolve_sketch_section(&section.sketch).ok()?)
        }
        (ProfileShape::Sketch, None) => return None,
        _ => None,
    };
    let wall_mm = wall_mm.filter(|t| positive(*t) && t * 2. < width_mm.min(*height_mm));
    let corner_radius_mm = corner_radius_mm
        .filter(|r| corner_radius_problem(*shape, *width_mm, *height_mm, wall_mm, *r).is_none());
    Some(ProfileSection {
        designation: designation.clone(),
        shape: *shape,
        width_mm: *width_mm,
        height_mm: *height_mm,
        wall_mm,
        stock_length_mm: *length_mm,
        custom,
        corner_radius_mm,
    })
}

/// A representative section of a family, for the type picker.
pub(crate) fn sample(shape: ProfileShape, hollow: bool) -> ProfileSection {
    let (w, h, wall, r) = match (shape, hollow) {
        (ProfileShape::Rectangular, true) => (40., 30., Some(3.), Some(5.)),
        (ProfileShape::Rectangular, false) => (40., 12., None, None),
        (ProfileShape::Round, true) => (40., 40., Some(4.), None),
        (ProfileShape::Round, false) => (40., 40., None, None),
        (ProfileShape::Angle, _) => (40., 40., Some(6.), None),
        (ProfileShape::ItemLine8_40x40 | ProfileShape::Sketch, _) => (40., 40., None, None),
    };
    ProfileSection {
        designation: String::new(),
        shape: if shape == ProfileShape::Sketch {
            ProfileShape::Rectangular
        } else {
            shape
        },
        width_mm: w,
        height_mm: h,
        wall_mm: wall,
        stock_length_mm: 6000.,
        custom: None,
        corner_radius_mm: r,
    }
}

/// Section outlines (outside first, then holes) in section millimetres.
pub(crate) fn loops(section: &ProfileSection) -> Vec<Vec<[f32; 2]>> {
    let outer = section
        .custom
        .as_ref()
        .map_or_else(|| section.outer(0.), |c| c.outer.clone());
    std::iter::once(outer)
        .chain(section.holes())
        .map(|boundary| {
            sample_section_loop(&boundary)
                .into_iter()
                .map(|[x, y]| [x as f32, y as f32])
                .collect()
        })
        .collect()
}

/// The section outline alone: the picture of a type.
pub(crate) fn art(section: &ProfileSection) -> Diagram {
    loops(section)
        .into_iter()
        .fold(Diagram::new(), |d, points| {
            d.path(&points, true, Line::Outline)
        })
}

/// Picture of a "based on sketch" family: a free outline with a sketch hint.
fn sketch_art() -> Diagram {
    Diagram::new()
        .path(
            &[
                [0., 0.],
                [40., 0.],
                [40., 14.],
                [26., 30.],
                [8., 30.],
                [0., 18.],
            ],
            true,
            Line::Outline,
        )
        .circle([18., 14.], 5., Line::Outline, false)
}

pub(crate) fn family_art(shape: ProfileShape, hollow: bool) -> Diagram {
    if shape == ProfileShape::Sketch {
        sketch_art()
    } else {
        art(&sample(shape, hollow))
    }
}

/// The editable fields a specification shows, in display order. Their keys
/// are exactly the callouts of [`diagram`].
pub(crate) fn keys(spec: &StockKind) -> Vec<Key> {
    let StockKind::Profile { shape, wall_mm, .. } = spec else {
        return vec![];
    };
    match shape {
        ProfileShape::Rectangular if wall_mm.is_some() => vec![WIDTH, HEIGHT, WALL, RADIUS],
        ProfileShape::Rectangular => vec![WIDTH, HEIGHT, RADIUS],
        ProfileShape::Round if wall_mm.is_some() => vec![DIAMETER, WALL],
        ProfileShape::Round => vec![DIAMETER],
        ProfileShape::Angle => vec![WIDTH, HEIGHT, WALL],
        ProfileShape::ItemLine8_40x40 | ProfileShape::Sketch => vec![],
    }
}

/// The dimensioned section: every editable field of the type is a callout.
pub(crate) fn diagram(spec: &StockKind) -> Diagram {
    let Some(section) = section(spec) else {
        return Diagram::new()
            .outline_rect(0., 0., 40., 30., Line::Hidden)
            .note([20., 15.], "Enter a positive size");
    };
    let StockKind::Profile {
        shape,
        wall_mm,
        corner_radius_mm,
        ..
    } = spec
    else {
        return Diagram::new();
    };
    let (w, h) = (section.width_mm as f32, section.height_mm as f32);
    let s = w.max(h);
    let keys = keys(spec);
    let mut d = art(&section);
    if keys.contains(&WIDTH) {
        d = d.dim(WIDTH, [0., 0.], [w, 0.], -14., Notation::Length, "Width");
    }
    if keys.contains(&DIAMETER) {
        d = d
            .line([-0.1 * s, h / 2.], [w + 0.1 * s, h / 2.], Line::Center)
            .line([w / 2., -0.1 * s], [w / 2., h + 0.1 * s], Line::Center)
            .dim(
                DIAMETER,
                [0., 0.],
                [w, 0.],
                -14.,
                Notation::Diameter,
                "Outside diameter",
            );
    }
    if keys.contains(&HEIGHT) {
        d = d.dim(HEIGHT, [w, 0.], [w, h], -14., Notation::Length, "Height");
    }
    if keys.contains(&WALL) {
        let t = wall_mm.unwrap_or(0.) as f32;
        let (at, name) = if *shape == ProfileShape::Angle {
            ([t / 2., h * 0.7], "Leg thickness")
        } else {
            ([t / 2., h * 0.55], "Wall")
        };
        d = d.leader(WALL, at, [-0.55 * s, h + 0.2 * s], Notation::Length, name);
    }
    if keys.contains(&RADIUS) {
        let r = corner_radius_mm.unwrap_or(0.) as f32;
        let c = std::f32::consts::FRAC_1_SQRT_2;
        let at = [w - r + r * c, h - r + r * c];
        d = d.leader(
            RADIUS,
            at,
            [w + 0.15 * s, h + 0.45 * s],
            Notation::Radius,
            "Corner radius",
        );
    }
    if *shape == ProfileShape::ItemLine8_40x40 {
        d = d.note([w / 2., -0.2 * s], "40 × 40 mm · four 8 mm slots");
    }
    d
}

/// The text a callout shows for its field.
pub(crate) fn callout_value(spec: &StockKind, key: Key) -> String {
    let StockKind::Profile {
        width_mm,
        height_mm,
        wall_mm,
        corner_radius_mm,
        ..
    } = spec
    else {
        return "?".into();
    };
    match key {
        WIDTH | DIAMETER => format!("{} mm", format_mm(*width_mm)),
        HEIGHT => format!("{} mm", format_mm(*height_mm)),
        WALL => wall_mm.map_or("?".into(), |t| format!("{} mm", format_mm(t))),
        RADIUS => corner_radius_mm.map_or("sharp".into(), |r| format!("{} mm", format_mm(r))),
        _ => "?".into(),
    }
}

/// The first field the specification's own validation refers to.
pub(crate) fn invalid_key(spec: &StockKind) -> Option<Key> {
    let StockKind::Profile {
        width_mm,
        height_mm,
        wall_mm,
        shape,
        corner_radius_mm,
        ..
    } = spec
    else {
        return None;
    };
    let positive = |v: f64| v.is_finite() && v > 0.;
    let keys = keys(spec);
    if !positive(*width_mm) {
        return keys.first().copied();
    }
    if !positive(*height_mm) && keys.contains(&HEIGHT) {
        return Some(HEIGHT);
    }
    if wall_mm.is_some_and(|t| !positive(t) || t * 2. > width_mm.min(*height_mm)) {
        return Some(WALL);
    }
    corner_radius_mm
        .filter(|r| corner_radius_problem(*shape, *width_mm, *height_mm, *wall_mm, *r).is_some())
        .map(|_| RADIUS)
}

/// Paint a section picture into a fixed-size cell.
pub(crate) fn picture(ui: &mut egui::Ui, size: egui::Vec2, art: &Diagram) -> egui::Response {
    let (rect, response) = ui.allocate_exact_size(size, egui::Sense::hover());
    if ui.is_rect_visible(rect) {
        ui.painter()
            .rect_filled(rect, 3., ui.visuals().extreme_bg_color);
        diagram::thumbnail(ui, rect, art);
    }
    response
}

/// One picker entry: the family's pictures (tube and bar where both exist),
/// a selectable label and the plain types it covers.
fn family_row(ui: &mut egui::Ui, shape: &mut ProfileShape, family: ProfileShape) -> bool {
    let (_, label, hint) = FAMILIES
        .iter()
        .find(|(s, ..)| *s == family)
        .copied()
        .expect("known family");
    let mut changed = false;
    ui.horizontal(|ui| {
        let size = egui::vec2(34., 26.);
        picture(ui, size, &family_art(family, true));
        if matches!(family, ProfileShape::Rectangular | ProfileShape::Round) {
            picture(ui, size, &family_art(family, false));
        }
        ui.vertical(|ui| {
            changed |= ui
                .selectable_value(shape, family, label)
                .on_hover_text(hint)
                .changed();
            ui.label(egui::RichText::new(hint).small().weak());
        });
    });
    changed
}

/// The type picker: a picture of the current section beside a dropdown of
/// pictured families. Returns whether the family changed.
pub(crate) fn type_picker(ui: &mut egui::Ui, spec: &StockKind, shape: &mut ProfileShape) -> bool {
    let before = *shape;
    ui.label("Section type");
    ui.horizontal(|ui| {
        let current = match section(spec) {
            Some(section) if section.custom.is_some() || *shape != ProfileShape::Sketch => {
                art(&section)
            }
            _ => family_art(*shape, true),
        };
        picture(ui, egui::vec2(40., 30.), &current);
        egui::ComboBox::from_id_salt("profile-shape")
            .selected_text(family_label(*shape))
            .width((ui.available_width() - 8.).clamp(80., 260.))
            .show_ui(ui, |ui| {
                ui.set_min_width(240.);
                for (family, ..) in FAMILIES {
                    family_row(ui, shape, family);
                }
            });
    });
    *shape != before
}

/// The current field being edited, for the label highlight.
fn label_text(ui: &egui::Ui, link: &DiagramLink, key: Key, label: &str) -> egui::RichText {
    let emphasis = link.emphasis();
    let text = egui::RichText::new(label);
    if emphasis.invalid == Some(key) {
        text.color(ui.visuals().error_fg_color).strong()
    } else if emphasis.active == Some(key) {
        text.color(ui.visuals().selection.stroke.color).strong()
    } else {
        text
    }
}

/// A millimetre field linked to its callout.
pub(crate) fn linked_number(
    ui: &mut egui::Ui,
    link: &mut DiagramLink,
    key: Key,
    label: &str,
    value: &mut f64,
    help: &str,
) -> bool {
    ui.horizontal_wrapped(|ui| {
        let label_response = ui.label(label_text(ui, link, key, label));
        let response = ui
            .add(
                egui::DragValue::new(value)
                    .speed(0.1)
                    .range(0.0..=f64::MAX)
                    .suffix(" mm"),
            )
            .labelled_by(label_response.id)
            .on_hover_text(help);
        link.track(key, &response);
        response.changed()
    })
    .inner
}

/// The stored link state of the stock form, kept between frames.
pub(crate) fn link_state(ui: &egui::Ui) -> DiagramLink {
    ui.data(|d| d.get_temp::<DiagramLink>(link_id()))
        .unwrap_or_default()
}
pub(crate) fn store_link(ui: &egui::Ui, link: DiagramLink) {
    ui.data_mut(|d| d.insert_temp(link_id(), link));
}
fn link_id() -> egui::Id {
    egui::Id::new("profile-stock-diagram-link")
}

/// Paint the plain description of a profile, if it has one.
pub(crate) fn description(spec: &StockKind) -> Option<String> {
    profile_description(spec)
}

/// The inside corner radius note for a rounded tube.
pub(crate) fn inside_radius_note(spec: &StockKind) -> Option<String> {
    let StockKind::Profile {
        shape: ProfileShape::Rectangular,
        wall_mm: Some(t),
        corner_radius_mm: Some(r),
        ..
    } = spec
    else {
        return None;
    };
    Some(match inner_corner_radius(Some(*r), Some(*t)) {
        Some(inner) => format!(
            "Inside radius R{} mm (corner radius − wall)",
            format_mm(inner)
        ),
        None => "Inside corners stay sharp (corner radius ≤ wall)".into(),
    })
}

/// A profile stock row for pickers in Model tools: its picture, plain
/// description and, as secondary text, its own name.
pub(crate) fn stock_row(
    ui: &mut egui::Ui,
    selected: bool,
    spec: &StockKind,
    name: &str,
) -> egui::Response {
    ui.horizontal(|ui| {
        match section(spec) {
            Some(section) => picture(ui, egui::vec2(30., 22.), &art(&section)),
            None => picture(ui, egui::vec2(30., 22.), &Diagram::new()),
        };
        ui.vertical(|ui| {
            let text = description(spec).unwrap_or_else(|| name.to_owned());
            let response = ui.selectable_label(selected, text);
            ui.label(egui::RichText::new(name).small().weak());
            response
        })
        .inner
    })
    .inner
}
