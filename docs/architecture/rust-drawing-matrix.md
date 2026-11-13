# Rust drawing capability matrix

This matrix pins every inventoried combination of a Rust body family or two-body scene with a drawing view. ✓ means the combination draws well-formed exact vectors. F*n* is an explicit unresolved finding: the drawing refuses with the listed Unsupported message. `crates/multicad-kernel/tests/drawing_inventory.rs` asserts every cell; regenerate the table only together with a reviewed change. Scenes place the second body 20 mm above the first and either overlapping in plan (2 mm offset) or separated (80 mm offset).

Views are taken through the body center. Axial looks along +Z, side along −Y, and oblique along −(1,1,1)/√3; sections cut through the center (the oblique section-only cut uses the oblique plane), finite depth retains 1 mm, the lateral window is 4 × 4 mm, and the detail circle has a 3 mm radius.

Analytic and OpenCascade-referenced regressions for the combinations added by CAD-346 are in `rust_bored_prism_drawing.rs`, `rust_gear_section_drawing.rs`, `rust_separated_scene_drawing.rs` and `rust_oblique_bore_drawing.rs`. The oblique bore family is the CAD-354 audit geometry with a 0.5 mm bore along the cube diagonal.

| Family or scene | ax | ax rev | ax hid | side | side hid | obl | obl hid | ax sec | ax only | ax depth | side sec | side only | side depth | side win | obl sec | obl only | ax det | side det |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| box | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| bored box | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | F1 | ✓ | ✓ | ✓ |
| oblique bore | F2 | F2 | F2 | F2 | F2 | F2 | F2 | F3 | ✓ | F3 | F3 | ✓ | F3 | F3 | F4 | ✓ | F2 | F2 |
| cylinder | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | F1 | ✓ | ✓ | ✓ |
| annulus | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | F5 | ✓ | F5 | F5 | F5 | ✓ | ✓ | ✓ |
| sphere | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| rounded slot | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | F1 | ✓ | F1 | ✓ | ✓ | ✓ | ✓ | F1 | F6 | ✓ | ✓ |
| spur gear | ✓ | ✓ | ✓ | F7 | F7 | F7 | F7 | ✓ | ✓ | ✓ | F8 | F8 | F8 | F8 | F8 | F8 | F7 | F7 |
| internal ring | ✓ | ✓ | ✓ | F9 | F9 | F9 | F9 | ✓ | ✓ | ✓ | F8 | F8 | F8 | F8 | F8 | F8 | F9 | F9 |
| helical gear | ✓ | ✓ | ✓ | F10 | F10 | F10 | F10 | ✓ | ✓ | ✓ | F10 | F11 | F10 | F10 | F10 | F11 | F10 | F10 |
| herringbone gear | ✓ | ✓ | ✓ | F10 | F10 | F10 | F10 | ✓ | ✓ | ✓ | F10 | F11 | F10 | F10 | F10 | F11 | F10 | F10 |
| bevel gear | ✓ | ✓ | ✓ | F12 | F12 | F12 | F12 | ✓ | ✓ | ✓ | F12 | F12 | F12 | F12 | F12 | F12 | F12 | F12 |
| za worm | ✓ | ✓ | ✓ | F13 | F13 | F13 | F13 | ✓ | ✓ | ✓ | F13 | F13 | F13 | F13 | F13 | F13 | F13 | F13 |
| worm wheel | ✓ | ✓ | ✓ | F13 | F13 | F13 | F13 | ✓ | ✓ | ✓ | F13 | F13 | F13 | F13 | F13 | F13 | F13 | F13 |
| box + cylinder overlapping | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| box + cylinder separated | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| box + sphere overlapping | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| box + sphere separated | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| box + spur gear overlapping | F14 | F14 | F14 | F7 | F7 | F14 | F14 | F14 | F14 | F14 | F8 | F8 | F8 | F8 | F14 | F14 | F14 | F7 |
| box + spur gear separated | ✓ | ✓ | ✓ | F7 | F7 | F7 | F7 | ✓ | ✓ | ✓ | F8 | F8 | F8 | F8 | F8 | F8 | F7 | F7 |
| cylinder + sphere overlapping | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| cylinder + sphere separated | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| cylinder + annulus overlapping | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | F1 | ✓ | ✓ | ✓ |
| cylinder + annulus separated | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | F1 | ✓ | ✓ | ✓ |
| box + internal ring overlapping | F14 | F14 | F14 | F9 | F9 | F14 | F14 | F14 | F14 | F14 | F8 | F8 | F8 | F8 | F14 | F14 | F14 | F9 |
| box + internal ring separated | ✓ | ✓ | ✓ | F9 | F9 | F9 | F9 | ✓ | ✓ | ✓ | F8 | F8 | F8 | F8 | F8 | F8 | F9 | F9 |
| spur gear + spur gear overlapping | F14 | F14 | F14 | F7 | F7 | F14 | F14 | F14 | F14 | F14 | F8 | F8 | F8 | F8 | F14 | F14 | F14 | F7 |
| spur gear + spur gear separated | ✓ | ✓ | ✓ | F7 | F7 | F7 | F7 | ✓ | ✓ | ✓ | F8 | F8 | F8 | F8 | F8 | F8 | F7 | F7 |
| box + bevel gear overlapping | F15 | F15 | F15 | F12 | F12 | F12 | F12 | F15 | F15 | F15 | F12 | F12 | F12 | F12 | F12 | F12 | F15 | F12 |
| box + bevel gear separated | ✓ | ✓ | ✓ | F12 | F12 | F12 | F12 | ✓ | ✓ | ✓ | F12 | F12 | F12 | F12 | F12 | F12 | F12 | F12 |
| box + za worm overlapping | F16 | F16 | F16 | F13 | F13 | F16 | F16 | F16 | F16 | F16 | F13 | F13 | F13 | F13 | F16 | F16 | F16 | F13 |
| box + za worm separated | ✓ | ✓ | ✓ | F13 | F13 | F13 | F13 | ✓ | ✓ | ✓ | F13 | F13 | F13 | F13 | F13 | F13 | F13 | F13 |
| box + herringbone gear overlapping | F17 | F17 | F17 | F10 | F10 | F17 | F17 | F17 | ✓ | F17 | F10 | F11 | F10 | F10 | F17 | F11 | F17 | F10 |
| box + herringbone gear separated | ✓ | ✓ | ✓ | F10 | F10 | F10 | F10 | ✓ | ✓ | ✓ | F10 | F11 | F10 | F10 | F10 | F11 | F10 | F10 |
| bored box + rounded slot overlapping | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| bored box + rounded slot separated | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |

## Findings

- **F1** (7 cells): Rust retained projection cannot yet clip this curved boundary
- **F2** (9 cells): Rust oblique bore drawings currently support sections of one solid or separated scenes
- **F3** (5 cells): Rust retained oblique bore sections require a view along the bore axis
- **F4** (1 cells): Rust retained oblique bore sections require the retained stock behind the cut face
- **F5** (4 cells): A Rust annulus section with geometry behind the cut requires a transverse plane and no lateral crop
- **F6** (1 cells): Rust machined sections require a stock-aligned cutting plane
- **F7** (24 cells): Rust spur gear drawing currently requires an uncropped axial view
- **F8** (42 cells): Rust spur and internal ring sections require an uncropped plane perpendicular to the gear axis
- **F9** (15 cells): Rust internal ring drawing currently requires an uncropped axial view
- **F10** (36 cells): Rust herringbone projection requires an uncropped axial view
- **F11** (8 cells): Rust herringbone sections require a plane perpendicular to the gear axis
- **F12** (35 cells): Rust bevel drawing requires an uncropped view along the gear axis
- **F13** (43 cells): Rust worm and worm-wheel drawings require an uncropped view along the gear axis
- **F14** (33 cells): Rust planar drawing projection requires straight edges and planar faces
- **F15** (7 cells): Rust bevel drawing supports one gear at a time; scene occlusion is unavailable
- **F16** (11 cells): Rust worm and worm-wheel drawings support one placed gear at a time
- **F17** (9 cells): Rust exact drawing projection for continuous herringbone gears is unavailable

358 of 648 combinations draw; 290 are explicit findings. No combination returns an invalid or partial result.
