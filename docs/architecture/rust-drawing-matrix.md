# Rust drawing capability matrix

This matrix pins every inventoried combination of a Rust body family or two-body scene with a drawing view. ✓ means the combination draws well-formed exact vectors. F*n* is an explicit unresolved finding: the drawing refuses with the listed Unsupported message. `crates/multicad-kernel/tests/drawing_inventory.rs` asserts every cell; regenerate the table only together with a reviewed change. Scenes place the second body 20 mm above the first and either overlapping in plan (2 mm offset) or separated (80 mm offset).

Views are taken through the body center. Axial looks along +Z, side along −Y, and oblique along −(1,1,1)/√3; sections cut through the center (the oblique section-only cut uses the oblique plane), finite depth retains 1 mm, the lateral window is 4 × 4 mm, and the detail circle has a 3 mm radius.

Analytic and OpenCascade-referenced regressions for the combinations added by CAD-346 are in `rust_bored_prism_drawing.rs`, `rust_gear_section_drawing.rs`, `rust_separated_scene_drawing.rs` and `rust_oblique_bore_drawing.rs`. The oblique bore family is the CAD-354 audit geometry with a 0.5 mm bore along the cube diagonal.

| Family or scene | ax | ax rev | ax hid | side | side hid | obl | obl hid | ax sec | ax only | ax depth | side sec | side only | side depth | side win | obl sec | obl only | ax det | side det |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| box | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| bored box | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| oblique bore | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| cylinder | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| annulus | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| sphere | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| rounded slot | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| spur gear | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| internal ring | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| helical gear | ✓ | ✓ | ✓ | F10 | F10 | F10 | F10 | ✓ | ✓ | ✓ | F10 | F11 | F10 | F10 | F10 | F11 | F10 | F10 |
| herringbone gear | ✓ | ✓ | ✓ | F10 | F10 | F10 | F10 | ✓ | ✓ | ✓ | F10 | F11 | F10 | F10 | F10 | F11 | F10 | F10 |
| bevel gear | ✓ | ✓ | ✓ | F12 | F12 | F12 | F12 | ✓ | ✓ | ✓ | F12 | F12 | F12 | F12 | F12 | F12 | F12 | F12 |
| za worm | ✓ | ✓ | ✓ | F13 | F13 | F13 | F13 | ✓ | ✓ | ✓ | F13 | F13 | F13 | F13 | F13 | F13 | F13 | F13 |
| worm wheel | ✓ | ✓ | ✓ | F13 | F13 | F13 | F13 | ✓ | ✓ | ✓ | F13 | F13 | F13 | F13 | F13 | F13 | F13 | F13 |
| box + cylinder overlapping | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| box + cylinder separated | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| box + sphere overlapping | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| box + sphere separated | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| box + spur gear overlapping | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| box + spur gear separated | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| cylinder + sphere overlapping | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| cylinder + sphere separated | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| cylinder + annulus overlapping | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| cylinder + annulus separated | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| box + internal ring overlapping | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| box + internal ring separated | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| spur gear + spur gear overlapping | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| spur gear + spur gear separated | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| box + bevel gear overlapping | F15 | F15 | F15 | F12 | F12 | F12 | F12 | F15 | F15 | F15 | F12 | F12 | F12 | F12 | F12 | F12 | F15 | F12 |
| box + bevel gear separated | ✓ | ✓ | ✓ | F12 | F12 | F12 | F12 | ✓ | ✓ | ✓ | F12 | F12 | F12 | F12 | F12 | F12 | F12 | F12 |
| box + za worm overlapping | F16 | F16 | F16 | F13 | F13 | F16 | F16 | F16 | F16 | F16 | F13 | F13 | F13 | F13 | F16 | F16 | F16 | F13 |
| box + za worm separated | ✓ | ✓ | ✓ | F13 | F13 | F13 | F13 | ✓ | ✓ | ✓ | F13 | F13 | F13 | F13 | F13 | F13 | F13 | F13 |
| box + herringbone gear overlapping | F17 | F17 | F17 | F10 | F10 | F17 | F17 | F17 | ✓ | F17 | F10 | F11 | F10 | F10 | F17 | F11 | F17 | F10 |
| box + herringbone gear separated | ✓ | ✓ | ✓ | F10 | F10 | F10 | F10 | ✓ | ✓ | ✓ | F10 | F11 | F10 | F10 | F10 | F11 | F10 | F10 |
| bored box + rounded slot overlapping | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| bored box + rounded slot separated | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |

## Findings

- **F10** (36 cells): Rust herringbone projection requires an uncropped axial view
- **F11** (8 cells): Rust herringbone sections require a plane perpendicular to the gear axis
- **F12** (35 cells): Rust bevel drawing requires an uncropped view along the gear axis
- **F13** (43 cells): Rust worm and worm-wheel drawings require an uncropped view along the gear axis
- **F15** (7 cells): Rust bevel drawing supports one gear at a time; scene occlusion is unavailable
- **F16** (11 cells): Rust worm and worm-wheel drawings support one placed gear at a time
- **F17** (9 cells): Rust exact drawing projection for continuous herringbone gears is unavailable

499 of 648 combinations draw; 149 are explicit findings. No combination returns an invalid or partial result.
