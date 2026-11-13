# Rust drawing capability matrix

This matrix pins every inventoried combination of a Rust body family or two-body scene with a drawing view. ✓ means the combination draws well-formed exact vectors. F*n* is an explicit unresolved finding: the drawing refuses with the listed Unsupported message. `crates/multicad-kernel/tests/drawing_inventory.rs` asserts every cell; regenerate the table only together with a reviewed change. Scenes place the second body 20 mm above the first and either overlapping in plan (2 mm offset) or separated (80 mm offset).

Views are taken through the body center. Axial looks along +Z, side along −Y, and oblique along −(1,1,1)/√3; sections cut through the center (the oblique section-only cut uses the oblique plane), finite depth retains 1 mm, the lateral window is 4 × 4 mm, and the detail circle has a 3 mm radius.

Analytic and OpenCascade-referenced regressions for the combinations added by CAD-346 are in `rust_bored_prism_drawing.rs`, `rust_gear_section_drawing.rs` and `rust_separated_scene_drawing.rs`.

| Family or scene | ax | ax rev | ax hid | side | side hid | obl | obl hid | ax sec | ax only | ax depth | side sec | side only | side depth | side win | obl sec | obl only | ax det | side det |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| box | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| bored box | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | F1 | ✓ | ✓ | ✓ |
| cylinder | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | F1 | ✓ | ✓ | ✓ |
| annulus | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | F2 | ✓ | F2 | F2 | F2 | ✓ | ✓ | ✓ |
| sphere | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| rounded slot | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | F1 | ✓ | F1 | ✓ | ✓ | ✓ | ✓ | F1 | F3 | ✓ | ✓ |
| spur gear | ✓ | ✓ | ✓ | F4 | F4 | F4 | F4 | ✓ | ✓ | ✓ | F5 | F5 | F5 | F5 | F5 | F5 | F4 | F4 |
| internal ring | ✓ | ✓ | ✓ | F6 | F6 | F6 | F6 | ✓ | ✓ | ✓ | F5 | F5 | F5 | F5 | F5 | F5 | F6 | F6 |
| helical gear | ✓ | ✓ | ✓ | F7 | F7 | F7 | F7 | ✓ | ✓ | ✓ | F7 | F8 | F7 | F7 | F7 | F8 | F7 | F7 |
| herringbone gear | ✓ | ✓ | ✓ | F7 | F7 | F7 | F7 | ✓ | ✓ | ✓ | F7 | F8 | F7 | F7 | F7 | F8 | F7 | F7 |
| bevel gear | ✓ | ✓ | ✓ | F9 | F9 | F9 | F9 | ✓ | ✓ | ✓ | F9 | F9 | F9 | F9 | F9 | F9 | F9 | F9 |
| za worm | ✓ | ✓ | ✓ | F10 | F10 | F10 | F10 | ✓ | ✓ | ✓ | F10 | F10 | F10 | F10 | F10 | F10 | F10 | F10 |
| worm wheel | ✓ | ✓ | ✓ | F10 | F10 | F10 | F10 | ✓ | ✓ | ✓ | F10 | F10 | F10 | F10 | F10 | F10 | F10 | F10 |
| box + cylinder overlapping | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| box + cylinder separated | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| box + sphere overlapping | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| box + sphere separated | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| box + spur gear overlapping | F11 | F11 | F11 | F4 | F4 | F11 | F11 | F11 | F11 | F11 | F5 | F5 | F5 | F5 | F11 | F11 | F11 | F4 |
| box + spur gear separated | ✓ | ✓ | ✓ | F4 | F4 | F4 | F4 | ✓ | ✓ | ✓ | F5 | F5 | F5 | F5 | F5 | F5 | F4 | F4 |
| cylinder + sphere overlapping | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| cylinder + sphere separated | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| cylinder + annulus overlapping | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | F1 | ✓ | ✓ | ✓ |
| cylinder + annulus separated | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | F1 | ✓ | ✓ | ✓ |
| box + internal ring overlapping | F11 | F11 | F11 | F6 | F6 | F11 | F11 | F11 | F11 | F11 | F5 | F5 | F5 | F5 | F11 | F11 | F11 | F6 |
| box + internal ring separated | ✓ | ✓ | ✓ | F6 | F6 | F6 | F6 | ✓ | ✓ | ✓ | F5 | F5 | F5 | F5 | F5 | F5 | F6 | F6 |
| spur gear + spur gear overlapping | F11 | F11 | F11 | F4 | F4 | F11 | F11 | F11 | F11 | F11 | F5 | F5 | F5 | F5 | F11 | F11 | F11 | F4 |
| spur gear + spur gear separated | ✓ | ✓ | ✓ | F4 | F4 | F4 | F4 | ✓ | ✓ | ✓ | F5 | F5 | F5 | F5 | F5 | F5 | F4 | F4 |
| box + bevel gear overlapping | F12 | F12 | F12 | F9 | F9 | F9 | F9 | F12 | F12 | F12 | F9 | F9 | F9 | F9 | F9 | F9 | F12 | F9 |
| box + bevel gear separated | ✓ | ✓ | ✓ | F9 | F9 | F9 | F9 | ✓ | ✓ | ✓ | F9 | F9 | F9 | F9 | F9 | F9 | F9 | F9 |
| box + za worm overlapping | F13 | F13 | F13 | F10 | F10 | F13 | F13 | F13 | F13 | F13 | F10 | F10 | F10 | F10 | F13 | F13 | F13 | F10 |
| box + za worm separated | ✓ | ✓ | ✓ | F10 | F10 | F10 | F10 | ✓ | ✓ | ✓ | F10 | F10 | F10 | F10 | F10 | F10 | F10 | F10 |
| box + herringbone gear overlapping | F14 | F14 | F14 | F7 | F7 | F14 | F14 | F14 | ✓ | F14 | F7 | F8 | F7 | F7 | F14 | F8 | F14 | F7 |
| box + herringbone gear separated | ✓ | ✓ | ✓ | F7 | F7 | F7 | F7 | ✓ | ✓ | ✓ | F7 | F8 | F7 | F7 | F7 | F8 | F7 | F7 |
| bored box + rounded slot overlapping | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| bored box + rounded slot separated | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |

## Findings

- **F1** (7 cells): Rust retained projection cannot yet clip this curved boundary
- **F2** (4 cells): A Rust annulus section with geometry behind the cut requires a transverse plane and no lateral crop
- **F3** (1 cells): Rust machined sections require a stock-aligned cutting plane
- **F4** (24 cells): Rust spur gear drawing currently requires an uncropped axial view
- **F5** (42 cells): Rust spur and internal ring sections require an uncropped plane perpendicular to the gear axis
- **F6** (15 cells): Rust internal ring drawing currently requires an uncropped axial view
- **F7** (36 cells): Rust herringbone projection requires an uncropped axial view
- **F8** (8 cells): Rust herringbone sections require a plane perpendicular to the gear axis
- **F9** (35 cells): Rust bevel drawing requires an uncropped view along the gear axis
- **F10** (43 cells): Rust worm and worm-wheel drawings require an uncropped view along the gear axis
- **F11** (33 cells): Rust planar drawing projection requires straight edges and planar faces
- **F12** (7 cells): Rust bevel drawing supports one gear at a time; scene occlusion is unavailable
- **F13** (11 cells): Rust worm and worm-wheel drawings support one placed gear at a time
- **F14** (9 cells): Rust exact drawing projection for continuous herringbone gears is unavailable

355 of 630 combinations draw; 275 are explicit findings. No combination returns an invalid or partial result.
