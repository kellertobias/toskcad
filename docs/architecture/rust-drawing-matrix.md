# Rust drawing capability matrix

This matrix pins every inventoried combination of a Rust body family or two-body scene with a drawing view. ✓ means the combination draws well-formed exact vectors. F*n* is an explicit unresolved finding: the drawing refuses with the listed Unsupported message. `crates/multicad-kernel/tests/drawing_inventory.rs` asserts every cell; regenerate the table only together with a reviewed change. Scenes place the second body 20 mm above the first and either overlapping in plan (2 mm offset) or separated (80 mm offset).

Views are taken through the body center. Axial looks along +Z, side along −Y, and oblique along −(1,1,1)/√3; sections cut through the center, finite depth retains 1 mm, the lateral window is 4 × 4 mm, and the detail circle has a 3 mm radius.

Analytic and OpenCascade-referenced regressions for the combinations added by CAD-346 are in `rust_bored_prism_drawing.rs`, `rust_gear_section_drawing.rs` and `rust_separated_scene_drawing.rs`.

| Family or scene | ax | ax rev | ax hid | side | side hid | obl | obl hid | ax sec | ax only | ax depth | side sec | side only | side depth | side win | obl sec | ax det | side det |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| box | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| bored box | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | F1 | ✓ | ✓ |
| cylinder | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | F1 | ✓ | ✓ |
| annulus | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | F2 | ✓ | F2 | F2 | F2 | ✓ | ✓ |
| sphere | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| rounded slot | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | F1 | ✓ | F1 | ✓ | ✓ | ✓ | ✓ | F1 | ✓ | ✓ |
| spur gear | ✓ | ✓ | ✓ | F3 | F3 | F3 | F3 | ✓ | ✓ | ✓ | F4 | F4 | F4 | F4 | F4 | F3 | F3 |
| internal ring | ✓ | ✓ | ✓ | F5 | F5 | F5 | F5 | ✓ | ✓ | ✓ | F4 | F4 | F4 | F4 | F4 | F5 | F5 |
| helical gear | ✓ | ✓ | ✓ | F6 | F6 | F6 | F6 | ✓ | ✓ | ✓ | F6 | F7 | F6 | F6 | F6 | F6 | F6 |
| herringbone gear | ✓ | ✓ | ✓ | F6 | F6 | F6 | F6 | ✓ | ✓ | ✓ | F6 | F7 | F6 | F6 | F6 | F6 | F6 |
| bevel gear | ✓ | ✓ | ✓ | F8 | F8 | F8 | F8 | ✓ | ✓ | ✓ | F8 | F8 | F8 | F8 | F8 | F8 | F8 |
| za worm | ✓ | ✓ | ✓ | F9 | F9 | F9 | F9 | ✓ | ✓ | ✓ | F9 | F9 | F9 | F9 | F9 | F9 | F9 |
| worm wheel | ✓ | ✓ | ✓ | F9 | F9 | F9 | F9 | ✓ | ✓ | ✓ | F9 | F9 | F9 | F9 | F9 | F9 | F9 |
| box + cylinder overlapping | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| box + cylinder separated | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| box + sphere overlapping | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| box + sphere separated | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| box + spur gear overlapping | F10 | F10 | F10 | F3 | F3 | F10 | F10 | F10 | F10 | F10 | F4 | F4 | F4 | F4 | F10 | F10 | F3 |
| box + spur gear separated | ✓ | ✓ | ✓ | F3 | F3 | F3 | F3 | ✓ | ✓ | ✓ | F4 | F4 | F4 | F4 | F4 | F3 | F3 |
| cylinder + sphere overlapping | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| cylinder + sphere separated | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| cylinder + annulus overlapping | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | F1 | ✓ | ✓ |
| cylinder + annulus separated | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | F1 | ✓ | ✓ |
| box + internal ring overlapping | F10 | F10 | F10 | F5 | F5 | F10 | F10 | F10 | F10 | F10 | F4 | F4 | F4 | F4 | F10 | F10 | F5 |
| box + internal ring separated | ✓ | ✓ | ✓ | F5 | F5 | F5 | F5 | ✓ | ✓ | ✓ | F4 | F4 | F4 | F4 | F4 | F5 | F5 |
| spur gear + spur gear overlapping | F10 | F10 | F10 | F3 | F3 | F10 | F10 | F10 | F10 | F10 | F4 | F4 | F4 | F4 | F10 | F10 | F3 |
| spur gear + spur gear separated | ✓ | ✓ | ✓ | F3 | F3 | F3 | F3 | ✓ | ✓ | ✓ | F4 | F4 | F4 | F4 | F4 | F3 | F3 |
| box + bevel gear overlapping | F11 | F11 | F11 | F8 | F8 | F8 | F8 | F11 | F11 | F11 | F8 | F8 | F8 | F8 | F8 | F11 | F8 |
| box + bevel gear separated | ✓ | ✓ | ✓ | F8 | F8 | F8 | F8 | ✓ | ✓ | ✓ | F8 | F8 | F8 | F8 | F8 | F8 | F8 |
| box + za worm overlapping | F12 | F12 | F12 | F9 | F9 | F12 | F12 | F12 | F12 | F12 | F9 | F9 | F9 | F9 | F12 | F12 | F9 |
| box + za worm separated | ✓ | ✓ | ✓ | F9 | F9 | F9 | F9 | ✓ | ✓ | ✓ | F9 | F9 | F9 | F9 | F9 | F9 | F9 |
| box + herringbone gear overlapping | F13 | F13 | F13 | F6 | F6 | F13 | F13 | F13 | ✓ | F13 | F6 | F7 | F6 | F6 | F13 | F13 | F6 |
| box + herringbone gear separated | ✓ | ✓ | ✓ | F6 | F6 | F6 | F6 | ✓ | ✓ | ✓ | F6 | F7 | F6 | F6 | F6 | F6 | F6 |
| bored box + rounded slot overlapping | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | F14 | ✓ | ✓ |
| bored box + rounded slot separated | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | F14 | ✓ | ✓ |

## Findings

- **F1** (7 cells): Rust retained projection cannot yet clip this curved boundary
- **F2** (4 cells): A Rust annulus section with geometry behind the cut requires a transverse plane and no lateral crop
- **F3** (24 cells): Rust spur gear drawing currently requires an uncropped axial view
- **F4** (37 cells): Rust spur and internal ring sections require an uncropped plane perpendicular to the gear axis
- **F5** (15 cells): Rust internal ring drawing currently requires an uncropped axial view
- **F6** (36 cells): Rust herringbone projection requires an uncropped axial view
- **F7** (4 cells): Rust herringbone sections require a plane perpendicular to the gear axis
- **F8** (32 cells): Rust bevel drawing requires an uncropped view along the gear axis
- **F9** (40 cells): Rust worm and worm-wheel drawings require an uncropped view along the gear axis
- **F10** (30 cells): Rust planar drawing projection requires straight edges and planar faces
- **F11** (7 cells): Rust bevel drawing supports one gear at a time; scene occlusion is unavailable
- **F12** (10 cells): Rust worm and worm-wheel drawings support one placed gear at a time
- **F13** (9 cells): Rust exact drawing projection for continuous herringbone gears is unavailable
- **F14** (2 cells): Rust machined sections require a stock-aligned cutting plane

338 of 595 combinations draw; 257 are explicit findings. No combination returns an invalid or partial result.
