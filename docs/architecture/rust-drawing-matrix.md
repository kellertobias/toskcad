# Rust drawing capability matrix

This matrix pins every inventoried combination of a Rust body family or two-body scene with a drawing view. ✓ means the combination draws well-formed exact vectors. F*n* is an explicit unresolved finding: the drawing refuses with the listed Unsupported message. `crates/multicad-kernel/tests/drawing_inventory.rs` asserts every cell; regenerate the table only together with a reviewed change. Scenes place the second body 20 mm above the first and either overlapping in plan (2 mm offset) or separated (80 mm offset).

Views are taken through the body center. Axial looks along +Z, side along −Y, and oblique along −(1,1,1)/√3; sections cut through the center (the oblique section-only cut uses the oblique plane), finite depth retains 1 mm, the lateral window is 4 × 4 mm, and the detail circle has a 3 mm radius.

Analytic and OpenCascade-referenced regressions for the combinations added by CAD-346 are in `rust_bored_prism_drawing.rs`, `rust_gear_section_drawing.rs`, `rust_separated_scene_drawing.rs` and `rust_oblique_bore_drawing.rs`. The oblique bore family is the CAD-354 audit geometry with a 0.5 mm bore along the cube diagonal.

Requests that the family providers refuse draw through the exact-curve scene engine in `crates/multicad-kernel/src/drawing/scene/` when every body of the scene has a scene member: right prisms over line, arc and involute profiles (boxes, cylinders, annuli, rounded slots, bored stock, spur gears with hubs, internal rings) and the box with an oblique bore. The engine clips each member's exact edges, outlines and section loops to the retained region, splits them where their projections meet, classifies every piece by the exact material along its sight line and publishes cut faces, window faces, detail crops and hatches uniformly. Its analytic regressions are `rust_scene_gear.rs`, `rust_scene_prism.rs` and `rust_scene_bore.rs`; `rust_scene_native.rs` and `rust_oblique_bore_drawing.rs` compare every class of view, section and overlapping scene with OpenCascade 7.9.3 (developer builds with `--features occt` only). Measured OpenCascade deviations on the STEP-exchanged gears, none of which widen the three-chord-tolerance bound: flank sections bridged by straight chords (sagitta up to 0.104 mm), rear root-fillet outlines and smooth flank seams published as visible behind the root cylinder, and silhouettes that fall on a face seam.

Spherical bevel gears are scene members too (`scene/bevel.rs`): every face is ruled (cones through the pitch apex over the root, fillet, meridian, spherical-involute and tip pieces of the tooth boundary, two back cones and the bore), so outlines are whole rulings, a sight line meets the solid where closed-form quadric roots or certified sign changes of the provider's boundary law occur, and a plane cuts each face strip between the crossings of its two rims. Uncropped axial views keep the dedicated bevel provider; every other bevel view, crop, section and scene draws through the engine. `rust_scene_bevel.rs` holds the analytic regressions and compares 60 view, section and scene cases with the independently constructed OpenCascade 7.9.3 bevel: curves agree both ways within 0.0122 mm (bound 0.03 mm). Measured OpenCascade deviations: it refuses oblique cuts and hidden-line lengthwise cuts of its own bevel ("Projected curve family is unsupported"), bridges transverse flank cuts with straight chords (sagitta up to 0.174 mm), draws tangent-continuous flank seams, and publishes up to 15.4 mm of rear tooth edges per oblique hidden-line view as visible although material lies in front of them (decided by marching the exact membership predicate in the member's unit tests). No Rust visible curve lies on a native hidden one. Edges of one body that pierce another interpenetrating body are not split at the piercing point by the scene engine; the inventoried scenes do not interpenetrate.

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
| bevel gear | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| thread end relief | F18 | F18 | F18 | F18 | F18 | F18 | F18 | F18 | F18 | F18 | F18 | F18 | F18 | F18 | F18 | F18 | F18 | F18 |
| rounded thread | F18 | F18 | F18 | F18 | F18 | F18 | F18 | F18 | F18 | F18 | F18 | F18 | F18 | F18 | F18 | F18 | F18 | F18 |
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
| box + bevel gear overlapping | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| box + bevel gear separated | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| box + za worm overlapping | F16 | F16 | F16 | F13 | F13 | F16 | F16 | F16 | F16 | F16 | F13 | F13 | F13 | F13 | F16 | F16 | F16 | F13 |
| box + za worm separated | ✓ | ✓ | ✓ | F13 | F13 | F13 | F13 | ✓ | ✓ | ✓ | F13 | F13 | F13 | F13 | F13 | F13 | F13 | F13 |
| box + herringbone gear overlapping | F17 | F17 | F17 | F10 | F10 | F17 | F17 | F17 | ✓ | F17 | F10 | F11 | F10 | F10 | F17 | F11 | F17 | F10 |
| box + herringbone gear separated | ✓ | ✓ | ✓ | F10 | F10 | F10 | F10 | ✓ | ✓ | ✓ | F10 | F11 | F10 | F10 | F10 | F11 | F10 | F10 |
| bored box + rounded slot overlapping | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| bored box + rounded slot separated | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |

## Findings

- **F10** (36 cells): Rust herringbone projection requires an uncropped axial view
- **F11** (8 cells): Rust herringbone sections require a plane perpendicular to the gear axis
- **F13** (43 cells): Rust worm and worm-wheel drawings require an uncropped view along the gear axis
- **F16** (11 cells): Rust worm and worm-wheel drawings support one placed gear at a time
- **F17** (9 cells): Rust exact drawing projection for continuous herringbone gears is unavailable
- **F18** (36 cells): Rust planar drawing projection requires straight edges and planar faces

541 of 684 combinations draw; 143 are explicit findings. No combination returns an invalid or partial result.
