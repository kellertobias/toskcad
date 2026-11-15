# Rust drawing capability matrix

This matrix pins every inventoried combination of a Rust body family or two-body scene with a drawing view. ✓ means the combination draws well-formed exact vectors. F*n* is an explicit unresolved finding: the drawing refuses with the listed Unsupported message. `crates/multicad-kernel/tests/drawing_inventory.rs` asserts every cell; regenerate the table only together with a reviewed change. Scenes place the second body 20 mm above the first and either overlapping in plan (2 mm offset) or separated (80 mm offset).

Views are taken through the body center. Axial looks along +Z, side along −Y, and oblique along −(1,1,1)/√3; sections cut through the center (the oblique section-only cut uses the oblique plane), finite depth retains 1 mm, the lateral window is 4 × 4 mm, and the detail circle has a 3 mm radius.

Analytic and OpenCascade-referenced regressions for the combinations added by CAD-346 are in `rust_bored_prism_drawing.rs`, `rust_gear_section_drawing.rs`, `rust_separated_scene_drawing.rs` and `rust_oblique_bore_drawing.rs`. The oblique bore family is the CAD-354 audit geometry with a 0.5 mm bore along the cube diagonal.

Requests that the family providers refuse draw through the exact-curve scene engine in `crates/multicad-kernel/src/drawing/scene/` when every body of the scene has a scene member: right prisms over line, arc and involute profiles (boxes, cylinders, annuli, rounded slots, bored stock, spur gears with hubs, internal rings) and the box with an oblique bore. The engine clips each member's exact edges, outlines and section loops to the retained region, splits them where their projections meet, classifies every piece by the exact material along its sight line and publishes cut faces, window faces, detail crops and hatches uniformly. Its analytic regressions are `rust_scene_gear.rs`, `rust_scene_prism.rs` and `rust_scene_bore.rs`; `rust_scene_native.rs` and `rust_oblique_bore_drawing.rs` compare every class of view, section and overlapping scene with OpenCascade 7.9.3 (developer builds with `--features occt` only). Measured OpenCascade deviations on the STEP-exchanged gears, none of which widen the three-chord-tolerance bound: flank sections bridged by straight chords (sagitta up to 0.104 mm), rear root-fillet outlines and smooth flank seams published as visible behind the root cylinder, and silhouettes that fall on a face seam.

Helical and continuous herringbone gears are scene members of their own (`scene/twisted.rs`): the uncropped axial views stay with the axial provider, and every other view, crop, section and overlapping scene is evaluated on the phase law of the transverse profile. Rims, vertex helices and the herringbone apex ring are exact edges; involute helicoids and coaxial cylinders have straight outlines and the radial flanks and fillets traced ones; sight intervals end on isolated roots of the flank surfaces and plane sections are traced on every flank. The member is certified against an independent exact membership predicate in its unit tests (interval ends within 1e-13 of the surface relative to the gear, outline normals within 2e-12 of perpendicular to the sight, section areas within 6e-4 of a scan-line quadrature), and `rust_scene_twisted.rs` certifies drawings analytically (the flank law of a section through the axis to 2e-14 mm, every occluding contour of the exact depth map within 0.025 mm of a visible curve, exact hiding behind a box) and against OpenCascade 7.9.3 for 84 view classes (Rust and native curves within 0.03 mm, section curves within 0.01 mm). Measured OpenCascade deviations: imported flank seams drawn as sharp edges, rear edges and a detail view's bore outline published as visible behind material, a stretch of a finite-depth far cut omitted on one tooth of the hubbed gears (accepted only where the exact predicate certifies the Rust curve), and the native error "Drawing curve has no finite trim domain" on some axial herringbone sections. A sight line or section plane that touches a flank along a whole curve (a view exactly along an involute ruling, a plane tangent to a flank) is refused with a typed error instead of a guessed result.

Spherical bevel gears are scene members too (`scene/bevel.rs`): every face is ruled (cones through the pitch apex over the root, fillet, meridian, spherical-involute and tip pieces of the tooth boundary, two back cones and the bore), so outlines are whole rulings, a sight line meets the solid where closed-form quadric roots or certified sign changes of the provider's boundary law occur, and a plane cuts each face strip between the crossings of its two rims. Uncropped axial views keep the dedicated bevel provider; every other bevel view, crop, section and scene draws through the engine. `rust_scene_bevel.rs` holds the analytic regressions and compares 60 view, section and scene cases with the independently constructed OpenCascade 7.9.3 bevel: curves agree both ways within 0.0122 mm (bound 0.03 mm). Measured OpenCascade deviations: it refuses oblique cuts and hidden-line lengthwise cuts of its own bevel ("Projected curve family is unsupported"), bridges transverse flank cuts with straight chords (sagitta up to 0.174 mm), draws tangent-continuous flank seams, and publishes up to 15.4 mm of rear tooth edges per oblique hidden-line view as visible although material lies in front of them (decided by marching the exact membership predicate in the member's unit tests). No Rust visible curve lies on a native hidden one. Edges of one body that pierce another interpenetrating body are not split at the piercing point by the scene engine; the inventoried scenes do not interpenetrate.

ZA worms and generated worm wheels are scene members as well (`scene/worm.rs`): both are slabs about their axis whose section at every height is star shaped, bounded by a closed chain of exact patches — the four helicoidal bands of the worm thread (crest, root and two ruled flanks per turn) and, per wheel sector, the two tip lands and the five generating patches of the provider (two conjugate flanks trimmed where they leave the blank, two tip-edge sweeps and the root torus). Edges are the patch junctions and both end-plane section chains; outlines are the traced zero sets of the exact normal against the sight line; a sight line ends on isolated roots of the radial margin against the provider's boundary radius, each interval decided by the provider's membership predicate; planes across the axis cut the exact section chain and every other plane is traced on each patch. Uncropped axial views without crops keep the axial worm provider. The member unit tests certify sight intervals against the provider's predicate in six views (both gears, with and without bore), outline normals within 1e-9 of perpendicular with no uncovered sign change of `n·d`, and closed, oriented sections whose areas agree with a rotated scan-line quadrature within 1e-3; the axial-plane cut of a worm reproduces its straight-sided profile. `rust_scene_worm.rs` certifies published section vertices within 5.1e-7 mm of the provider boundary (bound: twice the engine's 4e-7 mm weld), exact hiding behind a plate and crop bounds, and compares 108 single-gear and box-scene cases with OpenCascade 7.9.3: every native curve lies on a Rust curve within 0.011 mm outside the two worm axial details, and the 12 section-only cuts agree both ways within 0.0099 mm. Measured OpenCascade deviations on the STEP-exchanged worm: its hidden-line removal omits the helicoid flank silhouettes of side and oblique views and draws the bore generators as visible through the core, and its axial detail publishes the hidden root circle beyond the detail circle; the Rust-to-native role comparison of whole views is therefore reported but not asserted.

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
| helical gear | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| herringbone gear | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| bevel gear | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| za worm | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| worm wheel | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
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
| box + za worm overlapping | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| box + za worm separated | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| box + herringbone gear overlapping | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| box + herringbone gear separated | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| bored box + rounded slot overlapping | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| bored box + rounded slot separated | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |

## Findings

None.

All 648 combinations draw; no explicit findings remain. No combination returns an invalid or partial result.
