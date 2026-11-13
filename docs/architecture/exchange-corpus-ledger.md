# Exchange corpus ledger

This is the pinned pass/fail/unsupported ledger for the independent cross-CAD STEP corpus ([CAD-345](https://plans.tokenet.de/CAD/345)). The [manifest](../../crates/multicad-app/tests/fixtures/exchange/corpus/manifest.json) pins every fixture's producer, version, `FILE_NAME` strings, units, provenance, redistribution basis, SHA-256, geometry tags, expected Rust outcome, verdict and owner. Licensing and attribution are in the [corpus README](../../crates/multicad-app/tests/fixtures/exchange/corpus/README.md).

The machine-readable [ledger JSON](exchange-corpus-ledger.json) is the emitted output of the OpenCascade reference run. It records the measured source commit, tolerances, and every per-fixture measurement and timing. `corpus_manifest_pins_provenance_strata_and_owners` fails when the committed JSON and the manifest disagree on any fixture or outcome.

**Rerun after integration.** These results describe the CAD-345 branch, which is based on `dff6a0f8`. CAD-354 (oblique cylinder cuts) and CAD-346 (drawing) change import, Boolean and section behaviour. Rerun the commands below after both land and after every change to the reader or importer. Then update the manifest outcomes, verdicts and this ledger together.

## Result

| | Count |
| --- | ---: |
| Fixtures | 42 |
| Vendored, always run | 20 (15 authored, 5 KiCad library copies) |
| Repository fixtures, always run | 2 |
| Referenced from the installed KiCad 10.0.0 library by path and SHA-256 | 20 |
| **Pass** (full workflow) | **14** |
| **Unsupported**: an honest refusal of a family outside the contract | **18** |
| **Fail**: an incorrect refusal of a file that OpenCascade validates | **10** (6 invalid, 4 limit) |
| Unavailable on this host | 0 |
| Pass fixtures whose section is refused by drawing (CAD-346) | 4 |

A **pass** means that every check below succeeded. **Unsupported** means that the Rust importer refused a geometry family the [qualification matrix](rust-kernel-qualification.md) does not claim. **Fail** means that OpenCascade 7.9.3 reads the same source as a valid solid, so the Rust classification is wrong.

## Commands and artifacts

| Command | Scope | Artifact | Time (debug build, shared host) |
| --- | --- | --- | --- |
| `cargo test --locked -p multicad-app --test exchange_corpus` | Pinned ledger without the native reader, manifest and provenance checks, and the real executable over MCP for six fixtures | `target/cargo/tmp/exchange-corpus-ledger.json` | 9 s alone; up to 93 s under parallel load |
| `CARGO_TARGET_DIR=target/occt cargo test --locked -p multicad-app --features occt --test exchange_corpus corpus_matches_the_pinned_ledger` | The same ledger, plus OpenCascade 7.9.3 readback of every source and re-export and the deviation from OpenCascade's tessellation | `target/occt/tmp/exchange-corpus-ledger.json`. Set `MULTICAD_CORPUS_LEDGER` to write the committed JSON. | 162–243 s |
| `cargo test --locked -p multicad-app --test exchange_robustness` | Malformed, oversized, exhausting, healing, cancellation and damaged-project controls | `target/cargo/tmp/exchange-robustness-{malformed,healing,exhaustion}.json` | 13 s |
| `cargo test --locked -p multicad-app --bin multicad exchange_recovery` | Recovery copies and pinned library definitions | — | 0.2 s |
| `cargo test --locked -p multicad-app --test exchange_corpus -- --ignored kicad_library_sweep --nocapture` | Aggregate sweep over the whole KiCad library | `target/cargo/tmp/kicad-sweep.json` | 29.6 min wall (11.4 min CPU) under load |
| `node tools/verify-rust-runtime.mjs target/cargo/release/multicad` | Runtime contract of the compiled release executable | JSON report on stdout | < 1 s |

## Checks and tolerances

Each admitted fixture runs every check below. The tolerances are constants in the runner and are copied into the ledger JSON.

| Check | Tolerance |
| --- | --- |
| Analytic oracle: volume, bounds, face count and solid count | Volume within relative 1e-9; bounds within 1e-6 mm; counts exact |
| Independent Rust Part 21 readback of the source against the reconstructed solid | Volume within relative 1e-9; bounds within 1e-6 mm; members equal solids |
| Certified manufacturing mesh: each directed edge occurs once with its reverse, and the enclosed volume matches | Mesh volume within relative 2e-2; 0.05 mm linear deflection |
| Planar section at mid-height: exact line and arc lengths against the analytic perimeter where pinned | Relative 1e-9 |
| Material and section hatch: material, hatch, pinned-definition diagnostics and mass, after reopening | Exact; at least one hatch vector when the section succeeds |
| Canonical import: save the native project, delete the source, reopen, rebuild with Rust and compare the embedded source digest | Volume and bounds bit-identical; digest equal |
| Re-export, independent readback and reimport | Volume within relative 1e-9; bounds within 1e-6 mm; faces and solids exact |
| Geometric deviation of the readback mesh from the source mesh | ≤ 0.05 mm, searched within 0.2 mm |
| OpenCascade 7.9.3 source readback (`--features occt`) | Volume within relative 1e-6; bounds within 1e-3 mm |
| OpenCascade 7.9.3 re-export readback (`--features occt`) | Volume within relative 1e-6; bounds within 1e-3 mm; solids exact |
| Deviation between the OpenCascade tessellation and the Rust readback mesh (`--features occt`) | ≤ 0.1 mm, twice the deflection |
| A `fail` verdict (`--features occt`) | OpenCascade must accept the source as a valid solid |
| Real-executable MCP workflow | Same tolerances as above, plus an unchanged revision, project bytes and refused-file bytes after a refusal |

## Fixture ledger

| Fixture | Producer | Location | Tags | Rust import | Verdict | Owner |
| --- | --- | --- | --- | --- | --- | --- |
| `freecad-drilled-box` | FreeCAD | vendored | box, cylinder, trims-seams, placement, hole | pass | pass | section: CAD-346 |
| `occt-placed-box` | OpenCascade | vendored | box, placement | pass | pass | — |
| `occt-cylinder` | OpenCascade | vendored | cylinder, trims-seams, placement | pass | pass | — |
| `occt-torus` | OpenCascade | vendored | torus, revolution, trims-seams | pass | pass | — |
| `occt-box-void` | OpenCascade | vendored | box, voids | pass | pass | — |
| `occt-pocketed-plate` | OpenCascade | vendored | box, cylinder, trims-seams, hole, pocket | pass | pass | — |
| `occt-cone-frustum` | OpenCascade | vendored | cone, revolution | unsupported | unsupported | CAD-290 |
| `occt-sphere` | OpenCascade | vendored | sphere, revolution, trims-seams | invalid | fail | CAD-290 |
| `occt-revolved-vase` | OpenCascade | vendored | revolution, torus, cylinder | unsupported | unsupported | CAD-290 |
| `occt-spline-prism` | OpenCascade | vendored | spline, trims-seams | limit | fail | CAD-290 |
| `occt-two-bodies` | OpenCascade | vendored | multi-body, placement, assembly, box, cylinder | unsupported | unsupported | CAD-290 |
| `kicad-plate-holes` | KiCad | vendored | board, box, cylinder, hole, assembly | unsupported | unsupported | CAD-290 |
| `kicad-round-board` | KiCad | vendored | board, cylinder, hole, assembly | unsupported | unsupported | CAD-290 |
| `kicad-slotted-rounded-board` | KiCad | vendored | board, cylinder, trims-seams, pocket, assembly | unsupported | unsupported | CAD-290 |
| `kicad-board-with-part` | KiCad | vendored | board, multi-body, placement, assembly, cylinder | unsupported | unsupported | CAD-290 |
| `occt-annulus` | OpenCascade | repository | cylinder, hole, placement, trims-seams | pass | pass | — |
| `occt-regular-torus` | OpenCascade | repository | torus, revolution, placement | pass | pass | — |
| `kicad-vendored-Mounting_Wuerth_WA-SMST-3.3mm_H1mm_9774010960` | KiCad | vendored | cylinder, hole | pass | pass | section: CAD-346 |
| `kicad-vendored-L_Wuerth_HCM-7070` | KiCad | vendored | box, planar | pass | pass | — |
| `kicad-lib-Pin_D0.7mm_L6.5mm_W1.8mm_FlatFork` | FreeCAD | kicad-library | cylinder, trims-seams | pass | pass | section: CAD-346 |
| `kicad-lib-Heatsink_125x35x50mm_3xFixationM3` | FreeCAD | kicad-library | cylinder, hole, planar | pass | pass | section: CAD-346 |
| `kicad-lib-D_SOD-882` | KiCad | kicad-library | box, planar | pass | pass | — |
| `kicad-lib-LED_Rectangular_W3.0mm_H2.0mm` | FreeCAD | kicad-library | planar | pass | pass | — |
| `kicad-vendored-Mounting_Wuerth_WA-SMSE-ExternalM3_H5mm_9771050360` | KiCad | vendored | cone, revolution, torus | limit | fail | CAD-290 |
| `kicad-lib-RV_Disc_D9mm_W6.1mm_P5mm` | KiCad | kicad-library | revolution, torus | limit | fail | CAD-290 |
| `kicad-lib-Converter_ACDC_Hi-Link_HLK-20Mxx` | KiCad | kicad-library | sphere, torus | limit | fail | CAD-290 |
| `kicad-vendored-L_CommonModeChoke_TDK_ACM7060` | KiCad | vendored | planar | invalid | fail | CAD-290 |
| `kicad-lib-L_CommonModeChoke_Coilank_ACM2012` | KiCad | kicad-library | cylinder | invalid | fail | CAD-290 |
| `kicad-lib-SW_SPST_B3U-3000P-B` | KiCad | kicad-library | cone | invalid | fail | CAD-290 |
| `kicad-lib-LED_D3.0mm_Horizontal_O1.27mm_Z6.0mm` | FreeCAD | kicad-library | sphere | invalid | fail | CAD-290 |
| `kicad-lib-Buzzer_Mallory_AST1109MLTRQ` | KiCad | kicad-library | cylinder | invalid | fail | CAD-290 |
| `kicad-lib-Samtec_HPM-01-05-x-S_Straight_1x01_Pitch5.08mm` | KiCad | kicad-library | assembly, multi-body | unsupported | unsupported | CAD-290 |
| `kicad-lib-Mounting_Wuerth_WA-SMSI-M3_H1mm_9774010360` | KiCad | kicad-library | cone | unsupported | unsupported | CAD-290 |
| `kicad-lib-R_Box_L13.0mm_W4.0mm_P9.00mm` | KiCad | kicad-library | cylinder | unsupported | unsupported | CAD-290 |
| `kicad-lib-PinSocket_1x01_P2.54mm_Vertical` | KiCad | kicad-library | planar | unsupported | unsupported | CAD-290 |
| `kicad-lib-C_Axial_L5.1mm_D3.1mm_P10.00mm_Horizontal` | KiCad | kicad-library | revolution | unsupported | unsupported | CAD-290 |
| `kicad-vendored-WLP-4_0.83x0.83mm_P0.4mm` | KiCad | vendored | sphere | unsupported | unsupported | CAD-290 |
| `kicad-lib-D_SOD-923` | KiCad | kicad-library | spline | unsupported | unsupported | CAD-290 |
| `kicad-lib-R_LDR_5.2x5.2mm_P3.5mm_Horizontal` | FreeCAD | kicad-library | multi-body, surface-model, torus | unsupported | unsupported | CAD-290 |
| `kicad-lib-Converter_ACDC_MeanWell_IRM-60-xx_THT` | KiCad | kicad-library | assembly, multi-body, sphere, torus | unsupported | unsupported | CAD-290 |
| `kicad-lib-C_Radial_D6.3mm_H7.0mm_P2.50mm` | KiCad | kicad-library | revolution, spline, torus, voids | unsupported | unsupported | CAD-290 |
| `kicad-lib-L_Radial_D28.0mm_P29.20mm` | FreeCAD | kicad-library | multi-body, torus | unsupported | unsupported | CAD-290 |

## Admitted measurements

The deltas below come from the OpenCascade reference run.

- "Re-export Δ" compares the Rust readback of the re-export with the source import.
- "Mesh deviation" compares the source mesh with the readback mesh. Both use the same Rust tessellator, so the value is effectively zero.
- "OpenCascade mesh deviation" compares OpenCascade's own tessellation of the re-export with the Rust readback mesh.

| Admitted fixture | Volume mm³ | Faces | Mesh triangles | Section length mm | Re-export volume Δ mm³ | Bounds Δ mm | Mesh deviation mm | OpenCascade source Δ vol | OpenCascade re-export Δ vol | OpenCascade mesh deviation mm |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `freecad-drilled-box` | 4699.47 | 7 | 100 | refused (CAD-346) | 0 | 0 | 0 | 9.09495e-13 | 9.09495e-13 | 3.97205e-15 |
| `occt-placed-box` | 3456 | 6 | 12 | 80 | 0 | 1.77636e-15 | 1.77636e-15 | 5.09317e-10 | 2.72848e-12 | 1.77636e-15 |
| `occt-cylinder` | 1696.46 | 3 | 96 | 37.6991 | 0 | 0 | 0 | 3.41061e-12 | 3.63798e-12 | 0.0473118 |
| `occt-torus` | 2131.83 | 1 | 22464 | 150.796 | 0 | 0 | 0 | 5.00222e-12 | 5.00222e-12 | 0.032519 |
| `occt-box-void` | 5760 | 12 | 120 | 132 | 0 | 0 | 0 | 9.09495e-13 | 9.09495e-13 | 3.7943e-15 |
| `occt-pocketed-plate` | 9170.63 | 13 | 114 | 158.85 | 0 | 0 | 0 | 1.81899e-11 | 0 | 0.0335075 |
| `occt-annulus` | 117.626 | 4 | 200 | 18.5 | 4.26326e-14 | 8.88178e-16 | 1.98603e-15 | 3.97023e-10 | 3.12639e-13 | 0.0482975 |
| `occt-regular-torus` | 222.066 | 1 | 10816 | 25.9698 | 0 | 0 | 0 | 4.54747e-13 | 4.54747e-13 | 0.0215635 |
| `kicad-vendored-Mounting_Wuerth_WA-SMST-3.3mm_H1mm_9774010960` | 27.1434 | 6 | 252 | refused (CAD-346) | 0 | 0 | 0 | 6.39488e-14 | 6.39488e-14 | 0.0335075 |
| `kicad-vendored-L_Wuerth_HCM-7070` | 343.084 | 14 | 100 | 28 | 0 | 0 | 8.88178e-16 | 1.7053e-13 | 1.13687e-13 | 4.44089e-16 |
| `kicad-lib-Pin_D0.7mm_L6.5mm_W1.8mm_FlatFork` | 3.33226 | 15 | 212 | refused (CAD-346) | 0 | 0 | 0 | 1.33227e-15 | 4.44089e-16 | 0.00380829 |
| `kicad-lib-Heatsink_125x35x50mm_3xFixationM3` | 55195.2 | 73 | 1412 | refused (CAD-346) | 0 | 0 | 0 | 5.3551e-09 | 2.91038e-11 | 0.0100539 |
| `kicad-lib-D_SOD-882` | 0.2787 | 21 | 236 | 3.2 | 0 | 0 | 0 | 5.55112e-17 | 5.55112e-17 | 1.38778e-16 |
| `kicad-lib-LED_Rectangular_W3.0mm_H2.0mm` | 33.76 | 20 | 172 | 3.2 | 0 | 0 | 0 | 7.10543e-15 | 7.10543e-15 | 6.32925e-16 |

## Findings and owners

Every refusal or incorrect result is attributed in the manifest. Grouped:

| Finding | Fixtures | Owner |
| --- | --- | --- |
| **Incorrect: valid solids reported as invalid geometry.** OpenCascade validates them with BRepCheck. The cases are a complete OpenCascade sphere (the matrix lists complete spheres) and five KiCad models. In the sweep, the same false `Invalid STEP geometry` affects 298 library files. | `occt-sphere`, `kicad-vendored-L_CommonModeChoke_TDK_ACM7060`, `kicad-lib-L_CommonModeChoke_Coilank_ACM2012`, `kicad-lib-SW_SPST_B3U-3000P-B`, `kicad-lib-LED_D3.0mm_Horizontal_O1.27mm_Z6.0mm` | CAD-290 |
| **Incorrect: internal reconstruction mismatch.** "Rebuilt STEP extrusion differs from the imported analytic topology" on a model that OpenCascade validates. | `kicad-lib-Buzzer_Mallory_AST1109MLTRQ` | CAD-290 |
| **Incorrect: small inputs classified as resource limits.** A 9.6 kB spline prism and small torus and cone models report "STEP exceeds the Rust import budget". Every one of the 965 `Limit` results in the sweep contains `TOROIDAL_SURFACE`. | `occt-spline-prism`, `kicad-vendored-Mounting_Wuerth_WA-SMSE-ExternalM3_H5mm_9771050360`, `kicad-lib-RV_Disc_D9mm_W6.1mm_P5mm`, `kicad-lib-Converter_ACDC_Hi-Link_HLK-20Mxx` | CAD-290 |
| **Unsupported: STEP assembly product structure** (`NEXT_ASSEMBLY_USAGE_OCCURRENCE`, `CONTEXT_DEPENDENT_SHAPE_REPRESENTATION`, `ITEM_DEFINED_TRANSFORMATION`). Every KiCad board export, even a bare board, and OpenCascade's multi-body compound export are assemblies. | `occt-two-bodies`, all four `kicad-*board*` exports, `kicad-lib-Samtec_HPM-01-05-x-S_Straight_1x01_Pitch5.08mm` | CAD-290 |
| **Unsupported: general analytic families**: standalone cones, general revolutions, partial cylinders, non-orthogonal polyhedra, spherical caps, B-spline surfaces, surface models, toroidal blends, general `BREP_WITH_VOIDS` and other flat multi-solid files | `occt-cone-frustum`, `occt-revolved-vase` and 10 KiCad library entries | CAD-290 |
| **Section refused.** The Rust drawing projection refuses a horizontal section through imported solids with cylindrical faces: "requires straight edges and planar faces". OpenCascade blind-bore and pocket imports section correctly. | `freecad-drilled-box`, `kicad-vendored-Mounting_Wuerth_WA-SMST-3.3mm_H1mm_9774010960`, `kicad-lib-Pin_D0.7mm_L6.5mm_W1.8mm_FlatFork`, `kicad-lib-Heatsink_125x35x50mm_3xFixationM3` | CAD-346 |
| **Misclassification: malformed reference cycles** (an `EDGE_LOOP` or `VERTEX_POINT` that refers to itself) are refused as "unsupported analytic solid" rather than invalid. They are still refused atomically. | robustness `self-cycle.step`, `vertex-cycle.step` | CAD-290 |
| **Observation: `FILE_SCHEMA` is not checked.** A valid AP214 body under `CONFIG_CONTROL_DESIGN` or a nonsense schema name imports. Geometry, not header text, decides. This is pinned so that any change is deliberate. | robustness `file_schema_text_is_not_geometry_evidence` | CAD-290 (decision) |
| **Observation: STEP parsing has no interrupt.** Cancellation is checked before the read and again before commit, so a cancelled request never commits. A long parse cannot be interrupted midway. | robustness `cancelled_import_commits_nothing` | CAD-290 |

Oblique cylinder cuts are outside this corpus. During this work, the existing `multicad-mcp` test `rust_mcp_round_trips_a_step_body_and_refuses_an_unsupported_shape` failed at its fixture construction: "Oblique cylindrical subtraction requires a qualified curved Boolean provider". That failure is in the baseline, not caused by this branch, and is owned by CAD-354.

## Robustness ledger

Every case below is refused. Each refusal leaves the refused file, the valid source file and the previous valid document unchanged, for both a new import and an update of the existing import (`exchange_robustness.rs`).

| Input | Rust reader | Document/MCP path |
| --- | --- | --- |
| Empty; truncated at ¼, ½ and before the trailer; missing `DATA` | invalid | refused |
| Dangling reference; duplicate entity id; `NaN`; `1.E400` overflow | invalid | refused |
| Self-referencing `EDGE_LOOP` or `VERTEX_POINT` | unsupported (see findings) | refused |
| `1.E300` coordinate | limit | refused |
| Centimetre, unprefixed metre and conversion-based inch units | unsupported | refused |
| Binary garbage; non-UTF-8 byte | rejected before parsing | refused, "must contain UTF-8 text" |
| 40-level nesting; 70 KiB token | limit | refused |
| 1,048,577 entities (13.6 MB); 8,388,609 values (25.2 MB) | limit, after 5.6 s and 7.2 s | refused, "exceeds the 8 MiB import limit" |
| Exactly 8 MiB (padded valid file) | pass | **imported** |
| 8 MiB + 1 byte | — | refused, "exceeds the 8 MiB import limit" |
| Exactly 64 MiB | pass | refused by the 8 MiB document limit |
| 64 MiB + 1 byte | limit, also through `Solid::read_step` | refused |
| Cancelled request, new or update | — | `cancelled`; document unchanged |

Bounded healing: a vertex moved by 0.1 × the 1e-7 mm linear tolerance is accepted. The following are refused as invalid and never snapped, flipped or closed:

- a vertex moved by 10 × the tolerance or by 1 µm;
- a flipped face;
- a flipped edge use;
- a closed shell with a missing face.

Damaged projects were tested with the real executable (`--command --save`). The project was truncated at 8 bytes, 52 bytes, half its length and one byte short; bits were flipped in the magic, manifest, source and last byte; and one trailing byte was added. Each damaged project refuses to open and stays byte-identical. The last valid project still rebuilds its import without the original STEP file and produces the same parts.

Recovery copies were tested with FreeCAD, OpenCascade and KiCad library imports (`files/exchange_recovery_tests.rs`):

- Restoring a recovery copy rebuilds the same three solids without the sources, with the embedded digests equal to the foreign files.
- A recovery copy with a damaged payload is refused on restore.
- A recovery copy with a damaged manifest is excluded from the listing with a diagnostic.
- A stale `Replace` save and a non-project extension are refused, and the newer project bytes survive.
- Changing a pinned hatch definition reports `DefinitionChanged`. The recovery copy keeps the pinned original.

## KiCad library sweep

The sweep covered 7,236 STEP models in KiCad 10.0.0. Every file was classified; there were no panics or unclassified errors.

| Producer (`FILE_NAME`) | Files | Pass | Unsupported | Invalid | Limit |
| --- | ---: | ---: | ---: | ---: | ---: |
| KiCad | 5181 | 129 | 4480 | 53 | 519 |
| FreeCAD | 1031 | 44 | 654 | 51 | 282 |
| unstated | 796 | 0 | 453 | 179 | 164 |
| OpenCascade | 211 | 0 | 211 | 0 | 0 |
| Solid Edge / ST-DEVELOPER v16.7 | 14 | 0 | 0 | 14 | 0 |
| Autodesk Translation Framework v7.6.0.251 / ST-DEVELOPER v17.2 | 2 | 0 | 0 | 2 | 0 |
| Autodesk Translation Framework v13.14.0.145 / ST-DEVELOPER v20 | 1 | 0 | 1 | 0 | 0 |
| Total | 7236 | 173 | 5799 | 299 | 965 |

| Category (any surface or structure present) | Files | Pass | Unsupported | Invalid | Limit |
| --- | ---: | ---: | ---: | ---: | ---: |
| cylinder | 2902 | 127 | 2575 | 200 | 0 |
| torus | 1317 | 0 | 339 | 13 | 965 |
| spline | 1006 | 0 | 970 | 16 | 20 |
| sphere | 956 | 0 | 841 | 24 | 91 |
| planar | 866 | 46 | 769 | 51 | 0 |
| cone | 555 | 0 | 403 | 19 | 133 |
| multi-body | 537 | 0 | 522 | 15 | 0 |
| revolution | 532 | 0 | 400 | 0 | 132 |
| assembly | 382 | 0 | 368 | 14 | 0 |
| voids | 211 | 0 | 210 | 1 | 0 |
| surface-model | 2 | 0 | 2 | 0 | 0 |

The categories overlap: a file is counted once for each tag it carries. `planar` and `cylinder` are assigned only when no other tag applies. All 965 `limit` results carry the `torus` tag. Of the 299 `invalid` results, 298 report `Invalid STEP geometry` and one reports the extrusion reconstruction mismatch. The slowest single import took 4.5 s. The run took 29 min 34 s wall-clock and 11 min 23 s CPU in a debug build, under parallel load, at commit `4a243bfc`; that commit's import code is identical to `dff6a0f8`.

## Runtime contract

`cargo build --release --locked -p multicad-app` built the release executable in 13 min 33 s on the shared host. The production sources were those of `dff6a0f8`: CAD-345 changes only tests, fixtures, documentation and developer tooling. The executable is 52,470,160 bytes, with SHA-256 `51f8b2eaa950c416d0a9f1c93b4232232d84493ea4439d34f60debfa1bae071d`.

`node tools/verify-rust-runtime.mjs target/cargo/release/multicad` passes:

- `--diagnostics` reports `CAD kernel: MultiCAD Rust` and no OCCT source.
- `otool -L` lists only `/System/Library` frameworks, `libSystem`, `libc++`, `libobjc` and `libiconv`. No `TK*`, OCCT or OpenCascade library is linked.
- None of the OpenCascade runtime type names checked by the tool appear in the binary: `Standard_Transient`, `TopoDS_TShape`, `BRepCheck_Analyzer`, `STEPControl_Reader`, `BRepMesh_IncrementalMesh` and `Geom_BSplineSurface`.

As a control, the same tool refuses the `--features occt` debug executable. That executable identifies `OCCT 7.9.3 (statically linked)` and contains all six type names. The default debug executable passes.

OpenCascade appears only in test code (`#[cfg(feature = "occt")]`) and in the optional `target/occt` build. The unit tests for the tool run with `npm test`. In this worktree, `npm test` reports 3 of 5 files passing. The two release-contract files cannot load the `semver` devDependency because `node_modules` is not installed, and no installation was performed. `node --test tools/verify-rust-runtime.test.mjs` passes 3/3.

## Limitations

- These are debug-build test timings on a host shared with two other compiling agents. They are not performance claims.
- Mesh deviation between source and readback uses the same Rust tessellator, so it mainly proves identical geometry. The independent deviation comes from the OpenCascade tessellation in the reference run.
- The section check is one horizontal plane at mid-height. Analytic perimeters are pinned for five fixtures. KiCad section lengths are compared only between source and readback.
- Material persistence covers material, hatch, mass and pinned-definition diagnostics. Appearance and stock are covered elsewhere in the material tests.
- No FreeCAD, cadquery or OCP installation was available. FreeCAD-produced coverage comes from the existing fixture and from FreeCAD-authored KiCad library models.
- KiCad boards are always exported as assemblies, so they can pass only after the assembly finding is resolved.
