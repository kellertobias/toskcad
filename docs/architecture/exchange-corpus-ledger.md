# Exchange corpus ledger

This is the pinned pass/fail/unsupported ledger for the independent cross-CAD STEP corpus ([CAD-345](https://plans.tokenet.de/CAD/345)). The [manifest](../../crates/multicad-app/tests/fixtures/exchange/corpus/manifest.json) pins every fixture's producer, version, `FILE_NAME` strings, units, provenance, redistribution basis, SHA-256, geometry tags, expected Rust outcome, verdict and owner. Licensing and attribution are in the [corpus README](../../crates/multicad-app/tests/fixtures/exchange/corpus/README.md).

The machine-readable [ledger JSON](exchange-corpus-ledger.json) is the emitted output of the OpenCascade reference run. It records the measured source commit, tolerances, and every per-fixture measurement and timing. `corpus_manifest_pins_provenance_strata_and_owners` fails when the committed JSON and the manifest disagree on any fixture or outcome.

**Measured commit.** These results were measured at `1f1c4fe7` on the CAD-359 branch, which repairs the incorrect refusals that CAD-345 recorded at `dff6a0f8` (then 14 pass, 18 unsupported and 10 fail of 42). Rerun the commands below after every change to the reader, importer, drawing or Boolean code, and update the manifest outcomes, verdicts and this ledger together.

## Result

| | Count |
| --- | ---: |
| Fixtures | 46 |
| Vendored, always run | 20 (15 authored, 5 KiCad library copies) |
| Repository fixtures, always run | 2 |
| Referenced from the installed KiCad 10.0.0 library by path and SHA-256 | 24 |
| **Pass** (full workflow) | **22** |
| **Unsupported**: an honest refusal of a family outside the contract | **24** |
| **Fail**: an incorrect refusal of a file that OpenCascade validates | **0** |
| Unavailable on this host | 0 |
| Pass fixtures whose section is refused by drawing (CAD-346) | 6 |

A **pass** means that every check below succeeded. **Unsupported** means that the Rust importer refused a geometry family the [qualification matrix](rust-kernel-qualification.md) does not claim. **Fail** means that OpenCascade 7.9.3 reads the same source as a valid solid, so the Rust classification is wrong. No fixture is in that state at the measured commit.

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
| Analytic oracle: volume, bounds, face count, solid count and, where pinned, boundary area and planar face count | Volume and area within relative 1e-9; bounds within 1e-6 mm; counts exact |
| Independent Rust Part 21 readback of the source against the reconstructed solid | Volume within relative 1e-9; bounds within 1e-6 mm; members equal solids |
| Certified manufacturing mesh: each directed edge occurs once with its reverse, and the enclosed volume matches | Mesh volume within relative 2e-2; 0.05 mm linear deflection |
| Planar section at mid-height: exact line and arc lengths against the analytic perimeter where pinned | Relative 1e-9 |
| Material and section hatch: material, hatch, pinned-definition diagnostics and mass, after reopening | Exact; at least one hatch vector when the section succeeds |
| Canonical import: save the native project, delete the source, reopen, rebuild with Rust and compare the embedded source digest | Volume and bounds bit-identical; digest equal |
| Re-export, independent readback and reimport | Volume within relative 1e-9; bounds within 1e-6 mm; faces and solids exact |
| Geometric deviation of the readback mesh from the source mesh | ≤ 0.05 mm, searched within 0.2 mm |
| OpenCascade 7.9.3 source readback (`--features occt`) | Volume within relative 1e-6; bounds within 1e-3 mm |
| OpenCascade 7.9.3 source supports (`--features occt`): total boundary area and planar boundary area, which do not depend on how a surface is divided into faces | Relative 1e-6 |
| Deviation between OpenCascade's tessellation of the source and the Rust mesh of the import (`--features occt`) | ≤ 0.1 mm, twice the deflection |
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
| `occt-sphere` | OpenCascade | vendored | sphere, revolution, trims-seams | pass | pass | — |
| `occt-revolved-vase` | OpenCascade | vendored | revolution, torus, cylinder | unsupported | unsupported | CAD-290 |
| `occt-spline-prism` | OpenCascade | vendored | spline, trims-seams | pass | pass | section: CAD-346 |
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
| `kicad-vendored-Mounting_Wuerth_WA-SMSE-ExternalM3_H5mm_9771050360` | KiCad | vendored | cone, revolution, torus | unsupported | unsupported | CAD-360 |
| `kicad-lib-RV_Disc_D9mm_W6.1mm_P5mm` | KiCad | kicad-library | revolution, torus | unsupported | unsupported | CAD-360 |
| `kicad-lib-Converter_ACDC_Hi-Link_HLK-20Mxx` | KiCad | kicad-library | sphere, torus | unsupported | unsupported | CAD-360 |
| `kicad-vendored-L_CommonModeChoke_TDK_ACM7060` | KiCad | vendored | planar | pass | pass | — |
| `kicad-lib-L_CommonModeChoke_Coilank_ACM2012` | KiCad | kicad-library | cylinder | unsupported | unsupported | CAD-360 |
| `kicad-lib-SW_SPST_B3U-3000P-B` | KiCad | kicad-library | cone | unsupported | unsupported | CAD-360 |
| `kicad-lib-LED_D3.0mm_Horizontal_O1.27mm_Z6.0mm` | FreeCAD | kicad-library | sphere | unsupported | unsupported | CAD-360 |
| `kicad-lib-Buzzer_Mallory_AST1109MLTRQ` | KiCad | kicad-library | box, cylinder, hole, pocket | pass | pass | — |
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
| `kicad-lib-R_0603_1608Metric` | KiCad | kicad-library | cylinder, planar | pass | pass | section: CAD-346 |
| `kicad-lib-QFN-32-1EP_4x4mm_P0.4mm_EP2.9x2.9mm` | KiCad | kicad-library | box, planar | pass | pass | — |
| `kicad-lib-DFN-8-1EP_3x2mm_P0.5mm_EP1.3x1.5mm` | KiCad | kicad-library | box, planar | pass | pass | — |
| `kicad-lib-L_CommonModeChoke_Bourns_SRF1260` | KiCad | kicad-library | box, planar | pass | pass | — |

## Admitted measurements

The deltas below come from the OpenCascade reference run.

- "Re-export Δ" compares the Rust readback of the re-export with the source import.
- "Mesh deviation" compares the source mesh with the readback mesh. Both use the same Rust tessellator, so the value is effectively zero.
- "OpenCascade source mesh" compares OpenCascade's tessellation of the transmitted source with the Rust mesh of the import.
- "OpenCascade mesh deviation" compares OpenCascade's own tessellation of the re-export with the Rust readback mesh.
- "Faces" lists the Rust faces and, where OpenCascade keeps a finer division of the same surfaces, its count. The boundary area and the planar part of it agree with OpenCascade within relative 1e-6 for every row.

| Admitted fixture | Volume mm³ | Faces (OpenCascade) | Mesh triangles | Section length mm | Re-export volume Δ mm³ | Bounds Δ mm | Mesh deviation mm | OpenCascade source Δ vol | OpenCascade area Δ mm² | OpenCascade source mesh mm | OpenCascade re-export Δ vol | OpenCascade mesh deviation mm |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `freecad-drilled-box` | 4699.47 | 7 | 100 | refused (CAD-346) | 0 | 0 | 0 | 9.09495e-13 | 0 | 5.32907e-15 | 9.09495e-13 | 3.97205e-15 |
| `occt-placed-box` | 3456 | 6 | 12 | 80 | 0 | 1.77636e-15 | 1.77636e-15 | 5.09317e-10 | 9.28e-11 | 5.75523e-13 | 2.72848e-12 | 1.77636e-15 |
| `occt-cylinder` | 1696.46 | 3 | 96 | 37.6991 | 0 | 0 | 0 | 3.41061e-12 | 2.27e-13 | 0.0473118 | 3.63798e-12 | 0.0473118 |
| `occt-torus` | 2131.83 | 1 | 22464 | 150.796 | 0 | 0 | 0 | 5.00222e-12 | 2.27e-13 | 0.032519 | 5.00222e-12 | 0.032519 |
| `occt-box-void` | 5760 | 12 | 120 | 132 | 0 | 0 | 0 | 9.09495e-13 | 0 | 3.66205e-15 | 9.09495e-13 | 3.7943e-15 |
| `occt-pocketed-plate` | 9170.63 | 13 | 114 | 158.85 | 0 | 0 | 0 | 1.81899e-11 | 0 | 0.0335075 | 0 | 0.0335075 |
| `occt-sphere` | 523.599 | 1 | 960 | 31.4159 | 0 | 0 | 0 | 9.09495e-13 | 5.68e-14 | 0.0240217 | 9.09495e-13 | 0.0240217 |
| `occt-spline-prism` | 864 | 4 | 112 | refused (CAD-346) | 0 | 0 | 0 | 1.59162e-12 | 3.31e-05 | 0.036732 | 1.81899e-12 | 0.036732 |
| `occt-annulus` | 117.626 | 4 | 200 | 18.5 | 4.26326e-14 | 8.88178e-16 | 1.98603e-15 | 3.97023e-10 | 3.8e-10 | 0.04827 | 3.12639e-13 | 0.0482975 |
| `occt-regular-torus` | 222.066 | 1 | 10816 | 25.9698 | 0 | 0 | 0 | 4.54747e-13 | 2.84e-13 | 0.0215635 | 4.54747e-13 | 0.0215635 |
| `kicad-vendored-Mounting_Wuerth_WA-SMST-3.3mm_H1mm_9774010960` | 27.1434 | 6 | 252 | refused (CAD-346) | 0 | 0 | 0 | 6.39488e-14 | 1.42e-14 | 0.0335075 | 6.39488e-14 | 0.0335075 |
| `kicad-vendored-L_Wuerth_HCM-7070` | 343.084 | 14 (16) | 100 | 28 | 0 | 0 | 8.88178e-16 | 1.7053e-13 | 1.14e-13 | 8.88178e-16 | 1.13687e-13 | 4.44089e-16 |
| `kicad-lib-Pin_D0.7mm_L6.5mm_W1.8mm_FlatFork` | 3.33226 | 15 | 212 | refused (CAD-346) | 0 | 0 | 0 | 1.33227e-15 | 7.11e-15 | 0.00380829 | 4.44089e-16 | 0.00380829 |
| `kicad-lib-Heatsink_125x35x50mm_3xFixationM3` | 55195.2 | 73 | 1412 | refused (CAD-346) | 0 | 0 | 0 | 5.3551e-09 | 4.22e-10 | 0.0100539 | 2.91038e-11 | 0.0100539 |
| `kicad-lib-D_SOD-882` | 0.2787 | 21 | 236 | 3.2 | 0 | 0 | 0 | 5.55112e-17 | 8.88e-16 | 1.38778e-16 | 5.55112e-17 | 1.38778e-16 |
| `kicad-lib-LED_Rectangular_W3.0mm_H2.0mm` | 33.76 | 20 | 172 | 3.2 | 0 | 0 | 0 | 7.10543e-15 | 2.84e-14 | 5.23691e-16 | 7.10543e-15 | 6.32925e-16 |
| `kicad-vendored-L_CommonModeChoke_TDK_ACM7060` | 147 | 6 (15) | 12 | 26 | 0 | 0 | 0 | 5.68434e-14 | 5.68e-14 | 9.93014e-16 | 2.84217e-14 | 0 |
| `kicad-lib-Buzzer_Mallory_AST1109MLTRQ` | 168.929 | 24 (26) | 564 | 40 | 0 | 0 | 0 | 2.84217e-14 | 1.14e-13 | 0.00390921 | 0 | 0.00390921 |
| `kicad-lib-R_0603_1608Metric` | 0.537219 | 20 (26) | 228 | refused (CAD-346) | 0 | 0 | 0 | 2.22045e-16 | 0 | 0.000381795 | 2.22045e-16 | 0.000381795 |
| `kicad-lib-QFN-32-1EP_4x4mm_P0.4mm_EP2.9x2.9mm` | 9.6 | 6 (72) | 12 | 16 | 0 | 0 | 0 | 1.77636e-15 | 2.84e-14 | 3.25321e-15 | 0 | 0 |
| `kicad-lib-DFN-8-1EP_3x2mm_P0.5mm_EP1.3x1.5mm` | 3.6 | 6 (24) | 12 | 10 | 0 | 0 | 0 | 4.44089e-16 | 3.55e-15 | 8.67112e-16 | 4.44089e-16 | 0 |
| `kicad-lib-L_CommonModeChoke_Bourns_SRF1260` | 968.75 | 6 (15) | 12 | 50 | 0 | 0 | 0 | 0 | 2.27e-13 | 2.40178e-15 | 0 | 0 |

## Findings and owners

Every refusal is attributed in the manifest. The incorrect results that CAD-345 recorded were repaired or reclassified by [CAD-359](https://plans.tokenet.de/CAD/359):

| CAD-345 finding | Root cause | Result |
| --- | --- | --- |
| Valid solids reported as invalid geometry (298 library files). | Exchange writers divide one plane into several faces for pads, colours and markers. A circle or arc then joins two faces of one plane, which the family readers treated as malformed. Four narrower causes produced the rest: the verdict depended on which face was listed first, a round edge between two faces of one cylinder, a ring the profile certificate cannot resolve, and a straight spline edge trimmed by its vertices. | Faces that partition one plane, and the straight edges they split, are merged before family recognition when the transmitted shell is refused. `kicad-vendored-L_CommonModeChoke_TDK_ACM7060` and 22 more library files import; the others are refused as unsupported. `occt-sphere` imports: a vertex loop has no direction, so OpenCascade's reversed pole bound is the same face. |
| Internal reconstruction mismatch on a valid model. | The rebuilt drilled orthogonal shell merges coplanar source faces, but only the plain, pocketed and slotted variants were allowed to. | `kicad-lib-Buzzer_Mallory_AST1109MLTRQ` passes the full workflow. |
| Small inputs classified as resource limits (965 library files). | A toroidal face on a shell of more than five faces, and a spline shell of fewer than five, were reported as limits. All 965 files carry a toroidal face on a shell of at least nine faces. | They are refused as unsupported. `occt-spline-prism` imports: its cubic and chord section is certified by a convex control polygon and its rims by their master parameter lines. |
| Malformed reference cycles refused as unsupported. | A typed attribute that referred back to its own entity was reported as an unsupported entity kind. | A reference cycle below a solid is invalid for every entity kind. |

The open findings are:

| Finding | Fixtures | Owner |
| --- | --- | --- |
| **Unsupported: STEP assembly product structure** (`NEXT_ASSEMBLY_USAGE_OCCURRENCE`, `CONTEXT_DEPENDENT_SHAPE_REPRESENTATION`, `ITEM_DEFINED_TRANSFORMATION`). Every KiCad board export, even a bare board, and OpenCascade's multi-body compound export are assemblies. | `occt-two-bodies`, all four `kicad-*board*` exports, `kicad-lib-Samtec_HPM-01-05-x-S_Straight_1x01_Pitch5.08mm` | CAD-290 |
| **Unsupported: general analytic families**: standalone cones, general revolutions, partial cylinders, non-orthogonal polyhedra, spherical caps, B-spline surfaces, surface models, toroidal blends, general `BREP_WITH_VOIDS` and other flat multi-solid files | `occt-cone-frustum`, `occt-revolved-vase` and 10 KiCad library entries | CAD-290 |
| **Unsupported: general imported B-rep.** Valid shells that no analytic family describes: a revolved standoff with cones and a toroidal fillet, a varistor with `SURFACE_OF_REVOLUTION` faces, pins with spherical and toroidal blends, a hemispherical lens, and a box with notched corner pads. | `kicad-vendored-Mounting_Wuerth_WA-SMSE-ExternalM3_H5mm_9771050360`, `kicad-lib-RV_Disc_D9mm_W6.1mm_P5mm`, `kicad-lib-Converter_ACDC_Hi-Link_HLK-20Mxx`, `kicad-lib-LED_D3.0mm_Horizontal_O1.27mm_Z6.0mm`, `kicad-lib-L_CommonModeChoke_Coilank_ACM2012` | CAD-360 |
| **Unsupported: pinched non-manifold contact.** Four edges are each used by four faces in balanced directions. ISO 10303-42 requires every edge of a closed shell to be referenced exactly twice; OpenCascade tolerates the contact. The Rust kernel needs a two-manifold shell. | `kicad-lib-SW_SPST_B3U-3000P-B` | CAD-360 |
| **Section refused.** The Rust drawing projection refuses a horizontal section through imported solids with cylindrical faces ("requires straight edges and planar faces") and through spline sweeps ("requires an uncropped axial projection"). | `freecad-drilled-box`, `kicad-vendored-Mounting_Wuerth_WA-SMST-3.3mm_H1mm_9774010960`, `kicad-lib-Pin_D0.7mm_L6.5mm_W1.8mm_FlatFork`, `kicad-lib-Heatsink_125x35x50mm_3xFixationM3`, `kicad-lib-R_0603_1608Metric`, `occt-spline-prism` | CAD-346 |
| **Observation: `FILE_SCHEMA` is not checked.** A valid AP214 body under `CONFIG_CONTROL_DESIGN` or a nonsense schema name imports. Geometry, not header text, decides. This is pinned so that any change is deliberate. | robustness `file_schema_text_is_not_geometry_evidence` | CAD-290 (decision) |
| **Observation: STEP parsing has no interrupt.** Cancellation is checked before the read and again before commit, so a cancelled request never commits. A long parse cannot be interrupted midway. | robustness `cancelled_import_commits_nothing` | CAD-290 |

Merged partitions change the face count, not the boundary. The import of a partitioned model has one face for each plane, while OpenCascade keeps every transmitted face. The reference run therefore compares volume, bounds, total boundary area, planar boundary area and both tessellations, and the analytic oracles pin the merged face counts.

Oblique cylinder cuts are outside this corpus. During the CAD-345 work, the existing `multicad-mcp` test `rust_mcp_round_trips_a_step_body_and_refuses_an_unsupported_shape` failed at its fixture construction: "Oblique cylindrical subtraction requires a qualified curved Boolean provider". That failure is in the baseline and is owned by CAD-354.

## Robustness ledger

Every case below is refused. Each refusal leaves the refused file, the valid source file and the previous valid document unchanged, for both a new import and an update of the existing import (`exchange_robustness.rs`).

| Input | Rust reader | Document/MCP path |
| --- | --- | --- |
| Empty; truncated at ¼, ½ and before the trailer; missing `DATA` | invalid | refused |
| Dangling reference; duplicate entity id; `NaN`; `1.E400` overflow | invalid | refused |
| Self-referencing `EDGE_LOOP` or `VERTEX_POINT` | invalid | refused |
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

The sweep covered 7236 STEP models in KiCad 10.0.0. Every file was classified; there were no panics or unclassified errors.

| Producer (`FILE_NAME`) | Files | Pass | Unsupported | Invalid | Limit |
| --- | ---: | ---: | ---: | ---: | ---: |
| KiCad | 5181 | 174 | 5007 | 0 | 0 |
| FreeCAD | 1031 | 45 | 976 | 10 | 0 |
| unstated | 796 | 0 | 796 | 0 | 0 |
| OpenCascade | 211 | 0 | 211 | 0 | 0 |
| Solid Edge / ST-DEVELOPER v16.7 | 14 | 0 | 14 | 0 | 0 |
| Autodesk Translation Framework v7.6.0.251 / ST-DEVELOPER v17.2 | 2 | 0 | 2 | 0 | 0 |
| Autodesk Translation Framework v13.14.0.145 / ST-DEVELOPER v20 | 1 | 0 | 1 | 0 | 0 |
| Total | 7236 | 219 | 7007 | 10 | 0 |

| Category (any surface or structure present) | Files | Pass | Unsupported | Invalid | Limit |
| --- | ---: | ---: | ---: | ---: | ---: |
| cylinder | 2902 | 150 | 2748 | 4 | 0 |
| torus | 1317 | 0 | 1317 | 0 | 0 |
| spline | 1006 | 0 | 1006 | 0 | 0 |
| sphere | 956 | 0 | 955 | 1 | 0 |
| planar | 866 | 69 | 796 | 1 | 0 |
| cone | 555 | 0 | 551 | 4 | 0 |
| multi-body | 537 | 0 | 537 | 0 | 0 |
| revolution | 532 | 0 | 532 | 0 | 0 |
| assembly | 382 | 0 | 382 | 0 | 0 |
| voids | 211 | 0 | 210 | 1 | 0 |
| surface-model | 2 | 0 | 2 | 0 | 0 |

The categories overlap: a file is counted once for each tag it carries. `planar` and `cylinder` are assigned only when no other tag applies. The slowest single import took 1.5 s in a debug build at commit `1f1c4fe7`.

At `dff6a0f8` the same sweep gave 173 pass, 5,799 unsupported, 299 invalid and 965 limit. At `1f1c4fe7` it gives 219 pass, 7007 unsupported, 10 invalid and 0 limit. The 46 newly passing files are 23 boxes with partitioned faces, 22 chip resistors and fuses with partitioned ends, and the buzzer. All 965 former limits and 289 former invalid results are now either imports or unsupported refusals.

The 10 remaining invalid results are refusals that the tolerance and topology contract requires:

| Files | Measured cause |
| --- | --- |
| `Heatsink_AAVID_Extruded_531002B02100G`, `…531102B02100G`, `…531202B02100G` | Eight vertices lie 1.88e-7 mm off their lines. The file declares 2e-5 mm; the import tolerance is 1e-7 mm. |
| `ESP32-WROOM-32E` | 140 vertices lie up to 2.15e-7 mm off their lines (declared 2e-7 mm). |
| `R_Shunt_Ohmite_LVK12` | 22 vertices lie up to 2.94e-7 mm off their lines (declared 2e-5 mm). |
| `Relay_DPDT_FRT5_SMD` | One vertex lies 1.95e-3 mm off its line (declared 5e-4 mm). |
| `HDMI_Micro-D_Molex_46765-0xxx`, `-1xxx`, `-2xxx` | 308 vertices lie up to 2e-6 mm off their planes (declared 1e-5 or 2e-5 mm). |
| `Diode_Bridge_GeneSiC_KBPC_T` | Nine edge loops are not head to tail: one oriented edge in each has the wrong orientation flag. |

Healing is bounded to the 1e-7 mm linear tolerance, so the first nine are refused by design; honouring a file's declared uncertainty would be a product decision. The last file violates the ISO 10303-42 path rule and OpenCascade repairs it silently.

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
