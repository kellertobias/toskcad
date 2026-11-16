# Exchange corpus ledger

This is the pinned pass/fail/unsupported ledger for the independent cross-CAD STEP corpus ([CAD-345](https://plans.tokenet.de/CAD/345)). The [manifest](../../crates/multicad-app/tests/fixtures/exchange/corpus/manifest.json) pins every fixture's producer, version, `FILE_NAME` strings, units, provenance, redistribution basis, SHA-256, geometry tags, expected Rust outcome, verdict and owner. Licensing and attribution are in the [corpus README](../../crates/multicad-app/tests/fixtures/exchange/corpus/README.md).

The machine-readable [ledger JSON](exchange-corpus-ledger.json) is the emitted output of the OpenCascade reference run. It records the measured source commit, tolerances, and every per-fixture measurement and timing. `corpus_manifest_pins_provenance_strata_and_owners` fails when the committed JSON and the manifest disagree on any fixture or outcome.

**Measured commit.** These results were measured at `3cc58b11` on the CAD-361 branch `claude/step-assemblies`, which adds STEP product structure on top of the CAD-360 retained-boundary provider (merged at `1869dd5b`). The tracked changes during the run were this manifest's CAD-361 verdicts, the per-solid OpenCascade section reference and the multi-body MCP controls, all committed with this ledger. CAD-360 recorded 36 pass at `c0fddf10`. CAD-345 recorded 14 pass, 18 unsupported and 10 fail of 42 at `dff6a0f8`, and CAD-359 22 pass and 24 unsupported of 46 at `1f1c4fe7`. Rerun the commands below after every change to the reader, importer, drawing or Boolean code, and update the manifest outcomes, verdicts and this ledger together.

**Measured commit.** These results were measured at `b333779e` on the CAD-360 branch with a clean tree; since `c0fddf10` it reads horn tori, orients trim loops by their exact area, keeps a stated vertex nearer to its supports than their balance, refuses sliver faces and files consistent only beyond the 1e-5 mm incidence cap as unsupported, and publishes exact line and ellipse section edges. The pinned outcomes are unchanged. CAD-345 recorded 14 pass, 18 unsupported and 10 fail of 42 at `dff6a0f8`, CAD-359 22 pass and 24 unsupported of 46 at `1f1c4fe7`, and CAD-360 the same counts as now at `c0fddf10`. Rerun the commands below after every change to the reader, importer, drawing or Boolean code, and update the manifest outcomes, verdicts and this ledger together.

## Result

| | Count |
| --- | ---: |
| Fixtures | 46 |
| Vendored, always run | 20 (15 authored, 5 KiCad library copies) |
| Repository fixtures, always run | 2 |
| Referenced from the installed KiCad 10.0.0 library by path and SHA-256 | 24 |
| **Pass** (full workflow) | **42** (14 through the retained-boundary provider, CAD-360; 6 assemblies through it, CAD-361) |
| **Unsupported** or **invalid**: an honest refusal outside the contract | **2 unsupported, 2 invalid** |
| **Fail**: an incorrect refusal of a file that OpenCascade validates | **0** |
| Unavailable on this host | 0 |
| Pass fixtures whose section is refused by drawing | 0 (sections through the exact export, CAD-360) |

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
| `freecad-drilled-box` | FreeCAD | vendored | box, cylinder, trims-seams, placement, hole | pass | pass | — |
| `occt-placed-box` | OpenCascade | vendored | box, placement | pass | pass | — |
| `occt-cylinder` | OpenCascade | vendored | cylinder, trims-seams, placement | pass | pass | — |
| `occt-torus` | OpenCascade | vendored | torus, revolution, trims-seams | pass | pass | — |
| `occt-box-void` | OpenCascade | vendored | box, voids | pass | pass | — |
| `occt-pocketed-plate` | OpenCascade | vendored | box, cylinder, trims-seams, hole, pocket | pass | pass | — |
| `occt-cone-frustum` | OpenCascade | vendored | cone, revolution | pass (retained boundary) | pass | — |
| `occt-sphere` | OpenCascade | vendored | sphere, revolution, trims-seams | pass | pass | — |
| `occt-revolved-vase` | OpenCascade | vendored | revolution, torus, cylinder | pass (retained boundary) | pass | — |
| `occt-spline-prism` | OpenCascade | vendored | spline, trims-seams | pass | pass | — |
| `occt-two-bodies` | OpenCascade | vendored | multi-body, placement, assembly, box, cylinder | pass (assembly) | pass | CAD-361 (pass) |
| `kicad-plate-holes` | KiCad | vendored | board, box, cylinder, hole, assembly | pass (assembly) | pass | CAD-361 (pass) |
| `kicad-round-board` | KiCad | vendored | board, cylinder, hole, assembly | pass (assembly) | pass | CAD-361 (pass) |
| `kicad-slotted-rounded-board` | KiCad | vendored | board, cylinder, trims-seams, pocket, assembly | pass (assembly) | pass | CAD-361 (pass) |
| `kicad-board-with-part` | KiCad | vendored | board, multi-body, placement, assembly, cylinder | pass (assembly) | pass | CAD-361 (pass) |
| `occt-annulus` | OpenCascade | repository | cylinder, hole, placement, trims-seams | pass | pass | — |
| `occt-regular-torus` | OpenCascade | repository | torus, revolution, placement | pass | pass | — |
| `kicad-vendored-Mounting_Wuerth_WA-SMST-3.3mm_H1mm_9774010960` | KiCad | vendored | cylinder, hole | pass | pass | — |
| `kicad-vendored-L_Wuerth_HCM-7070` | KiCad | vendored | box, planar | pass | pass | — |
| `kicad-lib-Pin_D0.7mm_L6.5mm_W1.8mm_FlatFork` | FreeCAD | kicad-library | cylinder, trims-seams | pass | pass | — |
| `kicad-lib-Heatsink_125x35x50mm_3xFixationM3` | FreeCAD | kicad-library | cylinder, hole, planar | pass | pass | — |
| `kicad-lib-D_SOD-882` | KiCad | kicad-library | box, planar | pass | pass | — |
| `kicad-lib-LED_Rectangular_W3.0mm_H2.0mm` | FreeCAD | kicad-library | planar | pass | pass | — |
| `kicad-vendored-Mounting_Wuerth_WA-SMSE-ExternalM3_H5mm_9771050360` | KiCad | vendored | cone, revolution, torus | pass (retained boundary) | pass | — |
| `kicad-lib-RV_Disc_D9mm_W6.1mm_P5mm` | KiCad | kicad-library | revolution, torus | pass (retained boundary) | pass | — |
| `kicad-lib-Converter_ACDC_Hi-Link_HLK-20Mxx` | KiCad | kicad-library | sphere, torus | pass (retained boundary) | pass | — |
| `kicad-vendored-L_CommonModeChoke_TDK_ACM7060` | KiCad | vendored | planar | pass | pass | — |
| `kicad-lib-L_CommonModeChoke_Coilank_ACM2012` | KiCad | kicad-library | cylinder | pass (retained boundary) | pass | — |
| `kicad-lib-SW_SPST_B3U-3000P-B` | KiCad | kicad-library | cone | invalid | unsupported | CAD-360 |
| `kicad-lib-LED_D3.0mm_Horizontal_O1.27mm_Z6.0mm` | FreeCAD | kicad-library | sphere | pass (retained boundary) | pass | — |
| `kicad-lib-Buzzer_Mallory_AST1109MLTRQ` | KiCad | kicad-library | box, cylinder, hole, pocket | pass | pass | — |
| `kicad-lib-Samtec_HPM-01-05-x-S_Straight_1x01_Pitch5.08mm` | KiCad | kicad-library | assembly, multi-body | pass (assembly) | pass | CAD-361 (pass) |
| `kicad-lib-Mounting_Wuerth_WA-SMSI-M3_H1mm_9774010360` | KiCad | kicad-library | cone | pass (retained boundary) | pass | — |
| `kicad-lib-R_Box_L13.0mm_W4.0mm_P9.00mm` | KiCad | kicad-library | cylinder | pass (retained boundary) | pass | — |
| `kicad-lib-PinSocket_1x01_P2.54mm_Vertical` | KiCad | kicad-library | planar | pass (retained boundary) | pass | — |
| `kicad-lib-C_Axial_L5.1mm_D3.1mm_P10.00mm_Horizontal` | KiCad | kicad-library | revolution | pass (retained boundary) | pass | — |
| `kicad-vendored-WLP-4_0.83x0.83mm_P0.4mm` | KiCad | vendored | sphere | pass (retained boundary) | pass | — |
| `kicad-lib-D_SOD-923` | KiCad | kicad-library | spline | pass (retained boundary) | pass | — |
| `kicad-lib-R_LDR_5.2x5.2mm_P3.5mm_Horizontal` | FreeCAD | kicad-library | multi-body, surface-model, torus | unsupported | unsupported | CAD-363 |
| `kicad-lib-Converter_ACDC_MeanWell_IRM-60-xx_THT` | KiCad | kicad-library | assembly, multi-body, sphere, torus | unsupported | unsupported | CAD-361 |
| `kicad-lib-C_Radial_D6.3mm_H7.0mm_P2.50mm` | KiCad | kicad-library | revolution, spline, torus, voids | invalid | unsupported | CAD-360 |
| `kicad-lib-L_Radial_D28.0mm_P29.20mm` | FreeCAD | kicad-library | multi-body, torus | pass (retained boundary) | pass | — |
| `kicad-lib-R_0603_1608Metric` | KiCad | kicad-library | cylinder, planar | pass | pass | — |
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

| Admitted fixture | Volume mm³ | Faces (OpenCascade) | Mesh triangles | Section length mm | Re-export volume Δ mm³ | Bounds Δ mm | Mesh deviation mm | OpenCascade source Δ vol | OpenCascade area Δ mm² | OpenCascade source mesh mm | OpenCascade re-export Δ vol | OpenCascade mesh deviation mm | OpenCascade section Δ rel |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `freecad-drilled-box` | 4699.47 | 7 | 100 | 112.566 | 0 | 0 | 0 | 9.09495e-13 | 0 | 4.44089e-15 | 9.09495e-13 | 5.0243e-15 | 0 |
| `occt-placed-box` | 3456 | 6 | 12 | 80 | 0 | 1.77636e-15 | 1.77636e-15 | 5.09317e-10 | 9.28e-11 | 5.75523e-13 | 2.72848e-12 | 1.77636e-15 | 1.03029e-14 |
| `occt-cylinder` | 1696.46 | 3 | 96 | 37.6991 | 0 | 0 | 0 | 3.41061e-12 | 2.27e-13 | 0.0473038 | 3.63798e-12 | 0.0473038 | 0 |
| `occt-torus` | 2131.83 | 1 | 22464 | 150.796 | 0 | 0 | 0 | 5.00222e-12 | 2.27e-13 | 0.00976132 | 5.00222e-12 | 0.0233943 | 0 |
| `occt-box-void` | 5760 | 12 | 120 | 132 | 0 | 0 | 0 | 9.09495e-13 | 0 | 3.66205e-15 | 9.09495e-13 | 3.7943e-15 | 0 |
| `occt-pocketed-plate` | 9170.63 | 13 | 114 | 158.85 | 0 | 0 | 0 | 1.81899e-11 | 0 | 0.0334964 | 0 | 0.0334964 | 0 |
| `occt-cone-frustum` | 1172.86 | 3 | 188 | 37.6991 | 5.01814e-10 | 0 | 0.0468963 | 1.00499e-10 | 4.68e-11 | 0.0385173 | 1.00727e-10 | 0.0468904 | 0 |
| `occt-sphere` | 523.599 | 1 | 960 | 31.4159 | 0 | 0 | 0 | 9.09495e-13 | 5.68e-14 | 0.0475311 | 9.09495e-13 | 0.0475311 | 0 |
| `occt-revolved-vase` | 3196.45 | 5 | 2510 | 59.464 | 0 | 0 | 0 | 6.82121e-12 | 1.14e-12 | 0.0480756 | 6.82121e-12 | 0.0480756 | 1.12199e-05 |
| `occt-spline-prism` | 864 | 4 | 112 | 50.1524 | 0 | 0 | 0 | 1.59162e-12 | 3.31e-05 | 0.0366404 | 1.81899e-12 | 0.0366404 | 6.29293e-05 |
| `occt-annulus` | 117.626 | 4 | 200 | 18.5 | 4.26326e-14 | 8.88178e-16 | 1.98603e-15 | 3.97023e-10 | 3.8e-10 | 0.0482895 | 3.12639e-13 | 0.0482975 | 8.57452e-13 |
| `occt-regular-torus` | 222.066 | 1 | 10816 | 25.9698 | 0 | 0 | 0 | 4.54747e-13 | 2.84e-13 | 0.00870197 | 4.54747e-13 | 0.0151882 | 0.000785494 |
| `kicad-vendored-Mounting_Wuerth_WA-SMST-3.3mm_H1mm_9774010960` | 27.1434 | 6 | 252 | 23.5619 | 0 | 0 | 0 | 6.39488e-14 | 1.42e-14 | 0.0334964 | 6.39488e-14 | 0.0334964 | 0 |
| `kicad-vendored-L_Wuerth_HCM-7070` | 343.084 | 14 (16) | 100 | 28 | 0 | 0 | 8.88178e-16 | 1.7053e-13 | 1.14e-13 | 8.88178e-16 | 1.13687e-13 | 4.44089e-16 | 0 |
| `kicad-lib-Pin_D0.7mm_L6.5mm_W1.8mm_FlatFork` | 3.33226 | 15 | 212 | 3.63525 | 0 | 0 | 0 | 1.33227e-15 | 7.11e-15 | 0.00380829 | 4.44089e-16 | 0.00380829 | 2.07675e-15 |
| `kicad-lib-Heatsink_125x35x50mm_3xFixationM3` | 55195.2 | 73 | 1412 | 915.718 | 0 | 0 | 0 | 5.3551e-09 | 4.22e-10 | 0.0100539 | 2.91038e-11 | 0.0100539 | 2.48301e-16 |
| `kicad-lib-D_SOD-882` | 0.2787 | 21 | 236 | 3.2 | 0 | 0 | 0 | 5.55112e-17 | 8.88e-16 | 1.38778e-16 | 5.55112e-17 | 1.38778e-16 | 2.77556e-16 |
| `kicad-lib-LED_Rectangular_W3.0mm_H2.0mm` | 33.76 | 20 | 172 | 3.2 | 0 | 0 | 0 | 7.10543e-15 | 2.84e-14 | 5.23691e-16 | 7.10543e-15 | 6.32925e-16 | 0 |
| `kicad-vendored-Mounting_Wuerth_WA-SMSE-ExternalM3_H5mm_9771050360` | 180.852 | 11 | 2828 | 6.9115 | 0 | 0 | 0 | 5.96856e-13 | 5.68e-14 | 0.014441 | 5.96856e-13 | 0.014441 | 0 |
| `kicad-lib-RV_Disc_D9mm_W6.1mm_P5mm` | 298.9 | 11 | 16430 | 23.8005 | 0 | 0 | 0 | 2.36518e-07 | 1.71e-07 | 0.0495973 | 2.36518e-07 | 0.0495973 | 0.000379127 |
| `kicad-lib-Converter_ACDC_Hi-Link_HLK-20Mxx` | 45859.9 | 30 | 8330 | 181.571 | 0 | 0 | 0 | 3.63798e-11 | 1.82e-12 | 0.0032956 | 3.63798e-11 | 0.0032956 | 0 |
| `kicad-vendored-L_CommonModeChoke_TDK_ACM7060` | 147 | 6 (15) | 12 | 26 | 0 | 0 | 0 | 5.68434e-14 | 5.68e-14 | 9.93014e-16 | 2.84217e-14 | 0 | 0 |
| `kicad-lib-L_CommonModeChoke_Coilank_ACM2012` | 3.37691 | 27 | 300 | 6.4 | 8.88178e-16 | 0 | 0 | 2.22045e-15 | 1.78e-15 | 0.00077599 | 2.22045e-15 | 0.00077599 | 2.77556e-16 |
| `kicad-lib-LED_D3.0mm_Horizontal_O1.27mm_Z6.0mm` | 44.2398 | 33 | 1624 | 4.8 | 4.26326e-14 | 0 | 0 | 9.9476e-14 | 1.42e-14 | 0.0277456 | 9.9476e-14 | 0.0277456 | 1.85037e-16 |
| `kicad-lib-Buzzer_Mallory_AST1109MLTRQ` | 168.929 | 24 (26) | 564 | 40 | 0 | 0 | 0 | 2.84217e-14 | 1.14e-13 | 0.00390921 | 0 | 0.00390921 | 0 |
| `kicad-lib-Mounting_Wuerth_WA-SMSI-M3_H1mm_9774010360` | 35.9612 | 7 | 448 | 20.923 | 0 | 0 | 0 | 4.90274e-13 | 9.95e-14 | 0.014441 | 4.90274e-13 | 0.014441 | 0 |
| `kicad-lib-R_Box_L13.0mm_W4.0mm_P9.00mm` | 419.016 | 10 | 268 | 34 | 0 | 0 | 0 | 0 | 5.68e-14 | 0.00192174 | 0 | 0.00192174 | 0 |
| `kicad-lib-PinSocket_1x01_P2.54mm_Vertical` | 51.6409 | 36 | 116 | 12.76 | 0 | 0 | 0 | 4.9738e-14 | 1.42e-14 | 0 | 6.39488e-14 | 0 | 0 |
| `kicad-lib-C_Axial_L5.1mm_D3.1mm_P10.00mm_Horizontal` | 38.9307 | 15 | 4212 | 18.6987 | 0 | 0 | 0 | 1.27898e-13 | 7.11e-14 | 0.0105891 | 1.27898e-13 | 0.0105891 | 0 |
| `kicad-vendored-WLP-4_0.83x0.83mm_P0.4mm` | 0.249593 | 12 | 5912 | 3.312 | 5.61495e-14 | 0 | 0 | 1.47937e-14 | 4.44e-16 | 0.00412532 | 1.47937e-14 | 0.00412532 | 1.34085e-16 |
| `kicad-lib-D_SOD-923` | 0.174024 | 27 | 100 | 2.76364 | 0 | 0 | 0 | 3.84304e-13 | 2.87e-12 | 2.66604e-13 | 3.84304e-13 | 2.66604e-13 | 1.6069e-16 |
| `kicad-lib-L_Radial_D28.0mm_P29.20mm` | 17271.3 | 11 | 23866 | 95.5044 | 0 | 0 | 0 | 3.63798e-11 | 5e-12 | 0.0236861 | 3.63798e-11 | 0.0236861 | 8.92787e-16 |
| `kicad-lib-R_0603_1608Metric` | 0.537219 | 20 (26) | 228 | 4.8 | 0 | 0 | 0 | 2.22045e-16 | 0 | 0.000381795 | 2.22045e-16 | 0.000381795 | 0 |
| `kicad-lib-QFN-32-1EP_4x4mm_P0.4mm_EP2.9x2.9mm` | 9.6 | 6 (72) | 12 | 16 | 0 | 0 | 0 | 1.77636e-15 | 2.84e-14 | 3.25321e-15 | 0 | 0 | 0 |
| `kicad-lib-DFN-8-1EP_3x2mm_P0.5mm_EP1.3x1.5mm` | 3.6 | 6 (24) | 12 | 10 | 0 | 0 | 0 | 4.44089e-16 | 3.55e-15 | 8.67112e-16 | 4.44089e-16 | 0 | 0 |
| `kicad-lib-L_CommonModeChoke_Bourns_SRF1260` | 968.75 | 6 (15) | 12 | 50 | 0 | 0 | 0 | 0 | 2.27e-13 | 2.40178e-15 | 0 | 0 | 0 |

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
| **Resolved by CAD-361: STEP assembly product structure.** Next-assembly-usage occurrences with their context-dependent shape representations, item-defined transformations and rigid Cartesian operators, mapped items and identity shape-representation relationships are read into shared definitions, each validated once, and placed depth first. Mirroring or scaling operators, cycles, dangling or unplaced structure, unit-less or non-millimetre contexts and item placements on which ISO 10303-43 and OpenCascade disagree refuse with their kind. Re-export keeps one product per definition and one occurrence per instance. | `occt-two-bodies`, all four `kicad-*board*` exports, `kicad-lib-Samtec_HPM-01-05-x-S_Straight_1x01_Pitch5.08mm` | CAD-361 (pass) |
| **Unsupported: horn-torus blend in an assembly component.** The structure reads, but component `e50ef681-4cb7-11ee-8265-85d007bab819_part` (solid `#110`) has four `TOROIDAL_SURFACE` corner blends with major = minor = 0.125 mm, which the retained reader refuses; the whole import refuses rather than drop the component. OpenCascade reads 6 solids, 133479.60 mm³. | `kicad-lib-Converter_ACDC_MeanWell_IRM-60-xx_THT` | CAD-360 |
| **Resolved by CAD-360: general imported B-rep.** Cones, general revolutions, partial cylinders, polyhedra, spherical caps, B-spline surfaces and curves, surfaces of revolution and extrusion, toroidal blends and flat multi-solid files are imported with their own supports and trims retained. Volume, areas and centroids are integrated over the trimmed faces with the exact evaluators; the error bound is the measured quadrature convergence plus the linear tolerance swept over the boundary area. | `occt-cone-frustum`, `occt-revolved-vase`, both Wuerth standoffs, `RV_Disc`, `Hi-Link`, `Coilank`, `LED_D3.0mm`, `R_Box`, `PinSocket`, `C_Axial`, `WLP-4`, `D_SOD-923`, `L_Radial` | CAD-360 (pass) |
| **Invalid: non-manifold contact.** Four straight edges are each used by four faces. ISO 10303-42 requires every edge of a closed shell to be referenced exactly twice; OpenCascade tolerates the contact. The import reports it as invalid. | `kicad-lib-SW_SPST_B3U-3000P-B` | CAD-360 |
| **Invalid: inner-shell edge off its surface.** A degree-13 rational B-spline edge of the `BREP_WITH_VOIDS` inner shell lies up to 1.10e-5 mm from its cylinder, beyond the declared 2e-6 mm uncertainty and the 1e-5 mm incidence cap. | `kicad-lib-C_Radial_D6.3mm_H7.0mm_P2.50mm` | CAD-360 |
| **Unsupported: surface model.** A `SHELL_BASED_SURFACE_MODEL` with open shells sits beside the solids. There is no surface-body concept, and importing only the solids would drop part of the model. Surface bodies are deferred by product decision. | `kicad-lib-R_LDR_5.2x5.2mm_P3.5mm_Horizontal` | CAD-363 |
| **Resolved by CAD-360: sections through curved imports.** A plane section that a family projector refuses is cut through the solid's exact STEP export read back as a retained boundary, with its hatch. Retained sections publish exact lines and ellipses where the plane meets a plane, cylinder, cone nappe or sphere in one, and chords within the 0.01 mm drawing chord tolerance elsewhere. Every Rust section length agrees with OpenCascade's section of the same source within 9e-16 relative, except the toroidal faces' chords (7.9e-4 relative for `occt-regular-torus`). | the six former CAD-346 fixtures | CAD-360 (pass) |
| **Observation: `FILE_SCHEMA` is not checked.** A valid AP214 body under `CONFIG_CONTROL_DESIGN` or a nonsense schema name imports. Geometry, not header text, decides. This is pinned so that any change is deliberate. | robustness `file_schema_text_is_not_geometry_evidence` | CAD-290 (decision) |
| **Observation: STEP parsing has no interrupt.** Cancellation is checked before the read and again before commit, so a cancelled request never commits. A long parse cannot be interrupted midway. | robustness `cancelled_import_commits_nothing` | CAD-290 |

The OpenCascade reference mesh is drawn at 0.01 mm. At the Rust deflection, 0.05 mm, OpenCascade's mesh of `L_Radial`'s toroidal face lies 0.134 mm from a point of the exact surface, while the Rust mesh stays within 0.047 mm of the exact surface under dense barycentric sampling. The comparison bound is the sum of the two deflections, 0.06 mm.

A retained-boundary re-export of a cone frustum is read back by the cone family, which rebuilds it as a family body whose section drawing refuses curved faces; the source section of the retained import succeeds. Section lengths of retained imports are polyline lengths within the drawing chord tolerance, not exact circle lengths.

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

The sweep covered 7236 STEP models in KiCad 10.0.0 at `b333779e` (20.4 min wall, 16.9 min summed import time, debug build, host shared with OpenCascade and test builds); an earlier run at `6b26ab78` gave the same outcome for every file. Every file was classified; there were no panics or unclassified errors. A pass here means the import succeeded: volume, bounds and topology from the retained or recognised boundary. The full corpus workflow (mesh, section, persistence, export, OpenCascade) runs only for the pinned fixtures.

| Producer (`FILE_NAME`) | Files | Pass | Unsupported | Invalid | Limit |
| --- | ---: | ---: | ---: | ---: | ---: |
| KiCad | 5181 | 4642 | 214 | 325 | 0 |
| FreeCAD | 1031 | 832 | 79 | 120 | 0 |
| unstated | 796 | 796 | 0 | 0 | 0 |
| OpenCascade | 211 | 108 | 103 | 0 | 0 |
| Solid Edge / ST-DEVELOPER v16.7 | 14 | 0 | 14 | 0 | 0 |
| Autodesk Translation Framework v7.6.0.251 / ST-DEVELOPER v17.2 | 2 | 1 | 1 | 0 | 0 |
| Autodesk Translation Framework v13.14.0.145 / ST-DEVELOPER v20 | 1 | 0 | 1 | 0 | 0 |
| Total | 7236 | 6379 | 412 | 445 | 0 |

| Category (any surface or structure present) | Files | Pass | Unsupported | Invalid | Limit |
| --- | ---: | ---: | ---: | ---: | ---: |
| cylinder | 2902 | 2888 | 0 | 14 | 0 |
| torus | 1317 | 955 | 42 | 320 | 0 |
| spline | 1006 | 706 | 24 | 276 | 0 |
| sphere | 956 | 751 | 156 | 49 | 0 |
| planar | 866 | 853 | 0 | 13 | 0 |
| cone | 555 | 408 | 88 | 59 | 0 |
| multi-body | 537 | 117 | 385 | 35 | 0 |
| revolution | 532 | 264 | 37 | 231 | 0 |
| assembly | 382 | 0 | 382 | 0 | 0 |
| voids | 211 | 87 | 4 | 120 | 0 |
| surface-model | 2 | 0 | 2 | 0 | 0 |

At `dff6a0f8` the sweep gave 173 pass, 5,799 unsupported, 299 invalid and 965 limit; at `1f1c4fe7` (CAD-359) 219 pass, 7,007 unsupported, 10 invalid and 0 limit; at `5fc4f5df` (CAD-360) 6,115 pass, 1,115 unsupported, 6 invalid and 0 limit; at `046ad81a` 6,125 pass, 640 unsupported and 471 invalid, after the retained-boundary reader began to report boundaries that fail their own incidence as invalid. 382 unsupported files are assemblies (CAD-361).

Between `046ad81a` and this sweep, pass grew by 254 files, unsupported fell by 228 and invalid by 26, from horn tori, exact-area loop orientation, the stated-vertex fallback and the declared-uncertainty classification. Within that, three files with sliver faces (`ESP-WROOM-02`, `Conn_C14_Receptacle_RightAngle_Schurter_DD21.01xx`, `USB_C_Receptacle_GCT_USB4085`), which OpenCascade 7.9.3 accepts, moved from invalid to unsupported, and the three `L_Axial` Fastron 77A, MICC and SMCC inductors from invalid to pass: every stated vertex lies within 4.98e-6 mm of each of its supports (D = 5e-6 mm) and the evaluated edges within 1e-8 mm of their faces by the independent evaluator below. OpenCascade's own shape check rejects those three, for a reason not determined here, so their volumes have no OpenCascade comparison. A file is invalid only when it violates its own declared uncertainty or the STEP topology rules. The 445 invalid files fall in these classes. Distances are the first violating point the reader found (a lower bound of the largest one) against the file's declared uncertainty D. "Independent" is a separate evaluator of the STEP entities (de Boor for splines, closed forms for analytic surfaces), written for this check and sharing no code with the reader.

| Class | Files | Measured d/D | Independent evidence |
| --- | ---: | --- | --- |
| Edge or loop point off its face | 208 | 4 at 1.0–1.01, 66 at 1.01–1.5, 86 at 1.5–3, 40 at 3–10, 12 above 10 | 8 of 8 sampled at random (seed 360) and the 4 nearest 1.0 confirmed: largest edge-to-surface distance 3.7–13 D (`CP_Axial_L30.0mm_D18.0mm` 1.85e-5 mm with D 5e-6; `CP_Radial_Tantal_D7.0mm` 1.99e-5 mm with D 2e-6) |
| Vertex off its edge curve | 140 | 1 at 1.0–1.01, 35 at 1.01–1.5, 45 at 1.5–3, 15 at 3–10, 44 above 10 | 8 of 8 sampled confirmed; the one nearest 1.0, `L_Axial_L9.5mm_D4.0mm_P12.70mm`, at 5.0039e-6 mm with D 5e-6 |
| Quadrature point off its face | 38 | 3 at 1.0–1.01, 30 at 1.01–1.5, 5 at 1.5–3 | sampled files confirmed; the three nearest 1.0 at 1.0098 D (`C_Disc_D14.5mm_W5.0mm`, maximised along the edge), 1.126 D (`C_Disc_D7.5mm_W2.5mm`) and 1.77 D (`Sipeed-M1W`, surface of revolution) |
| Edge used by four faces | 48 | — | every file: one `EDGE_CURVE` entity is referenced by four faces of one `CLOSED_SHELL`, counted from the source text. ISO 10303-42 allows two. |
| Re-derived vertex beyond a support | 5 | — | every file: a stated vertex lies 3.2–6.6 D from one of its supports (`SW_SPST_TL3342`, `Resonator_Murata_CSTLSxxxX`, `Heatsink_AAVID_573300D00010G_TO-263`, `PQFP-256`, `MSOP-8-1EP`) |
| Loop not head to tail | 2 | — | `Fuseholder_Schurter_0031.7701` joins two edges 0.346 mm apart; `ST_HLGA-6` joins two distinct vertex entities at one point, against the `edge_loop` rule |
| Edge ends out of order | 2 | — | `L_CommonMode_Toroid_Vertical_L39.4mm` and `L41.9mm` (D 2e-4): spline edges whose vertices lie 5e-4 mm off their curves |
| Single-edge loop between two vertices | 1 | — | `BNC_Amphenol_B6252HB-NPP3G-50`: a circle edge from one vertex entity to another at the same point closes a loop alone |
| Family loops not head to tail | 1 | — | `Diode_Bridge_GeneSiC_KBPC_T`: 24 loop joints up to 2.158 mm apart |

OpenCascade loaded 398 of the 471 files invalid at `046ad81a` as valid (measured on that list before this change); its own checks heal or tolerate gaps within its tolerances, non-manifold edges and coincident vertex entities. The Rust reader does not heal: these files are invalid by their own declared terms, not because of a reader defect, and repairing them belongs to a healing step that is not implemented.

At `ca31f83c` (CAD-361, default build, 12.6 min) the sweep gave 6,462 pass, 258 unsupported, 516 invalid and 0 limit, with no panic or unclassified error. Of the 382 assembly files, 334 pass, 45 are refused as invalid and 3 as unsupported, all with the retained-boundary message; the per-category tables above are from `5fc4f5df`. Most of the rise in invalid results comes from `c0fddf10` (CAD-360), which reports an invalid retained boundary instead of the families' unsupported refusal. The sweep does not separate refusals of the product structure from refusals of a component's boundary.

Of the ten CAD-359 invalid results, `Heatsink_AAVID_Extruded_531002B02100G`, `…531102B02100G`, `…531202B02100G` and `R_Shunt_Ohmite_LVK12` import with their declared 2e-5 mm uncertainty (capped at 1e-5 mm, vertices re-derived on their supports). `ESP32-WROOM-32E` (vertices 2.15e-7 mm off, declared 2e-7 mm), `Relay_DPDT_FRT5_SMD` (1.95e-3 mm off, declared 5e-4 mm) and `Diode_Bridge_GeneSiC_KBPC_T` (loops not head to tail) stay invalid. The three `HDMI_Micro-D_Molex_46765` files were invalid at `5fc4f5df` because their cone-plane edges are `HYPERBOLA` curves; `9fbf5da1` reads them, and OpenCascade 7.9.3 agrees with their imports on volume within 1.3e-9 relative.

Sampled cross-check at `5fc4f5df`/`9a3fdda5`: 80 randomly drawn (seed 360) sweep passes, all retained-boundary imports, were compared with OpenCascade 7.9.3. Volume agreed within 1.2e-8 relative (median below 1e-14), bounds within 1e-7 mm, and every certified Rust mesh lay within 0.04 mm of OpenCascade's 0.01 mm tessellation. Three of the 80 manufacturing meshes are refused rather than produced: `USB_A_TE_292303-7_Horizontal` and `FanPinHeader_1x04_P2.54mm_Vertical` (a planar face whose ear clipping leaves a run of collinear boundary vertices uncovered) and `UQFN-20-1EP_4x4mm_P0.5mm_EP2.8x2.8mm` (certification: unresolved boundary). The first one panicked before `9a3fdda5`.

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
- KiCad boards are always exported as assemblies; they pass through the CAD-361 assembly reader. Analytic oracles are pinned for four of the five vendored assemblies. `kicad-slotted-rounded-board` has none: its transmitted corner arcs give 3333.1659232157 mm³, which OpenCascade confirms (3333.1659232183), 1.3e-9 relative from the 2100 + 25π mm² design outline times 1.53 mm, beyond the 1e-9 oracle tolerance.
- For a multi-body source the OpenCascade section reference cuts each transmitted solid on its own at the Rust plane. OpenCascade's section of the whole compound of `Samtec_HPM`, whose pin passes through its solid housing, measured 4.56 mm against the 27.42 mm sum of both bodies' cuts (22.86 mm housing, 4.56 mm pin); per solid it agrees exactly.

## Assembly reference (CAD-361)

`CARGO_TARGET_DIR=target/occt cargo test --locked -p multicad-app --features occt --test step_assembly -- --nocapture` compares every occurrence with OpenCascade 7.9.3. Each Rust world occurrence must match one OpenCascade solid (`STEPControl_Reader`) by volume within 1e-6 relative and bounds within 1e-3 mm, and each next-assembly-usage occurrence, found by its name or identifier, must place the same number of solids with the same volume and bounds in its parent's coordinates as the shape OpenCascade binds to it. The six passing corpus assemblies, a repeated component (one box definition placed twice), a nested placement (5, 6, 7) turned 90°, and a Cartesian-operator placement were checked, each also after its structured re-export: the largest volume gap was 7.7e-13 relative, the largest bounds gap 2.2e-16 mm, and every instance name and solid count agreed.
