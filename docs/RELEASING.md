# GitHub CI and desktop releases

## Repository and publication

GitHub is the primary repository at [kellertobias/TOSKCAD](https://github.com/kellertobias/TOSKCAD), as selected by the user on October 31, 2026. The release flow follows ToskLight: GitHub Actions tests and builds the application, then creates a `vVERSION` tag and hosted GitHub Release on the exact tested `main` commit. The publisher creates no release commit and leaves tracked manifests and changelogs unchanged. Forgejo mirroring is not part of this flow.

The release job alone has `contents: write` and uses the ephemeral GitHub repository token. Publication requires the repository variable `MULTICAD_RELEASES_ENABLED=true`. Unset or false keeps publication skipped while all checks and builds still run. Main/tag protection must permit this token to create release tags and assets. Tokens and signing keys never belong in source control.

## Version contract

| Commit | Release |
| --- | --- |
| `fix:` / `perf:` | Patch |
| `feat:` | Minor |
| `type!:` / `BREAKING CHANGE:` | Major |
| `docs:`, `test:`, `style:`, `refactor:`, `build:`, `ci:`, `chore:` without breaking intent | None |

Git tags are authoritative release versions. `Cargo.toml`/`Cargo.lock` retain the development workspace version; the private npm manifest has no product version and is never published. `MULTICAD_RELEASE_VERSION` and `MULTICAD_BUILD_COMMIT` stamp the built binary without changing tracked files. Archive names, embedded diagnostics and release manifests must agree with the selected version and tested SHA.

Like ToskLight, a history with no semantic tag uses a disposable `v0.0.0` analysis baseline on its root commit. It is removed before publication and never exported as a real release. The first feature selects `0.1.0`; a first patch selects `0.0.1`. A setup-only commit does not select a new version.

## Pipeline

`.github/workflows/release.yml` validates branch pushes, pull requests and manual builds. Tag pushes do not recursively trigger it.

- Metadata fetches full history and previews a release only for a main-branch push; other events use the latest reachable version or development version.
- Quality runs formatting, Clippy, actual release-engine contract tests and prototype checks.
- The reusable owned-kernel workflow tests and strictly lints the owned Rust exact-arithmetic and geometry crates on all three operating systems using the dedicated coverage profile. It gates publication alongside the workspace test archive runs. Configured jobs are separate from actual hosted-run evidence.
- Tests run under [cargo-nextest](https://nexte.st), pinned in `tools/native-test-shard.mjs`. CI installs that exact release as a prebuilt binary through `taiki-e/install-action`; the manual regression workflow still compiles it with `cargo install --locked`. nextest runs every test in its own process. Tests use the optimized `ci` Cargo profile (`--cargo-profile ci`) with debug assertions and overflow checks retained. Shipped executables use the normal `release` profile.
- **Compile tests** (one job per platform: Linux x64, macOS ARM64, Windows x64) builds every workspace test of the exact commit once, then splits the test binaries into four nextest archives, one per test runner (`node tools/native-test-shard.mjs archive DIRECTORY`). Each binary goes into exactly one archive, assigned longest first by the per-binary test seconds in `tools/test-weights.json` (test count for unlisted binaries). The weights only balance runtime; a stale table never drops a test. `partitions.json` records which binaries each archive holds. The job compiles with `CARGO_TARGET_DIR` set to the workspace's `target/` rather than the configured `target/cargo`: nextest stores and extracts the target directory as `target/`, and tests start `multicad` through its compile-time path. The same job runs the doctests, which nextest does not run, with `cargo test --doc`, plus the headless budget smoke in the `ci` profile (one run; not a performance result).
- **Test** runs one archive per runner, four runners per platform (`node tools/native-test-shard.mjs shard ARCHIVE`), so the four runners together run every workspace test exactly once. A single archive of all binaries would be about 7 GB compressed and 21 GB unpacked, too much for one hosted macOS runner; each split archive is about a quarter of that. The archive is extracted so its `target/` lands at the path it was built in, because tests that launch built executables use compile-time paths into it.
- **Build and package** builds, smoke-checks and archives the version-stamped executable on each platform in parallel with the tests, then produces versioned ZIPs, per-target identity records and SHA-256 checksums. With `MULTICAD_PACKAGE_TESTS=gated` it skips its own test run, because the test partitions and doctests of the same commit gate publication. Local packaging keeps the default and runs the full suite first.
- Jobs that compile restore compiled dependencies with `Swatinem/rust-cache`, keyed on the lockfile, toolchain and platform. Every job still checks out and builds the exact commit; the cache never supplies workspace crates or test results.
- `.github/workflows/native-regression.yml` keeps the named focused suites for diagnosing a failure by area. It runs only by hand (`workflow_dispatch`) and does not gate publication: every test it runs is part of the partitioned complete suite.
- Publication waits for every quality, package, test-archive and owned-kernel result and the publication activation variable. It verifies all three archives' version/SHA/target/checksum records, rejects missing/stale/extra files and writes the aggregate release manifest.
- The publisher verifies the checkout/current GitHub main SHA and the expected version again before tagging. It creates no branch commit. Hosted release assets are only the verified ZIPs, checksums and aggregate manifest; success/failure issue comments are disabled.

Packaging produces **native development builds**, with Parameters, constrained Sketch, initial modeling tools and native project save/recovery. Advanced workspaces and manufacturing exports remain under development; the local MCP server is included (`multicad --mcp`); current interactive evidence is recorded in the [modeling checkpoint](verification/native-modeling-2026-10-28.md). Release jobs package the MultiCAD Rust kernel in each runtime ZIP, with a matching source ZIP containing exact application source, vendored Cargo dependencies and offline rebuild instructions. Every normal package is checked for retired geometry archives, libraries and prototype assets, embedded version/commit/target identity, Rust kernel diagnostics, runtime library imports and statically linked OpenCascade code markers. Each packaging job also runs the archived executable's self-test (every bundled example opens) and a packaged MCP workflow: `multicad --mcp` opens a copy of the bundled drill-guide example, must report `geometry_backend: rust`, evaluates it and exports its largest part as a complete STEP file. There is one packaging mode. Aggregate verification repeats these archive checks for all three targets and requires both archives and their checksums; macOS uses ad-hoc signing without notarization. Node.js/npm are CI/developer release tools, not native runtime dependencies. The pinned release-engine/plugin versions match the inspected ToskLight lockfile.

## Local verification without publishing

```sh
npm ci --ignore-scripts --no-audit --no-fund
npm test
npm run release:dry-run
```

The preview uses disposable local clones/bare storage; it neither changes source refs/index nor contacts GitHub/Forgejo for publication. Semantic-release may log an internal push-permission dry run against that temporary local storage; it does not push this repository. Tests exercise patch/feature/breaking/non-release histories and an already-tagged retry, plus publication/asset rejection conditions.

From a clean committed source checkout, `npm run release:package -- 0.1.0` produces the current host's native runtime and corresponding-source archives under `artifacts/release/`. It runs native tests, builds and executes that binary, and validates its version/SHA/target before archiving. The application build uses the pinned Rust toolchain, and the tests need the pinned cargo-nextest (`cargo install cargo-nextest --version 0.9.152 --locked`). macOS also requires Apple Command Line Tools for the AppKit/Foundation window adapter. Ubuntu requires `libwayland-dev` and `libxkbcommon-dev`. This is local artifact generation, not a release. Remove/move older local output versions before aggregating a complete package set; verification deliberately rejects stale extras.

Completed releases are no-ops on repeated semantic analysis. An interrupted provider upload can leave a tag/draft requiring explicit inspection and repair of that same release; do not move/delete public tags or mint a different version to conceal it. Check the actual workflow run and hosted assets before reporting a release as available. CI builds and smoke tests are separate from interactive Windows/Linux testing.

## Project website publication

GitHub Pages serves `website/` at https://kellertobias.github.io/TOSKCAD/. The separate `.github/workflows/pages.yml` publishes a static artifact on website changes to `main` or a manual dispatch. Its deployment job has `pages: write` and `id-token: write`; application release activation remains independent.

The [download page](https://kellertobias.github.io/TOSKCAD/releases.html) reads published releases from the public GitHub API. It matches the versioned native archive names, links their corresponding sources and checksums, and displays an empty or error state when downloads are unavailable. JavaScript is required for the live list; a direct GitHub Releases link remains available.
