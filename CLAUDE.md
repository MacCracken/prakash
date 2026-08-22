# Prakash — Claude Code Instructions

## Project Identity

**Prakash** (Sanskrit: light/illumination) — Optics and light simulation — ray optics, wave optics, spectral math, lens geometry, PBR

- **Language**: Cyrius (ported from Rust; the `rust-old/` archive was removed in 2.2.3)
- **Type**: Cyrius library — `src/` modules bundled by `cyrius distlib` into `dist/prakash.cyr`
  (math-only) and `dist/prakash-ai.cyr` (adds the sandhi HTTP/TLS AI client)
- **License**: GPL-3.0
- **Toolchain pin**: cyrius 6.5.33 — single source of truth is `cyrius = "..."` in `cyrius.cyml`
- **Version**: SemVer, currently 2.2.x. `VERSION` is the single source; `cyrius.cyml`
  reads it via `${file:VERSION}` and the release workflow fails if the tag disagrees

⚠ **The Rust original is a FIDELITY reference, not a CORRECTNESS one.** 2.2.6 found six
wrong physics formulas, five of them character-for-character identical in the archive
(`git show 524a7aa^:rust-old/<path>.rs`). "It matches Rust" is never sufficient evidence
that a formula is right — check the literature (A&S, Coddington, CIE) instead.

## Consumers

soorat (PBR shading), kiran (lighting), ranga (lens effects)

## Development Process

### P(-1): Scaffold Hardening (before any new features)

0. Read roadmap, CHANGELOG, and open issues — know what was intended before auditing what was built
1. Test + benchmark sweep of existing code
2. Cleanliness check: `cyrius fmt <file> --check` (every `src/*.cyr tests/*.tcyr tests/*.bcyr examples/*.cyr`),
   `cyrius lint <file>` (⚠ exits 0 even with warnings — grep stdout for `^\s*warn `),
   `cyrius doc --check <file>`, `cyrius vet src/main.cyr`, `cyrius deny src/main.cyr`,
   `cyrius deps --verify`, and `cyrius audit` (the project sweep: fmt/lint/docs/tests/bench —
   it must exit **0**)
3. Get baseline benchmarks (`./scripts/bench-history.sh`)
4. Internal deep review — gaps, optimizations, security, logging/errors, docs
5. External research — domain completeness, missing capabilities, best practices, world-class accuracy
6. Cleanliness check — must be clean after review
7. Additional tests/benchmarks from findings
8. Post-review benchmarks — prove the wins
9. Repeat if heavy

### Work Loop / Working Loop (continuous)

1. Work phase — new features, roadmap items, bug fixes
2. Cleanliness check: `cyrius fmt <file> --check` (every `src/*.cyr tests/*.tcyr tests/*.bcyr examples/*.cyr`),
   `cyrius lint <file>` (⚠ exits 0 even with warnings — grep stdout for `^\s*warn `),
   `cyrius doc --check <file>`, `cyrius vet src/main.cyr`, `cyrius deny src/main.cyr`,
   `cyrius deps --verify`, and `cyrius audit` (the project sweep: fmt/lint/docs/tests/bench —
   it must exit **0**)
3. Test + benchmark additions for new code
4. Run benchmarks (`./scripts/bench-history.sh`)
5. Internal review — performance, memory, security, throughput, correctness
6. Cleanliness check — must be clean after audit
7. Deeper tests/benchmarks from audit observations
8. Run benchmarks again — prove the wins
9. If audit heavy → return to step 5
10. Documentation — update CHANGELOG, roadmap, docs
11. Version check — `VERSION` bumped; `cyrius.cyml` picks it up via `${file:VERSION}`
12. Regenerate the bundles — `cyrius distlib && cyrius distlib ai && ./scripts/sync-deps-sidecar.sh`,
    then confirm `git diff dist/` is what you expect. ⚠ `cyrius distlib` OVERWRITES the
    `.deps` sidecars with a backwards base-vs-profile inference; the sync script must run
    after it, every time (CI does exactly this and diffs)
13. Return to step 1

### Task Sizing

- **Low/Medium effort**: Batch freely — multiple items per work loop cycle
- **Large effort**: Small bites only — break into sub-tasks, verify each before moving to the next. Never batch large items together
- **If unsure**: Treat it as large. Smaller bites are always safer than overcommitting

### Refactoring

- Refactor when the code tells you to — duplication, unclear boundaries, performance bottlenecks
- Never refactor speculatively. Wait for the third instance before extracting an abstraction
- Refactoring is part of the work loop, not a separate phase. If a review (step 5) reveals structural issues, refactor before moving to step 6
- Every refactor must pass the same cleanliness + benchmark gates as new code

### Key Principles

- **Never skip benchmarks.** Numbers don't lie. The CSV history is the proof.
  ⚠ This host's run-to-run spread reaches 40% on sub-100 ns benchmarks. A whole-suite
  delta below ~10% is NOT a signal — put both versions in ONE binary and bench them
  against each other before claiming any movement.
- **Tests + benchmarks are the way.** Minimum 80%+ coverage target. ⚠ `cyrius coverage`
  reports *reference* coverage — whether a function is called at all. It is a floor,
  not a correctness proof.
- **A test that samples only the points a bug survives is worse than no test.** 2.2.6's
  six defects all passed a 5537-assertion suite: Fresnel was sampled at x = 0, 1, 10
  (two are fixed points of any approximation), every Seidel assertion wrapped its value
  in `f64_abs` so a sign error was invisible, and the doublet's entire coverage was a
  length and a name string. Pin VALUES against an independent reference, and pin the
  invariants (sign, monotonicity, bounds) separately.
- **Own the stack.** If an AGNOS library wraps an external one, depend on the AGNOS library.
  `lib/` is vendored (ganita, hisab, sakshi, sandhi, bayan, ...) — never edit it; a defect
  there is fixed upstream or worked around in `src/` with the reason written down.
- **No magic.** Every operation is measurable, auditable, traceable.
- **`#must_use`** on all pure functions (399 uses — the load-bearing attribute here).
- **`#derive(accessors)`** for struct field accessors rather than hand-written `load64`.
- ⚠ **`cycc` SILENTLY IGNORES UNKNOWN ATTRIBUTES** — `#definitely_not_real` compiles and
  lints clean, exactly like `#inline` does. Cyrius has no `#inline` and no enums, so the
  Rust-era `#[inline]` / `#[non_exhaustive]` rules are gone. Real attributes, as used by
  the vendored stdlib: `#must_use`, `#pure`, `#derive`, `#naked`, `#io`, `#host_only`,
  plus the `#ifdef` / `#ifndef` / `#define` / `#endif` preprocessor.
- **Write into a caller buffer, not a fresh allocation** — `fmt_*_buf` and
  `syscall(1, 1, ...)` over anything that allocates a temporary. (The Cyrius reading of
  the old "`write!` over `format!`".)
- **Borrow the pointer, don't copy the struct** — the bump allocator never frees, so a
  per-call scratch buffer is a leak. Prefer stack locals (`var f = 0; ... &f`).
- **Vec arena over hashmap** — when indices are known, direct access beats hashing.
- **Profile-gate optional deps** — `[lib]` vs `[lib.ai]` in `cyrius.cyml`; consumers pick
  ONE bundle, and `dist/prakash.cyr` must stay free of sandhi/TLS/net symbols (CI greps).
- **sakshi tracing on operations** — `_prk_trace("op")` at operation entry, level set by
  `prakash_set_log_level`.
- **Report, don't fabricate.** A function that cannot compute its answer returns a
  `PK_ERR_*` code through `err_out`, or a 0 sentinel — never a plausible-looking number.
  A null handle must be distinguishable from a computed zero.

## DO NOT
- **Do not commit or push** — the user handles all git operations (commit, push, tag)

- **NEVER use `gh` CLI** — use `curl` to GitHub API only
- Do not add unnecessary dependencies — keep it lean
- Do not edit anything under `lib/` — it is vendored by `cyrius deps`/`cyrius lib sync`
- Do not store through an `alloc()` result without checking it for 0 first, and do not
  index a buffer without a bounds check — Cyrius has no `unwrap`/`panic`, so the failure
  mode is a SIGSEGV that takes the process down. This is the 2.0.2 defect class.
- Do not skip benchmarks before claiming performance improvements
- Do not commit `build/` — but `cyrius.lock` IS committed and CI verifies it
  (`cyrius deps --verify`), and `dist/*.deps` are deliberately tracked
- Do not write temporary probes or scratch files under `src/` or `tests/` — put them in
  `/tmp` and run them with an absolute path from the repo root (includes resolve from
  there). Files left in-tree get swept into commits.

## Documentation Structure

```
Root files (required):
  README.md          — quick start, features, dependency stack, consumers, license
  CHANGELOG.md       — per-version changes (Added/Changed/Fixed/Removed)
  CLAUDE.md          — this file (development process, principles, DO NOTs)
  CONTRIBUTING.md    — fork, branch, run the cyrius gate, PR workflow (there is no Makefile)
  SECURITY.md        — supported versions, scope, reporting
  CODE_OF_CONDUCT.md — Contributor Covenant
  LICENSE            — GPL-3.0

docs/ (required):
  architecture/
    overview.md      — module map, data flow, consumers, dependency stack
    math.md          — (if applicable) mathematical reference for algorithms/formulas
  development/
    roadmap.md       — completed items, backlog, future features (demand-gated), v1.0 criteria

docs/ (when earned — not scaffolded empty):
  adr/
    NNN-title.md     — architectural decision records (when non-obvious choices are made)
  development/
    threat-model.md  — attack surface, mitigations (when security-relevant)
    dependency-watch.md — deps to monitor for updates/CVEs
  guides/
    usage.md         — patterns, philosophy, code examples
    testing.md       — test count, coverage, testing patterns

ADR format:
  # NNN — Title
  ## Status: Accepted/Superseded
  ## Context: Why this decision was needed
  ## Decision: What we chose
  ## Consequences: Trade-offs, what changes
```

## CHANGELOG Format

Follow [Keep a Changelog](https://keepachangelog.com/):

```markdown
# Changelog

## [Unreleased]
### Added — new features
### Changed — changes to existing features
### Fixed — bug fixes
### Removed — removed features
### Security — vulnerability fixes
### Performance — benchmark-proven improvements (include numbers)

## [X.Y.Z] - YYYY-MM-DD
### Added
- **module_name** — what was added and why
### Changed
- item: old behavior → new behavior
### Fixed
- issue description (root cause → fix)
### Performance
- benchmark_name: before → after (−XX%)
```

Rules:
- Every PR/commit that changes behavior gets a CHANGELOG entry
- Performance claims MUST include benchmark numbers
- Breaking changes get a **Breaking** section with migration guide
- Group by module when multiple changes in one release
- Link to ADR if a change was driven by an architectural decision
