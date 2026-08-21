# Contributing to prakash

Thank you for considering a contribution! This document covers the development
workflow, coding standards, and review process. prakash is written in
[Cyrius](https://github.com/MacCracken/cyrius).

## Development Setup

```bash
# Clone the repository
git clone https://github.com/MacCracken/prakash.git
cd prakash

# Install the Cyrius toolchain pinned in cyrius.cyml (the `cyrius = "..."` line).
# Use the upstream installer, which verifies the release checksum and signature
# and lays the toolchain out as ~/.cyrius/versions/<pin>/ — the layout `cyrius
# deps` resolves its stdlib snapshot from. CI does exactly this.
CYRIUS_VERSION="$(grep -oP '(?<=^cyrius = ")[^"]+' cyrius.cyml)"
curl -sSfL "https://raw.githubusercontent.com/MacCracken/cyrius/$CYRIUS_VERSION/scripts/install.sh" \
  | CYRIUS_VERSION="$CYRIUS_VERSION" sh

# Resolve the hisab git dependency (the stdlib subset is vendored under lib/).
cyrius deps
```

## Local CI gate

CI runs on every push (`.github/workflows/ci.yml`). Reproduce it locally before
pushing:

```bash
# Dependency hashes must match the committed cyrius.lock
cyrius deps --verify

# Lint + format (must be clean — 0 warnings, 0 drift).
# NOTE the glob includes examples/ — CI checks those too.
# ⚠ `cyrius fmt <f>` REWRITES IN PLACE as of cyrius 6.5.28 (it was stdout-only
# before). To fix drift run `cyrius fmt "$f"` — never `cyrius fmt "$f" > tmp`,
# which captures zero bytes and truncates the file.
for f in src/*.cyr tests/*.tcyr tests/*.bcyr examples/*.cyr; do
  cyrius lint "$f"          # fails on a `warn` line
  cyrius fmt "$f" --check   # fails on format drift (writes nothing in check mode)
done

# Vet the build entry
cyrius vet src/main.cyr

# Bundles must match lib/ (regenerate + commit if they drift).
# Sync the .deps sidecars after distlib: `cyrius distlib` gets them backwards
# for prakash's inverted two-bundle layout (base = narrow math-only, profile =
# wide + TLS). They are published metadata, not an enforced contract — see the
# script header for what was actually measured.
cyrius distlib && cyrius distlib ai
./scripts/sync-deps-sidecar.sh
git diff --quiet dist/ || echo "dist/ stale — commit it"

# Tests + benchmarks
for f in tests/*.tcyr; do cyrius test "$f"; done
cyrius bench tests/prakash.bcyr
```

## Pull Request Process

1. **Fork and branch** — create a feature branch from `main`.
2. **Keep commits focused** — one logical change per commit.
3. **Write tests** — new features require tests (a `.tcyr` suite); bug fixes
   require a regression test.
4. **Regenerate bundles** — if you touch a `[lib]`/`[lib.ai]` module, run
   `cyrius distlib && cyrius distlib ai && ./scripts/sync-deps-sidecar.sh` and
   commit the updated `dist/` (both `*.cyr` bundles **and** the `*.deps` sidecars).
5. **Run the local CI gate** (above) — lint, fmt, tests, and bench must be green.
6. **Open a PR** against `main` with a clear description of the change.

## Code Style

- Zero `cyrius lint` warnings; `cyrius fmt --check` clean (the linter counts
  **bytes** — keep lines ≤ 120).
- `#must_use` on pure functions, `#derive(accessors)` on structs. Error and
  status codes are module-level integer constants (`var PK_ERR_* = ...` in
  `src/error.cyr`): `0` (`PK_ERR_NONE`) = success, negative = error — there are
  no enums.
- No `unwrap`/`panic` equivalents in library code — return `PK_ERR_*` via an
  `err_out` pointer.
- Physics functions document their units.
- **Bit-fidelity to the Rust source** is the porting contract: encode `f64`
  constants as exact ratios (numerator/denominator < 2^53) or IEEE-754 hex bit
  patterns; replicate `powi` as square-and-multiply; preserve left-associative
  fold order. When in doubt, compare against the Rust original and add a test.
  It is no longer in-tree (removed in 2.2.3) — retrieve it with
  `git checkout 2.2.2 -- rust-old/`, or `git show 1.2.0:src/<file>.rs`.

## Testing Requirements

- All public API changes need unit tests in a `tests/*.tcyr` suite.
- Performance-sensitive changes: run `./scripts/bench-history.sh` before and
  after and cite the CSV numbers — **never claim a speedup without benchmarks.**
- Target: maintain or improve coverage (80%+) on every PR.

## Commit Messages

Use clear, imperative-mood commit messages:

```
add Sellmeier dispersion coefficients for sapphire

fix NaN handling in wavelength_to_rgb boundary check
```

## License

By contributing, you agree that your contributions will be licensed under the
same license as the project (see [LICENSE](LICENSE)).
