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

⚠ **CI runs on pushes to `main` and on PRs targeting `main`**
(`.github/workflows/ci.yml`) — **not** on a feature-branch push. Nothing runs
until you open the PR, so reproduce the gate locally first:

```bash
# Dependency hashes must match the committed cyrius.lock
cyrius deps --verify

# Lint + format (must be clean — 0 warnings, 0 drift).
# NOTE the glob includes examples/ — CI checks those too.
# ⚠ `cyrius fmt <f>` REWRITES IN PLACE as of cyrius 6.5.28 (it was stdout-only
# before). To fix drift run `cyrius fmt "$f"` — never `cyrius fmt "$f" > tmp`,
# which captures zero bytes and truncates the file.
# ⛔ `cyrius lint` EXITS 0 EVEN WITH WARNINGS — verified at 6.6.4, 6.6.6, 6.6.20, 6.7.9 and 6.7.18, a file with two
# `warn` lines still returns rc=0. Running it bare, as this block used to, gives a
# green local gate over code CI then rejects. CI greps stdout instead; so must you.
fail=0
for f in src/*.cyr tests/*.tcyr tests/*.bcyr examples/*.cyr; do
  out=$(cyrius lint "$f" 2>&1)
  if [ $? -ne 0 ] || printf '%s' "$out" | grep -qE '^[[:space:]]*warn '; then
    printf '%s\n' "$out" | grep -E '^[[:space:]]*warn '
    echo "lint: $f"; fail=1
  fi
  cyrius fmt "$f" --check >/dev/null 2>&1 || { echo "fmt drift: $f"; fail=1; }
  cyrius doc --check "$f" >/dev/null 2>&1 || { echo "undocumented fns: $f"; fail=1; }
done
[ $fail -eq 0 ] || echo "cleanliness gate FAILED"

# Vet the build entry, and enforce the project policies CI enforces
cyrius vet src/main.cyr
cyrius deny src/main.cyr

# ⚠ `#must_use` binds to the declaration that FOLLOWS it, so inserting a function
# above another silently steals its attribute and NO other gate reports it.
scripts/check-must-use.sh --selftest
scripts/check-must-use.sh

# Bundles must match lib/ (regenerate + commit if they drift).
# Sync the .deps sidecars after distlib: `cyrius distlib` gets them backwards
# for prakash's inverted two-bundle layout (base = narrow math-only, profile =
# wide + TLS). They are published metadata, not an enforced contract — see the
# script header for what was actually measured.
cyrius distlib && cyrius distlib ai
./scripts/sync-deps-sidecar.sh
git diff --quiet dist/ || echo "dist/ stale — commit it"

# The core bundle must link as a consumer links it: only dist/prakash.deps's
# stdlib leaves, NO hisab. ranga depends on this, and the examples cannot see it
# (they run inside this manifest, which declares hisab).
./scripts/check-consumer-link.sh

# Tests + benchmarks
for f in tests/*.tcyr; do cyrius test "$f"; done
cyrius bench tests/prakash.bcyr

# The core bundle must stay free of the TLS stack — CI greps for exactly this.
grep -nE '\bsandhi_|\btls_[a-z]|\bhttp_[a-z]|\btcp_|\bdynlib_|\brand_[a-z]' \
  dist/prakash.cyr && echo "TLS leaked into the math-only bundle"

# Or run the whole project sweep in one command — it must exit 0.
# ⚠ Exit 0 is not enough: at 6.7.18 `cyrius audit` exits 0 even when the bench fails
# to compile (2.15.0 hit it). Read its bench section for `=== 1 passed, 0 failed ===`.
# `cyrius bench` and `cyrius test` on their own do exit 1 on a compile error.
cyrius audit
```

⚠ **Counting assertions:** `cyrius test tests/` prints a per-suite line *and* a
suite-level line, so a naive `grep -E '^[0-9]+ passed'` sum overcounts by the
number of suites. Require the parenthesised total:

```bash
cyrius test tests/ | grep -E '^[0-9]+ passed, [0-9]+ failed \([0-9]+ total\)$' \
  | awk '{s+=$1} END {print s}'
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
  `err_out` pointer, or NaN from a scalar function, and never a plausible-looking
  number for an input the function cannot answer. `err_out` is a precondition (it
  must point to 8 writable bytes and is written first); a null handle, a null output
  buffer or a null input buffer is refused.
- Physics functions document their units.
- **Correctness comes from the literature, not from the Rust original.** Encode `f64`
  constants as exact ratios (numerator/denominator < 2^53) or IEEE-754 hex bit
  patterns, and pin a value against an independent reference (a closed form, a
  high-precision computation, a published table), never against the code's own
  output. The modules ported from Rust keep its arithmetic where it was right — but
  2.2.6 found six wrong formulas, five copied character for character from it, so
  "it matches Rust" is never evidence. The modules added since 2.6.0 (`wave_beam`,
  `ray_grin`, `wave_doe`, `ray_prt`, `wave_focus`, `lens_aberration`,
  `spectral_fluor`, and more) are not ports at all. The Rust source is retrievable
  with `git checkout 2.2.2 -- rust-old/`, or `git show 1.2.0:src/<file>.rs`.

## Testing Requirements

- All public API changes need unit tests in a `tests/*.tcyr` suite.
- Performance-sensitive changes: **never claim a speedup without benchmarks**, and
  never from two runs at different times — this kind of host moves 20% between
  identical runs. Put the old and new code in one bench binary and compare them there
  (with an A/A pair), or interleave runs of both trees and compare per-row medians.
  `./scripts/bench-history.sh` records the history row.
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
