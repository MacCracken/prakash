#!/usr/bin/env bash
# Rewrite dist/*.deps to state the truth for prakash's TWO bundles.
#
# WHY THIS EXISTS
# ---------------
# `cyrius distlib` writes a `dist/<lib>.deps` sidecar listing the stdlib folds a
# bundle needs in scope. It is published metadata describing the artifact.
#
# MEASURED SCOPE (re-measured at cyrius 6.5.33 for 2.2.4; do not overstate
# this): the sidecar is NOT an enforced contract. A consumer declaring fewer
# folds than it lists still resolves cleanly (exit 0), a bogus fold name in it
# raises no error, and in the git-dep flow the vendored set comes from prakash's
# own cyrius.cyml. We correct the sidecar because a published file should be
# true, not because a build breaks without it.
#
# ⛔ 2.10.1: THE REASON THIS SCRIPT EXISTS CHANGED UNDER IT. Through cyrius 6.6.17
# the generator emitted `[deps] stdlib` unioned with an include-scan for the base
# bundle and a pruned inference for profiles, which was backwards for prakash: the
# base (dist/prakash.cyr, math-only) is the NARROW bundle and the profile
# (dist/prakash-ai.cyr, + TLS) the WIDE one. 6.6.18 threw that away. A sidecar is
# now whatever a compile-verify fixpoint proves the bundle needs, and the bundle
# opens with a matching `# Requires` include block. That fixed the inversion and
# opened a new hole, because the verify runs with the package's NAMED deps in scope:
#
#   [deps.hisab]'s own requires block supplies math and ganita, and its
#   [deps.sakshi] supplies sakshi, so the generator derives
#     dist/prakash.deps     = bayan
#     dist/prakash-ai.deps  = sandhi bayan
#   "compile-verified", and FALSE for a consumer without hisab, which the README
#   promises works and ranga depends on. Measured at 6.6.20 and 6.7.9: the bundle with
#   only those leaves is ~800 errors (F64_ONE, ganita_*, sakshi_*), and
#   scripts/check-consumer-link.sh fails on the generated sidecar and passes on
#   this script's. With [deps.hisab] deleted from a scratch manifest the same
#   generator derives `math sakshi ganita string alloc vec fnptr result bayan str`
#   and warns on num_fft, so the verify is right about the question it is asked;
#   "does this build WITHOUT the optional dep" is not that question.
#
# ⚠ Two consequences of overriding it, both expected:
#   - `cyrius distlib --check` reads STALE after this script, because it compares
#     against the generated sidecar. CI runs distlib, this script and
#     `git diff dist/` instead; do the same.
#   - the `# Requires` block in each bundle still names only the generated leaves,
#     so `include "dist/prakash.cyr"` alone does not compile. Consumers declare
#     this sidecar's leaves, as the README says.
#
# Upstream issue (fix shipped v6.5.10, prakash's inverted layout not covered):
#   cyrius docs/development/issues/archived/2026-08-07-distlib-deps-sidecar-under-reports.md
# setu set the precedent of a lib owning dist/ and correcting its own sidecar.
#
# USAGE
#   ./scripts/sync-deps-sidecar.sh            # rewrite both sidecars
#   ./scripts/sync-deps-sidecar.sh --check    # exit 1 if either has drifted
#
# Run after every `cyrius distlib` / `cyrius distlib ai`.
set -euo pipefail

MANIFEST="cyrius.cyml"
CORE_DEPS="dist/prakash.deps"
AI_DEPS="dist/prakash-ai.deps"

# Folds pulled in ONLY by src/ai.cyr's sandhi HTTP client. These are exactly what
# separates the two bundles; everything else in `[deps] stdlib` is core-safe.
# Verified against the built bundles: dist/prakash.cyr contains zero references
# to any of them (see the ai-isolation check in CI).
AI_ONLY="net http tls async random fdlopen dynlib chrono sandhi"

# Parse the `[deps] stdlib = [ ... ]` array out of the manifest, in order.
read_declared() {
    awk '
        /^stdlib[[:space:]]*=[[:space:]]*\[/ { inarr = 1 }
        inarr {
            line = $0
            while (match(line, /"[^"]+"/)) {
                item = substr(line, RSTART + 1, RLENGTH - 2)
                if (item != "stdlib") print item
                line = substr(line, RSTART + RLENGTH)
            }
            if (index($0, "]") > 0 && !/^stdlib[[:space:]]*=[[:space:]]*\[[[:space:]]*$/) inarr = 0
        }
    ' "$MANIFEST"
}

is_ai_only() {
    local fold="$1" a
    for a in $AI_ONLY; do [ "$a" = "$fold" ] && return 0; done
    return 1
}

emit() {
    # $1 = "core" | "ai"
    printf '# cyrius dep sidecar - stdlib leaves this fold requires in scope.\n'
    printf '# Maintained by scripts/sync-deps-sidecar.sh (NOT by `cyrius distlib` —\n'
    printf '# see that script for why prakash overrides the generated sidecar).\n'
    local fold
    while IFS= read -r fold; do
        [ -z "$fold" ] && continue
        if [ "$1" = "core" ] && is_ai_only "$fold"; then continue; fi
        printf '%s\n' "$fold"
    done < <(read_declared)
}

declared_count=$(read_declared | wc -l | tr -d ' ')
if [ "$declared_count" -lt 5 ]; then
    echo "ERROR: parsed only $declared_count stdlib folds from $MANIFEST — refusing to write." >&2
    exit 1
fi

if [ "${1:-}" = "--check" ]; then
    rc=0
    for pair in "core:$CORE_DEPS" "ai:$AI_DEPS"; do
        profile="${pair%%:*}"; path="${pair#*:}"
        if ! diff -q <(emit "$profile") "$path" >/dev/null 2>&1; then
            echo "::error file=$path::deps sidecar drift — run ./scripts/sync-deps-sidecar.sh"
            diff <(emit "$profile") "$path" || true
            rc=1
        fi
    done
    [ $rc -eq 0 ] && echo "deps sidecars up to date (core $(( declared_count - $(echo $AI_ONLY | wc -w) )) folds, ai $declared_count)"
    exit $rc
fi

emit core > "$CORE_DEPS"
emit ai   > "$AI_DEPS"
echo "wrote $CORE_DEPS ($(grep -vc '^#' "$CORE_DEPS") folds — TLS stack excluded)"
echo "wrote $AI_DEPS ($(grep -vc '^#' "$AI_DEPS") folds — full stack)"
