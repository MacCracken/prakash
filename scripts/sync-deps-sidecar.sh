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
# own cyrius.cyml. The upstream issue below reports a hard-error behaviour that
# does NOT reproduce here. We correct the sidecar because a published file
# should be true, not because a build breaks without it.
#
# Since cyrius 6.5.10 the generator emits, for the BASE bundle, the declared
# `[deps] stdlib` unioned with an include-scan; PROFILE bundles keep a pruned
# inference "so a narrow bundle cannot over-report". That rule assumes the base
# bundle is the widest one and profiles are subsets of it. prakash is inverted:
#
#   [lib]     -> dist/prakash.cyr     math-only core, NO TLS   (the NARROW one)
#   [lib.ai]  -> dist/prakash-ai.cyr  core + AI client, + TLS  (the WIDE one)
#
# so the generator gets both backwards for us:
#   - dist/prakash.deps OVER-reports  — it inherits the whole sandhi/TLS stack
#     from `[deps] stdlib` (which must list it, because src/ai.cyr needs it),
#     even though the core bundle references none of it. That misdescribes the
#     math-only bundle to anyone reading it (and to any future tooling that
#     does treat it as authoritative).
#   - dist/prakash-ai.deps UNDER-reports — the pruned inference omits most of
#     what sandhi transitively needs.
#     ⚠ RE-MEASURED AT 6.5.33: through 6.5.20 this yielded literally
#     `syscalls io`. It now yields ten folds — `string alloc str vec math ganita
#     tagged bayan sandhi sakshi` — still missing `net http tls async random
#     fdlopen dynlib chrono` (the stack sandhi actually needs) as well as
#     `syscalls io fmt args assert result simd fnptr bench callback`. The
#     shortfall shrank; it did not close. The core over-reporting above
#     reproduces unchanged.
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
