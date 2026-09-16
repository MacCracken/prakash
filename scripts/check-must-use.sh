#!/usr/bin/env bash
# Guard against a function silently losing its #must_use.
#
# ⚠ WHY THIS EXISTS. `#must_use` binds to whatever declaration FOLLOWS it, so
# inserting a new function between an existing one and its `# doc` + `#must_use`
# block silently rebinds the attribute to the newcomer and leaves the original
# bare. It happened FOUR times in the 2.3.x cycle — 2.3.1
# (`atm_sky_radiance_single_scatter`), 2.3.4 (`rgb_to_json`), 2.3.6
# (`ray_snell_3d`), and once more caught in review — including twice AFTER the
# hazard was written into CLAUDE.md as a rule. A rule that is violated four times
# is not a control; this is.
#
# ⛔ NO OTHER GATE SEES IT. `cyrius lint` is clean either way. `cyrius audit`
# notices only if the DOC comment travelled with the attribute, and reports it as
# "undocumented", never as a lost attribute. The count alone is not enough either:
# add one function and remove another's marker and the total is unchanged — which
# is why this compares the SET OF NAMES, not the number.
#
# Usage:  scripts/check-must-use.sh [baseline-ref]      (default: HEAD)
#         scripts/check-must-use.sh --selftest
set -u

marked() {  # $1 = file contents on stdin; prints names whose previous line is #must_use
    awk '
        /^fn [A-Za-z_][A-Za-z0-9_]*[[:space:]]*\(/ {
            if (prev ~ /^[[:space:]]*#must_use[[:space:]]*$/) {
                name = $2; sub(/\(.*/, "", name); print name
            }
        }
        { prev = $0 }
    '
}

if [ "${1:-}" = "--selftest" ]; then
    tmp=$(mktemp); trap 'rm -f "$tmp"' EXIT
    printf '#must_use\nfn kept() { return 0; }\n\n# doc\n#must_use\nfn also() { return 0; }\n\nfn bare() { return 0; }\n' > "$tmp"
    got=$(marked < "$tmp" | sort | tr '\n' ' ')
    [ "$got" = "also kept " ] || { echo "SELFTEST FAILED: got '$got'"; exit 1; }
    # and the hazard itself: a function inserted into the gap must NOT inherit
    printf '# doc\n#must_use\nfn inserted() { return 0; }\nfn victim() { return 0; }\n' > "$tmp"
    got=$(marked < "$tmp" | sort | tr '\n' ' ')
    [ "$got" = "inserted " ] || { echo "SELFTEST FAILED (hazard): got '$got'"; exit 1; }
    echo "selftest ok"; exit 0
fi

BASE="${1:-HEAD}"
old=$(mktemp); new=$(mktemp); trap 'rm -f "$old" "$new"' EXIT

for f in $(git ls-files 'src/*.cyr'); do
    git show "$BASE:$f" 2>/dev/null | marked
done | sort -u > "$old"

for f in $(git ls-files 'src/*.cyr'); do
    marked < "$f"
done | sort -u > "$new"

lost=$(comm -23 "$old" "$new")
if [ -n "$lost" ]; then
    echo "ERROR: these functions had #must_use at $BASE and do not now:"
    echo "$lost" | sed 's/^/  /'
    echo
    echo "If a removal is intentional, say so in the CHANGELOG and re-baseline."
    echo "Otherwise a function was almost certainly inserted between one of them"
    echo "and its doc block — see the header of this script."
    exit 1
fi

echo "#must_use intact: $(wc -l < "$new") functions, 0 lost against $BASE"
