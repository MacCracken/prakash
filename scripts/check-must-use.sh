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

marked() {  # file contents on stdin; prints the names of #must_use-bearing fns
    # ⚠ The attribute binds ACROSS comments and blank lines, and this function got
    # that wrong when it shipped in 2.3.6. It matched only the immediately-preceding
    # line, so three functions whose `#must_use` sat above a doc comment —
    # pbr_fresnel_schlick, _ray_snell_3d_out, _pattern2d_new_uninit — were invisible
    # to the gate: it would not have noticed them losing the attribute, which is the
    # one thing it exists to catch. Verified against cycc in 2.4.0: a `#must_use`
    # separated from its fn by a comment, by several comments, or by a blank line
    # still warns on a discarded result, so all of those forms are live bindings and
    # must be tracked. Walk backwards over comments and blanks.
    awk '
        /^fn [A-Za-z_][A-Za-z0-9_]*[[:space:]]*\(/ {
            found = 0
            for (k = n; k >= 1; k--) {
                line = buf[k]
                if (line ~ /^[[:space:]]*#must_use[[:space:]]*$/) { found = 1; break }
                if (line ~ /^[[:space:]]*#/)  { continue }   # doc comment or other attr
                if (line ~ /^[[:space:]]*$/)  { continue }   # blank
                break                                        # any real code ends the run
            }
            if (found) { name = $2; sub(/\(.*/, "", name); print name }
        }
        { n++; buf[n] = $0 }
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
    # ⚠ The attribute binds ACROSS comments and blank lines — verified against cycc,
    # which warns on a discarded result in every one of these forms. This gate
    # missed all three until 2.4.0 and so could not have reported them going
    # missing. Each must be recognised.
    printf '#must_use\n# a doc comment\nfn after_comment() { return 0; }\n' > "$tmp"
    got=$(marked < "$tmp" | tr '\n' ' ')
    [ "$got" = "after_comment " ] || { echo "SELFTEST FAILED (comment): got '$got'"; exit 1; }
    printf '#must_use\n\nfn after_blank() { return 0; }\n' > "$tmp"
    got=$(marked < "$tmp" | tr '\n' ' ')
    [ "$got" = "after_blank " ] || { echo "SELFTEST FAILED (blank): got '$got'"; exit 1; }
    printf '#must_use\n# one\n# two\n\n# three\nfn after_many() { return 0; }\n' > "$tmp"
    got=$(marked < "$tmp" | tr '\n' ' ')
    [ "$got" = "after_many " ] || { echo "SELFTEST FAILED (many): got '$got'"; exit 1; }
    # a real statement between them DOES end the run — the attribute is not ours
    printf '#must_use\nvar x = 0;\nfn not_marked() { return 0; }\n' > "$tmp"
    got=$(marked < "$tmp" | tr '\n' ' ')
    [ "$got" = "" ] || { echo "SELFTEST FAILED (negative control): got '$got'"; exit 1; }
    echo "selftest ok"; exit 0
fi

BASE="${1:-HEAD}"
old=$(mktemp); new=$(mktemp); trap 'rm -f "$old" "$new"' EXIT

for f in $(git ls-files 'src/*.cyr'); do
    git show "$BASE:$f" 2>/dev/null | marked
done | sort -u > "$old"

# ⚠ The CURRENT side must include untracked-but-not-ignored files (2.6.0). With
# `git ls-files` alone a new module is invisible until it is staged, so the count
# printed below disagreed with CLAUDE.md's `grep -rhoE '^\s*#must_use' src/*.cyr`
# reconciliation (430 vs 461 for the new wave_beam.cyr), and a function MOVED into
# a new unstaged file would have been reported lost. The BASE side needs no change:
# an untracked file has nothing at the baseline to lose.
for f in $(git ls-files --cached --others --exclude-standard 'src/*.cyr'); do
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
