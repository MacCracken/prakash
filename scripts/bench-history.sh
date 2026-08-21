#!/usr/bin/env bash
# Run prakash's Cyrius benchmark suite, append results to a CSV history, and
# regenerate benchmarks.md with a 3-point trend table (baseline → mid → current).
#
# Usage:
#   ./scripts/bench-history.sh                     # default CSV: bench-history.csv
#   ./scripts/bench-history.sh results.csv         # custom CSV
#   ./scripts/bench-history.sh "" tests/foo.bcyr   # custom suite
#
# `cyrius bench` emits one line per benchmark:
#   <name>: <avg><unit> avg (min=<..> max=<..>) [<n> iters]
# We record the reported avg (per iter), normalized to nanoseconds.
#
# ⚠ THE INSTRUMENT CHANGED UNDER THIS SCRIPT AND THIS HEADER WAS LEFT ASSERTING
# THE OLD BEHAVIOUR. It used to say the harness's fixed per-call overhead
# "(~0.5-0.9us) is included ... the min in the raw output is the better absolute
# proxy". That has been false since cyrius 6.5.19, which taught the harness to
# CALIBRATE one clock read on the host and SUBTRACT it from every sample. The
# harness prints the figure it measured:
#   [timer floor 1.336us per clock read, measured; subtracted from every sample]
# and this very script has been echoing that banner while the header above it
# denied it.
#
# So the regime is now DERIVED, not declared: each row records `regime` and
# `floor_ns` read from whether the harness printed that line. A constant in a
# comment is exactly what went stale; a derived column cannot.
#
# ⚠ What `regime` can and cannot do, said plainly: it separates "floor
# subtracted" from "floor not subtracted". A future instrument change that kept
# printing the same banner would NOT flip it. This narrows the hole; it does not
# close it.
#
# Rows either side of that boundary are NOT comparable and every one of them
# says "avg", so the trend table below compares only within a single regime and
# reports how many rows it excluded rather than quietly dropping them.
set -euo pipefail

HISTORY_FILE="${1:-bench-history.csv}"
SUITE="${2:-tests/prakash.bcyr}"
BENCHMARKS_MD="benchmarks.md"
TIMESTAMP=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
COMMIT=$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")
BRANCH=$(git branch --show-current 2>/dev/null || echo "unknown")

if [ ! -f "$HISTORY_FILE" ]; then
    echo "timestamp,commit,branch,benchmark,estimate_ns,regime,floor_ns" > "$HISTORY_FILE"
fi

echo "╔══════════════════════════════════════════╗"
echo "║        prakash benchmark suite           ║"
echo "╠══════════════════════════════════════════╣"
echo "║  commit: $COMMIT"
echo "║  branch: $BRANCH"
echo "║  date:   $TIMESTAMP"
echo "║  suite:  $SUITE"
echo "╚══════════════════════════════════════════╝"
echo ""

normalize_to_ns() {
    local val="$1" unit="$2"
    case "$unit" in
        ps)    awk -v v="$val" 'BEGIN{printf "%.4f", v / 1000}' ;;
        ns)    echo "$val" ;;
        us|µs) awk -v v="$val" 'BEGIN{printf "%.4f", v * 1000}' ;;
        ms)    awk -v v="$val" 'BEGIN{printf "%.4f", v * 1000000}' ;;
        s)     awk -v v="$val" 'BEGIN{printf "%.4f", v * 1000000000}' ;;
        *)     echo "$val" ;;
    esac
}

BENCH_OUTPUT=$(cyrius bench "$SUITE" 2>&1 | sed 's/\x1b\[[0-9;]*m//g')
echo "$BENCH_OUTPUT"
echo ""

# Derive the measurement regime from the harness's own banner rather than from a
# hardcoded version check — see the header. The line looks like:
#   [timer floor 1.336us per clock read, measured; subtracted from every sample]
FLOOR_LINE=$(echo "$BENCH_OUTPUT" | grep -m1 -E 'timer floor .* subtracted from every sample' || true)
if [ -n "$FLOOR_LINE" ]; then
    REGIME="floor_subtracted"
    FLOOR_VAL=$(echo "$FLOOR_LINE" | grep -oE '[0-9]+(\.[0-9]+)?[[:space:]]*(ps|ns|µs|us|ms|s)' | head -n1)
    FLOOR_NUM=$(echo "$FLOOR_VAL" | grep -oE '[0-9]+(\.[0-9]+)?' | head -n1)
    FLOOR_UNIT=$(echo "$FLOOR_VAL" | grep -oE '(ps|ns|µs|us|ms|s)$')
    FLOOR_NS=$(normalize_to_ns "$FLOOR_NUM" "$FLOOR_UNIT")
else
    REGIME="floor_included"
    FLOOR_NS=""
fi
echo "regime: $REGIME (floor_ns=${FLOOR_NS:-n/a})"
echo ""


# Parse each "<name>: <avg><unit> avg (...)" line; the avg is the token before " avg".
while IFS= read -r line; do
    echo "$line" | grep -qE ':[[:space:]]*[0-9.]+[[:space:]]*(ps|ns|µs|us|ms|s)[[:space:]]+avg' || continue
    NAME=$(echo "$line" | sed -E 's/^[[:space:]]*//; s/:.*$//' | xargs)
    [ -z "$NAME" ] && continue
    TOKEN=$(echo "$line" | grep -oE '[0-9]+(\.[0-9]+)?[[:space:]]*(ps|ns|µs|us|ms|s)[[:space:]]+avg' | head -n1)
    VAL=$(echo "$TOKEN" | grep -oE '[0-9]+(\.[0-9]+)?' | head -n1)
    UNIT=$(echo "$TOKEN" | grep -oE '(ps|ns|µs|us|ms|s)[[:space:]]+avg' | grep -oE '(ps|ns|µs|us|ms|s)')
    NS=$(normalize_to_ns "$VAL" "$UNIT")
    echo "${TIMESTAMP},${COMMIT},${BRANCH},${NAME},${NS},${REGIME},${FLOOR_NS}" >> "$HISTORY_FILE"
done <<< "$BENCH_OUTPUT"

# 3-point trend table. Skip if Python missing — the CSV is the primary record.
command -v python3 >/dev/null 2>&1 || { echo "python3 absent — CSV updated, skipping benchmarks.md"; exit 0; }

python3 - "$HISTORY_FILE" "$BENCHMARKS_MD" "$COMMIT" <<'PYEOF'
import csv, sys
from collections import OrderedDict

history_file, md_file, cur_commit = sys.argv[1], sys.argv[2], sys.argv[3]
rows = list(csv.DictReader(open(history_file)))
if not rows:
    sys.exit(0)

# Rows either side of the timer-floor change are NOT comparable and both say
# "avg", so filter to a single regime before picking trend points. The regime is
# read from the row, which derived it from the harness banner (see the header).
# Rows predating the column read as "" and are classified once, below, by the
# only evidence left in the data: a floor-included run cannot have a sub-microsecond
# minimum, because the floor itself is ~1.3 us.
by_ts = OrderedDict()
for r in rows:
    by_ts.setdefault(r["timestamp"], []).append(r)

def regime_of(ts):
    rs = by_ts[ts]
    declared = {r.get("regime") or "" for r in rs} - {""}
    if declared:
        return sorted(declared)[0]
    vals = []
    for r in rs:
        try:
            vals.append(float(r["estimate_ns"]))
        except ValueError:
            pass
    if not vals:
        return "unknown"
    return "floor_included" if min(vals) > 800.0 else "floor_subtracted"

cur_regime = regime_of(list(by_ts)[-1])
all_stamps = list(by_ts)
stamps = [t for t in all_stamps if regime_of(t) == cur_regime]
excluded = [t for t in all_stamps if t not in stamps]
excluded_rows = sum(len(by_ts[t]) for t in excluded)
if excluded:
    print(f"trend: regime={cur_regime}; excluded {excluded_rows} rows "
          f"across {len(excluded)} run(s) measured in a different regime "
          f"({', '.join(sorted({regime_of(t) for t in excluded}))})")
rows = [r for r in rows if r["timestamp"] in stamps]

if len(stamps) >= 3:
    pick = [stamps[0], stamps[len(stamps)//2], stamps[-1]]
elif len(stamps) == 2:
    pick = [stamps[0], stamps[-1]]
else:
    pick = [stamps[0]]
seen = set(); pick = [t for t in pick if not (t in seen or seen.add(t))]

data, commits = {}, {}
for r in rows:
    if r["timestamp"] in pick:
        try:
            data.setdefault(r["benchmark"], {})[r["timestamp"]] = float(r["estimate_ns"])
        except ValueError:
            continue
        commits[r["timestamp"]] = r["commit"]

labels = []
for i, ts in enumerate(pick):
    if i == 0 and len(pick) > 1:
        labels.append(f"Baseline (`{commits[ts]}`)")
    elif i == len(pick) - 1:
        labels.append(f"Current (`{commits[ts]}`)")
    else:
        labels.append(f"Mid (`{commits[ts]}`)")

def fmt(ns):
    if ns is None: return "—"
    if ns >= 1e6: return f"{ns/1e6:.3f} ms"
    if ns >= 1e3: return f"{ns/1e3:.3f} µs"
    return f"{ns:.1f} ns"

lines = ["# prakash benchmarks", "",
         "Cyrius benchmark history for the ported optics modules. Generated by",
         "`scripts/bench-history.sh` from `tests/prakash.bcyr`. Values are the",
         "harness-reported avg per iter.", "",
         f"> All three points are measured in the **`{cur_regime}`** regime, so the Δ",
         "> column is like-for-like.", "",
         "> **⚠ This table deliberately does not span the timer-floor change.**",
         "> Cyrius 6.5.19 taught the bench harness to calibrate one clock read and",
         "> subtract it from every sample. Runs before that include the harness's",
         "> own ~1.32 µs call overhead in every number, which flattens all cheap",
         "> operations onto the floor — a row that \"improved\" by ~90% across that",
         "> boundary did not get faster, the overhead simply stopped being counted.",
         "> Both regimes label their numbers \"avg\", so they cannot be told apart",
         "> from the value alone; `bench-history.csv` carries a derived `regime`",
         "> column and this table filters on it.",
         (f"> **{excluded_rows} row(s) across {len(excluded)} earlier run(s) are excluded** "
          "from the table for that reason — they remain in the CSV.") if excluded
             else "> No rows were excluded — the whole history is in one regime.", "",
         "| Benchmark | " + " | ".join(labels) + " | Δ |", "|" + "---|" * (len(labels) + 2)]
for name in sorted(data):
    cells = [fmt(data[name].get(ts)) for ts in pick]
    first, last = data[name].get(pick[0]), data[name].get(pick[-1])
    delta = "—"
    if first and last and len(pick) > 1:
        pct = (last - first) / first * 100.0
        delta = f"{pct:+.1f}%"
    lines.append(f"| {name} | " + " | ".join(cells) + f" | {delta} |")
open(md_file, "w").write("\n".join(lines) + "\n")
print(f"benchmarks.md updated ({len(data)} benchmarks, {len(pick)} points)")
PYEOF
