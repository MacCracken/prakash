#!/usr/bin/env bash
# Prove dist/prakash.cyr links the way a CONSUMER links it: only the stdlib
# leaves dist/prakash.deps names, and NO hisab.
#
# ⚠ WHY THIS EXISTS. README.md promises that without hisab "the bundle still
# builds and every other entry point works; you get one `undefined function
# 'num_fft'` warning". ranga depends on exactly that — it vendors dist/prakash.cyr
# and deliberately omits hisab, whose premultiply_alpha / linear_to_srgb collide
# with its own — and until 2.5.1 nothing checked it: CI's examples run inside
# prakash's own manifest, which declares hisab, so they would pass with the
# promise broken. This compiles a consumer-shaped program straight through cycc
# — no `cyrius build`, so no dependency resolution against prakash's manifest —
# with explicit includes of the declared leaves, over the surface ranga calls and
# re-exports, and asserts:
#   1. it compiles, and the ONLY undefined function is num_fft — so hisab really
#      was absent (a build that silently pulled hisab in would prove nothing);
#   2. it runs, with ranga's own pinned values (D65 CCT 6503.46 K, etc.) and the
#      Xyz field layout that ranga's tests pin with raw load64;
#   3. calling one of the four FFT entry points WITHOUT hisab refuses to build,
#      and num_fft is the reason — not some other missing name.
#
# Usage:  scripts/check-consumer-link.sh [--root DIR] [path/to/cycc ...]
#   --root DIR  take the stdlib leaves from DIR/lib instead of prakash's lib/,
#               compiling from DIR as that consumer would. With
#               --root ~/Repos/ranga ~/.cyrius/versions/6.6.2/bin/cycc this is
#               ranga's own toolchain and its own vendored stdlib.
#   cycc ...    further compilers to run checks 1-2 with, besides the pinned one.
# CI runs it bare: the pinned cycc, prakash's own lib/.
set -u
cd "$(dirname "$0")/.." || exit 2
PRAKASH=$(pwd)

ROOT=$PRAKASH
if [ "${1:-}" = "--root" ]; then
    [ -n "${2:-}" ] || { echo "usage: $0 [--root DIR] [cycc ...]"; exit 2; }
    ROOT=$(cd "$2" && pwd) || exit 2
    shift 2
fi
CC=$(cyrius which) || { echo "::error::cyrius which failed"; exit 2; }
compilers=("$CC" "$@")

T=$(mktemp -d) || exit 2
trap 'rm -rf "$T"' EXIT
fail=0

leaves=$(grep -vE '^[[:space:]]*(#|$)' dist/prakash.deps)
for l in $leaves; do
    if [ ! -f "$ROOT/lib/$l.cyr" ]; then
        echo "::error::dist/prakash.deps names '$l' but $ROOT/lib/$l.cyr is missing (run cyrius deps)"
        exit 1
    fi
done
if [ "$ROOT" = "$PRAKASH" ]; then
    bundle='dist/prakash.cyr'
else
    bundle="$PRAKASH/dist/prakash.cyr"
fi

fixture() {  # $1 = extra statement placed before the verdict
    for l in $leaves; do echo "include \"lib/$l.cyr\""; done
    echo "include \"$bundle\""
    # Helpers carry a _ccl_ prefix: Cyrius globals are last-one-wins, and the
    # bundle has its own `_r` — a fixture must not quietly replace bundle code.
    cat <<'EOF'
fn _ccl_say(s) { syscall(1, 1, s, strlen(s)); return 0; }
fn _ccl_r(a, b) { return f64_div(f64_from(a), f64_from(b)); }
fn _ccl_off(a, b, tol) { return f64_gt(f64_abs(f64_sub(a, b)), tol); }
fn main() {
    alloc_init();
    var bad = 0;
    # The surface ranga's src/spectral.cyr and tests/spectral.tcyr call.
    var w = xyz_d65_white();
    var s = f64_add(f64_add(Xyz_x(w), Xyz_y(w)), Xyz_z(w));
    var cct = cct_from_xy(f64_div(Xyz_x(w), s), f64_div(Xyz_y(w), s));
    if (_ccl_off(cct, _ccl_r(650346, 100), _ccl_r(1, 50)) == 1) { _ccl_say("D65 CCT is not 6503.46 K\n"); bad = bad + 1; }
    var w50 = xyz_d50_white();
    var s50 = f64_add(f64_add(Xyz_x(w50), Xyz_y(w50)), Xyz_z(w50));
    var cct50 = cct_from_xy(f64_div(Xyz_x(w50), s50), f64_div(Xyz_y(w50), s50));
    if (_ccl_off(cct50, _ccl_r(500212, 100), _ccl_r(1, 50)) == 1) { _ccl_say("D50 CCT is not 5002.12 K\n"); bad = bad + 1; }
    # ranga/tests/spectral.tcyr pins this layout with raw load64 at 0/8/16.
    var c = xyz_new(_ccl_r(1, 4), _ccl_r(1, 2), _ccl_r(3, 4));
    if (load64(c) != Xyz_x(c)) { _ccl_say("Xyz.x is not at offset 0\n"); bad = bad + 1; }
    if (load64(c + 8) != Xyz_y(c)) { _ccl_say("Xyz.y is not at offset 8\n"); bad = bad + 1; }
    if (load64(c + 16) != Xyz_z(c)) { _ccl_say("Xyz.z is not at offset 16\n"); bad = bad + 1; }
    if (_ccl_off(Xyz_y(cie_cmf_at(f64_from(555))), F64_ONE, _ccl_r(1, 100)) == 1) { _ccl_say("ybar(555) is not 1\n"); bad = bad + 1; }
    var bb = spd_blackbody(f64_from(6500));
    if (Spd_len(bb) != 81) { _ccl_say("blackbody SPD is not 81 samples\n"); bad = bad + 1; }
    if (Spd_start_nm(bb) != f64_from(380)) { _ccl_say("blackbody SPD does not start at 380\n"); bad = bad + 1; }
    if (Spd_step_nm(bb) != f64_from(5)) { _ccl_say("blackbody SPD step is not 5\n"); bad = bad + 1; }
    if (spec_visible_min_nm() != f64_from(380)) { _ccl_say("visible min is not 380\n"); bad = bad + 1; }
    if (spec_visible_max_nm() != f64_from(780)) { _ccl_say("visible max is not 780\n"); bad = bad + 1; }
    var d = spd_to_xyz(illuminant_d65());
    var ds = f64_add(f64_add(Xyz_x(d), Xyz_y(d)), Xyz_z(d));
    var dcct = cct_from_xy(f64_div(Xyz_x(d), ds), f64_div(Xyz_y(d), ds));
    if (_ccl_off(dcct, f64_from(6504), f64_from(5)) == 1) { _ccl_say("D65 SPD does not integrate to ~6504 K\n"); bad = bad + 1; }
    # The names ranga-spectral re-exports without calling (src/spectral.cyr:9-17):
    # linking them proves they exist and need nothing beyond the declared leaves.
    if (&cie_1931_table == 0) { bad = bad + 1; }
    if (&color_rendering_index == 0) { bad = bad + 1; }
    if (&color_temperature_to_rgb == 0) { bad = bad + 1; }
    if (&illuminant_a == 0) { bad = bad + 1; }
    if (&illuminant_d50 == 0) { bad = bad + 1; }
    if (&illuminant_f2 == 0) { bad = bad + 1; }
    if (&illuminant_f11 == 0) { bad = bad + 1; }
    if (&linear_to_srgb_gamma == 0) { bad = bad + 1; }
    if (&planck_radiance == 0) { bad = bad + 1; }
    if (&srgb_gamma_to_linear == 0) { bad = bad + 1; }
    if (&wavelength_to_rgb == 0) { bad = bad + 1; }
    if (&wien_peak == 0) { bad = bad + 1; }
EOF
    echo "$1"
    cat <<'EOF'
    if (bad == 0) { _ccl_say("consumer surface ok\n"); }
    return bad;
}
var rc = main();
syscall(60, rc);
EOF
}

compile() {  # $1 = cycc, $2 = source, $3 = output, $4 = log; returns cycc's status
    local rc
    (cd "$ROOT" && CYRIUS_ALLOW_ABSOLUTE_INCLUDES=1 CYRIUS_NO_WARN_PIN_DRIFT=1 "$1" < "$2" > "$3" 2> "$4")
    rc=$?
    chmod +x "$3"
    return $rc
}

undefined_set() {  # distinct undefined-function names in a log, one per line
    grep -oE "undefined function '[A-Za-z0-9_]+'" "$1" | sort -u
}

# 1-2. The consumer surface, hisab absent, under each compiler.
fixture "" > "$T/surface.cyr"
for cc in "${compilers[@]}"; do
    if ! compile "$cc" "$T/surface.cyr" "$T/surface" "$T/surface.log"; then
        echo "::error::$cc cannot compile the core bundle without hisab"
        grep -iE 'error' "$T/surface.log" | head -5
        fail=1
        continue
    fi
    undef=$(undefined_set "$T/surface.log")
    if [ "$undef" != "undefined function 'num_fft'" ]; then
        echo "::error::$cc: expected num_fft as the ONLY undefined function (proof hisab was absent); got:"
        echo "${undef:-<none — hisab was linked, so this check proved nothing>}"
        fail=1
    fi
    if ! out=$("$T/surface"); then
        echo "::error::$cc: the consumer surface returned wrong values without hisab:"
        echo "$out"
        fail=1
    else
        echo "$cc [$ROOT/lib]: $out"
    fi
done

# 3. An FFT entry point without hisab must refuse to build, because of num_fft.
fixture "    if (psf_diffraction_limited(16, f64_from(1)) == 0) { bad = bad + 1; }" > "$T/fft.cyr"
if compile "$CC" "$T/fft.cyr" "$T/fft" "$T/fft.log"; then
    echo "::error::calling psf_diffraction_limited without hisab BUILT — a call to an undefined num_fft was linked"
    fail=1
elif ! grep -q "reachable undefined" "$T/fft.log"; then
    echo "::error::the FFT-without-hisab build failed, but not for the expected reason:"
    grep -iE 'error' "$T/fft.log" | head -5
    fail=1
elif [ "$(undefined_set "$T/fft.log")" != "undefined function 'num_fft'" ]; then
    echo "::error::the FFT-without-hisab build was refused, but not (only) for num_fft:"
    undefined_set "$T/fft.log"
    fail=1
else
    echo "FFT surface without hisab: refused at build time for num_fft, as documented"
fi

exit $fail
