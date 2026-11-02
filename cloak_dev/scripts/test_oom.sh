#!/bin/bash
# E16, allocation failure: run the largest fixture under a series of address-space caps
# (cloak_dev/tests/README.md).
#
# Usage (Debian terminal):
#   bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/test_oom.sh [BINARY]
#     BINARY   default ~/muscle-build/dev/bin/muscle
#
# The caps run from 16 MiB to 2 GiB in steps of sqrt(2), so they straddle the binary's need whatever
# its memory use; none is tuned to a binary. (MUSCLE itself starts at 8 MiB and segfaults at 4 MiB:
# `-version` under `ulimit -v`, 2026-11-02.) Every run must either exit 1 with the sentinel at the
# output path unchanged, or exit 0 with the expected output byte for byte, and at least one run of
# each kind is required. Stderr is not examined. Not run on the byte oracle (scripts/oracle_bin.sh),
# whose capped Python interpreter fails at startup with other exit statuses.
# Outputs go to a new directory under ~/muscle-runs/test_oom/. Exit status 0 only if every run passes.
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
DATA="$HERE/../tests/data"
BIN=${1:-$HOME/muscle-build/dev/bin/muscle}
[ -x "$BIN" ] || { echo "binary not found: $BIN" >&2; exit 2; }
bash "$HERE/check_fixtures.sh" || exit 2
IN="$DATA/PF00353full_n1226_s5.efa"; EXP="$DATA/expected/PF00353full_n1226_s5.m2.v030.fa"
CAPS="16384 23170 32768 46341 65536 92682 131072 185364 262144 370728 524288 741455 1048576 1482910 2097152"  # KiB, 16384 x 2^(i/2) rounded, i = 0..14
OUT="$HOME/muscle-runs/test_oom/$(date +%Y%m%d-%H%M%S)"
mkdir -p "$OUT/out" "$OUT/log"

echo "binary:  $BIN ($(sha256sum "$BIN" | cut -c1-12))"
echo "outputs: $OUT"
pass=0; fail=0; n0=0; n1=0
for cap in $CAPS; do
  out="$OUT/out/cap$cap.fa"; run="$OUT/run/cap$cap"; mkdir -p "$run"
  printf 'SENTINEL\n' > "$out"
  ( cd "$run" && ulimit -v "$cap" && timeout 120 "$BIN" -cloak "$IN" -output "$out" ) \
    > "$OUT/log/cap$cap.log" 2>&1
  rc=$?
  why=""
  if [ -n "$(ls -A "$run")" ]; then why="new files in the working directory"
  elif [ "$rc" = 1 ]; then
    [ "$(cat "$out" 2>/dev/null)" = "SENTINEL" ] || why="exit 1 but the output path was written"
  elif [ "$rc" = 0 ]; then
    cmp -s "$out" "$EXP" || why="exit 0 but the output differs from the expected output"
  else why="exit $rc"
  fi
  rm -rf "$run"
  if [ -z "$why" ]; then
    pass=$((pass+1)); [ "$rc" = 0 ] && n0=$((n0+1)) || n1=$((n1+1))
    printf "pass  cap %7s KiB  exit %s\n" "$cap" "$rc"
  else
    fail=$((fail+1)); printf "fail  cap %7s KiB  %s\n" "$cap" "$why"
  fi
done
rmdir "$OUT/run" 2>/dev/null
straddle=""
[ "$n0" -ge 1 ] && [ "$n1" -ge 1 ] || straddle="; the caps do not straddle the need (exit 0: $n0, exit 1: $n1)"
echo "caps $((pass+fail)): passed $pass (exit 0: $n0, exit 1: $n1), failed $fail$straddle" \
  | tee "$OUT/SUMMARY.txt"
[ "$fail" -eq 0 ] && [ -z "$straddle" ]
