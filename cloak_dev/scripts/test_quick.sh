#!/bin/bash
# Byte-compare CLOAK output with the expected outputs in cloak_dev/tests/data (23 cases).
#
# Usage (Debian terminal):
#   bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/test_quick.sh [options] [BINARY]
#     BINARY                default ~/muscle-build/dev/bin/muscle (from build.sh)
#     --check-baseline      also run the frozen baseline on the cases it finishes quickly, to confirm
#                           the expected outputs reproduce on this machine (about 1 min)
#     --check-baseline-slow also the two slow cases (about 4 min, up to 1.3 GB)
#
# scripts/check_fixtures.sh runs first (decision 009). Every run is capped at 4 GB of virtual
# memory (CLOAK_VMEM_KB overrides, in KiB or "unlimited", for sanitizer builds) and 120 s (600 s for the baseline on the slow cases). The cases with a baseline column of
# fast or slow keep every row (checked 2026-11-01), so decision 003's relation reduces to byte
# identity with the expected output for the baseline too. Outputs go to a new directory under
# ~/muscle-runs/test_quick/. Exit status 0 only if every case read from MANIFEST.tsv is compared and
# byte-identical.
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
DATA="$HERE/../tests/data"
BASE="$HOME/muscle-baseline/muscle-f3559af"
BIN="$HOME/muscle-build/dev/bin/muscle"
CHECK=none
for a in "$@"; do
  case "$a" in
    --check-baseline) CHECK=fast ;;
    --check-baseline-slow) CHECK=slow ;;
    -*) echo "unknown option $a" >&2; exit 2 ;;
    *) BIN="$a" ;;
  esac
done
[ -x "$BIN" ] || { echo "binary not found: $BIN" >&2; exit 2; }
bash "$HERE/check_fixtures.sh" || exit 2
OUT="$HOME/muscle-runs/test_quick/$(date +%Y%m%d-%H%M%S)"
mkdir -p "$OUT"

run_case() {  # label binary efa mincol expected [timeout] -> one result line, 0 if identical
  local label=$1 bin=$2 efa=$3 mc=$4 exp=$5 cap=${6:-120} name t0 t1 rc
  name="$label.$(basename "$efa" .efa).m$mc"
  t0=$(date +%s.%N)
  ( ulimit -v "${CLOAK_VMEM_KB:-4000000}"; timeout "$cap" "$bin" -cloak "$DATA/$efa" -mincol "$mc" -output "$OUT/$name.fa" > "$OUT/$name.log" 2>&1 )
  rc=$?
  t1=$(date +%s.%N)
  local secs; secs=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.2f", b-a}')
  if [ $rc -ne 0 ]; then printf "FAIL  %-9s %-26s mincol %s  exit %s (%s s)\n" "$label" "$efa" "$mc" "$rc" "$secs"; return 1; fi
  if cmp -s "$OUT/$name.fa" "$DATA/$exp"; then
    printf "ok    %-9s %-26s mincol %s  %6s s\n" "$label" "$efa" "$mc" "$secs"; return 0
  fi
  printf "FAIL  %-9s %-26s mincol %s  output differs from %s\n" "$label" "$efa" "$mc" "$exp"; return 1
}

echo "candidate: $BIN"
[ -f "$(dirname "$(dirname "$BIN")")/BUILD_INFO.txt" ] && sed 's/^/  /' "$(dirname "$(dirname "$BIN")")/BUILD_INFO.txt"
echo "outputs:   $OUT"
pass=0; fail=0; cases=0
while IFS=$'\t' read -r efa mc exp _src base || [ -n "$efa" ]; do
  base=${base%$'\r'}
  [ "$efa" = "ensemble" ] && continue
  cases=$((cases+1))
  if run_case candidate "$BIN" "$efa" "$mc" "$exp"; then pass=$((pass+1)); else fail=$((fail+1)); fi
  if { [ "$CHECK" = fast ] && [ "$base" = fast ]; } || { [ "$CHECK" = slow ] && [ "$base" != never ]; }; then
    cap=120; [ "$base" = slow ] && cap=600
    if run_case baseline "$BASE" "$efa" "$mc" "$exp" "$cap"; then pass=$((pass+1)); else fail=$((fail+1)); fi
  fi
done < "$DATA/MANIFEST.tsv"
want=$(($(grep -c . "$DATA/MANIFEST.tsv") - 1))
[ "$cases" -eq "$want" ] || { echo "read $cases cases, MANIFEST.tsv has $want" >&2; fail=$((fail+1)); }
echo "cases $cases: passed $pass, failed $fail" | tee "$OUT/SUMMARY.txt"
[ $fail -eq 0 ]
