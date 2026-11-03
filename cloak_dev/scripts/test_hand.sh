#!/bin/bash
# Run the hand fixtures (cloak_dev/tests/data/hand/MANIFEST.tsv) on a binary.
#
# Usage (Debian terminal):
#   bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/test_hand.sh [options] [BINARY]
#     BINARY            default ~/muscle-build/dev/bin/muscle
#     --data DIR        hand fixture directory (default cloak_dev/tests/data/hand)
#     --predict COLUMN  v010 or v020: instead of requiring every case to pass, require each case to
#                       pass or fail as that MANIFEST column predicts (the recorded red step)
#
# Every case runs in a fresh copy of the fixture directory, capped at 2 GB of virtual memory
# (CLOAK_VMEM_KB overrides, in KiB or "unlimited", for sanitizer builds), 100 MB of file size and
# 20 s. A case passes if the exit status is exactly the expected one and:
#   exit 0: the output equals the expected file byte for byte;
#   exit 1: a sentinel file at the output path is unchanged (nothing was written there).
# In both cases no file may appear in the working directory. Stderr is not examined.
# Outputs go to a new directory under ~/muscle-runs/test_hand/. Exit status 0 only if every case
# passes (or, with --predict, matches its prediction).
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
DATA="$HERE/../tests/data/hand"
BIN="$HOME/muscle-build/dev/bin/muscle"
PRED=""
while [ $# -gt 0 ]; do
  case "$1" in
    --data) DATA=$(cd "$2" && pwd); shift ;;
    --predict) PRED=$2; shift ;;
    -*) echo "unknown option $1" >&2; exit 2 ;;
    *) BIN=$1 ;;
  esac
  shift
done
[ -x "$BIN" ] || { echo "binary not found: $BIN" >&2; exit 2; }
[ -f "$DATA/MANIFEST.tsv" ] || { echo "no MANIFEST.tsv in $DATA" >&2; exit 2; }
case "$PRED" in ''|v010|v020) ;; *) echo "--predict takes v010 or v020" >&2; exit 2 ;; esac
if [ "$DATA" = "$(cd "$HERE/../tests/data/hand" 2>/dev/null && pwd)" ]; then
  bash "$HERE/check_fixtures.sh" || exit 2
fi
OUT="$HOME/muscle-runs/test_hand/$(date +%Y%m%d-%H%M%S)"
mkdir -p "$OUT/out" "$OUT/log"

echo "binary:  $BIN ($(sha256sum "$BIN" | cut -c1-12))"
echo "data:    $DATA"
echo "outputs: $OUT"
pass=0; fail=0; total=0; mismatch=0
while IFS=$'\t' read -r id args ex exp v010 v020 what || [ -n "$id" ]; do
  id=${id%$'\r'}; what=${what%$'\r'}
  [ "$id" = "id" ] && continue
  total=$((total+1))
  run="$OUT/run/$id"; mkdir -p "$run"; cp -r "$DATA/." "$run/"
  before=$(cd "$run" && find . -type f | sort)
  out="$OUT/out/$id.fa"
  [ "$ex" = 1 ] && printf 'SENTINEL\n' > "$out"
  # shellcheck disable=SC2086
  ( cd "$run" && ulimit -v "${CLOAK_VMEM_KB:-2000000}" && ulimit -f 100000 && timeout 20 "$BIN" ${args//\{out\}/$out} \
      > "$OUT/log/$id.log" 2>&1 )
  rc=$?
  after=$(cd "$run" && find . -type f | sort)
  why=""
  if [ "$rc" != "$ex" ]; then why="exit $rc, expected $ex"
  elif [ "$before" != "$after" ]; then why="new files in the working directory"
  elif [ "$ex" = 0 ] && ! cmp -s "$out" "$DATA/$exp"; then why="output differs from $exp"
  elif [ "$ex" = 1 ] && [ "$(cat "$out" 2>/dev/null)" != "SENTINEL" ]; then why="output path written"
  fi
  rm -rf "$run"
  if [ -z "$why" ]; then res=pass; pass=$((pass+1)); else res=fail; fail=$((fail+1)); fi
  if [ -n "$PRED" ]; then
    want=$v010; [ "$PRED" = v020 ] && want=$v020
    if [ "$res" = "$want" ]; then tag="as predicted"; else tag="NOT AS PREDICTED"; mismatch=$((mismatch+1)); fi
    printf "%-5s %-6s %-17s %s%s\n" "$res" "$id" "$tag" "$what" "${why:+ ($why)}"
  else
    printf "%-5s %-6s %s%s\n" "$res" "$id" "$what" "${why:+ ($why)}"
  fi
done < "$DATA/MANIFEST.tsv"
rmdir "$OUT/run" 2>/dev/null
lines=$(grep -c . "$DATA/MANIFEST.tsv")
[ "$total" -eq $((lines-1)) ] || { echo "read $total cases but MANIFEST has $((lines-1))" >&2; exit 2; }
if [ -n "$PRED" ]; then
  echo "cases $total: passed $pass, failed $fail; not as predicted ($PRED) $mismatch" | tee "$OUT/SUMMARY.txt"
  [ "$mismatch" -eq 0 ]
else
  echo "cases $total: passed $pass, failed $fail" | tee "$OUT/SUMMARY.txt"
  [ "$fail" -eq 0 ]
fi
