#!/bin/bash
# Error cases that need the shell rather than a hand MANIFEST row: failed writes, standard output
# and an empty -cloak value (cloak_dev/tests/README.md).
#
# Usage (Debian terminal):
#   bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/test_errors.sh [BINARY]
#     BINARY   default ~/muscle-build/dev/bin/muscle
#
# Every case runs in a fresh empty working directory, capped at 2 GB of virtual memory
# (CLOAK_VMEM_KB overrides, in KiB or "unlimited") and 20 s. A
# case passes if the exit status is exactly the expected one, the check below holds, and no file
# appears in the working directory. Stderr is not examined.
#   E20a  -output /dev/full, small output, so the error appears only at close: exit 1
#   E20b  -output - with standard output on /dev/full: exit 1
#   E20c  file-size limit below the output size, SIGXFSZ ignored, so a write fails midway:
#         exit 1 and no file at the output path
#   E21   -cloak '': exit 1, sentinel at the output path unchanged
#   H17   -output - redirected to a file: exit 0, the expected output byte for byte
# Outputs go to a new directory under ~/muscle-runs/test_errors/. Exit status 0 only if every case
# passes.
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
DATA="$HERE/../tests/data"
BIN=${1:-$HOME/muscle-build/dev/bin/muscle}
[ -x "$BIN" ] || { echo "binary not found: $BIN" >&2; exit 2; }
bash "$HERE/check_fixtures.sh" || exit 2
SMALL="$DATA/hand/H03.efa"; SMALL_EXP="$DATA/hand/expected/H03.fa"
BIG="$DATA/PF00745_n387_s5.efa"; BIG_EXP="$DATA/expected/PF00745_n387_s5.m2.fa"
LIMIT_KB=50
VMEM=${CLOAK_VMEM_KB:-2000000}   # unlimited for sanitizer builds (test_sanitize.sh)
[ "$(stat -c %s "$BIG_EXP")" -gt $((LIMIT_KB*1024)) ] ||
  { echo "E20c needs an expected output larger than $LIMIT_KB KiB" >&2; exit 2; }
OUT="$HOME/muscle-runs/test_errors/$(date +%Y%m%d-%H%M%S)"
mkdir -p "$OUT/out" "$OUT/log"

echo "binary:  $BIN ($(sha256sum "$BIN" | cut -c1-12))"
echo "outputs: $OUT"
pass=0; fail=0
start() { run="$OUT/run/$1"; mkdir -p "$run"; out="$OUT/out/$1.fa"; log="$OUT/log/$1.log"; }
finish() {  # ID WHAT EXPECTED_EXIT RC WHY (empty unless the case-specific check failed)
  local why=$5
  if [ "$4" != "$3" ]; then why="exit $4, expected $3"
  elif [ -n "$(ls -A "$run")" ]; then why="new files in the working directory"
  fi
  rm -rf "$run"
  if [ -z "$why" ]; then pass=$((pass+1)); printf "pass  %-5s %s\n" "$1" "$2"
  else fail=$((fail+1)); printf "fail  %-5s %s (%s)\n" "$1" "$2" "$why"; fi
}

start E20a
( cd "$run" && ulimit -v "$VMEM" && timeout 20 "$BIN" -cloak "$SMALL" -output /dev/full ) > "$log" 2>&1
finish E20a "-output /dev/full" 1 $? ""

start E20b
( cd "$run" && ulimit -v "$VMEM" && timeout 20 "$BIN" -cloak "$SMALL" -output - ) > /dev/full 2> "$log"
finish E20b "-output - on a full standard output" 1 $? ""

start E20c
( cd "$run" && trap '' XFSZ && ulimit -v "$VMEM" && ulimit -f $LIMIT_KB &&
  timeout 20 "$BIN" -cloak "$BIG" -output "$out" ) > "$log" 2>&1
rc=$?
why=""; [ -e "$out" ] && why="file left at the output path ($(stat -c %s "$out") bytes)"
finish E20c "write fails at a $LIMIT_KB KiB file-size limit" 1 $rc "$why"

start E21
printf 'SENTINEL\n' > "$out"
( cd "$run" && ulimit -v "$VMEM" && timeout 20 "$BIN" -cloak '' -output "$out" ) > "$log" 2>&1
rc=$?
why=""; [ "$(cat "$out" 2>/dev/null)" = "SENTINEL" ] || why="output path written"
finish E21 "-cloak with an empty value" 1 $rc "$why"

start H17
( cd "$run" && ulimit -v "$VMEM" && timeout 20 "$BIN" -cloak "$SMALL" -output - ) > "$out" 2> "$log"
rc=$?
why=""; cmp -s "$out" "$SMALL_EXP" || why="standard output differs from hand/expected/H03.fa"
finish H17 "-output - writes the output to standard output" 0 $rc "$why"

rmdir "$OUT/run" 2>/dev/null
total=$((pass+fail))
[ "$total" -eq 5 ] || { echo "ran $total cases, expected 5" >&2; exit 2; }
echo "cases $total: passed $pass, failed $fail" | tee "$OUT/SUMMARY.txt"
[ "$fail" -eq 0 ]
