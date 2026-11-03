#!/bin/bash
# One entry point for the test tiers (cloak_dev/tests/README.md).
#
# Usage (Debian terminal, conda env muscle-cloak):
#   bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/test.sh {edit|commit|release} [BINARY]
#     edit     after every build, about 1 min: fixture check, checker unit tests, hand cases,
#              test_errors.sh, and 30 tiny random cases through every oracle and relation
#     commit   before every commit, about 5 min: edit, plus test_oom.sh (not on the byte oracle),
#              test_quick.sh --check-baseline and 300 gen_cases + 100 gen_synth cases with
#              coverage asserted
#     release  gate A (scripts/gate_a.sh), about 3 h: commit, plus --check-baseline-slow, 2,000 +
#              800 more random cases, sanitizers (needs build_variant.sh asan of the same source),
#              warnings, mutants, performance against ddc9df1, digit-plane tagging and the
#              seen-set oracle; the first four are n/a on the byte oracle
#   BINARY     default ~/muscle-build/dev/bin/muscle; scripts/oracle_bin.sh tests the harness itself
# Prints one line per part; logs and outputs under ~/muscle-runs/test/<tier>-<time>/.
# Exit status 0 only if every part passes.
set -uo pipefail
TIER=${1:-}; BIN=${2:-$HOME/muscle-build/dev/bin/muscle}
HERE=$(cd "$(dirname "$0")" && pwd)
TOOLS="$HERE/../tools"
PY=${CLOAK_PY:-"$HOME/anaconda3/envs/muscle-cloak/bin/python"}
BASE="$HOME/muscle-baseline/muscle-f3559af"
REF="$HERE/../private/reference/cloak.py"
case "$TIER" in edit|commit|release) ;; *) echo "usage: test.sh {edit|commit|release} [BINARY]" >&2; exit 2 ;; esac
[ -e "$BIN" ] || { echo "binary not found: $BIN" >&2; exit 2; }
OUT="$HOME/muscle-runs/test/$TIER-$(date +%Y%m%d-%H%M%S)"
mkdir -p "$OUT"
ref=()
if [ -f "$REF" ]; then ref=(--cloakpy "$REF"); else echo "note: cloak.py not found, oracle O3 not run"; fi
nfail=0
part() {  # NAME COMMAND...: runs the command, logs it, prints pass or FAIL with the log's last line
  local name=$1 t0 rc res
  shift
  t0=$(date +%s)
  "$@" > "$OUT/$name.log" 2>&1
  rc=$?
  res=pass
  [ $rc -eq 0 ] || { res=FAIL; nfail=$((nfail+1)); }
  printf "%-5s %-12s %5ss  %s\n" "$res" "$name" "$(( $(date +%s) - t0 ))" \
    "$(grep -v '^ERROR conda' "$OUT/$name.log" | tail -1)"
}
ORACLE=0; [ "$(readlink -f "$BIN")" = "$(readlink -f "$HERE/oracle_bin.sh")" ] && ORACLE=1
na() { printf "n/a   %-12s %s\n" "$1" "$2"; }

echo "tier $TIER; binary $BIN ($(sha256sum "$BIN" | cut -c1-12)); outputs $OUT"
part fixtures bash "$HERE/check_fixtures.sh"
part checker "$PY" "$TOOLS/test_spec_check.py"
part hand bash "$HERE/test_hand.sh" "$BIN"
part errors bash "$HERE/test_errors.sh" "$BIN"
part tiny "$PY" "$TOOLS/run_suite.py" "$BIN" "$OUT/tiny" --schedules tiny --count 30 --relations 1 "${ref[@]}"
if [ "$TIER" != edit ]; then
  if [ "$ORACLE" = 1 ]; then na oom "the byte oracle has no memory behaviour of its own (test_oom.sh)"
  else part oom bash "$HERE/test_oom.sh" "$BIN"; fi
  part quick bash "$HERE/test_quick.sh" --check-baseline "$BIN"
  part random "$PY" "$TOOLS/run_suite.py" "$BIN" "$OUT/random" --count 60 --synth 50 --relations 3 \
    --baseline "$BASE" "${ref[@]}" --assert-coverage
fi
if [ "$TIER" = release ]; then
  part quick-slow bash "$HERE/test_quick.sh" --check-baseline-slow "$BIN"
  part random-big "$PY" "$TOOLS/run_suite.py" "$BIN" "$OUT/random-big" --count 400 --first 1001 \
    --synth 400 --relations 5 --baseline "$BASE" "${ref[@]}" --assert-coverage
  if [ "$ORACLE" = 1 ]; then
    for p in sanitizers warnings mutants performance; do na "$p" "tests the C++ source or binary, not the byte oracle"; done
  else
    part sanitizers bash "$HERE/test_sanitize.sh" "$HOME/muscle-build/asan/bin/muscle" "$BIN"
    part warnings bash "$HERE/warn_diff.sh"
    part mutants bash "$HERE/test_mutants.sh"
    part performance bash "$HERE/test_perf.sh" "$BIN"
  fi
  part tagging "$PY" "$TOOLS/tag_runs.py" "$BIN" "$OUT/tagging"
  part seenset bash "$HERE/test_seenset.sh" "$BIN"
fi
echo "tier $TIER: $nfail part(s) failed" | tee "$OUT/SUMMARY.txt"
[ $nfail -eq 0 ]
