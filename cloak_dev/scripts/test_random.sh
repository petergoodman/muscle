#!/bin/bash
# Byte-compare a candidate binary with the frozen baseline on random synthetic ensembles
# (cloak_dev/tools/gen_synth.py), at -mincol 2 and 3. For changes that must keep output identical
# (decision 003).
#
# Usage (Debian terminal, conda env muscle-cloak created from cloak_dev/environment.yml):
#   bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/test_random.sh [COUNT] [FIRST_SEED] [BINARY]
#     COUNT       ensembles per schedule (mixed and mild), default 200, so 4*COUNT comparisons
#     FIRST_SEED  default 1; the sandbox runs used seeds 1-400 (mixed) and 1001-1400 (mild)
#     BINARY      default ~/muscle-build/dev/bin/muscle
#
# The baseline runs capped at 2 GB and 60 s; a case where it times out is skipped and counted.
# Exit status 0 only if no comparison differs.
set -uo pipefail
COUNT=${1:-200}; FIRST=${2:-1}; BIN=${3:-$HOME/muscle-build/dev/bin/muscle}
HERE=$(cd "$(dirname "$0")" && pwd)
BASE="$HOME/muscle-baseline/muscle-f3559af"
PY=${CLOAK_PY:-"$HOME/anaconda3/bin/conda run --no-capture-output -n muscle-cloak python"}
OUT="$HOME/muscle-runs/test_random/$(date +%Y%m%d-%H%M%S)"
mkdir -p "$OUT/efa"
[ -x "$BIN" ] || { echo "binary not found: $BIN" >&2; exit 2; }

$PY "$HERE/../tools/gen_synth.py" --batch "$FIRST" "$COUNT" "$OUT/efa" mixed || exit 2
$PY "$HERE/../tools/gen_synth.py" --batch "$((FIRST + 1000))" "$COUNT" "$OUT/efa" mild || exit 2

tot=0; skip=0; bad=0; nonempty=0
for efa in "$OUT"/efa/*.efa; do
  for mc in 2 3; do
    rm -f "$OUT/b.fa" "$OUT/c.fa"
    ( ulimit -v 2000000; timeout 60 "$BASE" -cloak "$efa" -mincol $mc -output "$OUT/b.fa" > /dev/null 2>&1 ); eb=$?
    if [ $eb -eq 124 ] || [ $eb -ge 128 ]; then skip=$((skip+1)); continue; fi
    ( ulimit -v 2000000; timeout 60 "$BIN" -cloak "$efa" -mincol $mc -output "$OUT/c.fa" > /dev/null 2>&1 ); ec=$?
    tot=$((tot+1))
    if [ $eb -ne $ec ]; then bad=$((bad+1)); echo "exit status differs: $(basename "$efa") mincol $mc (baseline $eb, candidate $ec)"; continue; fi
    if [ -f "$OUT/b.fa" ]; then
      nonempty=$((nonempty+1))
      cmp -s "$OUT/b.fa" "$OUT/c.fa" || { bad=$((bad+1)); echo "output differs: $(basename "$efa") mincol $mc"; }
    elif [ -f "$OUT/c.fa" ]; then bad=$((bad+1)); echo "candidate wrote output, baseline did not: $(basename "$efa") mincol $mc"; fi
  done
done
echo "compared $tot (non-empty output $nonempty), skipped $skip (baseline over its cap), differing $bad" | tee "$OUT/SUMMARY.txt"
[ $bad -eq 0 ]
