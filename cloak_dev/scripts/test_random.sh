#!/bin/bash
# Random synthetic ensembles (cloak_dev/tools/gen_synth.py, schedules mixed and mild) at -mincol 2
# and 3, through every oracle of tools/run_suite.py: the definition checker, the byte oracle,
# decision 003's relation with the frozen baseline, and cloak.py's functions if they are present.
#
# Usage (Debian terminal, conda env muscle-cloak created from cloak_dev/environment.yml):
#   bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/test_random.sh [COUNT] [FIRST_SEED] [BINARY]
#     COUNT       ensembles per schedule (mixed and mild), default 200, so 4*COUNT comparisons
#     FIRST_SEED  default 1: seeds FIRST.. (mixed) and FIRST+1000.. (mild)
#     BINARY      default ~/muscle-build/dev/bin/muscle
#
# The baseline runs only on ensembles predict_recursion.py rates cheap, capped at 2 GB and 60 s; a
# case where it still exceeds a cap is a named skip (skips.txt), and any skip fails the run.
# Exit status 0 only if no comparison fails and nothing is skipped.
set -uo pipefail
COUNT=${1:-200}; FIRST=${2:-1}; BIN=${3:-$HOME/muscle-build/dev/bin/muscle}
HERE=$(cd "$(dirname "$0")" && pwd)
PY=${CLOAK_PY:-"$HOME/anaconda3/envs/muscle-cloak/bin/python"}
OUT="$HOME/muscle-runs/test_random/$(date +%Y%m%d-%H%M%S)"
[ -x "$BIN" ] || { echo "binary not found: $BIN" >&2; exit 2; }
REF="$HERE/../private/reference/cloak.py"
extra=(); [ -f "$REF" ] && extra=(--cloakpy "$REF")
mkdir -p "$(dirname "$OUT")"
cd "$HERE/../tools" && "$PY" run_suite.py "$BIN" "$OUT" --count 0 --synth "$COUNT" --first "$FIRST" \
  --baseline "$HOME/muscle-baseline/muscle-f3559af" "${extra[@]}"
