#!/bin/bash
# Gate B: a CLOAK binary on the users' own ensembles, run locally (PLAN.md, phase 1; where the data
# are and the pass lines in private/PROJECT_NOTES.md). The script names no data.
#
# Usage (Debian terminal): bash .../cloak_dev/scripts/gate_b.sh DATA LIST BINARY MODE
#   DATA    directory with ens/<name>.efa (an .efa with an .efa.tmp beside it is incomplete and is
#           left out) and cloak/<name>.fa, the reference outputs written by cloak-v0.1.0
#   LIST    file naming the ensembles that cloak-v0.1.0 could not finish, one per line
#   BINARY  the binary under test
#   MODE    bytes  every reference reproduced byte for byte (cloak-v0.2.0, decision 003)
#           strip  every reference equals the output without rows of at most one residue, and
#                  every output passes tools/spec_check.py --order with every sequence (0.3.0)
# Each run goes through tools/runmeasure.py with 3600 s and 5 GB of address space. results.tsv: name,
# kind (reference, failed, other), n, K, MSA-0 columns, output columns, status, wall s, peak RSS MB,
# result. In bytes mode scan.tsv also gives spec_check.py stats of each reference (rows in, rows
# out, rows missing, columns, columns below 2, all-gap columns). Outputs under
# ~/muscle-runs/gate_b/<time>/. Exit 0 only if every reference and every listed ensemble passes.
set -uo pipefail
[ $# -eq 4 ] || { sed -n '2,19p' "$0"; exit 2; }
DATA=$1; LIST=$2; BIN=$3; MODE=$4
HERE=$(cd "$(dirname "$0")" && pwd)
TOOLS="$HERE/../tools"
PY=${CLOAK_PY:-"$HOME/anaconda3/envs/muscle-cloak/bin/python"}
case "$MODE" in bytes|strip) ;; *) echo "MODE must be bytes or strip" >&2; exit 2 ;; esac
[ -x "$BIN" ] && [ -d "$DATA/ens" ] && [ -f "$LIST" ] || { echo "missing BINARY, DATA/ens or LIST" >&2; exit 2; }
OUT="$HOME/muscle-runs/gate_b/$(date +%Y%m%d-%H%M%S)-$MODE"
mkdir -p "$OUT/fa" "$OUT/log"
failed=$(tr -d '\r' < "$LIST")
{ echo "binary: $BIN"; echo "sha256: $(sha256sum "$BIN" | cut -d' ' -f1)"
  echo "version: $("$BIN" -version 2>&1 | head -1)"; echo "mode: $MODE"; echo "data: $DATA"; } > "$OUT/BINARY_INFO.txt"
cat "$OUT/BINARY_INFO.txt"
printf 'name\tkind\tn\tK\tcols0\tout_cols\tstatus\twall_s\trss_mb\tresult\n' > "$OUT/results.tsv"
[ "$MODE" = bytes ] && printf 'name\trows_in\trows_out\trows_missing\tcols\tcols_below_2\tcols_allgap\n' > "$OUT/scan.tsv"

dims() {  # EFA -> "n K cols0"
  awk '/^</ {k++; next} k == 1 && /^>/ {n++; next} k == 1 && n == 1 {w += length($0)}
       END {print n + 0, k + 0, w + 0}' "$1"
}
bad=0; incomplete=0
for efa in "$DATA"/ens/*.efa; do
  name=$(basename "$efa" .efa)
  if [ -e "$efa.tmp" ]; then incomplete=$((incomplete+1)); continue; fi
  ref="$DATA/cloak/$name.fa"
  if [ -f "$ref" ]; then kind=reference
  elif printf '%s\n' "$failed" | grep -qx "$name"; then kind=failed
  else kind=other; fi
  out="$OUT/fa/$name.fa"
  read -r status wall rss < <("$PY" "$TOOLS/runmeasure.py" 3600 5 -- "$BIN" -cloak "$efa" -output "$out" \
                               2> "$OUT/log/$name.err")
  outcols=0; [ -s "$out" ] && outcols=$(awk 'NR > 1 && /^>/ {exit} NR > 1 {w += length($0)} END {print w + 0}' "$out")
  if [ "$status" != OK ]; then result="FAIL:$status"
  elif [ ! -s "$out" ]; then result=FAIL:no-output
  elif [ "$MODE" = bytes ]; then
    if [ "$kind" = reference ]; then cmp -s "$out" "$ref" && result=IDENTICAL || result=FAIL:DIFFERENT
    else result=FINISHED; fi
  else
    result=CHECKED
    if ! "$PY" "$TOOLS/spec_check.py" check "$efa" "$out" --order > /dev/null 2>> "$OUT/log/$name.err"; then
      result=FAIL:spec_check
    elif [ "$kind" = reference ]; then
      "$PY" "$TOOLS/cloakio.py" strip "$out" "$OUT/fa/$name.strip"
      cmp -s "$OUT/fa/$name.strip" "$ref" && result=RELATION_OK || result=FAIL:relation
      rm -f "$OUT/fa/$name.strip"
    fi
  fi
  [ "$kind" = other ] || case "$result" in FAIL*) bad=$((bad+1)) ;; esac
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$kind" "$(dims "$efa" | tr ' ' '\t')" \
    "$outcols" "$status" "$wall" "$rss" "$result" | tee -a "$OUT/results.tsv"
  if [ "$MODE" = bytes ] && [ "$kind" = reference ]; then
    printf '%s\t%s\n' "$name" "$("$PY" "$TOOLS/spec_check.py" stats "$efa" "$ref")" >> "$OUT/scan.tsv"
  fi
done
summary=$(awk -F'\t' 'NR > 1 {c[$2 " " $10]++} END {for (k in c) printf "%s %d; ", k, c[k]}' "$OUT/results.tsv")
echo "gate B ($MODE): ${summary}incomplete ensembles left out $incomplete; failures $bad" | tee "$OUT/SUMMARY.txt"
[ $bad -eq 0 ]
