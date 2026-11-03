#!/bin/bash
# Run time and peak memory against a reference binary (tests/README.md; PLAN.md, phase 1, step 5).
#
# Usage (Debian terminal):
#   bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/test_perf.sh [CAND [REF]]
#     CAND  default ~/muscle-build/dev/bin/muscle
#     REF   default ~/muscle-build/ref-ddc9df1/bin/muscle (cloak-v0.2.0, its gate A binary)
#
# Gated inputs, written once under ~/muscle-runs/test_perf/inputs/ with fixed seeds:
#   scale  tools/gen_scale.py 1 1000 1500 16 0.3 0.3 10: 1,000 rows, 1,500 columns, 16 MSAs,
#          30% of rows re-placed within 10 columns in each MSA after the first
#   pf92   tests/data/PF00188_n92.efa
#   pairs  tools/gen_large.py pairs 1000 30: every MSA-0 column splits into 500 pairs
#   k30    tools/gen_scale.py 7 500 500 30 0.3 0.3 10: 30 MSAs
#   list   scale as a path list of 16 FASTA files
# Each gated input runs 3 times per binary, alternating, through tools/runmeasure.py (60 s, 3 GB of
# address space), and the medians are compared. Fail if, for any gated input, CAND's median time
# exceeds 1.3 x REF + 0.2 s or its median peak RSS 1.2 x REF + 20 MB, or if any run of either
# binary does not exit 0. These are engineering thresholds (private/PLAN_v030.md).
# Counted only: an n-doubling series (tools/gen_scale.py 2 N 500 16 0.3 0.3 10, N = 500 to 4,000)
# on CAND, and tools/fit_mem.py's fit over every CAND median. Outputs under
# ~/muscle-runs/test_perf/<time>/: perf.tsv (every run), fit.tsv. Exit 0 only if every gated input
# passes.
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
TOOLS="$HERE/../tools"; DATA="$HERE/../tests/data"
PY=${CLOAK_PY:-"$HOME/anaconda3/envs/muscle-cloak/bin/python"}
CAND=${1:-$HOME/muscle-build/dev/bin/muscle}
REF=${2:-$HOME/muscle-build/ref-ddc9df1/bin/muscle}
[ -x "$CAND" ] && [ -x "$REF" ] || { echo "binary not found: $CAND or $REF" >&2; exit 2; }
IN="$HOME/muscle-runs/test_perf/inputs"
OUT="$HOME/muscle-runs/test_perf/$(date +%Y%m%d-%H%M%S)"
mkdir -p "$IN" "$OUT/fa"

gen_scale() { [ -s "$IN/$1.efa" ] || "$PY" "$TOOLS/gen_scale.py" "${@:2}" "$IN/$1.efa" > /dev/null; }
gen_scale scale 1 1000 1500 16 0.3 0.3 10
gen_scale k30 7 500 500 30 0.3 0.3 10
[ -s "$IN/pairs.efa" ] || "$PY" "$TOOLS/gen_large.py" pairs "$IN/pairs.efa" 1000 30 > /dev/null
cp "$DATA/PF00188_n92.efa" "$IN/pf92.efa"
if [ ! -s "$IN/list.txt" ]; then
  mkdir -p "$IN/list"
  awk -v d="$IN/list" -v l="$IN/list.txt" '/^</ {f = sprintf("%s/m%02d.fa", d, k++); print f > l; next}
       {print > f}' "$IN/scale.efa"
fi
for s in 500 1000 2000 4000; do gen_scale "series$s" 2 "$s" 500 16 0.3 0.3 10; done

cells() {  # EFA -> sum over MSAs of rows x columns
  awk '{sub(/\r$/, "")} /^</ {c += n * w; n = 0; w = 0; next} /^>/ {n++; next} n == 1 {w += length($0)}
       END {print c + n * w}' "$1"
}
outcols() { awk '/^>/ {if (s != "") exit; next} {s = s $0} END {print length(s)}' "$1"; }
measure() {  # TAG INPUT BINARY NAME -> appends one line to perf.tsv, prints "status wall rss"
  local fa="$OUT/fa/$1.$4.fa" m
  rm -f "$fa"
  m=$("$PY" "$TOOLS/runmeasure.py" 60 3 -- "$3" -cloak "$2" -output "$fa" 2> "$OUT/fa/$1.$4.err")
  printf '%s\t%s\t%s\t%s\n' "$1" "$4" "$m" "$( [ -s "$fa" ] && outcols "$fa" || echo 0)" >> "$OUT/perf.tsv"
  echo "$m"
}
median() { sort -g | sed -n 2p; }

printf 'input\tbinary\tstatus\twall_s\tpeak_rss_mb\tout_cols\n' > "$OUT/perf.tsv"
printf 'n\tcells\tout_cols\tstatus\twall_s\tpeak_rss_mb\n' > "$OUT/fit.tsv"
echo "candidate: $CAND ($(sha256sum "$CAND" | cut -c1-12))"
echo "reference: $REF ($(sha256sum "$REF" | cut -c1-12))"
echo "outputs:   $OUT"
printf "%-6s %9s %9s %9s %9s  %s\n" input cand_s ref_s cand_MB ref_MB result
bad=0
for tag in scale pf92 pairs k30 list; do
  inp="$IN/$tag.efa"; [ "$tag" = list ] && inp="$IN/list.txt"
  for rep in 1 2 3; do measure "$tag" "$inp" "$CAND" cand > /dev/null; measure "$tag" "$inp" "$REF" ref > /dev/null; done
  ct=$(awk -v t="$tag" '$1 == t && $2 == "cand" {print $4}' "$OUT/perf.tsv" | median)
  rt=$(awk -v t="$tag" '$1 == t && $2 == "ref" {print $4}' "$OUT/perf.tsv" | median)
  cm=$(awk -v t="$tag" '$1 == t && $2 == "cand" {print $5}' "$OUT/perf.tsv" | median)
  rm_=$(awk -v t="$tag" '$1 == t && $2 == "ref" {print $5}' "$OUT/perf.tsv" | median)
  notok=$(awk -v t="$tag" '$1 == t && $3 != "OK"' "$OUT/perf.tsv" | wc -l)
  res=$(awk -v ct="$ct" -v rt="$rt" -v cm="$cm" -v rm="$rm_" -v no="$notok" 'BEGIN {
    r = ""; if (no > 0) r = r no " runs not OK; "
    if (ct > 1.3 * rt + 0.2) r = r "time over 1.3 x ref + 0.2 s; "
    if (cm > 1.2 * rm + 20) r = r "RSS over 1.2 x ref + 20 MB; "
    print (r == "" ? "pass" : "FAIL: " r) }')
  [ "$res" = pass ] || bad=$((bad+1))
  printf "%-6s %9s %9s %9s %9s  %s\n" "$tag" "$ct" "$rt" "$cm" "$rm_" "$res"
  [ "$tag" = list ] && continue
  oc=$(awk -v t="$tag" '$1 == t && $2 == "cand" {print $6; exit}' "$OUT/perf.tsv")
  st=$([ "$notok" -eq 0 ] && echo OK || echo FAIL)
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$(grep -c '^>' "$inp" | awk -v k="$(grep -c '^<' "$inp")" '{print $1 / k}')" \
    "$(cells "$inp")" "$oc" "$st" "$ct" "$cm" >> "$OUT/fit.tsv"
done
echo "series (candidate only, counted):"
for s in 500 1000 2000 4000; do
  inp="$IN/series$s.efa"
  read -r st wall rss <<< "$(measure "series$s" "$inp" "$CAND" cand)"
  oc=$(awk -v t="series$s" '$1 == t {print $6; exit}' "$OUT/perf.tsv")
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$s" "$(cells "$inp")" "$oc" "$st" "$wall" "$rss" >> "$OUT/fit.tsv"
  printf "  n %5s  %s s  %s MB  %s output columns  %s\n" "$s" "$wall" "$rss" "$oc" "$st"
done
"$PY" "$TOOLS/fit_mem.py" "$OUT/fit.tsv" | head -1 | sed 's/^/memory fit (counted): /'
echo "gated inputs failing: $bad" | tee "$OUT/SUMMARY.txt"
[ "$bad" -eq 0 ]
