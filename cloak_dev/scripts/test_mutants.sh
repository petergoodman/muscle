#!/bin/bash
# Mutation testing: every mutant of src/cloak.cpp that changes the output must be killed by the
# tests (tests/README.md).
#
# Usage (Debian terminal, conda off):
#   bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/test_mutants.sh [ONLY]
#     ONLY   run only the mutants whose names contain this text
#
# The mutants (tools/make_mutants.py) are made from ~/muscle-build/dev/src/cloak.cpp, the source of
# the current build, and each is built with build_variant.sh mutant (one compile of cloak.cpp).
# Each runs test.sh edit, then test.sh commit if edit passes, then, if both pass, test_sanitize.sh
# on its sanitizer build (build_variant.sh mutant-asan; build_variant.sh asan must have been run on
# the same source). A mutant is killed by the first part that fails; a mutant that does not compile
# is an error, not a kill. Pass: the control and the equivalent mutants survive, and every killable
# mutant is killed. results.tsv: name, kind, expected killer, result, killed by, failing parts.
# Outputs under ~/muscle-runs/test_mutants/<time>/. About an hour: keep the laptop awake.
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
TOOLS="$HERE/../tools"
PY=${CLOAK_PY:-"$HOME/anaconda3/envs/muscle-cloak/bin/python"}
MB="$HOME/muscle-build"
ONLY=${1:-}
src_of() { sed 's/^base   //' "$1/BUILD_INFO.txt" | grep -E '^(commit|src state):'; }
[ -x "$MB/asan/bin/muscle" ] && [ "$(src_of "$MB/asan")" = "$(src_of "$MB/dev")" ] ||
  { echo "the asan build is missing or not of dev's source: run build_variant.sh asan" >&2; exit 2; }
OUT="$HOME/muscle-runs/test_mutants/$(date +%Y%m%d-%H%M%S)"
mkdir -p "$OUT/log"
"$PY" "$TOOLS/make_mutants.py" "$MB/dev/src/cloak.cpp" "$OUT/src" || exit 2
echo "source:  $MB/dev ($(src_of "$MB/dev" | tr '\n' ' '))"
echo "outputs: $OUT"
printf 'name\tkind\texpected_killer\tresult\tkilled_by\tfailing_parts\n' > "$OUT/results.tsv"
bad=0; n=0
while IFS=$'\t' read -r name kind killer; do
  [ "$name" = name ] && continue
  [ -n "$ONLY" ] && [[ "$name" != *"$ONLY"* ]] && continue
  n=$((n+1)); cpp="$OUT/src/$name/cloak.cpp"; log="$OUT/log/$name"; mkdir -p "$log"
  by=""; parts=""
  if ! bash "$HERE/build_variant.sh" mutant "$name" "$cpp" > "$log/build.log" 2>&1; then by=BUILD_ERROR
  else
    bin="$MB/mutants/$name/muscle"
    for tier in edit commit; do
      if ! bash "$HERE/test.sh" "$tier" "$bin" > "$log/$tier.log" 2>&1; then
        by=$tier; parts=$(awk '$1 == "FAIL" {printf "%s%s", s, $2; s=","}' "$log/$tier.log"); break
      fi
    done
    if [ -z "$by" ]; then
      if ! bash "$HERE/build_variant.sh" mutant-asan "$name" "$cpp" > "$log/build_asan.log" 2>&1; then
        by=BUILD_ERROR
      elif ! bash "$HERE/test_sanitize.sh" "$MB/mutants/$name-asan/muscle" > "$log/sanitize.log" 2>&1; then
        by=sanitizers; parts=$(awk '/^[a-z]+ +cases/ && !/failed 0/ {printf "%s%s", s, $1; s=","}' "$log/sanitize.log")
      fi
    fi
  fi
  if [ "$by" = BUILD_ERROR ]; then res=ERROR
  elif [ "$kind" = killable ]; then [ -n "$by" ] && res=killed || res=SURVIVED
  else [ -z "$by" ] && res=survived || res="KILLED (expected to survive)"
  fi
  case "$res" in killed|survived) ;; *) bad=$((bad+1)) ;; esac
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$kind" "$killer" "$res" "${by:--}" "${parts:--}" >> "$OUT/results.tsv"
  printf "%-32s %-10s %-28s %s\n" "$name" "$kind" "$res${by:+ by $by}" "${parts}"
done < "$OUT/src/MUTANTS.tsv"
echo "mutants $n: not as required $bad" | tee "$OUT/SUMMARY.txt"
[ "$n" -gt 0 ] && [ "$bad" -eq 0 ]
