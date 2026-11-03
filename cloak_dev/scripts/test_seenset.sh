#!/bin/bash
# Oracle O5: the seen-set variant of f3559af (tools/make_seenset.py) on the fixtures that the frozen
# baseline cannot finish (tests/README.md).
#
# Usage (Debian terminal, conda off):
#   bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/test_seenset.sh [--control] [BINARY]
#     BINARY     default ~/muscle-build/dev/bin/muscle
#     --control  check the variant instead: on the MANIFEST cases whose baseline column is fast, its
#                output must equal the expected output byte for byte (those cases keep every row)
#
# The variant is f3559af's cloak.cpp, read from the build clone ~/muscle-build/dev, with the
# recursion memoised and the final de-duplication hashed, which keeps its output and makes the large
# fixtures finish (2026-10-30 audit); build_variant.sh seenset builds it. On every MANIFEST case
# whose baseline column is never, BINARY's output with its rows of at most one residue removed
# (tools/cloakio.py strip) must equal the variant's output byte for byte: decision 003's relation
# (the fixtures have no '.'). Every run goes through tools/runmeasure.py with 600 s and 3 GB.
# Outputs under ~/muscle-runs/test_seenset/<time>/. Exit 0 only if every case passes.
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
TOOLS="$HERE/../tools"; DATA="$HERE/../tests/data"
PY=${CLOAK_PY:-"$HOME/anaconda3/envs/muscle-cloak/bin/python"}
BIN="$HOME/muscle-build/dev/bin/muscle"; CONTROL=0
for a in "$@"; do case "$a" in --control) CONTROL=1 ;; -*) echo "unknown option $a" >&2; exit 2 ;; *) BIN=$a ;; esac; done
[ -x "$BIN" ] || { echo "binary not found: $BIN" >&2; exit 2; }
bash "$HERE/check_fixtures.sh" || exit 2
OUT="$HOME/muscle-runs/test_seenset/$(date +%Y%m%d-%H%M%S)"
mkdir -p "$OUT/fa"
git -C "$HOME/muscle-build/dev" show f3559af:src/cloak.cpp > "$OUT/cloak_f3559af.cpp" || exit 2
"$PY" "$TOOLS/make_seenset.py" "$OUT/cloak_f3559af.cpp" "$OUT/cloak_seenset.cpp" > /dev/null || exit 2
bash "$HERE/build_variant.sh" seenset "$OUT/cloak_seenset.cpp" > "$OUT/build.log" 2>&1 ||
  { tail -5 "$OUT/build.log"; exit 2; }
V="$HOME/muscle-build/seenset/bin/muscle"
want=never; [ "$CONTROL" = 1 ] && want=fast
echo "variant: $V ($(sha256sum "$V" | cut -c1-12)); mode: $([ "$CONTROL" = 1 ] && echo control || echo "S9 against $BIN")"
echo "outputs: $OUT"
pass=0; fail=0
while IFS=$'\t' read -r efa mc exp src base || [ -n "$efa" ]; do
  base=${base%$'\r'}
  [ "$efa" = ensemble ] && continue
  [ "$base" = "$want" ] || continue
  name="$(basename "$efa" .efa).m$mc"
  m=$("$PY" "$TOOLS/runmeasure.py" 600 3 -- "$V" -cloak "$DATA/$efa" -mincol "$mc" \
        -output "$OUT/fa/$name.seenset.fa" 2> "$OUT/fa/$name.seenset.err")
  why=""
  if [ "$(cut -f1 <<< "$m")" != OK ]; then why="variant: $m"
  elif [ "$CONTROL" = 1 ]; then
    cmp -s "$OUT/fa/$name.seenset.fa" "$DATA/$exp" || why="variant output differs from $exp"
  else
    c=$("$PY" "$TOOLS/runmeasure.py" 600 3 -- "$BIN" -cloak "$DATA/$efa" -mincol "$mc" \
          -output "$OUT/fa/$name.cand.fa" 2> "$OUT/fa/$name.cand.err")
    if [ "$(cut -f1 <<< "$c")" != OK ]; then why="binary: $c"
    else
      "$PY" "$TOOLS/cloakio.py" strip "$OUT/fa/$name.cand.fa" "$OUT/fa/$name.cand.strip.fa"
      cmp -s "$OUT/fa/$name.cand.strip.fa" "$OUT/fa/$name.seenset.fa" ||
        why="binary output without rows of at most one residue differs from the variant's"
    fi
  fi
  if [ -z "$why" ]; then pass=$((pass+1)); printf "pass  %-28s m%s  variant %s\n" "$efa" "$mc" "$m"
  else fail=$((fail+1)); printf "fail  %-28s m%s  (%s)\n" "$efa" "$mc" "$why"; fi
done < "$DATA/MANIFEST.tsv"
echo "cases $((pass+fail)): passed $pass, failed $fail" | tee "$OUT/SUMMARY.txt"
[ "$fail" -eq 0 ] && [ $((pass+fail)) -gt 0 ]
