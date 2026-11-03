#!/bin/bash
# Compiler warnings in cloak.cpp that a reference version does not have (tests/README.md).
#
# Usage (Debian terminal):
#   bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/warn_diff.sh [CAND_CPP [REF_CPP]]
#     CAND_CPP  default ~/muscle-build/dev/src/cloak.cpp (the source of the current build)
#     REF_CPP   default ddc9df1's src/cloak.cpp, read from the build clone ~/muscle-build/dev
#
# Each file is compiled (-c -o /dev/null, so that warnings that need the optimiser appear) with the
# Makefile's flags and ~/muscle-build/dev/src on the include path, by g++ and clang++, twice:
# -Wall -Wextra (gating) and with -Wconversion -Wsign-conversion added (counted only). A warning is
# compared as file name, message and the text of its source line, without line numbers, so code
# that moved is not new. Outputs go to a new directory under ~/muscle-runs/warn_diff/.
# Exit status 0 only if no gating warning is new.
set -uo pipefail
DEV="$HOME/muscle-build/dev"
CAND=${1:-$DEV/src/cloak.cpp}
OUT="$HOME/muscle-runs/warn_diff/$(date +%Y%m%d-%H%M%S)"
mkdir -p "$OUT/cand" "$OUT/ref"
tr -d '\r' < "$CAND" > "$OUT/cand/cloak.cpp" || exit 2
if [ -n "${2:-}" ]; then tr -d '\r' < "$2" > "$OUT/ref/cloak.cpp" || exit 2
else git -C "$DEV" show ddc9df1:src/cloak.cpp | tr -d '\r' > "$OUT/ref/cloak.cpp" || exit 2; fi
FLAGS="-ffast-math -march=native -O3 -DNDEBUG -fopenmp -I $DEV/src"
GATE="-Wall -Wextra"
INFO="-Wall -Wextra -Wconversion -Wsign-conversion"

keys() {  # LOG -> one key per warning: file|message|source line text
  grep -E '^[^ ]+:[0-9]+:[0-9]+: warning: ' "$1" | while IFS= read -r w; do
    f=${w%%:*}; rest=${w#*:}; line=${rest%%:*}
    src=$(sed -n "${line}p" "$f" 2>/dev/null | sed -E 's/^[[:space:]]+|[[:space:]]+$//g')
    printf '%s|%s|%s\n' "$(basename "$f")" "${w#*: warning: }" "$src"
  done | sort -u
}
echo "candidate: $CAND"; echo "outputs:   $OUT"
new_gate=0
for cxx in g++ clang++; do
  for set in gate info; do
    [ "$set" = gate ] && W=$GATE || W=$INFO
    for v in cand ref; do
      # shellcheck disable=SC2086
      $cxx $FLAGS $W -c "$OUT/$v/cloak.cpp" -o /dev/null > "$OUT/$cxx.$set.$v.log" 2>&1 ||
        { echo "$cxx failed to compile $v (see $OUT/$cxx.$set.$v.log)"; exit 2; }
      keys "$OUT/$cxx.$set.$v.log" > "$OUT/$cxx.$set.$v.keys"
    done
    comm -23 "$OUT/$cxx.$set.cand.keys" "$OUT/$cxx.$set.ref.keys" > "$OUT/$cxx.$set.new"
    n=$(wc -l < "$OUT/$cxx.$set.new")
    printf "%-8s %-4s warnings: candidate %3s, reference %3s, new %3s\n" "$cxx" "$set" \
      "$(wc -l < "$OUT/$cxx.$set.cand.keys")" "$(wc -l < "$OUT/$cxx.$set.ref.keys")" "$n"
    if [ "$set" = gate ] && [ "$n" -gt 0 ]; then new_gate=$((new_gate+n)); sed 's/^/  new: /' "$OUT/$cxx.$set.new"; fi
  done
done
echo "new gating warnings: $new_gate" | tee "$OUT/SUMMARY.txt"
[ "$new_gate" -eq 0 ]
