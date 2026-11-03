#!/bin/bash
# Run the hand cases, the error cases, the 23 fixtures and 30 tiny random cases on the sanitizer
# build (scripts/build_variant.sh asan), so that every AddressSanitizer, UBSan or
# _GLIBCXX_ASSERTIONS report fails a test through its exit status (tests/README.md).
#
# Usage (Debian terminal):
#   bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/test_sanitize.sh [ASAN_BINARY [BINARY]]
#     ASAN_BINARY  default ~/muscle-build/asan/bin/muscle
#     BINARY       if given, the sanitizer build must be of the same commit and source state as
#                  BINARY's build (the BUILD_INFO.txt beside each bin/), as gate A requires
#
# A sanitizer report exits 86 (an assertion aborts, 134), which no test accepts. UBSan is
# recoverable in this build only so that scripts/ubsan.supp can silence the null memcpy argument
# in MUSCLE's fasta2.cpp (2026-10-30 audit); halt_on_error stops on every other report. The tests'
# address-space caps are lifted (CLOAK_VMEM_KB=unlimited), since ASan reserves terabytes of
# virtual memory; hard_rss_limit_mb caps real memory instead. Leak checking is off (MUSCLE frees
# little at exit), and so is alloc_dealloc_mismatch (MUSCLE's own new/free pairs).
# Exit status 0 only if every part passes.
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
TOOLS="$HERE/../tools"
PY=${CLOAK_PY:-"$HOME/anaconda3/envs/muscle-cloak/bin/python"}
ASAN=${1:-$HOME/muscle-build/asan/bin/muscle}
[ -x "$ASAN" ] || { echo "sanitizer binary not found: $ASAN (build_variant.sh asan)" >&2; exit 2; }
source_of() {  # BINARY -> its build's commit and src state lines
  sed 's/^base   //' "$(dirname "$(dirname "$1")")/BUILD_INFO.txt" | grep -E '^(commit|src state):'
}
if [ -n "${2:-}" ] && [ "$(source_of "$ASAN")" != "$(source_of "$2")" ]; then
  echo "the sanitizer build is not of the binary's source:"; source_of "$ASAN"; source_of "$2"; exit 1
fi
export ASAN_OPTIONS=exitcode=86:detect_leaks=0:alloc_dealloc_mismatch=0:hard_rss_limit_mb=3000
export UBSAN_OPTIONS=halt_on_error=1:exitcode=86:print_stacktrace=1:suppressions=$HERE/ubsan.supp
export CLOAK_VMEM_KB=unlimited
OUT="$HOME/muscle-runs/test_sanitize/$(date +%Y%m%d-%H%M%S)"
mkdir -p "$OUT"
echo "sanitizer binary: $ASAN ($(sha256sum "$ASAN" | cut -c1-12)); outputs $OUT"
nfail=0
step() {  # NAME COMMAND...
  local name=$1; shift
  "$@" > "$OUT/$name.log" 2>&1 || nfail=$((nfail+1))
  printf "%-6s %s\n" "$name" "$(tail -1 "$OUT/$name.log")"
}
step hand bash "$HERE/test_hand.sh" "$ASAN"
step errors bash "$HERE/test_errors.sh" "$ASAN"
step quick bash "$HERE/test_quick.sh" "$ASAN"
step tiny "$PY" "$TOOLS/run_suite.py" "$ASAN" "$OUT/tiny" --schedules tiny --count 30 --relations 1
echo "sanitizer parts failed: $nfail" | tee "$OUT/SUMMARY.txt"
[ "$nfail" -eq 0 ]
