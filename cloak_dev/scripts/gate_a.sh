#!/bin/bash
# Gate A (PLAN.md, phase 1, step 6): the release tier on a binary built in full from a clean commit.
#
# Usage (Debian terminal, conda off; about 3 h, keep the laptop awake):
#   bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/gate_a.sh BINARY
#     BINARY  the release binary, ~/muscle-build/portable/bin/muscle (build.sh --full from a clean
#             commit, then build_variant.sh portable), or ~/muscle-build/dev/bin/muscle
#
# Refuses unless ~/muscle-build/dev/BUILD_INFO.txt records a full build of a clean commit that is
# the build clone's HEAD, and BINARY's own BUILD_INFO.txt (beside its bin/) has that commit, a clean
# source and BINARY's SHA-256. Then builds the sanitizer variant of the same source
# (build_variant.sh asan), runs test.sh release on BINARY, and writes
# ~/muscle-runs/gate_a/<time>/GATE_A.txt: the commit, BINARY's SHA-256 and BUILD_INFO, and the
# release tier's part lines and output directory. Exit status: the release tier's.
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
DEV="$HOME/muscle-build/dev"
BIN=${1:-}
[ -n "$BIN" ] && [ -x "$BIN" ] || { sed -n '2,15p' "$0"; exit 2; }
BIN=$(readlink -f "$BIN")
INFO="$(dirname "$(dirname "$BIN")")/BUILD_INFO.txt"
refuse() { echo "gate A refused: $*" >&2; exit 2; }
field() { sed 's/^base   //' "$1" | awk -v k="$2:" '$1 == k || ($1 " " $2) == k {sub(/^[^:]*: */, ""); print; exit}'; }
[ -f "$INFO" ] || refuse "no BUILD_INFO.txt beside $(dirname "$BIN")"
[ "$(field "$DEV/BUILD_INFO.txt" mode)" = full ] || refuse "the dev build is not a full build (build.sh --full)"
[ "$(field "$DEV/BUILD_INFO.txt" 'src state')" = "clean (equals the commit)" ] || refuse "the dev build has uncommitted source"
COMMIT=$(field "$DEV/BUILD_INFO.txt" commit)
[ "$COMMIT" = "$(git -C "$DEV" rev-parse HEAD)" ] || refuse "the dev build clone is not at the built commit"
[ "$(field "$INFO" commit)" = "$COMMIT" ] || refuse "$BIN is not built from $COMMIT"
[ "$(field "$INFO" 'src state')" = "clean (equals the commit)" ] || refuse "$BIN has uncommitted source"
SHA=$(sha256sum "$BIN" | cut -d' ' -f1)
[ "$(grep '^sha256:' "$INFO" | awk '{print $2}')" = "$SHA" ] || refuse "$BIN does not match the SHA-256 in $INFO"

OUT="$HOME/muscle-runs/gate_a/$(date +%Y%m%d-%H%M%S)"
mkdir -p "$OUT"
echo "gate A: $BIN ($SHA), commit $COMMIT; outputs $OUT"
bash "$HERE/build_variant.sh" asan > "$OUT/build_asan.log" 2>&1 || { tail -5 "$OUT/build_asan.log"; exit 1; }
bash "$HERE/test.sh" release "$BIN" 2>&1 | tee "$OUT/release.txt"
rc=${PIPESTATUS[0]}
{ echo "gate A, $(date -Is)"
  echo "commit:  $COMMIT"
  echo "binary:  $BIN"
  echo "sha256:  $SHA"
  echo "result:  $([ "$rc" -eq 0 ] && echo PASS || echo "FAIL (release tier exit $rc)")"
  echo; echo "BUILD_INFO:"; sed 's/^/  /' "$INFO"
  echo; echo "release tier:"; sed 's/^/  /' "$OUT/release.txt"
} > "$OUT/GATE_A.txt"
echo "gate A record: $OUT/GATE_A.txt"
exit "$rc"
