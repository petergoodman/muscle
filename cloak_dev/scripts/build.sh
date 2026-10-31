#!/bin/bash
# Build the Windows working copy's muscle in Debian, on the Linux filesystem.
#
# Usage (Debian terminal, conda off):  bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/build.sh [--full]
#
# ~/muscle-build/dev is a Linux clone of the working copy, checked out at the working copy's HEAD.
# The tracked files under src/ are mirrored into it with CRLF removed, rewriting only files whose
# content differs so that make stays incremental. Uncommitted edits are therefore built too, and
# BUILD_INFO.txt says so. --full (and the first run) does a clean build through vcxproj_make.py,
# exactly as build_linux.bash does, which also stamps the commit into `muscle -version`; an
# incremental build leaves that stamp at the last full build. A gate build is --full from a clean
# commit. Nothing is written outside ~/muscle-build/dev, and the working copy is only read.
set -euo pipefail
WIN="${MUSCLE_WIN:-/mnt/c/Users/peter/BioInformatics/muscle}"   # overridable for testing
B="${MUSCLE_BUILD_DIR:-$HOME/muscle-build/dev}"
FULL=0
[ "${1:-}" = "--full" ] && FULL=1

[ -d "$WIN/.git" ] || { echo "working copy not found: $WIN" >&2; exit 1; }
if [ ! -d "$B/.git" ]; then
  echo "first run: cloning $WIN into $B"
  git clone --quiet --no-checkout "$WIN" "$B"
  FULL=1
fi

HEAD=$(git -C "$WIN" rev-parse HEAD)            # read-only on the working copy's repository
git -C "$B" fetch --quiet origin
git -C "$B" checkout --quiet --force --detach "$HEAD"

# Mirror tracked src/ files (all are text; verified at f3559af) plus any new .cpp/.h files.
changed=0
cd "$B"
{ git ls-files -- src; (cd "$WIN" && find src -maxdepth 1 -type f \( -name '*.cpp' -o -name '*.h' \)); } \
  | sort -u | while IFS= read -r f; do
    [ -f "$WIN/$f" ] || continue
    tmp=$(mktemp "$B/.mirror.XXXXXX")
    tr -d '\r' < "$WIN/$f" > "$tmp"
    if cmp -s "$tmp" "$B/$f"; then rm -f "$tmp"; else mv "$tmp" "$B/$f"; echo "$f" >> "$B/.mirror.changed"; fi
  done
[ -f "$B/.mirror.changed" ] && changed=$(wc -l < "$B/.mirror.changed")
rm -f "$B/.mirror.changed"
DIRTY=$(git -C "$B" status --porcelain --untracked-files=no -- src)

cd "$B/src"
if [ "$FULL" -eq 1 ] || [ ! -f Makefile ]; then
  [ -f vcxproj_make.py ] || curl -fsSL https://raw.githubusercontent.com/rcedgar/vcxproj_make/806d016/vcxproj_make.py -o vcxproj_make.py
  mode=full
  /usr/bin/python3 vcxproj_make.py --openmp > build.log 2>&1 || { tail -20 build.log; exit 1; }
else
  mode=incremental
  make > build.log 2>&1 || { tail -20 build.log; exit 1; }
fi

BIN="$B/bin/muscle"
{
  echo "built:     $(date -Is)"
  echo "mode:      $mode"
  echo "commit:    $HEAD"
  if [ -n "$DIRTY" ]; then echo "src state: uncommitted changes"; echo "$DIRTY" | sed 's/^/           /'
  else echo "src state: clean (equals the commit)"; fi
  echo "sha256:    $(sha256sum "$BIN" | cut -d' ' -f1)"
  echo "version:   $("$BIN" -version 2>&1 | head -1)"
} > "$B/BUILD_INFO.txt"
cat "$B/BUILD_INFO.txt"
echo "binary:    $BIN   (files rewritten this run: $changed)"
