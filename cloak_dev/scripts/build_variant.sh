#!/bin/bash
# Build a variant of the binary in ~/muscle-build/dev, each in its own directory (tests/README.md).
#
# Usage (Debian terminal, conda off):
#   bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/build_variant.sh VARIANT [ARGS]
#     asan             ~/muscle-build/asan/bin/muscle: the whole tree with -O1 -g
#                      -fsanitize=address,undefined (UBSan recoverable, so scripts/ubsan.supp can
#                      silence the one known report in MUSCLE's fasta2.cpp) -D_GLIBCXX_ASSERTIONS,
#                      not static, not stripped (the 2026-10-30 audit's recipe)
#     mutant ID CPP    ~/muscle-build/mutants/ID/muscle: dev's tree with src/cloak.cpp replaced by CPP
#     mutant-asan ID CPP  the same on the asan tree: ~/muscle-build/mutants/ID-asan/muscle
#     seenset CPP      ~/muscle-build/seenset/bin/muscle: ~/muscle-build/ref-f3559af with
#                      src/cloak.cpp replaced by CPP (tools/make_seenset.py)
#     portable         ~/muscle-build/portable/bin/muscle: dev's build without -march=native, the
#                      release binary (PLAN.md, phase 1, step 6)
# asan and portable start from a fresh copy of ~/muscle-build/dev and rebuild every object; mutant
# and seenset reuse one copy per base and recompile only cloak.cpp. ~/muscle-build/dev and
# ~/muscle-build/ref-f3559af are only read. Each binary gets a BUILD_INFO.txt beside it: the base
# build's BUILD_INFO, the variant, the source of cloak.cpp and the SHA-256.
set -euo pipefail
MB="$HOME/muscle-build"
DEV="$MB/dev"
JOBS=8
VARIANT=${1:-}
[ -f "$DEV/BUILD_INFO.txt" ] || { echo "no build in $DEV (run build.sh)" >&2; exit 1; }

info() {  # BASE_DIR VARIANT_TEXT BIN OUT_FILE
  { sed 's/^/base   /' "$1/BUILD_INFO.txt"
    echo "variant:   $2"
    echo "built:     $(date -Is)"
    echo "sha256:    $(sha256sum "$3" | cut -d' ' -f1)"
  } > "$4"
  cat "$4"
}
fresh_copy() {  # SRC DST: DST becomes an exact copy of SRC
  rm -rf "$2"; cp -a "$1" "$2"
}
same_base() {  # BASE COPY: true if COPY was made from BASE's current build
  [ -f "$2/BASE_BUILD_INFO.txt" ] && cmp -s "$1/BUILD_INFO.txt" "$2/BASE_BUILD_INFO.txt"
}
swap_build() {  # TREE CPP MAKEFILE: replace cloak.cpp (CR removed) and rebuild incrementally
  tr -d '\r' < "$2" > "$1/src/cloak.cpp"
  ( cd "$1/src" && make -f "$3" -j"$JOBS" > variant_build.log 2>&1 ) ||
    { tail -20 "$1/src/variant_build.log"; exit 1; }
}

case "$VARIANT" in
asan)
  T="$MB/asan"
  fresh_copy "$DEV" "$T"
  cp "$DEV/BUILD_INFO.txt" "$T/BASE_BUILD_INFO.txt"
  cd "$T/src"
  sed -e 's#^CXXFLAGS := .*#CXXFLAGS := -O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined -D_GLIBCXX_ASSERTIONS -fopenmp -DNDEBUG#' \
      -e 's#^LDFLAGS := .*#LDFLAGS := -g -fsanitize=address,undefined -fopenmp#' \
      -e '/LDFLAGS += -static/d' -e '/^\tstrip /d' Makefile > Makefile.asan
  rm -rf o ../bin/muscle
  make -f Makefile.asan -j"$JOBS" > variant_build.log 2>&1 || { tail -20 variant_build.log; exit 1; }
  info "$DEV" "asan (Makefile.asan)" "$T/bin/muscle" "$T/BUILD_INFO.txt"
  ;;
mutant|mutant-asan)
  ID=${2:?mutant ID}; CPP=$(readlink -f "${3:?mutant cloak.cpp}")
  if [ "$VARIANT" = mutant ]; then BASE="$DEV"; T="$MB/mutants/tree"; MK=Makefile; OUTD="$MB/mutants/$ID"
  else BASE="$MB/asan"; T="$MB/mutants/tree-asan"; MK=Makefile.asan; OUTD="$MB/mutants/$ID-asan"
    [ -f "$BASE/BUILD_INFO.txt" ] || { echo "build the asan variant first" >&2; exit 1; }
  fi
  if ! same_base "$BASE" "$T"; then
    mkdir -p "$MB/mutants"; fresh_copy "$BASE" "$T"; cp "$BASE/BUILD_INFO.txt" "$T/BASE_BUILD_INFO.txt"
  fi
  swap_build "$T" "$CPP" "$MK"
  mkdir -p "$OUTD"; cp "$T/bin/muscle" "$OUTD/muscle"
  info "$BASE" "$VARIANT $ID, cloak.cpp from $CPP" "$OUTD/muscle" "$OUTD/BUILD_INFO.txt"
  ;;
seenset)
  CPP=$(readlink -f "${2:?seen-set cloak.cpp}")
  BASE="$MB/ref-f3559af"; T="$MB/seenset"
  [ -x "$BASE/bin/muscle" ] || { echo "no f3559af build in $BASE" >&2; exit 1; }
  [ -d "$T" ] || fresh_copy "$BASE" "$T"
  swap_build "$T" "$CPP" Makefile
  { echo "base:      $BASE (f3559af, $("$BASE/bin/muscle" -version 2>&1 | head -1))"
    echo "variant:   seenset, cloak.cpp from $CPP"
    echo "built:     $(date -Is)"
    echo "sha256:    $(sha256sum "$T/bin/muscle" | cut -d' ' -f1)"
  } > "$T/BUILD_INFO.txt"
  cat "$T/BUILD_INFO.txt"
  ;;
portable)
  T="$MB/portable"
  fresh_copy "$DEV" "$T"
  cd "$T/src"
  sed -i 's/ -march=native//' Makefile
  ! grep -q 'march=native' Makefile || { echo "-march=native still in the Makefile" >&2; exit 1; }
  rm -rf o ../bin/muscle
  make -j"$JOBS" > variant_build.log 2>&1 || { tail -20 variant_build.log; exit 1; }
  info "$DEV" "portable (-march=native removed from CXXFLAGS and LDFLAGS)" "$T/bin/muscle" "$T/BUILD_INFO.txt"
  ;;
*)
  sed -n '2,19p' "$0"; exit 2 ;;
esac
