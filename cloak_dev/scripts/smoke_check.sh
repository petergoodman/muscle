#!/bin/bash
# One-minute check that a released CLOAK binary works on the machine that will run it: its SHA-256
# is the released one, and it reproduces every fixture (scripts/test_quick.sh) and hand case
# (scripts/test_hand.sh) byte for byte, with the expected exit status. Needs only bash and
# coreutils: no Python, no conda, no compiler.
#
# Usage, from a clone of the fork at the release tag (git clone --branch cloak-vX.Y.Z ...):
#   bash cloak_dev/scripts/smoke_check.sh BINARY SHA256
# Prints PASS or FAIL. Outputs under ~/muscle-runs/.
set -uo pipefail
[ $# -eq 2 ] || { sed -n '2,10p' "$0"; exit 2; }
BIN=$1; WANT=$2
HERE=$(cd "$(dirname "$0")" && pwd)
got=$(sha256sum "$BIN" | cut -d' ' -f1)
if [ "$got" != "$WANT" ]; then echo "FAIL: SHA-256 is $got, expected $WANT"; exit 1; fi
ok=1
bash "$HERE/test_quick.sh" "$BIN" | tail -1 || ok=0
bash "$HERE/test_hand.sh" "$BIN" | tail -1 || ok=0
if [ $ok -eq 1 ]; then echo PASS; else echo FAIL; exit 1; fi
