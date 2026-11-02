#!/bin/bash
# The byte oracle (tools/oracle_direct.py) behind MUSCLE's command line, standing in for the binary
# when the test harness itself is tested: every test script must pass on it (the positive control).
# Usage: bash oracle_bin.sh -cloak IN -output OUT [-mincol M]   (Debian; conda env muscle-cloak)
exec "$HOME/anaconda3/envs/muscle-cloak/bin/python" "$(dirname "$(readlink -f "$0")")/../tools/oracle_direct.py" "$@"
