#!/bin/bash
# Verify that no existing test fixture has changed (decision 009). Every file under
# cloak_dev/tests/data except README.md, the MANIFEST.tsv files and SHA256SUMS itself must be
# listed in SHA256SUMS, and every listed file must match its checksum. A new fixture is added by
# appending its line (sha256sum PATH >> SHA256SUMS, from tests/data); an existing line is changed
# only with Peter's explicit approval, for a specification change recorded in DECISIONS.md.
#
# Usage (Debian or Git Bash): bash cloak_dev/scripts/check_fixtures.sh
# Exit status 0 only if every fixture is listed and unchanged. The test scripts run it first.
set -uo pipefail
cd "$(dirname "$0")/../tests/data" || exit 2
[ -f SHA256SUMS ] || { echo "check_fixtures: no SHA256SUMS" >&2; exit 2; }
bad=0
sha256sum --quiet -c SHA256SUMS || bad=1
unlisted=$(find . -type f ! -name README.md ! -name MANIFEST.tsv ! -name SHA256SUMS \
  | sed 's#^\./##' | sort | comm -23 - <(cut -c67- SHA256SUMS | sort))
if [ -n "$unlisted" ]; then
  echo "check_fixtures: files not in SHA256SUMS:"; echo "$unlisted" | sed 's/^/  /'; bad=1
fi
[ $bad -eq 0 ] && echo "check_fixtures: $(wc -l < SHA256SUMS) fixture files listed and unchanged"
exit $bad
