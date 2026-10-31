#!/bin/bash
# Fail if any term in cloak_dev/private/forbidden_terms.txt appears in a tracked file (staged or
# working-tree version) or in the message of a commit on HEAD that no origin branch contains
# (CLAUDE.md, "Git", "Public content"). Run before every commit and push.
#
# Usage (Git Bash, Windows git; never Linux git in this working copy):
#   bash cloak_dev/scripts/check_public.sh
#
# Terms are case-insensitive fixed strings, one per line; blank lines and lines starting with #
# are ignored. Exit status 0 only if nothing matches; 2 if there are no terms to check.
set -uo pipefail
cd "$(git rev-parse --show-toplevel)" || exit 2
TERMS=cloak_dev/private/forbidden_terms.txt
PAT=$(mktemp); trap 'rm -f "$PAT"' EXIT
[ -f "$TERMS" ] && tr -d '\r' < "$TERMS" | grep -v -e '^#' -e '^[[:space:]]*$' \
  | tr '[:upper:]' '[:lower:]' > "$PAT"
[ -s "$PAT" ] || { echo "no terms in $TERMS" >&2; exit 2; }

bad=0
echo "== staged files"
git grep --cached -I -i -n -F -f "$PAT" && bad=1
echo "== working-tree files"
git grep -I -i -n -F -f "$PAT" && bad=1
echo "== messages of commits not on origin"
# Lowercased instead of grep -i: Git Bash's grep 3.0 aborts on -i -F -f.
for c in $(git rev-list HEAD --not --remotes=origin); do
  git log -1 --format=%B "$c" | tr '[:upper:]' '[:lower:]' | grep -n -F -f "$PAT" \
    | sed "s/^/${c:0:9}:/" && bad=1
done
[ $bad -eq 0 ] && echo "clean" || echo "forbidden terms found"
exit $bad
