# Plan: CLOAK in MUSCLE

> **Status, 2026-10-30: phase 0 (setup) complete; phase 1 (hotfix cloak-v0.2.0): the code commit
> `ddc9df1` passed gate A; the project files are scrubbed for the public fork and staged, not
> committed; next is their commit, the push and gate B.** Peter owns this document.

This is the project's living core document: what we are building and why, the constraints, the
phases and their gates, open questions, where it stands, and the next step. It is read at the
start of every task and rewritten as the project moves. Decisions are cited by number from
`DECISIONS.md`; how CLOAK works and fails is in `ANALYSIS.md`. Project-specific details (users,
validation data, cluster paths, correspondence) stay out of the public repository, in
`private/PROJECT_NOTES.md` (gitignored).

## Goal

Make `muscle -cloak` usable on ensembles of up to a few thousand sequences, such as large Pfam
families, then turn it into a documented, tested command that can be offered upstream. CLOAK
stays a C++ command inside MUSCLE (decision 001).

## Constraints

1. Code changes stay in `src/cloak.cpp` and its three one-line hooks; Edgar's code is not edited.
2. `cloak-v0.2.x` output is byte-identical to `cloak-v0.1.0` (decision 003).
3. Every output-changing behaviour change in 0.3.0 is approved by Peter, and by the users whose
   results it changes, before it is implemented.
4. Validation data owned by others are read only; validation writes only to Peter's directories.
5. Expected outputs in `tests/data/` are never edited to make a test pass.

## Phases and gates

**Phase 0, setup (done 2026-10-29).** Fork synced (`dev`, `cloak`, tag `cloak-v0.1.0` at
`f3559af`), GitHub rulesets, frozen baseline binary, April 2025 smoke test identical.

**Phase 1, hotfix cloak-v0.2.0.** Direct grouping replaces the recursion (decision 004).
1. Peter reviews the diff and commits in two commits: `src/cloak.cpp`, then the project files.
2. Gate A (laptop, Debian): `build.sh --full` from the clean commit; `test_quick.sh
   --check-baseline` passes all cases; `test_random.sh 200` reports no difference.
3. Push `cloak`. Gate B (cluster; procedure and data in the private notes): every existing
   `cloak-v0.1.0` reference output is reproduced byte for byte, and every ensemble that
   `cloak-v0.1.0` could not finish completes within 16 GB and 1 h.
4. Tag `cloak-v0.2.0` on the gated commit, push the tag, and hand the release to its users with the
   binary path, commit, SHA-256, gate tables and the known issues (`ANALYSIS.md`, section 5; row
   dropping first).

**Phase 2, specification and tests.** A short specification of what CLOAK computes and of every
open question below, decided one by one; a slow checker written from the definition; tests for
each decided behaviour; automated tests on GitHub if Peter wants them. Gate: Peter approves the
specification, and the users whose outputs change approve those parts.

**Phase 3, revisions (0.3.0).** One change per commit, each with its test, against the
specification. Gate: all tests pass; gate B rerun with the expected differences only.

**Phase 4, release and upstream.** Tag `cloak-v0.3.0` and offer upstream a new branch from its
chosen base (`main` or `dev`, decision 002) that carries only the `src/` diff, with new commit
messages (decision 006).

## Open questions (to decide in phase 2)

1. Rows left with at most one residue: keep them all (as `cloak.py` does), drop them (as now), or
   an option.
2. `.` gaps: support them or reject them with a clear message (0.2.0 already stops).
3. Errors: report on stderr and exit nonzero for invalid options and unreadable input.
4. No surviving column: empty output file or an error.
5. Column order: keep the current order, or document it.
6. An optional map from output columns back to MSA 0 columns.
7. Usage text and README: default output name, `-mincol` meaning, `-efa_explode` syntax.
8. The log line that prints `userMincol`.

Questions for other people are kept in `private/PROJECT_NOTES.md`.

## Risks

- `-march=native`: a binary built on one CPU type may fail on another; cluster binaries are built
  on the cluster's compute nodes.
- The fixtures are subsampled Pfam ensembles, not production data; gate B covers that.
- The WSL VM has 6 GB: the baseline can exhaust it on unfamiliar input (`CLAUDE.md`, "Memory").

## Current state and next step

`ddc9df1` holds the 0.2.0 code and passed gate A on 2026-10-30 (`CHANGELOG.md`, top entry). The
project files are scrubbed for the public fork, pass `scripts/check_public.sh`, and are staged, not
committed; `cloak` is not yet pushed. Next step: Peter commits the project files, reruns
`check_public.sh`, and pushes `cloak` (permanent: force pushes are blocked), then gate B (phase 1,
step 3).
