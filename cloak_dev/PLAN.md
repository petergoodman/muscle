# Plan: CLOAK in MUSCLE

> **Status, 2026-11-01: phase 0 (setup) complete. Phase 1 is now the release `cloak-v0.3.0`: the
> speed rewrite of `ddc9df1` (gate A passed 2026-10-30) plus the corrections the reviews of
> 2026-10-30 and 2026-11-01 found, built test first. Steps 1 to 3 are done (records, test
> infrastructure, gate B0 passed on 2026-11-02); next is step 4, the code changes.** Peter owns
> this document.

This is the project's living core document: what we are building and why, the constraints, the
phases and their gates, open questions, where it stands, and the next step. It is read at the
start of every task and rewritten as the project moves. Decisions are cited by number from
`DECISIONS.md`; how CLOAK works and fails is in `ANALYSIS.md`. Project-specific details (users,
validation data, cluster paths, correspondence) stay out of the public repository, in
`private/PROJECT_NOTES.md` (gitignored).

## Goal

Make `muscle -cloak` usable on ensembles of up to a few thousand sequences, such as large Pfam
families, computing exactly the published algorithm (Mol Biol Evol 2026, doi
10.1093/molbev/msag182, Methods and Fig. 6), then turn it into a documented, tested command that
can be offered upstream. CLOAK stays a C++ command inside MUSCLE (decision 001).

## Constraints

1. Code changes stay in `src/cloak.cpp` and its three one-line hooks; Edgar's code is not edited.
2. Every release relates exactly to `cloak-v0.1.0` as decision 003 states, and is checked against
   the published definition.
3. Every output-changing behaviour change is approved by Peter before it is implemented, and the
   output must agree with the three references for what CLOAK computes: the published definition,
   `cloak.py`, and `cloak-v0.1.0` through decision 003's relation.
4. Validation data owned by others are read only; validation writes only to Peter's directories.
5. Existing test data are never changed to make a test pass; new test data may be added (decision
   009).

## Phases and gates

**Phase 0, setup (done 2026-10-29).** Fork synced (`dev`, `cloak`, tag `cloak-v0.1.0` at
`f3559af`), GitHub rulesets, frozen baseline binary, April 2025 smoke test identical.

**Phase 1, release cloak-v0.3.0.** Direct grouping (decision 004, in `ddc9df1`), every sequence
kept (007), `.` as a gap (008), every error exits 1 with no output, and output written with checked
I/O, each built test first: the failing tests are recorded before the code changes them.
1. Records for 0.3.0 committed; `cloak` pushed.
2. Test infrastructure: a checker written from the definition, independent oracles (including
   `cloak.py`'s own functions and a memoised build of `f3559af`), metamorphic relations, random
   generators with asserted coverage, hand fixtures (installed by Peter), and one entry point,
   `scripts/test.sh {edit|commit|release}`. Before any code change, the failing set on `ddc9df1`
   and `f3559af` is recorded and must match the prediction.
3. Gate B0 on `ddc9df1`, run locally on a read-only copy of the users' reference data (procedure
   in the private notes): every existing `cloak-v0.1.0` reference output reproduced byte for
   byte, and every ensemble that `cloak-v0.1.0` could not finish completed with non-empty output.
   On a pass, tag `cloak-v0.2.0` on `ddc9df1` (not released to users).
4. Code, one change per commit with its tests and records: rows and labels (007); errors and
   checked output; `.` gaps (008).
5. Quality pass: sanitizers, compiler warnings, mutation testing (every non-equivalent mutant
   killed), performance against `ddc9df1`.
6. Gate A: `build.sh --full` from a clean commit, then the portable build (no `-march=native`);
   `scripts/gate_a.sh` (release tier) passes on the portable binary with no skips.
7. Gate B1 on the portable binary, same data: every reference output meets decision 003's
   relation, every output passes the checker with all sequences present, and the ensembles that
   `cloak-v0.1.0` could not finish complete within 16 GB and 1 h (measured locally; anything above
   the local 5 GB cap is measured once on a larger machine).
8. Tag `cloak-v0.3.0` on the gated commit, push the tag, and hand the release to its users with
   the portable binary, its SHA-256, the commit, the gate tables, the known issues, and a
   one-minute smoke check (`scripts/smoke_check.sh`) to run once on their cluster.

**Phase 2, revisions (0.4.0).** The remaining open questions below, decided one by one, each with
its test, against the definition. Gate: all tiers pass; gate B rerun with the expected differences
only.

**Phase 3, release and upstream.** Tag the last revision and offer upstream a new branch from its
chosen base (`main` or `dev`, decision 002) that carries only the `src/` diff, with new commit
messages (decision 006).

## Open questions

1. Closed 2026-11-01: every sequence is kept (decision 007).
2. Closed 2026-11-01: `.` is a gap (decision 008).
3. Errors: 0.3.0 exits 1 with a message for invalid options, unreadable or empty input, allocation
   failure while loading, no surviving column and failed writes. Still open: `-mincol` values that
   MUSCLE's parser wraps (`ANALYSIS.md`, section 5, item 10).
4. No surviving column: exit 1 in 0.3.0; an empty output file instead remains possible.
5. Closed 2026-11-01: the column order is kept and documented (the paper does not specify it).
6. An optional map from output columns back to MSA 0 columns.
7. Usage text and README: default output name, `-mincol` meaning, `-efa_explode` syntax.
8. The log line that prints `userMincol`.
9. Memory: store classes sparsely if gate B shows the rows x output-columns term matters
   (`ANALYSIS.md`, section 7).

Questions for other people are kept in `private/PROJECT_NOTES.md`.

## Risks

- `-march=native`: a binary built on one CPU type may fail on another. The release binary is built
  without it, and its users run the smoke check once before use.
- The fixtures are subsampled Pfam ensembles, not production data; gate B covers that.
- The WSL VM has 6 GB: the baseline can exhaust it on unfamiliar input (`CLAUDE.md`, "Memory"), and
  a large ensemble may exceed the local cap in gate B.
- Peak memory has a rows x output-columns term (`ANALYSIS.md`, section 7).
- The baseline and the candidate share MUSCLE's loader, so loader defects are invisible to
  comparisons between them; the checker models the loader separately.
- The Windows build (`muscle.vcxproj`) is not tested.

## Current state and next step

`ddc9df1` holds the speed rewrite, passed gate A on 2026-10-30 and gate B0 on 2026-11-02; its binary
is kept at `~/muscle-build/ref-ddc9df1/`, and Peter tags it `cloak-v0.2.0`. The test infrastructure
is committed (`tests/README.md`), with the red step recorded in `CHANGELOG.md`. Next step: the rows
commit (step 4), then errors and checked output, then `.` gaps.
