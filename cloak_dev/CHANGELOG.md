# Changelog

Session history for the CLOAK work on the fork `petergoodman/muscle`, newest first. Each entry is
about ten lines: what was done, what failed and why, what is next. An entry that closes a gate
also holds its gate record. The next session starts from the top entry.

## 2026-11-03: release cloak-v0.3.0 (step 8)

- **Binary frozen:** the gated portable binary copied to `~/muscle-release/cloak-v0.3.0/`
  (read-only, with `BUILD_INFO.txt` and `.sha256`) and to a private share folder; both SHA-256
  `f461271d...a6fa467`.
- **Gate B rerun in bytes mode** on the frozen copy (`gate_b.sh ... bytes`,
  `~/muscle-runs/gate_b/20261103-101830-bytes/`): reference IDENTICAL 58, failed FINISHED 19,
  failures 0; all 77 outputs byte-identical to B1's. 0.3.0 reproduces every 0.1.0 reference exactly.
- **Public check:** `check_public.sh` covers HEAD's files and unpushed messages only, so the files
  of all 11 unpushed commits and the tag message were also scanned with its terms: clean.
- **Tag** `cloak-v0.3.0`, annotated, on `790cd12`, the gated commit; its message is the release
  notes (11 known issues), identical to the approved draft. The gate records are in this commit.
- **Rehearsal** from a Linux clone of the local tag
  (`~/muscle-runs/release-rehearsal/20261103-104113/`):
  - the share copy passes `sha256sum -c` and `smoke_check.sh` (23 of 23, 42 of 42, PASS); a
    relative BINARY path FAILs, as known issue 10 says;
  - a build from the tag with `build_linux.bash` reads `[790cd12]` (SHA-256 `397f2774...`), passes
    the smoke check, and `gate_b.sh ... bytes` gives IDENTICAL 58, FINISHED 19
    (`20261103-104551-bytes`), all 77 outputs identical to the release binary's.
- **Next:** Peter pushes `cloak` and the tag and sends the release to its users; a fresh GitHub
  clone is then checked; then the audit's harness fixes.

## 2026-11-02: pre-release audit of 790cd12

- **Result:** no defect in `src/cloak.cpp` or the portable binary (`f461271d...`). Every code
  change since `ddc9df1` is C1, C2 or C3; test data were only added; gate records match their logs.
- **Independent oracle** (written from the paper, no code shared with `tools/`; scripts and outputs
  in `~/muscle-runs/audit-20261102/`): 9,500 generated cases (5,000 general, 1,500 heavy, 1,500
  with `-log` on this binary and on the asan build), 0 failed; 2,531 also checked against `f3559af`
  and 2,707 against `cloak.py` run end to end. Fixtures and gate B ensembles byte-identical; gate
  B references identical to 0.1.0's without stripping. `test.sh commit` passes; the smoke check
  passes from a Linux clone.
- **Found:** `smoke_check.sh` FAILs falsely unless BINARY is an absolute path; relation runs ignore
  exit status, `check_fixtures.sh` lacks `--strict` and skips MANIFESTs, `gate_a.sh` does not pin the
  harness; `AI_DISCLOSURE.md` has no 0.3.0 entry.
- **Failed:** my first oracle expected output at `-mincol 1` (23 false failures); rerun clean.
- **Next:** step 8 (recommended: tag `790cd12`; users pass an absolute path to the smoke check),
  then the harness fixes.

## 2026-11-02: gate B1 passed on 790cd12

- **Gate record:** commit `790cd12`; portable release binary SHA-256 `f461271d...a6fa467` (the
  gate A1 binary); `gate_b.sh ... strip`; output `~/muscle-runs/gate_b/20261102-225557-strip/`.
- **Result:** every reference output meets decision 003's relation; every output passes
  `spec_check.py --order` with every sequence present; every ensemble 0.1.0 could not finish
  completes, all runs under 5 s and 1.4 GB. The smoke ensemble matches `cloak.py` (`d846718`,
  through `ref_cloakpy.py`). Tables and the reproducer comparison are in the private notes.
- **Failed:** the plan expected 5 incomplete ensembles left out; the script reports 0, correctly,
  since those `.efa.tmp` files have no `.efa` beside them. Background `wsl.exe` completion notices
  arrived before the WSL jobs had finished; results were read from the run directories.
- **Next:** the release (PLAN.md, step 8). First decide whether the tag goes on `790cd12`, whose
  binary was gated, or on the commit that carries the gate records.

## 2026-11-02: gate A1 passed on 790cd12

- **Gate record:**
  - commit `790cd12`;
  - `build.sh --full`: `src state: clean`, dev binary `4c3bb6a3d9337900159b5b98e15518bb371ea3789a951f45548940393218d6d2`;
  - portable release binary (`build_variant.sh portable`, `~/muscle-build/portable/bin/muscle`):
    SHA-256 `f461271d39d342c1c1f42164673238121cd3444f09e50023ec0caf4baa6fa467`;
  - output `~/muscle-runs/gate_a/20261102-190956/` (`GATE_A.txt`).
- **Result:** `gate_a.sh` passed the release tier, 0 of 16 parts failed, in 44 min:
  - random-big: 5,600 comparisons and 6,562 relation checks, 0 failed, 0 skipped;
  - mutants: 34, all as required;
  - performance within thresholds of `ddc9df1`.
- **Fixed in records:** the earlier mutant count (27 killable before the M07 relabel, not 28).
- **Next:** gate B1 (PLAN.md, step 7) on this portable binary.

## 2026-11-02: quality pass (Q1 to Q6)

- **Built:**
  - `build_variant.sh` (asan, mutant, mutant-asan, seenset, portable);
  - `test_sanitize.sh`, `warn_diff.sh`, `make_mutants.py` + `test_mutants.sh`, `test_perf.sh`,
    `tag_runs.py`, `make_seenset.py` + `test_seenset.sh`, `gate_a.sh`;
  - the release tier in `test.sh`, with no `todo` part left.

  Every part has a positive control that fails as required.
- **Run once on `e0eda41`** (numbers in `tests/README.md`):
  - sanitizers clean;
  - no new compiler warning;
  - performance within thresholds of `ddc9df1`;
  - tagging 23 of 23;
  - O5 4 of 4;
  - the portable build byte-identical on the hand and fixture cases.
- **Mutants:** all 27 killable mutants are killed (24 by the edit tier, 3 by the commit tier), and
  the control and 4 equivalents survive. The Q3 commit message says 28 and 25, which is a miscount.
- **Failed:** M07 (filtering at 2 in `create_output`), recorded as equivalent by the 2026-10-30
  audit, is killed by E15a, E15c and the tiny `-mincol 3` cases. Since C1, the no-survivor check
  reads `create_output`'s result, so that filter is not redundant (ANALYSIS 5.4). The shipped code
  is unaffected.
- **Decided** (Peter): M07 becomes killable (killers E15a and E15c), and X04 (skipping
  `remove_singletons`, the redundant filter) is added as an equivalent mutant. Rerun: M07 is
  killed, and X04 survives every tier.
- **Not run:** Valgrind, which is not installed and which Peter agreed to skip.
- **Next:** gate A1, run on the commit with the M07 change.

## 2026-11-02: `.` gaps (C3)

- **Code** (`src/cloak.cpp`, decision 008):
  - `convert_to_nums` and `produce_output_in_letters` use MUSCLE's `isgap`, so `.` is a gap.
  - The residue-count `Die` is an `asserta`, since the loader guarantees equal counts.
  - Comments and a header note are updated.
- **Red step:** the C2 failing set (5 hand cases, 20 tiny and 143 random comparisons, all `.`).
- **Green:** `test.sh commit` passes with no failure: hand 42 of 42; random 800 comparisons (680
  against `f3559af`, 289 against `cloak.py`) and 1,513 relation checks, 0 failed, 0 skipped.
- **Next:** the quality pass (PLAN.md, step 5).

## 2026-11-02: errors and checked output (C2)

- **Decided** (Peter, approving the C2 plan):
  - CLOAK writes its FASTA itself, because Edgar's stdio helpers leave a partial file and treat
    `-` as a file name. A failed write removes a partial regular file.
  - An allocation failure anywhere exits 1.
  - `test_oom.sh` is not run on the byte oracle.
- **Code**, in `src/cloak.cpp`:
  - `Die` for an empty `-cloak` value, for `-mincol` below 2, and for a path list naming no MSA.
  - `write_fasta` checks every write and the close; `-output -` still writes to stdout.
  - `cmd_cloak` turns `std::bad_alloc` into exit 1.
- **Tests:** `scripts/test_errors.sh` (E20a to E20c, E21, H17; edit tier) and `scripts/test_oom.sh`
  (E16; commit tier). `oracle_direct.py` now handles failed writes and stdout the same way.
- **Red step** on `20798ce` (`fed24aa7...`), as predicted:
  - E01 exits 0, E04 exits 139.
  - E20a to E20c and E21 exit 0; E20c leaves 51,200 bytes.
  - Memory caps of 16 to 45 MiB exit 0 with no output; 64 and 90 MiB exit 134.
- **Green:**
  - hand 37 of 42; errors 5 of 5; oom 15 of 15; quick 23 of 23;
  - tiny and random fail the same 20 and 143 comparisons as `20798ce`, all of them `.` cases;
  - the positive control passes.
- **Failed:** I edited `test.sh` while a run was executing it, and bash hit a syntax error after
  the run's last part. The results of the parts that had finished were kept and used.
- **Next:** C3 (`.` gaps), then the quality pass.

## 2026-11-02: test infrastructure, gate B0 on ddc9df1

- **Decided** (Peter). Agents may add test data but never change it (decision 009, guarded by
  `SHA256SUMS`). Constraint 3 now asks for agreement with the paper, `cloak.py` and 0.1.0, not
  user approval. Timing uses `tools/runmeasure.py`. A local pre-push hook refuses everything except
  `cloak`, `dev` and `cloak-v*` tags.
- **Written.**
  - Checker and loader model: `spec_check.py` with its tests, `cloakio.py`.
  - Oracles and cases: `ref_cloakpy.py` (`cloak.py`'s own functions), `oracle_direct.py` (0.3.0),
    `make_hand_fixtures.py` (42 hand cases, in `tests/data/hand/`), `gen_cases.py`, `run_suite.py`.
  - Scripts: `test.sh`, `test_hand.sh`, `gate_b.sh`, `smoke_check.sh`, `check_fixtures.sh`, and
    `test_quick.sh` and `test_random.sh` reworked. The design is in `tests/README.md`.
- **Verified** (numbers in `tests/README.md`).
  - Positive control: the harness on the byte oracle passes every hand and random case. That
    includes 680 random comparisons with `f3559af` through decision 003's relation and 289 with
    `cloak.py`.
  - `cloak.py` agrees with every applicable hand case.
- **Red step.** Hand cases: both binaries fail 16 of 42, as predicted. Random cases: 388 and 413
  of 800 comparisons fail, every one explained by short rows or `.`.
  `test.sh commit` on the oracle fails only on PF00353full_n1226_s5, whose MANIFEST row still
  expects 0.1.0's output.
- **Gate B0 passed** on `ddc9df1` (binary SHA-256 `c683e3b5...e2af1a`): every reference output
  reproduced byte for byte, and every previously failing ensemble finished. Tables are in the
  private notes.
- **Found.** An MSA with fewer rows than MSA 0 crashes MUSCLE's loader (`ANALYSIS.md`, section 5,
  item 12).
- **Rows commit (C1)**, in `src/cloak.cpp`:
  - Every sequence is written, under its exact label.
  - `Die` before output when no column survives.
  - With Peter's approval, the PF00353full_n1226_s5 MANIFEST row now points at its 0.3.0 expected
    output, which the C1 build reproduces byte for byte.
  - Results: `test_quick.sh` 23 of 23. Hand cases 35 of 42; the 7 still failing are the `.` cases
    and E01 and E04. The tiny random cases now fail only relation R3 (inserted `.`).
- **Next.**
  - Peter tags `cloak-v0.2.0` on `ddc9df1`.
  - The errors and checked-output commit (C2), then `.` gaps (C3), following the detailed plan in
    `private/PLAN_v030.md`.

## 2026-11-01: reviews, third audit, plan for cloak-v0.3.0

- **Reviews.** Read the 2026-10-30 audit of `ddc9df1` and the 2026-11-01 review against the
  published definition (both moved to `private/reviews/`), then audited again: the grouping core
  matches the definition, the code around it does not (`ANALYSIS.md`, section 5).
- **Found.**
  - Dropping rows after the column filter leaves columns with one residue or none (both binaries,
    two 4-row ensembles).
  - Writes are unchecked, so a failed write exits 0.
  - A leading `>` is stripped from labels.
  - Duplicate labels already stop in the loader (verified; that review item is closed).
- **Decided** (Peter):
  - Fix the seven review items and release as `cloak-v0.3.0`.
  - Every sequence is kept (decision 007), `.` is a gap (008), and decision 003 is rewritten.
  - Gate B runs locally on a copy of the users' data, and the release binary is portable.
  - Tag `ddc9df1` as `cloak-v0.2.0` after gate B0.
- **Records committed and `cloak` pushed:** DECISIONS 003, 004, 007, 008; PLAN; ANALYSIS sections
  1, 3, 4, 5, 6, 7; two CLAUDE.md lines. The detailed plan is in `private/PLAN_v030.md`. The gate A binary
  of `ddc9df1` is kept at `~/muscle-build/ref-ddc9df1/` (SHA-256 verified).
- **Measured:** a read-only scan of the 23 expected outputs. 22 keep every row, and PF00353full_n1226_s5
  loses 24 of 1,226. None has a column below its `-mincol`.
- **Next.** The test infrastructure (PLAN phase 1, step 2). The users' reference data must be
  copied to WSL for gate B0.

## 2026-10-30: public-repository scrub, gate A on ddc9df1

- **History.** Peter rebuilt `cloak` on `f3559af` with the 0.2.0 code alone (`8ee5bdc`), then its
  message was amended to drop a session link, giving `ddc9df1`; the tree is unchanged and `src/`
  equals that of `948ecec`. The earlier commits are on the local branch `backup-unscrubbed`, never
  pushed.
- **Scrubbed** the project files for the public fork: project-specific passages moved to
  `cloak_dev/private/PROJECT_NOTES.md` and the cluster jobs to `cloak_dev/private/slurm/`, both
  gitignored; added `scripts/check_public.sh`, the rule in `CLAUDE.md`, "Git", and decision 006.
  `check_public.sh` passes. Git Bash's grep 3.0 aborts on `grep -i -F -f`, so the script
  lowercases the messages instead.
- **Gate A passed on `ddc9df1`** (laptop, Debian 12, 2026-10-30).
  - Build: `build.sh --full`, `src state: clean (equals the commit)`, `-version`
    `muscle 5.3.linux64 [ddc9df1]`; binary SHA-256
    `c683e3b59e1b1f9663f75ff158d2a654312a4e63873e058649a17d7180e2af1a`.
  - `test_quick.sh --check-baseline`: `passed 40, failed 0` (23 candidate and 17 baseline cases);
    `~/muscle-runs/test_quick/20261030-165627/`.
  - `test_random.sh 200`: `compared 800 (non-empty output 345), skipped 0 (baseline over its cap),
    differing 0`; `~/muscle-runs/test_random/20261030-165651/`.
- **Earlier gate A run on `948ecec`**, 2026-10-29, with the same results; superseded because that
  commit is not published. Command lines were not logged, so the options are read from the
  outputs.
  - Build: `build.sh --full`, commit `948ecec`, `src state: clean (equals the commit)`, `-version`
    `muscle 5.3.linux64 [948ecec]`; binary SHA-256
    `947df45897709e9321105a2d31ea92f68bec7e0b34825168736067cde4847c7f` (rechecked 2026-10-30).
  - `test_quick.sh --check-baseline`: `passed 40, failed 0` (23 candidate and 17 baseline cases);
    `~/muscle-runs/test_quick/20261029-221054/`.
  - `test_random.sh 200`: `compared 800 (non-empty output 345), skipped 0 (baseline over its cap),
    differing 0` (seeds 1 to 200 mixed, 1001 to 1200 mild);
    `~/muscle-runs/test_random/20261029-221115/`.
- **Next.** Commit the project files, push `cloak` (permanent: force pushes are blocked), then
  gate B (`PLAN.md`, phase 1).

## 2026-10-29 (Cowork, second session): hotfix code, project files

- **Written** (uncommitted, for Peter's review): `CLAUDE.md`, `.claude/settings.json`,
  `.claude/.gitignore`, Peter's gitignored `.claude/settings.local.json`, and `cloak_dev/` (plan,
  decisions, analysis, scripts, tools, 23 fixture cases, cluster jobs).
- **Changed** `src/cloak.cpp`: `create_output` now groups each MSA-0 column's residues by their
  columns in MSAs 1..K-1 (decision 004); the recursion, its helpers and the duplicate-removal pass
  are removed. Rows whose residue counts differ between MSAs now stop with an error (decision 003).
- **Tested in the Cowork sandbox** (Ubuntu 24.04, g++ 13.3): byte-identical to `f3559af` on all
  fixtures and on `test_random.sh 400 1` (1,600 comparisons, none differing); the four fixtures
  `f3559af` cannot finish run in under 0.4 s (`ANALYSIS.md`, section 3). `build.sh` and
  `test_quick.sh` were exercised against a simulated CRLF working copy. Nothing yet run on the
  laptop.
- **Next.** Peter reviews and commits in two commits, then gate A (`PLAN.md`, phase 1).

## 2026-10-29: setup (Peter, steps 1 to 6)

- Debian toolchain present: g++ 12.2.0, GNU Make 4.3, git 2.39.5, Python 3.11, static `libc.a`,
  `libgomp.a`, `libstdc++.a`.
- Fork: cloned to `C:\Users\peter\BioInformatics\muscle`; `upstream` added with push disabled;
  `dev` (tracks `upstream/dev`) and `cloak` created at `f3559af`; tag `cloak-v0.1.0` at `f3559af`;
  `dev`, `cloak` and the tag pushed; default branch `cloak` (these verified from GitHub's public
  refs); rulesets on `cloak`, `dev` and `cloak-v*` (reported; not publicly visible). `main`
  untouched at `ee8f8b9`.
- Baseline: `~/muscle-baseline/muscle-f3559af`, built with `build_linux.bash` in a fresh clone of
  upstream at `f3559af`; `-version` reads `muscle 5.3.linux64 [f3559af]`; SHA-256 beside it.
- Smoke test: the April 2025 ensemble in `muscle_testing` gives output byte-identical to
  `cloak_test_output.fasta` (exit 0).

## 2026-10-27 to 2026-10-28 (Cowork): investigation

- Read `cloak.cpp` at `f3559af` and `cloak.py` at `d846718`; built `f3559af` in a sandbox and
  confirmed the reported diagnosis with exact counts: the recursion re-derives each output
  column along every order of splitting and stores every copy (`ANALYSIS.md`, sections 2 and 3).
- Measured on Pfam-derived ensembles: 20 sequences took 101 s and 0.6 GB; 92 sequences exhausted
  a 5 GB cap; a 1,226-sequence super5 ensemble is predicted to need about 25,500 GB.
- Found and tested the direct grouping, its output order, and the other problems in `ANALYSIS.md`,
  section 5. Failed approach: none; the early-exit patch works but stays slow on wide ensembles.
- Verified that the two CLOAK commits move onto upstream `main` (two one-line conflicts, identical
  output on three ensembles).
