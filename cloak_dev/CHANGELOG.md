# Changelog

Session history for the CLOAK work on the fork `petergoodman/muscle`, newest first. Each entry is
about ten lines: what was done, what failed and why, what is next. An entry that closes a gate
also holds its gate record. The next session starts from the top entry.

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
