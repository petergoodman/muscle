# MUSCLE fork: CLOAK

Fork of `rcedgar/muscle` (C++, built by Edgar's `vcxproj_make.py`). The project repairs and then
revises CLOAK, the `muscle -cloak` command, which keeps two residues in one output column only if
they share a column in every MSA of an ensemble (normally the 16 MSAs of `-stratified` or of
4 x `-super5 -perm all`). CLOAK is `src/cloak.cpp` plus one line each in `src/cmds.h`,
`src/myopts.h` and `src/muscle.vcxproj`. Goal, scope, phases and next step: `cloak_dev/PLAN.md`;
where it is silent, ask.

Project facts only. Personal working preferences live in `~/.claude/CLAUDE.md`.

## Read these before starting

Every session:

1. `cloak_dev/PLAN.md`: goal, constraints, phases and gates, open questions, current state, next
   step.
2. The top entry of `cloak_dev/CHANGELOG.md`: where the last session stopped.

Only when the task needs them:

3. `cloak_dev/DECISIONS.md`, the active rules: all of it before planning a change.
4. `cloak_dev/ANALYSIS.md`, before touching `src/cloak.cpp`: how CLOAK worked at `f3559af`, why it
   exhausted memory, why direct grouping gives identical output, and the other problems found.
5. `cloak_dev/private/PROJECT_NOTES.md` (gitignored), for anything involving other people, their
   data, the cluster or upstream: the project-specific details kept out of the public files.

Not read at the start: `cloak_dev/AI_DISCLOSURE.md` and older CHANGELOG entries.

## Environment

- **One working copy**, `C:\Users\peter\BioInformatics\muscle`, edited in VS Code and managed only
  with Windows git (`git.exe`). Edgar's `.gitattributes` (`* text=auto`) makes Windows check text
  out with CRLF, so Linux git in WSL sees false changes here: never run Linux git `status`, `diff`
  or any command that writes in this folder (`build.sh`'s read-only `rev-parse` and `clone` are
  fine). Do not create a second working copy.
- **Builds and runs happen in WSL2 Debian 12** (g++ 12.2, GNU Make 4.3, Python 3.11 at
  `/usr/bin/python3`), on the Linux filesystem, never in the working copy. `build.sh` keeps a
  Linux clone in `~/muscle-build/dev` for this (below); it is a build area, never edited.
- **Frozen baseline**: `~/muscle-baseline/muscle-f3559af` (SHA-256 in `muscle-f3559af.sha256`
  beside it), built 2026-10-29 from upstream `f3559af` with `build_linux.bash` in a fresh clone.
  `-version` reads `muscle 5.3.linux64 [f3559af]`. Never rebuilt, overwritten or moved.
- **Python** for `cloak_dev/tools/` runs only in the conda environment `muscle-cloak`
  (`cloak_dev/environment.yml`, standard library only), never in the Windows or Debian system
  Python: `~/anaconda3/bin/conda run --no-capture-output -n muscle-cloak python ...`. The one
  exception is muscle's own build, which runs `vcxproj_make.py` with `/usr/bin/python3` as upstream
  does.
- **Memory**: the WSL VM has 6 GB. Run the baseline only on the fixtures marked `fast` or `slow`
  in `cloak_dev/tests/data/MANIFEST.tsv`, or under a cap: `( ulimit -v 2000000; timeout 120 ... )`.
  On an unfamiliar ensemble, `cloak_dev/tools/predict_recursion.py` counts the baseline's work
  exactly without running it.
- **Long jobs need the laptop awake**: start
  `powershell.exe -NoProfile -ExecutionPolicy Bypass -File C:/Users/peter/.claude/keepawake/keepawake.ps1 PATTERN`
  in the background, as in the IQ-TREE project. Peter's gitignored `.claude/settings.local.json`
  allows it.

## Commands

```bash
# Debian terminal, any directory; conda off first (conda deactivate)
bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/build.sh            # incremental
bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/build.sh --full     # clean build
bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/test.sh edit       # after every build, ~1 min
bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/test.sh commit     # before every commit, ~5 min
bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/gate_a.sh BINARY   # gate A: release tier, ~3 h
bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/build_variant.sh portable  # release binary; also asan, mutant, seenset
bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/test_quick.sh       # 23 fixtures
bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/test_hand.sh        # 42 hand cases
bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/test_random.sh 200  # 800 random comparisons
bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/test_errors.sh      # failed writes, stdout
bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/test_oom.sh         # memory caps, ~6 s
# binary: ~/muscle-build/dev/bin/muscle; BUILD_INFO.txt beside bin/ records commit and src state
# test design, oracles, tiers and measurements: cloak_dev/tests/README.md
```

- `build.sh` checks out the working copy's HEAD in `~/muscle-build/dev`, mirrors the tracked
  `src/` files into it with CRLF removed (uncommitted edits included), and builds there.
  `--full` (and the first run) builds exactly as `build_linux.bash` does and stamps the commit into
  `-version`; an incremental build leaves the stamp at the last full build. A gate build is `--full`
  from a clean commit, and its `BUILD_INFO.txt` must say `src state: clean`.
- `test_quick.sh` byte-compares the candidate with the expected outputs; `--check-baseline` also
  reruns the baseline on the fast cases. `test_hand.sh` runs the hand cases (exact exit status;
  `--predict v010|v020` records the red step). `test_random.sh` and `test.sh` run random ensembles
  through every oracle (`tools/run_suite.py`). All exit nonzero on any failure or skip, and all
  start with `check_fixtures.sh`.
- Release-tier parts, each runnable alone: `test_sanitize.sh` (after `build_variant.sh asan`),
  `warn_diff.sh`, `test_mutants.sh` (~30 min), `test_perf.sh`, `tools/tag_runs.py`,
  `test_seenset.sh`. Variants build in their own `~/muscle-build/<variant>/`; `dev` is only read.
  `gate_a.sh` refuses a binary that is not from a `--full` build of a clean commit.
- `gate_b.sh DATA LIST BINARY bytes|strip` is gate B on the users' data (paths and pass lines in
  the private notes); `smoke_check.sh BINARY SHA256` is the users' one-minute check of a release.
- Timing and peak memory: `tools/runmeasure.py` (WSL has no `/usr/bin/time`).
- **Running WSL from Claude Code's Bash tool** (Git Bash): put multi-step commands in a script and
  run `MSYS_NO_PATHCONV=1 wsl.exe -d Debian -- bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/build.sh`.
  Without `MSYS_NO_PATHCONV=1`, Git Bash rewrites `/mnt/c/...`; `$VAR` inside
  `wsl.exe ... bash -c '...'` is expanded by the outer shell. Shells started this way do not
  activate conda.
- Run outputs go under `~/muscle-runs/`, never into the working copy.
- Cluster validation (gate B): `cloak_dev/private/slurm/README.md` (gitignored).

## Layout

- `src/cloak.cpp`: all of CLOAK. `cmd_cloak()` parses options, loads the ensemble (`Ensemble`,
  `ensemble.cpp`, which validates labels and ungapped sequences and sorts rows), calls
  `create_output()`, converts numbers to letters, drops columns under `-mincol` and rows with at
  most one residue, and writes FASTA.
- `cloak_dev/`: everything this fork adds besides CLOAK's code (the counterpart of `docs/agent/` in
  the IQ-TREE fork), so that a pull request to Edgar can leave it out in one step.
  - `scripts/` build and test drivers; `tools/` Python helpers; `tests/data/` fixtures, added but
    never changed (decision 009, `SHA256SUMS`, its `README.md`); `private/` project notes,
    staging and cluster jobs, gitignored.
- Not this project's: the rest of `src/`, `test_scripts/`, `test_data/`, `.github/` (Edgar's
  workflows, all manual-only).

## Conventions

- Code changes stay inside `src/cloak.cpp` and its three one-line hooks. Edgar's code is not
  edited.
- Errors: `Die()` (message to stderr, exit 1). `Log()` writes only to the file given by `-log`.
- Gaps are `-` and `.`, as in the `Ensemble` loader (decision 008, from the `.` commit of 0.3.0;
  earlier CLOAK numbered `.` as a residue, ANALYSIS.md); output gaps are `-`.
- Versions are tags `cloak-vMAJOR.MINOR.PATCH` on `cloak`. `cloak-v0.2.0` (`ddc9df1`) is the
  speed rewrite, byte-identical to `cloak-v0.1.0`; `cloak-v0.3.0` computes the published
  definition and relates to 0.1.0 exactly (decision 003).

## Gotchas

- The Makefile compiles with `-march=native -static`: a binary is built on the machine type that
  runs it. Cluster binaries are built on the cluster's compute nodes (`cloak_dev/private/slurm/`).
- `vcxproj_make.py` (downloaded by `build_linux.bash`, pinned commit `806d016`) always cleans and
  rebuilds everything and writes `gitver.txt` from `git log`, which is why builds need a git clone.
- Edgar's `cloak.cpp` starts with a UTF-8 byte-order mark; keep it.
- `-efa_explode` takes `-prefix` and `-suffix`, not `-output`.
- `muscle -cloak` accepts an EFA or a text file of MSA paths (CRLF line endings are fine).
- From the errors commit of 0.3.0, every error (an invalid `-mincol`, no surviving column, a failed
  write, running out of memory) exits 1 with a message on stderr and leaves no output file; before
  it, an invalid `-mincol` exits 0 with no output and no message unless `-log` is given.
  `-output -` writes to standard output.

## Git

- `origin` is the fork `petergoodman/muscle`. `upstream` is `rcedgar/muscle`, and its push URL is
  deliberately disabled (`DISABLED_NO_PUSH_TO_UPSTREAM`): do not re-enable it.
- **Never contact Edgar or anyone else.** No agent opens issues or pull requests, comments, or
  emails. Peter handles all communication; questions for others go to
  `cloak_dev/private/PROJECT_NOTES.md`.
- **Stage files by name**, never `git add -A` or `git add .`: local data folders in the working
  copy are kept out of git only by `.git/info/exclude`. A local pre-push hook (`.git/hooks/pre-push`,
  not tracked) refuses every push except `cloak`, `dev` and `cloak-v*` tags to `origin`, and runs
  `check_public.sh`.
- **Public content.** Tracked files and commit messages must contain nothing about collaborators,
  their projects or results, HPC paths, or private correspondence; that material lives only in
  `cloak_dev/private/` (gitignored, never pushed). Run `cloak_dev/scripts/check_public.sh` (Git
  Bash) before every commit and push; its terms are in `cloak_dev/private/forbidden_terms.txt`.
- Work happens on `cloak`. **Never commit to `dev`** (exact mirror of upstream `dev`, `f3559af` on
  2026-10-29) **or `main`** (stale, left as it was). GitHub rulesets block force pushes and
  deletion on `cloak` and `dev`, and updates or deletion of `cloak-v*` tags.
- `cloak-v0.1.0` = `f3559af`. Push a tag by name, never with `--tags`: the
  fetch from upstream also brings Edgar's release tags.
- To sync `dev`: `git switch dev; git merge --ff-only upstream/dev; git push origin dev;
  git switch cloak`. Bringing `dev` into `cloak` is a decision, not a routine step.
- Commit messages carry measurements and commands. Commits made with an agent carry a
  `Co-Authored-By` line. Pushes go through Git Credential Manager from Peter's machine.

## Records

Each kind of fact has one home; other files cite it.

| Home | Holds |
|---|---|
| `cloak_dev/PLAN.md` | goal, scope, constraints, phases and gates, open questions, risks, current state, next step |
| `cloak_dev/DECISIONS.md` | the active rules, edited in place, at most 12 lines per entry |
| `cloak_dev/ANALYSIS.md` | how CLOAK works and fails, with measurements; corrected in place, never prescriptive |
| `cloak_dev/CHANGELOG.md` | what each session did, what failed, what is next; gate records |
| `cloak_dev/AI_DISCLOSURE.md` | one entry per milestone |
| `cloak_dev/private/PROJECT_NOTES.md` | users, their data and results, cluster paths, correspondence, questions for others; never tracked |
| commit messages | measurements, commands, detail |
| this file | environment, commands, conventions, gotchas, git |

- **Decision bar.** A new decision only when Peter chooses between real alternatives and the choice
  constrains later code, tests or method. Corrections, measurements and status go to their homes
  above. A superseded decision is rewritten in place and its "Since" line updated; numbers are
  never reused.
- **Gate runs** use a `--full` build of a clean commit, and their records name the commit, the
  binary's SHA-256 and the output directory; a cluster gate's tables and paths go in the private
  notes.
