# AI disclosure

Record of substantive contributions by AI tools to this work, kept so that the Methods or
Acknowledgements of any resulting manuscript, its cover letter, and any release notes can state
them accurately. Newest first. One entry per milestone, at most 15 lines: tool and model,
operator, what it did, files touched. Commits made with an AI agent carry a `Co-Authored-By` line.
No AI-generated result in this project is accepted without the tests and gates in `PLAN.md`.

## 2026-10-30 to 2026-11-03: cloak-v0.3.0

**Tool:** Claude Code (VS Code extension) on Peter's laptop, building and running in WSL.
**Model:** `claude-opus-5-5` (Claude Opus 5.5, per the commits' `Co-Authored-By` lines).
**Operator:** Peter Goodman. All of this milestone's work, reviews included, ran in Claude Code.

**What it did.**
- **Reviews.** The 2026-10-30 audit of `ddc9df1`, the 2026-11-01 review against the paper (deleted
  rows, `.` as a residue, errors that exit 0, unchecked writes; `ANALYSIS.md`, section 5), and the
  2026-11-02 pre-release audit with an oracle written from the paper.
- **Records.** Drafted decisions 003, 007, 008 and 009, the 0.3.0 plan and the corrected analysis.
- **Tests.** Wrote the checker, oracles, random generators, hand fixtures (installed by Peter), the
  test tiers, the quality tooling and `gate_a.sh`, `gate_b.sh` and `smoke_check.sh`.
- **Code.** Wrote C1 to C3 in `src/cloak.cpp` (rows, errors and checked output, `.` gaps), test
  first.
- **Gates and release.** Ran gates B0, A1 and B1; drafted the release notes and the users'
  instructions.

**Files.** `src/cloak.cpp`; `CLAUDE.md`; `cloak_dev/` (scripts, tools, tests, records).

## 2026-10-27 to 2026-10-29: diagnosis, hotfix code and project setup

**Tool:** Claude in Cowork mode (Claude desktop app, cloud sandbox with read and write access to
folders Peter connected). **Model:** configured as `claude-opus-5-5`; the serving model may differ.
**Operator:** Peter Goodman.

**What it did.**
- **Diagnosis.** Read `cloak.cpp` (`f3559af`) and `cloak.py` (`d846718`), built `f3559af`, measured
  the failure on Pfam-derived ensembles and wrote an exact cost counter (`ANALYSIS.md`).
- **Code.** Wrote the direct-grouping `create_output` in `src/cloak.cpp` for cloak-v0.2.0.
- **Tests.** Generated the fixtures and their expected outputs (`tests/data/README.md` gives each
  output's source), the random-ensemble test and the cross-checking oracle.
- **Process.** Drafted the setup commands Peter ran, and wrote `CLAUDE.md`, `.claude/settings.json`
  and `cloak_dev/`.

**Files.** `src/cloak.cpp`; `CLAUDE.md`; `.claude/`; all of `cloak_dev/`.
