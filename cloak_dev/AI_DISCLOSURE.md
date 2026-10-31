# AI disclosure

Record of substantive contributions by AI tools to this work, kept so that the Methods or
Acknowledgements of any resulting manuscript, its cover letter, and any release notes can state
them accurately. Newest first. One entry per milestone, at most 15 lines: tool and model,
operator, what it did, files touched. Commits made with an AI agent carry a `Co-Authored-By` line.
No AI-generated result in this project is accepted without the tests and gates in `PLAN.md`.

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
