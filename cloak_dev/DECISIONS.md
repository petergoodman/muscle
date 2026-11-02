# Decisions

The active rules for this project's code, tests and process. Each entry states its rule as it
stands now; a superseded rule is rewritten in place and its "Since" line updated, and git keeps
the earlier wording. Numbers are permanent and never reused. An entry is added only when Peter
chooses between real alternatives and the choice constrains later work (`CLAUDE.md`, "Records").

Entry format: **Rule**, **Why**, **Not chosen**, **Since**; at most 12 lines.

### 001. CLOAK stays a C++ command inside MUSCLE

**Rule.** CLOAK is developed only as `muscle -cloak` (`src/cloak.cpp` and its three one-line hooks).
No separate implementation is maintained or released; Python in `cloak_dev/tools/` is test
tooling only.
**Why.** MUSCLE is permanent and widely available, CLOAK needs MUSCLE's 16-MSA ensembles anyway,
and the code is almost self-contained, so it can move with MUSCLE's code.
**Not chosen.** A standalone CLOAK (Python or C++); further work on `cloak.py` in
`phylowheeler/CLOAK`.
**Since** 2026-10-29.

### 002. Work on the fork's `cloak` branch, based on upstream `dev` at `f3559af`

**Rule.** All work is committed to `cloak` of `petergoodman/muscle`, which starts at `f3559af`
(`cloak-v0.1.0`). `dev` mirrors upstream `dev` and is never committed to; `main` is left as it
was; nothing is pushed upstream.
**Why.** Every difference from `cloak-v0.1.0` is then CLOAK's own change.
**Not chosen.** Basing `cloak` on upstream `main` (13 commits Edgar has not merged into `dev`,
including build fixes). The two CLOAK commits move onto it with two one-line conflicts and give
identical output (verified 2026-10-28), so this is revisited when the change set is offered
upstream (phase 4).
**Since** 2026-10-29.

### 003. Releases and their exact relation to cloak-v0.1.0

**Rule.** `cloak-v0.2.0` (`ddc9df1`, speed only) is byte-identical to `cloak-v0.1.0` wherever 0.1.0
handles the input correctly. `cloak-v0.3.0` computes the published definition (007, 008) and relates
to 0.1.0 exactly: for input with only `-` gaps and no label starting with `>`, its output minus the
rows with at most one residue is 0.1.0's output byte for byte (0.1.0 exits 1 if no row remains); its
output on input X equals its output on X with every `.` made `-`; where 0.1.0 exits 0 without output
or crashes, it exits 1. Pass lines: a checker written from the definition, plus these relations.
**Why.** Every 0.1.0 output stays explainable, and byte-level tests survive the intended changes.
**Not chosen.** Byte identity for 0.3.0 (keeps the row bug, 007); comparing column multisets only.
**Since** 2026-11-01 (first version 2026-10-29: 0.2.x changes speed and memory only).

### 004. Retained columns are found by direct grouping

**Rule.** For each column of MSA 0, residues are grouped by the tuple of columns they occupy in MSAs
1 to K-1; groups of at least `-mincol` residues are emitted in lexicographic tuple order, which is
the order of `f3559af`'s output (`ANALYSIS.md`, section 4).
**Why.** It computes the same columns in time about linear in the input. Peak memory also has a
term proportional to rows times output columns (`ANALYSIS.md`, section 7).
**Not chosen.** A per-column set of already expanded subsets, and an early exit from the
recursion. Both remove the duplicates and should be byte-identical, but both keep the
linear search for each residue's column in every MSA and the recursion's dense copies: the early
exit took 134 s against 3.9 s for direct grouping on a synthetic 1,000 x 1,500 ensemble.
**Since** 2026-10-29 (Why corrected 2026-11-01).

### 005. Python runs only in the conda environment `muscle-cloak`

**Rule.** Scripts in `cloak_dev/tools/` run in `muscle-cloak` (`cloak_dev/environment.yml`), never
in the Windows or Debian system Python. Exception: muscle's own build runs `vcxproj_make.py` with
`/usr/bin/python3`, as upstream does.
**Why.** Peter's standing rule that the system Python stays untouched, as in the IQ-TREE fork.
**Not chosen.** Debian's system Python, which would work for these standard-library scripts.
**Since** 2026-10-29.

### 006. Personal details stay in the fork; upstream receives only the code diff

**Rule.** The fork's tracked files may name Peter and describe his machine and runs, but hold
nothing about other people, their projects or results, cluster paths or private correspondence
(`CLAUDE.md`, "Git"). What is offered upstream is a new branch from upstream's chosen base carrying
only CLOAK's `src/` diff, with new commit messages; never `cloak` itself or its history.
**Why.** The fork stays Peter's working repository, and nothing personal or fork-specific reaches
`rcedgar/muscle`. The 0.2.0 code diff (`src/cloak.cpp` only) has no fork references (checked
2026-10-30).
**Not chosen.** Removing Peter's own details from the fork now; tracking only code information.
**Since** 2026-10-30.

### 007. Every input sequence is in the output

**Rule.** CLOAK writes every input sequence, under its label exactly as read, in label order, however
many residues it keeps; a sequence that keeps none is an all-gap row. Only columns are filtered, and
every output column holds at least `-mincol` residues.
**Why.** The published algorithm removes columns, never sequences, and `cloak.py` keeps every row.
Deleting rows after the column filter, as 0.1.0 does, also leaves columns with one residue or none
(`ANALYSIS.md`, section 5, item 1).
**Not chosen.** Deleting rows with at most one residue (0.1.0 and 0.2.0); an option to choose.
**Since** 2026-11-01.

### 008. Gaps are `-` and `.`

**Rule.** CLOAK treats exactly `-` and `.` as gaps, as MUSCLE's `isgap` and its ensemble loader do;
every other non-space character is a residue. Output gaps are written as `-`.
**Why.** The loader accepts an ensemble only if its sequences agree with both characters removed, so
numbering residues the same way makes residue p of a row the same residue in every alignment. Pfam
Stockholm alignments use `.` (the fixture recipe's `sto2fa.py` removes it); MUSCLE writes only `-`.
**Not chosen.** Numbering `.` as a residue (0.1.0, silently wrong); rejecting input that has `.`.
**Since** 2026-11-01.
