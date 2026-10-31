# Analysis of CLOAK at f3559af

How `src/cloak.cpp` worked at upstream `f3559af` (cloak-v0.1.0), why it exhausted memory, why
direct grouping gives identical output, and the other problems found. Written 2026-10-29 from the
investigation of 2026-10-27 to 2026-10-29 (`CHANGELOG.md`). This file describes and measures; what
to change is decided in `DECISIONS.md` and planned in `PLAN.md`.

Labels: (verified) means read in source or run during the investigation; (measured) gives the
number a kept script printed; (estimated) is arithmetic, not a run.

## 1. What CLOAK computes

Given K MSAs of the same sequences, two residues stay in one output column if and only if they
share a column in all K MSAs. "Share a column in MSA k" is an equivalence relation on residues, so
the intersection over k is one too, and its classes are the output columns. Every class lies
inside one column of MSA 0. Classes with fewer than `-mincol` residues (default 2) are dropped,
then rows (sequences) left with at most one residue are dropped, and the result is written as
FASTA with rows sorted by label. The set of classes does not depend on which MSA is MSA 0; only
the column order does.

## 2. How f3559af computes it (verified)

- `cmd_cloak` loads the ensemble through `Ensemble::FromFile` (labels and ungapped sequences
  checked, rows sorted to MSA 0's order, text uppercased), copies it to a 3-D char array sorted by
  label (`convertEnsembleTo3D`) and numbers residues per row from 1 (`convert_to_nums`; only `-`
  is a gap).
- `create_output` takes each column of MSA 0 and, through `create_indices`, finds each residue's
  column in every MSA by a linear scan of the row.
- If every MSA keeps the column's residues in one column (`divvy_is_required` false), the column is
  output whole. Otherwise `add_to_output_or_discard` splits it ("divvies") by each MSA i = 1..K-1
  in turn: each block of at least two residues that lies in one column of every MSA is output, and
  every other block is split again by the same procedure, starting from MSA 1.
- Each output entry is a vector of n `NumOrDash` (8 bytes each). Duplicates are removed only at
  the end, by comparing every entry with every distinct entry kept so far, then numbers are turned
  back into letters and the column and row filters applied.

## 3. Why it exhausts memory and time (verified, measured)

At every node of the recursion the loop tries all K-1 MSAs, not just one, as the next split. The
first MSA that splits a set already leads to every class inside it; the other branches reach the
same classes again along other orders of splitting, and every arrival is stored. The number of
paths grows roughly factorially with the number of MSAs that disagree on a column. The parameter
`seen` was apparently meant to prevent this but is never read or written; `cloak.py` has the same
unused `seen` (its line 176), so the port is faithful and the flaw is in the original design.

Cost therefore follows disagreement among the MSAs, not sequence count or width. An exact counter
(`tools/predict_recursion.py`, memoised on the residue set) reproduces the instrumented binary's
counts exactly on three ensembles (1,508 / 2,036; 88,402 / 114,763; 3,569,796 / 3,088,350 calls /
stored entries) and so predicts the cost of runs too large to make.

| Ensemble (fixture name) | f3559af | cloak-v0.2.0 candidate |
|---|---|---|
| PF00188_n20 (20 seqs, stratified) | 101 s, 0.60 GB; 3,569,796 calls; 3,088,544 stored entries for 340 distinct | 0.02 s |
| PF00188_n60 (60 seqs, stratified) | 120 s, 1.34 GB | 0.32 s |
| PF00188_n92_s5 (92 seqs, super5) | stopped by a 5 GB cap after 1,297 s (`std::bad_alloc`); predicted 1.78e7 entries, 13.8 GB | 0.12 s, 25 MB |
| PF00353full_n1226_s5 (1,226 seqs, 134 columns) | predicted 1.34e9 calls and 2.59e9 entries for 1,635 distinct, about 25,500 GB | 0.33 s, 94 MB |
| synthetic, 1,000 seqs x about 1,500 columns | not run | 3.9 s, 717 MB (the early-exit variant: 134 s) |

Times are from a 2-core cloud sandbox with other jobs running, so treat them as orders of
magnitude. Memory is peak RSS of the process.

Two secondary costs remain even without duplicates: the linear scan in `create_indices`,
O(K n c^2) for n rows and c columns (about 1e11 comparisons at n = c = 2,000, estimated), and the
dense length-n vectors built for every block at every node.

## 4. Direct grouping gives the same output (verified)

Key each residue of a column of MSA 0 by the tuple of its columns in MSAs 1..K-1. Residues with
equal keys form one class, by the definition in section 1.

- Every set the recursion outputs is a block of some MSA's split of its parent that lies in one
  column of every MSA, so it is contained in one class; and the parent, being an intersection of
  columns containing it, contains that class, so the block is the whole class.
- Every class of at least two residues is reached: at any set strictly larger than the class, some
  MSA splits the set, the block containing the class has at least two residues (so the size filter
  in `divvies_as_3darray` cannot drop it) and is strictly smaller, so the recursion ends at the
  class. A column that needs no split is a single class either way.
- Order: blocks are visited in `std::map` order, ascending column number. Within a block of MSA 1
  the recursion skips MSA 1 (it cannot split the block) and splits by MSA 2 in ascending order, and
  so on, so classes are first reached in lexicographic key order, and everything after the first
  splitting MSA at each node only repeats them. Keeping first occurrences (as the final
  de-duplication does) therefore yields lexicographic key order.

Tests (measured): byte-identical output to `f3559af` on 13 Pfam-derived ensembles and Peter's
April 2025 test, at `-mincol` 2, 3, 5 and 20 where run; `scripts/test_random.sh 400 1` (800
synthetic ensembles, `-mincol` 2 and 3, so 1,600 comparisons, 709 with non-empty output) with no
difference in output or exit status, and another 1,600 comparisons under a slightly different
size schedule likewise; a single-MSA EFA (identical); and an ensemble where nothing survives (both
exit 1, no output). On the four ensembles `f3559af` cannot
finish, the candidate, an early-exit variant of the recursion and `tools/oracle_direct.py` agree
byte for byte.

## 5. Other problems found (verified unless marked)

Status: (0.2.0) handled in the hotfix; (open) listed in `PLAN.md`, "Open questions".

1. **Rows are dropped.** After filtering columns, `cmd_cloak` drops every row with at most one
   residue left (`letterCount > 1`, independent of `-mincol`). `cloak.py` keeps all rows. In the
   PF00353full_n1226_s5 fixture 24 of 1,226 sequences disappear. (open)
2. **`.` gaps give wrong output without an error.** `Ensemble` treats `.` as a gap when it checks
   sequences, but `convert_to_nums` numbers `.` as a residue. Replacing `-` by `.` in one MSA of a
   20-sequence ensemble changed the output from 90 to 167 columns with exit 0. The 0.2.0 candidate
   stops with an error instead, because residue counts then differ between MSAs. MUSCLE itself
   writes `-`. (0.2.0 stops; open for proper support)
3. **Errors that exit 0.** A `-mincol` below 2 or a missing input path is reported only through
   `Log()`, which writes nothing unless `-log` is given, and `cmd_cloak` returns normally: exit 0,
   no output file. The `try`/`catch` around `E.FromFile` never fires, since `Die` exits. (open)
4. **Nothing survives.** A result with no column left stops with `Empty FASTA` (exit 1) rather than
   writing an empty file. (open)
5. **Log line.** The start message prints `userMincol` (0 when `-mincol` is absent) instead of the
   value used. (open)
6. **Documentation.** The default output name is `cloak_result.fasta`, while the README promises
   `{input}.cloak.fa`; `-efa_explode` takes `-prefix`/`-suffix`, not `-output`. (open)
7. **Single-threaded.** Columns of MSA 0 are independent, so a parallel loop is possible; with
   direct grouping it is unlikely to matter. (open)

## 6. Differences from cloak.py (phylowheeler/CLOAK d846718, verified)

The same columns, usually in a different order. `cloak.py` takes its reference alignment and MSA
order from directory listing order (`os.walk`), and orders blocks by first appearance (dictionary
insertion) rather than by column number. With both made deterministic, its output equals
`f3559af`'s row for row on three ensembles. Besides item 1 above, `cloak.py` keeps lower case,
appends to an existing result file on a rerun, fails if its temporary directory exists, and
collapses sequences whose headers share the text between `|` and the next `-`.
