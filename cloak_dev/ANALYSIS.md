# Analysis of CLOAK at f3559af

How `src/cloak.cpp` worked at upstream `f3559af` (cloak-v0.1.0), why it exhausted memory, why
direct grouping gives identical output, and the other problems found. Written 2026-10-29 from the
investigation of 2026-10-27 to 2026-10-29 (`CHANGELOG.md`); corrected and extended 2026-11-01 from
the reviews of 2026-10-30 and 2026-11-01. This file describes and measures; what
to change is decided in `DECISIONS.md` and planned in `PLAN.md`.

Labels: (verified) means read in source or run during the investigation; (measured) gives the
number a kept script printed; (estimated) is arithmetic, not a run.

## 1. What CLOAK computes

Given K MSAs of the same sequences, two residues stay in one output column if and only if they
share a column in all K MSAs. "Share a column in MSA k" is an equivalence relation on residues, so
the intersection over k is one too, and its classes are the output columns. Every class lies
inside one column of MSA 0. Classes with fewer than `-mincol` residues (default 2, the published
algorithm) are dropped, and the result is written as FASTA with rows sorted by label. The set of
classes does not depend on which MSA is MSA 0; only the column order does.

Up to `cloak-v0.2.0` the code also drops every row left with at most one residue (section 5, item
1) and numbers `.` as a residue (item 2). `cloak-v0.3.0` keeps every row and treats `.` as a gap
(decisions 007 and 008).

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
| PF00188_n20 (20 seqs, stratified) | 101 s, 0.60 GB; 3,569,796 calls; 3,088,544 stored entries (3,088,350 from the recursion, 194 unsplit columns), 340 distinct from the recursion | 0.02 s |
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

Assumptions: the loader has accepted the ensemble, n >= 1, and every row has the same number of
residues in every MSA (P2). The loader guarantees P2 when the code's gap characters are the
loader's (`-` and `.`); `f3559af` numbers `.`, so for it P2 holds only for input without `.`
(section 5, item 2). Under P2 each residue number p of row r names the same residue in every MSA.

- Every set the recursion works on is an intersection of columns and so a union of classes.
  `can_this_divvy_be_found_without_subsequent_divvies` holds exactly when a block lies in one column
  of every MSA (its check against MSA 0 is vacuous), so every set the recursion outputs is a whole
  class of at least two residues. A column of MSA 0 that no MSA splits is output whole: one class,
  possibly of one residue, which the `-mincol` filter removes.
- Every class of at least two residues is reached: at any set strictly larger than the class, some
  MSA splits the set, the block containing the class has at least two residues (so the size filter
  in `divvies_as_3darray` cannot drop it) and is strictly smaller, and the loop has no early exit,
  so by induction on set size the recursion ends at the class.
- Order: let i* be the first MSA that splits a set. For i < i* the only block equals the set and is
  skipped by the `toSearch == divvy` test; at i* the blocks are visited in ascending column, and
  each block's subtree outputs all its classes before the next block starts; later iterations only
  repeat classes. All classes of the set agree on MSAs 1..i*-1, so by induction first occurrences
  follow lexicographic key order.
- Output entries from different MSA-0 columns are distinct unless both are empty, so the final
  de-duplication keeps each column's classes in that order. Every residue maps to a letter other
  than `-`, so filtering by `-mincol` before conversion to letters removes the same columns as
  `remove_singletons` after it.

Tests (measured): byte-identical output to `f3559af` on 13 Pfam-derived ensembles and Peter's
April 2025 test, at `-mincol` 2, 3, 5 and 20 where run; `scripts/test_random.sh 400 1` (800
synthetic ensembles, `-mincol` 2 and 3, so 1,600 comparisons, 709 with non-empty output) with no
difference in output or exit status, and another 1,600 comparisons under a slightly different
size schedule likewise; a single-MSA EFA (identical); and an ensemble where nothing survives (both
exit 1, no output). On the four ensembles `f3559af` cannot
finish, the candidate, an early-exit variant of the recursion and `tools/oracle_direct.py` agree
byte for byte. Not kept, so not rerunnable: the `-mincol` 20 runs, the other size schedule, the
single-MSA EFA and the early-exit variant. `oracle_direct.py` uses the same key and sort as the
candidate, so it does not test the order argument. The 2026-10-30 audit built `f3559af` with a
per-column set of already expanded subsets and an order-preserving final de-duplication (which
removes only repeats); that build reproduced all 23 expected outputs byte for byte, including the
four that `f3559af` cannot finish.

## 5. Other problems found (verified unless marked)

Status: (0.3.0) fixed in `cloak-v0.3.0` (`PLAN.md`, phase 1); (open) listed in `PLAN.md`, "Open
questions"; (MUSCLE) in MUSCLE's own code, documented and not changed (constraint 1). Every item
holds for both `f3559af` and `ddc9df1` unless stated.

1. **Rows are dropped, and columns can then fall below `-mincol`.** After filtering columns,
   `cmd_cloak` drops every row with at most one residue left (`letterCount > 1`, independent of
   `-mincol`). The published algorithm and `cloak.py` keep all rows. In the PF00353full_n1226_s5
   fixture 24 of 1,226 sequences disappear. Because the column filter runs first, dropping a row
   can remove a residue that a column needed: on a 4-row ensemble (rows `A--`, `CD-`, `-E-`, `--F`
   in two identical MSAs) both binaries write only the row `CD`, two columns of one residue each,
   and on another 4-row ensemble they write an all-gap column (verified 2026-11-01; both inputs
   become hand fixtures). None of the 23 fixture outputs has such a column. (0.3.0, decision 007)
2. **`.` is a gap to the loader but a residue to CLOAK.** `Ensemble` accepts an ensemble only if
   its sequences agree with `-` and `.` removed, but `convert_to_nums` numbers `.`, so residue p of
   a row can be different residues in different MSAs. `f3559af` then writes wrong output with
   exit 0 (`.` gaps in MSA 1 of a 20-sequence ensemble: 90 columns became 167), reads out of
   bounds and can crash when MSA 0 has more `.` than another MSA, and is right by coincidence in
   some cases, such as `.` after a row's last residue in MSAs other than MSA 0. `ddc9df1` stops
   only when a row's residue count differs between MSAs: with equal counts it prints `.` as a
   residue as `f3559af` does, and it stops on cases that `f3559af` gets right (in the 2026-10-30
   audit, 242 of 1,926 unequal-count comparisons). MUSCLE writes only `-`.
   (0.3.0, decision 008)
3. **Errors that exit 0 or crash.** A `-mincol` below 2 is reported only through `Log()`, which
   writes nothing unless `-log` is given, and `cmd_cloak` returns normally: exit 0, no output file.
   An empty value (`-cloak ''`) does the same. An allocation failure while loading is caught by
   the `catch (...)` around `E.FromFile`: exit 0, no output, no message; after loading it is
   uncaught, exit 134. On PF00353full_n1226_s5 under `ulimit -v`, caps of 16 to 45 MiB give exit 0
   with no output and 64 and 90 MiB give exit 134 (`scripts/test_oom.sh` on `20798ce`,
   2026-11-02). A path list that names no MSA segfaults (exit 139). A missing input path and
   `-cloak` without a value already exit 1. (0.3.0: each exits 1 with a message, and an allocation
   failure anywhere is caught)
4. **Nothing survives.** A result with no column left stops with `Empty FASTA` (exit 1) rather than
   writing an empty file. Once every row is kept, an empty result would instead reach
   `MSA::ToFASTAFile` with no columns, where `(0 - 1)/60 + 1` wraps and the loop writes out of
   bounds, so 0.3.0 stops before building the output. (open: exit 1 kept in 0.3.0)
5. **Log line.** The start message prints `userMincol` (0 when `-mincol` is absent) instead of the
   value used. (open)
6. **Documentation.** The default output name is `cloak_result.fasta`, while the README promises
   `{input}.cloak.fa`; `-efa_explode` takes `-prefix`/`-suffix`, not `-output`. (open)
7. **Single-threaded.** Columns of MSA 0 are independent, so a parallel loop is possible; with
   direct grouping it is unlikely to matter. (open)
8. **Writes are not checked.** The output goes through MUSCLE's `TextFile`, whose `PutChar` and
   `PutString` check `fputc` and `fputs` only with `assert` (removed by `-DNDEBUG`) and whose
   `Close` ignores `fclose`'s result (`textfile.cpp`, `textfile.h`). A full disk or quota therefore
   leaves a truncated file with exit 0: under a 50 KiB file-size limit with `SIGXFSZ` ignored,
   `20798ce` exits 0 and leaves 51,200 bytes, and `-output /dev/full` exits 0
   (`scripts/test_errors.sh`, 2026-11-02). (0.3.0 writes the file itself, checks every write and
   the close or flush, and on failure removes a partial regular file and exits 1. MUSCLE's checked
   `WriteStdioFile` and `CloseStdioFile` `Die` without removing the partial file, and
   `CreateStdioFile` would create a file named `-`, which `TextFile` writes to standard output.)
9. **Labels lose a leading `>`.** `cmd_cloak` strips a leading `>` from every label, although the
   loader's labels never include the FASTA `>`, so an input label `>x` (line `>>x`) is written as
   `x`. (0.3.0, decision 007)
10. **`-mincol` overflow.** MUSCLE's `StrToUint` wraps without an error: `-mincol 4294967296` reads
    as 0 (so the default, 2) and `4294967297` as 1 (an error). (MUSCLE)
11. **Path lists drop some characters.** MUSCLE's FASTA reader drops `*`, `?` and `~` with a
    warning, which leaves the rows of an MSA unequal in length; the loader then stops, or, since
    `MSA::AppendSeq` checks lengths only on every fourth sequence, overflows the heap. EFA input
    keeps these characters (2026-10-30 audit). (MUSCLE)
12. **An MSA with fewer rows than MSA 0 crashes the loader.** `Ensemble::ToUpper` loops over MSA
    0's row count in every MSA before `SortMSAs` compares the counts, so it reads past the rows of
    a smaller MSA: an EFA with rows a, b in MSA 0 and a in MSA 1 exits 139 on both binaries
    (verified 2026-11-02). A later MSA with more rows reaches the check and exits 1. A truncated
    EFA can produce this. (MUSCLE)

## 6. Differences from cloak.py (phylowheeler/CLOAK d846718, verified)

The same columns, usually in a different order. `cloak.py` takes its reference alignment and MSA
order from directory listing order (`os.walk`), and orders blocks by first appearance (dictionary
insertion) rather than by column number. With both made deterministic, its output equals
`f3559af`'s row for row on three ensembles. Besides section 5, items 1 and 2 (`cloak.py` keeps
every row, as 0.3.0 does, and numbers `.` as `f3559af` does), `cloak.py` keeps lower case, appends
to an existing result file on a rerun, fails if its temporary directory exists, raises IndexError
on an all-gap column of its reference alignment, and collapses sequences whose headers share the
text between `|` and the next `-`.

## 7. Memory of direct grouping (measured, 2026-10-30 audit of `ddc9df1`; not rerun)

A least-squares fit over 33 runs (18 fixtures and 15 synthetic ensembles) gives peak RSS = 5.3 MB
+ 25.9 B per input cell (rows x columns, summed over MSAs) + 11.8 B per row x output column; every
run above 30 MB was predicted within 7%. The second term comes from storing each retained column as
a vector over all rows (8 bytes per entry) and copying it as letters, as `f3559af` also does. It
grows with rows twice, since output columns grow with rows: at 500 columns, 1,000, 2,000 and 4,000
rows took 249, 595 and 1,548 MB. A 2,500 x 2,500 x 16 synthetic ensemble with local disagreement
took 15.8 s and 2,698 MiB for 7,742 output columns; one whose classes are all pairs (n = 3,000,
45,000 output columns) stopped with `std::bad_alloc` under a 1.43 GiB cap. Storing each class as
its (row, residue) pairs would make the second term proportional to the residues written. (open)
