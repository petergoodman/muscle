# Tests

The tests decide whether a CLOAK binary computes what was published. Three references define that
(`PLAN.md`, constraint 3):
- the definition in the paper (Mol Biol Evol 2026, doi 10.1093/molbev/msag182, Methods and Fig. 6);
- `cloak.py` (phylowheeler/CLOAK `d846718`), the original implementation;
- `cloak-v0.1.0` (`f3559af`), through decision 003's relation: CLOAK's output without its rows of
  at most one residue equals 0.1.0's output on the same input with `.` made `-`.

Where they differ, the paper and `cloak.py` win: both keep every sequence (decision 007), and
`.` is a gap (decision 008). Every change is test first: a new test must fail on the code before
the change, and the failing set is recorded in `CHANGELOG.md`.

## Running them

```bash
# Debian terminal; conda env muscle-cloak (cloak_dev/environment.yml)
bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/test.sh edit     # after every build, about 1 min
bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/test.sh commit   # before every commit, about 5 min
bash /mnt/c/Users/peter/BioInformatics/muscle/cloak_dev/scripts/test.sh release  # gate A, hours
```

The binary defaults to `~/muscle-build/dev/bin/muscle`. Passing `scripts/oracle_bin.sh`, the byte
oracle behind MUSCLE's command line, tests the harness itself: every tier must pass on it, except
`test_oom.sh`, which `test.sh` reports as `n/a` there because a capped Python interpreter fails at
startup with other exit statuses.

| Tier | Parts |
|---|---|
| edit | `check_fixtures.sh`; `tools/test_spec_check.py`; `test_hand.sh`; `test_errors.sh`; 30 tiny random cases through every oracle and relation |
| commit | edit; `test_oom.sh`; `test_quick.sh --check-baseline`; `run_suite.py` on 300 `gen_cases.py` and 100 `gen_synth.py` cases, coverage asserted |
| release | commit; `--check-baseline-slow`; 2,000 + 800 more random cases; and, still to build (the tier fails until they exist): sanitizer builds, compiler warnings, mutation testing, performance against `ddc9df1`, digit-plane tagging, the memoised `f3559af` build on the large fixtures |

Gate B (`scripts/gate_b.sh`) runs a binary on the users' own ensembles and `cloak-v0.1.0`'s outputs
for them; `scripts/smoke_check.sh` is the one-minute check of a released binary on another machine.

## Oracles

| ID | Where | What it checks | Independence from `cloak.cpp` |
|---|---|---|---|
| O1 | `tools/spec_check.py` | the output meets the definition: every sequence, exact letters, columns equal to the classes (computed by tuple and, up to 200,000 pairs, pair by pair); exit 1 exactly when no class reaches `-mincol` | written from the paper's text |
| O2 | `tools/oracle_direct.py` | the output equals the byte oracle's, so order and format too | same algorithm, separate code |
| O3 | `tools/ref_cloakpy.py` | the same rows and columns as `cloak.py`'s own functions (sets; `-mincol` 2; inputs it can finish) | the original author's code |
| O4 | frozen baseline | decision 003's relation | MUSCLE's loader is shared |

All Python tools read input through `tools/cloakio.py`, a model of MUSCLE's loader (the EFA and
path-list readers, CR handling, upper case, label checks). The checker recovers which input
residue each output residue is from the forced order of the retained residues in each row, not from
its position in the output row. The handoff's earlier sketch made that mistake.

## Metamorphic relations (`tools/run_suite.py`)

R1 shuffled rows, R2 inserted all-gap columns, R3 `-` made `.`, R4 lower case, R5 CRLF and
wrapping, R6 path list, R10 a copied MSA, R14 reversed rows (output rows reversed), R15 the
ensemble concatenated with itself (output concatenated), R16 `-threads 8` and `-log`, R17 letters
renamed: same bytes. R8 rotation and R9 permutation of MSAs 1..K-1: same rows and columns as sets.
R11 a higher `-mincol`, R12 one MSA, R13 a deleted row: the output predicted from the original one.
R0 checks that the loader model's normalised ensemble gives the same output.

## Generated cases and coverage

`tools/gen_cases.py` derives MSAs 1..K-1 from MSA 0 by order-keeping operations (fission, fusion,
slide, rogue) under five schedules (fission, rogue, dots, kn, tiny), in three forms (EFA, CRLF EFA,
path list); `tools/gen_synth.py` keeps the earlier mixed and mild schedules. `run_suite.py
--assert-coverage` fails a run whose cases are too few in any class. The floors are about 60% of
the 2026-11-02 counts, per 300 cases:
- MSA-0 columns with two or more retained classes: 250 (`-mincol` 2) and 180 (3);
- such columns whose order needs an MSA after MSA 1: 200;
- cases with a sequence that keeps no residue: 120; one residue: 70;
- `.` in MSA 0: 15; `.` only in other MSAs: 10;
- one MSA: 10; 30 MSAs: 4; one sequence: 7;
- no surviving column: 50; path lists: 25; CRLF: 30.

## Hand fixtures

`tests/data/hand/` (42 cases; `tests/data/README.md`): tiny ensembles whose expected outputs were
written by hand. O1, O2 and O3 all agree with every one (`tools/make_hand_fixtures.py --verify
tests/data/hand private/reference/cloak.py`). `scripts/test_hand.sh --predict v010|v020` checks that
`f3559af` and `ddc9df1` fail exactly the predicted cases.

## Error cases outside the MANIFEST

`scripts/test_errors.sh` runs the cases a MANIFEST row cannot express, each with an exact exit
status:
- E20a, `-output /dev/full`: the error appears only at close.
- E20b, `-output -` on a full standard output.
- E20c, a 50 KiB file-size limit with `SIGXFSZ` ignored: no file may be left.
- E21, `-cloak ''`.
- H17, `-output -`: the expected bytes on standard output.

`scripts/test_oom.sh` (E16) runs PF00353full_n1226_s5 under 15 address-space caps, from 16 MiB to
2 GiB in steps of √2. Each run must either exit 1 with the output path untouched or exit 0 with
the expected bytes, and both kinds must occur. No cap is tuned to a binary.

## Quality tooling (release tier)

Variants are built by `scripts/build_variant.sh {asan|mutant ID CPP|mutant-asan ID CPP|seenset
CPP|portable}`, each in its own `~/muscle-build/` directory with a `BUILD_INFO.txt` naming its base
build.

- **Sanitizers.** `scripts/test_sanitize.sh` runs the hand, error, fixture and 30 tiny cases on the
  `asan` build (AddressSanitizer, UBSan, `_GLIBCXX_ASSERTIONS`; the 2026-10-30 audit's
  recipe). A report exits 86 (an assertion, 134), so the tests fail through exit status alone.
  `scripts/ubsan.supp` silences MUSCLE's known null `memcpy` argument in `fasta2.cpp`. Address-space
  caps are lifted with `CLOAK_VMEM_KB=unlimited`, which every test script honours. Valgrind is not
  installed in the WSL image and is not run.

## Measured (2026-11-02, laptop, WSL Debian 12)

| Run | Result |
|---|---|
| `test_spec_check.py` | 12 tests pass, including 12 corrupted outputs rejected |
| harness on `oracle_bin.sh` (positive control) | hand 42 of 42; `run_suite.py --count 60 --synth 50 --relations 3` with baseline and `cloak.py`: 800 comparisons, 680 against `f3559af`, 289 against `cloak.py`, 1,513 relation checks, 0 failed, 0 skipped, 141 s |
| red step, hand cases | `ddc9df1` and `f3559af` each fail 16 of 42, exactly as predicted |
| red step, random cases (same run) | `ddc9df1` fails 388 of 800 comparisons, `f3559af` 413; every failure is a case with a sequence keeping at most one residue, a `.` in the input, or a relation that inserts `.` or raises `-mincol` |
| red step for the errors commit, on `20798ce` | `test_errors.sh`: E20a, E20b, E20c and E21 exit 0 (E20c leaves 51,200 bytes), H17 passes; `test_oom.sh`: caps of 16 to 45 MiB exit 0 with no output, 64 and 90 MiB exit 134, 128 MiB and above pass |
| errors commit, `test.sh commit` | hand 37 of 42 (the 5 `.` cases fail); errors 5 of 5; oom 15 of 15 (6 exit 1, 9 exit 0; 5 s); quick 23 of 23; tiny and random fail the same 20 and 143 comparisons as `20798ce`, all `.` cases or relation R3 |
| `.` commit, `test.sh commit` | every part passes: hand 42 of 42; random 800 comparisons (680 against `f3559af`, 289 against `cloak.py`) and 1,513 relation checks, 0 failed, 0 skipped. The red step is the errors commit's failing set (5 hand cases, 20 tiny and 143 random comparisons, all `.`) |
| sanitizers, `e0eda41` | `build_variant.sh asan` 2.5 min; `test_sanitize.sh`: hand 42 of 42, errors 5 of 5, quick 23 of 23, tiny pass; no report; 43 s |
| sanitizers, positive control | mutant M09 (`n = 0` guard removed) on the asan tree fails E06 with exit 134 (`_GLIBCXX_ASSERTIONS`: `operator[]` out of range) |
| positive control, errors commit | `test.sh commit` on `oracle_bin.sh`: every part passes (errors 5 of 5), oom reported `n/a` |

The last row means that, wherever the known defects cannot act, both old binaries agree with the
definition, the byte oracle, `cloak.py` and every relation.

## Rules

- A failing test is evidence: never weaken a test or change existing test data to make it pass
  (decision 009; `scripts/check_fixtures.sh`).
- Tests assert exit status, presence or absence of the output, and output bytes or parsed structure;
  never stderr or log text.
- Errors must exit exactly 1: a crash (139), an abort (134) or a timeout is a failure.
- Every run is capped (`ulimit` or `tools/run_suite.py`'s limits); seeds are fixed by schedule and
  number; every skip is named and fails the run unless `--allow-skips`.
