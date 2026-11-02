# Test fixtures

Reference data for the test scripts in `cloak_dev/scripts/`. Never edit, regenerate or delete a
file here to make a test pass: a failing fixture is evidence. New fixtures may be added by anyone,
each with its line in `SHA256SUMS`; existing files and lines change only for a specification change
recorded in `DECISIONS.md`, with Peter's explicit approval (decision 009).
`scripts/check_fixtures.sh` checks this and runs first in every test script.

`cloak_dev/.gitattributes` marks this directory `-text`, so Windows checkout does not convert
line endings and the expected outputs stay byte-comparable.

## Cases

`MANIFEST.tsv` lists 23 cases: an ensemble, a `-mincol` value, the expected output, the source of
that output, and whether the frozen baseline (`~/muscle-baseline/muscle-f3559af`) can be run on it
routinely.

| source | meaning |
|---|---|
| `f3559af` | written by a build of upstream `f3559af` (cloak-v0.1.0) |
| `april2025` | written by Peter's April 2025 Windows build (same `cloak.cpp`), and reproduced byte for byte by `f3559af` builds in the Cowork sandbox and on Peter's laptop on 2026-10-29 |
| `independent` | `f3559af` cannot finish (memory). Written by the cloak-v0.2.0 candidate and byte-identical to two independent implementations: an early-exit variant of the `f3559af` recursion and `cloak_dev/tools/oracle_direct.py` |
| `v030` | the cloak-v0.3.0 output, where it differs from 0.1.0's (decision 007); see below |

`expected/PF00353full_n1226_s5.m2.v030.fa` (2026-11-02) is the cloak-v0.3.0 output, with all 1,226
rows, written by `tools/oracle_direct.py` before the code change and then reproduced by the C1
build. Removing its rows with at most one residue gives `expected/PF00353full_n1226_s5.m2.fa` (the
cloak-v0.1.0 and 0.2.0 output, kept) byte for byte, and `tools/spec_check.py check --order` accepts
it. Its MANIFEST row was changed with Peter's approval (decision 009). Every other case keeps every
row, so its 0.1.0 output is also its 0.3.0 output.

| baseline | meaning |
|---|---|
| `fast` | the baseline finishes in under 20 s and 400 MB; checked by `test_quick.sh --check-baseline` |
| `slow` | the baseline needs about 2 min and 0.6 to 1.3 GB; checked only by `--check-baseline-slow` |
| `never` | the baseline exhausts memory; never run it on these without a cap |

## Ensembles

All are 16-MSA ensembles in EFA format.

- `PF*.efa`: made on 2026-10-27 from Pfam seed alignments (`PFxxxxx_nN`, N sequences) or Pfam full
  alignments (`PFxxxxxfull_nN`), gaps removed, randomly subsampled, and realigned with muscle
  `f3559af`. Names ending `_s5` are 4 x `-super5 -perturb 0..3 -perm all`, the usual recipe for
  large families; the others are `-align -stratified`. Recipe and source URLs:
  `cloak_dev/tools/fixture_provenance/make_fixture_ensembles.sh`. Pfam data are CC0.
- `april2025_n49_s5.efa`: Peter's April 2025 test in `C:\Users\peter\BioInformatics\muscle_testing`
  (49 sequences; MSAs `abc.0` to `none.3`), concatenated in the order of its `cloak_test.txt`.

These ensembles were chosen to exercise the code, not to match production data: gate B (`PLAN.md`,
phase 1) is the test for that.

### Hand fixtures

`hand/MANIFEST.tsv` lists 42 cases on tiny ensembles written for cloak-v0.3.0. Their expected
outputs were written by hand from the definition (decisions 007 and 008), and
`tools/make_hand_fixtures.py --verify` confirms that `tools/spec_check.py` and
`tools/oracle_direct.py` agree with every one; the same script regenerates this directory byte
for byte. Columns: id; args (after the binary; `{out}` is the output path, other paths are relative
to `hand/`, the working directory of every run); exit; expected output (`-` for none); v010 and
v020, whether `f3559af` and `ddc9df1` are predicted to pass or fail the case; what the case
targets. `scripts/test_hand.sh` runs them. A case with exit 1 also requires that nothing is written
at the output path and no file appears in the working directory.
