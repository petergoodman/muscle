# Test fixtures

Reference data for `cloak_dev/scripts/test_quick.sh`. Never edit, regenerate or delete a file
here to make a test pass: a failing fixture is evidence. Changes to this directory are made by
Peter only, with the reason in the commit message.

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

