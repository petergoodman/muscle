"""Test oracle O3: the reference cloak.py's own functions, run on an ensemble without its file
handling. Test tooling only (decision 001).

cloak.py (phylowheeler/CLOAK d846718, kept privately in cloak_dev/private/reference/) is the
original implementation of the published algorithm. Its main() walks a directory and writes into
the working directory; this harness executes the file's text without its final `main()` call and
calls what main() calls, in main()'s order: convert_to_nums, create_output, main()'s de-duplication
loop, produce_output_in_letters, remove_singletons, transpose_output.

Input preparation: rows sorted by label (as main() does); letters uppercased (CLOAK's loader does,
cloak.py keeps case); '.' made '-' (decision 008; cloak.py numbers '.'); columns that are gaps in
every row removed from each alignment (cloak.py raises IndexError on an all-gap column of its
reference alignment; removing all-gap columns changes no class). An ensemble with no residue at
all counts as having no surviving column (cloak.py raises IndexError on it). -mincol is 2:
cloak.py has no other value.

cloak.py orders the parts of a split column by first appearance, and takes as reference whichever
file os.walk lists first, so its output is compared with CLOAK's as sets of columns and rows, not
as bytes. Its recursion is exponential in the disagreement between alignments: small inputs only.

Usage: ref_cloakpy.py CLOAK_PY IN OUT
Exit 0 with OUT written; 1 if no column survives (CLOAK exits 1 there too); 3 if cloak.py raised.
"""
import re
import sys
import traceback

from cloakio import LoadError, load, render


def load_reference(path):
    with open(path) as f:
        text, n = re.subn(r'^main\(\)\s*$', '', f.read(), flags=re.M)
    if n != 1:
        sys.exit(f'ref_cloakpy: expected exactly one main() call in {path}')
    ns = {'__name__': 'cloak_reference'}
    argv, sys.argv = sys.argv, ['cloak.py', 'unused']   # its top level reads sys.argv[1]
    try:
        exec(compile(text, path, 'exec'), ns)
    finally:
        sys.argv = argv
    return ns


def run(ns, msas):
    """msas[k][r]: rows in label order. Returns the output rows, or None if no column survives."""
    three_d = []
    for m in msas:
        rows = [row.replace(b'.', b'-') for row in m]
        keep = [c for c in range(len(rows[0])) if any(r[c] != 0x2D for r in rows)]
        three_d.append([[chr(r[c]) for c in keep] for r in rows])
    if not any(three_d[0]):
        return None
    out = ns['create_output'](ns['convert_to_nums'](three_d))
    dedup = []
    for a in out:
        if a not in dedup:
            dedup.append(a)
    letters = ns['produce_output_in_letters'](dedup, ns['transpose_output'](three_d[0], []))
    letters = ns['remove_singletons'](letters)
    if not letters:
        return None
    return [''.join(r).encode('latin-1') for r in ns['transpose_output'](letters, [])]


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    ns = load_reference(sys.argv[1])
    try:
        labels, msas = load(sys.argv[2])
    except LoadError as e:
        sys.exit(f'ref_cloakpy: invalid input: {e}')
    try:
        rows = run(ns, msas) if labels else None
    except Exception:
        traceback.print_exc()
        return 3
    if rows is None:
        print('ref_cloakpy: no column survives', file=sys.stderr)
        return 1
    with open(sys.argv[3], 'wb') as f:
        f.write(render(labels, rows))
    return 0


if __name__ == '__main__':
    sys.exit(main())
