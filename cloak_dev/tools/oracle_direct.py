"""Test oracle only: an independent re-implementation of CLOAK's output (cloak-v0.3.0), used to
cross-check the C++ byte for byte. It is not a CLOAK implementation for analyses (decision 001).

Computes, for each column of MSA 0, the groups of residues that occupy one column in every MSA
(key = column in MSAs 1..K-1), keeps groups of at least mincol residues in lexicographic key order,
and writes FASTA as cloak-v0.3.0 does: every row (decision 007), rows sorted by label, input
uppercased, '-' and '.' as gaps (decision 008), 60 characters per line; OUT '-' is standard
output. Input that CLOAK rejects, a result with no column, and a failed write exit 1 and leave no
file at OUT, as CLOAK does.

Usage: oracle_direct.py IN OUT [MINCOL]     (IN: EFA or path list; MINCOL as given to -mincol)
       oracle_direct.py -cloak IN -output OUT [-mincol M] [-threads T] [-log L]
The second form is MUSCLE's command line (-threads and -log ignored), so that
scripts/oracle_bin.sh can stand in for the binary when testing the test harness itself.
"""
import os
import sys

from cloakio import GAPS, LoadError, load, parse_mincol, render, ungap


def cloak(labels, msas, mincol):
    K, n = len(msas), len(labels)
    m0 = msas[0]
    # pos[k][r][p] = column of residue p (0-based) of row r in MSA k
    pos = [[[c for c, ch in enumerate(row) if ch not in GAPS] for row in m] for m in msas]
    out_cols = []                                        # each: list of (row, residue index)
    ptr = [0] * n                                        # next residue index per row in MSA 0
    for x in range(len(m0[0]) if n else 0):
        groups = {}
        for r in range(n):
            if m0[r][x] not in GAPS:
                p = ptr[r]
                ptr[r] += 1
                groups.setdefault(tuple(pos[k][r][p] for k in range(1, K)), []).append((r, p))
        out_cols += [groups[key] for key in sorted(groups) if len(groups[key]) >= mincol]
    if not out_cols:
        raise LoadError(f'no column with {mincol} or more residues')
    ungapped = [ungap(row) for row in m0]
    seqs = [bytearray(b'-' * len(out_cols)) for _ in range(n)]
    for j, g in enumerate(out_cols):
        for r, p in g:
            seqs[r][j] = ungapped[r][p]
    return render(labels, [bytes(s) for s in seqs])


def main():
    a = sys.argv[1:]
    if a[:1] == ['-cloak']:
        opts = dict(zip(a[::2], a[1::2]))
        if len(a) % 2 or set(opts) - {'-cloak', '-output', '-mincol', '-threads', '-log'}:
            sys.exit('oracle_direct: bad arguments')
        inp, out, mc = opts['-cloak'], opts.get('-output', 'cloak_result.fasta'), opts.get('-mincol')
    elif len(a) in (2, 3):
        inp, out, mc = a[0], a[1], a[2] if len(a) == 3 else None
    else:
        sys.exit(__doc__)
    try:
        data = cloak(*load(inp), parse_mincol(mc))
    except LoadError as e:
        sys.exit(f'oracle_direct: {e}')
    except MemoryError:
        sys.exit('oracle_direct: out of memory')
    try:
        if out == '-':  # standard output, written unbuffered so no flush at exit can fail later
            view = memoryview(data)
            while view:
                view = view[os.write(1, view):]
        else:
            with open(out, 'wb') as o:
                o.write(data)
    except OSError as e:
        if out != '-' and os.path.isfile(out):  # a failed write leaves no output, as CLOAK does
            os.remove(out)
        sys.exit(f'oracle_direct: cannot write {out}: {e.strerror}')


if __name__ == '__main__':
    main()
