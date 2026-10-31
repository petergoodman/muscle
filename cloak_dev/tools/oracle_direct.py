"""Test oracle only: an independent re-implementation of CLOAK's output, used to cross-check the C++
on inputs the f3559af binary cannot finish. It is not a CLOAK implementation for analyses
(decision 001).

Computes, for each column of MSA 0, the groups of residues that occupy one column in every MSA
(key = column in MSAs 1..K-1), keeps groups of at least mincol residues in lexicographic key order,
and writes FASTA the way cloak.cpp does: rows sorted by label, input uppercased, only '-' is a gap,
rows with at most one residue left out, 60 characters per line.

Usage: oracle_direct.py in.efa out.fa [mincol]
"""
import sys
from collections import defaultdict


def read_efa(path):
    msas, cur, lab = [], None, None
    for line in open(path):
        line = line.rstrip('\n')
        if not line:
            continue
        if line.startswith('<'):
            cur = {}
            msas.append(cur)
        elif line.startswith('>'):
            lab = line[1:]
            if lab in cur:
                sys.exit(f'duplicate label {lab}')
            cur[lab] = []
        else:
            cur[lab].append(''.join(line.split()).upper())
    return [{k: ''.join(v) for k, v in m.items()} for m in msas]


def main():
    efa, out = sys.argv[1], sys.argv[2]
    mincol = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    msas = read_efa(efa)
    labels = sorted(msas[0])
    for m in msas[1:]:
        if sorted(m) != labels:
            sys.exit('label sets differ')
    rows = [[m[l] for l in labels] for m in msas]          # rows[k][r] = aligned string
    n, K = len(labels), len(rows)
    # col[k][r][p] = column of residue p (0-based) of sequence r in MSA k
    col = []
    for k in range(K):
        ck = []
        for r in range(n):
            ck.append([c for c, ch in enumerate(rows[k][r]) if ch != '-'])
        col.append(ck)
    for k in range(1, K):
        for r in range(n):
            if len(col[k][r]) != len(col[0][r]):
                sys.exit(f'ungapped length differs: MSA {k} row {labels[r]}')

    out_cols = []                                        # each: list of (row, residue index)
    ptr = [0] * n                                        # next residue index per row in MSA 0
    for x in range(len(rows[0][0])):
        groups = defaultdict(list)
        for r in range(n):
            if rows[0][r][x] != '-':
                p = ptr[r]
                ptr[r] += 1
                key = tuple(col[k][r][p] for k in range(1, K))
                groups[key].append((r, p))
        for key in sorted(groups):
            g = groups[key]
            if len(g) >= mincol:
                out_cols.append(g)

    ungapped = [rows[0][r].replace('-', '') for r in range(n)]
    seqs = [['-'] * len(out_cols) for _ in range(n)]
    for j, g in enumerate(out_cols):
        for r, p in g:
            seqs[r][j] = ungapped[r][p]
    with open(out, 'w') as o:
        for r in range(n):
            s = ''.join(seqs[r])
            if len(s) - s.count('-') > 1:
                o.write(f'>{labels[r]}\n')
                for i in range(0, len(s), 60):
                    o.write(s[i:i + 60] + '\n')


if __name__ == '__main__':
    main()
