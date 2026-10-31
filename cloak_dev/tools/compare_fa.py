"""Compare two CLOAK outputs. Usage: compare_fa.py a.fa b.fa
Prints: n_rows_a n_rows_b n_cols_a n_cols_b rows_identical column_multiset_equal (labels+residues per column)"""
import sys
from collections import Counter


def read(path):
    d, lab = {}, None
    for line in open(path):
        line = line.rstrip('\n')
        if line.startswith('>'):
            lab = line[1:].strip()
            d[lab] = []
        elif line:
            d[lab].append(line.strip())
    return {k: ''.join(v) for k, v in d.items()}


a, b = read(sys.argv[1]), read(sys.argv[2])


def cols(d):
    labs = sorted(d)
    L = len(next(iter(d.values()))) if d else 0
    return Counter(tuple((l, d[l][j]) for l in labs if d[l][j] != '-') for j in range(L)), L


ca, la = cols(a)
cb, lb = cols(b)
rows_ident = (list(a.items()) == list(b.items()))
print(f'rows {len(a)} {len(b)}\tcols {la} {lb}\trows_identical={rows_ident}\tcolumn_multiset_equal={ca == cb}'
      f'\tlabels_equal={sorted(a) == sorted(b)}')
