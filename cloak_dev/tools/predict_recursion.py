"""Exactly count the work cloak.cpp (f3559af) would do on an EFA, without running it.

For each MSA-0 column whose residues are split by some other MSA, the recursion add_to_output_or_discard(P)
visits, for every MSA i>=1 and every block B of P's split by MSA i (|B| >= 2, B != P): emits B if B lies in
one column of every MSA, else recurses on B. The tree is re-walked along every path, but its value depends only
on the residue set P, so calls(P) and emits(P) can be memoised on P. Reports, per ensemble:
total recursive calls, total emitted entries (with duplicates), distinct emitted entries, and the estimated
bytes the C++ 'output' vector would hold (entries x (8*n + 40)).

Usage: predict_recursion.py in.efa [per_column_tsv]
"""
import sys
from functools import lru_cache
sys.setrecursionlimit(100000)


def read_efa(path):
    msas, cur, lab = [], None, None
    for line in open(path):
        line = line.rstrip('\n')
        if not line:
            continue
        if line.startswith('<'):
            cur = {}; msas.append(cur)
        elif line.startswith('>'):
            lab = line[1:]; cur[lab] = []
        else:
            cur[lab].append(''.join(line.split()).upper())
    return [{k: ''.join(v) for k, v in m.items()} for m in msas]


msas = read_efa(sys.argv[1])
labels = sorted(msas[0])
rows = [[m[l] for l in labels] for m in msas]
n, K = len(labels), len(rows)
col = [[[c for c, ch in enumerate(rows[k][r]) if ch != '-'] for r in range(n)] for k in range(K)]

total_calls = total_emits = total_unique = 0
worst = (0, None)
per_col = []
ptr = [0] * n
for x in range(len(rows[0][0])):
    res = []
    for r in range(n):
        if rows[0][r][x] != '-':
            res.append((r, ptr[r])); ptr[r] += 1
    if len(res) == 0:
        continue
    keys = {rp: tuple(col[k][rp[0]][rp[1]] for k in range(1, K)) for rp in res}
    splitters = [k for k in range(1, K) if len({keys[rp][k - 1] for rp in res}) > 1]
    if not splitters:
        continue

    def consistent(P):
        ks = {keys[rp] for rp in P}
        return len(ks) == 1

    @lru_cache(maxsize=None)
    def walk(P):  # P: frozenset of residues; returns (calls, emits)
        calls, emits = 1, 0
        for i in range(1, K):
            blocks = {}
            for rp in P:
                blocks.setdefault(keys[rp][i - 1], []).append(rp)
            for B in blocks.values():
                if len(B) < 2 or len(B) == len(P):
                    continue
                Bf = frozenset(B)
                if consistent(Bf):
                    emits += 1
                else:
                    c, e = walk(Bf)
                    calls += c; emits += e
        return calls, emits

    c, e = walk(frozenset(res))
    classes = {}
    for rp in res:
        classes.setdefault(keys[rp], 0); classes[keys[rp]] += 1
    u = sum(1 for v in classes.values() if v >= 2)
    total_calls += c; total_emits += e; total_unique += u
    per_col.append((x, len(res), len(splitters), u, c, e, walk.cache_info().currsize))
    if c > worst[0]:
        worst = (c, x)

entry_bytes = 8 * n + 40
print(f'n={n}\tK={K}\tcols0={len(rows[0][0])}\tsplit_cols={len(per_col)}\tcalls={total_calls:.4g}\t'
      f'emits={total_emits:.4g}\tunique={total_unique}\tdup_factor={total_emits / max(total_unique, 1):.3g}\t'
      f'est_output_GB={total_emits * entry_bytes / 1e9:.4g}\tworst_col_calls={worst[0]:.4g}')
if len(sys.argv) > 2:
    with open(sys.argv[2], 'w') as o:
        o.write('col\tresidues\tn_splitting_msas\tclasses_ge2\tcalls\temits\tdistinct_nodes\n')
        for t in per_col:
            o.write('\t'.join(map(str, t)) + '\n')
