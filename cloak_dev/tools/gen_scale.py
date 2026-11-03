"""Synthetic K-MSA ensembles (EFA) for the performance and memory tests of muscle -cloak (scripts/test_perf.sh).

MSA 0 places each row's residues at random columns of a common width (gap fraction G). In MSAs
1..K-1 each row is, with probability P, re-placed inside one random window of W columns (residue
order and count kept), which mimics the local disagreement of real aligners. All-gap columns are
removed and row order is shuffled per MSA.

Usage: gen_scale.py SEED N_SEQS WIDTH K G P W OUT.efa
Prints one line: seed n width k g p w msa0_cols total_cells
"""
import random
import sys

AA = 'ACDEFGHIKLMNPQRSTVWY'


def main():
    seed, n, width, K = (int(x) for x in sys.argv[1:5])
    g, p = float(sys.argv[5]), float(sys.argv[6])
    w = int(sys.argv[7])
    out = sys.argv[8]
    rng = random.Random(seed)
    base = []
    for _ in range(n):
        base.append(sorted(c for c in range(width) if rng.random() >= g))
    seqs = [''.join(rng.choice(AA) for _ in b) for b in base]
    cells = 0
    msa0_cols = 0
    with open(out, 'w', newline='\n') as o:
        for k in range(K):
            grid = [bytearray(b'-' * width) for _ in range(n)]
            for r in range(n):
                pos = base[r]
                if k > 0 and rng.random() < p and pos:
                    lo = rng.randrange(0, max(1, width - w + 1))
                    hi = min(width, lo + w)
                    inside = [i for i, c in enumerate(pos) if lo <= c < hi]
                    if inside:
                        new = sorted(rng.sample(range(lo, hi), len(inside)))
                        pos = list(pos)
                        for i, c in zip(inside, new):
                            pos[i] = c
                for c, ch in zip(pos, seqs[r]):
                    grid[r][c] = ord(ch)
            keep = [c for c in range(width) if any(grid[r][c] != 45 for r in range(n))]
            if k == 0:
                msa0_cols = len(keep)
            cells += n * len(keep)
            o.write(f'<m{k}\n')
            order = list(range(n))
            rng.shuffle(order)
            for r in order:
                row = grid[r]
                o.write(f'>s{r:05d}\n' + ''.join(chr(row[c]) for c in keep) + '\n')
    print(seed, n, width, K, g, p, w, msa0_cols, cells)


if __name__ == '__main__':
    if len(sys.argv) != 9:
        sys.exit(__doc__)
    main()
