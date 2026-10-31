"""Generate synthetic 16-MSA ensembles (EFA) of random protein sequences with controlled disagreement.

Each ensemble shares one base placement of every sequence's residues in a common width; in MSAs
1..15 each row is independently re-placed at random with probability p (residue order kept,
all-gap columns removed, row order shuffled per MSA as aligners do). The data have no biological
meaning: they exercise CLOAK's grouping logic on many disagreement patterns.

Single ensemble:  gen_synth.py SEED N_SEQS MAX_LEN P OUT.efa
Batch:            gen_synth.py --batch FIRST_SEED COUNT OUTDIR {mixed|mild}
  mixed: n = 3 + seed % 8,  max_len = 4 + seed % 9, p cycles 0.1, 0.3, 0.6, 1.0 (heavy disagreement)
  mild:  n = 4 + seed % 12, max_len = 4 + seed % 9, p cycles 0.03, 0.06, 0.1, 0.15
Batch files are named OUTDIR/<schedule>_<seed>.efa. Same seed and parameters, same file.
"""
import os
import random
import sys

AA = 'ACDEFGHIKLMNPQRSTVWY'


def make(seed, n, maxlen, p, out):
    rng = random.Random(seed)
    seqs = [''.join(rng.choice(AA) for _ in range(rng.randint(2, maxlen))) for _ in range(n)]
    width = maxlen + rng.randint(0, 4)

    def place(length):
        return sorted(rng.sample(range(width), length))

    base = [place(len(s)) for s in seqs]
    with open(out, 'w', newline='\n') as o:
        for k in range(16):
            pl = [place(len(s)) if (k > 0 and rng.random() < p) else base[r] for r, s in enumerate(seqs)]
            grid = [['-'] * width for _ in seqs]
            for r, s in enumerate(seqs):
                for c, ch in zip(pl[r], s):
                    grid[r][c] = ch
            keep = [c for c in range(width) if any(grid[r][c] != '-' for r in range(n))]
            o.write(f'<m{k}\n')
            order = list(range(n))
            rng.shuffle(order)
            for r in order:
                o.write(f'>s{r:03d}\n' + ''.join(grid[r][c] for c in keep) + '\n')


def params(seed, schedule):
    if schedule == 'mixed':
        return 3 + seed % 8, 4 + seed % 9, [0.1, 0.3, 0.6, 1.0][seed % 4]
    if schedule == 'mild':
        return 4 + seed % 12, 4 + seed % 9, [0.03, 0.06, 0.1, 0.15][seed % 4]
    sys.exit('schedule must be mixed or mild')


if __name__ == '__main__':
    a = sys.argv[1:]
    if a and a[0] == '--batch':
        first, count, outdir, schedule = int(a[1]), int(a[2]), a[3], a[4]
        os.makedirs(outdir, exist_ok=True)
        for seed in range(first, first + count):
            n, maxlen, p = params(seed, schedule)
            make(seed, n, maxlen, p, os.path.join(outdir, f'{schedule}_{seed}.efa'))
    elif len(a) == 5:
        make(int(a[0]), int(a[1]), int(a[2]), float(a[3]), a[4])
    else:
        sys.exit(__doc__)
