"""Write a large synthetic ensemble (EFA) for the performance and memory tests (scripts/test_perf.sh).
Usage:
  gen_large.py perturb OUT.efa N_ROWS N_COLS N_MSAS MOVE_FRAC SEED
      MSA 0: N_ROWS rows of N_COLS columns, each cell a residue with probability 0.8.
      MSA k >= 1: MSA 0 with about MOVE_FRAC of each row's residues moved one column into an
      adjacent gap (order preserved), plus 0-3 all-gap columns inserted; rows shuffled.
  gen_large.py pairs OUT.efa N_ROWS N_COLS
      Two MSAs. MSA 0: every row is N_COLS residues, no gaps. MSA 1: row r is shifted right by
      r // 2 columns, so every column of MSA 0 splits into N_ROWS/2 groups of 2 residues
      (the number of retained columns is N_COLS * N_ROWS / 2)."""
import random, sys

AA = 'ACDEFGHIKLMNPQRSTVWY'
mode, path = sys.argv[1], sys.argv[2]


def write(path, msas):
    with open(path, 'w', newline='\n') as f:
        for name, rows in msas:
            f.write('<' + name + '\n')
            for lab, seq in rows:
                f.write('>' + lab + '\n' + seq + '\n')


if mode == 'perturb':
    n, L, K, frac, seed = (int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]),
                           float(sys.argv[6]), int(sys.argv[7]))
    rng = random.Random(seed)
    labels = ['s%06d' % i for i in range(n)]
    base = [[rng.choice(AA) if rng.random() < 0.8 else '-' for _ in range(L)] for _ in range(n)]
    msas = [('m0', [(labels[r], ''.join(base[r])) for r in range(n)])]
    for k in range(1, K):
        rows = []
        inserts = sorted(rng.randrange(L + 1) for _ in range(rng.randrange(4)))
        for r in range(n):
            row = base[r][:]
            nres = sum(c != '-' for c in row)
            for _ in range(int(round(frac * nres))):
                i = rng.randrange(L)
                if row[i] == '-':
                    continue
                d = rng.choice((-1, 1))
                j = i + d
                if 0 <= j < L and row[j] == '-':
                    row[i], row[j] = row[j], row[i]
            for pos in reversed(inserts):
                row.insert(pos, '-')
            rows.append((labels[r], ''.join(row)))
        rng.shuffle(rows)
        msas.append(('m%d' % k, rows))
    write(path, msas)
elif mode == 'pairs':
    n, L = int(sys.argv[3]), int(sys.argv[4])
    rng = random.Random(1)
    labels = ['s%06d' % i for i in range(n)]
    seqs = [''.join(rng.choice(AA) for _ in range(L)) for _ in range(n)]
    half = (n + 1) // 2
    m0 = [(labels[r], seqs[r]) for r in range(n)]
    m1 = [(labels[r], '-' * (r // 2) + seqs[r] + '-' * (half - 1 - r // 2)) for r in range(n)]
    write(path, [('m0', m0), ('m1', m1)])
else:
    sys.exit('unknown mode')
print('wrote', path)
