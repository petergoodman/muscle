"""Random ensembles for the test driver (run_suite.py), with schedules aimed at the cases that
matter. Test tooling only (decision 001).

Each ensemble starts from one random placement of random sequences (MSA 0) and derives MSAs
1..K-1 by operations that keep every row's residue order, so each MSA is a valid alignment of the
same sequences (the operations of the 2026-10-30 audit's generators, re-expressed here):
  fission  split a column's residues into 2 to 4 adjacent new columns (multi-class MSA-0 columns)
  fusion   merge two adjacent columns that no row occupies twice
  slide    move one residue within its gap run
  rogue    place a whole row afresh
Schedules (parameters drawn from the seed):
  fission  n 4-40, length 3-12, K 2-16; mostly fission and slide
  rogue    n 3-12, length 0-4, K 2-16; rows placed afresh (rows keeping 0 or 1 residues)
  dots     a fission or rogue ensemble with '-' made '.' in MSA 0, in other MSAs or everywhere
           (leading, trailing, interior or random gaps)
  kn       K in {1, 2, 3, 4, 8, 16, 30}, n in {1, 2, 3, 5, 10, 40, 120, 200}, length 2-10
  tiny     K 1-5, n 1-7, length 0-4 over {A, C, D, E}: for the pairwise and cloak.py oracles
Letters: distinct letters and digits per row (each output residue identified by its letter),
protein, or {A, C}; some ensembles partly lower case. Labels: plain, or a hostile mix (spaces,
case, digits that sort differently as text). Forms: EFA; CRLF EFA wrapped at 7 with blanks; or a
path list of FASTA files (letters and digits only, as MUSCLE's FASTA reader keeps). All-gap columns
are kept in some MSAs. Same schedule and seed, same files.

Usage: gen_cases.py OUTDIR SCHEDULE FIRST_SEED COUNT
Writes OUTDIR/<schedule>_<seed>.efa, or OUTDIR/<schedule>_<seed>.txt listing
<schedule>_<seed>/m<k>.fa (paths relative to OUTDIR, the working directory for the run).
"""
import os
import random
import sys

ALNUM = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'
PROTEIN = 'ACDEFGHIKLMNPQRSTVWY'
HOSTILE = ['a9', 'a10', 'B', 'b', '_x', 's 1', 'S1', 'z', 'Z0', 'a', '1', 'x y']
SCHEDULES = ('fission', 'rogue', 'dots', 'kn', 'tiny')


def place(rng, length, width):
    return sorted(rng.sample(range(width), length))


def fission(rng, pos, width):
    used = sorted({c for p in pos for c in p})
    if not used:
        return width
    c = rng.choice(used)
    rows = [r for r, p in enumerate(pos) if c in p]
    b = min(len(rows), rng.randint(2, 4))
    if b < 2:
        return width
    block = {r: rng.randrange(b) for r in rows}
    for r, p in enumerate(pos):
        pos[r] = [x + b - 1 if x > c else (c + block[r] if x == c else x) for x in p]
    return width + b - 1


def fusion(rng, pos, width):
    if width < 2:
        return width
    c = rng.randrange(width - 1)
    if any(c in p and c + 1 in p for p in pos):
        return width
    for r, p in enumerate(pos):
        pos[r] = [x - 1 if x > c else x for x in p]
    return width - 1


def slide(rng, pos, width):
    cand = [(r, i) for r, p in enumerate(pos) for i in range(len(p))]
    if not cand:
        return width
    r, i = rng.choice(cand)
    p = pos[r]
    lo = p[i - 1] + 1 if i else 0
    hi = p[i + 1] - 1 if i + 1 < len(p) else width - 1
    p[i] = rng.randint(lo, hi)
    return width


def rogue(rng, pos, width):
    r = rng.randrange(len(pos))
    pos[r] = place(rng, len(pos[r]), width)
    return width


def derive(rng, base, width, ops, weights):
    pos = [list(p) for p in base]
    for _ in range(ops):
        width = rng.choices((fission, fusion, slide, rogue), weights)[0](rng, pos, width)
    return pos, width


def grid(rng, seqs, pos, width, drop_gap_cols):
    rows = [['-'] * width for _ in seqs]
    for r, (s, p) in enumerate(zip(seqs, pos)):
        for ch, c in zip(s, p):
            rows[r][c] = ch
    keep = [c for c in range(width) if not drop_gap_cols or any(row[c] != '-' for row in rows)]
    keep = keep or [0]
    return [''.join(row[c] for c in keep) for row in rows]


def dot(rng, row, where):
    res = [i for i, ch in enumerate(row) if ch != '-']
    first, last = (res[0], res[-1]) if res else (len(row), -1)
    out = list(row)
    for i, ch in enumerate(row):
        if ch == '-' and (where == 'all' or (where == 'leading' and i < first)
                          or (where == 'trailing' and i > last)
                          or (where == 'interior' and first < i < last)
                          or (where == 'random' and rng.random() < 0.5)):
            out[i] = '.'
    return ''.join(out)


def make(schedule, seed):
    rng = random.Random(f'{schedule}:{seed}')
    base_schedule = rng.choice(('fission', 'rogue')) if schedule == 'dots' else schedule
    if base_schedule == 'fission':
        K, n, lo, hi, weights = rng.randint(2, 16), rng.randint(4, 40), 3, 12, (6, 1, 3, 0.3)
    elif base_schedule == 'rogue':
        K, n, lo, hi, weights = rng.randint(2, 16), rng.randint(3, 12), 0, 4, (1, 1, 1, 6)
    elif base_schedule == 'kn':
        K, n, lo, hi, weights = (rng.choice((1, 2, 3, 4, 8, 16, 30)),
                                 rng.choice((1, 2, 3, 5, 10, 40, 120, 200)), 2, 10, (4, 1, 3, 1))
    else:
        K, n, lo, hi, weights = rng.randint(1, 5), rng.randint(1, 7), 0, 4, (3, 1, 2, 2)
    lengths = [rng.randint(lo, hi) for _ in range(n)]
    kind = 'tiny' if schedule == 'tiny' else rng.choices(('distinct', 'protein', 'lowent'), (5, 3, 2))[0]
    seqs = []
    for L in lengths:
        if kind == 'tiny':
            seqs.append(''.join(rng.choice('ACDE') for _ in range(L)))
        elif kind == 'distinct' and L <= len(ALNUM):
            seqs.append(''.join(rng.sample(ALNUM, L)))
        else:
            seqs.append(''.join(rng.choice(PROTEIN if kind != 'lowent' else 'AC') for _ in range(L)))
    width = max(lengths + [1]) + rng.randint(0, 4)
    base = [place(rng, L, width) for L in lengths]
    ops = max(1, n // 4) if base_schedule != 'rogue' else n
    msas = [grid(rng, seqs, base, width, rng.random() < 0.7)]
    for _ in range(1, K):
        pos, w = derive(rng, base, width, rng.randint(1, ops), weights)
        msas.append(grid(rng, seqs, pos, w, rng.random() < 0.7))
    if schedule == 'dots':
        which = rng.choice(('msa0', 'others', 'all'))
        where = rng.choice(('all', 'leading', 'trailing', 'interior', 'random'))
        for k, m in enumerate(msas):
            if which == 'all' or (which == 'msa0') == (k == 0):
                msas[k] = [dot(rng, row, where) for row in m]
    if rng.random() < 0.2:
        msas = [[row.lower() if rng.random() < 0.5 else row for row in m] for m in msas]
    if rng.random() < 0.3 and n <= len(HOSTILE):
        labels = rng.sample(HOSTILE, n)
    else:
        labels = [f's{r:03d}' for r in range(n)]
    form = rng.choices(('efa', 'crlf', 'pathlist'), (6, 2, 2))[0]
    return labels, msas, form


def write(outdir, schedule, seed):
    labels, msas, form = make(schedule, seed)
    name = f'{schedule}_{seed}'
    rng = random.Random(f'order:{schedule}:{seed}')
    orders = [rng.sample(range(len(labels)), len(labels)) for _ in msas]
    if form == 'pathlist':
        os.makedirs(os.path.join(outdir, name), exist_ok=True)
        with open(os.path.join(outdir, name + '.txt'), 'w', newline='\n') as f:
            for k, (m, order) in enumerate(zip(msas, orders)):
                f.write(f'{name}/m{k}.fa\n')
                with open(os.path.join(outdir, name, f'm{k}.fa'), 'w', newline='\n') as g:
                    g.write(''.join(f'>{labels[r]}\n{m[r]}\n' for r in order))
        return name + '.txt'
    nl, wrap = ('\r\n', 7) if form == 'crlf' else ('\n', 0)
    lines = []
    for k, (m, order) in enumerate(zip(msas, orders)):
        lines.append(f'<m{k}')
        for r in order:
            lines.append('>' + labels[r])
            row = m[r]
            lines += [' ' + row[i:i + wrap] + '\t' for i in range(0, len(row), wrap)] if wrap else [row]
    with open(os.path.join(outdir, name + '.efa'), 'w', newline='') as f:
        f.write(nl.join(lines) + nl)
    return name + '.efa'


def main():
    if len(sys.argv) != 5 or sys.argv[2] not in SCHEDULES:
        sys.exit(__doc__)
    outdir, schedule, first, count = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    os.makedirs(outdir, exist_ok=True)
    for seed in range(first, first + count):
        print(write(outdir, schedule, seed))


if __name__ == '__main__':
    main()
