"""Digit-plane tagging: which input residue each output residue really is (tests/README.md).

A row that repeats a letter (say AA) writes the same bytes whichever residue of the run CLOAK
picks, so byte comparisons cannot see a wrong pick. For each case of tests/data/MANIFEST.tsv this
tool writes P copies of the ensemble in which residue i (1-based) of every row becomes, in every
MSA, digit d of i - 1 in base 36 (0-9, A-Z), for d = 0..P-1, where P is the fewest planes that
number the longest row. Each copy goes through the binary. Every output must have the untagged
output's labels and gap mask (relation R17), and the digits at each output position give the
index of the residue written there. Each output row's indices must then increase, and the output
columns, as sets of (row, index), must be exactly the definition's classes with at least -mincol
residues (spec_check.expected), each once.

Usage: tag_runs.py BINARY OUTDIR [--only TEXT]   (cases whose ensemble name contains TEXT)
Every run is capped at 4 GB of address space and 300 s. Exit 0 only if every case passes.
"""
import argparse
import os
import resource
import subprocess
import sys

import cloakio
import spec_check

DIG = b'0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tests', 'data')


def run(binary, inp, out, mincol):
    def limits():
        resource.setrlimit(resource.RLIMIT_AS, (4 * 2**30,) * 2)
    if os.path.exists(out):
        os.remove(out)
    cmd = ['bash', binary] if binary.endswith('.sh') else [binary]   # the byte oracle, as in run_suite
    p = subprocess.run(cmd + ['-cloak', inp, '-output', out, '-mincol', mincol], preexec_fn=limits,
                       timeout=300, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if p.returncode != 0 or not os.path.exists(out):
        return None
    return [(lab, seq) for lab, seq, _ in cloakio.records(cloakio.slurp(out))]


def tagged(labels, msas, plane):
    """EFA bytes: the ensemble with residue i of every row replaced by digit `plane` of i - 1."""
    out = bytearray()
    for k, m in enumerate(msas):
        out += b'<m%d\n' % k
        for lab, row in zip(labels, m):
            t, i = bytearray(row), 0
            for c, ch in enumerate(row):
                if ch not in cloakio.GAPS:
                    t[c] = DIG[(i // 36 ** plane) % 36]
                    i += 1
            out += b'>' + lab + b'\n' + bytes(t) + b'\n'
    return bytes(out)


def check_case(binary, outdir, efa, mincol):
    inp = os.path.join(DATA, efa)
    labels, msas, classes = spec_check.expected(inp, cloakio.parse_mincol(mincol))
    longest = max(len(cloakio.ungap(row)) for row in msas[0])
    planes = 1
    while 36 ** planes < longest:
        planes += 1
    name = f'{os.path.splitext(efa)[0]}.m{mincol}'
    base = run(binary, inp, os.path.join(outdir, name + '.fa'), mincol)
    if base is None:
        return f'untagged run failed', planes
    if [lab for lab, _ in base] != labels:
        return 'output labels differ from the sorted input labels', planes
    idx = [[0 if ch != 0x2D else None for ch in seq] for _, seq in base]
    for d in range(planes):
        efa_d = os.path.join(outdir, f'{name}.plane{d}.efa')
        with open(efa_d, 'wb') as f:
            f.write(tagged(labels, msas, d))
        out = run(binary, efa_d, os.path.join(outdir, f'{name}.plane{d}.fa'), mincol)
        if out is None:
            return f'plane {d}: run failed', planes
        for r, ((lab, seq), (_, seq0)) in enumerate(zip(out, base)):
            if lab != labels[r] or len(seq) != len(seq0) or \
                    any((a == 0x2D) != (b == 0x2D) for a, b in zip(seq, seq0)):
                return f'plane {d}: labels or gap mask differ from the untagged output (R17)', planes
            for j, ch in enumerate(seq):
                if ch != 0x2D:
                    idx[r][j] += DIG.index(ch) * 36 ** d
    for r, row in enumerate(idx):
        got = [i for i in row if i is not None]
        if got != sorted(set(got)):
            return f'row {labels[r]!r}: residue indices do not increase', planes
    width = len(base[0][1]) if base else 0
    cols = [frozenset((r, idx[r][j] + 1) for r in range(len(idx)) if idx[r][j] is not None)
            for j in range(width)]
    if len(set(cols)) != len(cols):
        return 'two output columns hold the same residues', planes
    if set(cols) != set(classes.values()):
        missing = len(set(classes.values()) - set(cols))
        extra = len(set(cols) - set(classes.values()))
        return f'output columns differ from the classes ({missing} missing, {extra} not classes)', planes
    return None, planes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('binary')
    ap.add_argument('outdir')
    ap.add_argument('--only', default='')
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    binary = os.path.abspath(a.binary)
    n = bad = 0
    with open(os.path.join(DATA, 'MANIFEST.tsv'), newline='') as f:
        rows = [line.rstrip('\r\n').split('\t') for line in f if line.strip()]
    for efa, mincol, *_ in rows[1:]:
        if a.only not in efa:
            continue
        n += 1
        why, planes = check_case(binary, a.outdir, efa, mincol)
        bad += why is not None
        print(f'{"pass" if why is None else "FAIL"}  {efa} m{mincol}  {planes} plane(s)'
              + ('' if why is None else f'  ({why})'), flush=True)
    print(f'cases {n}: failed {bad}')
    return 1 if bad or n == 0 else 0


if __name__ == '__main__':
    sys.exit(main())
