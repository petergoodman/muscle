"""Checks CLOAK output against the published definition; written from the paper, not from cloak.cpp.

Definition (Mol Biol Evol 2026, doi 10.1093/molbev/msag182, Methods and Fig. 6; decisions 007 and
008): a residue is (row, i), i its 1-based index in the row's sequence with gaps ('-', '.') removed.
Two residues are in one output column iff they share a column in every input alignment; the output
columns are the classes of that relation with at least mincol residues; every input sequence is
written, in label order, with '-' for gaps and 60 characters per line.

Classes are computed twice: by grouping residues on their tuple of columns, and pair by pair from
each alignment's same-column relation, intersected over the alignments and checked for
transitivity (the second only up to PAIR_LIMIT pairs per alignment).

An output row shows which residues it keeps only through its letters. A correct output keeps, in
row r, exactly the residues R_r of the expected classes, in increasing order (a row has at most one
residue per column of alignment 0), so once the row's letters match, the k-th residue written is
the k-th index of R_r, and the output columns must then be exactly the expected classes. Column
order is reported; --order also requires CLOAK's order (alignment-0 column, then the columns in
alignments 1..K-1, lexicographically).

Usage:
  spec_check.py check IN OUT [--mincol M] [--order]   exit 0 if OUT is a correct output for IN
  spec_check.py expect IN [--mincol M]                 expected exit status, then columns or reason
  spec_check.py stats IN OUT [--mincol M]              rows in, rows out, rows missing, columns,
                                                       columns below mincol, all-gap columns
IN is an EFA or a path list (paths relative to the working directory). Exit 2: IN is outside the
model (cloakio.Undefined) or the arguments are wrong.
"""
import argparse
import sys
from collections import Counter, defaultdict

from cloakio import GAPS, WIDTH, LoadError, Undefined, load, parse_mincol, records, slurp, ungap

PAIR_LIMIT = 200_000


def residue_columns(msas):
    """col[(r, i)] = the residue's column in each alignment, in alignment order."""
    col = defaultdict(list)
    for m in msas:
        for r, row in enumerate(m):
            i = 0
            for c, ch in enumerate(row):
                if ch not in GAPS:
                    i += 1
                    col[(r, i)].append(c)
    return {a: tuple(v) for a, v in col.items()}


def classes_by_tuple(col):
    groups = defaultdict(set)
    for a, t in col.items():
        groups[t].add(a)
    return {t: frozenset(g) for t, g in groups.items()}


def classes_by_pairs(col, K):
    """R1 read pair by pair; returns None if an alignment has more than PAIR_LIMIT same-column pairs."""
    pairs = None
    for k in range(K):
        bycol = defaultdict(list)
        for a, t in col.items():
            bycol[t[k]].append(a)
        if sum(len(v) * (len(v) - 1) // 2 for v in bycol.values()) > PAIR_LIMIT:
            return None
        pk = {(a, b) for v in bycol.values() for a in v for b in v if a < b}
        pairs = pk if pairs is None else pairs & pk
    partner = defaultdict(set)
    for a, b in pairs or ():
        partner[a].add(b)
        partner[b].add(a)
    classes = set()
    for a in col:
        cls = frozenset(partner[a] | {a})
        if any(frozenset(partner[b] | {b}) != cls for b in cls):
            raise AssertionError('same-column relation is not transitive')
        classes.add(cls)
    return classes


def expected(inp, mincol):
    """(labels, msas, {column tuple: class}) for the classes with at least mincol residues."""
    labels, msas = load(inp)
    col = residue_columns(msas)
    by_tuple = classes_by_tuple(col)
    by_pairs = classes_by_pairs(col, len(msas))
    if by_pairs is not None and by_pairs != set(by_tuple.values()):
        raise AssertionError('tuple and pair formulations disagree')
    return labels, msas, {t: c for t, c in by_tuple.items() if len(c) >= mincol}


def format_problems(data):
    """S7: '>label' lines, then sequence lines of 60 (the last of a record 1 to 60), LF, final LF."""
    if b'\r' in data:
        return [], ['output contains CR']
    if not data.endswith(b'\n'):
        return [], ['output does not end with LF']
    recs, problems = [], []
    for s in data[:-1].split(b'\n'):
        if s.startswith(b'>'):
            recs.append((s[1:], []))
        elif not recs:
            return [], ['output does not start with >']
        else:
            recs[-1][1].append(s)
    for lab, lines in recs:
        if not lines or any(not 0 < len(s) <= WIDTH for s in lines) \
                or any(len(s) != WIDTH for s in lines[:-1]):
            problems.append(f'record {lab!r} is not wrapped at {WIDTH}')
    return [(lab, b''.join(lines)) for lab, lines in recs], problems


def check(inp, out, mincol, order=False):
    """Returns (problems, note)."""
    try:
        labels, msas, want = expected(inp, mincol)
    except LoadError as e:
        return [f'CLOAK must exit 1 without output: invalid input ({e})'], ''
    if not want:
        return [f'CLOAK must exit 1 without output: no class of {mincol} or more residues'], ''
    recs, problems = format_problems(slurp(out))
    if problems or not recs:
        return problems or ['empty output'], ''
    if [lab for lab, _ in recs] != labels:
        missing = sorted(set(labels) - {lab for lab, _ in recs})
        return [f'labels differ from the input labels in sorted order; missing {missing[:5]}'], ''
    rows = [row for _, row in recs]
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        return ['rows differ in length'], ''
    kept = defaultdict(list)
    for cls in want.values():
        for r, i in cls:
            kept[r].append(i)
    columns = [set() for _ in range(width)]
    for r, row in enumerate(rows):
        idx = sorted(kept[r])
        seq = ungap(msas[0][r])
        if row.replace(b'-', b'') != bytes(seq[i - 1] for i in idx):
            problems.append(f'row {labels[r]!r}: letters differ from the residues it must keep')
            continue
        k = 0
        for j, ch in enumerate(row):
            if ch != 0x2D:
                columns[j].add((r, idx[k]))
                k += 1
    if problems:
        return problems, ''
    got = Counter(frozenset(c) for c in columns)
    if got != Counter(want.values()):
        empty = got.get(frozenset(), 0)
        extra = sum((got - Counter(want.values())).values())
        missing = sum((Counter(want.values()) - got).values())
        return [f'columns differ from the classes: {extra} unexpected ({empty} all-gap), '
                f'{missing} missing'], ''
    canonical = [want[t] for t in sorted(want)] == [frozenset(c) for c in columns]
    if order and not canonical:
        return ['column order differs from CLOAK order'], ''
    return [], 'order: ' + ('CLOAK order' if canonical else 'other')


def stats(inp, out, mincol):
    labels, _ = load(inp)
    recs = records(slurp(out))
    rows = [seq for _, seq, _ in recs]
    width = len(rows[0]) if rows else 0
    counts = [sum(row[j] != 0x2D for row in rows) for j in range(width)]
    missing = len(set(labels) - {lab for lab, _, _ in recs})
    return [len(labels), len(rows), missing, width, sum(c < mincol for c in counts),
            counts.count(0)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('mode', choices=['check', 'expect', 'stats'])
    ap.add_argument('inp')
    ap.add_argument('out', nargs='?')
    ap.add_argument('--mincol')
    ap.add_argument('--order', action='store_true')
    a = ap.parse_args()
    if (a.out is None) != (a.mode == 'expect'):
        ap.error('check and stats take IN and OUT; expect takes IN only')
    try:
        if a.mode == 'expect':
            try:
                mincol = parse_mincol(a.mincol)
                _, _, want = expected(a.inp, mincol)
            except LoadError as e:
                print(f'1\t{e}')
                return 0
            print(f'0\t{len(want)}' if want else f'1\tno class of {mincol} or more residues')
            return 0
        mincol = parse_mincol(a.mincol) if a.mode == 'check' else int(a.mincol or 2)
        if a.mode == 'stats':
            print('\t'.join(map(str, stats(a.inp, a.out, mincol))))
            return 0
        problems, note = check(a.inp, a.out, mincol, a.order)
    except Undefined as e:
        print(f'undefined input: {e}', file=sys.stderr)
        return 2
    except LoadError as e:
        print(f'CLOAK must exit 1 without output: {e}', file=sys.stderr)
        return 1
    for p in problems:
        print(p, file=sys.stderr)
    if note:
        print(note)
    return 1 if problems else 0


if __name__ == '__main__':
    sys.exit(main())
