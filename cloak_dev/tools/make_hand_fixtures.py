"""Writes the hand fixtures for cloak-v0.3.0 (inputs, hand-written expected outputs, MANIFEST.tsv)
into a staging directory, for Peter to review and install as cloak_dev/tests/data/hand/.

The expected outputs below are written by hand from the definition (decisions 007, 008), not
computed by any tool. --verify then requires the definition checker (spec_check.py) and the byte
oracle (oracle_direct.py) to agree with every one of them.

MANIFEST.tsv columns: id; args (after the binary; {out} is the output path, other paths are
relative to the hand directory, which is the working directory of every run); exit; expected
output ('-' if none); v010 and v020: whether f3559af and ddc9df1 are predicted to pass or fail the
case (the recorded red step); what the case targets.

Usage: make_hand_fixtures.py OUTDIR                     write the fixtures (OUTDIR must not exist)
       make_hand_fixtures.py --verify DIR [CLOAK_PY]    check them against spec_check.py (O1) and
                                                        oracle_direct.py (O2), and with CLOAK_PY
                                                        against cloak.py's functions (O3, -mincol 2)
"""
import os
import sys
import tempfile


def efa(*msas, crlf=False, wrap=0, noise=False):
    """msas: lists of (label, row). wrap splits rows into lines of that length; noise adds blanks."""
    nl = '\r\n' if crlf else '\n'
    out = []
    for k, m in enumerate(msas):
        out.append(f'<m{k}')
        for lab, row in m:
            out.append('>' + lab)
            parts = [row[i:i + wrap] for i in range(0, len(row), wrap)] if wrap else [row]
            out += [(' ' + p + '\t') if noise else p for p in parts]
    return nl.join(out) + nl


def fasta(m, crlf=False):
    nl = '\r\n' if crlf else '\n'
    return ''.join(f'>{lab}{nl}{row}{nl}' for lab, row in m)


def out(*rows):
    """Expected output: rows in label order, 60 characters per line."""
    return ''.join(f'>{lab}\n' + ''.join(row[i:i + 60] + '\n' for i in range(0, len(row), 60))
                   for lab, row in rows)


def rows(spec):
    """'a=A- b=C-' -> [('a', 'A-'), ('b', 'C-')]"""
    return [tuple(x.split('=', 1)) for x in spec.split()]


H01 = [rows('a=A b=C c=D d=E'), rows('a=A- b=C- c=-D d=-E')]
H02 = [rows('a=AV b=CW c=DY d=EH e=FK'),
       rows('e=-FK c=D-Y a=-AV d=E-H b=-CW'),
       rows('b=C-W d=-EH a=A-V e=F-K c=-DY')]
H03 = [rows('p=AL q=CM r=DN s=EP t=FQ u=GR'),
       rows('p=A-L q=C-M r=-DN s=-EP t=F-Q u=G-R'),
       rows('p=--AL q=--CM r=D--N s=E--P t=-F-Q u=-G-R')]
H08 = [rows('a=.AC b=GAC c=-AC'), rows('a=A.C b=GAC c=A-C')]
H03_OUT = out(*rows('p=-A-L q=-C-M r=--DN s=--EP t=F--Q u=G--R'))
H08_OUT = out(*rows('a=AC b=-C c=AC'))

FILES = {
    'H01.efa': efa(*H01),
    'H02.efa': efa(*H02),
    'H03.efa': efa(*H03),
    'H04.efa': efa(rows('a=K b=L c=M'), rows('a=K- b=L- c=-M')),
    'H05.efa': efa(rows('a=A-- b=CD- c=-E- d=--F'), rows('a=A-- b=CD- c=-E- d=--F')),
    'H06.efa': efa(rows('a=A-- b=-DE c=C-- x=-FG'), rows('a=A-- b=-DE c=C-- x=-FG')),
    'H07.efa': efa(rows('a=AC- b=AC- c=---'), rows('a=A-C b=A-C c=---')),
    'H08.efa': efa(*H08),
    'H09.efa': efa(rows('a=AC- b=AC-'), rows('a=AC. b=AC-')),
    'H10.efa': efa(rows('a=A. b=C. c=D. d=E.'), rows('a=A. b=C. c=.D d=.E')),
    'H11.efa': efa([('b', 'ac-t'), ('a', '-cgt')]),
    'H12.efa': efa([('a9', 'AW'), ('s 1', 'AW'), ('B', 'AW'), ('a10', 'AW'), ('_x', 'AW')]),
    'H13.efa': efa([('a', '*1'), ('b', '?~')]),
    'H14.efa': efa([('x', 'A' * 61), ('y', 'A' * 61)]),
    'H15a.efa': efa(*H03, crlf=True, wrap=2, noise=True),
    'H15b.efa': efa(*H08, crlf=True, wrap=2),
    'H15c.txt': ''.join(f'H15c/m{k}.fa\r\n' for k in range(3)),
    'H15d.txt': ''.join(f'H15d/m{k}.fa\n' for k in range(2)),
    'H16.efa': efa([('>x', 'AW'), ('y', 'AW')]),
    'E03.efa': '',
    'E04.txt': '\r',
    'E05.txt': '\n',
    'E06.txt': 'E06/m0.fa\nE06/m1.fa\n',
    'E06/m0.fa': '',
    'E06/m1.fa': '',
    'E07.efa': '<m0\n',
    'E08.efa': '<m0\n<m1\n>a\nAC\n>b\nAC\n',
    'E09.efa': efa(rows('a=AC b=AC'), rows('a=AC c=AC')),
    'E10.efa': efa(rows('a=AC'), rows('a=AC b=AC')),
    'E11.efa': efa(rows('a=AC b=AC'), rows('a=AG b=AC')),
    'E12.efa': efa(rows('a=AC a=AC'), rows('a=AC a=AC')),
    'E13.efa': efa(rows('a=AC b=AC'), rows('a=AC a=AC')),
    'E14.efa': efa(rows('a=AC b=A')),
    'E15b.efa': efa(rows('a=AC'), rows('a=A-C')),
    'E19.fa': fasta(rows('a=AC b=AC')),
    'expected/H01.fa': out(*rows('a=A- b=C- c=-D d=-E')),
    'expected/H02.fa': out(*rows('a=-AV b=-CW c=D-Y d=E-H e=-FK')),
    'expected/H03.fa': H03_OUT,
    'expected/H04.fa': out(*rows('a=K b=L c=-')),
    'expected/H05.fa': out(*rows('a=A- b=CD c=-E d=--')),
    'expected/H06.fa': out(*rows('a=A-- b=-DE c=C-- x=-FG')),
    'expected/H07.fa': out(*rows('a=AC b=AC c=--')),
    'expected/H08.fa': H08_OUT,
    'expected/H08m3.fa': out(*rows('a=C b=C c=C')),
    'expected/H09.fa': out(*rows('a=AC b=AC')),
    'expected/H11.fa': out(*rows('a=CT b=CT')),
    'expected/H12.fa': out(('B', 'AW'), ('_x', 'AW'), ('a10', 'AW'), ('a9', 'AW'), ('s 1', 'AW')),
    'expected/H13.fa': out(*rows('a=*1 b=?~')),
    'expected/H14.fa': out(('x', 'A' * 61), ('y', 'A' * 61)),
    'expected/H16.fa': out(('>x', 'AW'), ('y', 'AW')),
}
for k, m in enumerate(H03):
    FILES[f'H15c/m{k}.fa'] = fasta([(lab, row.lower()) for lab, row in m], crlf=True)
for k, m in enumerate(H08):
    FILES[f'H15d/m{k}.fa'] = fasta(m)

RUN = '-cloak {} -output {{out}}'
CASES = [  # id, args, exit, expected, v010, v020, what
    ('H01', RUN.format('H01.efa'), 0, 'expected/H01.fa', 'fail', 'fail', 'divvying, -mincol absent'),
    ('H01m0', RUN.format('H01.efa') + ' -mincol 0', 0, 'expected/H01.fa', 'fail', 'fail', '-mincol 0 means 2'),
    ('H01m2', RUN.format('H01.efa') + ' -mincol 2', 0, 'expected/H01.fa', 'fail', 'fail', '-mincol 2'),
    ('H02', RUN.format('H02.efa'), 0, 'expected/H02.fa', 'pass', 'pass', 'column order within an MSA-0 column'),
    ('H03', RUN.format('H03.efa'), 0, 'expected/H03.fa', 'pass', 'pass', 'order tie broken by a later MSA'),
    ('H04', RUN.format('H04.efa'), 0, 'expected/H04.fa', 'fail', 'fail', 'rows: every row keeps at most one residue'),
    ('H05', RUN.format('H05.efa'), 0, 'expected/H05.fa', 'fail', 'fail', 'rows: dropping left 1-residue columns'),
    ('H06', RUN.format('H06.efa'), 0, 'expected/H06.fa', 'fail', 'fail', 'rows: dropping left an all-gap column'),
    ('H07', RUN.format('H07.efa'), 0, 'expected/H07.fa', 'fail', 'fail', 'all-gap row and columns'),
    ('H08', RUN.format('H08.efa'), 0, 'expected/H08.fa', 'fail', 'fail', "'.' gaps, equal counts per row"),
    ('H08m3', RUN.format('H08.efa') + ' -mincol 3', 0, 'expected/H08m3.fa', 'fail', 'fail', "'.' gaps, -mincol 3"),
    ('H09', RUN.format('H09.efa'), 0, 'expected/H09.fa', 'pass', 'fail', "trailing '.' in MSA 1"),
    ('H10', RUN.format('H10.efa'), 0, 'expected/H01.fa', 'fail', 'fail', "'.' gaps everywhere, all-'.' column"),
    ('H11', RUN.format('H11.efa'), 0, 'expected/H11.fa', 'pass', 'pass', 'K = 1, lower case, label sort'),
    ('H12', RUN.format('H12.efa'), 0, 'expected/H12.fa', 'pass', 'pass', 'bytewise label order'),
    ('H13', RUN.format('H13.efa'), 0, 'expected/H13.fa', 'pass', 'pass', 'non-letter residues'),
    ('H14', RUN.format('H14.efa'), 0, 'expected/H14.fa', 'pass', 'pass', 'wrapping at 60'),
    ('H15a', RUN.format('H15a.efa'), 0, 'expected/H03.fa', 'pass', 'pass', 'CRLF, wrapped EFA with blanks'),
    ('H15b', RUN.format('H15b.efa'), 0, 'expected/H08.fa', 'fail', 'fail', "CRLF, wrapped EFA with '.'"),
    ('H15c', RUN.format('H15c.txt'), 0, 'expected/H03.fa', 'pass', 'pass', 'CRLF path list, lower case'),
    ('H15d', RUN.format('H15d.txt'), 0, 'expected/H08.fa', 'fail', 'fail', "path list with '.'"),
    ('H16', RUN.format('H16.efa'), 0, 'expected/H16.fa', 'fail', 'fail', "label starting with '>'"),
    ('E01', RUN.format('H03.efa') + ' -mincol 1', 1, '-', 'fail', 'fail', '-mincol 1'),
    ('E02', RUN.format('missing.efa'), 1, '-', 'pass', 'pass', 'missing input'),
    ('E03', RUN.format('E03.efa'), 1, '-', 'pass', 'pass', 'empty input file'),
    ('E04', RUN.format('E04.txt'), 1, '-', 'fail', 'fail', 'path list naming no MSA (CR only)'),
    ('E05', RUN.format('E05.txt'), 1, '-', 'pass', 'pass', 'path list with an empty line'),
    ('E06', RUN.format('E06.txt'), 1, '-', 'fail', 'pass', 'FASTA files with no records'),
    ('E07', RUN.format('E07.efa'), 1, '-', 'pass', 'pass', 'EFA with a name and no MSA'),
    ('E08', RUN.format('E08.efa'), 1, '-', 'pass', 'pass', 'EFA with two names, one MSA'),
    ('E09', RUN.format('E09.efa'), 1, '-', 'pass', 'pass', 'label sets differ'),
    ('E10', RUN.format('E10.efa'), 1, '-', 'pass', 'pass', 'MSA 1 has more rows than MSA 0'),
    ('E11', RUN.format('E11.efa'), 1, '-', 'pass', 'pass', 'ungapped sequences differ'),
    ('E12', RUN.format('E12.efa'), 1, '-', 'pass', 'pass', 'duplicate label in MSA 0'),
    ('E13', RUN.format('E13.efa'), 1, '-', 'pass', 'pass', 'duplicate label in MSA 1 only'),
    ('E14', RUN.format('E14.efa'), 1, '-', 'pass', 'pass', 'rows of unequal length'),
    ('E15a', RUN.format('H01.efa') + ' -mincol 3', 1, '-', 'pass', 'pass', 'no column survives'),
    ('E15b', RUN.format('E15b.efa'), 1, '-', 'pass', 'pass', 'one sequence'),
    ('E15c', RUN.format('H03.efa') + ' -mincol 4294967295', 1, '-', 'pass', 'pass', 'largest -mincol'),
    ('E17', '-cloak H03.efa -output nodir/out.fa', 1, '-', 'pass', 'pass', 'output directory missing'),
    ('E18', '-cloak', 1, '-', 'pass', 'pass', '-cloak without a value'),
    ('E19', RUN.format('E19.fa'), 1, '-', 'pass', 'pass', 'plain FASTA given as input'),
]


def write(outdir):
    os.makedirs(outdir)
    for name, text in FILES.items():
        path = os.path.join(outdir, name)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'wb') as f:
            f.write(text.encode())
    with open(os.path.join(outdir, 'MANIFEST.tsv'), 'wb') as f:
        f.write(b'id\targs\texit\texpected\tv010\tv020\twhat\n')
        for c in CASES:
            f.write('\t'.join(map(str, c)).encode() + b'\n')
    print(f'{outdir}: {len(FILES)} files, {len(CASES)} cases')


def verify(d, cloak_py=None):
    """Each case's input against spec_check.py (O1), oracle_direct.py (O2) and cloak.py (O3)."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import spec_check
    import ref_cloakpy
    from cloakio import LoadError, load, parse_mincol, render, slurp
    from oracle_direct import cloak
    ns = ref_cloakpy.load_reference(os.path.abspath(cloak_py)) if cloak_py else None
    os.chdir(d)
    bad = 0
    for cid, args, ex, exp, *_ in CASES:
        a = args.split()
        if cid in ('E17', 'E18'):
            continue            # argument errors, not input errors: only the binary can show them
        inp = a[a.index('-cloak') + 1]
        mc = a[a.index('-mincol') + 1] if '-mincol' in a else None
        try:
            data = cloak(*load(inp), parse_mincol(mc))
        except LoadError:
            data = None
        o2 = (data is None) if ex else (data == slurp(exp))
        if ex:
            try:
                o1 = not spec_check.expected(inp, parse_mincol(mc))[2]
            except LoadError:
                o1 = True
        else:
            o1 = not spec_check.check(inp, exp, parse_mincol(mc), order=True)[0]
        o3 = '-'
        if ns is not None and mc in (None, '0', '2'):
            try:
                labels, msas = load(inp)
            except LoadError:
                o3 = 'n/a'           # invalid input: cloak.py has no defined behaviour
            else:
                rows = ref_cloakpy.run(ns, msas) if labels else None
                if rows is None:
                    ok = ex == 1
                else:
                    tmp = os.path.join(tempfile.gettempdir(), f'cloak_o3_{os.getpid()}.fa')
                    with open(tmp, 'wb') as f:
                        f.write(render(labels, rows))
                    ok = ex == 0 and not spec_check.check(inp, tmp, 2)[0]
                    os.remove(tmp)
                o3 = 'ok' if ok else 'DISAGREES'
                bad += not ok
        if not (o1 and o2):
            bad += 1
        print(f'{cid:6} exit {ex}  O1 {"ok" if o1 else "DISAGREES"}  O2 {"ok" if o2 else "DISAGREES"}'
              f'  O3 {o3}')
    print(f'{len(CASES) - 2} cases checked, {bad} disagreements')
    return 1 if bad else 0


if __name__ == '__main__':
    if len(sys.argv) in (3, 4) and sys.argv[1] == '--verify':
        sys.exit(verify(sys.argv[2], sys.argv[3] if len(sys.argv) == 4 else None))
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    write(sys.argv[1])
