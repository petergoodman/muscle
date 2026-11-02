"""Runs random ensembles through a CLOAK binary against the test oracles: the random part of the
commit and release tiers (scripts/test.sh). Test tooling only (decision 001).

Cases: gen_cases.py schedules and gen_synth.py's mixed and mild, each at -mincol 2 and 3. For each:
  O1  spec_check: CLOAK exits 1 exactly when no class reaches mincol, and otherwise writes an
      output that meets the definition, in CLOAK's column order;
  O2  oracle_direct: the output equals the byte oracle's;
  O4  (--baseline B; EFA cases the baseline can finish, pre-screened with predict_recursion.py)
      decision 003's relation: the output without rows of at most one residue equals the
      baseline's output on the input with '.' made '-' (the baseline exits 1 if nothing remains);
  O3  (--cloakpy P; at most 60 residues; mincol 2) cloak.py's functions give the same rows and
      columns, compared as sets.
On --relations N cases (mincol 2) it checks the metamorphic relations of tests/README.md:
R1 row shuffle, R2 all-gap columns, R3 '-' to '.', R4 lower case, R5 CRLF and wrapping, R6 path
list, R8 rotation, R9 permutation of MSAs 1..K-1, R10 a copied MSA, R11 higher mincol, R12 K = 1,
R13 a deleted row, R14 reversed rows, R15 concatenation, R16 -threads and -log, R17 letters renamed.
A case the baseline or cloak.py cannot finish within its cap is skipped by name; any skip fails the
run unless --allow-skips. Coverage counters are printed; --assert-coverage fails the run below
minima scaled from tests/README.md (per 300 gen_cases cases).

Usage: run_suite.py BINARY OUTDIR [--count N] [--first SEED] [--schedules fission,rogue,...]
         [--synth N] [--baseline B] [--cloakpy P] [--relations N] [--allow-skips]
         [--assert-coverage] [--expect-fail]
--expect-fail inverts the exit status (for recording the red step on an old binary).
"""
import argparse
import os
import random
import resource
import signal
import subprocess
import sys
import time
from collections import Counter, defaultdict

import cloakio
import gen_cases
import gen_synth
import oracle_direct
import spec_check
from cloakio import GAPS, LoadError, load, render, strip_rows

HERE = os.path.dirname(os.path.abspath(__file__))
# Per 300 cases: about 60% of what --count 60 --synth 50 produced on 2026-11-02 (scaled from 400
# cases), so that a generator change that stops producing these cases fails the run.
MINIMA = {'multi2': 250, 'multi3': 180, 'tie': 200, 'row0': 120, 'row1': 70, 'dot0': 15, 'dotk': 10,
          'K1': 10, 'K30': 4, 'n1': 7, 'empty': 50, 'pathlist': 25, 'crlf': 30}


def run(binary, args, cwd, cap_gb=2, timeout=60):
    """Exit status (128 + signal for a signal, 'timeout') and wall time."""
    def limits():
        resource.setrlimit(resource.RLIMIT_AS, (int(cap_gb * 2**30),) * 2)
        resource.setrlimit(resource.RLIMIT_FSIZE, (200 * 2**20,) * 2)
    t0 = time.time()
    cmd = (['bash', binary] if binary.endswith('.sh') else [binary]) + args
    try:
        p = subprocess.run(cmd, cwd=cwd, preexec_fn=limits, timeout=timeout,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        rc = p.returncode if p.returncode >= 0 else 128 - p.returncode
    except subprocess.TimeoutExpired:
        rc = 'timeout'
    return rc, time.time() - t0


class Suite:
    def __init__(self, a):
        self.a, self.binary = a, os.path.abspath(a.binary)
        self.out = os.path.abspath(a.outdir)
        self.cases_dir = os.path.join(self.out, 'cases')
        self.tmp = os.path.join(self.out, 'tmp')
        os.makedirs(self.cases_dir)
        os.makedirs(self.tmp)
        os.chdir(self.cases_dir)      # path lists name their files relative to the cases directory
        self.fail, self.skips, self.n = [], [], Counter()
        self.cov = Counter()
        self.ns = None
        if a.cloakpy:
            import ref_cloakpy
            self.ns = ref_cloakpy.load_reference(a.cloakpy)

    def cloak(self, binary, inp, mincol=None, extra=(), cwd=None, cap_gb=2, timeout=60):
        """Runs binary on inp; returns (exit status, output bytes or None)."""
        outp = os.path.join(self.tmp, 'out.fa')
        if os.path.exists(outp):
            os.remove(outp)
        args = ['-cloak', inp, '-output', outp] + (['-mincol', str(mincol)] if mincol else [])
        rc, _ = run(binary, args + list(extra), cwd or self.cases_dir, cap_gb, timeout)
        return rc, cloakio.slurp(outp) if os.path.exists(outp) else None

    def bad(self, case, what):
        self.fail.append(f'{case}: {what}')
        print(f'FAIL {case}: {what}', flush=True)

    def write_efa(self, name, labels, msas, nl='\n'):
        path = os.path.join(self.tmp, name)
        lines = []
        for k, m in enumerate(msas):
            lines.append(f'<m{k}')
            for lab, row in zip(labels, m):
                lines += ['>' + lab.decode('latin-1'), row.decode('latin-1')]
        with open(path, 'wb') as f:
            f.write((nl.join(lines) + nl).encode('latin-1'))
        return path

    # ---- one case
    def case(self, inp, rel_n):
        path = os.path.join(self.cases_dir, inp)
        try:
            labels, msas = load(path)
        except cloakio.Undefined as e:
            self.bad(inp, f'generator produced undefined input ({e})')
            return
        raw = cloakio.slurp(path)
        self.coverage(inp, labels, msas, raw)
        self.pred = None
        for mincol in (2, 3):
            cid = f'{inp} m{mincol}'
            self.n['comparisons'] += 1
            rc, out = self.cloak(self.binary, inp, mincol)
            try:
                _, _, want = spec_check.expected(path, mincol)
            except LoadError:
                want = {}
            if mincol == 3:
                self.cov['multi3'] += multi(want)
            exp_rc = 0 if want else 1
            if rc != exp_rc:
                self.bad(cid, f'O1: exit {rc}, expected {exp_rc}')
                continue
            o1_ok = True
            if rc == 0:
                outp = os.path.join(self.tmp, 'out.fa')
                problems = spec_check.check(path, outp, mincol, order=True)[0]
                if problems:
                    self.bad(cid, 'O1: ' + '; '.join(problems))
                    o1_ok = False
            try:
                data = oracle_direct.cloak(labels, msas, mincol)
            except LoadError:
                data = None
            if data != out:
                self.bad(cid, 'O2: output differs from oracle_direct.py')
            if self.a.baseline and inp.endswith('.efa'):
                self.baseline(cid, path, out, labels, msas)
            if self.ns is not None and mincol == 2:
                self.cloakpy(cid, labels, msas, out)
            if rel_n and mincol == 2 and rc == 0 and o1_ok:   # predictions need a correct output
                self.relations(cid, labels, msas, out)

    def calls(self, labels, msas):
        """Recursive calls f3559af and cloak.py would make (predict_recursion.py), cached per case."""
        if self.pred is None:
            p = self.write_efa('pred.efa', labels, [[r.replace(b'.', b'-') for r in m] for m in msas])
            try:
                r = subprocess.run([sys.executable, os.path.join(HERE, 'predict_recursion.py'), p],
                                   capture_output=True, text=True, timeout=20)
                self.pred = float(r.stdout.split('calls=')[1].split()[0])
            except (subprocess.TimeoutExpired, IndexError, ValueError):
                self.pred = float('inf')
        return self.pred

    def baseline(self, cid, path, out, labels, msas):
        if self.calls(labels, msas) > 1e6:
            self.n['baseline not eligible'] += 1
            return
        dashed = os.path.join(self.tmp, 'dashed.efa')
        with open(dashed, 'wb') as f:
            f.write(cloakio.dots_to_dashes(cloakio.slurp(path).replace(b'\r', b'')))
        mincol = cid.split()[-1][1:]
        rc, b = self.cloak(self.a.baseline, dashed, mincol, cwd=self.tmp)
        if rc == 'timeout' or (isinstance(rc, int) and rc >= 128):
            self.skips.append(f'{cid}: baseline {rc}')
            return
        self.n['baseline compared'] += 1
        stripped = strip_rows(out) if out else b''
        if stripped and (rc, b) != (0, stripped):
            self.bad(cid, f'O4: baseline (exit {rc}) differs from the output minus short rows')
        elif not stripped and rc != 1:
            self.bad(cid, f'O4: baseline exit {rc}, expected 1 (no row would remain)')

    def cloakpy(self, cid, labels, msas, out):
        if sum(len(row) - sum(row.count(g) for g in GAPS) for row in msas[0]) > 60 \
                or (labels and self.calls(labels, msas) > 2e4):
            self.n['cloak.py not eligible'] += 1
            return
        import ref_cloakpy

        def alarm(*_):
            raise TimeoutError
        signal.signal(signal.SIGALRM, alarm)
        signal.alarm(20)
        try:
            rows = ref_cloakpy.run(self.ns, msas) if labels else None
        except (RecursionError, TimeoutError) as e:
            self.skips.append(f'{cid}: cloak.py {type(e).__name__}')
            return
        except Exception as e:
            self.bad(cid, f'O3: cloak.py raised {type(e).__name__}')
            return
        finally:
            signal.alarm(0)
        self.n['cloak.py compared'] += 1
        if (rows is None) != (out is None):
            self.bad(cid, 'O3: cloak.py and CLOAK disagree on whether any column survives')
        elif rows is not None and not same_sets(render(labels, rows), out):
            self.bad(cid, 'O3: columns or rows differ from cloak.py')

    # ---- metamorphic relations, on the normalised ensemble (labels sorted, upper case)
    def relations(self, cid, labels, msas, out0):
        rng = random.Random(cid)
        K = len(msas)
        self.n['relation cases'] += 1

        def same(rel, path, expect, mincol=2, extra=(), cwd=None, sets=False):
            self.n['relations'] += 1
            rc, out = self.cloak(self.binary, path, mincol, extra, cwd=cwd)
            ok = (out is None and expect is None) if (out is None or expect is None) else \
                (same_sets(out, expect) if sets else out == expect)
            if not ok:
                self.bad(cid, f'{rel}: exit {rc}, output differs from the predicted one')

        base = self.write_efa('base.efa', labels, msas)
        same('R0 normalised input', base, out0)
        order =[rng.sample(range(len(labels)), len(labels)) for _ in msas]
        lines = []
        for k, m in enumerate(msas):
            lines.append(b'<m%d' % k)
            for r in order[k]:
                lines += [b'>' + labels[r], m[r]]
        p = os.path.join(self.tmp, 'r1.efa')
        with open(p, 'wb') as f:
            f.write(b'\n'.join(lines) + b'\n')
        same('R1 row shuffle', p, out0)
        gapcols = [gap_columns(rng, m) for m in msas]
        same('R2 all-gap columns', self.write_efa('r2.efa', labels, gapcols), out0)
        dotted = [[bytes(c if c != 0x2D or rng.random() < 0.5 else 0x2E for c in row) for row in m]
                  for m in msas]
        same('R3 dots', self.write_efa('r3.efa', labels, dotted), out0)
        lower = [[row.lower() if rng.random() < 0.5 else row for row in m] for m in msas]
        same('R4 lower case', self.write_efa('r4.efa', labels, lower), out0)
        p = os.path.join(self.tmp, 'r5.efa')
        lines = []
        for k, m in enumerate(msas):
            lines.append(b'<m%d' % k)
            for lab, row in zip(labels, m):
                lines += [b'>' + lab] + [b' ' + row[i:i + 5] + b'\t' for i in range(0, len(row), 5)]
        with open(p, 'wb') as f:
            f.write(b'\r\n'.join(lines) + b'\r\n')
        same('R5 CRLF and wrapping', p, out0)
        if all(all(c in GAPS or chr(c).isalnum() for c in row) for m in msas for row in m) \
                and all(len(row) for m in msas for row in m):
            d = os.path.join(self.tmp, 'r6')
            os.makedirs(d, exist_ok=True)
            with open(os.path.join(d, 'list.txt'), 'wb') as f:
                for k, m in enumerate(msas):
                    f.write(b'm%d.fa\n' % k)
                    with open(os.path.join(d, f'm{k}.fa'), 'wb') as g:
                        g.write(b''.join(b'>' + lab + b'\n' + row + b'\n' for lab, row in zip(labels, m)))
            same('R6 path list', 'list.txt', out0, cwd=d)
        j = rng.randrange(K)
        same('R8 rotation', self.write_efa('r8.efa', labels, msas[j:] + msas[:j]), out0, sets=True)
        rest = msas[1:]
        rng.shuffle(rest)
        same('R9 permuted MSAs', self.write_efa('r9.efa', labels, msas[:1] + rest), out0, sets=True)
        j = rng.randrange(K)
        same('R10 copied MSA', self.write_efa('r10.efa', labels, msas[:j + 1] + [msas[j]] + msas[j + 1:]),
             out0)
        same('R11 mincol 3', base, drop_columns(out0, 3), mincol=3)
        rows0 = [row.replace(b'.', b'-') for row in msas[0]]
        same('R12 K = 1', self.write_efa('r12.efa', labels, msas[:1]),
             drop_columns(render(labels, rows0), 2))
        if len(labels) > 1:
            r = rng.randrange(len(labels))
            keep = [i for i in range(len(labels)) if i != r]
            recs = cloakio.records(out0)
            expect = drop_columns(render([recs[i][0] for i in keep], [recs[i][1] for i in keep]), 2)
            same('R13 deleted row', self.write_efa('r13.efa', [labels[i] for i in keep],
                                                   [[m[i] for i in keep] for m in msas]), expect)
        recs = cloakio.records(out0)
        same('R14 reversed rows', self.write_efa('r14.efa', labels, [[row[::-1] for row in m] for m in msas]),
             render([lab for lab, _, _ in recs], [seq[::-1] for _, seq, _ in recs]))
        same('R15 concatenation', self.write_efa('r15.efa', labels, [[row + row for row in m] for m in msas]),
             render([lab for lab, _, _ in recs], [seq + seq for _, seq, _ in recs]))
        same('R16 -threads 8 -log', base, out0, extra=('-threads', '8', '-log', os.path.join(self.tmp, 'log')))
        table = bytes.maketrans(b'ACDEFGHIKLMNPQRSTVWY', b'CDEFGHIKLMNPQRSTVWYA')
        same('R17 renamed letters', self.write_efa('r17.efa', labels, [[row.translate(table) for row in m] for m in msas]),
             render([lab for lab, _, _ in recs], [seq.translate(table) for _, seq, _ in recs]))

    def coverage(self, inp, labels, msas, raw):
        _, _, want = spec_check.expected(os.path.join(self.cases_dir, inp), 2)
        self.cov['cases'] += 1
        self.cov['multi2'] += multi(want)
        by0 = defaultdict(list)
        for t in want:
            by0[t[0]].append(t)
        self.cov['tie'] += sum(1 for ts in by0.values() if len(ts) != len({t[1:2] for t in ts}))
        kept = Counter(r for cls in want.values() for r, _ in cls)
        self.cov['row0'] += any(kept[r] == 0 for r in range(len(labels)))
        self.cov['row1'] += any(kept[r] == 1 for r in range(len(labels)))
        self.cov['dot0'] += any(0x2E in row for row in msas[0])
        self.cov['dotk'] += (not any(0x2E in row for row in msas[0])) and any(0x2E in row for m in msas[1:] for row in m)
        self.cov['K1'] += len(msas) == 1
        self.cov['K30'] += len(msas) == 30
        self.cov['n1'] += len(labels) == 1
        self.cov['empty'] += not want
        self.cov['pathlist'] += inp.endswith('.txt')
        self.cov['crlf'] += b'\r' in raw


def multi(want):
    return sum(1 for v in Counter(t[0] for t in want).values() if v > 1)


def gap_columns(rng, m):
    """Insert 1 to 3 all-gap columns at random places."""
    rows = [bytearray(r) for r in m]
    for _ in range(rng.randint(1, 3)):
        c = rng.randint(0, len(rows[0]) if rows else 0)
        for r in rows:
            r[c:c] = b'-'
    return [bytes(r) for r in rows]


def drop_columns(fasta, mincol):
    """The FASTA with columns of fewer than mincol residues removed (all rows kept); None if none
    remain."""
    if fasta is None:
        return None
    recs = cloakio.records(fasta)
    rows = [seq for _, seq, _ in recs]
    keep = [j for j in range(len(rows[0])) if sum(r[j] != 0x2D for r in rows) >= mincol] if rows else []
    if not keep:
        return None
    return render([lab for lab, _, _ in recs], [bytes(r[j] for j in keep) for r in rows])


def same_sets(a, b):
    """Same labels in the same order and the same multiset of columns (read down the rows)."""
    if a is None or b is None:
        return a is b
    ra, rb = cloakio.records(a), cloakio.records(b)
    if [x[0] for x in ra] != [x[0] for x in rb]:
        return False
    cols = lambda recs: Counter(bytes(s[j] for _, s, _ in recs) for j in range(len(recs[0][1]))) if recs else Counter()
    return cols(ra) == cols(rb)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('binary')
    ap.add_argument('outdir')
    ap.add_argument('--count', type=int, default=60, help='cases per gen_cases schedule')
    ap.add_argument('--first', type=int, default=1)
    ap.add_argument('--schedules', default=','.join(gen_cases.SCHEDULES))
    ap.add_argument('--synth', type=int, default=0, help='gen_synth.py cases per schedule (mixed, mild)')
    ap.add_argument('--baseline')
    ap.add_argument('--cloakpy')
    ap.add_argument('--relations', type=int, default=0, help='run the relations on every Nth case')
    ap.add_argument('--allow-skips', action='store_true')
    ap.add_argument('--assert-coverage', action='store_true')
    ap.add_argument('--expect-fail', action='store_true')
    a = ap.parse_args()
    s = Suite(a)
    names = []
    for sched in filter(None, a.schedules.split(',')):
        names += [gen_cases.write(s.cases_dir, sched, seed) for seed in range(a.first, a.first + a.count)]
    for sched, first in (('mixed', a.first), ('mild', a.first + 1000)):
        for seed in range(first, first + a.synth):
            name = f'{sched}_{seed}.efa'
            gen_synth.make(seed, *gen_synth.params(seed, sched), os.path.join(s.cases_dir, name))
            names.append(name)
    print(f'binary {a.binary}; {len(names)} cases in {s.cases_dir}', flush=True)
    for i, name in enumerate(names):
        s.case(name, a.relations and i % a.relations == 0)
    with open(os.path.join(s.out, 'skips.txt'), 'w') as f:
        f.write(''.join(x + '\n' for x in s.skips))
    counts = ', '.join(f'{k} {v}' for k, v in sorted(s.n.items()))
    print(f'{counts}; failed {len(s.fail)}; skipped {len(s.skips)}')
    print('coverage: ' + ', '.join(f'{k} {s.cov[k]}' for k in ['cases'] + list(MINIMA)))
    ok = not s.fail and (a.allow_skips or not s.skips)
    if a.assert_coverage:
        scale = s.cov['cases'] / 300
        low = [k for k, v in MINIMA.items() if s.cov[k] < v * scale]
        if low:
            print('coverage below minima: ' + ', '.join(f'{k} {s.cov[k]} < {MINIMA[k] * scale:.0f}' for k in low))
            ok = False
    with open(os.path.join(s.out, 'SUMMARY.txt'), 'w') as f:
        f.write(f'{counts}; failed {len(s.fail)}; skipped {len(s.skips)}\n')
    return (0 if not ok else 1) if a.expect_fail else (0 if ok else 1)


if __name__ == '__main__':
    sys.exit(main())
