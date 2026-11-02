"""Tests of the test tools themselves: the definition checker must accept correct outputs and reject
each deliberate corruption, agree with the byte oracle on random ensembles, and the loader model
must follow MUSCLE's rules.

Usage: python test_spec_check.py      (standard library unittest; run from any directory)
"""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cloakio                                   # noqa: E402
import gen_synth                                 # noqa: E402
import make_hand_fixtures as hand                # noqa: E402
import oracle_direct                             # noqa: E402
import spec_check                                # noqa: E402


def rows_to_out(spec):
    return hand.out(*hand.rows(spec)).encode()


class Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def put(self, name, data):
        path = os.path.join(self.dir, name)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'wb') as f:
            f.write(data if isinstance(data, bytes) else data.encode())
        return path

    def problems(self, case, out, mincol=2, order=False):
        return spec_check.check(self.put(case, hand.FILES[case]), self.put('out.fa', out),
                                mincol, order)[0]


class CheckerAcceptsAndRejects(Base):
    H03 = hand.H03_OUT.encode()

    def test_correct_outputs_pass(self):
        for case, exp, mc in [('H01.efa', 'expected/H01.fa', 2), ('H03.efa', 'expected/H03.fa', 2),
                              ('H05.efa', 'expected/H05.fa', 2), ('H06.efa', 'expected/H06.fa', 2),
                              ('H07.efa', 'expected/H07.fa', 2), ('H08.efa', 'expected/H08.fa', 2),
                              ('H08.efa', 'expected/H08m3.fa', 3), ('H14.efa', 'expected/H14.fa', 2),
                              ('H16.efa', 'expected/H16.fa', 2)]:
            with self.subTest(case=case, mincol=mc):
                self.assertEqual(self.problems(case, hand.FILES[exp], mc, order=True), [])

    def test_corruptions_fail(self):
        bad = {
            'swap residues within a column': 'p=-C-L q=-A-M r=--DN s=--EP t=F--Q u=G--R',
            'move a residue to another column': 'p=-A-L q=-C-M r=--DN s=--EP t=-F-Q u=G--R',
            'drop a row': 'p=-A-L q=-C-M r=--DN t=F--Q u=G--R',
            'drop a column': 'p=-AL q=-CM r=--N s=--P t=F-Q u=G-R',
            'duplicate a column': 'p=-A-LL q=-C-MM r=--DNN s=--EPP t=F--QQ u=G--RR',
            'split a class': 'p=--A-L q=--C-M r=---DN s=---EP t=F---Q u=-G--R',
            'merge two classes': 'p=-AL q=-CM r=D-N s=E-P t=F-Q u=G-R',
            'change a letter': 'p=-A-K q=-C-M r=--DN s=--EP t=F--Q u=G--R',
            "add a '.'": 'p=.A-L q=-C-M r=--DN s=--EP t=F--Q u=G--R',
            'lower case': 'p=-a-L q=-C-M r=--DN s=--EP t=F--Q u=G--R',
            'rows out of label order': 'q=-C-M p=-A-L r=--DN s=--EP t=F--Q u=G--R',
            'all-gap column added': 'p=-A-L- q=-C-M- r=--DN- s=--EP- t=F--Q- u=G--R-',
        }
        for what, spec in bad.items():
            with self.subTest(what):
                self.assertNotEqual(self.problems('H03.efa', rows_to_out(spec)), [])

    def test_format_corruptions_fail(self):
        a61 = b'A' * 61
        for what, out in {'re-wrapped at 61': b'>x\n' + a61 + b'\n>y\n' + a61 + b'\n',
                          'CRLF': self.H03.replace(b'\n', b'\r\n'),
                          'no final LF': self.H03[:-1],
                          'blank line': self.H03.replace(b'\n>q', b'\n\n>q')}.items():
            with self.subTest(what):
                case = 'H14.efa' if what.startswith('re-wrapped') else 'H03.efa'
                self.assertNotEqual(self.problems(case, out), [])

    def test_column_order_only_with_order(self):
        swapped = rows_to_out('p=A--L q=C--M r=--DN s=--EP t=-F-Q u=-G-R')
        self.assertEqual(self.problems('H03.efa', swapped), [])
        self.assertNotEqual(self.problems('H03.efa', swapped, order=True), [])

    def test_old_row_dropping_output_fails(self):
        self.assertNotEqual(self.problems('H05.efa', rows_to_out('b=CD')), [])
        self.assertNotEqual(self.problems('H06.efa', rows_to_out('b=-DE x=-FG')), [])


class OraclesAgree(Base):
    def test_random_ensembles(self):
        n_ok = n_empty = 0
        for schedule, first in (('mixed', 1), ('mild', 1001)):
            for seed in range(first, first + 40):
                efa = os.path.join(self.dir, f'{schedule}_{seed}.efa')
                gen_synth.make(seed, *gen_synth.params(seed, schedule), efa)
                for mincol in (2, 3):
                    with self.subTest(efa=os.path.basename(efa), mincol=mincol):
                        _, _, want = spec_check.expected(efa, mincol)   # raises if formulations differ
                        try:
                            data = oracle_direct.cloak(*cloakio.load(efa), mincol)
                        except cloakio.LoadError:
                            self.assertEqual(want, {})
                            n_empty += 1
                            continue
                        out = self.put('o.fa', data)
                        self.assertEqual(spec_check.check(efa, out, mincol, order=True)[0], [])
                        n_ok += 1
        self.assertGreater(n_ok, 40)
        self.assertGreater(n_empty, 0)


class LoaderModel(Base):
    def load_case(self, name):
        for f, text in hand.FILES.items():
            self.put(f, text)
        return cloakio.load(name, cwd=self.dir)

    def test_rejected_inputs(self):
        for name in ('missing.efa', 'E03.efa', 'E04.txt', 'E05.txt', 'E07.efa', 'E08.efa',
                     'E09.efa', 'E10.efa', 'E11.efa', 'E12.efa', 'E13.efa', 'E14.efa', 'E19.fa'):
            with self.subTest(name):
                with self.assertRaises(cloakio.LoadError):
                    self.load_case(name)

    def test_forms_give_the_same_ensemble(self):
        self.assertEqual(self.load_case('H15a.efa'), self.load_case('H03.efa'))
        self.assertEqual(self.load_case('H15c.txt'), self.load_case('H03.efa'))
        self.assertEqual(self.load_case('H15b.efa'), self.load_case('H15d.txt'))

    def test_fewer_rows_in_a_later_msa_is_undefined(self):
        with self.assertRaises(cloakio.Undefined):
            cloakio.load(self.put('few.efa', hand.efa(hand.rows('a=AC b=AC'), hand.rows('a=AC'))))

    def test_no_records_gives_no_rows(self):
        self.assertEqual(self.load_case('E06.txt'), ([], [[], []]))

    def test_lines_and_mincol(self):
        self.assertEqual(cloakio.lines_of(b''), [])
        self.assertEqual(cloakio.lines_of(b'\r'), [])
        self.assertEqual(cloakio.lines_of(b'\n'), [b''])
        self.assertEqual(cloakio.lines_of(b'a\r\nb'), [b'a', b'b'])
        self.assertEqual([cloakio.parse_mincol(x) for x in (None, '0', '2', '7', '4294967296')],
                         [2, 2, 2, 7, 2])
        for x in ('1', '4294967297', '-1', '+3', 'abc', ''):
            with self.subTest(x):
                with self.assertRaises(cloakio.LoadError):
                    cloakio.parse_mincol(x)

    def test_strip_rows(self):
        full = rows_to_out('a=A- b=CD c=-E d=--')
        self.assertEqual(cloakio.strip_rows(full), rows_to_out('b=CD'))
        self.assertEqual(cloakio.strip_rows(b'>a\nA\n'), b'')
        self.assertEqual(cloakio.dots_to_dashes(b'<m.0\n>a.b\n.A.\n'), b'<m.0\n>a.b\n-A-\n')


if __name__ == '__main__':
    unittest.main()
