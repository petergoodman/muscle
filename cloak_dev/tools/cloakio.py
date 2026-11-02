"""Python model of MUSCLE's ensemble loader and of CLOAK's FASTA output, shared by the test tools.

Test tooling only (decision 001). Loading follows MUSCLE's own code, so that the tools see the
ensemble the binary sees:
  - The first byte decides the form: '<' is an EFA, anything else a list of FASTA paths
    (Ensemble::FromFile). An empty file stops.
  - EFA (Ensemble::FromEFA, MSA::FromStrings): every CR is deleted and lines split at LF; '<name'
    starts an MSA, '>label' a row whose label is the rest of the line; every other byte of a
    sequence line except whitespace is kept; rows of an MSA must be of equal length.
  - Path list (Ensemble::FromMSAPaths, GetFastaSeq): one path per line, relative to the working
    directory. In each FASTA file a label runs to the first CR or LF; in sequence data, gaps, ASCII
    letters (uppercased) and digits are kept, every other byte is dropped; a record that keeps
    nothing is skipped.
  - Then ASCII letters are uppercased, labels must be unique and the same in every MSA, and the
    sequences with gaps removed must agree (Ensemble::SetDerived).
Gaps are '-' and '.' (MUSCLE's isgap; decision 008). LoadError stands for the binary stopping with
exit 1; Undefined for input on which MUSCLE's loader has undefined behaviour (unequal row lengths in
a path-list FASTA, an MSA with fewer rows than MSA 0, NUL bytes), which tests must not use.

Command line (transforms for the test scripts):
  cloakio.py strip IN.fa OUT.fa   drop records with at most one residue (decision 003's relation)
  cloakio.py dots IN OUT          turn '.' into '-' in the sequence lines of an EFA or FASTA file
"""
import os
import sys

GAPS = b'-.'
SPACE = b' \t\n\v\f\r'
WIDTH = 60


class LoadError(Exception):
    pass


class Undefined(Exception):
    pass


def lines_of(data):
    """ReadLineStdioFile: every CR deleted, split at LF, a last line without LF kept if not empty."""
    out = data.replace(b'\r', b'').split(b'\n')
    if out[-1] == b'':
        out.pop()
    return out


def _from_strings(strings):
    """MSA::FromStrings and FromStrings2: one EFA block to [(label, row)]."""
    if not strings:
        raise LoadError('MSA with no data')
    labels, seqs, cur = [], [], bytearray()
    for s in strings:
        if s[:1] == b'>':
            if labels:
                seqs.append(bytes(cur))
            labels.append(s[1:])
            cur = bytearray()
        else:
            cur.extend(c for c in s if c not in SPACE)
    seqs.append(bytes(cur))
    if len(seqs) != len(labels):
        raise LoadError('labels and sequences differ in number')
    if len({len(q) for q in seqs}) > 1:
        raise LoadError('MSA not aligned')
    return list(zip(labels, seqs))


def _read_efa(data):
    names, blocks, cur = [], [], []
    for s in lines_of(data):
        if s[:1] == b'<':
            if cur:
                blocks.append(cur)
                cur = []
            names.append(s[1:])
        else:
            cur.append(s)
    blocks.append(cur)
    msas = [_from_strings(b) for b in blocks]
    if len(msas) != len(names):
        raise LoadError('MSA count differs from name count')
    return msas


def _read_fasta(data):
    """GetFastaSeq, record by record."""
    recs, i, n = [], 0, len(data)
    while i < n:
        if data[i] != 0x3E:
            raise LoadError('expected > to start a FASTA record')
        j = i + 1
        while j < n and data[j] not in (10, 13):
            j += 1
        if j == n:
            raise LoadError('end of file in FASTA label')
        label, i, prev, seq = data[i + 1:j], j + 1, 10, bytearray()
        while i < n:
            c = data[i]
            if c == 0x3E:
                if prev in (10, 13):
                    break
                raise LoadError('> inside FASTA sequence data')
            i += 1
            if c in SPACE:
                pass
            elif c in GAPS or 0x30 <= c <= 0x39 or 0x41 <= c <= 0x5A:
                seq.append(c)
            elif 0x61 <= c <= 0x7A:
                seq.append(c - 32)
            else:
                continue   # dropped with a warning; does not count as the previous character
            prev = c
        if seq:
            recs.append((label, bytes(seq)))
    if len({len(q) for _, q in recs}) > 1:
        raise Undefined('rows of unequal length in a path-list FASTA')
    return recs


def slurp(path):
    with open(path, 'rb') as f:
        return f.read()


def _read(path):
    try:
        data = slurp(path)
    except OSError as e:
        raise LoadError(f'cannot open {path}: {e.strerror}')
    if b'\0' in data:
        raise Undefined(f'NUL byte in {path}')
    return data


def load(path, cwd='.'):
    """Returns (labels, msas): labels sorted bytewise, msas[k][r] the uppercased row of labels[r]."""
    data = _read(os.path.join(cwd, path))
    if not data:
        raise LoadError('empty input file')
    if data[0] == 0x3C:
        msas = _read_efa(data)
    else:
        paths = lines_of(data)
        if not paths:
            raise LoadError('path list names no MSA')
        msas = [_read_fasta(_read(os.path.join(cwd, os.fsdecode(p)))) for p in paths]
    if any(len(m) < len(msas[0]) for m in msas[1:]):
        # Ensemble::ToUpper loops over MSA 0's row count in every MSA, before any count check.
        raise Undefined('an MSA has fewer rows than MSA 0')
    dicts = []
    for m in msas:
        d = {lab: row.upper() for lab, row in m}
        if len(d) != len(m):
            raise LoadError('duplicate label')
        dicts.append(d)
    msas = dicts
    labels = sorted(msas[0])
    for m in msas[1:]:
        if sorted(m) != labels:
            raise LoadError('label sets differ')
    for lab in labels:
        if len({ungap(m[lab]) for m in msas}) > 1:
            raise LoadError(f'sequences differ: {lab!r}')
    return labels, [[m[lab] for lab in labels] for m in msas]


def ungap(row):
    return bytes(c for c in row if c not in GAPS)


def parse_mincol(arg):
    """-mincol through MUSCLE's StrToUint (wraps modulo 2**32), then cmd_cloak's default and check."""
    if arg is None:
        return 2
    if not arg.isascii() or not arg.isdigit():
        raise LoadError(f'invalid integer {arg!r}')
    v = int(arg) % 2**32
    if v == 1:
        raise LoadError('-mincol must be at least 2')
    return v or 2


def render(labels, rows):
    """FASTA as CLOAK writes it: '>label', then the row in lines of 60."""
    out = bytearray()
    for lab, row in zip(labels, rows):
        out += b'>' + lab + b'\n'
        for i in range(0, len(row), WIDTH):
            out += row[i:i + WIDTH] + b'\n'
    return bytes(out)


def records(data):
    """Splits FASTA bytes into (label, sequence, raw record bytes)."""
    recs = []
    for chunk in data.split(b'\n>'):
        if not recs and not chunk.startswith(b'>'):
            raise ValueError('FASTA does not start with >')
        raw = chunk if not recs else b'>' + chunk
        head, _, body = raw[1:].partition(b'\n')
        recs.append([head, body.replace(b'\n', b''), raw])
    for r in recs[:-1]:
        r[2] += b'\n'
    return [tuple(r) for r in recs]


def strip_rows(data):
    """Decision 003's relation: the records that keep at least two residues, bytes unchanged."""
    return b''.join(raw for _, seq, raw in records(data) if len(seq) - seq.count(b'-') >= 2)


def dots_to_dashes(data):
    """'.' made '-' in every line that is neither an MSA name nor a label."""
    return b'\n'.join(s if s[:1] in (b'<', b'>') else s.replace(b'.', b'-')
                      for s in data.split(b'\n'))


def main(argv):
    if len(argv) != 4 or argv[1] not in ('strip', 'dots'):
        sys.exit(__doc__)
    data = slurp(argv[2])
    out = strip_rows(data) if argv[1] == 'strip' else dots_to_dashes(data)
    with open(argv[3], 'wb') as f:
        f.write(out)


if __name__ == '__main__':
    main(sys.argv)
