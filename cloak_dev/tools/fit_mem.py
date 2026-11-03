"""Least-squares fit of peak RSS = c0 + a*cells + b*(n*out_cols) over measurement TSVs (ANALYSIS.md,
section 7; scripts/test_perf.sh writes such a TSV).

cells = sum over MSAs of rows x columns; out_cols = columns of the CLOAK output.
Usage: fit_mem.py TSV...   (columns named n, cells, out_cols, peak_rss_mb, status; rows whose status
is not OK are ignored). Prints the coefficients and each point's relative error.
"""
import csv
import sys


def solve3(A, y):
    # Normal equations by Gaussian elimination (3 unknowns).
    M = [[sum(A[i][r] * A[i][c] for i in range(len(A))) for c in range(3)] +
         [sum(A[i][r] * y[i] for i in range(len(A)))] for r in range(3)]
    for p in range(3):
        piv = max(range(p, 3), key=lambda r: abs(M[r][p]))
        M[p], M[piv] = M[piv], M[p]
        for r in range(3):
            if r != p:
                f = M[r][p] / M[p][p]
                M[r] = [a - f * b for a, b in zip(M[r], M[p])]
    return [M[r][3] / M[r][r] for r in range(3)]


def main():
    rows = []
    for path in sys.argv[1:]:
        with open(path) as f:
            for d in csv.DictReader(f, delimiter='\t'):
                if d['status'] != 'OK':
                    continue
                rows.append((path.split('/')[-1], int(d['n']), int(d['cells']), int(d['out_cols']),
                             float(d['peak_rss_mb']) * 1048576))
    if len(rows) < 3:
        sys.exit(f'fit_mem: {len(rows)} usable points, at least 3 needed')
    A = [[1.0, c, n * oc] for _, n, c, oc, _ in rows]
    y = [r[4] for r in rows]
    c0, a, b = solve3(A, y)
    print(f'points={len(rows)}  c0={c0 / 1e6:.1f} MB  a={a:.2f} B/cell  b={b:.2f} B/(row*out_col)')
    worst = 0
    for (src, n, c, oc, m), row in zip(rows, A):
        pred = c0 + a * c + b * n * oc
        err = (pred - m) / m
        worst = max(worst, abs(err))
        print(f'{src:12s} n={n:5d} cells={c:10d} out_cols={oc:6d} rss={m / 1e6:8.1f} MB '
              f'pred={pred / 1e6:8.1f} MB err={err:+.1%}')
    print(f'worst relative error {worst:.1%}')


if __name__ == '__main__':
    main()
