"""Convert gzipped Stockholm to unaligned FASTA (gaps removed, uppercase). Usage: sto2fa.py in.sto.gz out.fa"""
import gzip, sys, collections
seqs = collections.OrderedDict()
with gzip.open(sys.argv[1], 'rt', errors='replace') as f:
    for line in f:
        if not line.strip() or line.startswith('#') or line.startswith('//'):
            continue
        name, s = line.split()[:2]
        seqs[name] = seqs.get(name, '') + s
with open(sys.argv[2], 'w') as o:
    for n, s in seqs.items():
        u = ''.join(c for c in s if c not in '.-').upper()
        o.write(f'>{n}\n{u}\n')
print(sys.argv[2], len(seqs), 'seqs')
