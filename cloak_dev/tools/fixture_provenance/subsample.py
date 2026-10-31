"""Random subsample of a FASTA. Usage: subsample.py in.fa n seed out.fa"""
import sys, random
recs=[]; name=None
for line in open(sys.argv[1]):
    line=line.rstrip('\n')
    if line.startswith('>'): recs.append([line[1:], '']); 
    else: recs[-1][1]+=line
n=int(sys.argv[2]); seed=int(sys.argv[3]); random.seed(seed)
sub = recs if n>=len(recs) else random.sample(recs, n)
with open(sys.argv[4],'w') as o:
    for a,b in sub: o.write(f'>{a}\n{b}\n')
