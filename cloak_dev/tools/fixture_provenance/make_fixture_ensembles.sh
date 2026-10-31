#!/bin/bash
# Recipe that produced the Pfam-derived test ensembles in cloak_dev/tests/data (2026-10-27, in a
# Claude Cowork sandbox: Ubuntu 24.04, g++ 13.3, 2 cores). Kept so the fixtures are reproducible;
# rerunning it is not part of any test. muscle -align and -super5 are not guaranteed to reproduce
# the same MSAs on another machine or thread count, so regenerated ensembles may differ.
#
# Usage: make_fixture_ensembles.sh MUSCLE_F3559AF WORKDIR
#   MUSCLE_F3559AF  muscle binary built from rcedgar/muscle f3559af
#   WORKDIR         empty directory; raw/, fa/, sub/, ens/ are created inside
#
# Sources: Pfam seed alignments for PF26345, PF00188, PF00353, PF00745, PF17862 and Pfam full
# alignments for PF00188 and PF00353, from the InterPro API (Pfam release not recorded):
#   https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/<PF>/?annotation=alignment:seed
#   https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/<PF>/?annotation=alignment:full
# Pfam data are CC0.
set -euo pipefail
M=$1; W=$2; T=$(cd "$(dirname "$0")" && pwd)
PY="${PY:-python3}"
mkdir -p "$W"/{raw,fa,sub,ens}; cd "$W"

for pf in PF26345 PF00188 PF00353 PF00745 PF17862; do
  curl -fsS -o raw/$pf.seed.sto.gz "https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/$pf/?annotation=alignment:seed"
  $PY "$T/sto2fa.py" raw/$pf.seed.sto.gz fa/$pf.fa
done
for pf in PF00188 PF00353; do
  curl -fsS -o raw/$pf.full.sto.gz "https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/$pf/?annotation=alignment:full"
  $PY "$T/sto2fa.py" raw/$pf.full.sto.gz fa/${pf}_full.fa
done

strat() {  # family n seed: stratified ensemble (16 MSAs) of a random subsample
  $PY "$T/subsample.py" fa/$1.fa $2 $3 sub/${1}_n$2.fa
  $M -align sub/${1}_n$2.fa -stratified -output ens/${1}_n$2.efa -threads 2
}
super5() {  # input_fasta out_name: 4 x super5 -perturb 0..3 -perm all, the recipe for large families
  out=ens/$2.efa; : > "$out"
  for p in 0 1 2 3; do
    $M -super5 "$1" -perturb $p -perm all -output tmp_$2_@.$p.afa -threads 2
    for g in tmp_$2_*.$p.afa; do echo "<$g" >> "$out"; cat "$g" >> "$out"; rm "$g"; done
  done
}

strat PF26345 9 1;  strat PF26345 20 1; strat PF00188 20 1; strat PF00353 40 1
strat PF00745 40 1; strat PF17862 40 1; strat PF00188 60 1; strat PF00745 80 1
strat PF00353 80 1; strat PF17862 80 1; strat PF00188 92 1
for x in "PF17862 532" "PF00353 205" "PF00745 387" "PF00188 92"; do
  set -- $x; $PY "$T/subsample.py" fa/$1.fa $2 1 sub/${1}_n$2.fa; super5 sub/${1}_n$2.fa ${1}_n$2_s5
done
$PY "$T/subsample.py" fa/PF00188_full.fa 203 7 sub/PF00188full_n203.fa;   super5 sub/PF00188full_n203.fa PF00188full_n203_s5
$PY "$T/subsample.py" fa/PF00353_full.fa 1226 7 sub/PF00353full_n1226.fa; super5 sub/PF00353full_n1226.fa PF00353full_n1226_s5
