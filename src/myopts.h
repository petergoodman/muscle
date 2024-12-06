#ifndef MY_VERSION
#define MY_VERSION	"5.3"
#endif

#define PROGRAM_NAME	"muscle"

////////////////////
// Commands
#define C(x)	STR_OPT(x)
#include "cmds.h"
////////////////////

STR_OPT(log)
STR_OPT(output)
STR_OPT(output1)
STR_OPT(output2)
STR_OPT(output3)
STR_OPT(output4)
STR_OPT(query)
STR_OPT(input)
STR_OPT(input2)
STR_OPT(joinprefix)
STR_OPT(joinpaths)
STR_OPT(linkage)
STR_OPT(joins)
STR_OPT(tsvout)
STR_OPT(html)
STR_OPT(jalview)
STR_OPT(distmxin)
STR_OPT(guidetreein)
STR_OPT(guidetreeout)
STR_OPT(prefix)
STR_OPT(suffix)
STR_OPT(nodes)
STR_OPT(label)
STR_OPT(labels2)
STR_OPT(savedir)
STR_OPT(db)
STR_OPT(label1)
STR_OPT(label2)
STR_OPT(subtreeout)
STR_OPT(supertreeout)
STR_OPT(refmsa)
STR_OPT(ref)
STR_OPT(refdir)
STR_OPT(indir)
STR_OPT(testdir)
STR_OPT(testdir1);
STR_OPT(testdir2);
STR_OPT(outdir)
STR_OPT(hmmin)
STR_OPT(hmmout)
STR_OPT(report)
STR_OPT(accalnout)
STR_OPT(perm)
STR_OPT(calnout)
STR_OPT(anchor_letter)
STR_OPT(centroids)
STR_OPT(substmx)
STR_OPT(fev)
STR_OPT(gridspec)
STR_OPT(spatterspec)
STR_OPT(kmerdist)

UNS_OPT(threads)
UNS_OPT(consiters)
UNS_OPT(refineiters)
UNS_OPT(randseed)
UNS_OPT(paircount)
UNS_OPT(n)
UNS_OPT(splitcount)
UNS_OPT(maxcoarse)
UNS_OPT(perturb)
UNS_OPT(replicates)
UNS_OPT(maxcols)
UNS_OPT(minsuper)
//UNS_OPT(muloccs)
//UNS_OPT(normaafreqs)
UNS_OPT(minpctid)
UNS_OPT(maxpctid)
UNS_OPT(maxiters)
UNS_OPT(maxfailiters)
UNS_OPT(triesperiter)
UNS_OPT(blosumpct)
UNS_OPT(blosumparamset)
UNS_OPT(warmup_pct)
UNS_OPT(treeiters)
UNS_OPT(shrub_size)

FLT_OPT(min_cons_pct)
FLT_OPT(max_gap_fract)
FLT_OPT(max_gap_fract_row)
FLT_OPT(minea)
FLT_OPT(super6_maxpd1)
FLT_OPT(super5_minea1)
FLT_OPT(super4_minea1)
FLT_OPT(super4_minea2)
FLT_OPT(pctid)
FLT_OPT(perturb_var)
FLT_OPT(minconf)
FLT_OPT(maxpd)
FLT_OPT(shrink)

FLT_OPT(gapopen)
FLT_OPT(gapext)
FLT_OPT(termgapopen)
FLT_OPT(termgapext)
FLT_OPT(center)

FLT_OPT(s_is)
FLT_OPT(s_il)
FLT_OPT(m_is)
FLT_OPT(m_il)
FLT_OPT(is_is)
FLT_OPT(il_il)

FLAG_OPT(quiet)
FLAG_OPT(compilerinfo)
FLAG_OPT(right)
FLAG_OPT(scaledist)
FLAG_OPT(eadist)
FLAG_OPT(force_super4)
FLAG_OPT(force_probcons)
FLAG_OPT(allpairs)
FLAG_OPT(nt)
FLAG_OPT(amino)
FLAG_OPT(accs)
FLAG_OPT(verbose)
FLAG_OPT(basename)
FLAG_OPT(intsuffix)
FLAG_OPT(stratified)
FLAG_OPT(diversified)
FLAG_OPT(randomchaintree)
FLAG_OPT(input_order)
FLAG_OPT(tree_order)
FLAG_OPT(muscle3_randomorder)
FLAG_OPT(confseq1)
FLAG_OPT(q2)
FLAG_OPT(missingtestseqok)
FLAG_OPT(missingtestfileok)
FLAG_OPT(bysequence)
FLAG_OPT(reseek)
FLAG_OPT(mega)
FLAG_OPT(squeeze)

#undef FLAG_OPT
#undef UNS_OPT
#undef FLT_OPT
#undef STR_OPT
