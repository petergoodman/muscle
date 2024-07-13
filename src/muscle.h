#pragma once

#define _CRTDBG_MAP_ALLOC

#if	DEBUG && !_DEBUG
#define _DEBUG	1
#endif

#if	_DEBUG && !DEBUG
#define DEBUG	1
#endif

#if	_MSC_VER
#define TIMING	0
#endif

#ifdef	_MSC_VER	// Miscrosoft compiler
#pragma warning(disable : 4800)	// int-bool conversion
#pragma warning(disable : 4996)	// deprecated names like strdup, isatty.
#define brk(x)	if (x) __debugbreak()
#endif

#include "myutils.h"
#include "types.h"
#include "multisequence.h"
#include "textfile.h"
#include "mysparsemx.h"
#include "scoretype.h"
#include "treeperm.h"
#include "pairhmm.h"
#include "alpha.h"
#include "msa.h"
#include "mpcflat.h"
#include "kmerscan.h"
#include "alpha3.h"
#include "cachemem3.h"
#include "profile3.h"
#include "mx.h"
#include "mega.h"
#ifndef _WIN32
#define stricmp strcasecmp
#define strnicmp strncasecmp
#define	_snprintf snprintf
#define _fsopen(name, mode, share)	fopen((name), (mode))
#endif

const double VERY_NEGATIVE_DOUBLE = -9e29;
const float VERY_NEGATIVE_FLOAT = (float) -9e29;

void CalcEADistMx(FILE *f, MultiSequence* sequences,
  vector<vector<float> > &DistMx, vector<MySparseMx * > *SparsePostVec = 0);
void PermuteTree(const Tree &InputTree,
  Tree &TreeABC, Tree &TreeACB, Tree &TreeBCA,
  vector<string> &LabelsA, vector<string> &LabelsB, vector<string> &LabelsC);
void PermTree(Tree &InputTree, TREEPERM TP);
void StringsToFile(const string &FileName, const vector<string> &v);
void MakeReplicateFileName_N(const string &Pattern, uint N, string &FileName);
void MakeReplicateFileName(const string &Pattern, TREEPERM TP,
  uint PerturbSeed, string &FileName);
MultiSequence &LoadGlobalInputMS(const string &FileName);
MultiSequence &GetGlobalInputMS();
void ShowGlobalInputSeqStats();
double GetGlobalMSMeanSeqLength();
uint GetGlobalMSSeqCount();
uint GetGSICount();
uint GetAssertSameSeqsOkCount();
const Sequence &GetGlobalInputSeq(uint GSI);
const string &GetGlobalInputSeqLabel(uint GSI);
void ClearGlobalInputMS();
void CharVecToStr(const vector<char> &Vec, string &Str);
void LogAln(const Sequence &X, const Sequence &Y, const string &PathXY);
void LogAln(const string &X, const string &Y, const string &PathXY);
void ReadStringsFromFile(const string &FileName,
  vector<string> &Strings);
void GetGuideTreeJoinOrder(const Tree &GuideTree,
  const map<string, uint> &LabelToIndex,
  vector<uint> &Indexes1, vector<uint> &Indexes2);
void ValidateJoinOrder(const vector<uint> &Indexes1,
  const vector<uint> &Indexes2);

void _AssertSameLabels(const char *File, uint Line,
  const MultiSequence &MS);
void _AssertSameSeqs(const char *File, uint Line, 
  const MultiSequence &MS1, const MultiSequence &MS2);
void _AssertSameSeqsVec(const char *File, uint Line, 
  const MultiSequence &MS, vector<MultiSequence *> &v);
void _AssertSameSeqsJoin(const char *File, uint Line, 
  const MultiSequence &MS1, const MultiSequence &MS2, const MultiSequence &MS12);
void _AssertSeqsEqInput(const char *File, uint Line,
  const MultiSequence &MS);

void _AssertSeqsEq(const char *FileName, uint LineNr,
  const MultiSequence &MSA1, const MultiSequence &MSA2);
#define AssertSeqsEq(MSA1, MSA2)	_AssertSeqsEq(__FILE__, __LINE__, MSA1, MSA2)

#define AssertSeqsEqInput(MS)		_AssertSeqsEqInput(__FILE__, __LINE__, MS)
#define AssertSameLabels(MS)		_AssertSameLabels(__FILE__, __LINE__, MS)
#define AssertSameSeqs(MS1, MS2)	_AssertSameSeqs(__FILE__, __LINE__, MS1, MS2)
#define AssertSameSeqsVec(MS, v)	_AssertSameSeqsVec(__FILE__, __LINE__, MS, v)
#define AssertSameSeqsJoin(MS1, MS2, MS12)	_AssertSameSeqsJoin(__FILE__, __LINE__, MS1, MS2, MS12)

void LogFlatMx(const string &Name, const float *Flat, uint LX, uint LY);
void LogFlatMxs(const string &Name, const float *Flat, uint LX, uint LY);
void LogFlatMx1(const string &Name, const float *Flat, uint LX, uint LY);

float AlignMSAsFlat(const string &aProgressStr,
  const MultiSequence &MSA1, const MultiSequence &MSA2,
  uint TargetPairCount, string &Path);

void InitProbcons();
void AlignMSAsByPath(const MultiSequence &MSA1, const MultiSequence &MSA2,
  const string &Path, MultiSequence &MSA12);

void CalcFwdFlat(const byte *X, uint LX, const byte *Y, uint LY, float *Flat);
void CalcBwdFlat(const byte *X, uint LX, const byte *Y, uint LY, float *Flat);

void CalcFwdFlat_mega(const Mega &M, uint ProfileIdxX, uint ProfileIdxY, float *Flat);
void CalcBwdFlat_mega(const Mega &M,uint ProfileIdxX, uint ProfileIdxY, float *Flat);

void CalcPostFlat(const float *FlatFwd, const float *FlatBwd,
  uint LX, uint LY, float *Post);
float CalcAlnFlat(const float *Post, uint LX, uint LY,
  float *DPRows, char *TB, string &Path);
void CalcPosteriorFlat3(const MultiSequence &MSA1,
  const MultiSequence &MSA2,
  const vector<uint> &SeqIndexes1,
  const vector<uint> &SeqIndexes2,
  const vector<MySparseMx *> &SparseMxs,
  float *Flat);
float AlignPairFlat(const Sequence *Seq1, const Sequence *Seq2, string &Path);
float AlignPairFlat_SparsePost(const Sequence *Seq1, const Sequence *Seq2, string &Path,MySparseMx *SparsePost);
float AlignPairFlat_mega(const Mega * M, string &Path, uint index_X, uint index_Y);
float AlignPairFlat_mega_SparsePost(const Mega * M,string &Path, MySparseMx *SparsePost, uint index_X, uint index_Y);

void GetAllPairs(uint SeqCount,
  vector<uint> &SeqIndexes1, vector<uint> &SeqIndexes2);
void GetAllPairs(uint SeqCount1, uint SeqCount2,
  vector<uint> &SeqIndexes1, vector<uint> &SeqIndexes2);
void GetPairs(uint SeqCount1, uint SeqCount2, uint TargetPairCount,
  vector<uint> &SeqIndexes1, vector<uint> &SeqIndexes2);
float GetPostPairsAlignedFlat(const string &aProgressStr,
  const MultiSequence &MSA1, const MultiSequence &MSA2,
  const vector<uint> &SeqIndexes1, const vector<uint> &SeqIndexes2, 
  vector<MySparseMx *> &SparsePosts);
void SeqToFasta(FILE *f, const byte *Seq, unsigned L, const char *Label);
void SeqToFasta(FILE *f, const char *Seq, unsigned L, const char *Label);

float NWSmall3(CacheMem3 &CM, const Profile3 &ProfA,
  const Profile3 &ProfB, string &Path);
void AlignTwoMSAsGivenPath(const MultiSequence &msaA,
  const MultiSequence &msaB, const string &Path,
  MultiSequence &msa2);
void AlignTwoProfsGivenPath(const Profile3 &ProfA, float WeightA,
  const Profile3 &ProfB, float WeightB,
  const Mx2020 &SubstMx_Letter, float GapOpen,
  const string &Path, Profile3 &ProfAB);
void GetKimuraDistMx(const MultiSequence &MSA,
  vector<vector<float> > &DistMx);
void GetKimuraDistMx_Viterbi(const MultiSequence &MS,
  vector<vector<float> > &DistMx);
class PathInfo;
void LogAln(const byte *X, uint LX, const byte *Y, uint LY, const PathInfo &PI);
bool GetNextEnumGrid(const vector<uint> &Sizes, vector<uint> &Indexes);

void GetSubstMx_Letter_Blosum(uint PctId, float Mx[20][20]);
void ReadSubstMx_Letter_FromFile(const string &FileName, float Mx[20][20]);
void GetGapParams_Blosum(uint PctId, uint n, float *ptrOpen, float *ptrCenter);
float ScoreProfPos2(const ProfPos3 &PPA, const ProfPos3 &PPB);
double hscore(const double *xs, const double *ys, uint N, double X);

typedef bool (*ptr_GetMSAColIsAligned)(const MSA &Aln, uint Col);
uint GetOverlap(uint Lo1, uint Hi1, uint Lo2, uint Hi2);
