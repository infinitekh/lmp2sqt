#ifndef __lmp2sqt_h__ 
#define __lmp2sqt_h__ 
/*
 * =====================================================================================
 *
 *       Filename:  lmp2sqt.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2015년 11월 27일 15시 31분 18초
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

#include "common.h"
typedef struct {
	real R, I;
} Cmplx;
typedef struct { 
	real x,y,z;
} VecR3;
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "snapshot.h"
#define ALLOC(type) type*  alloc_ ## type(size_t n) { \
	return (type *) malloc(sizeof(type)*n); \
}
#define Min(a,b) ( (a<b)?(a):(b))
ALLOC(double);
ALLOC(real);
ALLOC(int);
ALLOC(Cmplx);
ALLOC(VecR3);
#define CSet(a, x, y) a.R = x, a.I = y
#define CAdd(a, b, c) a.R = b.R + c.R, a.I = b.I + c.I
#define CSub(a, b, c) a.R = b.R - c.R, a.I = b.I - c.I
#define CMul(a, b, c) a.R = b.R * c.R - b.I * c.I, a.I = b.R * c.I + b.I * c.R
#define Sqr(x)     ( (x)*(x))
#define CHAR_MINUS '-'
#define CHAR_ZERO  '0'’
#define NameVal(x)                     \
	if (! strncmp (bp, #x, strlen (#x))) { \
		bp += strlen (#x);                   \
		x = strtod (bp, &bp);                \
	}



#define CSet(a,x,y) 																		\
	a.R = x,																							\
	a.I = y
#define BUFF_LEN 1024
char *header[]= {"cur-long", "cur-trans", "densty"},
		 *txtCorr = "space-time corr";
void* unused_pointer;
int ununused_value;
#define AllocMem(a, n, t)                       \
	a = (t *) malloc ((n) * sizeof (t))     
#define AllocMem2(a, n1, n2, t)                        \
	AllocMem (a, n1, t *);                               \
AllocMem (a[0], (n1) * (n2), t);                     \
for (ununused_value = 1; ununused_value < n1; ununused_value ++) a[ununused_value] = a[ununused_value - 1] + n2;

void ZeroDiffuse ();
void InitDiffuse ();
void PrintDiffuse (FILE *fp);
void EvalDiffuse (Snapshot * snap);
void AccumDiffuse (int nCol);

void ZeroSpacetimeCorr ();
void InitSpacetimeCorr ();
void EvalOtherInformation ();
void PrintSpacetimeCorr (FILE *fp);
void EvalSpacetimeCorr (Snapshot * snap);
void AllocArray();
void Alloc_more(int);
void AccumSpacetimeCorr (int nMol);
real kVal, deltaT;

typedef struct {
	real **acfST, *orgST;
	VecR3 *orgR, *rTrue;
	// VecR3 *orgVel; real *acfVel;
	real *rrDiffuse;
	int count;
	int countDiff;
} TBuf;

TBuf *tBuf;
real **avAcfST, *valST, **valDqt, **valGammaQT;
int countCorrAv, limitCorrAv, nBuffCorr, nFunCorr, nValCorr;
real *rrDiffuseAv;
int countDiffuseAv, limitDiffuseAv, nBuffDiffuse, nValDiffuse;


#endif
