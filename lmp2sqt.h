/*!
 *    \file  lmp2sqt.h
 *   \brief  
 *
 *  Evaluate intermediate scattering function from lammps dump files.
 *    		\f$F_s(q,t)  \f$	(1.1)  self intermediate scattering function 
 *    														TAG name ptl1
 *
 *  \author  KIM Hyeok (kh), ekh0324@gmail.com
 *
 *  \internal
 *       Created:  2017년 05월 29일
 *      Revision:  none
 *      Compiler:  gcc
 *  Organization:  Konkuk University
 *     Copyright:  Copyright (c) 2017, KIM Hyeok
 *
 *  This source code is released for free distribution under the terms of the
 *  GNU General Public License as published by the Free Software Foundation.
 */

#ifndef __lmp2sqt_h__ 
#define __lmp2sqt_h__ 

/* #####   HEADER FILE INCLUDES   ################################################### */
#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "snapshot.h"
/* #####   EXPORTED TYPE DEFINITIONS   ############################################## */
typedef struct {
	real R, I;
} Cmplx;
/*!
 *  \brief  struct for complex real value
 */
typedef struct { 
	real x,y,z;
} VecR3;
typedef struct {
	real **acfFcol, *orgFcol;
	real **acfFself, *orgFself;
	VecR3 *orgR, *rTrue;
	// VecR3 *orgVel; real *acfVel;
	real *rrMSD;
	real *rrMQD;
	int **DrTable;
	int count;
	int countDiff;
} TBuf;


/* #####   EXPORTED FUNCTION DECLARATIONS   ######################################### */
#define ALLOC(type) type*  alloc_ ## type(size_t n) { \
	return (type *) malloc(sizeof(type)*n); \
}
/*!
 *  \def ALLOC(type) 
 *  \brief  macro for  basic allocation function named alloc_(type)
 */
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
#define TAG 1
#define DIMEN  (3)                                /*!< \brief 3d - > 3 */
#define DOF  (2*(1+3+3))                        /*!< \brief 2x(density+velocity+magnet) */
#define FDOF (DIMEN*DOF)                          /*!< \brief memory for correlator */
#define AVDOF  5                                //!< (vel,mag)*(long,trans) + density
int n_FDOF = FDOF;
int n_DOF = DOF;
int n_AVDOF = AVDOF;
enum  { DEN=6, V_X=0,V_Y=1,V_Z=2, M_X=3,M_Y=4,M_Z=5,
	V_LONG=0,V_TRANS=1,M_LONG=2,M_TRANS=3, AV_DEN=4} ;


#define CSet(a,x,y) 																		\
	a.R = x,																							\
	a.I = y
#define BUFF_LEN 1024
/* char *header[]= {"cur-long", "cur-trans", "densty"},
 * 		 *txtCorr = "space-time corr";
 */
void* unused_pointer;
int ununused_value;
long long ll_mem_size=0;
#define AllocMem(a, n, t)                       \
	a = (t *) malloc ((n) * sizeof (t));      \
  ll_mem_size += n* sizeof(t);
#define AllocMem2(a, n1, n2, t)                        \
	AllocMem (a, n1, t *);                               \
AllocMem (a[0], (n1) * (n2), t);                     \
for (ununused_value = 1; ununused_value < n1; ununused_value ++) a[ununused_value] = a[ununused_value - 1] + n2;
/*!
 *  \def AllocMem(a, n,  t)  
 *  @param[out] a     array of type t
 *  @param[in] n    array size
 *  @param      t     type 
 */
/*!
 *  \def AllocMem2(a, n1, n2, t)  
 *  @param[out] a     2 dimensitional array of type t
 *  @param[in] n1,n2   n1 x n2 array size
 *  @param      t     type 
 */
/* #####   EXPORTED FUNCTION DECLARATIONS   ######################################### */
void Init_reciprocal_space(Snapshot * snap);

void ZeroSpacetimeCorr ();
void InitSpacetimeCorr ();
void EvalOtherInformation ();
void PrintSpacetimeCorr (FILE *fp);
void EvalSpacetimeCorr (Snapshot * snap);
void AllocArray();
void Alloc_more();
void AccumSpacetimeCorr ();

/* #####   EXPORTED DATA TYPES   #################################################### */
real kVal, deltaT, rVal, g_Vol;
real L ;                                        /*!< \brief box length */
int nPtls;


TBuf *tBuf;
real **avAcfFcol, *valFcol, **valDqt, **valGammaQT ;
real **avAcfFself, *valFself;

	real **avDrTable;
	real *factorDr;
int countCorrAv, limitCorrAv, nCBuffer, nCSpatial, nCTime;
real *rrMSDAv;
real *rrMQDAv;
real *rrDt;


#endif
