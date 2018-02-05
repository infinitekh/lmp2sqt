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
 *       Created:  2017- 05- 29
 *      Revision:  1
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
#include <omp.h>
/* #####   EXPORTED TYPE DEFINITIONS   ############################################## */
//typedef struct {
//	real R, I;
//} Cmplx;
/*!
 * move to common 
 *  \brief  struct for complex real value
 */
//typedef struct { 
//	real xx,xy,xz;
//	real yx,yy,yz;
//	real zx,zy,zz;
//} Rank2R3;
//typedef struct { 
//	real x,y,z;
//} VecR3;
typedef struct {
	real **F_qq2, *org_rho_q1;
	real **F_s_qq2, **F_d_qq2, **org_rho_s_q1 , **org_rho_d_q1  ;
	VecR3 *orgR, *rTrue; 
	real *rho_q1;
	real **rho_s_q1, **rho_d_q1;
	real *rho_s_q1_temp;
	Rank2R3 orgSumVR;
	Rank2R3 sumVR_ct;
	// VecR3 *orgVel; real *acfVel;
	Rank2R3 *rrMSR2_VR;
	VecR3 *rrMSR1_R;
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
#define N_AXIS  (3+4+6)                                /*!< \brief 3d - > 3 */
#define DOF  (2*(1))                        /*!< \brief 2x(density+velocity+magnet) */
#define FDOF (N_AXIS*DOF)                          /*!< \brief memory for correlator */
#define AVDOF  1                                 //!< (vel,mag)*(long,trans) + density
int n_FDOF = FDOF;
int n_DOF = DOF;
int n_AVDOF = AVDOF;
enum  { DEN=6, V_X=0,V_Y=1,V_Z=2, M_X=3,M_Y=4,M_Z=5,
	V_LONG=0,V_TRANS=1,M_LONG=2,M_TRANS=3, AV_DEN=4} ;
enum {VXC =0, VXS, VYC, VYS, VZC, VZS,
	MXC =6, MXS, MYC, MYS, MZC, MZS, 
	OneC=12, OneS} ; // Those are not used.


#define CSet(a,x,y) 																		\
	a.R = x,																							\
	a.I = y
#define BUFF_LEN 1024
/* char *header[]= {"cur-long", "cur-trans", "densty"},
 * 		 *txtCorr = "space-time corr";
 */
void* unused_pointer;
int ununused_value;
long long int ll_mem_size=0;
int ErrorAllocMem=0; 
#define AllocMem(a, n, t)                       \
	a = (t *) malloc ((n) * sizeof (t));      \
	if (a == NULL ) ErrorAllocMem=1;          \
  ll_mem_size += (long long)n* sizeof(t);
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
/* #####   EXPORTED DATA TYPES   #################################################### */
// local Data
typedef struct {
	TBuf *tBuf;
	int flag_alloc,flag_alloc_more;
	Snapshot* snap;
} MakeSqtClass;
MakeSqtClass* classSqt;
// global data
real kVal, deltaT, rVal, g_Vol,mass;
real L ;                                        /*!< \brief box length */
int nPtls;
/*!
 *  \brief  for Intermediate scattering function <rho(q,t)rho(-q,0)>
 */
real **avF_qq2,  **valDqt, **valGammaQT ;
real **avF_s_qq2, **avF_d_qq2;


real *factorDr, *radius;
int countCorrAv, limitCorrAv, nCBuffer, nCSpatial, nCTime;



/*!
 *  \brief  for van Hove function (globall
 */
real *rrDt;
real **avDrTable;
/*!
 *  \brief  Accumulate and average value ( globall)
 */
real *rrMSR2_VR_Av_offdig;    
real *rrMSR2_VR_Av_dig;
VecR3 *rrMSR1_R_Av;
real *rrMQDAv;
real *rrMSDAv;
Rank2R3 *rrMSR2_VR_Av;


/*!
 *  \brief  openmp locker
 */
int nthreads, threadID;
omp_lock_t write_lock,read_lock;
/* #####   EXPORTED FUNCTION DECLARATIONS   ######################################### */
void Init_reciprocal_space(Snapshot*);

void ZeroAvSpacetimeCorr ();
void InitSpacetimeCorr (MakeSqtClass* );
void EvalOtherInformation ();
void PrintSpacetimeCorr (FILE *fp);
void EvalSpacetimeCorr (MakeSqtClass*);
void AllocArray(MakeSqtClass* );
void AllocMemCheck ();
void Alloc_more(MakeSqtClass* );
void AccumSpacetimeCorr (MakeSqtClass*); // __threadsafe__

#endif
