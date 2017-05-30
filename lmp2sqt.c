/*!
 *    \file  lmp2sqt.c
 *   \brief  
 *
 *  
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

#include "lmp2sqt.h"
#include "snapshot.h"

#include <assert.h>

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include	<unistd.h>
#define M_PI       3.14159265358979323846
#include<math.h>
#include"snapshot.h"


#define DIM 3

typedef enum {N_I, N_R} VType;
typedef struct {
	char *vName;
	void *vPtr;
	VType vType;
	int vLen, vStatus;
} NameList;
#define NameI(x)                      \
{#x, &x, N_I, sizeof (x) / sizeof (int)}
#define NameR(x)                       \
{#x, &x, N_R, sizeof (x) / sizeof (real)}

typedef struct {
	void *vPtr;
	VType vType;
	int vLen;
} ValList;
#define ValI(x)    \
{&x, N_I, sizeof (x) / sizeof (int)}
#define ValR(x)    \
{&x, N_R, sizeof (x) / sizeof (real)}

#define NP_I        \
	((int *) (nameList[k].vPtr) + j)
#define NP_R         \
	((real *) (nameList[k].vPtr) + j)

NameList nameList[] = {
	NameR   (rVal),
	NameR   (kVal),
	NameR   (deltaT),
	NameI  	(limitCorrAv),
	NameI   (nBuffCorr),       // number of simul. time seq
	NameI   (nFunCorr),        // number of spatial seq
	NameI   (nValCorr)         // number of time seq
};

void PrintNameList (FILE *fp);
int GetNameList (int argc, char **argv);

int main(int argc, char** argv) {
	/*!
	 *  \brief  main 함수. 설명이 필요없다. 
	 *
	 *  \param   argc 
	 *  \param   argv 
	 */
	char filename[100];
	int n_snap;
	if(argc <2) {
		perror("#run inputfilename");
		return 1;
	}
	GetNameList(argc,argv);


	AllocArray();
	InitSpacetimeCorr();

	strcpy( filename,argv[1]);
	FILE* fp = fopen( filename ,"r");
	Snapshot* snap;
	n_snap = 0;	

	// kVal value have be changed because reciprocal information
	snap = read_dump(fp);
	Init_reciprocal_space(snap);
	rewind(fp);
	PrintNameList(stdout);

	while(1) {
		snap =	read_dump(fp);

		if (snap == NULL)
			break;

		EvalSpacetimeCorr(snap);

		free_Snapshot(snap);
		n_snap++; 
	}

	if (n_snap <5){
		perror("The # of snap is too small(<5)!!");
		return 23;
	}

	return 0;
}

void AccumSpacetimeCorr ()
	/*!
	 *  \brief  계산된 현재 시간의 SpaceTime correlation을 누적한다. 
	 */
{
	real fac;
	int j,  nb, nr, n;
	for (nb = 0; nb < nBuffCorr; nb ++) {
		if (tBuf[nb].count == nValCorr) {
			// S(q,t), M(q,t) part
			for (j = 0; j < 3 * nFunCorr; j ++) {
				for (n = 0; n < nValCorr; n ++)
					avAcfFcol[j][n] += tBuf[nb].acfFcol[j][n];
			}
			// Diffuse Part
			for (j = 0; j < nValCorr; j ++) {
				rrMSDAv[j] += tBuf[nb].rrMSD[j];
				rrMQDAv[j] += tBuf[nb].rrMQD[j];
				for ( nr=0; nr<nFunCorr; nr++) 
					avDrTable[nr][j] += tBuf[nb].DrTable[nr][j];
			}
			// buffer nb reset
			tBuf[nb].count = 0;
			++ countCorrAv;
			if (countCorrAv == limitCorrAv) {
				for (j = 0; j < 3 * nFunCorr; j ++) {
					for (n = 0; n < nValCorr; n ++)
						avAcfFcol[j][n] /= 3. * nPtls * limitCorrAv;
				}
//				fac = 1./ ( DIM * 2 * nPtls * deltaT * limitCorrAv); 
/*-----------------------------------------------------------------------------
*   rrMSDAv -> mean square displacemnt 
*   rrMQDAv -> mean quadropole displacemnt 
*-----------------------------------------------------------------------------*/
//				fac = 1./ ( DIM * 2 * nPtls * deltaT * limitCorrAv); 
				fac = 1./ ( nPtls *  limitCorrAv); 
				for (j = 1; j < nValCorr; j ++) {
					rrMSDAv[j] *= fac;
					rrMQDAv[j] *= fac;
					for ( nr=0; nr<nFunCorr; nr++) 
						avDrTable[nr][j] *= factorDr[nr];
				}
				PrintSpacetimeCorr (stdout);
				ZeroSpacetimeCorr ();
			}
		}
	}
}


void InitSpacetimeCorr ()
	/*!
	 *  \brief  프로그램 초기에 시간 평균을 낼 수 있도록 index를 부여하는 과정
	 */
{
	int nb;
	if (nBuffCorr > nValCorr) {
		fputs("Error nBuffCorr> nValCorr", stderr);
		exit(1);
	}

	for (nb = 0; nb < nBuffCorr; nb ++){
		tBuf[nb].count = - nb * nValCorr / nBuffCorr;
		tBuf[nb].countDiff = - nb * nValCorr / nBuffCorr;
	}
	ZeroSpacetimeCorr ();
}
void ZeroSpacetimeCorr ()
	/*!
	 *  \brief  출력 후 또는, 프로그램 시작시 평균 계산을 위한 메모리를 
	 *  			0값으로 초기화
	 */
{
	int j, n, nr;
	countCorrAv = 0;
	for (j = 0; j < 3 * nFunCorr; j ++) {
		for (n = 0; n < nValCorr; n ++) avAcfFcol[j][n] = 0.;
	}

	countDiffuseAv = 0;
	for (j = 0; j < nValCorr; j ++) {
		rrMSDAv[j] = 0.;
		rrMQDAv[j] = 0.;
	}
	for (j = 0; j < nValCorr; j ++) 
		for (nr = 0; nr < nFunCorr; nr ++) avDrTable[nr][j]= 0.;
}
void EvalOtherInformation () 
	/*!
	 *  \brief  \f$ F(q,t) \f$를 출력전에 미분해서 data를 뽑아낸다. 
	 *
	 */
{                            // this evaluation yield analysis.c
#define Fqt_FIX_q avAcfFcol[3*(j) +2] 
	int j,  n,  ppT, pT, cT, nT, nnT;
	extern real kVal;
	real kVal2 = kVal*kVal;
	n=0; nnT = n+2; nT = n+1; cT = n; {  //Forward O(h^2)
		for (j = 0; j < nFunCorr; j ++) {
			valGammaQT[j][n]=		 (-(Fqt_FIX_q[nnT]) +4.*(Fqt_FIX_q[nT]) -3.*(Fqt_FIX_q[cT]) )/ (2.0* deltaT*Fqt_FIX_q[cT]);
			valDqt [j][n] = - valGammaQT[j][n] / (kVal2*j*j) ;
		}
	}
	for (n = 1; n < nValCorr-1; n ++) {     /* centerd O(h^2) */
		pT = n-1; nT = n+1; cT = n;
		for (j = 0; j < nFunCorr; j ++) {
			valGammaQT[j][n] = ( (Fqt_FIX_q[nT]) -(Fqt_FIX_q[pT]) )/ (2.0* deltaT*Fqt_FIX_q[cT]);
			valDqt [j][n] = - valGammaQT[j][n] / (kVal2*j*j) ;
		}
	}
	n= nValCorr-1; ppT = n-2; pT = n-1; cT = n; { /* Backward O(h^2) */
		for (j = 0; j < nFunCorr; j ++) {
			valGammaQT[j][n] = (+(3.*Fqt_FIX_q[cT]) -4.*(Fqt_FIX_q[pT]) +(Fqt_FIX_q[ppT]) )/ (2.0* deltaT*Fqt_FIX_q[cT]);
			valDqt [j][n] = - valGammaQT[j][n] / (kVal2*j*j) ;
		}
	}

}


void PrintSpacetimeCorr (FILE *fp)
	/*!
	 *  \brief  결과를 출력하는 함수
	 *
	 *  \param  fp output file descriptor
	 */
{
	extern real kVal;
	real tVal;
	int j, k, n, k2, nr;
//	char *header[] = {"cur-long", "cur-trans", "density", "vanHove-self"};
	char *header[] = {"cur-long", "cur-trans", "density", "vanHove-self","self-long", "self-trans", "self"};
	fprintf (fp, "space-time corr\n");
	//for (k = 0; k < 3; k ++) {
	for (k2 = 0; k2 < sizeof(header)/ sizeof(char*); k2 ++) {
		/* 		fprintf (fp, "%s", header[k]);
		 * 		for (j = 0; j < nFunCorr; j ++)
		 * 			fprintf (fp, " %7.3f", kVal*(j+1));
		 * 		fprintf (fp, "\n");
		 */

		//    EvalOtherInformation ();
		fprintf (fp, "# %s %7.3f %7.3f %7.3f\n", header[k2] , kVal, 1.0*deltaT, rVal);
		switch ( k2) {
			case 0: case 1: case 2: 
				/*-----------------------------------------------------------------------------
				 *  avAcfFcol[3*i+k][j] -> F(q_i,t_j) k=0 longi k=1 tranv k=2 density
				 *-----------------------------------------------------------------------------*/
				k= k2%3;
				for (n = 0; n < nValCorr; n ++) {
					/* 			deltaT = n *1. * deltaT;
					 * 			fprintf (fp, "%7.3f", deltaT);
					 */
					for (j = 0; j < nFunCorr; j ++){
						fprintf (fp, " %8.4e", avAcfFcol[3 * j + k][n]);
					}
					fprintf (fp, "\n");
				} 
				break;
			case 3:                /* gamma_qt */
				//        fprintf (fp, "#van Hove function\n");

				for (j = 0; j < nValCorr; j ++) {
					for ( nr=0; nr<nFunCorr; nr++)  {
						fprintf (fp, " %8.4e", avDrTable[nr][j] );
					}
					fprintf (fp, "\n");
				}
				break;
			case 4: case 5: case 6: 
				/*-----------------------------------------------------------------------------
				 *  avAcfFself[3*i+k][j] -> F(q_i,t_j) k=0 longi k=1 tranv k=2 density
				 *-----------------------------------------------------------------------------*/

				k= (k2-4)%3;
				for (n = 0; n < nValCorr; n ++) {
					/* 			deltaT = n *1. * deltaT;
					 * 			fprintf (fp, "%7.3f", deltaT);
					 */
					for (j = 0; j < nFunCorr; j ++){
						fprintf (fp, " %8.4e", avAcfFself[3 * j + k][n]);
					}
					fprintf (fp, "\n");
				} 
				break;
		}
		fprintf (fp, "\n");
	}

	//  char filename1[100] ="Dq00.info" ;
	//  char filename2[100] ="Ft00.info" ;
	char filename1[100];
	char filename2[100];
	char filename3[100];
	int nfile = 0;
		sprintf(filename1, "Dt%02d.info",nfile);
		sprintf(filename2, "vanHove%02d.info",nfile);
		sprintf(filename3, "SSF%02d.info",nfile);
		//printf( "access(%s) -> return %d", filename1, access(filename1,F_OK));
	while( 0 == access(filename1,F_OK) ) {
/* 		fprintf(stderr, "Files are  exist at least . (%02d) \n", nfile);
 * 		sleep(1);
 */
		nfile++;
		sprintf(filename1, "Dt%02d.info",nfile);
		sprintf(filename2, "vanHove%02d.info",nfile);
		sprintf(filename3, "SSF%02d.info",nfile);
	}

	/* 	FILE* fp_Dq = fopen(filename1,"w");
	 * 	fprintf (fp_Dq, "# dt = %7.3f\n", deltaT);
	 * 	for (j = 0; j < nFunCorr; j ++) {
	 * 		fprintf (fp_Dq, "%8.4f" , j*kVal );
	 * 		for (n = 1; n < nValCorr; n ++) {   
	 * 			fprintf (fp_Dq, " %8.4e" ,  valDqt[j][n]);
	 * 		}
	 * 		fprintf (fp_Dq, "\n");
	 * 	}
	 * 	fclose(fp_Dq);
	 */

	/* 	FILE* fp_Ft = fopen(filename2,"w");
	 * 	fprintf (fp_Ft, "# dq = %7.3e\n", kVal);
	 * 	for (n = 0; n < nValCorr; n ++) {   
	 * 		fprintf (fp_Ft, "%8.4f" , n*deltaT );
	 * 		for (j = 0; j < nFunCorr; j ++) {
	 * 			fprintf (fp_Ft, " %8.4e" ,  avAcfFcol[(3*j)+2][n]/avAcfFcol[(3*j)+2][0]);
	 * 		}
	 * 		fprintf (fp_Ft, "\n");
	 * 	}
	 * 	fclose(fp_Ft);
	 */

	FILE* fp_SSF = fopen(filename3,"w");
	n=0;
	for (j = 0; j < nFunCorr; j ++) {
		fprintf (fp_SSF, "%8.4f" " %8.4e""\n" , (j+1)*kVal , avAcfFcol[(3*j)+2][0]);
	}
	fclose(fp_SSF);

	//  fprintf (fp_SSF, "# dq = %7.3e\n", kVal);




	FILE* fp_Dt = fopen(filename1,"w");
	fprintf (fp_Dt, "#time MSD diffusion\n");

	real fac = 1./( 2.* deltaT * DIM * 2);
	j=0;
	rrDt[j] = fac*(-rrMSDAv[j+2]  +4.*rrMSDAv[j+1] - 3.* rrMSDAv[j]);
	for ( j = 1; j < nValCorr-1; j += 1 ) {
		rrDt[j] = fac*(rrMSDAv[j+1]  -rrMSDAv[j-1] );
	}
	j=nValCorr-1;
	rrDt[j] = fac*(rrMSDAv[j-2]  -4.*rrMSDAv[j-1] + 3.* rrMSDAv[j]);


	for ( j = 0; j < nValCorr; j += 1 ) {
		tVal = j * deltaT;
		fprintf (fp_Dt, "%8.4f %8.4e %8.4e %8.4e\n", tVal, rrMSDAv[j] , rrDt[j], rrMQDAv[j]); 
	}
	fclose(fp_Dt); 

	/*-----------------------------------------------------------------------------
	 *  van Hove function part
	 *-----------------------------------------------------------------------------*/
	/* 	fprintf (fp_Gr, "#van Hove function\n");
	 * 	FILE* fp_Gr = fopen(filename2,"w");
	 * 	
	 * 	for ( nr=0; nr<nFunCorr; nr++)  {
	 * 		for (j = 0; j < nValCorr; j ++) {
	 * 			fprintf (fp_Gr, " %8.4e", avDrTable[nr][j] );
	 * 		}
	 * 		fprintf (fp_Gr, "\n");
	 * 	}
	 * 	fprintf (fp_Gr, "\n");
	 * 	fclose(fp_Gr);
	 */

}

void EvalSpacetimeCorr(Snapshot* snap)
	/*!
	 *  \brief  space time correlation을 계산한다. 
	 *  				q space value는 x, y, z 방향 세개의 방향으로 
	 *  				longitudinal version, translational version과 가장 기본적인 방향성분 없는 density
	 *  				version 3개를 구함.
	 *						기계 편의적으로 코드가 짜져 있어서 사람이 보기에 별로 직관적이지 못해서 고칠 예정이고 
	 *			거기다가 참조한 기본 코드에서 magnetization에 대한 version으로 바꾸면서 많이 복잡해지고 좋지 
	 *			않아짐. 그리고 개인적으로 속도 vector도 다시 정보로 가져올 것이기 때문에 바뀔 야정 
	 *
	 * 			PREV. $M_T(q,t)$ $M_L(q,t)$ 
	 *			Todo. 
	 *  \param  Snapshot* Snapshot 포인터 
	 */
{
	real b, c, c0, c1, c2, s, s1, s2, w;
	extern real kVal;
	int j, k, m, n, nb, nc, ni, nv, nr;
	real r[3], mu[3];          // L[3];

	VecR3 dr;
	real deltaR2;
	int i_Dr;
	real L = snap->box.xhigh- snap->box.xlow;
	g_Vol  = L*L*L;
	nPtls = snap->n_atoms;
	static int first_run = 0;
	if (first_run ==0 ) 
		Alloc_more();
	first_run = 1;
	atom* col_i;
	// zero initalize current time value
	for (j = 0; j < 24 * nFunCorr; j ++) valFcol[j] = 0.;
	// we assume L0=L1 = L2 
	/* 	L[0] = snap->box.xhigh- snap->box.xlow;
	 * 	L[1] = snap->box.yhigh- snap->box.ylow;
	 * 	L[2] = snap->box.zhigh- snap->box.zlow;
	 */
	kVal = 2.*M_PI / L;

	//  kVal = 2. * M_PI / nFunCorr;
	//  Begin Calculatie current time information
	/*-----------------------------------------------------------------------------
	 *  Direct calculate  rho(q)
	 *-----------------------------------------------------------------------------*/
	for (n=0; n<nPtls; n++) {
		col_i = &(snap->atoms[n]);
		/* 		r[0] = col_i->x; r[0] = r[0] - L* floor(r[0]/L)- L/2.;  
		 * 		r[1] = col_i->y; r[1] = r[1] - L* floor(r[1]/L)- L/2.;  
		 * 		r[2] = col_i->z; r[2] = r[2] - L* floor(r[2]/L)- L/2.;  
		 */
		r[0] = col_i->x; 
		r[1] = col_i->y; 
		r[2] = col_i->z; 
		mu[0] = col_i->mux; 
		mu[1] = col_i->muy; 
		mu[2] = col_i->muz;
		j = 0;
		// 
		for (k = 0; k < 3; k ++) {
			for (m = 0; m < nFunCorr; m ++) {
				if (m == 0) {
					b = kVal * r[k];
					c = cos (b);
					s = sin (b);
					c0 = c;
				} else if (m == 1) {
					c1 = c;
					s1 = s;
					c = 2. * c0 * c1 - 1.;
					s = 2. * c0 * s1;
				} else {
					c2 = c1;
					s2 = s1;
					c1 = c;
					s1 = s;
					c = 2. * c0 * c1 - c2;
					s = 2. * c0 * s1 - s2;
				}
				valFcol[j ++] += mu[0] * c;
				valFcol[j ++] += mu[0] * s;
				valFcol[j ++] += mu[1] * c;
				valFcol[j ++] += mu[1] * s;
				valFcol[j ++] += mu[2] * c;
				valFcol[j ++] += mu[2] * s;
				valFcol[j ++] += c;
				valFcol[j ++] += s;
			}
		}
		if (n == 0 ) {
			for (j = 0; j < 24 * nFunCorr; j ++) valFself[j] = valFcol[j];
//			memcpy(valFself, valFcol,sizeof(real)*24*nFunCorr);
		}
  }                                             /* for loop : n<nPtls */
	// End Calculate Current time value
	// Begin Two time corrlation function
	for (nb = 0; nb < nBuffCorr; nb ++) {
		if (tBuf[nb].count == 0) {
			/*-----------------------------------------------------------------------------
			 *   t_w information 
			 *-----------------------------------------------------------------------------*/
			for (n=0; n<nPtls; n++) {
				col_i = &(snap->atoms[n]);
				tBuf[nb].orgR[n].x = col_i->x;
				tBuf[nb].orgR[n].y = col_i->y;
				tBuf[nb].orgR[n].z = col_i->z;
			}
//			memcpy(tBuf[nb].orgFcol ,  valFcol,sizeof(real)*24*nFunCorr);
//			memcpy(tBuf[nb].orgFself, valFself,sizeof(real)*24*nFunCorr);
			for (j = 0; j < 24 * nFunCorr; j ++){
				tBuf[nb].orgFcol[j] = valFcol[j];
				tBuf[nb].orgFself[j] = valFself[j];
			}
		}                        // End   buffer count ==0

		if (tBuf[nb].count >= 0) {
			/*-----------------------------------------------------------------------------
			 *  Get Delta_r(t)
			 *-----------------------------------------------------------------------------*/
			ni = tBuf[nb].count;
			/*-----------------------------------------------------------------------------
			 *  Zero initializing
			 *-----------------------------------------------------------------------------*/
			tBuf[nb].rrMSD[ni]= 0.;
			tBuf[nb].rrMQD[ni]= 0.;
			for ( nr=0; nr<nFunCorr; nr++) 
				tBuf[nb].DrTable[nr][ni] =0;

			/*-----------------------------------------------------------------------------
			 *  Calculation at this time
			 *-----------------------------------------------------------------------------*/

			for (n=0; n<nPtls; n++) {
				col_i = &(snap->atoms[n]);
				dr.x =  col_i->x-tBuf[nb].orgR[n].x ;
				dr.y =  col_i->y-tBuf[nb].orgR[n].y ;
				dr.z =  col_i->z-tBuf[nb].orgR[n].z ;
				deltaR2 = (dr.x*dr.x+dr.y*dr.y+dr.z*dr.z);

				tBuf[nb].rrMSD[ni] += deltaR2;
				tBuf[nb].rrMQD[ni] += deltaR2*deltaR2;

				i_Dr    = (int) (sqrt(deltaR2)/rVal);
				if (i_Dr<nFunCorr)
					tBuf[nb].DrTable[i_Dr][ni] ++;
			}

			/*-----------------------------------------------------------------------------
			 *  two time correlation 
			 *-----------------------------------------------------------------------------*/
			//acfFcol 0 initialization
			for (j = 0; j < 3 * nFunCorr; j ++) {
				tBuf[nb].acfFcol[j][ni] = 0.;
				tBuf[nb].acfFself[j][ni] = 0.;
			}
			// add AcfFcol
			for (j = 0,k = 0; k < 3; k ++) { // 3 loop
				for (m = 0; m < nFunCorr; m ++) {
					for (nc = 0; nc < 4; nc ++) {  
						nv = 3 * m + 2;
						if (nc < 3) {    /*  */
							//              w = Sqr (kVal * (m + 1));
							w=1.0;
							/*-----------------------------------------------------------------------------
							 *  nv %3 = >> 0 - longitunidal, 1 - transverse, 2 - density 
							 *-----------------------------------------------------------------------------*/
							-- nv;         // transverse  3*m +1
							if (nc == k) -- nv; // longitudinal 3*m
							else w = 0.5;
							//              else w *= 0.5;
						} else w = 1.;   // density   3*m+2
						// cos(q*r(t)) cos(q*r(t_w) +sin sin
						tBuf[nb].acfFself[nv][ni] +=
							w * (valFself[j] * tBuf[nb].orgFself[j] +
									valFself[j + 1] * tBuf[nb].orgFself[j + 1]);
						tBuf[nb].acfFcol[nv][ni] +=
							w * (valFcol[j] * tBuf[nb].orgFcol[j] +
									valFcol[j + 1] * tBuf[nb].orgFcol[j + 1]);
						j += 2;
					}  // for nc, total j+=8
					assert ( j%8 ==0);
				}    // for m , total j+=8*nFunCorr
				assert (j%(8*nFunCorr)==0);
			}      // for k , total j+=3*8*nFunCorr
			assert( j== 3 * 8 *nFunCorr);
		}                        // End buffer count >=0
		++ tBuf[nb].count;
	}
	AccumSpacetimeCorr (nPtls );
}

void AllocArray ()
	/*!
	 *  \brief   이름그대로 memory 할다함. 
	 *       		대상은 valFself   self  intermediate scattering function과
	 *     				valFcol     collective intermediate scatterring function 
	 *     	등.
	 */
{
	int nb;
	AllocMem (valFself, 3 * 8 * nFunCorr, real);
	AllocMem2 (avAcfFself, 3 * nFunCorr, nValCorr, real);

	AllocMem (valFcol, 3 * 8 * nFunCorr, real);
	AllocMem2 (avAcfFcol, 3 * nFunCorr, nValCorr, real);

	AllocMem2 (valDqt,  nFunCorr, nValCorr, real);
	AllocMem2 (valGammaQT,  nFunCorr, nValCorr, real);
	AllocMem (tBuf, nBuffCorr, TBuf);
	for (nb = 0; nb < nBuffCorr; nb ++) {
		AllocMem (tBuf[nb].orgFself, 24 * nFunCorr, real);
		AllocMem2 (tBuf[nb].acfFself, 3 * nFunCorr, nValCorr, real);

		AllocMem (tBuf[nb].orgFcol, 24 * nFunCorr, real);
		AllocMem2 (tBuf[nb].acfFcol, 3 * nFunCorr, nValCorr, real);
	}
	// AllocArray for Diffuse ()
	AllocMem (rrMSDAv, nValCorr, real);
	AllocMem (rrDt, nValCorr, real);
	AllocMem (rrMQDAv, nValCorr, real);
	AllocMem2 (avDrTable, nFunCorr,nValCorr, real);
}
void Alloc_more () {
	/*!
	 *  \brief  Alloc_more 
	 *
	 */
	int nb,nr; real rho0, shell_Vol;
	for (nb = 0; nb < nBuffCorr; nb ++) {
		AllocMem (tBuf[nb].orgR, nPtls, VecR3);
		AllocMem (tBuf[nb].rrMSD, nValCorr, real);
		AllocMem (tBuf[nb].rrMQD, nValCorr, real);
		AllocMem2 (tBuf[nb].DrTable, nFunCorr,nValCorr, int);
	}
	AllocMem (factorDr, nFunCorr, real);

	rho0 = nPtls/g_Vol;
	for (nr = 0; nr < nFunCorr; nr ++) {
		if (nr ==0) shell_Vol = 4*M_PI /3. * pow(rVal,3);
		else shell_Vol = 4*M_PI * pow(rVal,3)* (nr*nr + 1./12.);

		factorDr[nr] = 1./( pow(rho0,2) * g_Vol *shell_Vol*limitCorrAv);
		/* 		printf("rho0=%.2e, Vol=%.2e, shell_Vol=%.2e, factorDr=%.2e\n", 
		 * 				rho0,g_Vol,shell_Vol,factorDr[nr]);
		 */

	}

}


int GetNameList (int argc, char **argv)
	/* 
	 * ===  FUNCTION  ======================================================================
	 *         Name:  GetNameList
	 *  Description:  Rapaport 책에서 가져온 코드로, 
	 *  						Name  value 
	 *  						Name2 value 
	 *  				형식으로 되어 있는 초기값을 불러오는데 사용됨. 
	 *  				지원하는 형태는 integer와 real 값을 받아다가 초기값 배정함. 
	 *  				잘 작동함. 
	 * =====================================================================================
	 */
{
	int  j, k, match, ok;
	char buff[80], *token;
	FILE *fp;
	strcpy (buff, argv[0]);
	strcat (buff, ".in");
	if ((fp = fopen (buff, "r")) == 0)  {
		fp = fopen(buff, "w");
		for (k = 0; k < sizeof (nameList) / sizeof (NameList); k ++) {
			fprintf (fp, "%s\t", nameList[k].vName);
			if (strlen (nameList[k].vName) < 8) fprintf (fp, "\t");
			for (j = 0; j < nameList[k].vLen; j ++) {
				switch (nameList[k].vType) {
					case N_I:
						fprintf (fp, "%d ", *NP_I);
						break;
					case N_R:
						fprintf (fp, "%#g ", *NP_R);
						break;
				}
				fprintf (fp, "\n");
			}
		}
		fprintf (fp, "----\n");
		fclose(fp);
		exit (1);
	}

	for (k = 0; k < sizeof (nameList) / sizeof (NameList); k ++)
		nameList[k].vStatus = 0;
	ok = 1;
	while (1) {
		fgets (buff, 80, fp);
		if (feof (fp)) break;
		token = strtok (buff, " \t\n");
		if (! token) break;
		match = 0;
		for (k = 0; k < sizeof (nameList) / sizeof (NameList); k ++) {
			if (strcmp (token, nameList[k].vName) == 0) {
				match = 1;
				if (nameList[k].vStatus == 0) {
					nameList[k].vStatus = 1;
					for (j = 0; j < nameList[k].vLen; j ++) {
						token = strtok (NULL, ", \t\n");
						if (token) {
							switch (nameList[k].vType) {
								case N_I:
									*NP_I = atol (token);
									break;
								case N_R:
									*NP_R = atof (token);
									break;
							}
						} else {
							nameList[k].vStatus = 2;
							ok = 0;
						}
					}
					token = strtok (NULL, ", \t\n");
					if (token) {
						nameList[k].vStatus = 3;
						ok = 0;
					}
					break;
				} else {
					nameList[k].vStatus = 4;
					ok = 0;
				}
			}
		}
		if (! match) ok = 0;
	}
	fclose (fp);

	if(nBuffCorr > nValCorr ) nBuffCorr = nValCorr;

	for (k = 0; k < sizeof (nameList) / sizeof (NameList); k ++) {
		if (nameList[k].vStatus != 1) ok = 0;
	}
	return (ok);
}
void PrintNameList (FILE *fp)
	/*!
	 *  \brief 초기값을 출력하는 함수 getNameList의 짝함수이다.   
	 *
	 *  \param  fp  FILE* file descriptor
	 */
{
	int j, k;
	fprintf (fp, "NameList -- data\n");
	for (k = 0; k < sizeof (nameList) / sizeof (NameList); k ++) {
		fprintf (fp, "%s\t", nameList[k].vName);
		if (strlen (nameList[k].vName) < 8) fprintf (fp, "\t");
		if (nameList[k].vStatus > 0) {
			for (j = 0; j < nameList[k].vLen; j ++) {
				switch (nameList[k].vType) {
					case N_I:
						fprintf (fp, "%d ", *NP_I);
						break;
					case N_R:
						fprintf (fp, "%#g ", *NP_R);
						break;
				}
			}
		}
		switch (nameList[k].vStatus) {
			case 0:
				fprintf (fp, "** no data");
				break;
			case 1:
				break;
			case 2:
				fprintf (fp, "** missing data");
				break;
			case 3:
				fprintf (fp, "** extra data");
				break;
			case 4:
				fprintf (fp, "** multiply defined");
				break;
		}
		fprintf (fp, "\n");
	}
	fprintf (fp, "----\n");
}

void Init_reciprocal_space(Snapshot * snap) {
	/*!
	 *  \brief  뭔가 실패한 시도. 관심 가질 필요 없음. 
	 *           원래는 delta_k 값이 2pi/L보다 작을시에 simulation 박스를 복사하여 
	 *           replica를 포함하는 공간에 대해서 계산하는 식으로 하려고 하였으나 
	 *           복소평면 상에서 그려보았을 때 전부 더하면 0이 되는 관계로 
	 *           더이상 구할 수 없고 아무 의미 없음이 확인하었다. 
	 *           periodic boundary라 아닐 때는 상관없이 좋은 결과를 낼 수 있을 것이다. 
	 *           고로 원래 목적과 달리 delta_k를 2pi/L * n(정수)로 맞추도록 한다.
	 *  \param  snap Snaptshot* 스냅샷 포인터
	 */
	extern real kVal;
	int n_mul;
	real L[3];
	// zero initalize current time value
	// we assume L0=L1 = L2 
	L[0] = snap->box.xhigh- snap->box.xlow;
	L[1] = snap->box.yhigh- snap->box.ylow;
	L[2] = snap->box.zhigh- snap->box.zlow;

	/* 	for (k = 0; k < sizeof (nameList) / sizeof (NameList); k ++) {
	 * 		if ( strcmp(vName, nameList[k].vName)== 0 )  {
	 * 			j=0;
	 * 			p_kVal = NP_R;
	 * 		}
	 * 	}
	 * 	printf( "kVal %p kValp %p\n", &kVal, p_kVal);
	 */

	n_mul = kVal/ (2.*M_PI/ L[0] );
	if (n_mul <=0) n_mul =1;
	kVal = (2.*M_PI/L[0]) * n_mul;
}
