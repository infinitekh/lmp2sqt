/*!
 *    \file  lmp2sqt.c
 *   \brief  
 *
 *  \author  KIM Hyeok (kh), ekh0324@gmail.com
 *
 *  \internal
 *       Created:  2017- 05- 29
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
#include <stdbool.h>
#define M_PI       3.14159265358979323846
#include<math.h>
#include"snapshot.h"


#define DIM 3
#define flagSelf 1

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
	/*!
	 *  \brief  system information 
	 */
	NameR   (kVal),
	NameR   (deltaT),
	NameR   (mass),
	/*!
	 *  \brief input parameter for evaluation 
	 */
	NameR   (rVal),
	NameI  	(limitCorrAv),
	NameI   (nCBuffer),       // number of simul. time seq
	NameI   (nCSpatial),        // number of spatial seq
	NameI   (nCTime)         // number of time seq
};
typedef real coordi3[3];
#define OOS2 (0.7071067811865475)
#define OOS3 (0.5773502691896258)
coordi3 coordi_list[13]= {
	{ +1,+0,+0},
	{ +0,+1,+0},
	{ +0,+0,+1}, // basic
	{ +OOS2,+OOS2,+0},
	{ +OOS2,-OOS2,+0},
	{ +OOS2,+0,+OOS2},
	{ +OOS2,+0,-OOS2},
	{ +0,+OOS2,+OOS2},
	{ +0,+OOS2,-OOS2}, // 2axis combination
	{ +OOS3,+OOS3,+OOS3},
	{ +OOS3,+OOS3,-OOS3},
	{ +OOS3,-OOS3,+OOS3},
	{ +OOS3,-OOS3,-OOS3} // 3axis combination
};


char inputFilename[100]= "in.lmp2sqt";
void UpdateNameList ();
void PrintNameList2File (FILE *fp);
int GetNameList (int argc, char **argv);
int Number_call_Print =0;
int flag_global_alloc =0 ;
int flag_global_alloc_more =0 ;
void PrintSpacetimeCorr_binary ( FILE*);

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
		perror("#run inputfilename inputfilename2 ifile3 ifile4...");
		return 1;
	}
	GetNameList(argc,argv);
	int opt_num=1;
	int full_n_snaps=0, full_n_snaps_index=0,progress=-1;
	bool files_on [argc+3];
	long int num_data [argc+3];
	limitCorrAv = 0;

	omp_init_lock(&write_lock);
	omp_init_lock(&read_lock);
	

	for( opt_num = 1;   opt_num < argc; opt_num++)  {
		strcpy( filename,argv[opt_num]);
		FILE* fp = fopen( filename ,"r");
		if( fp == NULL) {
			files_on [opt_num] = false;
			fprintf(stderr,"Can`t open file (%s)!!\n", filename);
			continue;
		}

		// kVal value have be changed because reciprocal information

		n_snap = 0;	
		while(1) {
			bool check =	read_dump_OnlyCheck(fp);

			if (check == false )
				break;

			n_snap++; 
		}
		// if ( flag_Max_eval)
		if (n_snap <5){
			fprintf(stderr,"The # of snap is too small(<5)\n"
					"We would  not use 	this file(%s)!!\n", filename);
			files_on [opt_num] = false;
			//			return 23;
		}
		else if ( floor((n_snap - nCTime) /(nCTime/ nCBuffer)) <1) {
			fprintf(stderr,"The # of snap is too small()\n"
					     "it dont make a infomation\n"
					"We would  not use 	this file(%s)!!\n", filename);
			files_on [opt_num] = false;
		}
		else {
			num_data [opt_num ] = n_snap;
			if (nCBuffer < 1) {
				nCBuffer =1;
			}
			
			limitCorrAv += floor((n_snap - nCTime) /(nCTime/ nCBuffer));
			full_n_snaps += n_snap;
			files_on [opt_num] = true;
		}
		fclose (fp);
	}

	PrintNameList2File(stderr);
	UpdateNameList ();
	int num_files =0;
	for( opt_num = 1;   opt_num < argc; opt_num++)  {
		if (files_on[opt_num] == true) num_files ++;
	}
	nthreads = omp_get_max_threads();
	nthreads = (nthreads< num_files)? nthreads:num_files;
	nthreads = 1;
	
	real * r_done_works;
	int * i_done_works;
	AllocMem(classSqt,nthreads,MakeSqtClass);
	AllocMem(r_done_works,nthreads,real);
	AllocMem(i_done_works,nthreads,int);
	
	for ( int tID =0; tID < nthreads; tID++) {
		classSqt[tID].flag_alloc = 0;
		classSqt[tID].flag_alloc_more = 0;
		i_done_works[tID]= 0;
	}


//#pragma omp parallel for    schedule(dynamic)
	for( opt_num = 1;   opt_num < argc; opt_num++)  {
		if (files_on[opt_num] == false) continue;
		strcpy( filename,argv[opt_num]);
		FILE* fp = fopen( filename ,"r");
		Snapshot* snap, *firstSnap;
		int threadID = omp_get_thread_num();
		MakeSqtClass* cl_sqt = & classSqt[threadID];
		InitSpacetimeCorr(cl_sqt);
		/*!
		 *  \brief  Start Calculation.
		 */
		omp_set_lock(&read_lock);
		rewind(fp);  
			firstSnap = read_dump(fp);
		Init_reciprocal_space(firstSnap);
		rewind(fp);
		omp_unset_lock(&read_lock);
		fprintf(stderr,"FULL_SNAPS = %5d\n",full_n_snaps);
		for (int ns = 0 ; ns < num_data[opt_num] ; ns++ ) {
			/* 		while(1) {}
			*/
			omp_set_lock(&read_lock);
			snap =	read_dump(fp);
			omp_unset_lock(&read_lock);
			if (snap == NULL)
				break;
			cl_sqt->snap = snap;
			EvalSpacetimeCorr(cl_sqt);

			free_Snapshot(snap);
#pragma omp atomic
			full_n_snaps_index++;

			i_done_works[threadID]++;

			int new_progress =(1000.0*full_n_snaps_index/ full_n_snaps);
			r_done_works[threadID] =
				(1000.0*i_done_works[threadID]/ full_n_snaps);
			/* 			fprintf(stderr,"FULL_SNAPS = %5d, newprogress %d\n",full_n_snaps_index,
			 * 					new_progress);
			 */
			if (new_progress != progress ) {
				progress = new_progress;
				fprintf(stderr, "\r %4.1f%%(", progress*.1);
				for (int k =0 ; k<nthreads; k++) 
					fprintf(stderr, "-tid:%d:%4.1f%%", k,r_done_works[k]*.1);
				fprintf(stderr, ")");
				//				fflush(stderr);
			}
		}
		//		opt_num ++;
		fclose (fp);
		fprintf(stderr, "\nEnd : file : %s\n", filename);
	}


	if ( Number_call_Print ==0 ) {
		fprintf(stderr, "limit corr  = %d, countCorrAv = %d\n",
				limitCorrAv,countCorrAv);
		FILE* output = fopen("full_results.info", "w");
		limitCorrAv =countCorrAv;
		PrintSpacetimeCorr(output);
		fclose(output);

	}

	omp_destroy_lock(&write_lock);
	omp_destroy_lock(&read_lock);
	return 0;
}

void AccumSpacetimeCorr (MakeSqtClass* cl_sqt) // __thread_safe__
	/*!
	 *  \brief  계산된 현재 시간의 SpaceTime correlation을 누적한다. 
	 */
{
	int k,  nb, nr, n, nt;
	TBuf* tBuf = cl_sqt->tBuf;
	TBuf* pt;
	for (nb = 0; nb < nCBuffer; nb ++) {
		if (tBuf[nb].count == nCTime) {
			omp_set_lock(&write_lock);
			// check!! that  data is full
			// S(q,t), M(q,t) part
			pt = &tBuf[nb];
			for (k = 0; k < AVDOF * nCSpatial; k ++) {
				for (n = 0; n < nCTime; n ++) {
					avF_qq2[k][n] += pt->F_qq2[k][n];
					avF_s_qq2[k][n] += pt->F_s_qq2[k][n];
					avF_d_qq2[k][n] += pt->F_d_qq2[k][n];
				}
			}
			// Diffuse Part
			for (nt = 0; nt < nCTime; nt ++) {
				rrMSDAv[nt] += pt->rrMSD[nt];
				rrMQDAv[nt] += pt->rrMQD[nt];
			  real_tensor_increase_r1_r1(&rrMSR1_R_Av[nt], 
					 &pt->rrMSR1_R[nt]);
			  real_tensor_increase_r2_r2(&rrMSR2_VR_Av[nt], 
					 &pt->rrMSR2_VR[nt]);
				for ( nr=0; nr<nCSpatial; nr++) {
					avDrTable[nr][nt] += pt->DrTable[nr][nt];
				}
			}
			// buffer nb reset
			pt->count = 0;
			++ countCorrAv;
			if (countCorrAv == limitCorrAv) {
				FILE* fout = fopen( "output.binary", "w+");
				FILE* fout_text = fopen( "output.txt", "w+");
				PrintSpacetimeCorr (fout_text);
				PrintSpacetimeCorr_binary( fout);
				fclose(fout);
				fclose(fout_text);
				ZeroAvSpacetimeCorr ();  //FIXIT   아무의미없음 출력에서 교체합시다. 
			}
			omp_unset_lock(&write_lock);
		} //   if tBuf[nb].count is full
	}   // for all buffer
}

void InitSpacetimeCorr (MakeSqtClass* cl_sqt)
	/*!
	 *  \brief  프로그램 초기에 시간 평균을 낼 수 있도록 index를 부여하는 과정
	 */
{
	if (cl_sqt->flag_alloc == 0 ) {
		omp_set_lock(&write_lock);
		AllocArray(cl_sqt);
		cl_sqt->flag_alloc = 1;
		flag_global_alloc =1;
		AllocMemCheck();
		omp_unset_lock(&write_lock);
	}
	if (nCBuffer > nCTime) {
		fputs("Error nCBuffer> nCTime\n", stderr);
		exit(1);
	}
	TBuf* tBuf = cl_sqt->tBuf;

	for (int nb = 0; nb < nCBuffer; nb ++){
		tBuf[nb].count = - nb * nCTime / nCBuffer;
		tBuf[nb].countDiff = - nb * nCTime / nCBuffer;
	}
}

void ZeroAvSpacetimeCorr ()
	/*!
	 *  \brief  출력 후 또는, 프로그램 시작시 평균 계산을 위한 메모리를 
	 *  			0값으로 초기화
	 */
{
	int nk,nt,  nr;
	countCorrAv = 0;
	for (nk = 0; nk < AVDOF * nCSpatial; nk ++) {
		for (nt = 0; nt < nCTime; nt ++) {
			avF_qq2[nk][nt] = 0.;
			avF_s_qq2[nk][nt] = 0.;
			avF_d_qq2[nk][nt] = 0.;
		}
	}
	for (nt = 0; nt < nCTime; nt ++) {
		rrMSDAv[nt] = 0.; rrMQDAv[nt] = 0.;
		real_tensor_zero_r1( &rrMSR1_R_Av[nt] );
		real_tensor_zero_r2(	&rrMSR2_VR_Av [nt] ) ;
		rrMSR2_VR_Av_offdig[nt] = 0.;
		rrMSR2_VR_Av_dig [nt]   = 0.;
	}
	for (nt = 0; nt < nCTime; nt ++) {
		for (nr = 0; nr < nCSpatial; nr ++)  {
			avDrTable[nr][nt]= 0.;
		}
	}
}

void EvalOtherInformation () 
	/*!
	 *  \brief  \f$ F(q,t) \f$를 출력전에 미분해서 data를 뽑아낸다. 
	 *
	 */
{                            // this evaluation yield analysis.c
#define Fqt_FIX_q avF_qq2[AVDOF*(nk) +AV_DEN] 
	int nk,  n,  ppT, pT, cT, nT, nnT;
	real kVal2 = kVal*kVal, q_sqr;
	n=0; nnT = n+2; nT = n+1; cT = n; 
	{  //Forward O(h^2)
		for (nk = 0; nk < nCSpatial; nk ++) {
			valGammaQT[nk][n]=		 (-(Fqt_FIX_q[nnT]) +4.*(Fqt_FIX_q[nT]) -3.*(Fqt_FIX_q[cT]) )/ (2.0* deltaT*Fqt_FIX_q[cT]);
			q_sqr = kVal2 * (nk+1)*(nk+1);
			valDqt [nk][n] = - valGammaQT[nk][n] / q_sqr ;
		}
	}
	for (n = 1; n < nCTime-1; n ++) {     /* centerd O(h^2) */
		pT = n-1; nT = n+1; cT = n;
		for (nk = 0; nk < nCSpatial; nk ++) {
			valGammaQT[nk][n] = ( (Fqt_FIX_q[nT]) -(Fqt_FIX_q[pT]) )/ (2.0* deltaT*Fqt_FIX_q[cT]);
			q_sqr = kVal2 * (nk+1)*(nk+1);
			valDqt [nk][n] = - valGammaQT[nk][n] / q_sqr ;
		}
	}
	n= nCTime-1; ppT = n-2; pT = n-1; cT = n; { /* Backward O(h^2) */
		for (nk = 0; nk < nCSpatial; nk ++) {
			valGammaQT[nk][n] = (+(3.*Fqt_FIX_q[cT]) -4.*(Fqt_FIX_q[pT]) +(Fqt_FIX_q[ppT]) )/ (2.0* deltaT*Fqt_FIX_q[cT]);
			q_sqr = kVal2 * (nk+1)*(nk+1);
			valDqt [nk][n] = - valGammaQT[nk][n] / q_sqr ;
		}
	}
}

void prePrintProcess () 
{
	real scale_factor = 1./(13.*nPtls*countCorrAv);
#pragma omp parallel for
	for (int nr = 0; nr < AVDOF * nCSpatial; nr ++) {
		for (int nt = 0; nt < nCTime; nt ++){
			avF_qq2[nr][nt] *= scale_factor;
			avF_s_qq2[nr][nt] *= scale_factor;
			avF_d_qq2[nr][nt] *= 0.5*scale_factor;
		}
	}
	//				fac = 1./ ( DIM * 2 * nPtls * deltaT * limitCorrAv); 
	/*-----------------------------------------------------------------------------
	 *   rrMSDAv -> mean square displacemnt 
	 *   rrMQDAv -> mean quadropole displacemnt 
	 *-----------------------------------------------------------------------------*/
	//				fac = 1./ ( DIM * 2 * nPtls * deltaT * limitCorrAv); 
	scale_factor = 1./ ( nPtls *  countCorrAv); 
	real factor_dig = 1./(3.*countCorrAv * g_Vol);
	real factor_offdig = 1./(6.*countCorrAv * g_Vol);

#pragma omp parallel for
	for (int nt = 1; nt < nCTime; nt ++) {
		rrMSDAv[nt] *= scale_factor;
		rrMQDAv[nt] *= scale_factor;
		real_tensor_product_r1_r0r1(&rrMSR1_R_Av[nt], scale_factor,
				&rrMSR1_R_Av[nt]);
		
		/*!
		 *  \brief  if all mass of particles is same value
		 */
		real_tensor_product_r2_r0r2(&rrMSR2_VR_Av[nt]
				, (.5*mass*mass),&rrMSR2_VR_Av[nt]);
		rrMSR2_VR_Av_dig[nt] = 
			factor_dig*	real_tensor_sum_dig_r2(&rrMSR2_VR_Av[nt]);
		rrMSR2_VR_Av_offdig[nt] = 
			factor_offdig*  real_tensor_sum_offdig_r2(&rrMSR2_VR_Av[nt]);

		for ( int nr=0; nr<nCSpatial; nr++) {
			avDrTable[nr][nt] *= factorDr[nr];
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

	prePrintProcess ();
	extern real kVal;
	int  nType,  k2, nr;
	Number_call_Print ++;
	//	char *header[] = {"cur-long", "cur-trans", "density", "vanHove-self"};
	char *header[] = {
		"full-density"   	 ,                        // 0
		"self-density"			,                       // 1
		"cross-density" 		,                       // 2
		"self-vanHove"                              // 3
	};
	fprintf (fp, "%s\n",txtCorr);
	//for (nType = 0; k < 3; k ++) {
	for (k2 = 0; k2 < sizeof(header)/ sizeof(char*); k2 ++) {
		/* 		fprintf (fp, "%s", header[nType]);
		 * 		for (j = 0; j < nCSpatial; j ++)
		 * 			fprintf (fp, " %7.3f", kVal*(j+1));
		 * 		fprintf (fp, "\n");
		 */

		//    EvalOtherInformation ();
		fprintf (fp, "# %s %7.3f %7.3f %7.3f\n", header[k2] , kVal, 1.0*deltaT, rVal);
		switch ( k2) {
			case 0: 
				/*!-----------------------------------------------------------------------------
				 *  avF_qq2[AVDOF*i+nType][k] -> F(q_i,t_k) 
				 *-----------------------------------------------------------------------------*/
				nType= 0;
				for (int nt = 0; nt < nCTime; nt ++) {
					/* 			deltaT = n *1. * deltaT;
					 * 			fprintf (fp, "%7.3f", deltaT);
					 */
					for (int nk = 0; nk < nCSpatial; nk ++){
						fprintf (fp, " %8.4e", avF_qq2[AVDOF * nk + nType][nt]);
					}
					fprintf (fp, "\n");
				} 
				break;
/*-----------------------------------------------------------------------------
 *  avF_s_qq2[3*i+nType][j] -> F_s(q_i,t_j) 
 *-----------------------------------------------------------------------------*/
			case 1: 
				nType =  0;
				for (int nt = 0; nt < nCTime; nt ++) {
					/* 			deltaT = n *1. * deltaT;
					 * 			fprintf (fp, "%7.3f", deltaT);
					 */
					for (int nr = 0; nr < nCSpatial; nr ++){
						fprintf (fp, " %8.4e", avF_s_qq2[AVDOF * nr + nType][nt]);
					}
					fprintf (fp, "\n");
				} 
				break;
/*-----------------------------------------------------------------------------
 *  avF_d_qq2[3*i+nType][j] -> F_d(q_i,t_j) 
 *  magnetic
 *-----------------------------------------------------------------------------*/
			case 2:
				nType = 0;
				for (int nt = 0; nt < nCTime; nt ++) {
					/* 			deltaT = n *1. * deltaT;
					 * 			fprintf (fp, "%7.3f", deltaT);
					 */
					for (int nr = 0; nr < nCSpatial; nr ++){
						fprintf (fp, " %8.4e", avF_d_qq2[AVDOF * nr + nType][nt]);
					}
					fprintf (fp, "\n");
				} 
				break;
			case 3: 
				//        fprintf (fp, "#van Hove function\n");

				for (int nt = 0; nt < nCTime; nt ++) {
					for ( nr=0; nr<nCSpatial; nr++)  {
						fprintf (fp, " %8.4e", avDrTable[nr][nt] );
					}
					fprintf (fp, "\n");
				}
				break;
		}
		fprintf (fp, "\n");
	}
	void PrintEtc();
	PrintEtc ();
}
void PrintSpacetimeCorr_binary (FILE *fp)
	/*!
	 *  \brief  결과를 출력하는 함수
	 *
	 *  \param  fp output file descriptor
	 */
{

	extern real kVal;
	int  nType,  k2, nr;
//	Number_call_Print ++;
	//	char *header[] = {"cur-long", "cur-trans", "density", "vanHove-self"};
	char *header[] = {
		"full-density"   	 ,                        // 0
		"self-density"			,                       // 1
		"cross-density" 		,                       // 2
		"self-vanHove"                              // 3
	};
	int nTypes = sizeof(header)/ sizeof(char*);
	fwrite (txtCorr,sizeof(char),strlen(txtCorr),fp);
	fwrite (&nTypes,sizeof(int) ,1,fp);
//	fwrite (fp, "%s\n",txtCorr);
	//for (nType = 0; k < 3; k ++) {
	for (k2 = 0; k2 < nTypes; k2 ++) {
		/* 		fprintf (fp, "%s", header[nType]);
		 * 		for (j = 0; j < nCSpatial; j ++)
		 * 			fprintf (fp, " %7.3f", kVal*(j+1));
		 * 		fprintf (fp, "\n");
		 */

		//    EvalOtherInformation ();
		real col2 =  kVal; 
		real col3 = 1.0*deltaT; 
		real col4 = rVal; 
		fwrite(header[k2],sizeof(char), strlen(header[k2]), fp);

		fwrite(&col2, sizeof(real),1,fp);
		fwrite(&col3, sizeof(real),1,fp);
		fwrite(&col4, sizeof(real),1,fp);

/* 		fprintf (fp, "# %s %7.3f %7.3f %7.3f\n", header[k2] , kVal, 1.0*deltaT, rVal);
 */
		switch ( k2) {
			case 0: 
				/*!-----------------------------------------------------------------------------
				 *  avF_qq2[AVDOF*i+nType][k] -> F(q_i,t_k) 
				 *-----------------------------------------------------------------------------*/
				nType= 0;
				for (int nt = 0; nt < nCTime; nt ++) {
					/* 			deltaT = n *1. * deltaT;
					 * 			fprintf (fp, "%7.3f", deltaT);
					 */
					for (int nk = 0; nk < nCSpatial; nk ++){
						fwrite( &(avF_qq2[AVDOF * nk + nType][nt]), sizeof(real),1,fp);
					}
/* 					for (int nk = 0; nk < nCSpatial; nk ++){
 * 						fprintf (fp, " %8.4e", avF_qq2[AVDOF * nk + nType][nt]);
 * 					}
 * 					fprintf (fp, "\n");
 */
				} 
				break;
/*-----------------------------------------------------------------------------
 *  avF_s_qq2[3*i+nType][j] -> F_s(q_i,t_j) 
 *-----------------------------------------------------------------------------*/
			case 1: 
				nType =  0;
				for (int nt = 0; nt < nCTime; nt ++) {
					/* 			deltaT = n *1. * deltaT;
					 * 			fprintf (fp, "%7.3f", deltaT);
					 */
					for (int nk = 0; nk < nCSpatial; nk ++){
						fwrite( &(avF_s_qq2[AVDOF * nk + nType][nt]), sizeof(real),1,fp);
					}
/* 					for (int nr = 0; nr < nCSpatial; nr ++){
 * 						fprintf (fp, " %8.4e", avF_s_qq2[AVDOF * nr + nType][nt]);
 * 					}
 * 					fprintf (fp, "\n");
 */
				} 
				break;
/*-----------------------------------------------------------------------------
 *  avF_d_qq2[3*i+nType][j] -> F_d(q_i,t_j) 
 *  magnetic
 *-----------------------------------------------------------------------------*/
			case 2:
				nType = 0;
				for (int nt = 0; nt < nCTime; nt ++) {
					/* 			deltaT = n *1. * deltaT;
					 * 			fprintf (fp, "%7.3f", deltaT);
					 */
					for (int nk = 0; nk < nCSpatial; nk ++){
						fwrite( &(avF_d_qq2[AVDOF * nk + nType][nt]), sizeof(real),1,fp);
					}
/* 					for (int nr = 0; nr < nCSpatial; nr ++){
 * 						fprintf (fp, " %8.4e", avF_d_qq2[AVDOF * nr + nType][nt]);
 * 					}
 * 					fprintf (fp, "\n");
 */
				} 
				break;
			case 3: 
				//        fprintf (fp, "#van Hove function\n");

				for (int nt = 0; nt < nCTime; nt ++) {
					for (int nk = 0; nk < nCSpatial; nk ++){
						fwrite( &(avDrTable[nr][nt]), sizeof(real),1,fp);
					}
/* 					for ( nr=0; nr<nCSpatial; nr++)  {
 * 						fprintf (fp, " %8.4e", avDrTable[nr][nt] );
 * 					}
 * 					fprintf (fp, "\n");
 */
				}
				break;
		}
/* 		fprintf (fp, "\n");
 */
	}
/* 	void PrintEtc();
 * 	PrintEtc ();
 */
}
void PrintEtc () {

	//  char filename1[100] ="Dq00.info" ;
	//  char filename2[100] ="Ft00.info" ;
	char filename1[100];
	char filename2[100];
	char filename3[100];
	int nfile = 0;
	sprintf(filename1, "Dt%03d.info",nfile);
	sprintf(filename2, "vanHove%03d.info",nfile);
	sprintf(filename3, "SSF%03d.info",nfile);
	//printf( "access(%s) -> return %d", filename1, access(filename1,F_OK));
	//		
	//		{
	//			while( 0 == access(filename1,F_OK) ) {
	//				/* 		fprintf(stderr, "Files are  exist at least . (%03d) \n", nfile);
	//				 * 		sleep(1);
	//				 */
	//				nfile++;
	//				sprintf(filename1, "Dt%03d.info",nfile);
	//				sprintf(filename2, "vanHove%03d.info",nfile);
	//				sprintf(filename3, "SSF%03d.info",nfile);
	//			}
	//		}

	/* 	FILE* fp_Dq = fopen(filename1,"w");
	 * 	fprintf (fp_Dq, "# dt = %7.3f\n", deltaT);
	 * 	for (j = 0; j < nCSpatial; j ++) {
	 * 		fprintf (fp_Dq, "%8.4f" , j*kVal );
	 * 		for (n = 1; n < nCTime; n ++) {   
	 * 			fprintf (fp_Dq, " %8.4e" ,  valDqt[j][n]);
	 * 		}
	 * 		fprintf (fp_Dq, "\n");
	 * 	}
	 * 	fclose(fp_Dq);
	 */

	/* 	FILE* fp_Ft = fopen(filename2,"w");
	 * 	fprintf (fp_Ft, "# dq = %7.3e\n", kVal);
	 * 	for (n = 0; n < nCTime; n ++) {   
	 * 		fprintf (fp_Ft, "%8.4f" , n*deltaT );
	 * 		for (j = 0; j < nCSpatial; j ++) {
	 * 			fprintf (fp_Ft, " %8.4e" ,  avF_qq2[(3*j)+2][n]/avF_qq2[(3*j)+2][0]);
	 * 		}
	 * 		fprintf (fp_Ft, "\n");
	 * 	}
	 * 	fclose(fp_Ft);
	 */

	FILE* fp_SSF = fopen(filename3,"w");
	for (int nr = 0; nr < nCSpatial; nr ++) {
		fprintf (fp_SSF, "%8.4f" " %8.4e""\n" , (nr+1)*kVal , 
				//				avF_qq2[(AVDOF*nr)+AV_DEN][0]);
						avF_qq2[(AVDOF*nr)+0][0]);
	}
	fclose(fp_SSF);

	//  fprintf (fp_SSF, "# dq = %7.3e\n", kVal);




	FILE* fp_Dt = fopen(filename1,"w");
	fprintf (fp_Dt, "#time MSD msdx msdy msdz D(t) MQD  MSVR_dig MSVR_offdig xy yx zy yz xz zx\n");

	real fac = 1./( 2.* deltaT * DIM * 2);
	int nr=0;
	rrDt[nr] = fac*(-rrMSDAv[nr+2]  +4.*rrMSDAv[nr+1] - 3.* rrMSDAv[nr]);
	for ( nr = 1; nr < nCTime-1; nr += 1 ) {
		rrDt[nr] = fac*(rrMSDAv[nr+1]  -rrMSDAv[nr-1] );
	}
	nr=nCTime-1;
	rrDt[nr] = fac*(rrMSDAv[nr-2]  -4.*rrMSDAv[nr-1] + 3.* rrMSDAv[nr]);


	for ( int  nt = 0; nt < nCTime; nt += 1 ) {
		real tVal = nt * deltaT;
		fprintf (fp_Dt, "%8.4f %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e \n", 
				tVal, rrMSDAv[nt] , 
				rrMSR1_R_Av[nt].x, rrMSR1_R_Av[nt].y, rrMSR1_R_Av[nt].z,
				rrDt[nt], rrMQDAv[nt],
				rrMSR2_VR_Av_dig[nt], rrMSR2_VR_Av_offdig[nt]
				, rrMSR2_VR_Av[nt].xy
				, rrMSR2_VR_Av[nt].yx
				, rrMSR2_VR_Av[nt].zy
				, rrMSR2_VR_Av[nt].yz
				, rrMSR2_VR_Av[nt].xz
				, rrMSR2_VR_Av[nt].zx
				); 
	}
	fclose(fp_Dt); 

	/*-----------------------------------------------------------------------------
	 *  van Hove function part
	 *-----------------------------------------------------------------------------*/
	/* 	fprintf (fp_Gr, "#van Hove function\n");
	 * 	FILE* fp_Gr = fopen(filename2,"w");
	 * 	
	 * 	for ( nr=0; nr<nCSpatial; nr++)  {
	 * 		for (j = 0; j < nCTime; j ++) {
	 * 			fprintf (fp_Gr, " %8.4e", avDrTable[nr][j] );
	 * 		}
	 * 		fprintf (fp_Gr, "\n");
	 * 	}
	 * 	fprintf (fp_Gr, "\n");
	 * 	fclose(fp_Gr);
	 */

}

void ZeroOneTimeCorr(MakeSqtClass* cl_sqt)
{
	TBuf* tBuf = cl_sqt->tBuf;
	real *	rho_q1  = tBuf->rho_q1 ;
	real **	rho_s_q1  = tBuf->rho_s_q1 ;
	real **	rho_d_q1 = tBuf->rho_d_q1;
	for (int j = 0; j < FDOF * nCSpatial; j ++) {
		rho_q1[j] = 0.;
	}
	if ( flagSelf ) {
		for (int n=0; n<nPtls; n++) {
			for (int j = 0; j < FDOF * nCSpatial; j ++) {
				rho_s_q1[n][j] = 0.;
				rho_d_q1[n][j] = 0.;
			}
		}
	}
	real_tensor_zero_r2(&tBuf->sumVR_ct);
}
void EvalOneTimeSumVR(MakeSqtClass* cl_sqt) 
{
	Rank2R3 VR;  
	VecR3  vecr3,vel;
	atom* col_i;
	TBuf* tBuf = cl_sqt->tBuf;
	Snapshot* snap = cl_sqt->snap;
	for (int n=0; n<nPtls; n++) {
		col_i = &(snap->atoms[n]);

		vecr3.x = col_i->x;
		vecr3.y = col_i->y;
		vecr3.z = col_i->z;
		vel.x = col_i->vx;
		vel.y = col_i->vy;
		vel.z = col_i->vz;

		real_tensor_product_r2_r1r1 (& VR, &vel, &vecr3);
		real_tensor_increase_r2_r2(&tBuf->sumVR_ct, &VR);
	}
}

void EvalOneTimeKspace(MakeSqtClass* cl_sqt)
{
	real r[3], v[3],mu[3];
	Snapshot* snap = cl_sqt->snap;
	TBuf* tBuf = cl_sqt->tBuf;
	real *	rho_q1  = tBuf->rho_q1 ;
	real *	rho_s_q1_temp  = tBuf->rho_s_q1_temp ;
	real **	rho_s_q1  = tBuf->rho_s_q1 ;
	real **	rho_d_q1 = tBuf->rho_d_q1;
	atom* col_i;

/*-----------------------------------------------------------------------------
 *  Direct calculate  rho(q)
 *-----------------------------------------------------------------------------*/
	for (int n=0; n<nPtls; n++) {
		col_i = &(snap->atoms[n]);

		r[0]  =  col_i->x;   r[1] = col_i->y;   r[2]  = col_i->z; 
		v[0]  =  col_i->vx;  v[1] = col_i->vy;  v[2]  = col_i->vz;
		mu[0] = col_i->mux; mu[1] = col_i->muy; mu[2] = col_i->muz;
		// 
#pragma omp parallel  for 
		for (int k = 0; k < N_AXIS; k ++) {          
			real b,c,s,c0,c1,s1,c2,s2;
			coordi3* c3=&coordi_list[k];
			for (int m = 0; m < nCSpatial; m ++) {  
			  int 	marker = (nCSpatial* DOF)*k + DOF*m;
				// Because the time of integer calculation  is small, 
				// change the code more explicity
#ifndef SLOW_FULL_MATH
				if (m == 0) {
					b = kVal * ( r[0]*(*c3)[0] + r[1]*(*c3)[1] + r[2]*(*c3)[2]);
					c = cos (b);
					s = sin (b);
					c0 = c;
				} else if (m == 1) {
					c1 = c;
					s1 = s;
					c = 2. * c0 * c1 - 1.; //cos(2x)=2cos(x)^2-1
					s = 2. * c0 * s1;      // sin(2x)=2sin(x)cos(x)
				} else {
					c2 = c1;
					s2 = s1;
					c1 = c;
					s1 = s;
					c = 2. * c0 * c1 - c2; // Check true!!
					s = 2. * c0 * s1 - s2; // Check true!!
				}
#else
				b = kVal * r[k];
				c = cos( (m+1)*b);
				s = sin( (m+1)*b);
#endif
				rho_s_q1_temp[marker +0] = c;
				rho_s_q1_temp[marker +1] = s;
			}// loop spatial slice
		}  // for DIMEN
		   //	memcpy(rho_s_q1, rho_q1,sizeof(real)*24*nCSpatial);
		for(int nk=0; nk< FDOF * nCSpatial; nk++ ) {
			rho_q1 [ nk] += rho_s_q1_temp[nk];
		}
		if ( flagSelf ) {
			memcpy(rho_s_q1[n], rho_s_q1_temp, sizeof(real)*FDOF*nCSpatial);
		}
	} /* for loop : n<nPtls */

	if ( flagSelf ) {
		for(int nk=0; nk< FDOF * nCSpatial; nk++ ) {
			for (int n=0; n<nPtls; n++) {
				rho_d_q1[n] [ nk] = rho_q1[nk] - rho_s_q1 [n][nk];
			}
		}
	}
}

void EvalOneTimeCorr(MakeSqtClass* cl_sqt)
	/*!
	 * 
	 *  \brief  one time Correlation을 계산한다. 
	 * q space value는 x, y, z 방향 세개의 방향으로 
	 * longitudinal version, translational version과 가장 기본적인 방향성분 없는 density
	 */
{
	void ZeroOneTimeCorr(MakeSqtClass* cl_sqt);
	void EvalOneTimeSumVR(MakeSqtClass* cl_sqt) ;
	void EvalOneTimeKspace(MakeSqtClass* cl_sqt);

	ZeroOneTimeCorr(cl_sqt);

	EvalOneTimeSumVR(cl_sqt);
	EvalOneTimeKspace(cl_sqt);
}

void SetWaitedTimeCorr(MakeSqtClass* cl_sqt, TBuf* tBuf_tw)
{
	/*-----------------------------------------------
	 *   t_w information 
	 *-----------------------------------------------*/
	TBuf* tBuf = cl_sqt->tBuf;
	Snapshot* snap = cl_sqt->snap;
	real *	rho_q1  = tBuf->rho_q1 ;
	real **	rho_s_q1  = tBuf->rho_s_q1 ;
	real **	rho_d_q1 = tBuf->rho_d_q1;

	real_tensor_copy_r2r2(& tBuf_tw->orgSumVR, &tBuf->sumVR_ct);
	for (int n=0; n<nPtls; n++) {
		tBuf_tw->orgR[n].x = snap->atoms[n].x;
		tBuf_tw->orgR[n].y = snap->atoms[n].y;
		tBuf_tw->orgR[n].z = snap->atoms[n].z;
	}

	for (int j = 0; j < FDOF * nCSpatial; j ++){
		tBuf_tw->org_rho_q1[j] = rho_q1[j];
		if ( flagSelf ) {
			for (int n=0; n<nPtls; n++) {
				tBuf_tw->org_rho_s_q1[n][j] = rho_s_q1[n][j];
				tBuf_tw->org_rho_d_q1[n][j] = rho_d_q1[n][j];
			}  // for n
		} // if flagSelf
	}   // for j
}
void InitTwoTimeCorr (MakeSqtClass* cl_sqt, TBuf* tBuf_tw, int subtime)
{
	/*------------------------------
	 *  Zero initializing
	 *-----------------------------*/
	tBuf_tw->rrMSD[subtime]= 0.;
	tBuf_tw->rrMQD[subtime]= 0.;
	real_tensor_zero_r1 (&tBuf_tw->rrMSR1_R[subtime]);
	real_tensor_zero_r2 (&tBuf_tw->rrMSR2_VR[subtime]);

	for (int  nr=0; nr<nCSpatial; nr++) {
		tBuf_tw->DrTable[nr][subtime] =0;
	}

			//F_qq2 0  KSpace
	for (int k = 0; k < AVDOF * nCSpatial; k ++) {
		tBuf_tw->F_qq2[k][subtime] = 0.;
		tBuf_tw->F_s_qq2[k][subtime] = 0.;
		tBuf_tw->F_d_qq2[k][subtime] = 0.;
	}
}
void EvalTwoTimeEach(MakeSqtClass* cl_sqt, TBuf* tBuf_tw, int subtime)
{
	Snapshot* snap = cl_sqt->snap;
	for (int n=0; n<nPtls; n++) {
		VecR3 dr;
		real dx2,dy2,dz2,dr2;
		atom* col_i = &(snap->atoms[n]);

		dr.x =  col_i->x-tBuf_tw->orgR[n].x ;
		dr.y =  col_i->y-tBuf_tw->orgR[n].y ;
		dr.z =  col_i->z-tBuf_tw->orgR[n].z ;

		dx2 = dr.x*dr.x; 
		dy2 = dr.y*dr.y; 
		dz2 = dr.z*dr.z; 
		dr2 = dx2 + dy2 + dz2;

		int  i_Dr    = floor (sqrt(dr2)/rVal);
		if (i_Dr<nCSpatial) tBuf_tw->DrTable[i_Dr][subtime] ++;

		tBuf_tw->rrMSD[subtime] += dr2;
		tBuf_tw->rrMSR1_R[subtime].x += dx2;
		tBuf_tw->rrMSR1_R[subtime].y += dy2;
		tBuf_tw->rrMSR1_R[subtime].z += dz2;
		tBuf_tw->rrMQD[subtime] += dr2*dr2;
	}
}
void EvalTwoTimeCollective(MakeSqtClass* cl_sqt, TBuf* tBuf_tw, int subtime)
{
	Rank2R3 subVR,sqVR;
	TBuf* tBuf = cl_sqt->tBuf;
	real_tensor_sub_r2_r2r2(&subVR, &tBuf->sumVR_ct, &tBuf_tw->orgSumVR);
	real_tensor_product_r2_r2r2 (& sqVR, & subVR, & subVR);

	real_tensor_increase_r2_r2(&tBuf_tw->rrMSR2_VR[subtime],
			&sqVR);
}

void EvalTwoTimeKSpace(MakeSqtClass* cl_sqt, TBuf* tBuf_tw, int subtime)
{
	TBuf* tBuf = cl_sqt->tBuf;
	real *	rho_q1  = tBuf->rho_q1 ;
	real **	rho_s_q1  = tBuf->rho_s_q1 ;
	real **	rho_d_q1 = tBuf->rho_d_q1;
	for (int axis_b = 0; axis_b < N_AXIS; axis_b ++) { // 3 loop
#pragma omp parallel for 
		for (int nk = 0; nk < nCSpatial; nk ++) {
			const int avMarker = nk*AVDOF;
			const int marker = (nCSpatial* DOF)*axis_b + DOF*nk;
			int nc= 0;
//			for (int nc = 0; nc < 7; nc ++) {  //DOF/2 = 7
				int nav;
				real w;
				int ncos = marker + 2*nc;
				int nsin = marker + 2*nc +1;
/* 				if (nc < 3) {
 * 					int axis_a = nc;
 * 					if (axis_a == axis_b) {
 * 						w = 1.0;
 * 						nav = avMarker +V_LONG ;
 * 					}
 * 					else {
 * 						w = 0.5;    //   
 * 						nav = avMarker +V_TRANS ;
 * 					}
 * 					//              else w *= 0.5;
 * 				}
 * 				else if (nc<6) {
 * 					int axis_a = nc -3;
 * 					//              w = Sqr (kVal * (m + 1));
 * 					if (axis_a == axis_b) { // longitudinal
 * 						w = 1.0;
 * 						nav = avMarker + M_LONG;
 * 					}
 * 					else { //trasverse
 * 						w = 0.5;    //   
 * 						nav = avMarker +M_TRANS ;
 * 					}
 * 					//              else w *= 0.5;
 * 				}
 * 				else if (nc==6){
 * 					w = 1.;  
 * 					nav = avMarker + AV_DEN;
 * 				};   // density   3*m+4
 */
				w = 1.;
				nav = avMarker;
				// cos(q*r(t)) cos(q*r(t_w) +sin sin
				if (flagSelf ) {
					for (int n=0; n<nPtls; n++) {
						tBuf_tw->F_s_qq2[nav][subtime] +=
							w * (rho_s_q1[n][ncos] * tBuf_tw->org_rho_s_q1[n][ncos] +
									rho_s_q1[n][nsin] * tBuf_tw->org_rho_s_q1[n][nsin]);
						tBuf_tw->F_d_qq2[nav][subtime] +=
							w * (rho_d_q1[n][ncos] * tBuf_tw->org_rho_s_q1[n][ncos] +
									rho_d_q1[n][nsin] * tBuf_tw->org_rho_s_q1[n][nsin])+
							w * (rho_s_q1[n][ncos] * tBuf_tw->org_rho_d_q1[n][ncos] +
									rho_s_q1[n][nsin] * tBuf_tw->org_rho_d_q1[n][nsin]);
					}
				}

				tBuf_tw->F_qq2[nav][subtime] +=
					w * (rho_q1[ncos] * tBuf_tw->org_rho_q1[ncos] +
							rho_q1[nsin] * tBuf_tw->org_rho_q1[nsin]);
//			}  // for nc, 
		}    // for nk , 
	} 
}
void EvalTwoTimeCorr(MakeSqtClass* cl_sqt, TBuf* tBuf_tw, int subtime)
{
	InitTwoTimeCorr(cl_sqt, tBuf_tw, subtime);
	EvalTwoTimeEach(cl_sqt, tBuf_tw, subtime);
	EvalTwoTimeCollective(cl_sqt, tBuf_tw, subtime);
	EvalTwoTimeKSpace(cl_sqt, tBuf_tw, subtime);
}
void EvalSpacetimeCorr(MakeSqtClass* cl_sqt)
	/*!
	 * 
	 *  \brief  space time correlation을 계산한다. 
	 * q space value는 x, y, z 방향 세개의 방향으로 
	 * longitudinal version, translational version과 가장 기본적인 방향성분 없는 density
	 * version 3개를 구함.
	 * PREV. $M_T(q,t)$ $M_L(q,t)$ 
	 * Todo. 
	 *  \param  Snapshot* Snapshot 포인터 
	 */
{
	extern real kVal;
	void EvalOneTimeCorr(MakeSqtClass* cl_sqt);
	TBuf* tBuf = cl_sqt->tBuf;
	Snapshot* snap = cl_sqt->snap;

	L = snap->box.xhigh- snap->box.xlow;
	g_Vol  = L*L*L;
	nPtls = snap->n_atoms;
	if (cl_sqt->flag_alloc_more ==0 ) {
		omp_set_lock(&write_lock);
		Alloc_more(cl_sqt);
		cl_sqt->flag_alloc_more =1;
		flag_global_alloc_more = 1;
		AllocMemCheck ();
		omp_unset_lock(&write_lock);
	}

	kVal = 2.*M_PI / L;

	EvalOneTimeCorr(cl_sqt);

	// End Calculate Current time value
	// Begin Two time corrlation function
	for (int nb = 0; nb < nCBuffer; nb ++) {
		if (tBuf[nb].count == 0) {
			SetWaitedTimeCorr(cl_sqt, &tBuf[nb]);
		}     // End   buffer count ==0

		if (tBuf[nb].count >= 0) {
			EvalTwoTimeCorr(cl_sqt,&tBuf[nb],tBuf[nb].count);
		}                        // End buffer count >=0
		++ tBuf[nb].count;
	}
	AccumSpacetimeCorr (cl_sqt);
}
void AllocMemCheck ()
{
	if (ErrorAllocMem == 1) {
		printf("Reserving memory Error!!!!!!\n");
		exit(1);
	}
}
void AllocArray (MakeSqtClass* cl_sqt)
	/*!
	 *  \brief   이름그대로 memory 할다함. 
	 *     				rho_q1 functions of q      
	 *     	등.
	 */
{
	int nb;

	if (flag_global_alloc ==0 ) {
		AllocMem2 (avF_s_qq2, AVDOF * nCSpatial, nCTime, real);
		AllocMem2 (avF_d_qq2, AVDOF * nCSpatial, nCTime, real);
		AllocMem2 (avF_qq2,  AVDOF * nCSpatial, nCTime, real);

		AllocMem2 (valDqt,  nCSpatial, nCTime, real);
		AllocMem2 (valGammaQT,  nCSpatial, nCTime, real);
	}
	AllocMem (cl_sqt->tBuf, nCBuffer, TBuf);
	TBuf* tBuf = cl_sqt->tBuf;
	AllocMem (tBuf->rho_q1, FDOF * nCSpatial, real);
	for (nb = 0; nb < nCBuffer; nb ++) {
		AllocMem (tBuf[nb].org_rho_q1, FDOF * nCSpatial, real);

		AllocMem2 (tBuf[nb].F_s_qq2, AVDOF * nCSpatial, nCTime, real);
		AllocMem2 (tBuf[nb].F_d_qq2, AVDOF * nCSpatial, nCTime, real);
		AllocMem2 (tBuf[nb].F_qq2, AVDOF * nCSpatial, nCTime, real);
	}
	/*!
	 *  \brief  Memory for Green-Kubo formula
	 */
	// AllocArray for Diffuse ()
	if (flag_global_alloc ==0 ) {
		AllocMem (rrMSDAv, nCTime, real);
		AllocMem (rrMSR1_R_Av , nCTime, VecR3);
		AllocMem (rrMQDAv, nCTime, real);
		// AllocArray for shear viscosity
		// 				 (diffusion of momentum)
		AllocMem (rrMSR2_VR_Av, nCTime, Rank2R3);
		AllocMem (rrMSR2_VR_Av_dig, nCTime, real);
		AllocMem (rrMSR2_VR_Av_offdig, nCTime, real);
		AllocMem (rrDt, nCTime, real);
		AllocMem2 (avDrTable, nCSpatial,nCTime, real);

		ZeroAvSpacetimeCorr();
	}

	fprintf(stderr, "Reserving memory on heap via AllocMem : %lld GB\n",  ll_mem_size/1000ll/1000ll/1000ll);
}
void Alloc_more (MakeSqtClass* cl_sqt) 
{
	/*!
	 *  \brief  Alloc_more 
	 *          Allocing   using nPtls  is post-process
	 *
	 */
	int nb,nr; real rho0, shell_Vol;

	TBuf* tBuf = cl_sqt->tBuf;

	AllocMem (tBuf->rho_s_q1_temp, FDOF * nCSpatial, real);

	if (flagSelf ) {
		AllocMem (tBuf->rho_s_q1, nPtls, real*);
		AllocMem (tBuf->rho_d_q1, nPtls, real*);

		for (int natom=0; natom <nPtls ; natom++) {
			AllocMem (tBuf->rho_s_q1[natom], FDOF * nCSpatial, real);
			AllocMem (tBuf->rho_d_q1[natom], FDOF * nCSpatial, real);
		}
	fprintf(stderr, "Reserving memory on heap via AllocMem : %lld GB\n",  ll_mem_size/1000ll/1000ll/1000ll);
	}
	for (nb = 0; nb < nCBuffer; nb ++) {
		AllocMem (tBuf[nb].orgR, nPtls, VecR3);
		AllocMem (tBuf[nb].rrMSD, nCTime, real);
		AllocMem (tBuf[nb].rrMQD, nCTime, real);
		AllocMem (tBuf[nb].rrMSR1_R, nCTime, VecR3);
		AllocMem (tBuf[nb].rrMSR2_VR, nCTime, Rank2R3);
		AllocMem2 (tBuf[nb].DrTable, nCSpatial,nCTime, int);

		if (flagSelf) {
			AllocMem (tBuf[nb].org_rho_s_q1, nPtls, real*);
			AllocMem (tBuf[nb].org_rho_d_q1, nPtls, real*);

			for (int natom=0; natom <nPtls ; natom++) {
				AllocMem (tBuf[nb].org_rho_s_q1[natom], FDOF * nCSpatial, real);
				AllocMem (tBuf[nb].org_rho_d_q1[natom], FDOF * nCSpatial, real);
			}
		}
	}
	fprintf(stderr, "Reserving memory on heap via AllocMem : %lld GB\n",  ll_mem_size/1000ll/1000ll/1000ll);
	if (flag_global_alloc_more ==0 ) {
		AllocMem (factorDr, nCSpatial, real);
		AllocMem (radius, nCSpatial, real);

		rho0 = nPtls/g_Vol;
		for (nr = 0; nr < nCSpatial; nr ++) {
			if (nr ==0) {
				shell_Vol = 4*M_PI /3. * pow(rVal,3);
			}
			else{
				shell_Vol = (4./3.)*M_PI * 
					( pow( (nr+1)*rVal,3)-pow(nr*rVal,3)) ; 
			}
			//		else shell_Vol = 4*M_PI * pow(rVal,3)* (nr*nr + 1./12.);
			// else 부분 확실히 해야함 최근에 다룬적 있음. 

			radius  [nr] = (nr+.5) * rVal;
			factorDr[nr] = 1./( pow(rho0,2) * g_Vol *shell_Vol*limitCorrAv);
			/* 		printf("rho0=%.2e, Vol=%.2e, shell_Vol=%.2e, factorDr=%.2e\n", 
			 * 				rho0,g_Vol,shell_Vol,factorDr[nr]);
			 */
		}
	} // if flag_global_alloc_more   

	fprintf(stderr, "Reserving memory on heap via AllocMem : %lld GB\n",  ll_mem_size/1000ll/1000ll/1000ll);
}

int GetNameList (int argc, char **argv)
	/*!
	 *  \brief    from book of rapaport 
	 *  	ex)	  input.file	
	 *  	Name  value 
	 *    Name2 value 
	 *  	 value type : real(double type) int 		 
	 *  			 Well defined.
	 */
{
	int  j, k, match, ok;
	char buff[100], *token;
	FILE *fp;
	strcpy (buff, inputFilename);
	//	strcpy (buff, argv[0]);
	//	strcat (buff, ".in");
	if ((fp = fopen (buff, "r")) == 0)  {
		fp = fopen(buff, "w");
		for (k = 0; k < sizeof (nameList) / sizeof (NameList); k ++) {
			fprintf (fp, "%s\t", nameList[k].vName);
			if (strlen (nameList[k].vName) < 8) fprintf (fp, "\t");
			for (j = 0; j < nameList[k].vLen; j ++) {
				switch (nameList[k].vType) {
					case N_I:
						fprintf (fp, "%d ", 0);
						//						fprintf (fp, "%d ", *NP_I);
						break;
					case N_R:
						fprintf (fp, "%#g ", 0.00);
						//						fprintf (fp, "%#g ", *NP_R);
						break;
				}
				fprintf (fp, "\n");
			}
		}
		fprintf (fp, "----\n");
		fclose(fp);
		printf("GetDataError\n");
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

	if(nCBuffer > nCTime ) nCBuffer = nCTime;

	for (k = 0; k < sizeof (nameList) / sizeof (NameList); k ++) {
		if (nameList[k].vStatus != 1) ok = 0;
	}
	return (ok);
}
void UpdateNameList ()
	/*!
	 *  \brief 초기값을 출력하는 함수 getNameList의 짝함수이다.   
	 *
	 *  \param  fp  FILE* file descriptor
	 */
{
	char buff[100];
	FILE* fp;
	strcpy (buff, inputFilename);
	fp = fopen( buff, "w");
	PrintNameList2File(fp);
	fclose(fp);
}
void PrintNameList2File (FILE *fp)
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
	 *          
	 *  고로 원래 목적과 달리 delta_k를 2pi/L * n(정수)로 맞추도록 한다.
	 *  \param  snap Snaptshot* 스냅샷 포인터
	 */
	extern real kVal;
	real new_dk;
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

	n_mul = round(kVal/ (2.*M_PI/ L[0] ));
	if (n_mul <=0) n_mul =1;
	new_dk = (2.*M_PI/L[0]) * n_mul;
	fprintf(stderr, "Update for input dk param: %f -> %f \n"
			, kVal, new_dk);
	kVal = new_dk;
}
