/*!
 *    \file  lmp2sqt.c
 *   \brief  
 *
 *  \author  KIM Hyeok (kh), ekh0324@gmail.com
 *
 *  \internal
 *       Created:  2017-05-29
 *      Revision:  2019-04-12
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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include	<unistd.h>
#include <stdbool.h>
#include <getopt.h>
#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif
#include <math.h>
#include "snapshot.h"
#include <time.h>


#define DIM 3
//#define flagSelf 1

typedef enum {N_I, N_R} VType;
typedef struct {
	char *vName;
	void *vPtr;
	VType vType;
	size_t vLen ;
	int vStatus;
} NameList;
#define NameI(x)                      \
{#x, &x, N_I, sizeof (x) / sizeof (int), 0}
#define NameR(x)                       \
{#x, &x, N_R, sizeof (x) / sizeof (real),0}

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

/*!
 *  \brief  Data List
 *  kVal
 *  deltaT
 *  mass
 *  rVal
 *  limitCorrAv
 *	limitCorrAv 
 *	nCBuffer  = number of simul. time seq
 *	nCSpatial = number of spatial seq
 *	nCTime    = number of time seq
 *	nCSkip    = number of time seq
 *
 */
#define DATALIST { \
	NameR   (kVal),   \
	NameR   (deltaT),   \
	NameR   (mass),   \
	NameR   (rVal),   \
	NameI  	(limitCorrAv),\
	NameI   (nCBuffer),   \
	NameI   (nCSpatial),   \
	NameI   (nCTime),   \
	NameI   (nCSkip)   \
}



NameList nameList[] = DATALIST;


bool flag_f=false;
bool flag_F=false;
bool flag_t=false;
bool flag_s=false;
bool flag_d=false;

bool flag_magnet=false;
bool flag_velocity=false;

char inputFilename[100]= "in.lmp2sqt";
void UpdateNameList ();
void PrintNameList2File (FILE *fp);
int GetNameList ();
int Number_call_Print =0;
int flag_global_alloc =0 ;
int flag_global_alloc_more =0 ;
void PrintSpacetimeCorr_binary ( FILE*);
char datetime_data[100] ;
int main(int argc, char** argv) {
	char filename[100];
	int n_snap;
	int opt_num, opt;
	
	int full_n_snaps=0, full_n_snaps_index=0,progress=-1;
	bool files_on [argc+3];
	long int num_data [argc+3];
	limitCorrAv = 0;

	time_t rawtime;
	struct tm* timeinfo;
	rawtime = time(NULL);
	timeinfo   =  localtime(&rawtime);
//	strftime(datetime_data,100,"%y%d%m_%H%M", timeinfo);
	strftime(datetime_data,100,"%y%d%m", timeinfo);
	char help_message[2000];
	sprintf(help_message,"%s -f -t -s fn_Dt filename 2 ... \n"
"	options: \n"
"		-f) collective full calculation required long time and large memory \n"
"		-F) collective self calculation required long time and large memory \n"
"			space-time, twotime, stress etc...\n"
"		-t) only twotime correlation \n"
"		-s) only twotime correlation + stress tensor \n"
"		-V) velocity correlation \n"
"		-M) magnet correlation \n"
"		-d) Get Amount of dynamic memory allocation \n"
"   -i) argment filename for input variable \n"
			, argv[0]);
		
	
	if( argc <2) {
		puts(help_message);
		return EXIT_FAILURE;
	}
	while (1) {
		opt = getopt (argc,argv, "tfFsVMi:d");
		if (opt == -1) break;
		switch (opt) {
			case 's' :
				puts("flag_s on : stress tensor calculation on");
				puts("require flag_t : also on");
				flag_t = true;
				flag_s = true; break;
			case 't' :
				puts("flag_t on : two time correlation");
				flag_t = true; break;
			case 'f' :
				puts("flag_f on : collective space time correlation ");
				flag_f = true; break;
			case 'F' :
				puts("flag_F on : self-collective space time correlation ");
				puts("flag_f on : collective space time correlation ");
				flag_f= flag_F = true; break;
			case 'V' :
				puts("flag_velocity on : velocity correlation on");
				flag_velocity = true; break;
			case 'd' :
				puts("flag_d on : get amount of dynamic memory allocation");
				flag_d = true; break;
			case 'M' :
				puts("flag_magnet on : magnet correlation on");
				flag_magnet = true; break;
			case 'i' :
				strcpy(inputFilename, optarg);
				printf("intput filename : %s\n" ,inputFilename);
				break;
			case 'h':
			case '?':
			default:
				puts(help_message);
			return EXIT_FAILURE;
			
		}
	}

	if(optind ==argc ) {
		puts(help_message);
		return EXIT_FAILURE;
	}
	
	int runable = 0;
	GetNameList();
	for( opt_num = optind;   opt_num < argc; opt_num++)  {
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
			bool check =	ReadDumpForOnlyCheck(fp);

			if (check == false )
				break;

			n_snap++; 
		}
		// if ( flag_Max_eval)
		if (n_snap <5){
			fprintf(stderr,"The # of snap is too small(<5) \n"
					", or this file is not valid \n"
					"We would  not use 	this file(%s)!!\n", filename);
			files_on [opt_num] = false;
			//			return 23;
		}
		else if ( floor((n_snap - nCTime) /(nCTime/ nCBuffer)) <1) {
			fprintf(stderr,"The # of snap is too small()\n"
					     "it dont make a infomation\n"
							 "floor((n_snap - nCTime) /(nCTime/ nCBuffer)) <1\n"
							 "floor((%d - %d) /(%d/ %d)) <1\n"
					"We would  not use 	this file(%s)!!\n", 
					n_snap,nCTime, nCTime, nCBuffer,
					filename);
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
			runable += 1;
		}
		fclose (fp);
	}

	if ( runable ==0 ) {
		puts("Check your input options or datafile length!!");
		exit(1);
	}
	PrintNameList2File(stderr);
	UpdateNameList ();
	int num_files =0;
	for( opt_num = optind;   opt_num < argc; opt_num++)  {
		if (files_on[opt_num] == true) num_files ++;
	}
  //   One process serially use each file.
	
	real  r_done_works;
	int  i_done_works;

	AllocMem(classSqt,1,MakeSqtClass);
	
	classSqt->flag_alloc = 0;
	classSqt->flag_alloc_more = 0;
	i_done_works= 0;


// recoding without openmp that is not suitable. 
	for( opt_num = optind;   opt_num < argc; opt_num++)  {
		if (files_on[opt_num] == false) continue;
		strcpy( filename,argv[opt_num]);
		FILE* fp = fopen( filename ,"r");
		Snapshot* snap, *firstSnap;

		MakeSqtClass* cl_sqt =  classSqt;
		InitSpacetimeCorr(cl_sqt);
		/*!
		 *  \brief  Start Calculation.
		 */
		rewind(fp);  
			firstSnap = ReadDump(fp);
		Init_reciprocal_space(firstSnap);
		rewind(fp);

		fprintf(stderr,"FULL_SNAPS = %5d\n",full_n_snaps);
		for (int ns = 0 ; ns < num_data[opt_num] ; ns++ ) {
			/* 		while(1) {}
			*/
			snap =	ReadDump(fp);
			if (snap == NULL)
				break;
			cl_sqt->snap = snap;
			
			EvalSpacetimeCorr(cl_sqt);

			FreeSnapshot(snap);
//#pragma omp atomic
			full_n_snaps_index++;

			i_done_works++;

			int new_progress =(10000.0*full_n_snaps_index/ full_n_snaps);
			r_done_works =
				(10000.0*i_done_works/ full_n_snaps);
			/* 			fprintf(stderr,"FULL_SNAPS = %5d, newprogress %d\n",full_n_snaps_index,
			 * 					new_progress);
			 */
			if (new_progress != progress ) {
				progress = new_progress;
				fprintf(stderr, "\r %4.2f%%(", progress*.01);
				fprintf(stderr, "main:%4.2f%%", r_done_works*.01);
				fflush(stderr);
			}
		}
		//		opt_num ++;
		fclose (fp);
		fprintf(stderr, "\nEnd : file : %s\n", filename);
	}


	if ( Number_call_Print ==0 ) {
		fprintf(stderr, "limit corr  = %d, countCorrAv = %d\n",
				limitCorrAv,countCorrAv);
		limitCorrAv =countCorrAv;
		MakeSqtClass* cl_sqt =  classSqt;
		PrintProcess(cl_sqt);
	}
  
	return 0;
}

void PrintProcess(MakeSqtClass* cl_sqt) 
{
	FILE* fout_log  = fopen( "output.log", "w+");
	fputs("Call PrintProcess : \n", fout_log);
	
	if(NULL== cl_sqt){
		fprintf(stderr,"FILE : %s,  FUNCTION : %s, LINE : %d, cl_sqt is NULL\n", __FILE__, __func__ , __LINE__);
		exit(1);
	}

	void prePrintProcess () ;
	void PrintEtc();

	prePrintProcess ();

	PrintEtc ();
	if (flag_f == true) {
		char* fout_filename = "output.binary";
		FILE* fout = fopen( fout_filename, "w+");
		FILE* fout_text = fopen( "output.txt", "w+");
		fprintf(fout_log, "%s %ld %d \n", fout_filename,ftell(fout), limitCorrAv);
		fprintf(stderr, "%s %ld %d \n", fout_filename,ftell(fout), limitCorrAv);
		PrintSpacetimeCorr (fout_text);
		PrintSpacetimeCorr_binary( fout);
		fclose(fout_text);
		fclose(fout);
	}
	fclose(fout_log);
	ZeroAvSpacetimeCorr ();  //FIXIT   아무의미없음 출력에서 교체합시다. 
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
			//omp_set_lock(&write_lock);
			// check!! that  data is full
			// S(q,t), M(q,t) part
			pt = &tBuf[nb];
			if (flag_f == true) {
				for (k = 0; k < AVDOF * nCSpatial; k ++) {
					for (n = 0; n < nCTime; n ++) {
						avF_qq2[k][n] += pt->F_qq2[k][n];
						StdDevF_qq2[k][n] += pt->F_qq2[k][n]*pt->F_qq2[k][n];
						if (flag_F == true ) {
							avF_s_qq2[k][n] += pt->F_s_qq2[k][n];
							StdDevF_s_qq2[k][n] += pt->F_s_qq2[k][n]*pt->F_s_qq2[k][n];

							avF_d_qq2[k][n] += pt->F_d_qq2[k][n];
							StdDevF_d_qq2[k][n] += pt->F_d_qq2[k][n]*pt->F_d_qq2[k][n];
						}
						if(flag_velocity==true){
							avC2_v_rho[k][n] += pt->C2_v_rho[k][n];
							avC2_rho_v[k][n] += pt->C2_rho_v[k][n];
							avC2_v_v[k][n] += pt->C2_v_v[k][n];
						}
						if(flag_velocity && flag_magnet){
							avC2_v_mu[k][n] += pt->C2_v_mu[k][n];
							avC2_mu_v[k][n] += pt->C2_mu_v[k][n];
						}
						if(flag_magnet==true){
							avC2_mu_rho[k][n] += pt->C2_mu_rho[k][n];
							avC2_rho_mu[k][n] += pt->C2_rho_mu[k][n];
							avC2_mu_mu[k][n] += pt->C2_mu_mu[k][n];
						}
					}
				}
				for (nt = 0; nt < nCTime; nt ++) {
					for ( nr=0; nr<nCSpatial; nr++) {
						avDrTable[nr][nt] += pt->DrTable[nr][nt];
					}
				}
			}
			// Diffuse Part
			if (flag_t == true) {
#pragma omp parallel for
				for (nt = 0; nt < nCTime; nt ++) {
					rrMSDAv[nt] += pt->rrMSD[nt];
					rrMSDCMAv[nt] += pt->rrMSDCM[nt];
					rrMQDAv[nt] += pt->rrMQD[nt];
					if (flag_velocity == true) {
						rrCvvAv[nt] += pt->rrCvv[nt];
						rrCvcmvcmAv[nt] += pt->rrCvcmvcm[nt];
					}
					if (flag_s == true) {
						real_tensor_increase_r1_r1(&rrMSR1_R_Av[nt], 
								&pt->rrMSR1_R[nt]);
						real_tensor_increase_r2_r2(&rrMSR2_VR_Av[nt], 
								&pt->rrMSR2_VR[nt]);
					}
				}
			}
			// buffer nb reset
			pt->count = 0;
			++ countCorrAv;
			if (countCorrAv == limitCorrAv) {
				PrintProcess (cl_sqt);
			}
			//omp_unset_lock(&write_lock);
		} //   if tBuf[nb].count is full
	}   // for all buffer
}

void InitSpacetimeCorr (MakeSqtClass* cl_sqt)
	/*!
	 *  \brief  프로그램 초기에 시간 평균을 낼 수 있도록 index를 부여하는 과정
	 */
{
	cl_sqt->nSkip =0;
	if (cl_sqt->flag_alloc == 0 ) {
		//omp_set_lock(&write_lock);
		if(flag_d == true ) {
			AmountAllocArray(cl_sqt);
		}
		AllocArray(cl_sqt);
		cl_sqt->flag_alloc = 1;
		flag_global_alloc =1;
		AllocMemCheck();
		//omp_unset_lock(&write_lock);
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
	if (flag_f == true) {
		for (nk = 0; nk < AVDOF * nCSpatial; nk ++) {
			for (nt = 0; nt < nCTime; nt ++) {
				avF_qq2[nk][nt] = 0.;
				StdDevF_qq2[nk][nt] = 0.;
				if (flag_F == true) {
					avF_s_qq2[nk][nt] = 0.;
					avF_d_qq2[nk][nt] = 0.;
					StdDevF_s_qq2[nk][nt] = 0.;
					StdDevF_d_qq2[nk][nt] = 0.;
				}
				if (flag_velocity == true) {
					avC2_v_rho [nk][nt] = 0. ;
					avC2_rho_v [nk][nt] = 0. ;
					avC2_v_v [nk][nt] = 0. ;
				}
				if (flag_velocity && flag_magnet ) {
					avC2_mu_v [nk][nt] = 0. ;
					avC2_v_mu [nk][nt] = 0. ;
				}
				if (flag_magnet == true) {
					avC2_mu_rho [nk][nt] = 0. ;
					avC2_rho_mu [nk][nt] = 0. ;
					avC2_mu_mu [nk][nt] = 0. ;
				}
			}
		}
		for (nt = 0; nt < nCTime; nt ++) {
			for (nr = 0; nr < nCSpatial; nr ++)  {
				avDrTable[nr][nt]= 0.;
			}
		}
	}
	if (flag_t == true) {
#pragma omp parallel for
		for (nt = 0; nt < nCTime; nt ++) {
			rrMSDAv[nt] = 0.; rrMQDAv[nt] = 0.;
			rrMSDCMAv[nt] = 0.;
			if (flag_velocity == true) {
				rrCvvAv[nt] = 0.; 
				rrCvcmvcmAv[nt] = 0.; 
			}
			if (flag_magnet == true) {
				rrCmmAv[nt] = 0.;
			}

			if (flag_s == true) {
				real_tensor_zero_r1( &rrMSR1_R_Av[nt] );
				real_tensor_zero_r2(	&rrMSR2_VR_Av [nt] ) ;
				rrMSR2_VR_Av_offdig[nt] = 0.;
				rrMSR2_VR_Av_dig [nt]   = 0.;
			}
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
	real inverseNptls = 1./nPtls;
	real inverseNptlsSq = 1./(nPtls*nPtls);
	real average_number = 1./(3.0*countCorrAv);
	if (flag_f == true) {
#pragma omp parallel for
		for (int nr = 0; nr < AVDOF * nCSpatial; nr ++) {
			for (int nt = 0; nt < nCTime; nt ++){
				real average = (avF_qq2[nr][nt] * inverseNptls) * average_number ;
				real averageSq = (StdDevF_qq2[nr][nt] * inverseNptlsSq) * average_number;
				avF_qq2[nr][nt] = average;
				StdDevF_qq2[nr][nt] = sqrt(   (averageSq - average*average) );
				//				ErrF_qq2[nr][nt] = sqrt( average_number * (averageSq - average*average) );
				if (flag_F == true) {
					/* 					avF_s_qq2[nr][nt] *= scale_factor;
					 * 					avF_d_qq2[nr][nt] *= 1.*scale_factor;
					 * 					StdDevF_s_qq2[nr][nt] *= scale_factor;
					 * 					StdDevF_d_qq2[nr][nt] *= 1. *scale_factor;
					 */

					average = (avF_s_qq2[nr][nt] * inverseNptls ) * average_number;
					averageSq = (StdDevF_s_qq2[nr][nt] * inverseNptlsSq )* average_number;
					avF_s_qq2[nr][nt] = average;
					StdDevF_s_qq2[nr][nt] = sqrt(  (averageSq - average*average) );
					//					ErrF_s_qq2[nr][nt] = sqrt( average_number * (averageSq - average*average) );

					average = (avF_d_qq2[nr][nt] * inverseNptls ) * average_number;
					averageSq = (StdDevF_d_qq2[nr][nt] * inverseNptlsSq) * average_number;
					avF_d_qq2[nr][nt] = average;
					StdDevF_d_qq2[nr][nt] = sqrt(  1. * (averageSq - average*average) );
					//					ErrF_d_qq2[nr][nt] = sqrt( average_number * 1. * (averageSq - average*average) );
				}
				if (flag_velocity == true) {
					avC2_v_rho [nr][nt] *= inverseNptls * average_number ;
					avC2_rho_v [nr][nt] *= inverseNptls * average_number ;
					avC2_v_v [nr][nt]   *= inverseNptls * average_number ;
				}
				if (flag_velocity && flag_magnet ) {
					avC2_mu_v [nr][nt] *= inverseNptls * average_number ;
					avC2_v_mu [nr][nt] *= inverseNptls * average_number ;
				}
				if (flag_magnet == true) {
					avC2_mu_rho [nr][nt] *= inverseNptls * average_number ;
					avC2_rho_mu [nr][nt] *= inverseNptls * average_number ;
					avC2_mu_mu [nr][nt] *= inverseNptls * average_number ;
				}
			}
		}
#pragma omp parallel for
		for (int nt = 1; nt < nCTime; nt ++) {
			for ( int nr=0; nr<nCSpatial; nr++) {
				avDrTable[nr][nt] *= factorDr[nr];
			}
		}
	}
	//				fac = 1./ ( DIM * 2 * nPtls * deltaT * limitCorrAv); 
	/*-----------------------------------------------------------------------------
	 *   rrMSDAv -> mean square displacemnt 
	 *   rrMQDAv -> mean quadropole displacemnt 
	 *-----------------------------------------------------------------------------*/
	//				fac = 1./ ( DIM * 2 * nPtls * deltaT * limitCorrAv); 
	real scale_countAv = 1./countCorrAv;
	real scale_4stress = .5*mass*mass/(countCorrAv* g_Vol);
	real scale_factor = 1./ ( nPtls *  countCorrAv); 
	real factor_Cvv = 1./(nPtls* countCorrAv*3.);
	real factor_msdcm = 1./( countCorrAv);
	real factor_Cvcmvcm = 1./( countCorrAv*3.);

	if (flag_t == true) {
#pragma omp parallel for
		for (int nt = 0; nt < nCTime; nt ++) {
			rrMSDAv[nt] *= scale_factor;
			rrMSDCMAv[nt] *= factor_msdcm;
			rrMQDAv[nt] *= scale_factor;
			if (flag_velocity == true ) {
				rrCvvAv[nt] *= factor_Cvv;
				rrCvcmvcmAv[nt] *= factor_Cvcmvcm;
			}
			/*!
			 *  \brief  if all mass of particles is same value
			 */


			if (flag_s == true) {
				real_tensor_product_r1_r0r1(&rrMSR1_R_Av[nt], scale_factor,
						&rrMSR1_R_Av[nt]);
				real_tensor_product_r2_r0r2(&rrMSR2_VR_Av[nt]
						, scale_4stress,&rrMSR2_VR_Av[nt]);

				rrMSR2_VR_Av_dig[nt] = 
					real_tensor_avg_dig_r2(&rrMSR2_VR_Av[nt]);
				rrMSR2_VR_Av_offdig[nt] = 
					real_tensor_avg_offdig_r2(&rrMSR2_VR_Av[nt]);

			}
		}
	}
	Number_call_Print ++;
}

void PrintSpacetimeCorr (FILE *fp)
	/*!
	 *  \brief  결과를 출력하는 함수
	 *
	 *  \param  fp output file descriptor
	 */
{

	extern real kVal;
	size_t k2;
	int  nType,   nr;
	//	char *header[] = {"cur-long", "cur-trans", "density", "vanHove-self"};
	char *header[] = {
		"full-density",  // 0
		"self-density",  // 1
		"cross-density", // 2
		"self-vanHove"   // 3
	};

	fprintf (fp, "%s\n",txtCorr);
	for (k2 = 0; k2 < sizeof(header)/ sizeof(char*); k2 ++) {

		//    EvalOtherInformation ();
		fprintf (fp, "# %s %7.3f %7.3f %7.3f\n", 
				header[k2] , kVal, 
				(nCSkip+1)*deltaT, 
				rVal);
		switch ( k2) {
			case 0: 
/*!-----------------------------------------------------------------------------
*  avF_qq2[AVDOF*i+nType][k] -> F(q_i,t_k) 
*-----------------------------------------------------------------------------*/
				nType= 0;
				for (int nt = 0; nt < nCTime; nt ++) {
					real time_d = nt *1. * deltaT;
					fprintf (fp, "%7.3f", time_d);

					for (int nk = 0; nk < nCSpatial; nk ++){
						fprintf (fp, " %8.4e %8.4e", avF_qq2[AVDOF * nk + nType][nt], StdDevF_qq2[AVDOF * nk + nType][nt]);
					}
					fprintf (fp, "\n");
				} 
				break;
				/*-----------------------------------------------------------------------------
				 *  avF_s_qq2[3*i+nType][j] -> F_s(q_i,t_j) 
				 *-----------------------------------------------------------------------------*/
			case 1: 
				if ( flag_F == true) {
					nType =  0;
					for (int nt = 0; nt < nCTime; nt ++) {
						real time_d = nt *1. * deltaT;
						fprintf (fp, "%7.3f", time_d);
						for (int nk = 0; nk < nCSpatial; nk ++){
							fprintf (fp, " %8.4e %8.4e", avF_s_qq2[AVDOF * nk + nType][nt],StdDevF_s_qq2[AVDOF * nk + nType][nt]);
						}
						fprintf (fp, "\n");
					} 
				}
				break;
/*-----------------------------------------------------------------------------
*  avF_d_qq2[3*i+nType][j] -> F_d(q_i,t_j) 
*  magnetic
*-----------------------------------------------------------------------------*/
			case 2:
				if (flag_F==true) {
					nType = 0;
					for (int nt = 0; nt < nCTime; nt ++) {
						real time_d = nt *1. * deltaT;
						fprintf (fp, "%7.3f", time_d);
						for (int nk = 0; nk < nCSpatial; nk ++){
							fprintf (fp, " %8.4e %8.4e", avF_d_qq2[AVDOF * nk + nType][nt], StdDevF_d_qq2[AVDOF * nk + nType][nt]);
						}
						fprintf (fp, "\n");
					} 
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
}

void PrintSpacetimeCorr_binary (FILE *fp)
	/*!
	 *  \brief  결과를 binary형태로 출력함
	 *
	 *  \param  fp output file descriptor
	 */
{

	extern real kVal;
	int  nType,  k2;
  int nr  __attribute__((unused));
//	Number_call_Print ++;
	//	char *header[] = {"cur-long", "cur-trans", "density", "vanHove-self"};
	char *header[] = {
		"full-density"   	 ,                        // 0
		"self-density"			,                       // 1
		"cross-density" 		,                       // 2
		"self-vanHove"                              // 3
	};
	int nTypes = sizeof(header)/ sizeof(char*);
	int header_txtCorr = strlen(txtCorr);

	fwrite (&header_txtCorr,sizeof(int),1,fp);
	fwrite (txtCorr,1,header_txtCorr,fp);

	fwrite (&nTypes,sizeof(int) ,1,fp);
	for (k2 = 0; k2 < nTypes; k2 ++) {

		real col2 =  kVal; 
		real col3 = 1.0*deltaT*(nCSkip+1); 
		real col4 = rVal; 
		int DataYes = 1;
		int DataNo = 0;
		int length = strlen(header[k2]);
		fwrite(&length, sizeof(int),1,fp);
		fwrite(header[k2],sizeof(char), length, fp);

		fwrite(&col2, sizeof(real),1,fp);
		fwrite(&col3, sizeof(real),1,fp);
		fwrite(&col4, sizeof(real),1,fp);

		switch ( k2) {
			case 0: 
/*!-----------------------------------------------------------------------------
*  avF_qq2[AVDOF*i+nType][k] -> F(q_i,t_k) 
*-----------------------------------------------------------------------------*/
				nType= 0;
				fwrite(&DataYes, sizeof(int),1,fp);
				for (int nt = 0; nt < nCTime; nt ++) {
					for (int nk = 0; nk < nCSpatial; nk ++){
						fwrite( &(avF_qq2[AVDOF * nk + nType][nt]), sizeof(real),1,fp);
					}
				} 
				break;
/*-----------------------------------------------------------------------------
 *  avF_s_qq2[3*i+nType][j] -> F_s(q_i,t_j) 
 *-----------------------------------------------------------------------------*/
			case 1: 
				if (flag_F==true) {
					nType =  0;
					fwrite(&DataYes, sizeof(int),1,fp);
					for (int nt = 0; nt < nCTime; nt ++) {
						for (int nk = 0; nk < nCSpatial; nk ++){
							fwrite( &(avF_s_qq2[AVDOF * nk + nType][nt]), sizeof(real),1,fp);
						}
					} 
				}
				else fwrite(&DataNo, sizeof(int),1,fp);
				break;
/*-----------------------------------------------------------------------------
 *  avF_d_qq2[3*i+nType][j] -> F_d(q_i,t_j) 
 *  magnetic
 *-----------------------------------------------------------------------------*/
			case 2:
				if (flag_F==true) {
					nType = 0;
					fwrite(&DataYes, sizeof(int),1,fp);
					for (int nt = 0; nt < nCTime; nt ++) {
						for (int nk = 0; nk < nCSpatial; nk ++){
							fwrite( &(avF_d_qq2[AVDOF * nk + nType][nt]), sizeof(real),1,fp);
						}
					}
				} 
				else fwrite(&DataNo, sizeof(int),1,fp);
				break;
			case 3: 
//        fprintf (fp, "#van Hove function\n");
				fwrite(&DataYes, sizeof(int),1,fp);
				for (int nt = 0; nt < nCTime; nt ++) {
					for (int nk = 0; nk < nCSpatial; nk ++){
						fwrite( &(avDrTable[nk][nt]), sizeof(real),1,fp);
					}
				}
				break;
		}
	}
}
void PrintEtc () {

	//  char fn_Dt[100] ="Dq00.info" ;
	//  char fn_vanHove[100] ="Ft00.info" ;
	char fn_Dt[200];
	char fn_vanHove[200];
	char fn_SSF[200];

	char fn_C2_v_rho[200];
	char fn_C2_rho_v[200];
	char fn_C2_v_v[200];

	char fn_C2_mu_v[200];
	char fn_C2_v_mu[200];

	char fn_C2_mu_rho[200];
	char fn_C2_rho_mu[200];
	char fn_C2_mu_mu[200];

	char filename_stress[200];
	int nfile = 0;
	do {
		sprintf(fn_Dt, "Dt%03d.info.%s",nfile,datetime_data);
		sprintf(filename_stress, "Stress%03d.info.%s",nfile,datetime_data);
		sprintf(fn_vanHove, "vanHove%03d.info.%s",nfile,datetime_data);
		sprintf(fn_SSF, "SSF%03d.info.%s",nfile,datetime_data);
		sprintf(fn_C2_v_v, "C2_v_v%03d.info.%s",nfile,datetime_data);
		sprintf(fn_C2_v_rho, "C2_v_rho%03d.info.%s",nfile,datetime_data);
		sprintf(fn_C2_rho_v, "C2_rho_v%03d.info.%s",nfile,datetime_data);

		sprintf(fn_C2_mu_v, "C2_mu_v%03d.info.%s",nfile,datetime_data);
		sprintf(fn_C2_v_mu, "C2_v_mu%03d.info.%s",nfile,datetime_data);

		sprintf(fn_C2_mu_rho, "C2_mu_rho%03d.info.%s",nfile,datetime_data);
		sprintf(fn_C2_mu_mu, "C2_mu_mu%03d.info.%s",nfile,datetime_data);
		sprintf(fn_C2_rho_mu, "C2_rho_mu%03d.info.%s",nfile,datetime_data);
		nfile++;
	} while( 0 == access(fn_Dt,F_OK) ) ;
	/* 	FILE* fp_Dq = fopen(fn_Dt,"w");
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

	/* 	FILE* fp_Ft = fopen(fn_vanHove,"w");
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

	if (flag_f == true) {
		FILE* fp_SSF = fopen(fn_SSF,"w");
		FILE *fp_C2_mu_mu, *fp_C2_v_v;
		FILE *fp_C2_mu_rho, *fp_C2_v_rho;
		FILE *fp_C2_rho_mu, *fp_C2_rho_v;
		FILE *fp_C2_v_mu, *fp_C2_mu_v;
		if (flag_velocity == true) {
			fp_C2_v_rho = fopen(fn_C2_v_rho,"w");
			fp_C2_rho_v = fopen(fn_C2_rho_v,"w");
			fp_C2_v_v   = fopen(fn_C2_v_v,"w");
		}
		if (flag_velocity && flag_magnet) {
			fp_C2_v_mu = fopen(fn_C2_v_mu,"w");
			fp_C2_mu_v = fopen(fn_C2_mu_v,"w");
		}
		if (flag_magnet == true) {
			fp_C2_mu_rho = fopen(fn_C2_mu_rho,"w");
			fp_C2_rho_mu = fopen(fn_C2_rho_mu,"w");
			fp_C2_mu_mu   = fopen(fn_C2_mu_mu,"w");
		}
		for (int nt = 0; nt < nCTime; nt ++) {
			int nType= 0;
			real time_d = nt *1. * deltaT;
			if (flag_velocity == true) {
				fprintf (fp_C2_v_v, "%7.3f", time_d);
				fprintf (fp_C2_rho_v, "%7.3f", time_d);
				fprintf (fp_C2_v_rho, "%7.3f", time_d);
			}
			if (flag_velocity && flag_magnet ) {
				fprintf (fp_C2_v_mu, "%7.3f", time_d);
				fprintf (fp_C2_mu_v, "%7.3f", time_d);
			}
			if (flag_magnet == true) {
				fprintf (fp_C2_mu_mu, "%7.3f", time_d);
				fprintf (fp_C2_rho_mu, "%7.3f", time_d);
				fprintf (fp_C2_mu_rho, "%7.3f", time_d);
			}
			for (int nk = 0; nk < nCSpatial; nk ++){
				if (flag_velocity == true) {
					fprintf (fp_C2_v_v, " %8.4e", avC2_v_v[AVDOF * nk + nType][nt]);
					fprintf (fp_C2_rho_v, " %8.4e", avC2_rho_v[AVDOF * nk + nType][nt]);
					fprintf (fp_C2_v_rho, " %8.4e", avC2_v_rho[AVDOF * nk + nType][nt]);
				}
				if (flag_velocity && flag_magnet ) {
					fprintf (fp_C2_v_mu, " %8.4e", avC2_v_mu[AVDOF * nk + nType][nt]);
					fprintf (fp_C2_mu_v, " %8.4e", avC2_mu_v[AVDOF * nk + nType][nt]);
				}
				if (flag_magnet == true) {
					fprintf (fp_C2_mu_mu, " %8.4e", avC2_mu_mu[AVDOF * nk + nType][nt]);
					fprintf (fp_C2_rho_mu, " %8.4e", avC2_rho_mu[AVDOF * nk + nType][nt]);
					fprintf (fp_C2_mu_rho, " %8.4e", avC2_mu_rho[AVDOF * nk + nType][nt]);
				}
			}
			if (flag_velocity == true) {
				fprintf (fp_C2_v_v, "\n");
				fprintf (fp_C2_rho_v, "\n");
				fprintf (fp_C2_v_rho, "\n");
			}
			if (flag_velocity && flag_magnet ) {
				fprintf (fp_C2_v_mu, "\n");
				fprintf (fp_C2_mu_v, "\n");
			}
			if (flag_magnet == true) {
				fprintf (fp_C2_mu_mu, "\n");
				fprintf (fp_C2_rho_mu, "\n");
				fprintf (fp_C2_mu_rho, "\n");
			}
		} 
		for (int nk = 0; nk < nCSpatial; nk ++){
			real value1 = avF_qq2[(AVDOF*nk)+0][0];
			real value2 = sqrt( StdDevF_qq2[(AVDOF*nk)+0][0]  - value1*value1 );
			fprintf (fp_SSF, "%8.4f" " %8.4g" " %8.4g""\n" , (nk+1)*kVal , 
					//				avF_qq2[(AVDOF*nr)+AV_DEN][0]);
							value1, value2);
		}
		fclose(fp_SSF);
	}

	//  fprintf (fp_SSF, "# dq = %7.3e\n", kVal);
	if (flag_s == true) {
		FILE* fp_stress = fopen(filename_stress,"w");
		fprintf (fp_stress, "#time MSVR_dig MSVR_offdig xy yx zy yz xz zx\n");
		fprintf (stderr, "print out stress file\n");


		for ( int  nt = 0; nt < nCTime; nt += 1 ) {
			real tVal = nt * deltaT* (nCSkip+1);
			fprintf (fp_stress, "%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g \n", 
					tVal,  
					rrMSR2_VR_Av_dig[nt], rrMSR2_VR_Av_offdig[nt]
					, rrMSR2_VR_Av[nt].xy
					, rrMSR2_VR_Av[nt].yx
					, rrMSR2_VR_Av[nt].zy
					, rrMSR2_VR_Av[nt].yz
					, rrMSR2_VR_Av[nt].xz
					, rrMSR2_VR_Av[nt].zx
					); 
		}
		fclose(fp_stress); 
	}


	if (flag_t == true) {
		FILE* fp_Dt = fopen(fn_Dt,"w");
		fprintf (fp_Dt, "#time MSD msdx msdy msdz D(t) anotherform MQD");
		if (flag_velocity==true) {
			fprintf (fp_Dt, " Cvv Cvcmvcm");
		}
		if (flag_magnet==true) {
			fprintf (fp_Dt, " Cmm");
		}
		fprintf ( fp_Dt, " MSDCM\n");
		fprintf (stderr, "print time correlation \n");

		real fac = 1./( 2.* deltaT* (nCSkip+1) * DIM * 2);
		int nr=0;
		rrDt[nr] = fac*(-rrMSDAv[nr+2]  +4.*rrMSDAv[nr+1] - 3.* rrMSDAv[nr]);

#pragma omp parallel for
		for ( nr = 1; nr < nCTime-1; nr += 1 ) {
			rrDt[nr] = fac*(rrMSDAv[nr+1]  -rrMSDAv[nr-1] );
		}
		nr=nCTime-1;
		rrDt[nr] = fac*(rrMSDAv[nr-2]  -4.*rrMSDAv[nr-1] + 3.* rrMSDAv[nr]);


		for ( int  nt = 0; nt < nCTime; nt += 1 ) {
			real tVal = nt * deltaT* (nCSkip+1);
			real anotherDt=0;
			if( nt !=0 ) { 
				anotherDt =   rrMSDAv[nt]/(  tVal * DIM * 2);
			}
//		fprintf (fp_Dt, "#time MSD msdx msdy msdz D(t) anotherform MQD");
			fprintf (fp_Dt, "%8.4f %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g", 
					tVal, rrMSDAv[nt] , 
					rrMSR1_R_Av[nt].x, rrMSR1_R_Av[nt].y, rrMSR1_R_Av[nt].z,
					rrDt[nt], anotherDt,rrMQDAv[nt]
					); 
			if (flag_velocity==true) {
				fprintf (fp_Dt, " %8.4f %8.4g", 
						rrCvvAv[nt] ,  rrCvcmvcmAv[nt]);
			}
			if (flag_magnet==true) {
				fprintf (fp_Dt, " %8.4f", rrCmmAv[nt]);
			}
			fprintf (fp_Dt, " %8.4f\n", rrMSDCMAv[nt]);
		}
		fclose(fp_Dt); 
	}

	/*-----------------------------------------------------------------------------
	 *  van Hove function part
	 *-----------------------------------------------------------------------------*/
	/* 	fprintf (fp_Gr, "#van Hove function\n");
	 * 	FILE* fp_Gr = fopen(fn_vanHove,"w");
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
	real *	kvel_q1  = tBuf->kvel_q1 ;
	real *	kmu_q1  = tBuf->kmu_q1 ;
	real **	rho_s_q1  = tBuf->rho_s_q1 ;
	real **	rho_d_q1 = tBuf->rho_d_q1;
	for (int j = 0; j < FDOF * nCSpatial; j ++) {
		rho_q1[j] = 0.;
		if (flag_velocity == true ) {
			kvel_q1[j] = 0.;
		}
		if (flag_magnet == true ) {
			kmu_q1[j] = 0.;
		}
	}
	if ( flag_F == true ) {
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
	Atom* col_i;
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
	real *	kvel_q1  = tBuf->kvel_q1 ;
	real *	kmu_q1  = tBuf->kmu_q1 ;
	real *	rho_s_q1_temp  = tBuf->rho_s_q1_temp ;
	real **	rho_s_q1  = tBuf->rho_s_q1 ;
	real **	rho_d_q1 = tBuf->rho_d_q1;
	Atom* col_i;

/*-----------------------------------------------------------------------------
 *  Direct calculate  rho(q)
 *  FIXIT : (~v)(q)   (~mu)(q) 
 *-----------------------------------------------------------------------------*/
	for (int n=0; n<nPtls; n++) {
		col_i = &(snap->atoms[n]);

		r[0]  =  col_i->x;   r[1] = col_i->y;   r[2]  = col_i->z; 
		v[0]  =  col_i->vx;  v[1] = col_i->vy;  v[2]  = col_i->vz;
		mu[0] = col_i->mux; mu[1] = col_i->muy; mu[2] = col_i->muz;
		// 
#pragma omp parallel  for 
		for (int k = 0; k < N_AXIS; k ++) {          
			real b,c,s,c0,c1,s1,c2,s2, kv,km;
			
			for (int m = 0; m < nCSpatial; m ++) {  
			  const int 	markerR = (nCSpatial* DOF)*k + DOF*m +0;
			  const int 	markerI = (nCSpatial* DOF)*k + DOF*m +1;
				const int qn = m+1;
				// Because the time of integer calculation  is small, 
				// change the code more explicity
#ifndef SLOW_FULL_MATH
				if (m == 0) {
//					b = kVal * ( r[0]*(*c3)[0] + r[1]*(*c3)[1] + r[2]*(*c3)[2]);
					b = kVal * r[k]; 
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
					c = 2. * c0 * c1 - c2; // cos(nx)=2*cos((n-1)x)*cos(x) - cos((n-2)x)
					s = 2. * c0 * s1 - s2; // sin(nx)=2*cos((n-1)x)*sin(x) - sin((n-2)x)
				}
#else
				b = kVal * r[k];
				c = cos( qn*b);
				s = sin( qn*b);
#endif
				if ( flag_velocity==true ) // velocity is defined 
				{
					kv = qn * kVal*  v[k];
					kvel_q1[markerR] += - kv * s;
					kvel_q1[markerI] += + kv * c;
				}
				if ( flag_magnet == true ) // magnetization  is defined 
				{
					km = qn * kVal*  mu[k];
					kmu_q1[markerR] += - km * s;
					kmu_q1[markerI] += + km * c;
				}
				rho_s_q1_temp[markerR] = c;
				rho_s_q1_temp[markerI] = s;
			}// loop spatial slice
		}  // for DIMEN
		   //	memcpy(rho_s_q1, rho_q1,sizeof(real)*24*nCSpatial);
		for(int nk=0; nk< FDOF * nCSpatial; nk++ ) {
			rho_q1 [ nk] += rho_s_q1_temp[nk];
		}
		if ( flag_F == true ) {
			memcpy(rho_s_q1[n], rho_s_q1_temp, sizeof(real)*FDOF*nCSpatial);
		}
	} /* for loop : n<nPtls */

	if ( flag_F == true ) {
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


	if (flag_f == true) {
		ZeroOneTimeCorr(cl_sqt);
		EvalOneTimeKspace(cl_sqt);
	}

	if (flag_s == true) {
		EvalOneTimeSumVR(cl_sqt);
	}
}

void SetWaitedTimeCorr(MakeSqtClass* cl_sqt, TBuf* tBuf_tw)
{
	/*-----------------------------------------------
	 *   t_w information 
	 *-----------------------------------------------*/
	TBuf* tBuf = cl_sqt->tBuf;
	Snapshot* snap = cl_sqt->snap;
	real *	rho_q1  = tBuf->rho_q1 ;
	real *	kvel_q1  = tBuf->kvel_q1 ;
	real *	kmu_q1  = tBuf->kmu_q1 ;
	real **	rho_s_q1  = tBuf->rho_s_q1 ;
	real **	rho_d_q1 = tBuf->rho_d_q1;

	if (flag_s == true) {
		real_tensor_copy_r2r2(& tBuf_tw->orgSumVR, &tBuf->sumVR_ct);
	}

	if (flag_t == true) {
#pragma omp parallel for
		for (int n=0; n<nPtls; n++) {
			tBuf_tw->orgR[n].x = snap->atoms[n].x;
			tBuf_tw->orgR[n].y = snap->atoms[n].y;
			tBuf_tw->orgR[n].z = snap->atoms[n].z;
			if (flag_velocity == true ) {
				tBuf_tw->orgV[n].x = snap->atoms[n].vx;
				tBuf_tw->orgV[n].y = snap->atoms[n].vy;
				tBuf_tw->orgV[n].z = snap->atoms[n].vz;
			}
			if (flag_magnet == true ) {
				tBuf_tw->orgMu[n].x = snap->atoms[n].mux;
				tBuf_tw->orgMu[n].y = snap->atoms[n].muy;
				tBuf_tw->orgMu[n].z = snap->atoms[n].muz;
			}
		}
	}

	if (flag_f == true) {
//#pragma omp parallel for 
		for (int j = 0; j < FDOF * nCSpatial; j ++){
			tBuf_tw->org_rho_q1[j] = rho_q1[j];
			if (flag_velocity == true) {
				tBuf_tw->org_kvel_q1[j] = kvel_q1[j];
			}
			if (flag_magnet == true) {
				tBuf_tw->org_kmu_q1[j] = kmu_q1[j];
			}
			if ( flag_F == true ) {
				for (int n=0; n<nPtls; n++) {
					tBuf_tw->org_rho_s_q1[n][j] = rho_s_q1[n][j];
					tBuf_tw->org_rho_d_q1[n][j] = rho_d_q1[n][j];
				}  // for n
			} // if flagSelf
		}   // for j
	}
}
void InitTwoTimeCorr (MakeSqtClass* cl_sqt, TBuf* tBuf_tw, int subtime)
{
	/*------------------------------
	 *  Zero initializing
	 *-----------------------------*/
	if (cl_sqt == NULL ){
		fprintf(stderr,"FILE : %s,  FUNCTION : %s, LINE : %d, cl_sqt is NULL\n", __FILE__, __func__ , __LINE__);
		exit(1);
	}
	if (flag_t == true) {
		tBuf_tw->rrMSD[subtime]= 0.;
		tBuf_tw->rrMQD[subtime]= 0.;
		if (flag_velocity == true) {
			tBuf_tw->rrCvv[subtime]= 0.;
			tBuf_tw->rrCvcmvcm[subtime]= 0.;
		}
		if (flag_magnet == true) {
			tBuf_tw->rrCmm[subtime]= 0.;
		}
		if (flag_s == true) {
			real_tensor_zero_r1 (&tBuf_tw->rrMSR1_R[subtime]);
			real_tensor_zero_r2 (&tBuf_tw->rrMSR2_VR[subtime]);
		}
	}

	if (flag_f == true) {
		for (int  nr=0; nr<nCSpatial; nr++) {
			tBuf_tw->DrTable[nr][subtime] =0;
		}
		//F_qq2 0  KSpace
		for (int k = 0; k < AVDOF * nCSpatial; k ++) {
			tBuf_tw->F_qq2[k][subtime] = 0.;
			if (flag_velocity == true) {
				tBuf_tw->C2_v_rho[k][subtime] = 0.;
				tBuf_tw->C2_rho_v[k][subtime] = 0.;
				tBuf_tw->C2_v_v[k][subtime]   = 0.;
			}
			if (flag_velocity && flag_magnet ) {
				tBuf_tw->C2_v_mu[k][subtime] = 0.;
				tBuf_tw->C2_mu_v[k][subtime] = 0.;
			}
			if (flag_magnet == true) {
				tBuf_tw->C2_mu_rho[k][subtime] = 0.;
				tBuf_tw->C2_rho_mu[k][subtime] = 0.;
				tBuf_tw->C2_mu_mu[k][subtime]  = 0.;
			}
			if (flag_F == true) {
				tBuf_tw->F_s_qq2[k][subtime] = 0.;
				tBuf_tw->F_d_qq2[k][subtime] = 0.;
			}
		}
	}
}
void EvalTwoTimeEach(MakeSqtClass* cl_sqt, TBuf* tBuf_tw, int subtime)
{
	if (flag_t == true) {
		Snapshot* snap = cl_sqt->snap;
		VecR3 sum_ri= {0,0,0};
		VecR3 sum_rj= {0,0,0};// t_w
		VecR3 sum_vi= {0,0,0};
		VecR3 sum_vj= {0,0,0};// t_w
		real Cvcmvcm=0, displacement_cm;
#pragma omp parallel for
		for (int n=0; n<nPtls; n++) {
			VecR3 dr;
			real dx2,dy2,dz2,dr2,Cvv,Cmm;
			Atom* col_i = &(snap->atoms[n]);
			VecR3* pos_j = &(tBuf_tw->orgR[n]);
			VecR3* vel_j = &(tBuf_tw->orgV[n]);
			VecR3* mu_j;
			if (flag_magnet == true) {
				mu_j = &(tBuf_tw->orgMu[n]);
			}

			sum_ri.x += col_i->x;
			sum_ri.y += col_i->y;
			sum_ri.z += col_i->z;

			sum_rj.x += pos_j->x;
			sum_rj.y += pos_j->y;
			sum_rj.z += pos_j->z;

			if (flag_velocity == true) {
				sum_vi.x += col_i->vx;
				sum_vi.y += col_i->vy;
				sum_vi.z += col_i->vz;

				sum_vj.x += vel_j->x;
				sum_vj.y += vel_j->y;
				sum_vj.z += vel_j->z;
			}

			dr.x =  col_i->x-pos_j->x ;
			dr.y =  col_i->y-pos_j->y ;
			dr.z =  col_i->z-pos_j->z ;

			if (flag_velocity == true) {
				Cvv = col_i->vx * vel_j->x;
				Cvv += col_i->vy * vel_j->y;
				Cvv += col_i->vz * vel_j->z;
			}

			if (flag_magnet == true) {
				Cmm = col_i->mux * mu_j->x;
				Cmm += col_i->muy * mu_j->y;
				Cmm += col_i->muz * mu_j->z;
			}

			dx2 = dr.x*dr.x; 
			dy2 = dr.y*dr.y; 
			dz2 = dr.z*dr.z; 
			dr2 = dx2 + dy2 + dz2;

			if (flag_f == true) {
				int  i_Dr    = floor (sqrt(dr2)/rVal);
				if (i_Dr<nCSpatial) tBuf_tw->DrTable[i_Dr][subtime] ++;
			}

			tBuf_tw->rrMSD[subtime] += dr2;
			tBuf_tw->rrMSR1_R[subtime].x += dx2;
			tBuf_tw->rrMSR1_R[subtime].y += dy2;
			tBuf_tw->rrMSR1_R[subtime].z += dz2;
			tBuf_tw->rrMQD[subtime] += dr2*dr2;
			if (flag_velocity == true) {
				tBuf_tw->rrCvv[subtime] += Cvv;
			}
			if (flag_magnet == true) {
				tBuf_tw->rrCmm[subtime] += Cmm;
			}

		} // for  n  in nPtls
		
		displacement_cm = (sum_ri.x - sum_rj.x)*(sum_ri.x - sum_rj.x);
		displacement_cm+= (sum_ri.y - sum_rj.y)*(sum_ri.y - sum_rj.y);
		displacement_cm+= (sum_ri.z - sum_rj.z)*(sum_ri.z - sum_rj.z);
		tBuf_tw->rrMSDCM[subtime] += displacement_cm/(nPtls*nPtls);

		if (flag_velocity == true) {
			Cvcmvcm = sum_vi.x * sum_vj.x;
			Cvcmvcm += sum_vi.y * sum_vj.y;
			Cvcmvcm += sum_vi.z * sum_vj.z;
			tBuf_tw->rrCvcmvcm[subtime] += Cvcmvcm/(nPtls*nPtls);
		}

	}
}
void EvalTwoTimeCollective(MakeSqtClass* cl_sqt, TBuf* tBuf_tw, int subtime)
{
	if (flag_s == true) {
		Rank2R3 subVR,sqVR;
		TBuf* tBuf = cl_sqt->tBuf;
		real_tensor_sub_r2_r2r2(
				&subVR, &tBuf->sumVR_ct, 
				&tBuf_tw->orgSumVR);
		real_tensor_product_r2_r2r2 (
				& sqVR, & subVR, & subVR);

		real_tensor_increase_r2_r2(
				&tBuf_tw->rrMSR2_VR[subtime], &sqVR);
	}
}

void EvalTwoTimeKSpace(MakeSqtClass* cl_sqt, TBuf* tBuf_tw, int subtime)
{
	if (flag_f == true) {
		TBuf* tBuf = cl_sqt->tBuf;
		real *	rho_q1  = tBuf->rho_q1 ;
		real *	kvel_q1  = tBuf->kvel_q1 ;
		real *	kmu_q1  = tBuf->kmu_q1 ;
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
				int markerR = marker + 2*nc;
				int markerI = marker + 2*nc +1;
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
				if ( flag_F ) {
					for (int n=0; n<nPtls; n++) {
						tBuf_tw->F_s_qq2[nav][subtime] +=
							w * (rho_s_q1[n][markerR] * tBuf_tw->org_rho_s_q1[n][markerR] +
									rho_s_q1[n][markerI] * tBuf_tw->org_rho_s_q1[n][markerI]);
						tBuf_tw->F_d_qq2[nav][subtime] +=
							w * (rho_d_q1[n][markerR] * tBuf_tw->org_rho_s_q1[n][markerR] +
									rho_d_q1[n][markerI] * tBuf_tw->org_rho_s_q1[n][markerI])+
							w * (rho_s_q1[n][markerR] * tBuf_tw->org_rho_d_q1[n][markerR] +
									rho_s_q1[n][markerI] * tBuf_tw->org_rho_d_q1[n][markerI]);
					}
				}
				if( flag_velocity == true) 
				{
					tBuf_tw->C2_v_v[nav][subtime] +=
						w * (kvel_q1[markerR] * tBuf_tw->org_kvel_q1[markerR] +
								kvel_q1[markerI] * tBuf_tw->org_kvel_q1[markerI]);

					tBuf_tw->C2_rho_v[nav][subtime] +=
						w * (kvel_q1[markerR] * tBuf_tw->org_rho_q1[markerR] +
								kvel_q1[markerI] * tBuf_tw->org_rho_q1[markerI]);

					tBuf_tw->C2_v_rho[nav][subtime] +=
						w * (rho_q1[markerR] * tBuf_tw->org_kvel_q1[markerR] +
								rho_q1[markerI] * tBuf_tw->org_kvel_q1[markerI]);
				}
				if( flag_velocity&&flag_magnet) 
				{
					tBuf_tw->C2_mu_v[nav][subtime] +=
						w * (kvel_q1[markerR] * tBuf_tw->org_kmu_q1[markerR] +
								kvel_q1[markerI] * tBuf_tw->org_kmu_q1[markerI]);
					tBuf_tw->C2_v_mu[nav][subtime] +=
						w * (kmu_q1[markerR] * tBuf_tw->org_kvel_q1[markerR] +
								kmu_q1[markerI] * tBuf_tw->org_kvel_q1[markerI]);
				}
				if( flag_magnet == true) 
				{
					tBuf_tw->C2_rho_mu[nav][subtime] +=
						w * (kmu_q1[markerR] * tBuf_tw->org_rho_q1[markerR] +
								kmu_q1[markerI] * tBuf_tw->org_rho_q1[markerI]);
					tBuf_tw->C2_mu_rho[nav][subtime] +=
						w * (rho_q1[markerR] * tBuf_tw->org_kmu_q1[markerR] +
								rho_q1[markerI] * tBuf_tw->org_kmu_q1[markerI]);
					tBuf_tw->C2_mu_mu[nav][subtime] +=
						w * (kmu_q1[markerR] * tBuf_tw->org_kmu_q1[markerR] +
								kmu_q1[markerI] * tBuf_tw->org_kmu_q1[markerI]);
/* 					printf("j = %d,rho(t) = %g+%gi,\t = %g+%gi\n",nav,
 * 							kmu_q1[markerR],
 * 							kmu_q1[markerI],
 * 							tBuf_tw->org_kmu_q1[markerR],
 * 							tBuf_tw->org_kmu_q1[markerI]
 * 							);
 */
				}
				tBuf_tw->F_qq2[nav][subtime] +=
					w * (rho_q1[markerR] * tBuf_tw->org_rho_q1[markerR] +
							rho_q1[markerI] * tBuf_tw->org_rho_q1[markerI]);
				//			}  // for nc, 
			}    // for nk , 
		} 
	} // if flag_f
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
	if( cl_sqt->nSkip < nCSkip) {
		cl_sqt->nSkip ++;
		return;
	}
	void EvalOneTimeCorr(MakeSqtClass* cl_sqt);
	TBuf* tBuf = cl_sqt->tBuf;
	Snapshot* snap = cl_sqt->snap;

	L = snap->Box.xhigh- snap->Box.xlow;
	g_Vol  = L*L*L;
	nPtls = snap->NumAtoms;
	if (cl_sqt->flag_alloc_more ==0 ) {
		//omp_set_lock(&write_lock);
		Alloc_more(cl_sqt);
		cl_sqt->flag_alloc_more =1;
		flag_global_alloc_more = 1;
		AllocMemCheck ();
		//omp_unset_lock(&write_lock);
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
	cl_sqt->nSkip =0;
}
void AllocMemCheck ()
{
	if (ErrorAllocMem == 1) {
		printf("Reserving memory Error!!!!!!\n");
		exit(1);
	}
}
void AmountAllocArray(MakeSqtClass * cl_sqt) 
{
#define AAM2 AmountAllocMem2
#define AAM AmountAllocMem

	long long int mem = 0;
	long long int part_mem = 0;
	if (flag_f == true) {
		if (flag_F == true) {
			mem += AAM2 (avF_s_qq2, AVDOF * nCSpatial, nCTime, real);
			mem += AAM2 (avF_d_qq2, AVDOF * nCSpatial, nCTime, real);
			mem += AAM2 (StdDevF_s_qq2, AVDOF * nCSpatial, nCTime, real);
			mem += AAM2 (StdDevF_d_qq2, AVDOF * nCSpatial, nCTime, real);
		}

		mem += AAM2 (avF_qq2,  AVDOF * nCSpatial, nCTime, real);
		mem += AAM2 (StdDevF_qq2,  AVDOF * nCSpatial, nCTime, real);
		if (flag_velocity == true) {
			mem += AAM2 (avC2_v_rho,  AVDOF * nCSpatial, nCTime, real);
		}
		if (flag_magnet == true) {
			mem += AAM2 (avC2_mu_rho,  AVDOF * nCSpatial, nCTime, real);
			mem += AAM2 (avC2_mu_mu,  AVDOF * nCSpatial, nCTime, real);
		}
		mem += AAM2 (valDqt,  nCSpatial, nCTime, real);
		mem += AAM2 (valGammaQT,  nCSpatial, nCTime, real);
	}

	part_mem += AAM (cl_sqt->tBuf, nCBuffer, TBuf);
	TBuf* tBuf = cl_sqt->tBuf;
	if (flag_f == true) {
		part_mem += AAM (tBuf->rho_q1, FDOF * nCSpatial, real);
		if (flag_velocity == true) {
			part_mem += AAM (tBuf->kvel_q1, FDOF * nCSpatial, real);
		}
		if (flag_magnet == true) {
			part_mem += AAM (tBuf->kmu_q1, FDOF * nCSpatial, real);
		}
		for (int nb = 0; nb < nCBuffer; nb ++) {
			TBuf* b = &tBuf[nb];
			part_mem += AAM (b->org_rho_q1, FDOF * nCSpatial, real);
			if (flag_velocity == true) {
				part_mem += AAM (b->org_kvel_q1, FDOF * nCSpatial, real);
				part_mem += AAM2 (b->C2_v_rho, AVDOF * nCSpatial, nCTime, real);
				part_mem += AAM2 (b->C2_rho_v, AVDOF * nCSpatial, nCTime, real);
				part_mem += AAM2 (b->C2_v_v, AVDOF * nCSpatial, nCTime, real);
			}
			if(flag_velocity&&flag_magnet) {
				part_mem += AAM2 (b->C2_mu_v, AVDOF * nCSpatial, nCTime, real);
				part_mem += AAM2 (b->C2_v_mu, AVDOF * nCSpatial, nCTime, real);
			}
			if (flag_magnet == true) {
				part_mem += AAM (b->org_kmu_q1, FDOF * nCSpatial, real);
				part_mem += AAM2 (b->C2_mu_rho, AVDOF * nCSpatial, nCTime, real);
				part_mem += AAM2 (b->C2_rho_mu, AVDOF * nCSpatial, nCTime, real);
				part_mem += AAM2 (b->C2_mu_mu, AVDOF * nCSpatial, nCTime, real);
			}

			if (flag_F == true) {
				part_mem += AAM2 (b->F_s_qq2, AVDOF * nCSpatial, nCTime, real);
				part_mem += AAM2 (b->F_d_qq2, AVDOF * nCSpatial, nCTime, real);
			}
			part_mem += AAM2 (b->F_qq2, AVDOF * nCSpatial, nCTime, real);
		}
	}

	if (flag_t == true) {
		mem += AAM (rrMSDAv, nCTime, real);
		mem += AAM (rrMSDCMAv, nCTime, real);
		if (flag_velocity == true) {
			mem += AAM (rrCvvAv, nCTime, real);
			mem += AAM (rrCvcmvcmAv, nCTime, real);
		}
		if (flag_magnet == true) {
			mem += AAM (rrCmmAv, nCTime, real);
		}
		mem += AAM (rrMQDAv, nCTime, real);
		mem += AAM (rrMSR1_R_Av , nCTime, VecR3);
		mem += AAM (rrDt, nCTime, real);
		// AllocArray for shear viscosity
		// 				 (diffusion of momentum)
		if (flag_s == true) {
			mem += AAM (rrMSR2_VR_Av, nCTime, Rank2R3);
			mem += AAM (rrMSR2_VR_Av_dig, nCTime, real);
			mem += AAM (rrMSR2_VR_Av_offdig, nCTime, real);
		}
	}
	if (flag_f == true) {
		mem += AAM2 (avDrTable, nCSpatial,nCTime, real);
	}


	if (flag_f == true) {
		part_mem += AAM (tBuf->rho_s_q1_temp, 
				FDOF * nCSpatial, real);

		if (flag_F== true ) {
			part_mem += AAM (tBuf->rho_s_q1, nPtls, real*);
			part_mem += AAM (tBuf->rho_d_q1, nPtls, real*);

			for (int natom=0; natom <nPtls ; natom++) {
				part_mem += AAM (tBuf->rho_s_q1[natom], FDOF * nCSpatial, real);
				part_mem += AAM (tBuf->rho_d_q1[natom], FDOF * nCSpatial, real);
			}
		}
	}
	for (int nb = 0; nb < nCBuffer; nb ++) {
		if (flag_t == true) {
			part_mem += AAM (tBuf[nb].orgR, nPtls, VecR3);
			if (flag_magnet == true) {
				part_mem += AAM (tBuf[nb].orgMu, nPtls, VecR3);
				part_mem += AAM (tBuf[nb].rrCmm, nCTime, real);
			}
			if (flag_velocity == true) {
				part_mem += AAM (tBuf[nb].orgV, nPtls, VecR3);
				part_mem += AAM (tBuf[nb].rrCvv, nCTime, real);
				part_mem += AAM (tBuf[nb].rrCvcmvcm, nCTime, real);
			}
			part_mem += AAM (tBuf[nb].rrMSD, nCTime, real);
			part_mem += AAM (tBuf[nb].rrMSDCM, nCTime, real);
			part_mem += AAM (tBuf[nb].rrMQD, nCTime, real);
			part_mem += AAM (tBuf[nb].rrMSR1_R, nCTime, VecR3);
			if (flag_s == true) {
				part_mem += AAM (tBuf[nb].rrMSR2_VR, nCTime, Rank2R3);
			}
		}

		if (flag_f == true) {
			part_mem += AAM2 (tBuf[nb].DrTable, nCSpatial,nCTime, int);
			if (flag_F) {
				part_mem += AAM (tBuf[nb].org_rho_s_q1, nPtls, real*);
				part_mem += AAM (tBuf[nb].org_rho_d_q1, nPtls, real*);

				for (int natom=0; natom <nPtls ; natom++) {
					part_mem += AAM (tBuf[nb].org_rho_s_q1[natom], FDOF * nCSpatial, real);
					part_mem += AAM (tBuf[nb].org_rho_d_q1[natom], FDOF * nCSpatial, real);
				}
			}
		}
	}
	
	mem += AAM (factorDr, nCSpatial, real);
	mem += AAM (radius, nCSpatial, real);


	
	mem += part_mem;
	unsigned int GiB,MiB,KiB,Byte;
	Byte = mem % 1024; mem /= 1024;
	KiB = mem % 1024; mem /= 1024;
	MiB = mem % 1024; 
	GiB = mem/1024;
	printf("memory Consume : ");
	if (GiB >0 ) printf("%dGiB ",GiB);
	if (MiB >0 ) printf("%dMiB ",MiB);
	if (KiB >0 ) printf("%dKiB ",KiB);
	printf("%dByte\n ",Byte);
	puts("End fake Job!!");
	

	exit(1);

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
		if (flag_f == true) {
			if (flag_F == true) {
				AllocMem2 (avF_s_qq2, AVDOF * nCSpatial, nCTime, real);
				AllocMem2 (avF_d_qq2, AVDOF * nCSpatial, nCTime, real);
				AllocMem2 (StdDevF_s_qq2, AVDOF * nCSpatial, nCTime, real);
				AllocMem2 (StdDevF_d_qq2, AVDOF * nCSpatial, nCTime, real);
			}
			AllocMem2 (avF_qq2,  AVDOF * nCSpatial, nCTime, real);
			AllocMem2 (StdDevF_qq2,  AVDOF * nCSpatial, nCTime, real);
			if (flag_velocity == true) {
				AllocMem2 (avC2_v_rho,  AVDOF * nCSpatial, nCTime, real);
				AllocMem2 (avC2_rho_v,  AVDOF * nCSpatial, nCTime, real);
				AllocMem2 (avC2_v_v,  AVDOF * nCSpatial, nCTime, real);
			}
			if (flag_velocity && flag_magnet ) {
				AllocMem2 (avC2_mu_v,  AVDOF * nCSpatial, nCTime, real);
				AllocMem2 (avC2_v_mu,  AVDOF * nCSpatial, nCTime, real);
			}
			if (flag_magnet == true) {
				AllocMem2 (avC2_mu_rho,  AVDOF * nCSpatial, nCTime, real);
				AllocMem2 (avC2_rho_mu,  AVDOF * nCSpatial, nCTime, real);
				AllocMem2 (avC2_mu_mu,  AVDOF * nCSpatial, nCTime, real);
			}
			AllocMem2 (valDqt,  nCSpatial, nCTime, real);
			AllocMem2 (valGammaQT,  nCSpatial, nCTime, real);
		}
	}
	AllocMem (cl_sqt->tBuf, nCBuffer, TBuf);
	TBuf* tBuf = cl_sqt->tBuf;
	if (flag_f == true) {
		for (nb = 0; nb < nCBuffer; nb ++) {
			TBuf* b = &tBuf[nb];
			AllocMem (b->rho_q1, FDOF * nCSpatial, real);
			AllocMem (b->org_rho_q1, FDOF * nCSpatial, real);
			AllocMem2 (b->F_qq2, AVDOF * nCSpatial, nCTime, real);
			if (flag_velocity == true) {
				AllocMem (b->kvel_q1, FDOF * nCSpatial, real);
				AllocMem (b->org_kvel_q1, FDOF * nCSpatial, real);
				AllocMem2(b->C2_v_rho, AVDOF * nCSpatial, nCTime, real);
				AllocMem2(b->C2_rho_v, AVDOF * nCSpatial, nCTime, real);
				AllocMem2(b->C2_v_v, AVDOF * nCSpatial, nCTime, real);
			}
			if(flag_velocity && flag_magnet){
				AllocMem2(b->C2_mu_v, AVDOF * nCSpatial, nCTime, real);
				AllocMem2(b->C2_v_mu, AVDOF * nCSpatial, nCTime, real);
			}
			if (flag_magnet == true) {
				AllocMem (b->kmu_q1, FDOF * nCSpatial, real);
				AllocMem (b->org_kmu_q1, FDOF * nCSpatial, real);
				AllocMem2 (b->C2_mu_mu, AVDOF * nCSpatial, nCTime, real);
				AllocMem2 (b->C2_mu_rho, AVDOF * nCSpatial, nCTime, real);
				AllocMem2 (b->C2_rho_mu, AVDOF * nCSpatial, nCTime, real);
			}

			if (flag_F == true) {
				AllocMem2 (b->F_s_qq2, AVDOF * nCSpatial, nCTime, real);
				AllocMem2 (b->F_d_qq2, AVDOF * nCSpatial, nCTime, real);
			}
		}
	}
	/*!
	 *  \brief  Memory for Green-Kubo formula
	 */
	// AllocArray for Diffuse ()
	if (flag_global_alloc ==0 ) {
		if (flag_t == true) {
			AllocMem (rrMSDAv, nCTime, real);
			AllocMem (rrMSDCMAv, nCTime, real);
			if (flag_velocity == true) {
				AllocMem (rrCvvAv, nCTime, real);
				AllocMem (rrCvcmvcmAv, nCTime, real);
			}
			if (flag_magnet == true) {
				AllocMem (rrCmmAv, nCTime, real);
			}
			AllocMem (rrMQDAv, nCTime, real);
			AllocMem (rrMSR1_R_Av , nCTime, VecR3);
			AllocMem (rrDt, nCTime, real);
			// AllocArray for shear viscosity
			// 				 (diffusion of momentum)
			if (flag_s == true) {
				AllocMem (rrMSR2_VR_Av, nCTime, Rank2R3);
				AllocMem (rrMSR2_VR_Av_dig, nCTime, real);
				AllocMem (rrMSR2_VR_Av_offdig, nCTime, real);
			}
		}
		if (flag_f == true) {
			AllocMem2 (avDrTable, nCSpatial,nCTime, real);
		}

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

	if (flag_f == true) {
		AllocMem (tBuf->rho_s_q1_temp, 
				FDOF * nCSpatial, real);

		if (flag_F== true ) {
			AllocMem (tBuf->rho_s_q1, nPtls, real*);
			AllocMem (tBuf->rho_d_q1, nPtls, real*);

			for (int natom=0; natom <nPtls ; natom++) {
				AllocMem (tBuf->rho_s_q1[natom], FDOF * nCSpatial, real);
				AllocMem (tBuf->rho_d_q1[natom], FDOF * nCSpatial, real);
			}
			fprintf(stderr, "(AllocMore)Reserving memory on heap via AllocMem : %lld GB\n",  ll_mem_size/1000ll/1000ll/1000ll);
		}
	}
	for (nb = 0; nb < nCBuffer; nb ++) {
		if (flag_t == true) {
			AllocMem (tBuf[nb].orgR, nPtls, VecR3);
			if (flag_magnet == true) {
				AllocMem (tBuf[nb].orgMu, nPtls, VecR3);
				AllocMem (tBuf[nb].rrCmm, nCTime, real);
			}
			if (flag_velocity == true) {
				AllocMem (tBuf[nb].orgV, nPtls, VecR3);
				AllocMem (tBuf[nb].rrCvv, nCTime, real);
				AllocMem (tBuf[nb].rrCvcmvcm, nCTime, real);
			}
			AllocMem (tBuf[nb].rrMSD, nCTime, real);
			AllocMem (tBuf[nb].rrMSDCM, nCTime, real);
			AllocMem (tBuf[nb].rrMQD, nCTime, real);
			AllocMem (tBuf[nb].rrMSR1_R, nCTime, VecR3);
			if (flag_s == true) {
				AllocMem (tBuf[nb].rrMSR2_VR, nCTime, Rank2R3);
			}
		}

		if (flag_f == true) {
			AllocMem2 (tBuf[nb].DrTable, nCSpatial,nCTime, int);
			if (flag_F) {
				AllocMem (tBuf[nb].org_rho_s_q1, nPtls, real*);
				AllocMem (tBuf[nb].org_rho_d_q1, nPtls, real*);

				for (int natom=0; natom <nPtls ; natom++) {
					AllocMem (tBuf[nb].org_rho_s_q1[natom], FDOF * nCSpatial, real);
					AllocMem (tBuf[nb].org_rho_d_q1[natom], FDOF * nCSpatial, real);
				}
			}
		}
	}
	fprintf(stderr, "(AllocMore)Reserving memory on heap via AllocMem : %lld GB\n",  ll_mem_size/1000ll/1000ll/1000ll);
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

	fprintf(stderr, "(AllocMore)Reserving memory on heap via AllocMem : %lld GB\n",  ll_mem_size/1000ll/1000ll/1000ll);
}

int GetNameList ()
	/*!
	 *  \brief    from book of rapaport 
	 *  	ex)	  input.file	
	 *  	Name  value 
	 *    Name2 value 
	 *  	 value type : real(double type) int 		 
	 *  			 Well defined.
	 */
{
	size_t  j, k;
	int match, ok;
	char buff[100], *token;
	FILE *fp;
	strcpy (buff, inputFilename);
	//	strcpy (buff, argv[0]);
	//	strcat (buff, ".in");
	//default value
	nCSkip=0;

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
	size_t j, k;
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
	L[0] = snap->Box.xhigh- snap->Box.xlow;
	L[1] = snap->Box.yhigh- snap->Box.ylow;
	L[2] = snap->Box.zhigh- snap->Box.zlow;

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
