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

char inputFilename[100]= "in.lmp2sqt";
void UpdateNameList ();
void PrintNameList2File (FILE *fp);
int GetNameList (int argc, char **argv);
int Number_call_Print =0;

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
	int opt_fileMax=1;
	bool files_on [argc+3];
	limitCorrAv = 0;

	for( opt_num = 1;   opt_num < argc; opt_num++)  {
		strcpy( filename,argv[opt_num]);
		FILE* fp = fopen( filename ,"r");
		if( fp == NULL) {
			files_on [opt_num] = false;
			fprintf(stderr,"Can`t open file (%s)!!\n", filename);
			continue;
		}
		Snapshot* snap, *firstSnap;

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
		else {
			if (nCBuffer < 1) {
				nCBuffer =1;
			}
			limitCorrAv += floor((n_snap - nCTime) /(nCTime/ nCBuffer));
			files_on [opt_num] = true;
		}
		fclose (fp);
	}

	PrintNameList2File(stdout);
	UpdateNameList ();


	AllocArray();
	ZeroSpacetimeCorr ();
	for( opt_num = 1;   opt_num < argc; opt_num++)  {
		if (files_on[opt_num] == false) continue;
		strcpy( filename,argv[opt_num]);
		FILE* fp = fopen( filename ,"r");
		Snapshot* snap, *firstSnap;
		n_snap = 0;	
		InitSpacetimeCorr();
		/*!
		 *  \brief  Start Calculation.
		 */
		rewind(fp);  n_snap =0;
		firstSnap = read_dump(fp);
		Init_reciprocal_space(firstSnap);
		rewind(fp);
		while(1) {
			snap =	read_dump(fp);

			if (snap == NULL)
				break;

			EvalSpacetimeCorr(snap);

			free_Snapshot(snap);
			n_snap++; 
		}
		opt_num ++;
		fclose (fp);
		fprintf(stderr, "End : file : %s\n", filename);
	}

	if ( Number_call_Print ==0 ) {
		fprintf(stderr, "limit corr  = %d, countCorrAv = %d\n",
				limitCorrAv,countCorrAv);
		FILE* output = fopen("full_results.info", "w");
		PrintSpacetimeCorr(output);
		fclose(output);

	}
	Number_call_Print ++;

	return 0;
}

void AccumSpacetimeCorr ()
	/*!
	 *  \brief  계산된 현재 시간의 SpaceTime correlation을 누적한다. 
	 */
{
	real fac;
	int j,  nb, nr, n;
	for (nb = 0; nb < nCBuffer; nb ++) {
		if (tBuf[nb].count == nCTime) {
			// check!! that  data is full
			// S(q,t), M(q,t) part
			for (j = 0; j < AVDOF * nCSpatial; j ++) {
				for (n = 0; n < nCTime; n ++){
					avF_qq2[j][n] += tBuf[nb].F_qq2[j][n];
					avF_s_qq2[j][n] += tBuf[nb].F_s_qq2[j][n];
					avF_d_qq2[j][n] += tBuf[nb].F_d_qq2[j][n];
				}
			}
			// Diffuse Part
			for (j = 0; j < nCTime; j ++) {
				rrMSDAv[j] += tBuf[nb].rrMSD[j];
				rrMQDAv[j] += tBuf[nb].rrMQD[j];
			  real_tensor_increase_r2_r2(&rrMSR2_VR_Av[j], 
					 &tBuf[nb].rrMSR2_VR[j]);
				for ( nr=0; nr<nCSpatial; nr++) {
					avDrTable[nr][j] += tBuf[nb].DrTable[nr][j];
				}
			}
			// buffer nb reset
			tBuf[nb].count = 0;
			++ countCorrAv;
			if (countCorrAv == limitCorrAv) {
				PrintSpacetimeCorr (stdout);
			}
		}
	}
}

void InitSpacetimeCorr ()
	/*!
	 *  \brief  프로그램 초기에 시간 평균을 낼 수 있도록 index를 부여하는 과정
	 */
{
	if (nCBuffer > nCTime) {
		fputs("Error nCBuffer> nCTime\n", stderr);
		exit(1);
	}

	for (int nb = 0; nb < nCBuffer; nb ++){
		tBuf[nb].count = - nb * nCTime / nCBuffer;
		tBuf[nb].countDiff = - nb * nCTime / nCBuffer;
	}
	ZeroSpacetimeCorr ();
}

void ZeroSpacetimeCorr ()
	/*!
	 *  \brief  출력 후 또는, 프로그램 시작시 평균 계산을 위한 메모리를 
	 *  			0값으로 초기화
	 */
{
	int j, n, nr,d1,d2;
	countCorrAv = 0;
	for (j = 0; j < AVDOF * nCSpatial; j ++) {
		for (n = 0; n < nCTime; n ++) {
			avF_qq2[j][n] = 0.;
			avF_s_qq2[j][n] = 0.;
			avF_d_qq2[j][n] = 0.;
		}
	}

	for (j = 0; j < nCTime; j ++) {
		rrMSDAv[j] = 0.;
		rrMQDAv[j] = 0.;
		real_tensor_zero_r2(	&rrMSR2_VR_Av [j] ) ;
		rrMSR2_VR_Av_offdig[j] = 0.;
		rrMSR2_VR_Av_dig [j]   = 0.;
	}
	for (j = 0; j < nCTime; j ++) 
		for (nr = 0; nr < nCSpatial; nr ++) avDrTable[nr][j]= 0.;
}

void EvalOtherInformation () 
	/*!
	 *  \brief  \f$ F(q,t) \f$를 출력전에 미분해서 data를 뽑아낸다. 
	 *
	 */
{                            // this evaluation yield analysis.c
#define Fqt_FIX_q avF_qq2[AVDOF*(j) +AV_DEN] 
	int j,  n,  ppT, pT, cT, nT, nnT;
	extern real kVal;
	real kVal2 = kVal*kVal;
	n=0; nnT = n+2; nT = n+1; cT = n; {  //Forward O(h^2)
		for (j = 0; j < nCSpatial; j ++) {
			valGammaQT[j][n]=		 (-(Fqt_FIX_q[nnT]) +4.*(Fqt_FIX_q[nT]) -3.*(Fqt_FIX_q[cT]) )/ (2.0* deltaT*Fqt_FIX_q[cT]);
			valDqt [j][n] = - valGammaQT[j][n] / (kVal2*j*j) ;
		}
	}
	for (n = 1; n < nCTime-1; n ++) {     /* centerd O(h^2) */
		pT = n-1; nT = n+1; cT = n;
		for (j = 0; j < nCSpatial; j ++) {
			valGammaQT[j][n] = ( (Fqt_FIX_q[nT]) -(Fqt_FIX_q[pT]) )/ (2.0* deltaT*Fqt_FIX_q[cT]);
			valDqt [j][n] = - valGammaQT[j][n] / (kVal2*j*j) ;
		}
	}
	n= nCTime-1; ppT = n-2; pT = n-1; cT = n; { /* Backward O(h^2) */
		for (j = 0; j < nCSpatial; j ++) {
			valGammaQT[j][n] = (+(3.*Fqt_FIX_q[cT]) -4.*(Fqt_FIX_q[pT]) +(Fqt_FIX_q[ppT]) )/ (2.0* deltaT*Fqt_FIX_q[cT]);
			valDqt [j][n] = - valGammaQT[j][n] / (kVal2*j*j) ;
		}
	}
}

void prePrintProcess () 
{
	real scale_factor = 1./(3.*nPtls*countCorrAv);
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

	for (int nt = 1; nt < nCTime; nt ++) {
		rrMSDAv[nt] *= scale_factor;
		rrMQDAv[nt] *= scale_factor;
		
		/*!
		 *  \brief  if all mass is same value
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
	real tVal;
	int j, nType, n, k2, nr;
	Number_call_Print ++;
	//	char *header[] = {"cur-long", "cur-trans", "density", "vanHove-self"};
	char *header[] = {"cur-long", "cur-trans", "mag-long", "mag-trans","density", "vanHove-self",
		"self-long", "self-trans", "self",
		"cross-long", "cross-trans", "cross"};
	fprintf (fp, "space-time corr\n");
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
			case 0: case 1:  case 2: case 3: case 4:
				/*!-----------------------------------------------------------------------------
				 *  avF_qq2[AVDOF*i+nType][j] -> F(q_i,t_j) 
				 *  nType=0 longi  - cur
				 *  nType=1 tranv  - cur
				 *  nType=2 longi  - magnetic
				 *  nType=3 tranv  - maagnetic
				 *  nType=4 density - density correlator
				 *-----------------------------------------------------------------------------*/
				nType= k2;
				for (int nt = 0; nt < nCTime; nt ++) {
					/* 			deltaT = n *1. * deltaT;
					 * 			fprintf (fp, "%7.3f", deltaT);
					 */
					for (int nr = 0; nr < nCSpatial; nr ++){
						fprintf (fp, " %8.4e", avF_qq2[AVDOF * nr + nType][nt]);
					}
					fprintf (fp, "\n");
				} 
				break;
			case 5:              
				//        fprintf (fp, "#van Hove function\n");

				for (int nt = 0; nt < nCTime; nt ++) {
					for ( nr=0; nr<nCSpatial; nr++)  {
						fprintf (fp, " %8.4e", avDrTable[nr][nt] );
					}
					fprintf (fp, "\n");
				}
				break;
			case 6: case 7: case 8: 
				/*-----------------------------------------------------------------------------
				 *  avF_s_qq2[3*i+nType][j] -> F_s(q_i,t_j) k=0 longi k=1 tranv k=2 density
				 *-----------------------------------------------------------------------------*/
				nType= (k2-6)%3;
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
			case 9: case 10: case 11: 
				/*-----------------------------------------------------------------------------
				 *  avF_s_qq2[3*i+nType][j] -> F_s(q_i,t_j) k=0 longi k=1 tranv k=2 density
				 *-----------------------------------------------------------------------------*/
				nType= (k2-9)%3;
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
		}
		fprintf (fp, "\n");
	}
	void PrintEtc();
	PrintEtc ();
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
				avF_qq2[(AVDOF*nr)+AV_DEN][0]);
	}
	fclose(fp_SSF);

	//  fprintf (fp_SSF, "# dq = %7.3e\n", kVal);




	FILE* fp_Dt = fopen(filename1,"w");
	fprintf (fp_Dt, "#time MSD diffusion MeanQuad-Disp MSVR_dig MSVR_offdig\n");

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
		fprintf (fp_Dt, "%8.4f %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e \n", 
				tVal, rrMSDAv[nt] , rrDt[nt], rrMQDAv[nt],
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

	ZeroSpacetimeCorr ();
}

void ZeroOneTimeCorr(Snapshot* snap)
{
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
	real_tensor_zero_r2(&sumVR_ct);
}
void EvalOneTimeSumVR(Snapshot* snap) 
{
	Rank2R3 VR;  
	VecR3 dr, vecr3,vel;
	for (int n=0; n<nPtls; n++) {
		atom* col_i;

		col_i = &(snap->atoms[n]);

		vecr3.x = col_i->x;
		vecr3.y = col_i->y;
		vecr3.z = col_i->z;
		vel.x = col_i->vx;
		vel.y = col_i->vy;
		vel.z = col_i->vz;

		real_tensor_product_r2_r1r1 (& VR, &vel, &vecr3);
		real_tensor_increase_r2_r2(&sumVR_ct, &VR);
	}
}

void EvalOneTimeKspace(Snapshot* snap)
{
	real r[3], v[3],mu[3];

	/*-----------------------------------------------------------------------------
	 *  Direct calculate  rho(q)
	 *-----------------------------------------------------------------------------*/
	for (int n=0; n<nPtls; n++) {
		atom* col_i;
		col_i = &(snap->atoms[n]);
		/* 		r[0] = col_i->x; r[0] = r[0] - L* floor(r[0]/L)- L/2.;  
		 * 		r[1] = col_i->y; r[1] = r[1] - L* floor(r[1]/L)- L/2.;  
		 * 		r[2] = col_i->z; r[2] = r[2] - L* floor(r[2]/L)- L/2.;  
		 */
		r[0]  =  col_i->x;   r[1] = col_i->y;   r[2]  = col_i->z; 
		v[0]  =  col_i->vx;  v[1] = col_i->vy;  v[2]  = col_i->vz;
		mu[0] = col_i->mux; mu[1] = col_i->muy; mu[2] = col_i->muz;
		int j = 0;
		real b,c,s,c0,c1,s1,c2,s2;
		// 
		for (int k = 0; k < DIM; k ++) {          
			for (int m = 0; m < nCSpatial; m ++) {  
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
				rho_s_q1_temp[j ++] = v[0] * c;
				rho_s_q1_temp[j ++] = v[0] * s;
				rho_s_q1_temp[j ++] = v[1] * c;
				rho_s_q1_temp[j ++] = v[1] * s;
				rho_s_q1_temp[j ++] = v[2] * c;
				rho_s_q1_temp[j ++] = v[2] * s;
				rho_s_q1_temp[j ++] = mu[0] * c;
				rho_s_q1_temp[j ++] = mu[0] * s;
				rho_s_q1_temp[j ++] = mu[1] * c;
				rho_s_q1_temp[j ++] = mu[1] * s;
				rho_s_q1_temp[j ++] = mu[2] * c;
				rho_s_q1_temp[j ++] = mu[2] * s;
				rho_s_q1_temp[j ++] = c;
				rho_s_q1_temp[j ++] = s;
			}
		}
		//			memcpy(rho_s_q1, rho_q1,sizeof(real)*24*nCSpatial);
		for(j=0; j< FDOF * nCSpatial; j++ ) {
			rho_q1 [ j] += rho_s_q1_temp[j];
		}
		if ( flagSelf ) {
			memcpy(rho_s_q1[n], rho_s_q1_temp, sizeof(real)*FDOF*nCSpatial);
		}
	} /* for loop : n<nPtls */

	if ( flagSelf ) {
		for(int j=0; j< FDOF * nCSpatial; j++ ) {
			for (int n=0; n<nPtls; n++) {
				rho_d_q1[n] [ j] = rho_q1[j] - rho_s_q1 [n][j];
			}
		}
	}
}

void EvalOneTimeCorr(Snapshot* snap)
	/*!
	 * 
	 *  \brief  one time Correlation을 계산한다. 
	 * q space value는 x, y, z 방향 세개의 방향으로 
	 * longitudinal version, translational version과 가장 기본적인 방향성분 없는 density
	 */
{
	void ZeroOneTimeCorr(Snapshot* snap);
	void EvalOneTimeSumVR(Snapshot* snap) ;
	void EvalOneTimeKspace(Snapshot* snap);

	ZeroOneTimeCorr(snap);

	EvalOneTimeSumVR(snap);
	EvalOneTimeKspace(snap);
}

void SetWaitedTimeCorr(Snapshot* snap, TBuf* tBuf_tw)
{
	/*-----------------------------------------------
	 *   t_w information 
	 *-----------------------------------------------*/
	real_tensor_copy_r2r2(& tBuf_tw->orgSumVR, &sumVR_ct);
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
void InitTwoTimeCorr (Snapshot* snap, TBuf* tBuf_tw, int subtime)
{
	/*------------------------------
	 *  Zero initializing
	 *-----------------------------*/
	tBuf_tw->rrMSD[subtime]= 0.;
	tBuf_tw->rrMQD[subtime]= 0.;
	real_tensor_zero_r2 (&tBuf_tw->rrMSR2_VR[subtime]);

	for (int  nr=0; nr<nCSpatial; nr++) {
		tBuf_tw->DrTable[nr][subtime] =0;
	}

			//F_qq2 0  KSpace
	for (int j = 0; j < AVDOF * nCSpatial; j ++) {
		tBuf_tw->F_qq2[j][subtime] = 0.;
		tBuf_tw->F_s_qq2[j][subtime] = 0.;
		tBuf_tw->F_d_qq2[j][subtime] = 0.;
	}
}
void EvalTwoTimeEach(Snapshot* snap, TBuf* tBuf_tw, int subtime)
{
	for (int n=0; n<nPtls; n++) {
		VecR3 dr, vecr3,vel;
		atom* col_i = &(snap->atoms[n]);

		dr.x =  col_i->x-tBuf_tw->orgR[n].x ;
		dr.y =  col_i->y-tBuf_tw->orgR[n].y ;
		dr.z =  col_i->z-tBuf_tw->orgR[n].z ;

		real deltaR2 = (dr.x*dr.x+dr.y*dr.y+dr.z*dr.z);

		int  i_Dr    = floor (sqrt(deltaR2)/rVal);
		if (i_Dr<nCSpatial) tBuf_tw->DrTable[i_Dr][subtime] ++;

		tBuf_tw->rrMSD[subtime] += deltaR2;
		tBuf_tw->rrMQD[subtime] += deltaR2*deltaR2;
	}
}
void EvalTwoTimeCollective(Snapshot* snap, TBuf* tBuf_tw, int subtime)
{
	real_tensor_sub_r2_r2r2(&subVR, &sumVR_ct, &tBuf_tw->orgSumVR);
	real_tensor_product_r2_r2r2 (& sqVR, & subVR, & subVR);

	real_tensor_increase_r2_r2(&tBuf_tw->rrMSR2_VR[subtime],
			&sqVR);

}

void EvalTwoTimeKSpace(Snapshot* snap, TBuf* tBuf_tw, int subtime)
{
	real w;
	int nv;
	for (int j=0,k = 0; k < DIM; k ++) { // 3 loop
		for (int m = 0; m < nCSpatial; m ++) {
			const int diffMarker = m*FDOF;
			for (int nc = 0; nc < 7; nc ++) {  
				//-----------------------------------------------
				//   n_c = 0 1 2 vx vy vz   3 4 5 mx my mz 6 density
				//-----------------------------------------------
				//-----------------------------------------------
				//   n_v = 0 1  v_long v_trans    2 3 m_long m_trans 4  density
				//--------------------------------------------
				if (nc < 3) {    /*   */
					int axis = nc;
					if (axis == k) {
						w = 1.0;
						nv = diffMarker +V_LONG ;
					}
					else {
						w = 0.5;    //   
						nv = diffMarker +V_TRANS ;
					}
					//              else w *= 0.5;
				}
				else if (nc<6) {
					int axis = nc -3;
					//              w = Sqr (kVal * (m + 1));
					if (axis == k) {
						w = 1.0;
						nv = diffMarker + M_LONG;
					}
					else {
						w = 0.5;    //   
						nv = diffMarker +M_TRANS ;
					}
					//              else w *= 0.5;
				}
				else {
					w = 1.;  
					nv = diffMarker + AV_DEN;
				};   // density   3*m+4
				// cos(q*r(t)) cos(q*r(t_w) +sin sin
				if (flagSelf ) {
					for (int n=0; n<nPtls; n++) {
						tBuf_tw->F_s_qq2[nv][subtime] +=
							w * (rho_s_q1[n][j] * tBuf_tw->org_rho_s_q1[n][j] +
									rho_s_q1[n][j + 1] * tBuf_tw->org_rho_s_q1[n][j + 1]);
						tBuf_tw->F_d_qq2[nv][subtime] +=
							w * (rho_d_q1[n][j] * tBuf_tw->org_rho_s_q1[n][j] +
									rho_d_q1[n][j + 1] * tBuf_tw->org_rho_s_q1[n][j + 1])+
							w * (rho_s_q1[n][j] * tBuf_tw->org_rho_d_q1[n][j] +
									rho_s_q1[n][j + 1] * tBuf_tw->org_rho_d_q1[n][j + 1]);
					}
				}

				tBuf_tw->F_qq2[nv][subtime] +=
					w * (rho_q1[j] * tBuf_tw->org_rho_q1[j] +
							rho_q1[j + 1] * tBuf_tw->org_rho_q1[j + 1]);
				j += 2;
			}  // for nc, total j+=8
			assert ( j% DOF ==0);
		}    // for m , total j+=8*nCSpatial
		assert (j%(DOF*nCSpatial)==0);
	} 
}
void EvalTwoTimeCorr(Snapshot* snap, TBuf* tBuf_tw, int subtime)
{
	InitTwoTimeCorr(snap, tBuf_tw, subtime);
	EvalTwoTimeEach(snap, tBuf_tw, subtime);
	EvalTwoTimeCollective(snap, tBuf_tw, subtime);
	EvalTwoTimeKSpace(snap, tBuf_tw, subtime);

}
void EvalSpacetimeCorr(Snapshot* snap)
	/*!
	 * 
	 *  \brief  space time correlation을 계산한다. 
	 * q space value는 x, y, z 방향 세개의 방향으로 
	 * longitudinal version, translational version과 가장 기본적인 방향성분 없는 density
	 * version 3개를 구함.
	 *기계 편의적으로 코드가 짜져 있어서 사람이 보기에 별로 직관적이지 못해서 고칠 예정이고 
	 *거기다가 참조한 기본 코드에서 magnetization에 대한 version으로 바꾸면서 많이 복잡해지고 좋지 
	 *않아짐. 그리고 개인적으로 속도 vector도 다시 정보로 가져올 것이기 때문에 바뀔 야정 
	 *
	 * PREV. $M_T(q,t)$ $M_L(q,t)$ 
	 * Todo. 
	 *  \param  Snapshot* Snapshot 포인터 
	 */
{
	extern real kVal;
	void EvalOneTimeCorr(Snapshot* snap);

	L = snap->box.xhigh- snap->box.xlow;
	g_Vol  = L*L*L;
	nPtls = snap->n_atoms;
	static int first_run = 0;
	if (first_run ==0 ) {
		Alloc_more();
		first_run++;
	}

	kVal = 2.*M_PI / L;

	EvalOneTimeCorr(snap);

	// End Calculate Current time value
	// Begin Two time corrlation function
	for (int nb = 0; nb < nCBuffer; nb ++) {
		if (tBuf[nb].count == 0) {
			SetWaitedTimeCorr(snap, &tBuf[nb]);
		}     // End   buffer count ==0

		if (tBuf[nb].count >= 0) {
			EvalTwoTimeCorr(snap,&tBuf[nb],tBuf[nb].count);
		}                        // End buffer count >=0
		++ tBuf[nb].count;
	}
	AccumSpacetimeCorr ();
}
void AllocMemCheck ()
{
	if (ErrorAllocMem == 1) {
		printf("Reserving memory Error!!!!!!\n");
		exit(1);
	}
}
void AllocArray ()
	/*!
	 *  \brief   이름그대로 memory 할다함. 
	 *     				rho_q1 functions of q      
	 *     	등.
	 */
{
	int nb, natom;


	AllocMem (rho_q1, FDOF * nCSpatial, real);

	AllocMem2 (avF_s_qq2, AVDOF * nCSpatial, nCTime, real);
	AllocMem2 (avF_d_qq2, AVDOF * nCSpatial, nCTime, real);
	AllocMem2 (avF_qq2,  AVDOF * nCSpatial, nCTime, real);

	AllocMem2 (valDqt,  nCSpatial, nCTime, real);
	AllocMem2 (valGammaQT,  nCSpatial, nCTime, real);
	AllocMem (tBuf, nCBuffer, TBuf);
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
	AllocMem (rrMSDAv, nCTime, real);
	AllocMem (rrMQDAv, nCTime, real);
	// AllocArray for shear viscosity
	// 				 (diffusion of momentum)
	AllocMem (rrMSR2_VR_Av, nCTime, Rank2R3);
	AllocMem (rrMSR2_VR_Av_dig, nCTime, real);
	AllocMem (rrMSR2_VR_Av_offdig, nCTime, real);
	AllocMem (rrDt, nCTime, real);
	AllocMem2 (avDrTable, nCSpatial,nCTime, real);

	fprintf(stderr, "Reserving memory on heap via AllocMem : %d mb\n", (int) ll_mem_size/1000/1000);
}
void Alloc_more () {
	/*!
	 *  \brief  Alloc_more 
	 *          Allocing   using nPtls  is post-process
	 *
	 */
	int nb,nr; real rho0, shell_Vol;
	AllocMem (rho_s_q1_temp, FDOF * nCSpatial, real);

	if (flagSelf ) {
		AllocMem (rho_s_q1, nPtls, real*);
		AllocMem (rho_d_q1, nPtls, real*);

		for (int natom=0; natom <nPtls ; natom++) {
			AllocMem (rho_s_q1[natom], FDOF * nCSpatial, real);
			AllocMem (rho_d_q1[natom], FDOF * nCSpatial, real);
		}
		fprintf(stderr, "Reserving memory on heap via AllocMem : %d mb\n", (int) ll_mem_size/1000/1000);
	}
	for (nb = 0; nb < nCBuffer; nb ++) {
		AllocMem (tBuf[nb].orgR, nPtls, VecR3);
		AllocMem (tBuf[nb].rrMSD, nCTime, real);
		AllocMem (tBuf[nb].rrMQD, nCTime, real);
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
	fprintf(stderr, "Reserving memory on heap via AllocMem : %d mb\n", (int) ll_mem_size/1000/1000);
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
	fprintf(stderr, "Reserving memory on heap via AllocMem : %d mb\n", (int) ll_mem_size/1000/1000);
}

int GetNameList (int argc, char **argv)
	/*!
	 *  \brief    Rapaport 책에서 가져온 코드로, 
	 *  						Name  value 
	 *  						Name2 value 
	 *  				형식으로 되어 있는 초기값을 불러오는데 사용됨. 
	 *  				지원하는 형태는 integer와 real 값을 받아다가 초기값 배정함. 
	 *  				잘 작동함. 
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

	n_mul = kVal/ (2.*M_PI/ L[0] );
	if (n_mul <=0) n_mul =1;
	new_dk = (2.*M_PI/L[0]) * n_mul;
	fprintf(stderr, "Update for input dk param: %f -> %f \n"
			, kVal, new_dk);
	kVal = new_dk;
}
