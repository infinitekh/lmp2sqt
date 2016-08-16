/*
 * =====================================================================================
 *
 *       Filename:  lmp2sqt.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2015년 11월 27일 16시 04분 59초
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

#include "lmp2sqt_GrKu.h"
#include "snapshot.h"



#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include	<unistd.h>
#define M_PI       3.14159265358979323846
#include<math.h>
#include"snapshot.h"

#define FORW_DIFF(A,y)  ( -A[(y)+2] +4. *A[(y)+1] -3.* A[(y) ]  )   
#define CENT_DIFF(A,y)  ( A[(y) +1] A[(y)-1]  )   
#define BACK_DIFF(A,y)  ( A[(y) -2] A[(y)-1] A[(y)] )   

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
//	NameR    (kVal),
	NameR    (deltaT),
	NameI   		(limitCorrAv),
	NameI 	(nBuffCorr), // number of simul. time seq
	NameI 	(nFunCorr),  // number of spatial seq
	NameI 	(nValCorr)   // number of time seq
};

void PrintNameList (FILE *fp);
int GetNameList (int argc, char **argv);

int main(int argc, char** argv) {
	char filename[100];
	int lowcut;
	int highcut;
	int n_snap;
	int n;
	int i,id,bin;
	int iii,jjj,kkk;
	atom* atoms, *ppi;
		FILE* fp_out;
	if(argc <2) {
		perror("#run inputfilename");
		return 1;
	}
	GetNameList(argc,argv);
	PrintNameList(stdout);

	int maxAtom =0;

	AllocArray();
	InitSpacetimeCorr();
	InitDiffuse ();
	

	strcpy( filename,argv[1]);
	FILE* fp = fopen( filename ,"r");
	Snapshot* snap;
	n_snap = 0;	
	while(1) {
		snap =	read_dump(fp);
		
		if (snap == NULL)
			break;

		EvalSpacetimeCorr(snap);
		EvalDiffuse (snap);


		
		free_Snapshot(snap);
		n_snap++;
		/* 			free(snap->atoms);
		 * 			free(snap);
		 */
	}

	if (n_snap <5){
		perror("The # of snap is too small(<5)!!");
		return 23;
	}

	return 0;
}

void AccumSpacetimeCorr ( int nCol)
{
	int j, n, nb;
	for (nb = 0; nb < nBuffCorr; nb ++) {
		if (tBuf[nb].count == nValCorr) {
			for (j = 0; j < 3 * nFunCorr; j ++) {
				for (n = 0; n < nValCorr; n ++)
					avAcfST[j][n] += tBuf[nb].acfST[j][n];
			}
			tBuf[nb].count = 0;
			++ countCorrAv;
			if (countCorrAv == limitCorrAv) {
				for (j = 0; j < 3 * nFunCorr; j ++) {
					for (n = 0; n < nValCorr; n ++)
						avAcfST[j][n] /= 3. * nCol * limitCorrAv;
				}
				PrintSpacetimeCorr (stdout);
				ZeroSpacetimeCorr ();
			}
		}
	}
}


void InitSpacetimeCorr ()
{
	int nb;
	if (nBuffCorr > nValCorr) {
		fputs("Error nBuffCorr> nValCorr", stderr);
		exit(1);
	}
		
	for (nb = 0; nb < nBuffCorr; nb ++)
		tBuf[nb].count = - nb * nValCorr / nBuffCorr;
	ZeroSpacetimeCorr ();
}
void ZeroSpacetimeCorr ()
{
	int j, n;
	countCorrAv = 0;
	for (j = 0; j < 3 * nFunCorr; j ++) {
		for (n = 0; n < nValCorr; n ++) avAcfST[j][n] = 0.;
	}
}
void EvalOtherInformation () 
{
#define Fqt_FIX_q avAcfST[3*(j) +2] 
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
{
	extern real kVal;
	int j, k, n, k2, nr;
	char *header[] = {"cur-long", "cur-trans", "density", 
		  "gamma_qt", "Dqt", "vanHove-self"};
	fprintf (fp, "space-time corr\n");
	//for (k = 0; k < 3; k ++) {
	for (k2 = 0; k2 < sizeof(header)/ sizeof(char*); k2 ++) {
		k= k2%3;
		/* 		fprintf (fp, "%s", header[k]);
		 * 		for (j = 0; j < nFunCorr; j ++)
		 * 			fprintf (fp, " %7.3f", kVal*(j+1));
		 * 		fprintf (fp, "\n");
		 */

		EvalOtherInformation ();
		fprintf (fp, "# %s %7.3f %7.3f %7.3f\n", header[k2] , kVal, 1.0*deltaT, rVal);
		switch ( k2) {
			case 0: case 1: case 2: 
				/*-----------------------------------------------------------------------------
				 *  avAcfST[3*i+k][j] -> F(q_i,t_j) k=0 longi k=1 tranv k=2 density
				 *-----------------------------------------------------------------------------*/
				for (n = 0; n < nValCorr; n ++) {
					/* 			deltaT = n *1. * deltaT;
					 * 			fprintf (fp, "%7.3f", deltaT);
					 */
					for (j = 0; j < nFunCorr; j ++)
						fprintf (fp, " %8.4e", avAcfST[3 * j + k][n]);
					fprintf (fp, "\n");
				} 
				break;
			case 3:                                   /* gamma_qt */
				for (n = 0; n < nValCorr; n ++) {   
					for (j = 0; j < nFunCorr; j ++)
						fprintf (fp, " %8.4e", valGammaQT[j][n]);
					fprintf (fp, "\n");
				}
				break ;
			case 4:                                   /* D_qt before lim t->0 */
				for (n = 0; n < nValCorr; n ++) {   
					for (j = 0; j < nFunCorr; j ++)
						fprintf (fp, " %8.4e", valDqt[j][n]);
					fprintf (fp, "\n");
				}
				break;
			case 5:                                   /* Van Hove before lim G_s(r,t) */
//				fprintf (fp, "#van Hove function\n");
	
				for ( nr=0; nr<nFunCorr; nr++)  {
					for (j = 0; j < nValCorr; j ++) {
						fprintf (fp, " %8.4e", avDrTable[nr][j] *factorDr[nr]);
					}
					fprintf (fp, "\n");
				}
				break;
		}
		fprintf (fp, "\n");
	}

	char filename1[100] ="Dq00.info" ;
	char filename2[100] ="Ft00.info" ;
	char filename3[100] ="SSF00.info" ;
	int nfile = 0;
	while( 0 == (access(filename1,F_OK))+(access(filename2,F_OK))+(access(filename3,F_OK))) {
		nfile++;
		sprintf(filename1, "Dq%02d.info",nfile);
		sprintf(filename2, "Ft%02d.info",nfile);
		sprintf(filename3, "SSF%02d.info",nfile);
	}
	FILE* fp_Dq = fopen(filename1,"w");
	FILE* fp_Ft = fopen(filename2,"w");
	FILE* fp_SSF = fopen(filename3,"w");

	fprintf (fp_Dq, "# dt = %7.3f\n", deltaT);
	fprintf (fp_Ft, "# dq = %7.3e\n", kVal);
//	fprintf (fp_SSF, "# dq = %7.3e\n", kVal);
	for (j = 0; j < nFunCorr; j ++) {
		fprintf (fp_Dq, "%8.4f" , j*kVal );
		for (n = 1; n < nValCorr; n ++) {   
			fprintf (fp_Dq, " %8.4e" ,  valDqt[j][n]);
		}
		fprintf (fp_Dq, "\n");
	}
	for (n = 0; n < nValCorr; n ++) {   
		fprintf (fp_Ft, "%8.4f" , n*deltaT );
		for (j = 0; j < nFunCorr; j ++) {
			fprintf (fp_Ft, " %8.4e" ,  avAcfST[(3*j)+2][n]/avAcfST[(3*j)+2][0]);
		}
		fprintf (fp_Ft, "\n");
	}
	
	n=0;
	for (j = 0; j < nFunCorr; j ++) {
		fprintf (fp_SSF, "%8.4f" " %8.4e""\n" , (j+1)*kVal , avAcfST[(3*j)+2][0]);
	}

	fclose(fp_Dq);
	fclose(fp_Ft);
	fclose(fp_SSF);

}

void EvalSpacetimeCorr(Snapshot* snap)
{
	real b, c, c0, c1, c2, s, s1, s2, w;
	extern real kVal;
	int j, k, m, n, nb, nc, ni, nv;
	real r[3], mu[3], L[3];
	// zero initalize current time value
	for (j = 0; j < 24 * nFunCorr; j ++) valST[j] = 0.;
// we assume L0=L1 = L2	
	L[0] = snap->box.xhigh- snap->box.xlow;
	L[1] = snap->box.yhigh- snap->box.ylow;
	L[2] = snap->box.zhigh- snap->box.zlow;
	kVal = 2.*M_PI / L[0];

	int nCol = snap->n_atoms;
	atom* col_i;
//	kVal = 2. * M_PI / nFunCorr;
//	Begin Calculatie current time information
	for (n=0; n<nCol; n++) {
		col_i = &(snap->atoms[n]);
		r[0] = col_i->x; r[0] = r[0] - L[0]* floor(r[0]/L[0])- L[0]/2.;  
		r[1] = col_i->y; r[1] = r[1] - L[1]* floor(r[1]/L[1])- L[1]/2.;  
		r[2] = col_i->z; r[2] = r[2] - L[2]* floor(r[2]/L[2])- L[2]/2.;  
		mu[0] = col_i->mux; mu[1]=col_i->muy; mu[2]=col_i->muz;
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
				valST[j ++] += mu[0] * c;
				valST[j ++] += mu[0] * s;
				valST[j ++] += mu[1] * c;
				valST[j ++] += mu[1] * s;
				valST[j ++] += mu[2] * c;
				valST[j ++] += mu[2] * s;
				valST[j ++] += c;
				valST[j ++] += s;
			}
		}
	}
	// End Calculate Current time value
	// Begin Two time corrlation function
	for (nb = 0; nb < nBuffCorr; nb ++) {
		if (tBuf[nb].count == 0) {
			for (j = 0; j < 24 * nFunCorr; j ++)
				tBuf[nb].orgST[j] = valST[j];
		}
		if (tBuf[nb].count >= 0) {
			for (j = 0; j < 3 * nFunCorr; j ++)
				tBuf[nb].acfST[j][tBuf[nb].count] = 0.;
			j = 0;
			for (k = 0; k < 3; k ++) {
				for (m = 0; m < nFunCorr; m ++) {
					for (nc = 0; nc < 4; nc ++) {
						nv = 3 * m + 2;
						if (nc < 3) {                       /*  */
							//							w = Sqr (kVal * (m + 1));
							w=1.0;
							-- nv;              // transverse  3*m +1
							if (nc == k) -- nv; // longitudinal 3*m
							else w *= 0.5;
						} else w = 1.;        // density   3*m+2
						tBuf[nb].acfST[nv][tBuf[nb].count] +=
							w * (valST[j] * tBuf[nb].orgST[j] +
									valST[j + 1] * tBuf[nb].orgST[j + 1]);
						j += 2;
					}
				}
			}
		}
		++ tBuf[nb].count;
	}
	AccumSpacetimeCorr (nCol );
}

void AllocArray ()
{
	int nb;
	AllocMem (valST, 24 * nFunCorr, real);
	AllocMem2 (avAcfST, 3 * nFunCorr, nValCorr, real);
	AllocMem2 (valDqt,  nFunCorr, nValCorr, real);
	AllocMem2 (valGammaQT,  nFunCorr, nValCorr, real);
	AllocMem (tBuf, nBuffCorr, TBuf);
	for (nb = 0; nb < nBuffCorr; nb ++) {
		AllocMem (tBuf[nb].orgST, 24 * nFunCorr, real);
		AllocMem2 (tBuf[nb].acfST, 3 * nFunCorr, nValCorr, real);
	}
	// AllocArray for Diffuse ()
	AllocMem (rrDiffuseAv, nValCorr, real);
	AllocMem2 (avDrTable, nFunCorr,nValCorr, real);
}
void Alloc_more (int  nCol) {
	int nb,nr;
	for (nb = 0; nb < nBuffCorr; nb ++) {
		AllocMem (tBuf[nb].orgR, nCol, VecR3);
		AllocMem (tBuf[nb].rTrue, nCol,VecR3);
		AllocMem (tBuf[nb].rrDiffuse, nValCorr, real);
		AllocMem2 (tBuf[nb].DrTable, nFunCorr,nValCorr, int);
	}
	AllocMem (factorDr, nFunCorr, real);
	for (nr = 0; nr < nFunCorr; nr ++) {
		if (nr ==0) factorDr[nr] =0;
		else factorDr[nr] = 1.*g_Vol/( nCol*nCol*4.*M_PI* ( nr*nr) * rVal*rVal*rVal*limitCorrAv );
	}
	
}


int GetNameList (int argc, char **argv)
{
	int id, j, k, match, ok;
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
	for (k = 0; k < sizeof (nameList) / sizeof (NameList); k ++) {
		if (nameList[k].vStatus != 1) ok = 0;
	}
	return (ok);
}
void PrintNameList (FILE *fp)
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





/*-----------------------------------------------------------------------------
 *  Diffusion part
 *
 *-----------------------------------------------------------------------------*/
void InitDiffuse () {
	int nb;
	for (nb = 0; nb < nBuffCorr; nb ++)
		tBuf[nb].countDiff = - nb * nValCorr / nBuffCorr;

	ZeroDiffuse ();

}

void ZeroDiffuse () {
	int j,nr;
	countDiffuseAv = 0;
	for (j = 0; j < nValCorr; j ++) rrDiffuseAv[j] = 0.;

	for (j = 0; j < nValCorr; j ++) 
		for (nr = 0; nr < nFunCorr; nr ++) avDrTable[nr][j]= 0.;

}

void PrintDiffuse (FILE *fp){
	real tVal;
	int j,nr;

	fprintf (fp, "#diffusion\n");

	for ( j = 0; j < nValCorr; j += 1 ) {
		tVal = j * deltaT;
		fprintf (fp, "%8.4f %8.4f\n", tVal, rrDiffuseAv[j] ); 
	}
	/*-----------------------------------------------------------------------------
	 *  van Hove function part
	 *-----------------------------------------------------------------------------*/
	fprintf (fp, "#van Hove function\n");
	
	for ( nr=0; nr<nFunCorr; nr++)  {
		for (j = 0; j < nValCorr; j ++) {
			fprintf (fp, " %8.4e", avDrTable[nr][j] *factorDr[nr]);
		}
		fprintf (fp, "\n");
	}
	fprintf (fp, "\n");
}

void EvalDiffuse (Snapshot * snap){
	VecR3 dr;
	int n, nb, ni, nr;
	real deltaR;
	int i_Dr;
	real L = snap->box.xhigh- snap->box.xlow;
	g_Vol  = L*L*L;
	int nCol = snap->n_atoms;
	static int first_run = 0;
	if (first_run ==0 ) 
		Alloc_more(nCol);
	first_run = 1;
	atom* col_i;

	for (nb = 0; nb < nBuffCorr; nb ++) {
		if (tBuf[nb].countDiff == 0) {
			for (n=0; n<nCol; n++) {
				col_i = &(snap->atoms[n]);
				tBuf[nb].orgR[n].x = col_i->x;
				tBuf[nb].orgR[n].y = col_i->y;
				tBuf[nb].orgR[n].z = col_i->z;
				tBuf[nb].rTrue[n].x = col_i->x;
				tBuf[nb].rTrue[n].y = col_i->y;
				tBuf[nb].rTrue[n].z = col_i->z;
			}
		}

		if (tBuf[nb].countDiff >= 0) {
			ni = tBuf[nb].countDiff;
			/*-----------------------------------------------------------------------------
			 *  Zero initializing
			 *-----------------------------------------------------------------------------*/
			tBuf[nb].rrDiffuse[ni]= 0.;
			for ( nr=0; nr<nFunCorr; nr++) 
				tBuf[nb].DrTable[nr][ni] =0;

			for (n=0; n<nCol; n++) {
				col_i = &(snap->atoms[n]);
				dr.x = tBuf[nb].orgR[n].x - col_i->x;
				dr.y = tBuf[nb].orgR[n].y - col_i->y;
				dr.z = tBuf[nb].orgR[n].z - col_i->z;
				deltaR = sqrt(dr.x*dr.x+dr.y*dr.y+dr.z*dr.z);
				tBuf[nb].rrDiffuse[ni] += deltaR;

				i_Dr    = (int) (deltaR/rVal);
				if (i_Dr<nFunCorr)
					tBuf[nb].DrTable[i_Dr][ni] ++;

			}
		}
		++ tBuf[nb].countDiff;
	}
	AccumDiffuse (nCol) ;
}

#define DIM 3
void AccumDiffuse (int nCol){
	real fac;
	int j,  nb, nr;
	for (nb = 0; nb < nBuffCorr; nb ++) {
		if (tBuf[nb].countDiff == nValCorr) {
			for (j = 0; j < nValCorr; j ++) {
				rrDiffuseAv[j] += tBuf[nb].rrDiffuse[j];
				for ( nr=0; nr<nFunCorr; nr++) 
					avDrTable[nr][j] += tBuf[nb].DrTable[nr][j];
			}
			tBuf[nb].countDiff = 0;
			++ countDiffuseAv;
			if (countDiffuseAv == limitCorrAv) {
				fac = 1./ ( DIM * 2 * nCol * deltaT * limitCorrAv); 
				for (j = 1; j < nValCorr; j ++)
					rrDiffuseAv[j] *= fac/ j;
				PrintDiffuse (stdout);
				ZeroDiffuse ();
			}
		}
	}

}

