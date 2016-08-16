/*
 * =====================================================================================
 *
 *       Filename:  analysis.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2016년 07월 26일 14시 24분 35초
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ph.D. Candidate KIM Hyeok (kh), ekh0324@gmail.com
 *        Company:  Konkuk University
 *
 * =====================================================================================
 */


typedef double real;
typedef struct {
	real R, I;
} Cmplx;
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "snapshot.h"
#define ALLOC(type) type*  alloc_ ## type(size_t n) { \
	  return (type *) malloc(sizeof(type)*n); \
}
#define Min(a,b) ( (a<b)?(a):(b))
#define CSet(a, x, y) a.R = x, a.I = y
#define CAdd(a, b, c) a.R = b.R + c.R, a.I = b.I + c.I
#define CSub(a, b, c) a.R = b.R - c.R, a.I = b.I - c.I
#define CMul(a, b, c) a.R = b.R * c.R - b.I * c.I, a.I = b.R * c.I + b.I * c.R
#define Sqr(x)     ( (x)*(x))
#define CHAR_MINUS '-'
#define CHAR_ZERO  '0'’
#define NameString(x)                     \
if (! strncmp (bp, #x, strlen (#x))) { \
	bp += strlen (#x);                   \
	x = strtok (bp, " ");                \
}
#define NameVal(x)                     \
if (! strncmp (bp, #x, strlen (#x))) { \
	bp += strlen (#x);                   \
	x = strtod (bp, &bp);                \
}

#define BUFF_LEN 1024
ALLOC(double); ALLOC(real); ALLOC(int); ALLOC(Cmplx);

char *header[] = {"cur-long", "cur-trans", "density", "gamma_qt", "Dqt", "vanHove-self"},
		 *txtCorr = "space-time corr";
int nDataTypes = sizeof(header)/sizeof(char*);
void PrintHelp ( char *pName);
void FftComplex (Cmplx *a, int size);
int main ( int argc, char **argv)
{
	Cmplx *work;
	real *corrSum[nDataTypes], *corrSumSq[nDataTypes], damp, deltaT, deltaTCorr,
			 omegaMax, tMax, w,x;
	int doFourier,doWindow,j,k,n,nData,nFunCorr,nSet,nSetSkip,
			nv, nValCorr;
	char *bp, *fName, buff[BUFF_LEN], *lmpFileName;
	FILE *fp;
	Snapshot* pSnap;
/* 	int nDataTypes = sizeof(header);
 * 	int tempa      = sizeof(char*);
 * 	printf( "%d %d \n ", nDataTypes, tempa);
 * 	exit(1);
 */

	n = 1;
	if (-- argc <1 || ! strcmp (argv[1], "-h")) PrintHelp (argv[0]);
	doFourier =1;
	doWindow =0;
	nSetSkip = 1;
	while (-- argc >= 0) {
		if (! strcmp (argv[n], "-t")) doFourier =0;
		else if (! strcmp (argv[n], "-w")) doWindow =1;
		else if (! strcmp (argv[n], "-s")) nSetSkip = atoi (argv[n]+2);
		else {
			fName = argv[n];
			break;
		}
		++ n;
	}

	if (argc >0) PrintHelp (argv[0]);
	omegaMax = 10.;
	tMax = 5.;
	if((fp = fopen (fName, "r")) == 0) {
		printf ("no file\n");
		exit (0);
	}

	while (1) {
		bp = fgets (buff, BUFF_LEN, fp);
		if (*bp == CHAR_MINUS) break;
		NameVal (deltaT);
		NameVal (nFunCorr);
		NameVal (nValCorr);
//		NameString (lmpFileName);
	}
	deltaTCorr =  deltaT;
	for (j = 0; j < nDataTypes; j ++) {
		corrSum[j] = alloc_real(nFunCorr * nValCorr);

		corrSumSq[j] = alloc_real( nFunCorr * nValCorr);
		for (n =0; n < nFunCorr* nValCorr; n++) {
			corrSum [j][n] = 0.;
			corrSumSq [j][n] = 0.;
		}
	}

	// The Preamble
	work = alloc_Cmplx( 2 * (nValCorr -1));
	nData =0;
	nSet =0;
	while (1) {
		//		if ( NULL == (pSnap = read_dump (fp ) )) break;
		if (! (bp = fgets (buff, BUFF_LEN, fp))) break;
		if (! strncmp (bp, txtCorr, strlen (txtCorr))) {
			++ nSet;
			if (nSet < nSetSkip) continue;
			++ nData;
			for ( j =0; j < nDataTypes; j++) {
				bp = fgets (buff, BUFF_LEN, fp); // header types check(not completed)
				for ( n =0; n<nValCorr; n ++) {
					bp = fgets (buff, BUFF_LEN, fp);

					for ( k = 0; k < nFunCorr; k += 1 ) {
						w = strtod (bp, &bp);
						corrSum[j][k * nValCorr + n] += w;
						corrSumSq[j][k * nValCorr + n] += Sqr(w);
					}
//					bp = fgets (buff, BUFF_LEN, fp);// null line
				}
			}
		}
	}
	fclose (fp);
	printf ("%d\n", nData);
 
	for ( j = 0; j < nDataTypes; j += 1 ) {
		for ( n = 0; n < nFunCorr*nValCorr; n += 1 ) {
			corrSum[j][n] /= nData;
			corrSumSq[j][n] = sqrt ( corrSumSq[j][n]/nData-Sqr(corrSum[j][n]));
		}
	}
	if (doFourier) {
		
		for ( j = 0; j < 3; j += 1 ) {
			for ( k = 0; k < nFunCorr; k += 1 ) {
				for ( n = 0; n < nValCorr; n += 1 ) {
					if (doWindow) damp = (nValCorr - n ) / (nValCorr +.5);
					else damp = 1.;
					CSet (work[n], corrSum[j][k * nValCorr +n] * damp,0.);
				}

				for ( n = nValCorr; n < 2*(nValCorr-1); n += 1 ) 
					work[n] = work[2* (nValCorr -1 ) - n ];
				FftComplex (work, 2 * nValCorr -2);

				for ( n = 0; n < nValCorr; n += 1 ) 
					corrSum[j][k * nValCorr +n] = work[n].R;
			}
		}
		omegaMax= Min(omegaMax, M_PI / deltaTCorr);
		nv = nValCorr * omegaMax / (M_PI / deltaTCorr);
	} else {
		
		for ( j = 0; j < nDataTypes; j += 1 ) {
			for ( k = 0; k < nFunCorr; k += 1 ) {
				for ( n = 1; n < nValCorr; n += 1 ) 
					corrSum[j][k * nValCorr +n] /= corrSum[j][k*nValCorr];
				corrSum[j][k * nValCorr] = 1.;
			}
		}
		tMax = Min ( tMax, (nValCorr -1) * deltaTCorr);
		nv = nValCorr * tMax / ( (nValCorr - 1) / deltaTCorr);
	}

	for ( j = 0; j < nDataTypes; j += 1 ) {
		printf("%s\n", header[j]);
		for (n=0; n < nv; n++) {
			if (doFourier) x = n * omegaMax / nv;
			else x = n * deltaTCorr;
			printf ( "%9.4f", x);
			
			for ( k = 0; k < nFunCorr; k += 1 ) 
				printf (" %9.4f", corrSum[j][k * nValCorr +n]);
			printf ("\n");
		}
	}
}

void PrintHelp ( char *pName)
{
	printf ("Usage: %s [-t{time_corr}] [-s{skip}n] [-w{window}]"
			 " input-file \n" , pName);
	exit(0);
}

void FftComplex (Cmplx *a, int size)
{
	Cmplx t, w, wo;
	real theta;
	int i, j, k, n;
	k = 0;
	for (i = 0; i < size; i ++) {
		if (i < k) {
			t = a[i];
			a[i] = a[k];
			a[k] = t;
		}
		n = size / 2;
		while (n >= 1 && k >= n) {
			k -= n;
			n /= 2;
		}
		k += n;
	}
	for (n = 1; n < size; n *= 2) {
		theta = M_PI / n;
		CSet (wo, cos (theta) - 1., sin (theta));
		CSet (w, 1., 0.);
		for (k = 0; k < n; k ++) {
			for (i = k; i < size; i += 2 * n) {
				j = i + n;
				CMul (t, w, a[j]);
				CSub (a[j], a[i], t);
				CAdd (a[i], a[i], t);
			}
			CMul (t, w, wo);
			CAdd (w, w, t);
		}
	}
}
