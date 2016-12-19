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


#include "common.h"
typedef struct {
	real R, I;
} Cmplx;
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "snapshot.h"
#include <errno.h>
#include <getopt.h>
#include <unistd.h>


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

#define BUFF_LEN 4096
ALLOC(double); ALLOC(real); ALLOC(int); ALLOC(Cmplx);
int  type[]    = { 1,   1,  1,  0 };
int  scail[]    = { 1,   1,  1,  0 };
char *header[] = {"cur-long", "cur-trans", "density", "vanHove-self"},
		 *txtCorr = "space-time corr";
int nDataTypes = sizeof(header)/sizeof(char*);
void PrintHelp ( char *pName);
void FftComplex (Cmplx *a, int size);


static int verbose_flag;
int main ( int argc, char **argv)
{
	Cmplx *work;
	real *corrSum[nDataTypes], *corrSumSq[nDataTypes], *corrSumErr[nDataTypes],
	damp, deltaT, deltaTCorr,
	omegaMax, tMax, w,x,  kVal, kVal2, qVal, qVal2;
	real valGamma, valDq, valSq;
	int doFourier,doWindow;
	int j,k,n,nT,nnT,nnnT,cT, pT,ppT, pppT;
	int nData,nFunCorr,nSet,nSetSkip,
			nv, nValCorr, NValDiff;
	char *bp, *fName, buff[BUFF_LEN], *lmpFileName, output_filename[BUFF_LEN];
	FILE *input, *output;
	Snapshot* pSnap;
	/* 	int nDataTypes = sizeof(header);
	 * 	int tempa      = sizeof(char*);
	 * 	printf( "%d %d \n ", nDataTypes, tempa);
	 * 	exit(1);
	 */


	/*-----------------------------------------------------------------------------
	 *  Argument Check!! And open a file.
	 *-----------------------------------------------------------------------------*/
	n = 1;
	if (-- argc <1 || ! strcmp (argv[1], "-h")) PrintHelp (argv[0]);
	doFourier =0;
	doWindow =0;
	nSetSkip = 1;
	/* 	while (-- argc >= 0) {
	 * 		if (! strcmp (argv[n], "-t")) doFourier =0;
	 * 		else if (! strcmp (argv[n], "-w")) doWindow =1;
	 * 		else if (! strcmp (argv[n], "-s")) nSetSkip = atoi (argv[n]+2);
	 * 		else {
	 * 			fName = argv[n];
	 * 			break;
	 * 		}
	 * 		++ n;
	 * 	}
	 */
	int   param_opt;
	opterr   = 0;
	NValDiff = 4;
	static struct option long_options[] = 
	{
		{"verbose", no_argument,       &verbose_flag, 1},
		{"time", no_argument, 0, 't'},
		{"fourier", no_argument, 0, 'f'},
		{"window", no_argument, 0, 'w'},
		{"nskip", required_argument, 0, 's'},
		{"ndiff", required_argument, 0, 'd'},
		{"help", no_argument, 0, 'h'},
		{0,0,0,0}
	};
	int option_index=0;

	while( -1 !=( param_opt = getopt_long( argc, argv, "ftws:d:h", long_options, &option_index)))
	{
		switch( param_opt)
		{
			case  't'   :  
				doFourier =0;
				break;
			case  'f'   :  
				doFourier =1;
				break;
			case  'w'   :  
				doWindow =1;
				break;
			case  's'   :  
				nSetSkip = atoi(optarg);
				break;
			case  'd'   :  
				NValDiff = atoi(optarg);
				break;
			case  '?'   :  
				printf( "알 수 없는 옵션: %cn", optopt);
			case 'h'    :
				PrintHelp(argv[0]);
				break;
		}
	}
	if (optind >0) PrintHelp (argv[0]);
	//	fName = argv[1];
	fName = argv[optind];

	omegaMax = 10.;
	tMax = 100.;

	if(!strcmp(fName,"-")) {
		input = stdin;
	} else {
		input = fopen(fName,"r");
		if (NULL == input) {
			fprintf(stderr, "Unable to open '%s': %s\n",
					fName, strerror(errno));
			exit(EXIT_FAILURE);
		}
	}
	//	if((input = fopen (fName, "r")) == 0) {
	//		printf ("no file : %s\n", fName);
	//		fp
	//		exit (0);
	//	}

	while (1) {
		bp = fgets (buff, BUFF_LEN, input);
		if (*bp == CHAR_MINUS) break;
		NameVal (deltaT);
		NameVal (nFunCorr);
		NameVal (nValCorr);
		NameVal (kVal);
		//		NameString (lmpFileName);
	}
	deltaTCorr =  deltaT;
	kVal2 = kVal*kVal;

	/*-----------------------------------------------------------------------------
	 *   Alloc memory and initialization
	 *-----------------------------------------------------------------------------*/
	for (j = 0; j < nDataTypes; j ++) {
		corrSum[j] = alloc_real(nFunCorr * nValCorr);

		corrSumSq[j] = alloc_real( nFunCorr * nValCorr);
		corrSumErr[j] = alloc_real( nFunCorr * nValCorr);
		for (n =0; n < nFunCorr* nValCorr; n++) {
			corrSum [j][n] = 0.;
			corrSumSq [j][n] = 0.;
			corrSumErr [j][n] = 0.;
		}
	}

	// The Preamble
	if (doFourier)
		work = alloc_Cmplx( 2 * (nValCorr -1));
	nData =0;
	nSet =0;
	while (1) {
		//		if ( NULL == (pSnap = read_dump (input ) )) break;
		if (! (bp = fgets (buff, BUFF_LEN, input))) break;
		if (! strncmp (bp, txtCorr, strlen (txtCorr))) {
			++ nSet;
			if (nSet < nSetSkip) continue;
			++ nData;
			for ( j =0; j < nDataTypes; j++) {
				bp = fgets (buff, BUFF_LEN, input); // header types check(not completed)
				printf("## header check : %s \n", buff);

				for ( n =0; n<nValCorr; n ++) {
					bp = fgets (buff, BUFF_LEN, input);
#ifndef NDEBUG
					puts(buff);
#endif
					for ( k = 0; k < nFunCorr; k += 1 ) {
						w = strtod (bp, &bp);
#ifndef NDEBUG
						printf(" %8.4f",w);
#endif
						corrSum[j][k * nValCorr + n] += w;
						corrSumSq[j][k * nValCorr + n] += Sqr(w);
					}
#ifndef NDEBUG
					puts("\n");
#endif
				}
				bp = fgets (buff, BUFF_LEN, input);// null line
				printf("## endline null string default : %s \n", buff);
			}
		}
	}
	fclose (input);
	printf ("%d\n", nData);

	for ( j = 0; j < nDataTypes; j += 1 ) {
		for ( n = 0; n < nFunCorr*nValCorr; n += 1 ) {
			corrSum[j][n] /= nData;
			corrSumSq[j][n] = sqrt ( corrSumSq[j][n]/nData-Sqr(corrSum[j][n]));
			corrSumErr[j][n] =  corrSumSq[j][n]/ sqrt(nData);
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
	}   //  print A(q,omega) or A(r,omega)
	else {   // print A(q,t), or A(r,t)
		for ( j = 0; j < nDataTypes; j += 1 ) {
			if (scail[j] ){
				/*-----------------------------------------------------------------------------
				 *  D(q) ~ S(q) ~ or density and  spin longi,transverse
				 *-----------------------------------------------------------------------------*/
				sprintf( output_filename, "%s00.out", header[j] );
				FILE *output = fopen( output_filename, "w");
				fprintf(output,"#kVal Sq Gamma Dq\n"
						"#deltaT = %9.4f\n", deltaT);
				int nValDiffMin=2;
				for (k = 0; k < nFunCorr; k ++) {
					qVal = (k+1)*kVal; qVal2=qVal*qVal;

					//					cT= k*nValCorr+NValDiff; 
					//					for (cT = k*nValCorr+nValDiffMin; cT < (k+1)*nValCorr-2; cT++) {
					//					for (cT = k*nValCorr+nValDiffMin; cT < (k)*nValCorr+6; cT++) {
					//						if ( corrSum[j][cT] < 0.2) break;
					//					}
					cT = k*nValCorr+NValDiff;
					nnT = cT+2; nT = cT+1; nnnT=cT+3;   //Forward first Derivative
					ppT = cT-2; pT = cT-1; pppT=cT-3;   //central first Derivative
					valSq    =   corrSum[j][k*nValCorr];
#define Xm3 (log(corrSum[j][pppT]))
#define Xm2 (log(corrSum[j][ppT]))
#define Xm1 (log(corrSum[j][pT]) )
#define Xp3 (log(corrSum[j][nnnT]))
#define Xp2 (log(corrSum[j][nnT]))
#define Xp1 (log(corrSum[j][nT]) )
#define X   (log(corrSum[j][cT]) )
					//					valGamma =	 (-(corrSum[j][nnT]) +4.*(corrSum[j][nT]) -3.*(corrSum[j][cT]) )/ (2.0* deltaT*corrSum[j][cT]);
					//					valGamma =   ( (-Xp2+4.*Xp1 -3.* X ) / (2.*deltaT )) ;  // Forward O(h^2) first Derivative
					//					valGamma =   ( (Xp1 - X ) / deltaT ) ;  // Forward O(h) first Derivative
					//					valGamma = (-11.*X + 18.*Xp1-9.*Xp2+2.*Xp3)/ ( 6.*deltaT);// F O(h^3)
					valGamma =  (-Xp2+ 8.*Xp1 - 8.*Xm1 + Xm2)/ (12.*deltaT); // Cetral O(h^4)
					//					valGamma  =    
					valDq    = - valGamma / (qVal2) ;
					fprintf( output, "%8.4e %8.4e %8.4e %8.4e %8.4e\n", qVal, valSq,valGamma,valDq, corrSumErr[j][cT]); 
				}
				fclose(output);


				//      scaling by function of  t=0
				for ( k = 0; k < nFunCorr; k += 1 ) {
					for ( n = 1; n < nValCorr; n += 1 ) 
						corrSum[j][k * nValCorr +n] /= corrSum[j][k*nValCorr];
					corrSum[j][k * nValCorr] = 1.;
				}
				}
			}
			tMax = Min ( tMax, (nValCorr -1) * deltaTCorr);
			nv = nValCorr * tMax /  (nValCorr - 1) / deltaTCorr;
			printf("nv = %d, tMax = %f, nValCorr = %d, deltaTCorr = %f\n" , nv ,tMax, nValCorr, deltaTCorr);
		} // else end

		for ( j = 0; j < nDataTypes; j += 1 ) {
			printf("%s\n", header[j]);
			for (n=0; n < nv; n++) {
				if (doFourier) x = n * omegaMax / nv;
				else x = n * deltaTCorr;
				printf ( "%9.4e", x);

				for ( k = 0; k < nFunCorr; k += 1 ) 
					printf (" %9.4e", corrSum[j][k * nValCorr +n]);
				printf ("\n");
			}
		}
	}

	void PrintHelp ( char *pName)
	{
		printf ("Usage: %s [-t{time_corr}] [-s{skip}n] [-w{window}]"
				" input-file \n" 
				"\t--time -t    : time space(default) \n"
				"\t--fourier -f : omega space  \n"
				"\t--window -w  : do windows \n"
				"\t--nskip -s (with option int): Skip early data\n"
				"\t--ndiff -d (with option int): \n"
				"\t--help -h    : usage of this function\n"
				"\t--verbose -v : equaivalent with help \n"
				" if you want to use stdin, you should used -  \n"
				, pName);
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
