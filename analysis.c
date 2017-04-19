/*
 * =====================================================================================
 *
 *       Filename:  analysis.cpp
 *
 *    Description:  
 *
 *        Version:  1.1 (Full Dqt)
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
void PrintHelp ( char *pName,int);
void FftComplex (Cmplx *a, int size);
real finite_diff1_log(real* diff, real* record,int size, real dt);
real finite_diff1(real* diff, real* record, int size, real dt);
real finite_diff2_log(real* diff2,real* diff, real* record, int size, real dt);
real finite_diff2(real* diff2,real* diff, real* record, int size, real dt);


static int verbose_flag;
int main ( int argc, char **argv)
{
	Cmplx *work;
	real *corrSum[nDataTypes], *corrSumSq[nDataTypes], *corrSumErr[nDataTypes],
	*Fqt[nDataTypes],
	*Dqt[nDataTypes],*Hqt[nDataTypes],*corrSumD1[nDataTypes],damp, deltaT, deltaTCorr,
	*GammaQT[nDataTypes],omegaMax, tMax, w,x,  kVal, kVal2, qVal, qVal2;
	real valGamma, valDq, valSq, Fq0;
	int doFourier,doWindow;
	int j,k,n,nT,nnT,nnnT,cT, pT,ppT, pppT;
	int nData,nFunCorr,nSet,nSetSkip,
			nv, nValCorr, NValDiff;
	real fFqtUnderLimit= exp(-2.5);
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
	if (-- argc <1 || ! strcmp (argv[1], "-h")) PrintHelp (argv[0], __LINE__);
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
		{"ulimit", required_argument, 0, 'u'},
		{"help", no_argument, 0, 'h'},
		{0,0,0,0}
	};
	int option_index=0;

	while( -1 !=( param_opt = getopt_long( argc, argv, "ftws:d:u:h", long_options, &option_index)))
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
			case  'u'   :  
				fFqtUnderLimit = atof(optarg);
				break;
			case  '?'   :  
				printf( "알 수 없는 옵션: %cn", optopt);
			case 'h'    :
				PrintHelp(argv[0], __LINE__);
				break;
		}
	}
	if (optind <1) PrintHelp (argv[0], __LINE__);
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
		corrSumD1[j] = alloc_real(nFunCorr * nValCorr);
		Dqt[j] = alloc_real(nFunCorr * nValCorr);
		Fqt[j] = alloc_real(nFunCorr * nValCorr);
		Hqt[j] = alloc_real(nFunCorr * nValCorr);
		GammaQT[j] = alloc_real(nFunCorr * nValCorr);

		corrSumSq[j] = alloc_real( nFunCorr * nValCorr);
		corrSumErr[j] = alloc_real( nFunCorr * nValCorr);
		for (n =0; n < nFunCorr* nValCorr; n++) {
			corrSum [j][n] = 0;
			corrSumD1 [j][n] = NAN;
			Dqt [j][n] = NAN;
			Hqt [j][n] = NAN;
			Fqt [j][n] = NAN;
			GammaQT [j][n] = NAN;
			corrSumSq [j][n] = 0;
			corrSumErr [j][n] = 0;
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
#ifndef NDEBUG
				printf("## header check : %s \n", buff);
#endif

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
/* 						sleep(1);
 * 						printf("why is this value nan? %f %f  \n",corrSum[j][k*nValCorr +n], w);
 */
					}
#ifndef NDEBUG
					puts("\n");
#endif
				}
				bp = fgets (buff, BUFF_LEN, input);// null line
#ifndef NDEBUG
				printf("## endline null string default : %s \n", buff);
#endif
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
				int nValDiffTemp=NValDiff;
				for (k = 0; k < nFunCorr; k ++) {
					qVal = (k+1)*kVal; qVal2=qVal*qVal;
					Fq0 = corrSum[j] [k*nValCorr];

/*-----------------------------------------------------------------------------
*      get slope on nValDiff time or where is return of function lower than 0.2
*-----------------------------------------------------------------------------*/
#define Xm3 ((corrSum[j][pppT]))
#define Xm2 ((corrSum[j][ppT]))
#define Xm1 ((corrSum[j][pT]) )
#define Xp3 ((corrSum[j][nnnT]))
#define Xp2 ((corrSum[j][nnT]))
#define Xp1 ((corrSum[j][nT]) )
#define X   ((corrSum[j][cT]) )
#define Xlogm3 (log(corrSum[j][pppT]))
#define Xlogm2 (log(corrSum[j][ppT]))
#define Xlogm1 (log(corrSum[j][pT]) )
#define Xlogp3 (log(corrSum[j][nnnT]))
#define Xlogp2 (log(corrSum[j][nnT]))
#define Xlogp1 (log(corrSum[j][nT]) )
#define Xlog   (log(corrSum[j][cT]) )
					cT = k*nValCorr;
					nnT = cT+2; nT = cT+1;  //Forward Records
					corrSumD1[j] [cT] =   ( (-Xp2+4.*Xp1 -3.* X ) / (2.*deltaT )) ;  // Forward O(h^2) first Derivative
					GammaQT[j] [cT] =   ( (-Xlogp2+4.*Xlogp1 -3.* Xlog ) / (2.*deltaT )) ;  // Forward O(h^2) first Derivative
					Dqt [j] [cT] = - GammaQT[j][cT]/qVal2;
					Hqt [j] [cT] =  Dqt[j][cT] * Fq0;


					for (cT= k*nValCorr+1; cT < (k+1)*nValCorr-1; cT++) {
						 nT = cT+1;   //Forward Records
						 pT = cT-1;   //Backward Records 

//						if( corrSum[j][nnT] <0.05 ) break; // Fail to do second central derivative
						if( corrSum[j][nT] <0.1 ) break; // Fail to do second central derivative
						corrSumD1[j] [cT] =   ( (Xp1- Xm1) / (2.*deltaT )) ;  // Central O(h^2)  Derivative
//						corrSumD1[j] [cT] =   ( (-Xp2+8.*Xp1-8.*Xm1+Xm2 ) / (12.*deltaT )) ;  // Central O(h^4) Derivative
						GammaQT[j] [cT] =   ( (Xlogp1 - Xlogm1 ) / (2.*deltaT )) ;  // Central O(h^2)  Derivative
//						GammaQT[j] [cT] =   ( (-Xlogp2+8.*Xlogp1-8.*Xlogm1+Xlogm2 ) / (12.*deltaT )) ;  // Central O(h^4) Derivative
						Dqt [j] [cT] = - GammaQT[j][cT]/qVal2;
						Hqt [j] [cT] =  Dqt[j][cT] * Fq0;
					}
					corrSumD1[j] [cT] =   ( (Xp1- Xm1 ) / (2.*deltaT )) ;  // Central O(h^2) first Derivative
					GammaQT[j] [cT] =   ( (Xlogp1- Xlogm1 ) / (2.*deltaT )) ;  // Central O(h^2) first Derivative
					Dqt [j] [cT] = - GammaQT[j][cT]/qVal2;
					Hqt [j] [cT] =  Dqt[j][cT] * Fq0;

					ppT = cT-2; pT=cT -1;
					corrSumD1[j] [cT] =   ( (+Xm2-4.*Xm1 +3.* X ) / (2.*deltaT )) ;  // Forward O(h^2) first Derivative
					GammaQT[j] [cT] =   ( (+Xlogm2-4.*Xlogm1 +3.* Xlog ) / (2.*deltaT )) ;  // Forward O(h^2) first Derivative
					Dqt [j] [cT] = - GammaQT[j][cT]/qVal2;
					Hqt [j] [cT] =  Dqt[j][cT] * Fq0;

					cT = cT < k*nValCorr +NValDiff  ? cT: k*nValCorr+NValDiff;
					fprintf(output, "%lf %le %le %le %le %le %le \n", qVal, corrSum[j][k*nValCorr], corrSumD1[j][k*nValCorr], 
							corrSum[j][cT],corrSumD1[j][cT], 
							GammaQT[j][cT], Dqt[j][cT]);

				}
				fclose(output);


				//      scaling by function of  t=0
				for ( k = 0; k < nFunCorr; k += 1 ) {
					for ( n = 1; n < nValCorr; n += 1 ){ 
						Fqt [j][k*nValCorr + n] = corrSum [j][k*nValCorr + n] ;
						corrSum[j][k * nValCorr +n] /= corrSum[j][k*nValCorr];
					}
					Fqt [j][k*nValCorr ] = corrSum [j][k*nValCorr ] ;
					corrSum[j][k * nValCorr] = 1.;
				}
				}
			}
			tMax = Min ( tMax, (nValCorr -1) * deltaTCorr);
			nv = nValCorr * tMax /  (nValCorr - 1) / deltaTCorr;
			printf("nv = %d, tMax = %f, nValCorr = %d, deltaTCorr = %f\n" , nv ,tMax, nValCorr, deltaTCorr);
		} // else end

		char fNGammaQT[100] ;
		char fNDqt[100] ;
		char fNHqt[100] ;
		char fNFqt[100] ;
		char fNFqt1[100];

		for ( j = 0; j < nDataTypes; j += 1 ) {
			
			sprintf(fNGammaQT,"GammaQT.%s.info", header[j]);
			sprintf(fNDqt, "Dqt.%s.info", header[j]);
			sprintf(fNHqt, "Hqt.%s.info", header[j]);
			sprintf(fNFqt, "Fqt.%s.info", header[j]);
			sprintf(fNFqt1, "Fqt1.%s.info", header[j]);

			FILE* fGammaQT = fopen(fNGammaQT, "w");
			FILE* fDqt = fopen(fNDqt, "w");
			FILE* fHqt = fopen(fNHqt, "w");
			FILE* fFqt = fopen(fNFqt, "w");
			FILE* fFqt1 = fopen(fNFqt1, "w");

/* 			printf("%s\n", header[j]);
 * 			fprintf(fGammaQT,"%s\n", header[j]);
 * 			fprintf(fDqt,"%s\n", header[j]);
 * 			fprintf(fHqt,"%s\n", header[j]);
 * 			fprintf(fFqt,"%s\n", header[j]);
 * 			fprintf(fFqt1,"%s\n", header[j]);
 */
  		fprintf(fGammaQT,"%d", nFunCorr);
  		fprintf(fDqt,"%d", nFunCorr);
  		fprintf(fHqt,"%d", nFunCorr);
  		fprintf(fFqt,"%d", nFunCorr);
  		fprintf(fFqt1,"%d", nFunCorr);
			for ( k = 0; k < nFunCorr; k += 1 ) {
				x= (k+1)*kVal; 
				fprintf (fGammaQT, "%9.4e", x);
				fprintf (fDqt, "%9.4e", x);
				fprintf (fHqt, "%9.4e", x);
				fprintf (fFqt, "%9.4e", x);
				fprintf (fFqt1, "%9.4e", x);
			}

			for (n=0; n < nv; n++) {
				if (doFourier) x = n * omegaMax / nv;
				else x = n * deltaTCorr;
				printf ( "%9.4e", x);
				fprintf (fGammaQT, "%9.4e", x);
				fprintf (fDqt, "%9.4e", x);
				fprintf (fHqt, "%9.4e", x);
				fprintf (fFqt, "%9.4e", x);
				fprintf (fFqt1, "%9.4e", x);

				for ( k = 0; k < nFunCorr; k += 1 ) {
					printf (" %9.4e", corrSum[j][k * nValCorr +n]);
					fprintf (fGammaQT," %9.4e", GammaQT[j][k * nValCorr +n]);
					fprintf (fDqt," %9.4e", Dqt[j][k * nValCorr +n]);
					fprintf (fHqt," %9.4e", Hqt[j][k * nValCorr +n]);
					fprintf (fFqt," %9.4e", Fqt[j][k * nValCorr +n]);
					fprintf (fFqt1," %9.4e", corrSumD1[j][k * nValCorr +n]);
				}
				printf ("\n");
				fprintf ( fGammaQT,"\n");
				fprintf (fDqt,"\n");
				fprintf (fHqt,"\n");
				fprintf (fFqt,"\n");
				fprintf (fFqt1,"\n");
			}
			fclose(fGammaQT);
			fclose(fDqt );
			fclose(fHqt );
			fclose(fFqt );
			fclose(fFqt1);
			sleep(1);
		}
	}

	void PrintHelp ( char *pName, int linenumber)
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
				"Error line number %d\n"
				, pName, linenumber);
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
real finite_diff1(real* diff, real* record,int size, real dt)
{ /// p previous n next
	real fac = 1./(2.*dt); int i;
	if (size > 3){
		diff[0] = fac* ( -record[2] + 4.*record[1] - 3.*record[0]);
		diff[size-1] = fac* ( record[size-3] - 4.*record[size-2] + 3.*record[size-1]);
	}
	else { diff[0] = 0; diff[size-1] = 0; }
	
	for ( i = 1; i < size-1; i += 1 ){
		diff[i] = fac* (record[i+1] - record[i-1]);
	}
}
real finite_diff1_log(real* diff, real* record,int size, real dt)
{ /// p previous n next
#define LOGREC(x) log(record[x])
	real fac = 1./(2.*dt); int i;
	if (size > 3){
		diff[0] = fac* ( -LOGREC(2) + 4.*LOGREC(1) - 3.*LOGREC(0));
		diff[size-1] = fac* ( +LOGREC(size-3) - 4.*LOGREC(size-2) + 3.*LOGREC(size-1));
	}
	else { diff[0] = 0; diff[size-1] = 0; }
	
	for ( i = 1; i < size-1; i += 1 ){
		diff[i] = fac* (LOGREC(i+1) - LOGREC(i-1));
	}
#undef LOGREC
}
real finite_diff2(real* diff2,real* diff, real* record,int size, real dt)
{ /// p previous n next
	real fac = 1./(2.*dt); int i;
	real fac2 = 1./(dt*dt);
	if (size > 4){
		diff[0] = fac* ( -record[2] + 4.*record[1] - 3.*record[0]);
		diff[size-1] = fac* ( record[size-3] - 4.*record[size-2] + 3.*record[size-1]);

		diff2[0] = fac2* ( -record[3] + 4.*record[2] - 5.*record[1]+ 2.*record[0]);
		diff2[size-1] =fac2*(+record[size-4]-4.*record[size-3]+5.*record[size-2]-2.*record[size-1]);
	}
	else { diff[0] = 0; diff[size-1] = 0; }
	
	for ( i = 1; i < size-1; i += 1 ){
		diff[i] = fac* (record[i+1] - record[i-1]);
		diff2[i] = fac2* (record[i+1]- 2.*record[i] + record[i-1]);
	}
}
real finite_diff2_log(real* diff2,real* diff, real* record,int size, real dt)
{ /// p previous n next
#define LOGREC(x) log(record[x])
	real fac = 1./(2.*dt); int i;
	real fac2 = 1./(dt*dt);
	if (size > 4){
		diff[0] = fac* ( -LOGREC(2) + 4.*LOGREC(1) - 3.*LOGREC(0));
		diff[size-1] = fac* ( +LOGREC(size-3) - 4.*LOGREC(size-2) + 3.*LOGREC(size-1));

		diff2[0] = fac2* ( -LOGREC(3) + 4.*LOGREC(2) - 5.*LOGREC(1)+ 2.*LOGREC(0));
		diff2[size-1] =fac2*(+LOGREC(size-4)-4.*LOGREC(size-3)+5.*LOGREC(size-2)-2.*LOGREC(size-1));
	}
	else { diff[0] = 0; diff[size-1] = 0; }
	
	for ( i = 1; i < size-1; i += 1 ){
		diff[i] = fac* (LOGREC(i+1) - LOGREC(i-1));
		diff2[i] = fac2* (LOGREC(i+1)- 2.*LOGREC(i) + LOGREC(i-1));
	}
#undef LOGREC
}
