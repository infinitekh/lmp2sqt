/*
 * =====================================================================================
 *
 *       Filename:  analysis_binary.c
 *
 *    Description:  
 *
 *        Version:  1.4 (Full Dqt)
 *        Created:  2018년 03월 06일 (화) 오후 07시 02분 42초
 *       Revision:  2018년 10월 30일 (화) 오후 02시 42분 34초
 *       Compiler:  gcc
 *
 *         Author:  Ph.D. Candidate KIM Hyeok (kh), ekh0324@gmail.com
 *        Company:  Konkuk University
 *
 * =====================================================================================
 */


#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "snapshot.h"
#include <errno.h>
#include <getopt.h>
#include <unistd.h>
#include "common.h"


#define ALLOC_pointer(type) type**  alloc_pointer_ ## type(size_t n) { \
	  return (type **) malloc(sizeof(type*)*n); \
}
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
#define NameVals(x,n)                     \
if (! strncmp (bp, #x, strlen (#x))) { \
	bp += strlen (#x);                   \
	for (int __nnn=0; ___nnn< n ; ___nnn++ ) \
		x[___nnn]  = strtod (bp, &bp);                \
}

#define BUFF_LEN 50000
ALLOC(double); ALLOC(real); ALLOC(int); ALLOC(Cmplx);
ALLOC_pointer(real); ALLOC_pointer(int);
int  type[]    = { 1,   1,  1,  0 };
int  scail[]    = { 1,   1,  1,  0 };
char *header[] = {
		"full-density"   	 ,                        // 0
		"self-density"			,                       // 1
		"cross-density" 		,                       // 2
		"self-vanHove"                              // 3
		};
int header_flag_calc[] = {
	1,
	1,
	1,
	1 };
int header_flag_more[] = {
	1,
	1,
	1,
	0 };

/* char *header[] = {
 * 		"full-cur-long"  	 ,                        // 0
 * 		"full-cur-trans" 	 ,                        // 1
 * 		"full-mag-long"  	 ,                        // 2
 * 		"full-mag-trans" 	 ,                        // 3
 * 		"full-density"   	 ,                        // 4
 * 		"self-cur-long"			,                       // 5
 * 		"self-cur-trans"		,                       // 6
 * 		"self-mag-long"			,                       // 7
 * 		"self-mag-trans"		,                       // 8
 * 		"self-density"			,                       // 9
 * 		"cross-cur-long"		,                       // 10
 * 		"cross-cur-trans"		,                       // 11
 * 		"cross-mag-long"		,                       // 12
 * 		"cross-mag-trans"		,                       // 13
 * 		"cross-density" 		,                       // 14
 * 		"self-vanHove"                              // 15
 * 		};
 * int header_flag_calc[] = {
 * 	0,0,1,1,1,
 * 	0,0,1,1,1,
 * 	0,0,1,1,1,
 * 	1 };
 * int header_flag_more[] = {
 * 	0,0,1,1,1,
 * 	0,0,1,1,1,
 * 	0,0,1,1,1,
 * 	0 };
 */
int nDataTypes = sizeof(header)/sizeof(char*);
void PrintHelp ( char *pName,int);
void FftComplex (Cmplx *a, int size);
real finite_diff1_log(real* diff, real* record,int size, real dt);
real finite_diff1(real* diff, real* record, int size, real dt);
real finite_diff2_log(real* diff2,real* diff, real* record, int size, real dt);
real finite_diff2(real* diff2,real* diff, real* record, int size, real dt);

int omegaMax,doFourier,doWindow;
int nData,nCSpatial,nSet,nSetSkip,
		nv, nCTime, NValDiff;
real deltaT,deltaTCorr,kVal,kVal2,rVal;
void Print_R2_data ( FILE*, real*);

static int verbose_flag;
extern char* txtCorr;



/*!
 *  \brief  important 
 */
long int g_weight =0;
real **corrSum, **corrSumSq, **corrSumErr,
		 **Fqt, **Dqt, **Dqt_c,**Hqt,**Hqt_c,**Fqt1,damp, 
		 **GammaQT, tMax, w,x, qVal, qVal2;
Cmplx *work;
real valGamma, valDq, valSq, Fq0;
int j,k,n,nT,nnT,nnnT,cT, pT,ppT, pppT, endT, beginT,endT;
real fFqtUnderLimit= exp(-2.5);
int flag_ul;
char *bp, *fName, buff[BUFF_LEN], *lmpFileName, output_filename[BUFF_LEN];
FILE *input, *output,*input_text;
Snapshot* pSnap;

char char_do_fourier[2];


//-----------------------------------------------------------------------------
//  Define functions
//-----------------------------------------------------------------------------
void init_memory();
void parser_option(int argc, char** argv);
void import_settings(char*);
void load_raw_data();
void scaling_data();
void do_fourier();
void do_not_fourier();
void print_output ();

int main ( int argc, char **argv)
{
	/* 	int nDataTypes = sizeof(header);
	 * 	int tempa      = sizeof(char*);
	 * 	printf( "%d %d \n ", nDataTypes, tempa);
	 * 	exit(1);
	 */


	/*-------------------------------------
	 *  Argument Check!! And open a file.
	 *-------------------------------------*/
	parser_option(argc, argv);

	omegaMax = 10.;
	tMax = 100.;

//	char inputFilename[100]= "in.lmp2sqt";
	import_settings("in.lmp2sqt");


/*---------------------------------------------
*   Alloc memory and initialization
*----------------------------------------------*/
	init_memory();
	load_raw_data();

	scaling_data();

	char_do_fourier[1] = '\0';
	if (doFourier) {
		do_fourier();
		char_do_fourier[0] = 'w';

	}   //  print A(q,omega) or A(r,omega)
	else {   // print A(q,t), or A(r,t)
		do_not_fourier ();
		char_do_fourier[0] = 't';
	} // else end

	print_output () ;
}

void Print_R2_data ( FILE* fp, real* datas)
{
	real x; int nr, n ;
	fprintf(fp,"#%d", nCSpatial);
	for (nr = 0; nr < nCSpatial; nr += 1 ) {
		x= (nr+1)*kVal; 
		printf("### kVal = %d * %f= %f\n", nr+1,kVal,x);
		fprintf (fp, " %9.4e", x);
	}
	fprintf (fp, "\n");

	for (n=0; n < nv; n++) {
		if (doFourier) x = n * omegaMax / nv;
		else x = n * deltaTCorr;
		printf ( "%9.4e", x);
		fprintf (fp, "%9.4e", x);

		for ( nr = 0; nr < nCSpatial; nr += 1 ) {
			printf (" %9.4e", datas[nr * nCTime +n]);
			fprintf (fp," %9.4e", datas[nr * nCTime +n]);
		}
		fprintf (fp,"\n");
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
			"\t--ulimit -r (with option float:default->exp(-2.5)): \n"
			"\t--help -h    : usage of this function\n"
			"\t--verbose -v : equaivalent with help \n"
			"\t--text  : text style data file \n"
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
void init_memory() 
{
	corrSum = alloc_pointer_real(nDataTypes);
	Fqt1 = alloc_pointer_real(nDataTypes);
	Dqt = alloc_pointer_real(nDataTypes);
	Dqt_c = alloc_pointer_real(nDataTypes);
	Fqt = alloc_pointer_real(nDataTypes);
	Hqt = alloc_pointer_real(nDataTypes);
	Hqt_c = alloc_pointer_real(nDataTypes);
	GammaQT = alloc_pointer_real(nDataTypes);

	corrSumSq = alloc_pointer_real(nDataTypes);
	corrSumErr = alloc_pointer_real(nDataTypes);

	for (int j = 0; j < nDataTypes; j ++) {
		corrSum[j] = alloc_real(nCSpatial * nCTime);
		Fqt1[j] = alloc_real(nCSpatial * nCTime);
		Dqt[j] = alloc_real(nCSpatial * nCTime);
		Dqt_c[j] = alloc_real(nCSpatial * nCTime);
		Fqt[j] = alloc_real(nCSpatial * nCTime);
		Hqt[j] = alloc_real(nCSpatial * nCTime);
		Hqt_c[j] = alloc_real(nCSpatial * nCTime);
		GammaQT[j] = alloc_real(nCSpatial * nCTime);

		corrSumSq[j] = alloc_real( nCSpatial * nCTime);
		corrSumErr[j] = alloc_real( nCSpatial * nCTime);
		for (int n =0; n < nCSpatial* nCTime; n++) {
			corrSum [j][n] = 0;
			corrSumSq [j][n] = 0;
			corrSumErr [j][n] = 0;
			Fqt [j][n] = 0.;
			Fqt1 [j][n] = NAN;
			Dqt [j][n] = NAN;
			Dqt_c [j][n] = NAN;
			Hqt [j][n] = NAN;
			Hqt_c [j][n] = NAN;
			GammaQT [j][n] = NAN;
		}
	}
}
void parser_option(int argc, char** argv)
{
	n = 1;
	if (-- argc <1 || ! strcmp (argv[1], "-h")) PrintHelp (argv[0], __LINE__);
	doFourier =0;
	doWindow =0;
	nSetSkip = -1;
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
		{"text", no_argument, 0, 1  },
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
}

void import_settings(char* inputFilename)
{

	input = fopen( inputFilename, "r");
	input_text = fopen( inputFilename, "r");
	while (1) {
		bp = fgets (buff, BUFF_LEN, input);
		if (*bp == CHAR_MINUS) break;
		NameVal (deltaT);
		NameVal (nCSpatial);
		NameVal (nCTime);
		NameVal (kVal);
		//		NameString (lmpFileName);
	}
	deltaTCorr =  deltaT;
	kVal2 = kVal*kVal;
	fclose(input);

	if(!strcmp(fName,"-")) {
		input = stdin;
	} else {
		input = fopen(fName,"rb");
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

}
void load_raw_data_a_pack(real* rbuffer, int weight)
{
	//		if ( NULL == (pSnap = read_dump (input ) )) break;
	int suc_count = 0, full_count = strlen(txtCorr) ;
	// To reach end of file, break!!
  fread(full_count,sizeof(int),1,input);
	suc_count = fread (bp,sizeof(char) ,full_count,input);
	if ( full_count != suc_count)  return;
	int nTypes;

	if (! strncmp (bp, txtCorr, strlen (txtCorr))) {

		fread (&nTypes,sizeof(int) ,1,input);
		for ( j =0; j < nTypes; j++) {
			real col2 ; 
			real col3 ; 
			real col4 ; 

			/* 				char str_temp[100] = "# ";
			 * 				strncpy(str_temp +2 , header[j], strlen(header[j]));
			 */
			char str_temp[100];
			int strlength = strlen(header[j]);

			fread(strlength,sizeof(int),1,input);
			fread(str_temp, sizeof(char), strlength, input);

			fread(&col2, sizeof(real),1,input);
			fread(&col3, sizeof(real),1,input);
			fread(&col4, sizeof(real),1,input);

			int ret_cmp = strncmp (header[j], str_temp,  strlen (header[j]));
#ifndef NDEBUG
			printf( "str_temp : %s\n"
					"col2 : %f\n"
					"col3 : %f\n"
					"col4 : %f\n", str_temp,col2,col3,col4);
			printf("ret_cmp %d, %s -- %s   %d\n", ret_cmp,
					header[j],str_temp,(int)strlen(str_temp));

#endif

			// header types check(not completed)
			if ( !ret_cmp )  {
				fread(rbuffer, sizeof(real), (nCTime*nCSpatial), input);
				for ( n =0; n<nCTime; n ++) {
					//						bp = fgets (buff, BUFF_LEN, input);
					for ( k = 0; k < nCSpatial; k += 1 ) {
						w = weight*rbuffer[ n*nCSpatial + k ]; // transpose
						corrSum[j][k * nCTime + n] += w;
						corrSumSq[j][k * nCTime + n] += Sqr(w);
					}
				}
			}
		}
		g_weight +=weight;
	}
}
void load_raw_data_text()// Need Check!!!!!
{
	real* rbuffer = alloc_real(nCSpatial * nCTime);
	// The Preamble
	if (doFourier)
		work = alloc_Cmplx( 2 * (nCTime -1));
	nData =0;
	nSet =0;
	while (1) {
		//		if ( NULL == (pSnap = read_dump (input ) )) break;

		// To reach end of file, break!!
		if (! (bp = fgets (buff, BUFF_LEN, input_text))) break;
		int nTypes;
		if (! strncmp (bp, txtCorr, strlen (txtCorr))) {
			++ nSet;
			if (nSet < nSetSkip) continue;
			++ nData;
			fread (&nTypes,sizeof(int) ,1,input);
			// if (nTypes > nDataTypes) break; // error;

			for ( j =0; j < nTypes; j++) {
				real col2 ; 
				real col3 ; 
				real col4 ; 

				/* 				char str_temp[100] = "# ";
				 * 				strncpy(str_temp +2 , header[j], strlen(header[j]));
				 */
				char str_temp[100];
				/* 				char str_temp[100] = "# ";
				 * 				strncpy(str_temp +2 , header[j], strlen(header[j]));
				 */
				bp = fgets (buff, BUFF_LEN, input_text); 

				size_t strlength = strlen(header[j]);
				fread(str_temp, sizeof(char), strlength, input_text);
				str_temp[strlength]='\0';
				fread(&col2, sizeof(real),1,input_text);
				fread(&col3, sizeof(real),1,input_text);
				fread(&col4, sizeof(real),1,input_text);

				int ret_cmp = strncmp (header[j], str_temp,  strlen (header[j]));
#ifndef NDEBUG
				printf( "str_temp : %s\n"
						"col2 : %f\n"
						"col3 : %f\n"
						"col4 : %f\n", str_temp,col2,col3,col4);
				printf("ret_cmp %d, %s -- %s   %d\n", ret_cmp,
						header[j],str_temp,(int)strlen(str_temp));

#endif

				// header types check(not completed)
				if ( !ret_cmp )  {
					fread(rbuffer, sizeof(real), (nCTime*nCSpatial), input_text);
					if ( !ret_cmp )  {
						for ( n =0; n<nCTime; n ++) {
							bp = fgets (buff, BUFF_LEN, input_text);
#ifndef NDEBUG
							printf("TIME%5d : %s\n",n,buff);
#endif
							for ( k = 0; k < nCSpatial; k += 1 ) {
								w = strtod (bp, &bp);
#ifndef NDEBUG
								//							printf("LINE %4d : %8.4f",__LINE__,w);
#endif
								corrSum[j][k * nCTime + n] += w;
								corrSumSq[j][k * nCTime + n] += Sqr(w);
								/* 						sleep(1);
								 * 						printf("why is this value nan? %f %f  \n",corrSum[j][k*nCTime +n], w);
								 */
							}
#ifndef NDEBUG
							puts("\n");
#endif
						}
					}
				}
			}
		}
	}
	fclose (input_text);
	printf ("%d\n", nData);
}
void load_raw_data()
{
	real* rbuffer = alloc_real(nCSpatial * nCTime);
	// The Preamble
	if (doFourier)
		work = alloc_Cmplx( 2 * (nCTime -1));
	nData =0;
	nSet =0;
	while (1) {
		//		if ( NULL == (pSnap = read_dump (input ) )) break;

		int suc_count = 0, full_count = strlen(txtCorr) ;
		// To reach end of file, break!!
		suc_count = fread (bp,sizeof(char) ,full_count,input);
		if ( full_count != suc_count)  break;
		int nTypes;

		if (! strncmp (bp, txtCorr, strlen (txtCorr))) {
			++ nSet;
			if (nSet < nSetSkip) continue;
			++ nData;

			fread (&nTypes,sizeof(int) ,1,input);
			for ( j =0; j < nTypes; j++) {
				real col2 ; 
				real col3 ; 
				real col4 ; 

/* 				char str_temp[100] = "# ";
 * 				strncpy(str_temp +2 , header[j], strlen(header[j]));
 */
				char str_temp[100];
				size_t strlength = strlen(header[j]);
				fread(str_temp, sizeof(char), strlength, input);

				fread(&col2, sizeof(real),1,input);
				fread(&col3, sizeof(real),1,input);
				fread(&col4, sizeof(real),1,input);

				kVal = col2;
				deltaT = col3;
				rVal = col4;

				int ret_cmp = strncmp (header[j], str_temp,  strlen (header[j]));
#ifndef NDEBUG
				printf( "str_temp : %s\n"
						"col2 : %f\n"
						"col3 : %f\n"
						"col4 : %f\n", str_temp,col2,col3,col4);
				printf("ret_cmp %d, %s -- %s   %d\n", ret_cmp,
						header[j],str_temp,(int)strlen(str_temp));

#endif

				// header types check(not completed)
				if ( !ret_cmp )  {
					fread(rbuffer, sizeof(real), (nCTime*nCSpatial), input);
					for ( n =0; n<nCTime; n ++) {
//						bp = fgets (buff, BUFF_LEN, input);
						for ( k = 0; k < nCSpatial; k += 1 ) {
							w = rbuffer[ n*nCSpatial + k ]; // transpose
							corrSum[j][k * nCTime + n] += w;
							corrSumSq[j][k * nCTime + n] += Sqr(w);
						}
					}
				}
			}
		}
	}
	fclose (input);
	printf ("%d\n", nData);
}
void scaling_data()
{
	if (nData>1) {
		for ( j = 0; j < nDataTypes; j += 1 ) {
			for ( n = 0; n < nCSpatial*nCTime; n += 1 ) {
				corrSum[j][n] /= nData;
				corrSumSq[j][n] = sqrt ( corrSumSq[j][n]/nData-Sqr(corrSum[j][n]));
				corrSumErr[j][n] =  corrSumSq[j][n]/ sqrt(nData);
			}
		}
	}
}
void do_fourier () 
{
	for ( j = 0; j < 3; j += 1 ) {
		for ( k = 0; k < nCSpatial; k += 1 ) {
			for ( n = 0; n < nCTime; n += 1 ) {
				if (doWindow) damp = (nCTime - n ) / (nCTime +.5);
				else damp = 1.;
				CSet (work[n], corrSum[j][k * nCTime +n] * damp,0.);
			}

			for ( n = nCTime; n < 2*(nCTime-1); n += 1 ) 
				work[n] = work[2* (nCTime -1 ) - n ];
			FftComplex (work, 2 * nCTime -2);

			for ( n = 0; n < nCTime; n += 1 ) 
				corrSum[j][k * nCTime +n] = work[n].R;
		}
	}
	omegaMax= Min(omegaMax, M_PI / deltaTCorr);
	nv = nCTime * omegaMax / (M_PI / deltaTCorr);
}
void do_not_fourier ()
{
	for ( j = 0; j < nDataTypes; j += 1 ) {
		if (header_flag_calc[j] ){
			/*-------------------------------------------------------
			 *  D(q) ~ S(q) ~ or density and  spin longi,transverse
			 *--------------------------------------------------------*/
			sprintf( output_filename, "%s00.out", header[j] );
			FILE *output = fopen( output_filename, "w");
			fprintf(output,"#kVal F0 F1_0 Gamma Dq0\n"
					"#deltaT = %9.4f\n", deltaT);
			int nValDiffTemp=NValDiff;
			for (k = 0; k < nCSpatial; k ++) {
				qVal = (k+1)*kVal; qVal2=qVal*qVal;
				Fq0 = corrSum[j] [k*nCTime];
				double lnFq0 = log(Fq0);

				/*--------------------------------------------
				 *      get slope on nValDiff time or
				 *     where is return of function lower than 0.2
				 *--------------------------------------------*/
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
				cT = beginT = k*nCTime;
				endT = (k+1)*nCTime-1; // biggest time 
first:
				{
					nnT = cT+2; nT = cT+1;  //Forward Records
					Fqt1[j] [cT] =   ( (-Xp2+4.*Xp1 -3.* X ) / (2.*deltaT )) ;  // Forward O(h^2) first Derivative
					GammaQT[j] [cT] =   ( (-Xlogp2+4.*Xlogp1 -3.* Xlog ) / (2.*deltaT )) ;  // Forward O(h^2) first Derivative
					Dqt [j] [cT] = - GammaQT[j][cT]/qVal2;
					Dqt_c [j] [cT] = Dqt [j] [cT];
					Hqt [j] [cT] =  Dqt[j][cT] * Fq0;
					Hqt_c [j] [cT] =  Dqt_c[j][cT] * Fq0;
				}
				flag_ul = 0;
middle:
				for (cT= beginT+1; cT < endT; cT++) {
					if( corrSum[j][nT] <fFqtUnderLimit ) {
						flag_ul = 1;
						break; // Fail to do second central derivative
					}
					nT = cT+1;   //Forward Records
					pT = cT-1;   //Backward Records 
					// Central O(h^2)  Derivative
					Fqt1[j] [cT] =   ( (Xp1- Xm1) / (2.*deltaT )) ;  
					GammaQT[j] [cT] =   ( (Xlogp1 - Xlogm1 ) / (2.*deltaT )) ;  	
					Dqt [j] [cT] = - GammaQT[j][cT]/qVal2;
					Dqt_c [j] [cT] = - (Xlog - lnFq0)/(qVal2*deltaT*(cT-beginT)) ;
					Hqt [j] [cT] =  Dqt[j][cT] * Fq0;
					Hqt_c [j] [cT] =  Dqt_c[j][cT] * Fq0;
				}
last:
				if( flag_ul == 0 ) {
					cT = endT;
					{
						// Forward O(h^2) first Derivative  Last term
						ppT = cT-2; pT=cT -1;
						Fqt1[j] [cT] =   ( (+Xm2-4.*Xm1 +3.* X ) / (2.*deltaT )) ;  
						GammaQT[j] [cT] =   ( (+Xlogm2-4.*Xlogm1 +3.* Xlog ) / (2.*deltaT )) ;  
						Dqt [j] [cT] = - GammaQT[j][cT]/qVal2;
						Dqt_c [j] [cT] = - (Xlog - lnFq0)/(qVal2*deltaT*(cT-beginT)) ;
						Hqt [j] [cT] =  Dqt[j][cT] * Fq0;
						Hqt_c [j] [cT] =  Dqt_c[j][cT] * Fq0;
					}
				}

				cT = cT < k*nCTime +NValDiff  ? cT: k*nCTime+NValDiff;
				fprintf(output, "%lf %le %le %le %le \n", qVal, 
						//							corrSum[j][k*nCTime], Fqt1[j][k*nCTime], 
						corrSum[j][cT],Fqt1[j][cT], 
						GammaQT[j][cT], Dqt[j][cT]);

			}

			fclose(output);

			//      scaling by function of  t=0
			for ( k = 0; k < nCSpatial; k += 1 ) {
#pragma omp parallel for
				for ( n = 1; n < nCTime; n += 1 ){ 
					Fqt [j][k*nCTime + n] = corrSum [j][k*nCTime + n] ;
					corrSum[j][k * nCTime +n] /= corrSum[j][k*nCTime];
				}
				Fqt [j][k*nCTime ] = corrSum [j][k*nCTime ] ;
				corrSum[j][k * nCTime] = 1.;
			}
		}
	}
	tMax = Min ( tMax, (nCTime -1) * deltaTCorr);
	nv = nCTime * tMax /  (nCTime - 1) / deltaTCorr;
	printf("nv = %d, tMax = %f, nCTime = %d, deltaTCorr = %f\n" , nv ,tMax, nCTime, deltaTCorr);
}
void print_output () 
{
	char fNGammaQT[100] ;
	char fNDqt[100] ;
	char fNDqt_c[100] ;
	char fNHqt[100] ;
	char fNHqt_c[100] ;
	char fNFqt[100] ;
	char fNFqt_scaled[100] ;
	char fNFqt1[100];

	for ( j = 0; j < nDataTypes; j += 1 ) {
		if ( header_flag_calc[j] ==1 ) {

			sprintf(fNFqt_scaled, "Fq%s_scaled.%s.info",char_do_fourier, header[j]);
			FILE* fFqt_scaled = fopen(fNFqt_scaled, "w");
			Print_R2_data(fFqt_scaled, corrSum[j]);
			fclose(fFqt_scaled );

			sprintf(fNFqt, "Fq%s.%s.info",char_do_fourier, header[j]);
			FILE* fFqt = fopen(fNFqt, "w");
			Print_R2_data(fFqt,  Fqt[j]);
			fclose(fFqt );

/* 		corrSum[j] = alloc_real(nCSpatial * nCTime);
 * 		Fqt1[j] = alloc_real(nCSpatial * nCTime);
 * 		Dqt[j] = alloc_real(nCSpatial * nCTime);
 * 		Fqt[j] = alloc_real(nCSpatial * nCTime);
 * 		Hqt[j] = alloc_real(nCSpatial * nCTime);
 * 		GammaQT[j] = alloc_real(nCSpatial * nCTime);
 */

			if ( header_flag_more[j] ==1 ) {
				sprintf(fNGammaQT,"GammaQ%s.%s.info", char_do_fourier,header[j]);
				FILE* fGammaQT = fopen(fNGammaQT, "w");
				Print_R2_data(fGammaQT, GammaQT[j]);
				fclose(fGammaQT);

				sprintf(fNDqt_c, "Dq%s_c.%s.info", char_do_fourier, header[j]);
				FILE* fDqt_c = fopen(fNDqt_c, "w");
				Print_R2_data(fDqt_c, Dqt_c[j]);
				fclose(fDqt_c );

				sprintf(fNDqt, "Dq%s.%s.info", char_do_fourier, header[j]);
				FILE* fDqt = fopen(fNDqt, "w");
				Print_R2_data(fDqt, Dqt[j]);
				fclose(fDqt );

				sprintf(fNHqt_c, "D0Hq%s_c.%s.info", char_do_fourier, header[j]);
				FILE* fHqt_c = fopen(fNHqt_c, "w");
				Print_R2_data(fHqt_c, Hqt_c[j]);
				fclose(fHqt_c );

				sprintf(fNHqt, "D0Hq%s.%s.info", char_do_fourier, header[j]);
				FILE* fHqt = fopen(fNHqt, "w");
				Print_R2_data(fHqt, Hqt[j]);
				fclose(fHqt );

				sprintf(fNFqt1, "Fq%s1.%s.info", char_do_fourier, header[j]);
				FILE* fFqt1 = fopen(fNFqt1, "w");
				Print_R2_data(fFqt1,  Fqt1[j]);
				fclose(fFqt1);
			}
		}
	}
}
