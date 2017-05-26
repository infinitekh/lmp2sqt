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
#include <unistd.h>
#include <errno.h>


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
#define CHAR_ZERO  '0'
#define CHAR_SKIP  '#'
#define BUFF_LEN   1024
#define VERSION    2

typedef struct { real x,  y; int n ;} data;
ALLOC(data);
ALLOC(real);

size_t nlines (char *fn)
{
	if (!fn) return 0;

	size_t n = 0, noeol = 0;
	char buf[BUFF_LEN] = "";
	FILE *fp = fopen (fn, "r");

	if (!fp) return 0;
#ifndef NDEBUG
	fprintf(stderr,"cat file : %s \n", fn);
#endif

	while (fgets (buf, BUFF_LEN, fp)) {
#ifndef NDEBUG
			fprintf(stderr,"%s",buf);
#endif
		if( strchr(buf, '#' )) continue;
		noeol = 0;
		if (!strchr (buf, '\n')) {
			noeol = 1;  /* noeol flag for last line */
			continue;
		}
		n++;
	}

	if (noeol) n++;     /* check if noeol, add 1 */

	fseek(fp,0,SEEK_SET);

#ifndef NDEBUG
	printf("cat file END: %s \n", fn);
#endif
	fclose (fp);

	return n;
}
int compare_real( const void * x, const void * y)
{
	real sub2   = (*(real*) x - *(real*) y);
	int ret;
	if (sub2>0.0)  ret=1;
	else if (sub2<0.0)  ret=-1;
	else ret =0;

	return ret;
}

int main ( int argc, char **argv)
{

	char* fName;
	FILE* input;
	real* stats_data;
	real* stats_data_err;
	real* rawdata;
	int lines;
	char* bp;
	char buff[BUFF_LEN];

	int k,i,j,nColumn=2;

	int param_opt;

	opterr   = 0;

	while( -1 != ( param_opt = getopt( argc, argv, "n:h")) )
	{
		switch( param_opt)
		{
			case  'n'   :  
				nColumn = atoi(optarg);
				break;
			case  '?'   :  
				printf( "알 수 없는 옵션: %c\n", optopt);

			case 'h'    :
				printf("%s -n {nColumn=2} file\n"
						"Version : %d\n", argv[0], VERSION);
				exit(1);
				break;
		}
	}

	if (argc <2){
		printf("argc %d\n", argc);
				printf("%s -n {nColumn=2} file\n", argv[0]);
		exit(1);
	}

	for (i=optind; i<argc; i++) {
		fprintf(stderr,"%d %s\n", i,argv[i]);
	}

	fName = argv[optind];


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

	lines=nlines( fName);

	stats_data = alloc_real( (lines+1)*nColumn);
	stats_data_err = alloc_real( (lines+1)*nColumn);
	rawdata = alloc_real( (lines+1)*nColumn);
	/*-----------------------------------------------------------------------------
	 *   Alloc memory and initialization
	 *-----------------------------------------------------------------------------*/
	real w;
	int nMaxRawdata=0;
	i=0;
	j=0;
	// The Preamble
	/* 	while (1) {
	*/

	while(  bp = fgets (buff, BUFF_LEN, input) ) {
		j++;
		if('#' == bp[0]) {
			fputs("catch ## preamble\n",stderr);
			fputs(bp,stderr);
			continue;	
		};
		if( (strchr(bp, '\n' ) - bp) < 6)  continue;
		if( strlen(bp) <6) continue;
		for ( k = 0; k < nColumn; k += 1 ) {
			w = strtod (bp, &bp);
#ifndef NDEBUG
			fprintf(stderr," %8.4f",w);
#endif
			rawdata[i] = w;
			i++;
		}
#ifndef NDEBUG
		fputs("\n",stderr);
#endif
		if(j==lines)
			break;
	}
	nMaxRawdata =i;


	j=0;
#ifndef NDEBUG
	while (j<nMaxRawdata) {
		if (j%nColumn ==0)	{
			printf("\n%8.4f ", rawdata[j]);
		}
		else {
			printf("%8.4e ", rawdata[j]);
		}
		j++;
	}
	fputs("\n DataCheck\n\nqsort\n",stderr);
#endif
	/* 	}
	*/
	int nRecord = nMaxRawdata/ nColumn;
	qsort( rawdata, nRecord, sizeof(real)*nColumn, compare_real);

#ifndef NDEBUG
	j=0;
	while (j<nMaxRawdata) {
		if (j%nColumn ==0)	{
			printf("\n%8.4f ", rawdata[j]);
		}
		else {
			printf("%8.4e ", rawdata[j]);
		}
		j++;
	}
	puts("\n qsort END");
#endif

	real record[nColumn], record2[nColumn],x,pre_x,y,y2;

	
	 j=0; int nDup=0;
	 // First data;
	 {
		 record[0] = x= rawdata[0];


		 for (k=1; k<nColumn; k++){
			 record[k] = rawdata[i+k];
		 }

		 i=nColumn; nDup =1;
		 pre_x = x;
	 }
	while ( i<nMaxRawdata) {

		 x= rawdata[i];

		if( pre_x == x) {
			nDup ++;
			for (k=1; k<nColumn; k++){
				y = rawdata[i+k];
				record[k] +=y; 
				record2[k] += y*y;
			}
		} // END if 
		else { 
			// stats record
			{
				stats_data[j+0] = record[0];
				for (k=1; k<nColumn; k++){
					y = record[k]/ (real)nDup;
					y2 = record2[k]/ (real)nDup;
					stats_data[j+k] = y;
					stats_data_err[j+k] = sqrt( (y2 - y*y )  /(real)nDup );
				}
			}
			// new record
			{
				record[0] = x;
				for (k=1; k<nColumn; k++){
					y = rawdata[i+k];
					record[k] =y; 
					record2[k] = y*y;
				}
				pre_x = x;
				nDup=1;
			}
			j+= nColumn;
		} // END else
		// Next data
		i+=nColumn;  
	} // END WHILE i<nMaxRawdata
	// last line stats record
	{
		stats_data[j+0] = record[0];
		for (k=1; k<nColumn; k++){
			y = record[k]/ (real)nDup;
			y2 = record2[k]/ (real)nDup;
			stats_data[j+k] = y;
			stats_data_err[j+k] = sqrt( (y2 - y*y )  /(real)nDup );
		}
	}

	nMaxRawdata = j+nColumn;

	
	printf("## nData %d nMaxRawdata  %d", nColumn-1, nMaxRawdata/nColumn);
	j=0;
	while (j<nMaxRawdata) {
		if (j%nColumn ==0)	{
			printf("\n%8.4f ", stats_data[j]);
		}
		else {
			printf("%8.4e %8.4e ", stats_data[j], stats_data_err[j] );
//			printf("%8.4e ", stats_data[j]);
		}
		j++;
	}
	
	puts("");



	return 0;
	
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
