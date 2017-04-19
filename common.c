/*
 * =====================================================================================
 *
 *       Filename:  common.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2017년 04월 19일 11시 18분 46초
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ph.D. Candidate KIM Hyeok (kh), ekh0324@gmail.com
 *        Company:  Konkuk University
 *
 * =====================================================================================
 */

#include "common.h"
	int fwrite_matrix( FILE *fout, real **m, int xsize, int ysize, real *rt, real *ct)
{
	int j;
	int status;
	real length = ysize;

	if ((status = fwrite((char *) &length, sizeof(real), 1, fout)) != 1) {
		fprintf(stderr, "fwrite 1 returned %d\n", status);
		return (0);
	}
	fwrite((char *) ct, sizeof(real), ysize, fout);
	for (j = 0; j < xsize; j++) {
		fwrite((char *) &rt[j], sizeof(real), 1, fout);
		fwrite((char *) (m[j]), sizeof(real), ysize, fout);
	}

	return (1);
}
