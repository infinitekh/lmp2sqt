#pragma once
/*
 * =====================================================================================
 *
 *       Filename:  common.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2016년 11월 24일 12시 47분 34초
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ph.D. Candidate KIM Hyeok (kh), ekh0324@gmail.com
 *        Company:  Konkuk University
 *
 * =====================================================================================
 */
#include <stdio.h>
typedef double real;
int fwrite_matrix( FILE *fout, real **m, int xsize, int ysize, real *rt, real *ct);
