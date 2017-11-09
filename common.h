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


typedef struct {
	real R, I;
} Cmplx;
/*!
 *  \brief  struct for complex real value
 */
typedef struct { 
	real xx,xy,xz;
	real yx,yy,yz;
	real zx,zy,zz;
} Rank2R3;
typedef struct { 
	real x,y,z;
} VecR3;

void real_tensor_copy_r1r1(VecR3*,VecR3*);
void real_tensor_copy_r2r2(Rank2R3*,Rank2R3*);
void real_tensor_product_r2_r1r1( Rank2R3* ,
		VecR3 *r1a, VecR3* ) ;
void real_tensor_product_r2_r2r2( Rank2R3* ,
		Rank2R3 *, Rank2R3* ) ;
void real_tensor_product_r2_r0r2( Rank2R3* ,
		real , Rank2R3*) ;
void real_tensor_add_r2_r2r2( Rank2R3* ,
		Rank2R3*, Rank2R3* ) ;
void real_tensor_sub_r2_r2r2( Rank2R3* ,
		Rank2R3 *, Rank2R3* ) ;
void real_tensor_zero_r2( Rank2R3* );
void real_tensor_zero_r1( VecR3* );
real real_tensor_sum_dig_r2(Rank2R3*);
real real_tensor_sum_offdig_r2(Rank2R3*);
