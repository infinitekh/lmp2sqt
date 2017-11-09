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

void real_tensor_copy_r1r1(VecR3*  lvalue,VecR3* rvalue) {
	lvalue->x = rvalue->x;
	lvalue->y = rvalue->y;
	lvalue->z = rvalue->z;
}
void real_tensor_copy_r2r2(Rank2R3*  lvalue,Rank2R3* rvalue) {
	lvalue->xx = rvalue->xx;
	lvalue->xy = rvalue->xy;
	lvalue->xz = rvalue->xz;
	lvalue->yx = rvalue->yx;
	lvalue->yy = rvalue->yy;
	lvalue->yz = rvalue->yz;
	lvalue->zx = rvalue->zx;
	lvalue->zy = rvalue->zy;
	lvalue->zz = rvalue->zz;
}
void real_tensor_product_r2_r1r1( Rank2R3* r2,
		VecR3 *r1a, VecR3* r1b) {
	r2->xx = r1a->x * r1b->x;
	r2->xy = r1a->x * r1b->y;
	r2->xz = r1a->x * r1b->z;
	r2->yx = r1a->y * r1b->x;
	r2->yy = r1a->y * r1b->y;
	r2->yz = r1a->y * r1b->z;
	r2->zx = r1a->z * r1b->x;
	r2->zy = r1a->z * r1b->y;
	r2->zz = r1a->z * r1b->z;
}
void real_tensor_product_r2_r0r2( Rank2R3* r2,
		real a, 
		Rank2R3* r2b) {
	r2->xx = a * r2b->xx;
	r2->xy = a * r2b->xy;
	r2->xz = a * r2b->xz;
	r2->yx = a * r2b->yx;
	r2->yy = a * r2b->yy;
	r2->yz = a * r2b->yz;
	r2->zx = a * r2b->zx;
	r2->zy = a * r2b->zy;
	r2->zz = a * r2b->zz;
}
void real_tensor_product_r2_r2r2( Rank2R3* r2,
		Rank2R3 * r2a, 
		Rank2R3* r2b) {
	r2->xx = r2a->xx * r2b->xx;
	r2->xy = r2a->xy * r2b->xy;
	r2->xz = r2a->xz * r2b->xz;
	r2->yx = r2a->yx * r2b->yx;
	r2->yy = r2a->yy * r2b->yy;
	r2->yz = r2a->yz * r2b->yz;
	r2->zx = r2a->zx * r2b->zx;
	r2->zy = r2a->zy * r2b->zy;
	r2->zz = r2a->zz * r2b->zz;
}
void real_tensor_sub_r2_r2r2( Rank2R3* r2,
		Rank2R3 * r2a, 
		Rank2R3* r2b) {
	r2->xx = r2a->xx - r2b->xx;
	r2->xy = r2a->xy - r2b->xy;
	r2->xz = r2a->xz - r2b->xz;
	r2->yx = r2a->yx - r2b->yx;
	r2->yy = r2a->yy - r2b->yy;
	r2->yz = r2a->yz - r2b->yz;
	r2->zx = r2a->zx - r2b->zx;
	r2->zy = r2a->zy - r2b->zy;
	r2->zz = r2a->zz - r2b->zz;
}

void real_tensor_add_r2_r2r2( Rank2R3* r2,
		Rank2R3* r2a, 
		Rank2R3* r2b) {
	r2->xx = r2a->xx + r2b->xx;
	r2->xy = r2a->xy + r2b->xy;
	r2->xz = r2a->xz + r2b->xz;
	r2->yx = r2a->yx + r2b->yx;
	r2->yy = r2a->yy + r2b->yy;
	r2->yz = r2a->yz + r2b->yz;
	r2->zx = r2a->zx + r2b->zx;
	r2->zy = r2a->zy + r2b->zy;
	r2->zz = r2a->zz + r2b->zz;
}

void real_tensor_zero_r2( Rank2R3* r2) {
	r2->xx = 0.;
	r2->xy = 0.;
	r2->xz = 0.;
	r2->yx = 0.;
	r2->yy = 0.;
	r2->yz = 0.;
	r2->zx = 0.;
	r2->zy = 0.;
	r2->zz = 0.;
}
real real_tensor_sum_offdig_r2( Rank2R3* r2) {
	return r2->xy +r2->yz + r2->zx +    r2->yx + r2->zy + r2->zx;
}
real real_tensor_sum_dig_r2( Rank2R3* r2) {
	return r2->xx + r2->yy + r2->zz;
}
void real_tensor_zero_r1( VecR3* r1) {
	r1->x =0.;
	r1->y =0.;
	r1->z =0.;
}
