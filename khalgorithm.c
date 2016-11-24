/*
 * =====================================================================================
 *
 *       Filename:  khalgorithm.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2016년 11월 10일 13시 57분 16초
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ph.D. Candidate KIM Hyeok (kh), ekh0324@gmail.com
 *        Company:  Konkuk University
 *
 * =====================================================================================
 */
#include <math.h>



typedef enum { SUCCESS, FAILED,} SUCCESS_TYPE;

SUCCESS_TYPE get_autocorelation_from_seq ( int width, double* seq, int max_diff, double* autoFunction, double* autoError) {
	// Init
	int count[max_diff] ;
	int ii,jj, diff ;
	double mul_seq_ii_jj;

	for ( diff = 0; diff < max_diff; diff += 1 ){
		count[diff]        = 0;
		autoFunction[diff] = 0.; 
		autoError[diff]    = 0.;
	}

// Evaluation
	
	for ( ii = 0; ii < width; ii += 1 )
		for ( jj = 0; jj <= ii; jj += 1 )
		{
			diff = ii-jj;
			if(diff<max_diff)
			{
				count[diff]++;
				mul_seq_ii_jj   = seq[ii]*seq[jj];
				autoFunction[diff] +=mul_seq_ii_jj;
				autoError[diff]    +=mul_seq_ii_jj*mul_seq_ii_jj;
			}
		}
// averaging
 	
	for ( diff = 0; diff < max_diff; diff += 1 ){
		autoFunction[diff] /=  count[diff]; 
		autoError[diff]    =  sqrt( (autoError[diff]/count[diff] 
				- autoFunction[diff]*autoFunction[diff]) /count[diff]);  
	}

	return SUCCESS ;
}

