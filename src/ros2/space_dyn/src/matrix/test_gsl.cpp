#include <stdio.h>

#include "matrix.h"
#include "vector.h"

int
main (void)
{
    int i,j;
    /*
    double a_data[] = { 0.18, 0.60, 0.57, 0.96,
			0.41, 0.24, 0.99, 0.58,
			0.14, 0.30, 0.97, 0.66,
			0.51, 0.13, 0.19, 0.85 };
    */
    double a[] = { 0.18, 0.60, 0.57, 0.96,
		   0.41, 0.24, 0.99, 0.58,
		   0.14, 0.30, 0.97, 0.66,
		   0.14, 0.30, 0.97, 0.66 };
    
    double *u = matrix_get(4,4);
    double *s = vector_get(4);
    double *v = matrix_get(4,4);
    
    matrix_svd( 4, 4, a, u, s, v );
    
    matrix_print( 4, 4, u );
    vector_print( 4, s );
    matrix_print( 4, 4, v );
    
    matrix_pinv( 4,4, a, a );
    
    matrix_print(4,4,a);
    
    /*
    double b_data[] = { 1.0, 2.0, 3.0, 4.0 };
    
    gsl_pinv_matrix_mxn( 4,4,a_data,a );
    
    print_matrix(4,a);
    
    
    
    inv_matrix(4,a_data,a);
    print_matrix(4,a);
    */
    return 0;
    
}

