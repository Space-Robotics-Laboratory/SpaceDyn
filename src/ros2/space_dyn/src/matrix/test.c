// g++ test.c libmatrix.a

#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include "vector.h"

int main()
{
    int i;
    double *a, *b, *c, *d;
    
    a = get_matrix( 3 );
    b = get_matrix_mxn( 6, 3 );
    c = get_matrix_mxn( 6, 3 );
    d = get_matrix_mxn( 3, 6 );
    
    for( i=0 ; i<3*3 ; i++ )
	a[i] = i;
    for( i=0 ; i<3*6 ; i++ )
	b[i] = 1.0;
    
    for( i=0 ; i<3*6 ; i++ )
	c[i] = i;
    //print_matrix( 3, a );
    //print_matrix_mxn( 6, 3, b );
    
//    mult_matrix_mxm_nxm( 3, 6, a, b, c );
    
    //print_matrix_mxn( 6, 3, c );
    
    trans_matrix_mxn( 6, 3, c, d );
    print_matrix_mxn( 6, 3, c );
    print_matrix_mxn( 3, 6, d );
    
}
