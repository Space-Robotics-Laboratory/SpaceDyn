#include <stdio.h>

#include "matrix.h"
#include "vector.h"

int main()
{
    int i,m,n;
    
    m = 3;
    n = 2;
    
    double *a;
     // *a2,*b,*c;
//     double *v,*v2;
    double *l,*u;
    
    
    a = matrix_get( m, n );
    l = matrix_get( m, n );
    u = matrix_get( m, n );
    
    
    for(i=0;i<m*n;i++)
	a[i] = i;

//     a[0] = 100;
//     a[1] = 30;
//     a[3] = 1;
//     a[4] = 1;
//     a[5] = 12;
    
    matrix_print( m, n, a );
    
    matrix_LU( m, a, l, u );
     double d;
     d = matrix_det( m, a );
     printf("%lf\n", d );
     matrix_print( m, n,  l );
     matrix_print( m, n,  u );
    
    return 0;
}
