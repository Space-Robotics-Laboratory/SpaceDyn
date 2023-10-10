//_/_/_/_/ vector library _/_/_/_//

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "vector.h"
using namespace std;

//_/_/_/ get vector /_/_/_//
double *vector_get( int n )
{
    return ( double * )new double[n];
}

/*//_/_/_/ get vector /_/_/_//
=== old expression with C langurage
double *vector_get( int n )
{
    return ( double * ) malloc( sizeof( double ) * n );
}
*/
//_/_/_/ copy vector /_/_/_//
void vector_cpy( int n, double *a, double *u )
{
    int i;
    
    for ( i = 0 ; i < n ; i++ )
    {
	u[ i ] = a[ i ];
    }
}

//_/_/ print vector _/_//
void vector_print( int n, double *a )
{
    int	i;
    
    for ( i=0 ; i<n ; i++ )
        cout << a[i] << " " ;
    cout << endl;
}


//_/_/ print vector _/_//
void vector_print_int( int n, int *a )
{
    int	i;
    
    for ( i=0 ; i<n ; i++ )
        cout << a[i] << " " ;
    cout << endl;
}


//_/_/_/ addition /_/_/_//
void vector_add( int n, double *a, double *b, double *c )
{
    int i;
    
    double *A, *B;
    A = vector_get( n );
    B = vector_get( n );

    vector_cpy( n, a, A );
    vector_cpy( n, b, B );
    
    for( i=0 ; i<n ; i++ )
    {
	c[i] = a[i] + b[i];
    }
    
    delete [] A;
    delete [] B;    
    
}


//_/_/_/ subtraction /_/_/_//
void vector_sub( int n, double *a, double *b, double *c )
{
    int i;
    
    double *A, *B;
    A = vector_get( n );
    B = vector_get( n );

    vector_cpy( n, a, A );
    vector_cpy( n, b, B );
    
    for( i=0 ; i<n ; i++ )
    {
	c[i] = A[i] - B[i];
    }

    delete [] A;
    delete [] B;    
    
}


//_/_/_/ return an inner product of matrix a & b /_/_/_//
double vector_inner( int n, double *a, double *b )
{
    int i;
    double sum;
    double *tmp = vector_get( n );
    
    sum = 0.0;
    
    for( i=0 ; i<n ; i++ )
    {
	tmp[i] = a[i]*b[i];
	sum +=tmp[i];
    }
  
   delete [] tmp;
    
   return sum;
}


//_/_/_/ return a cross product of matrix a(3x1) & b(3x1) /_/_/_//
void vector_cross3( double *a, double *b, double *c )
{
	
    double *A, *B;
    A = vector_get( 3 );
    B = vector_get( 3 );

    vector_cpy( 3, a, A );
    vector_cpy( 3, b, B );
    
    c[0] = A[1]*B[2]-A[2]*B[1];
    c[1] = A[2]*B[0]-A[0]*B[2];
    c[2] = A[0]*B[1]-A[1]*B[0];
    
    delete [] A;
    delete [] B;
    
}


//_/_/_/ return a tilde matrix form vector3 /_/_/_//
void vector_tilde3( double *a, double *b )
{

    double *A;
    A = vector_get( 3 );
    vector_cpy( 3, a, A );

    b[0] = 0;
    b[1] = -A[2];
    b[2] = A[1];
    b[3] = A[2];
    b[4] = 0;
    b[5] = -A[0];
    b[6] = -A[1];
    b[7] = A[0];
    b[8] = 0;

    delete [] A;

}

void vector_z( int n, double *a )
{
    int i;
    
    for( i=0 ; i<n ; i++ )
	a[i] = 0.0;
}

void vector_i( int n, double *a )
{
    int i;
    
    for( i=0 ; i<n ; i++ )
	a[i] = 1.0;
}

void vector_scale( int n, double s, double *a, double *b )
{
    int i;
    
    double *A = vector_get(n);
    
    vector_cpy(n,a,A);
    
    for( i=0 ; i<n ; i++ )
	b[i] = s*A[i];
    
    delete [] A;
}
