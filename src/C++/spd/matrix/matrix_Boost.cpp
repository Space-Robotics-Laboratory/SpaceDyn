//_/_/_/_/ matrix library _/_/_/_//

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_linalg.h>

#include "matrix.h"
#include "vector.h"

static void gsl_matrix_svd( int m, int n, double *a, double *u, double *s, double *v );
static void gsl_matrix_pinv( int m, int n, double *a, double *u );


//_/_/_/ get matirx /_/_/_//
double *matrix_get( int m, int n )
{
    return ( double * )new double[m * n];
//    return new double[m * n];
}

//_/_/_/ get matirx /_/_/_//
/* === old c expression ===
double *matrix_get( int m, int n )
{
    return ( double * ) malloc( sizeof( double ) * m * n );
}
*/

//_/_/_/ copy matrix /_/_/_//
void matrix_cpy( int m, int n, double *a, double *u )
{
    int max, i;
    
    max = m * n ;
    for ( i = 0 ; i < max ; i++ )
    {
	u[ i ] = a[ i ];
    }
}

//_/_/_/ set unit matrix /_/_/_//
void matrix_I( int n, double *l )
{
    int i, j ;
    
    for ( i = 0 ; i < n ; i++ )
    {
	for ( j = 0 ; j < n ; j++ )
	{
	    l[ i*n + j ] = 0;
	}
    }
    for ( i = 0 ; i < n ; i++ )
    {
	l[ i*n + i ] = 1;
    }
}

//_/_/_/ set zero matrix /_/_/_//
void matrix_Z( int m, int n, double *l )
{
    int max, i; 
    
    max = n * m ;
  
    for ( i = 0 ; i < max ; i++ )
    {
	l[ i ] = 0.0;
    }
}


//_/_/ print matrix _/_//
void matrix_print( int m, int n, double *a )
{
    int	i,j;
    
    for ( i=0 ; i<m ; i++ )
    {
	for ( j=0 ; j<n ; j++ )
	    printf("%5f ", a[ i*n + j ]);
	printf("\n");
    } 
    printf("\n");
}


//_/_/ multiplication mxn x nxl _/_//
void matrix_mult( int m, int n, int l, double *a ,double *b, double *c )
{
    int	i,j,k;
    double s;
    
    double *A = matrix_get( m, n );
    double *B = matrix_get( n, l );
    
    matrix_cpy( m, n, a, A );
    matrix_cpy( n, l, b, B );
    
    for ( i=0 ; i<m ; i++ )
    {
	for ( j=0 ; j<l ; j++ )
	{
	    s=0.0;
	    
	    for ( k=0 ; k<n ; k++ )
		s += A[ i*n + k ] * B[ k*l + j ];
	    
	    c[ i*l + j ] = s;
	}
    }
    
    free(A);
    free(B);
}


//_/_/ addition _/_//
void matrix_add( int m, int n, double *a ,double *b, double *c )
{
    int	i,j;
    
    double *A = matrix_get( m, n );
    double *B = matrix_get( m, n );
    
    matrix_cpy( m, n, a, A );
    matrix_cpy( m, n, b, B );
    
    for ( i=0 ; i<m ; i++ )
    {
	for ( j=0 ; j<n ; j++ )
	{
	    c[ i*n + j ] = A[ i*n + j ] + B[ i*n + j ];
	}
    }
    
    free(A);
    free(B);
}


//_/_/ subtraction _/_//
void matrix_sub( int m, int n, double *a ,double *b, double *c )
{
    int	i,j;
    
    double *A = matrix_get( m, n );
    double *B = matrix_get( m, n );
    
    matrix_cpy( m, n, a, A );
    matrix_cpy( m, n, b, B );
    
    for ( i=0 ; i<m ; i++ )
    {
	for ( j=0 ; j<n ; j++ )
	{
	    c[ i*n + j ] = A[ i*n + j ] - B[ i*n + j ];
	}
    }
    
    free(A);
    free(B);
}


//_/_/ transposed matrix _/_//
void matrix_trans( int m, int n, double *a ,double *b )
{
    int	i,j;
  
    double *A = matrix_get( m, n );
    
    matrix_cpy( m, n, a, A );
    
    for( i=0 ; i<m ; i++ )
    {
	for( j=0 ; j<n ; j++ )
	{
	    b[ j*m + i ] = A[ i*n + j ];
	}
    }
    
    free(A);
}


//_/_/ round matrix _/_//
void matrix_round( int m, int n, int r, double *a )
{
    int	i;
    double p;
    
    p = pow( (double)10, r );
    
    for( i=0 ; i<m*n ; i++ )
	a[i] = ( rint( a[i]*p ) )/p;
    
}

//_/_/ scale _/_//
void matrix_scale( int m, int n, double s, double *a, double *b )
{
    int	i;
    
    double *A = matrix_get(m,n);
    
    matrix_cpy(m,n,a,A);
    
    for( i=0 ; i<m*n ; i++ )
	b[i] = s*A[i];
    
    free(A);
}

/* //_/_/ extract submatrix /_/_/_/ */
/* void matrix_submat_ext( int m, int n, ) */
/* { */

/* } */

/* //_/_/ copy in submatrix /_/_/_/ */
/* void matrix_submat_copy( int m, int n, ) */
/* { */

/* } */
//_/_/   Matrix Multiplication   _/_//
//_/_/ using pointer calculation _/_//
void mxm1( int n, double *a , double *b , double *c )
{
    int    i,j,k,kk ;
    double s ;
    double *ai,*bj,*cj ;
    
    for ( j=0 , bj=b , cj=c ; j<n ; ++j , ++bj , ++cj )
    {
	for ( i=0, ai=a ; i<n ; ++i , ai+=n )
	{
	    s=0.0;
	    for ( k=0 , kk=0 ; k<n ; ++k , kk+=n )
		s+=ai[k]*bj[kk];
	    cj[i*n]=s;
	}
    }
}

//_/_/   Inverse Matrix   _/_//
void matrix_inv( int n, double *a, double *b)
{
    int i,j,k;
    double p,q;
    
    double *A = matrix_get( n, n );
    matrix_cpy( n, n, a, A );
    matrix_I( n, b );
    
    for ( k=0 ; k<n ; ++k )
    {
	p=A[k*n+k];
	
	for ( j=0 ; j<n ; ++j )
	{
	    b[k*n+j] /= p;
	    A[k*n+j] /= p;
	}
	
	for ( i=0 ; i<n ; ++i )
	{
	    if ( i!= k )
	    {
		q = A[ i*n + k ];
		
		for ( j=0 ; j<n ; ++j ){
		    A[ i*n + j ] -= q*A[ k*n + j ];
		    b[ i*n + j ] -= q*b[ k*n + j ];
		}
	    }
	}
    }
    free(A);
}


//_/_/   Pseudoinverse Matrix   _/_//
void matrix_pinv( int m, int n, double *a, double *b )
{
    
    double *A = matrix_get(m,n);
    matrix_cpy(m,n,a,A);
    
    if(m>=n)
    {
	gsl_matrix_pinv( m, n, A, b );
    }
    else if(m<n)
    {
	matrix_trans(m,n,A,A);
	gsl_matrix_pinv(n,m,A,A);
	matrix_trans(m,n,A,b);
    }
    
    free(A);
    
}


//_/_/_/ LU decomposition /_/_/_//
void matrix_LU( int n, double *a, double *l, double *u )
{
    int p, i, j;
    double c;
    
    matrix_cpy( n, n, (double *)a, (double *)u ); /* 行列 Uに 行列 Aをコピーする */
    matrix_I( n, (double *)l );				  /* 行列 Lを単位行列にする		 */
    
    for ( p = 0 ; p < n ; p++ ) /* 掃き出し操作を行いながら */
    {							/* 行列を L,Uに分解する	   */
	for ( i = p + 1 ; i < n ; i++ ) 
	{
	    c = u[ i*n + p ] / u[ p*n + p ];
	    l[ i*n + p ] = c;
	    for ( j = p ; j < n ; j++ )
	    {
		u[ i*n + j ] -= c * u[ p*n + j ];
	    }
	}
    }
}

//_/_/_/ return the product of diagonal elements of matrix *a  /_/_/_//
double matrix_mlting( int n, double *a )
{
    int i;
    double ret;
    
    ret = 1;
    for ( i = 0 ; i < n ; i++ )
    {
	ret *= a[ i*n + i ];
    }
    return ret;
}

//_/_/_/ get determinant /_/_/_// 
double matrix_det( int n, double *a )
{
    double *L, *U ;
    double ret;
    
    L = matrix_get( n, n );
    U = matrix_get( n, n );
    
    matrix_LU( n, a, L, U );		  /* LU分解を行う */
    ret = matrix_mlting( n, U ); /* 対角要素の積を求める */
    
    free(L);
    free(U);
    
    return ret;
}


//_/_/_/ singlar value decomposition /_/_/_/_//
void matrix_svd( int m, int n, double *a, double *u, double *s, double *v )
{
    gsl_matrix_svd( m, n, a, u, s, v );
}


void gsl_matrix_svd( int m, int n, double *a, double *u, double *s, double *v )
{
    int i, j;
    
    gsl_matrix *A = gsl_matrix_alloc (m,n);
    gsl_matrix *V = gsl_matrix_alloc (n,n);
    gsl_vector *sigma = gsl_vector_alloc (n);
    gsl_vector *work = gsl_vector_alloc (n);
    
    
    for(i=0;i<m;i++)
	for(j=0;j<n;j++)
	    gsl_matrix_set(A,i,j,a[i*n+j]); 
    
    
    gsl_linalg_SV_decomp(A, V, sigma, work);
    
    
    matrix_Z(m,n,u);
    for(i=0;i<m;i++)
	for(j=0;j<n;j++)
	    u[i*n+j] = gsl_matrix_get(A,i,j);
    
    matrix_Z(n,n,v);
    for(i=0;i<n;i++)
	for(j=0;j<n;j++)
	    v[i*n+j] = gsl_matrix_get(V,i,j); 
    
    vector_z(n,s);
    for(i=0;i<n;i++)
	s[i] = gsl_vector_get(sigma,i); 
    
    
    gsl_matrix_free(A);
    gsl_matrix_free(V);
    gsl_vector_free(sigma);
    gsl_vector_free(work);
    
}

void gsl_matrix_pinv( int m, int n, double *a, double *u )
{
    int i;
    
    double *U = matrix_get(m,n);
    double *V = matrix_get(n,n);
    double *S = vector_get(n);
    double *S2 = matrix_get(n,n);
    
    matrix_svd( m, n, a, U, S, V );
    
    matrix_Z(n,n,S2);
    for(i=0;i<n;i++)
	if(S[i]>1e-10)
	    S2[i*n+i] = 1.0/S[i];
    
    matrix_mult(n,n,n,V,S2,S2);
    matrix_trans(m,n,U,U);
    matrix_mult(n,n,m,S2,U,u);
    
    free(U);
    free(V);
    free(S);
    free(S2);
    
}

//_/_/_/ extract a row vector from a matrix /_/_/_//
void matrix_ext_row( int m, int n, int i, double *A, double *row){

  int k;

  for( k=0; k<n; k++ ){
    row[k] = A[n*(i-1)+k];
  }
}

//_/_/_/ extract a column vector from a matrix /_/_/_//
void matrix_ext_col( int m, int n, int j, double *A, double *colmun){

  int k;

  for( k=0; k<m; k++ ){
    colmun[k] = A[n*k+(j-1)];
  }
}

//_/_/_/ copy a vector into i-th row of a matrix /_/_/_//
void matrix_cpy_row( int m, int n, int i, double *row, double *A){

  int k;

  for( k=0; k<n; k++ ){
    A[n*(i-1)+k] = row[k];
  }
}

//_/_/_/ copy a column into j-th column of a matrix /_/_/_//
void matrix_cpy_col( int m, int n, int j, double *column, double *A){

  int k;

  for( k=0; k<m; k++ ){
    A[n*k+(j-1)] = column[k];
  }
}

//_/_/_/ extract a (ixj) submatrix from a matrix /_/_/_//
void matrix_ext_sub( int m, int n, int i_s, int i_f, int j_s, int j_f, double *A, double *ans ){

  int i, j;
  int row, col;

  row = i_f - (i_s-1);
  col = j_f - (j_s-1);

  for( i=(i_s-1); i<i_f; i++){
    for( j=(j_s-1); j<j_f; j++ ){
      ans[col*(i-(i_s-1))+(j-(j_s-1))] = A[n*i+j];
    }
  }
}


//_/_/_/ copy a (ixj) submatrix into a matrix /_/_/_//
void matrix_cpy_sub( int m, int n, int i_s, int i_f, int j_s, int j_f, double *A, double *ans ){

  int i, j;
  int row, col;

  row = i_f - (i_s-1);
  col = j_f - (j_s-1);

  for( i=(i_s-1); i<i_f; i++){
    for( j=(j_s-1); j<j_f; j++ ){
      ans[n*i+j] = A[col*(i-(i_s-1))+(j-(j_s-1))];
    }
  }
}
