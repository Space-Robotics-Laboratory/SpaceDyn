//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
//
// Function : rref( int m, int n, double *A, double *ans )
//            calculate reduced row echelon form.
//            
//            ans = [size(A)]
//            
// s.abiko [2007.8]
//
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
#include <iostream>
#include <cmath>
using namespace std;

#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/rot.h"
#include "../include/spn.h"
#include "../include/spd.h"

void rref(int m, int n, double *A, double *ans)
{

  int i, j, k, l, a;
  double tmp, p;

  i = 0;
  j = 0;
  a = 0;

  while( ( i < m ) && ( j < n ) ){
    p = 0.0;
    tmp = 0.0;

    for( k=i; k<m; k++ ){
      tmp = fabs(A[k*n+j]);
      if( tmp > p ){
	p = tmp;
	a = k;
      }
    }
    
    if( p <= 1e-12 ){
      for( k=i; k<m; k++){
	A[k*n+j] = 0;
      }
      j++;
    }
    else{
      // Swap i-th and k-th rows.
      for( k=j; k<n; k++){
	tmp = A[n*i+k];
	A[n*i+k] = A[n*a+k];
	A[n*a+k] = tmp;
      }

     // Divide the pivot row by the pivot element.
      tmp = A[n*i+j];
      for( k=j; k<n; k++ )
	A[n*i+k] = A[n*i+k]/tmp;
      
      // Subtract multiples of the pivot row from all the other rows.
      for( l=0; l<i; l++ ){
	tmp = A[n*l+j];
	for( k=j; k<n; k++ ){
	  A[n*l+k] = A[n*l+k] - tmp*A[n*i+k];
	}
      }
       
      for( l=i+1; l<m; l++ ){
	tmp = A[n*l+j];
	for( k=j; k<n; k++ ){
	  A[n*l+k] = A[n*l+k] - tmp*A[n*i+k];
	}
      }
      i++;
      j++;	
    }
  }

  matrix_cpy( m, n, A, ans );
}

// === EOF ===
