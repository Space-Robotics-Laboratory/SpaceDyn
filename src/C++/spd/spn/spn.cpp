//_/_/_/_/ Spatial Notation library _/_/_/_//

#include <iostream>
#include <cstdlib>
#include <cmath>
using namespace std;

#include <gsl/gsl_linalg.h>
#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/common.h"
#include "../include/rot.h"
#include "../include/spn.h"

#define _F(i,j) cF[(i)*6+j]
#define _V(i,j) cV[(i)*6+j]
#define _A(i,j) A[(i)*6+j]
#define _M(i,j) M[(i)*6+j]
#define _I(i,j) I[(i)*3+j]

//_/_/_/ crossF /_/_/_// */
void crossF( double *f, double *cF )
{
  int i;
  double *cV; 
   
  cV = matrix_get(6,6); // initialization
  
  crossM( f, cV );
  for(i=0;i<6*6;i++)
    cV[i] = -1*cV[i];
  
  matrix_trans( 6, 6, cV, cF );
  
  delete [] cV;
  //vcross = -crossM(v)'; 
 } 

//_/_/_/ crossM /_/_/_//
void crossM( double *v, double *cV )
{
  matrix_Z( 6, 6, cV );

  _V(0,1) = -v[5];
  _V(0,2) =  v[4];
  _V(0,4) = -v[2];
  _V(0,5) =  v[1];

  _V(1,0) =  v[5];
  _V(1,2) = -v[3];
  _V(1,3) =  v[2];
  _V(1,5) = -v[0];

  _V(2,0) = -v[4];
  _V(2,1) =  v[3];
  _V(2,3) = -v[1];
  _V(2,4) =  v[0];

  _V(3,4) = -v[5];
  _V(3,5) =  v[4];
  
  _V(4,3) =  v[5];
  _V(4,5) = -v[3];

  _V(5,3) = -v[4];
  _V(5,4) =  v[3];

  /* cV = [  0    -v(6)  v(5)   0    -v(3)  v(2) ;
    	     v(6)  0    -v(4)   v(3)  0    -v(1) ;
	    -v(5)  v(4)  0     -v(2)  v(1)  0    ;
	     0     0     0      0    -v(6)  v(5) ;
	     0     0     0      v(6)  0    -v(4) ;
	     0     0     0     -v(5)  v(4)  0    ]; 
  */
}

/* //_/_/_/ spI /_/_/_// */
void spI ( double m, double *r, double *I, double *M )
{
  
  matrix_Z( 6, 6, M );
  
  //  m*eye(3)
  _M(0,0) = m;
  _M(1,1) = m;
  _M(2,2) = m;
  
  // m*tilde(c)'
  _M(0,4) =  m*r[2];
  _M(0,5) = -m*r[1];
  _M(1,3) = -m*r[2];
  _M(1,5) =  m*r[0];
  _M(2,3) =  m*r[1];
  _M(2,4) = -m*r[0];

  // m*tilde(c)
  _M(3,1) = -m*r[2];
  _M(3,2) =  m*r[1];
  _M(4,0) =  m*r[2];
  _M(4,2) = -m*r[0];
  _M(5,0) = -m*r[1];
  _M(5,1) =  m*r[0];


  // I+m*tilde(c)*tilde(c)'
  // m+tilde(c)*tilde(c)'
  _M(3,3) = I[0] + m*( r[1]*r[1] + r[2]*r[2] );
  _M(3,4) = I[1] - m*r[0]*r[1];
  _M(3,5) = I[2] - m*r[0]*r[2];
  _M(4,3) = I[3] - m*r[0]*r[1];
  _M(4,4) = I[4] + m*( r[0]*r[0] + r[2]*r[2] );
  _M(4,5) = I[5] - m*r[1]*r[2];
  _M(5,3) = I[6] - m*r[0]*r[2];
  _M(5,4) = I[7] - m*r[1]*r[2];
  _M(5,5) = I[8] + m*( r[0]*r[0] + r[1]*r[1] );

  /* RBmci  MF6 rigid-body inertia tensor for a body with
     mass m, centre of mass at c, and (3x3) rotational inertia about CoM of I.

     rbi = [ m*eye(3), m*tilde(c)'; 
             m*tilde(c), I+m*tilde(c)*tilde(c)' ]; */
}

/* //_/_/_/ XrotBB /_/_/_// */
void XrotBB( double *Q, double *A )
{
  int i,j;
  double *R ,*C;
  
  C = matrix_get(3,3);
  R = matrix_get(3,3);
  
  matrix_Z(6,6,A);
  
  rpy2dc(Q, C);
  matrix_trans( 3, 3, C, R );
  
  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      _A(i,j) = R[3*i+j];
    }
  }

  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      _A(i+3,j+3) = R[3*i+j];
    }
  }

  /* A = [    R       zeros(3,3);
           zeros(3,3) R ] */
  delete [] R;
  delete [] C;
  
}

/* //_/_/_/ XrotAA /_/_/_// */
void XrotAA( double *R, double *A )
{
  int i,j;

  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      _A(i,j) = R[3*i+j];
    }
  }

  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      _A(i+3,j+3) = R[3*i+j];
    }
  }


  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      _A(i,j+3) = 0;
    }
  }


  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      _A(i+3,j) = 0;
    }
  }

  /* A = [    R       zeros(3,3);
           zeros(3,3) R ] */
}

/* //_/_/_/ Xtrans /_/_/_// */
void Xtrans( double *r, double *A )
{

  matrix_I( 6, A );

  _A(3,1) = -r[2];
  _A(3,2) =  r[1];
  _A(4,0) =  r[2];
  _A(4,2) = -r[0];
  _A(5,0) = -r[1];
  _A(5,1) =  r[0];

  /* Xtrans  MM6 coordinate transform from 3D translation vector.

  X = [  1     0     0    0  0  0 ;
         0     1     0    0  0  0 ;
         0     0     1    0  0  0 ;
         0    -r(3)  r(2) 1  0  0 ;
         r(3)  0    -r(1) 0  1  0 ;
        -r(2)  r(1)  0    0  0  1 ]; */
}

/* //_/_/_/ Xrotx /_/_/_// */
void Xrotx( double h, double *A )
{
  double c, s;

  c = cos(h);
  s = sin(h);

  _A(0,0) =  1.0;
  _A(1,1) =  c;
  _A(1,2) =  s;
  _A(2,1) = -s;
  _A(2,2) =  c;
  _A(3,3) =  1.0;
  _A(4,4) =  c;
  _A(4,5) =  s;
  _A(5,4) = -s;
  _A(5,5) =  c;
  
  /* Xrotx  MM6 coordinate transform from X-axis rotation.
     c = cos(h);  s = sin(h);

     X = [ 1  0  0  0  0  0 ;
           0  c  s  0  0  0 ;
           0 -s  c  0  0  0 ;
           0  0  0  1  0  0 ;
           0  0  0  0  c  s ;
           0  0  0  0 -s  c ]; */
  
}

/* //_/_/_/ Xroty /_/_/_// */
void Xroty( double h, double *A )
{
  double c, s;

  c = cos(h);
  s = sin(h);

  _A(0,0) =  c;
  _A(0,2) = -s;
  _A(1,1) =  1.0;
  _A(2,0) =  s;
  _A(2,2) =  c;
  _A(3,3) =  c;
  _A(3,5) = -s;
  _A(4,4) =  1.0;
  _A(5,3) =  s;
  _A(5,5) =  c;
  
  /* Xroty  MM6 coordinate transform from Y-axis rotation.
  c = cos(h);  s = sin(h);

  X = [ c  0 -s  0  0  0 ;
        0  1  0  0  0  0 ;
        s  0  c  0  0  0 ;
        0  0  0  c  0 -s ;
        0  0  0  0  1  0 ;
        0  0  0  s  0  c ]; */

}

/* //_/_/_/ Xrotz /_/_/_// */
void Xrotz( double h, double *A )
{
  double c, s;

  c = cos(h);
  s = sin(h);

  _A(0,0) =  c;
  _A(0,1) =  s;
  _A(1,0) = -s;
  _A(1,1) =  c;
  _A(2,2) =  1.0;
  _A(3,3) =  c;
  _A(3,4) =  s;
  _A(4,3) = -s;
  _A(4,4) =  c;
  _A(5,5) =  1.0;

  /* Xrotz  MM6 coordinate transform from Z-axis rotation.
  c = cos(h);  s = sin(h);

  X = [  c  s  0  0  0  0 ;
        -s  c  0  0  0  0 ;
         0  0  1  0  0  0 ;
         0  0  0  c  s  0 ;
         0  0  0 -s  c  0 ;
         0  0  0  0  0  1 ]; */
}

/* //_/_/_/ i_tilde /_/_/_// */
void i_tilde( double *A, double *B )
{
  matrix_Z(3,1,B);

  B[0] = A[7]; // (3,2)
  B[1] = A[2]; // (1,3)
  B[2] = A[3]; // (2,1)
  
}
// --- EOF ---
