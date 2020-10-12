//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//
// Function : calc_GJb( MODEL, ans )
//            calculate Generalized jacobian matrix related to the base
//
//            ans = (6 x 6) matrix
//
// s.abiko [2008.3]
//
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/rot.h"
#include "../include/spn.h"
#include "../include/spd.h"


void calc_GJb( MODEL &m, double *GJb){

  matrix_Z( 6, m.LINKNUM-1, GJb );      
  
  int i, j, k;

  Matrix6 *IC;
  IC = new Matrix6[m.LINKNUM];
  
  double *tmp, *tmp2, *tmp3;
  double *Xup_T;
  double *IC_i;
  double *IC_s;

  // initialization
  IC_i = matrix_get(6, 6);
  IC_s = matrix_get(6, 1);
  Xup_T = matrix_get(6, 6);
  tmp = matrix_get(6, 6);
  tmp2 = matrix_get(6, 1);
  tmp3 = matrix_get(6, 6);
  
  // copy the Inertia Matrix in joint space
  for(i=0;i<m.LINKNUM;i++){
    matrix_cpy( 6, 6, m.Isp[i], IC[i] );
  }

  if( m.LINKNUM != 1 ){
    for( i=m.LINKNUM-1; i>0; i-- ){

      // ---- calculate the composite inertia matrix ----
      matrix_trans( 6, 6, m.Xup[i], Xup_T );
      matrix_mult( 6, 6, 6, IC[i], Xup_T, IC_i );
      matrix_mult( 6, 6, 6, m.Xup[i], IC_i, IC_i );
      matrix_add( 6, 6, IC[m.BB[i]], IC_i, IC[m.BB[i]] );
    
      // ---- calculate the generalized Jacobian for the base ----
      j = i;
      matrix_I( 6, tmp );
      while( j > 0 ){
      	matrix_mult( 6, 6, 6, m.Xup[j], tmp, tmp );
	j = m.BB[j];
      }
     
      matrix_mult( 6, 6, 1, IC[i], m.S[i], tmp2);
      matrix_mult( 6, 6, 1, tmp, tmp2, IC_s);
      matrix_scale( 6, 1, -1, IC_s, IC_s );
       
      // copy the IC_s in GJb(:,i-1)
      for(k=0;k<6;k++){
      	GJb[(m.LINKNUM-1)*k + (i-1)] = IC_s[k];
      }
   }
  }
  
  matrix_inv( 6, IC[0], tmp3 );
  matrix_mult( 6, 6, m.LINKNUM-1, tmp3, GJb, GJb );
  matrix_mult( 6, 6, m.LINKNUM-1, m.Xup[0], GJb, GJb );
  
  // ---- clear memories ----
  delete [] IC;

  delete [] tmp;
  delete [] tmp2;
  delete [] tmp3;
  delete [] Xup_T;
  delete [] IC_i;
  delete [] IC_s;

}

// === EOF ===
