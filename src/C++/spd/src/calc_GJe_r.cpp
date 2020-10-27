//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//
// Function : calc_GJe_r( MODEL, e_num, ans )
//            calculate generalized jacobian matrix for selected end-effector in the recursive calculation
//
//            ans = (6 x LINKNUM-1) matrix
//
// s.abiko [2008.3]
//
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/rot.h"
#include "../include/spn.h"
#include "../include/spd.h"


void calc_GJe_r(MODEL &m, int e_num, double *GJe){

  matrix_Z( 6, m.LINKNUM-1, GJe );      

  int i, j, k;

  Matrix6 *IC;
  IC = new Matrix6[m.LINKNUM];
  
  double *tmp, *tmp2, *tmp3, *tmp4;
  double *Xup_T, *Xup_AT;
  double *IC_i;
  double *IC_s;
  double *GJb_tmp, *GJe_tmp;
  double *tmp_ORI;
  int *joints;
  
  int p = 0;
  int j_num = 0;
  
  // initialization
  IC_i = matrix_get(6, 6);
  IC_s = matrix_get(6, 1);
  Xup_T = matrix_get(6, 6);
  Xup_AT = matrix_get(6, 6);
  tmp = matrix_get(6, 6);
  tmp2 = matrix_get(6, 1);
  tmp3 = matrix_get(6, 6);
  tmp4 = matrix_get(6, m.LINKNUM-1);
  GJb_tmp = matrix_get(6, m.LINKNUM-1);
  GJe_tmp = matrix_get(6, m.LINKNUM-1);
  tmp_ORI = matrix_get(6, 6);
  
  joints = new int[m.LINKNUM]; // joints[0] invalid

  if( m.LINKNUM != 1 ){
    for( i=1; i<m.LINKNUM; i++ ){
      if(m.EE[i] == e_num){
	p = i;
	j = i;
	
	while(j != 0){ // to size memory
	  j_num++;
	  j = m.BB[j];
	}
	
	j = i;
	for( k=j_num; k>0; k--){
	  joints[k] = j;
	  j= m.BB[j];
	}
      }
    }


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
      	GJb_tmp[(m.LINKNUM-1)*k + (i-1)] = IC_s[k];
      }
   }
  }
  
  matrix_inv( 6, IC[0], tmp3 );
  matrix_mult( 6, 6, m.LINKNUM-1, tmp3, GJb_tmp, GJb_tmp );

  matrix_I( 6, Xup_AT );
  
  for( i=j_num;i>0;i-- ){
  
 	matrix_trans( 6, 6, m.Xup[joints[i]], tmp );
 	matrix_mult( 6, 6, 6, Xup_AT, tmp, Xup_AT );
	
	matrix_I( 6, Xup_T );
	if( i != j_num ){
		for( k=joints[i+1];k<=j_num;k++){
		 	matrix_trans( 6, 6, m.Xup[joints[k]], tmp );
			matrix_mult( 6, 6, 6, tmp, Xup_T, Xup_T );
		}
	}
	
	matrix_mult( 6, 6, 1, Xup_T, m.S[joints[i]], tmp2 );
	for( k=0;k<6;k++){
		GJe_tmp[ (m.LINKNUM-1)*k + joints[i]-1 ] = tmp2[k];
	}
  }  
  matrix_mult(6, 6, m.LINKNUM-1, Xup_AT, GJb_tmp, tmp4 );
  matrix_add( 6, m.LINKNUM-1, GJe_tmp, tmp4, GJe );
 
  matrix_trans( 6, 6, m.jXe[p], tmp );
  matrix_mult( 6, 6, m.LINKNUM-1, tmp, GJe, GJe );
  
  f_kin_e( m, e_num );
 
  matrix_Z( 6, 6, tmp_ORI );
  matrix_cpy_sub( 6, 6, 1, 3, 1, 3, m.ORI_e[e_num], tmp_ORI );
  matrix_cpy_sub( 6, 6, 4, 6, 4, 6, m.ORI_e[e_num], tmp_ORI );
  
  matrix_mult( 6, 6, m.LINKNUM-1, tmp_ORI, GJe, GJe );
 
  }

   // ---- clear memories ----
  delete [] IC;

  delete [] tmp;
  delete [] tmp2;
  delete [] tmp3;
  delete [] tmp4;
  delete [] Xup_T;
  delete [] IC_i;
  delete [] IC_s;

  delete [] Xup_AT;
  delete [] GJb_tmp;
  delete [] GJe_tmp;
  delete [] tmp_ORI;

  delete [] joints;

}

// === EOF ===
