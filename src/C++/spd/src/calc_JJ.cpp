//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
//
// Function : calc_JJ( MODEL &m, double **ans )
//            calculate jacobian matrix for each joint w.r.t. Inertia Frame
//            
//            ans = [i][6x( LINKNUM-1 )]
//                = 6 x (LINKNUM-1) Jacobian matrix of Link[i]
//            
// s.abiko [2007.5]
//
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/rot.h"
#include "../include/spn.h"
#include "../include/spd.h"

void calc_JJ(MODEL &m, double **ans){
 
  int i, j, k;
  double *tmp, *tmp1, *tmp2, *tmp3;
  double *I_AA_T, *Xup_k_T;
  Matrix6 *I_AA_i; 

  double ***Xup_k;

  tmp  = new Matrix6;
  tmp1 = new Matrix3;
  tmp2 = new Matrix6;
  tmp3 = new Vector6;

  I_AA_i = new Matrix6[m.LINKNUM-1];
  I_AA_T = new Matrix6;

  Xup_k_T = new Matrix6;


  Xup_k = new double**[m.LINKNUM-1];
  for( i=0; i<m.LINKNUM-1; i++ )
    Xup_k[i] = new double*[m.LINKNUM-1];
  for( i=0; i<m.LINKNUM-1; i++ ){
    for( j=0; i<m.LINKNUM-1; i++ ){
      Xup_k[i][j] = new Matrix6;
      matrix_Z( 6, 6, Xup_k[i][j] );
    }
  }

  // if multibody
  if(m.LINKNUM != 1){

    calc_Xup_I( m ); // w.r.t. Inertia frame

    for( i=1; i<m.LINKNUM; i++){

      for(j=0;j<3;j++){
	for(k=0;k<3;k++){
	  tmp1[3*j+k] = m.Xup_I[i][(6*j)+k];
	}
      }
      // set I_AA_i
      XrotAA(tmp1, I_AA_i[i-1]);
    }


    // initialization of JJ
    for( i=1; i<m.LINKNUM; i++)
      matrix_Z(6,(m.LINKNUM-1), ans[i]);

    // itialization of Xup_k
    Xup_k = new double**[m.LINKNUM-1];
    for( i=0; i<m.LINKNUM-1; i++ )
      Xup_k[i] = new double*[m.LINKNUM-1];
    for( i=0; i<m.LINKNUM-1; i++ ){
      for( j=0; j<m.LINKNUM-1; j++ ){     
	Xup_k[i][j] = new Matrix6;
	matrix_Z(6, 6, Xup_k[i][j]);
      }
    }

    for(i=1; i<m.LINKNUM; i++){
      j = i;
      
      matrix_I( 6, tmp );
      
      while( j > 0 ){

	matrix_mult(6,6,6, tmp, m.jXc[i], Xup_k[i-1][j-1]);
 	matrix_mult(6,6,6, tmp, m.Xup[j], tmp);

	j = m.BB[j];
      }
      
      for( j=1; j<(i+1); j++ ){
	
	matrix_trans( 6, 6, I_AA_i[i-1], I_AA_T );
	matrix_mult( 6, 6, 6, Xup_k[i-1][j-1], I_AA_T, Xup_k[i-1][j-1] );
	matrix_trans( 6, 6, Xup_k[i-1][j-1], Xup_k_T );
	
 	matrix_mult( 6, 6, 1, Xup_k_T, m.S[j], tmp3 );
	
	for( k=0; k<6; k++ )
	  ans[i][(j-1)+(m.LINKNUM-1)*k] = tmp3[k];

      }
      //matrix_print(6,m.LINKNUM-1, ans[i]);
    }
  }
    
    delete [] tmp;
    delete [] tmp1;
    delete [] tmp2;
    delete [] tmp3;
    
    delete [] I_AA_i;
    delete [] I_AA_T;
    
    delete [] Xup_k;
    delete [] Xup_k_T;
}

// === EOF === 
