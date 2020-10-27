//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//
// Function : calc_Je( MODEL, e_num, ans )
//            calculate jacobian matrix for selected end-effector
//
//            ans = (6 x LINKNUM-1) matrix
//
// s.abiko [2007.5]
//
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/rot.h"
#include "../include/spn.h"
#include "../include/spd.h"

void calc_Je(MODEL &m, int e_num, double *ans){
  
  int i, j, k;
  int *joints;
  double *I_AA_e;
  double *JJ;

  Matrix6 Xup_k_T;

  double ***Xup_k; // don't confuse due to many pointers
  double *tmp, *tmp1, *tmp2, *tmp3;

  int p = 0;
  int j_num = 0;
  int j_num_k = 0;

  I_AA_e = new Matrix6;
  matrix_Z( 6, 6, I_AA_e );

  tmp  = new Matrix6;
  tmp1 = new Matrix6;
  tmp2 = new Vector6;
  tmp3 = new Matrix3;

  JJ = new double[6*(m.LINKNUM-1)];

  joints = new int[m.LINKNUM];

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

    // itialization of Xup_k
    Xup_k = new double**[j_num];

    for( i=0; i<j_num; i++ )
      Xup_k[i] = new double*[j_num];
    for( i=0; i<j_num; i++ ){
      for( j=0; j<j_num; j++ ){     
	Xup_k[i][j] = new Matrix6;
	matrix_Z(6, 6, Xup_k[i][j]);
      }
    }
    // initialization of JJ
    matrix_Z(6, m.LINKNUM-1, JJ);    

    calc_Xup_I( m ); // w.r.t. Inertia frame
    matrix_mult( 6, 6, 6, m.Xup_I[p], m.jXe[p], tmp );
    // set I_AA_e
    for( i=0; i<3; i++)
      for( j=0; j<3; j++)
	tmp3[3*i+j] = tmp[6*i+j];

    XrotAA( tmp3, I_AA_e );  
    
    matrix_I( 6, tmp );
    j = p;
    j_num_k = j_num;
    while( j > 0 ){
      matrix_mult( 6, 6, 6, tmp, m.jXe[p], Xup_k[j_num-1][j_num_k-1]);
    //  matrix_mult( 6, 6, 6, tmp, m.Xup[j], tmp );
      matrix_mult( 6, 6, 6, m.Xup[j], tmp, tmp );
      
      j = m.BB[j];
      j_num_k--;
    }
    
    for( j=0; j<j_num; j++ ){	  
      matrix_trans( 6, 6, Xup_k[j_num-1][j], Xup_k_T );	  
      matrix_mult( 6, 6, 1, Xup_k_T, m.S[joints[j+1]], tmp2 );

      for( k=0;k<6;k++)
	JJ[joints[j+1] + k*(m.LINKNUM-1)-1 ] = tmp2[k];
    }                
    
    matrix_mult( 6, 6, m.LINKNUM-1, I_AA_e , JJ, ans);
//    matrix_print( 6, m.LINKNUM-1, ans);
    
  }
  else{
    cerr << "!!The system has only one regid body!!" << endl;
  }

  delete [] I_AA_e;
  delete [] joints;
  delete [] JJ;

  delete [] tmp;
  delete [] tmp1;
  delete [] tmp2;
  delete [] tmp3; 
  
  for( i=0; i<j_num; i++ ){
    for( j=0; j<j_num; j++ )
      	delete [] Xup_k[i][j];
  }
  for(i=0; i<j_num; i++ )
	delete [] Xup_k[i];
  delete [] Xup_k;

}
// === EOF ===
