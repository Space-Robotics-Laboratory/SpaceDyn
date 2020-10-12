//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
//
// Function : f_kin_e( MODEL, e_num )
//            calculate position and orientation of the selected end-effector
//
//
// Output : m.ORI_e  3x3
//        : m.POS_e  3x1 
//
// s.abiko [2007.5]
//
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/rot.h"
#include "../include/spn.h"
#include "../include/spd.h"


void f_kin_e(MODEL &m, int e_num){

  int i, j, k;
  double *I_Xup_e, *tmp, *tmp1;

  I_Xup_e = new Matrix6;
  matrix_Z(6,6,I_Xup_e);

  tmp = new Matrix3;
  tmp1 = new Matrix3;

  matrix_Z(3,3,tmp);
  matrix_Z(3,3,tmp1);

  if( m.LINKNUM != 1 ){
    calc_Xup_I(m);
    
    for( i=1; i<m.LINKNUM; i++){
      if(m.EE[i] == e_num){

	matrix_mult( 6, 6, 6, m.Xup_I[i], m.jXe[i], I_Xup_e);

	for(j=0;j<3;j++){
	  for(k=0;k<3;k++){
	    m.ORI_e[m.EE[i]][3*j+k] = I_Xup_e[6*j+k]; // copy term for orientation
	    tmp[3*j+k] = I_Xup_e[6*(j+3)+k]; // copy term for position
	  }
	}
	//matrix_print(6,6,I_Xup_e);

	matrix_trans( 3, 3, m.ORI_e[m.EE[i]], tmp1 );
	matrix_mult( 3, 3, 3, tmp, tmp1, tmp );
	i_tilde( tmp, m.POS_e[m.EE[i]] );

// 	matrix_print(3,3,m.ORI_e[m.EE[i]]);
// 	matrix_print(3,1,m.POS_e[m.EE[i]]);
      }
    }
  }
  else{
    cerr << " The model is one rigid body system." << endl;
  }

  delete [] tmp;
  delete [] tmp1;
  delete [] I_Xup_e;

}

// === EOF ===
