//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
//
// Function : f_kin_j( MODEL )
//            calculate position and orientation of each link w.r.t. Inertia Frame
//
//
// Output : ORI_j[i] 3x3 for Link[i]
//          POS_j[i] 3x1 for Link[i]
//
// s.abiko [2007.5]
//
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/rot.h"
#include "../include/spn.h"
#include "../include/spd.h"


//void f_kin_j(MODEL &m, double *ans_oj, double *ans_pj){
void f_kin_j(MODEL &m){
  
  int i, j, k;
  double *tmp, *tmp1;

  tmp = new Matrix3;
  tmp1 = new Matrix3;

  matrix_Z(3,3,tmp);
  matrix_Z(3,3,tmp1);

  if(m.LINKNUM != 1){
    calc_Xup_I( m );

    for(i=1;i<m.LINKNUM;i++){

      for(j=0;j<3;j++){
	for(k=0;k<3;k++){
	  m.ORI_j[i][3*j+k] = m.Xup_I[i][6*j+k]; // copy term for orientation
	  tmp1[3*j+k] = m.Xup_I[i][6*(j+3)+k]; // copy term for position
	}
      }
      matrix_trans(3,3,m.ORI_j[i], tmp); //ORI_j transpose
      
      matrix_mult(3,3,3, tmp1, tmp, tmp); //POS_j
      i_tilde(tmp,m.POS_j[i]);

    }
  }
  else{
    cerr << "The system is a rigid body system.";
  }
  
  delete [] tmp;
  delete [] tmp1;

}

// === EOF ===
