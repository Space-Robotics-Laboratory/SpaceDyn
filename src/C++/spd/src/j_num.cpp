//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
//
// Function : j_num( MODEL, num_e, joints, j_num, e_joint )
//            calculate the corresponding joints for the selected end-effector 
//
// s.abiko [2007.5]
//
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/rot.h"
#include "../include/spn.h"
#include "../include/spd.h"

// --- joint connectivity in terms of each arm ----
void j_num(MODEL &m, int num_e, int *joints, int j_num, int e_joint){

  int i, j, k;

  j = 0;
  j_num = 0;
  e_joint = 0;

  for( i=1; i<m.LINKNUM; i++ ){
    if(m.EE[i] == num_e){
      e_joint = i;
      j = i;

      while(j != 0){ // to size memory
	j_num++;
	j = m.BB[j];
	cout << j_num << endl;
      }

      j = i;
      for( k=j_num; k>0; k--){
	joints[k] = j;
   	cout << k <<" " << joints[k] << endl;
	j= m.BB[j];
      }
    }
  }
}

// === EOF ===
