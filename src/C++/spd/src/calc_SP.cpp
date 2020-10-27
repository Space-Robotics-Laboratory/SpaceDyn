//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
//
// Function : calc_SP( MODEL )
//            calculate position and orientation of each joint w.r.t. each joint
//
// Output : Xup[i]           
//
// s.abiko [2007.5]
//
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/rot.h"
#include "../include/spn.h"
#include "../include/spd.h"

void calc_SP( MODEL &m )
{
  int i;
  double* Xj;
  double* tmp;

  Xj = new Matrix6;
  tmp = new Vector3;

  // ---- calculate Xup ----
  // --- base ---
  XrotAA( m.A0, m.Xup[0]);
 
  // --- multi-body system ---
  if( m.LINKNUM != 1 ){
    for( i=1; i<m.LINKNUM; i++ ){
      matrix_Z(6,6,Xj);
      matrix_Z(3,1,tmp);
      
      if( m.J_type[i] != 1 ){ // revolute joint
	  Xrotz( m.q[i], Xj );
	  matrix_trans( 6, 6, Xj, Xj );
      }
      else{                // prismatic joint
	tmp[2] = m.q[i]; 
	Xtrans( tmp, Xj );
      }
      matrix_mult( 6, 6, 6, m.Xsp[i], Xj, m.Xup[i] );
    }    
  }

  delete [] Xj;
  delete [] tmp;
}

// === EOF === 
