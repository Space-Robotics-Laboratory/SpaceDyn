//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
//
// Function : calc_Xup_I( MODEL )
//            calculate position and orientation of each joint w.r.t. Inertia Frame
//
// Output : m.Xup_I[i]
//
// s.abiko [2007.5]
//
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/rot.h"
#include "../include/spn.h"
#include "../include/spd.h"

// ---- Transfer Matrix w.r.t. Inertia Frame ----
void calc_Xup_I(MODEL &m){

  int i;
  double *tmp;
  
  tmp = new Matrix6;

  Xtrans(m.POS0, tmp);
  matrix_mult( 6, 6, 6, tmp, m.Xup[0], m.Xup_I[0] );

  if( m.LINKNUM != 1 ){
    for(i=1;i<m.LINKNUM;i++)
      matrix_mult( 6, 6, 6, m.Xup_I[m.BB[i]], m.Xup[i], m.Xup_I[i]);
  }

  delete [] tmp;
}

// === EOF ===
