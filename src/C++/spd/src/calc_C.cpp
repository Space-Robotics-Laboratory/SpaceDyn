//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
//
// Function : calc_C( MODEL, Gravity, CC )
//            numerically calculate non-linear velocity dependent term 
//
// Output : CC ( 6 + (LINKNUM-1) ) x 1
//
// s.abiko [2007.5]
//
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/rot.h"
#include "../include/spn.h"
#include "../include/spd.h"

void calc_C( MODEL &m, double *Gravity, double *CC){
  
  int i;

  double *vd0_tmp, *wd0_tmp, *qdd_tmp;
  double *Fe_tmp[m.LINKNUM], *Te_tmp[m.LINKNUM];

  vd0_tmp = matrix_get(3,1);
  wd0_tmp = matrix_get(3,1);
  qdd_tmp = matrix_get( m.LINKNUM, 1 );
  
  for(i=0;i<m.LINKNUM;i++){
    Fe_tmp[i] = matrix_get(3,1);
    Te_tmp[i] = matrix_get(3,1);
  }

  // copy the original data into tmp values
  matrix_cpy( 3, 1, m.vd0, vd0_tmp);
  matrix_cpy( 3, 1, m.wd0, wd0_tmp);
  matrix_cpy( m.LINKNUM, 1, m.qdd, qdd_tmp);
  for(i=1;i<m.LINKNUM;i++){
    matrix_cpy( 3, 1, m.Fe[i], Fe_tmp[i]);
    matrix_cpy( 3, 1, m.Te[i], Te_tmp[i]);
  }
  
  // set the original data "0"
  matrix_Z( 3, 1, m.vd0 );
  matrix_Z( 3, 1, m.wd0 );
  matrix_Z( m.LINKNUM, 1, m.qdd );
  for(i=1;i<m.LINKNUM;i++){
    matrix_Z( 3, 1, m.Fe[i]);
    matrix_Z( 3, 1, m.Te[i]);
  }

  i_dyn(m, Gravity, CC);

  // reset the original data
  matrix_cpy( 3, 1, vd0_tmp, m.vd0 );
  matrix_cpy( 3, 1, wd0_tmp, m.wd0 );
  matrix_cpy( m.LINKNUM, 1, qdd_tmp, m.qdd );
  for(i=1;i<m.LINKNUM;i++){
    matrix_cpy( 3, 1, Fe_tmp[i], m.Fe[i]);
    matrix_cpy( 3, 1, Te_tmp[i], m.Te[i]);
  }

  delete [] vd0_tmp;
  delete [] wd0_tmp;
  delete [] qdd_tmp;

  for(i=0;i<m.LINKNUM;i++){
    delete [] Fe_tmp[i];
    delete [] Te_tmp[i];
  }
  
}

// === EOF ===
