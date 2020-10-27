//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
//
// Function : f_dyn_rk( MODEL )
//            Runge-Kutta Integrator for fixed step
//            calculate forward dynamics
//
// s.abiko [2007.5]
//
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/rot.h"
#include "../include/spn.h"
#include "../include/spd.h"

//#include <QDebug>
// --------------------------
void f_dyn_rk_ffg(MODEL &m, double *Gravity, double step){

  double recip = 0.166666666666666;
  double recip4 = 0.25;
  
  double *k1_q, *k1_qd, *k1_POS0, *k1_Qtn0, *k1_v0, *k1_w0;
  double *k2_q, *k2_qd, *k2_POS0, *k2_Qtn0, *k2_v0, *k2_w0;
  double *k1_q_tmp, *k1_qd_tmp, *k1_POS0_tmp, *k1_Qtn0_tmp, *k1_v0_tmp, *k1_w0_tmp;
  double *k2_q_tmp, *k2_qd_tmp, *k2_POS0_tmp, *k2_Qtn0_tmp, *k2_v0_tmp, *k2_w0_tmp;

  double *k3_q, *k3_qd, *k3_POS0, *k3_Qtn0, *k3_v0, *k3_w0;
  double *k4_q, *k4_qd, *k4_POS0, *k4_Qtn0, *k4_v0, *k4_w0;
  
  double *qdd_ori, *qd_ori, *q_ori; 
  double *vd0_ori, *wd0_ori, *v0_ori, *w0_ori, *POS0_ori, *Qtn0_ori, *A0_ori;
  
  double *qdd_tmp, *vd0_tmp, *wd0_tmp;
   
  // ---- memory allocation ----
  k1_q = matrix_get( m.LINKNUM, 1 );
  k2_q = matrix_get( m.LINKNUM, 1 );
  k3_q = matrix_get( m.LINKNUM, 1 );
  k4_q = matrix_get( m.LINKNUM, 1 );

  k1_q_tmp = matrix_get( m.LINKNUM, 1 );
  k2_q_tmp = matrix_get( m.LINKNUM, 1 );

  k1_qd = matrix_get( m.LINKNUM, 1 );
  k2_qd = matrix_get( m.LINKNUM, 1 );
  k3_qd = matrix_get( m.LINKNUM, 1 );
  k4_qd = matrix_get( m.LINKNUM, 1 );

  k1_qd_tmp = matrix_get( m.LINKNUM, 1 );
  k2_qd_tmp = matrix_get( m.LINKNUM, 1 );

  k1_POS0 = matrix_get( 3, 1 );
  k2_POS0 = matrix_get( 3, 1 );
  k3_POS0 = matrix_get( 3, 1 );
  k4_POS0 = matrix_get( 3, 1 );

  k1_POS0_tmp = matrix_get( 3, 1 );
  k2_POS0_tmp = matrix_get( 3, 1 );

  k1_Qtn0 = matrix_get( 4, 1 );
  k2_Qtn0 = matrix_get( 4, 1 );
  k3_Qtn0 = matrix_get( 4, 1 );
  k4_Qtn0 = matrix_get( 4, 1 );

  k1_Qtn0_tmp = matrix_get( 4, 1 );
  k2_Qtn0_tmp = matrix_get( 4, 1 );

  k1_v0 = matrix_get( 3, 1 );
  k2_v0 = matrix_get( 3, 1 );
  k3_v0 = matrix_get( 3, 1 );
  k4_v0 = matrix_get( 3, 1 );

  k1_v0_tmp = matrix_get( 3, 1 );
  k2_v0_tmp = matrix_get( 3, 1 );

  k1_w0 = matrix_get( 3, 1);
  k2_w0 = matrix_get( 3, 1);
  k3_w0 = matrix_get( 3, 1);
  k4_w0 = matrix_get( 3, 1 );

  k1_w0_tmp = matrix_get( 3, 1);
  k2_w0_tmp = matrix_get( 3, 1);
  
  qdd_ori = matrix_get( m.LINKNUM, 1 );
  qd_ori = matrix_get( m.LINKNUM, 1 );
  q_ori = matrix_get( m.LINKNUM, 1 );
  
  vd0_ori = matrix_get( 3, 1 );
  wd0_ori = matrix_get( 3, 1 );
  v0_ori = matrix_get( 3, 1 );
  w0_ori = matrix_get( 3, 1 );
  POS0_ori = matrix_get( 3, 1 );
  Qtn0_ori = matrix_get( 4, 1 );
  A0_ori = matrix_get( 3, 3 );
  
  qdd_tmp = matrix_get( m.LINKNUM, 1 );
  vd0_tmp = matrix_get( 3, 1 );
  wd0_tmp = matrix_get( 3, 1 );
  
  matrix_Z( m.LINKNUM, 1, qdd_tmp );
  matrix_Z( 3 , 1, vd0_tmp );
  matrix_Z( 3 , 1, wd0_tmp );
  
  // ---- copy original data ----
  matrix_cpy( 3, 1, m.vd0,  vd0_ori );
  matrix_cpy( 3, 1, m.wd0,  wd0_ori );
  matrix_cpy( 3, 1, m.v0,   v0_ori );
  matrix_cpy( 3, 1, m.w0,   w0_ori );
  matrix_cpy( 3, 1, m.POS0, POS0_ori );
  matrix_cpy( 4, 1, m.Qtn0, Qtn0_ori );
  matrix_cpy( 3, 3, m.A0, A0_ori );
  
  m.q[0] = 0; // invalid, to avoid computing problems
  m.qd[0] = 0;
  m.qdd[0] = 0;
  matrix_cpy( m.LINKNUM, 1, m.qdd, qdd_ori );
  matrix_cpy( m.LINKNUM, 1, m.qd, qd_ori );
  matrix_cpy( m.LINKNUM, 1, m.q, q_ori );
  
  // ================================================
  m.q[0] = 0; // invalid, to avoid computing problems
  m.qd[0] = 0;
  m.qdd[0] = 0;
  
  // 1st step
  w2dQtn( m.w0, m.Qtn0, m.dQtn0 );

  f_dyn_ffg( m, Gravity, m.vd0, m.wd0, m.qdd );
  
  // each joint 
  if( m.LINKNUM != 1){
    matrix_scale( m.LINKNUM, 1, step, m.qd, k1_q );
    matrix_scale( m.LINKNUM, 1, step, m.qdd, k1_qd );

    matrix_scale( m.LINKNUM, 1, 0.5, k1_q, k1_q_tmp );
    matrix_scale( m.LINKNUM, 1, 0.5, k1_qd, k1_qd_tmp );

    matrix_add( m.LINKNUM, 1, q_ori, k1_q_tmp, m.q );
    matrix_add( m.LINKNUM, 1, qd_ori, k1_qd_tmp, m.qd );
    
    // copy the acceleration
    matrix_cpy( m.LINKNUM, 1, m.qdd, qdd_tmp );
  }
  

  // base position and orientation (quartanion)
  matrix_scale( 3, 1, step, m.v0, k1_POS0 );
  matrix_scale( 4, 1, step, m.dQtn0, k1_Qtn0 );

  matrix_scale( 3, 1, 0.5, k1_POS0, k1_POS0_tmp );
  matrix_scale( 4, 1, 0.5, k1_Qtn0, k1_Qtn0_tmp );

  matrix_add( 3, 1, POS0_ori, k1_POS0_tmp, m.POS0 );
  matrix_add( 4, 1, Qtn0_ori, k1_Qtn0_tmp, m.Qtn0 );

  qtn2dc( m.Qtn0, m.A0 );
  matrix_trans( 3, 3, m.A0, m.A0 );

  // base linear and angular velocity
  matrix_scale( 3, 1, step, m.vd0, k1_v0 );
  matrix_scale( 3, 1, step, m.wd0, k1_w0 );

  matrix_scale( 3, 1, 0.5, k1_v0, k1_v0_tmp );
  matrix_scale( 3, 1, 0.5, k1_w0, k1_w0_tmp );

  matrix_add( 3, 1, v0_ori, k1_v0_tmp, m.v0 );
  matrix_add( 3, 1, w0_ori, k1_w0_tmp, m.w0 );

  // copy the acceleration0
  matrix_cpy( 3, 1, m.vd0, vd0_tmp);
  matrix_cpy( 3, 1, m.wd0, wd0_tmp);
  
  // ================================================
  m.q[0] = 0; // invalid, to avoid computing problems
  m.qd[0] = 0;
  m.qdd[0] = 0;

  // 2nd step
  w2dQtn( m.w0, m.Qtn0, m.dQtn0 );
  f_dyn_ffg( m, Gravity, m.vd0, m.wd0, m.qdd );
  
  // each joint 
  if( m.LINKNUM != 1){
    matrix_scale( m.LINKNUM, 1, step, m.qd, k2_q );
    matrix_scale( m.LINKNUM, 1, step, m.qdd, k2_qd );

    matrix_scale( m.LINKNUM, 1, 0.5, k2_q, k2_q_tmp );
    matrix_scale( m.LINKNUM, 1, 0.5, k2_qd, k2_qd_tmp );

    matrix_add( m.LINKNUM, 1, q_ori, k2_q, m.q );
    matrix_add( m.LINKNUM, 1, qd_ori, k2_qd, m.qd );
    
    // copy the acceleration
    matrix_add( m.LINKNUM, 1, qdd_tmp, m.qdd, qdd_tmp );

  }

  // base position and orientation (quartanion)
  matrix_scale( 3, 1, step, m.v0, k2_POS0 );
  matrix_scale( 4, 1, step, m.dQtn0, k2_Qtn0 );

  matrix_scale( 3, 1, 0.5, k2_POS0, k2_POS0_tmp );
  matrix_scale( 4, 1, 0.5, k2_Qtn0, k2_Qtn0_tmp );

  matrix_add( 3, 1, POS0_ori, k2_POS0_tmp, m.POS0 );
  matrix_add( 4, 1, Qtn0_ori, k2_Qtn0_tmp, m.Qtn0 );

  qtn2dc( m.Qtn0, m.A0 );
  matrix_trans( 3, 3, m.A0, m.A0 );

  // base linear and angular velocity
  matrix_scale( 3, 1, step, m.vd0, k2_v0 );
  matrix_scale( 3, 1, step, m.wd0, k2_w0 );

  matrix_scale( 3, 1, 0.5, k2_v0, k2_v0_tmp );
  matrix_scale( 3, 1, 0.5, k2_w0, k2_w0_tmp );

  matrix_add( 3, 1, v0_ori, k2_v0_tmp, m.v0 );
  matrix_add( 3, 1, w0_ori, k2_w0_tmp, m.w0 );
    
  // copy the acceleration0
  matrix_add( 3, 1, vd0_tmp, m.vd0, vd0_tmp);
  matrix_add( 3, 1, wd0_tmp, m.wd0, wd0_tmp);
    
  // ================================================
  m.q[0] = 0; // invalid, to avoid computing problems
  m.qd[0] = 0;
  m.qdd[0] = 0;

  // 3nd step
  w2dQtn( m.w0, m.Qtn0, m.dQtn0 );
  f_dyn_ffg( m, Gravity, m.vd0, m.wd0, m.qdd );
	
  // each joint 
  if( m.LINKNUM != 1){
    matrix_scale( m.LINKNUM, 1, step, m.qd, k3_q );
    matrix_scale( m.LINKNUM, 1, step, m.qdd, k3_qd );

    matrix_add( m.LINKNUM, 1, q_ori, k3_q, m.q );
    matrix_add( m.LINKNUM, 1, qd_ori, k3_qd, m.qd );
    
    // copy the acceleration
    matrix_add( m.LINKNUM, 1, qdd_tmp, m.qdd, qdd_tmp );

  }
  
  // base position and orientation (quartanion)
  matrix_scale( 3, 1, step, m.v0, k3_POS0 );
  matrix_scale( 4, 1, step, m.dQtn0, k3_Qtn0 );

  matrix_add( 3, 1, POS0_ori, k3_POS0, m.POS0 );
  matrix_add( 4, 1, Qtn0_ori, k3_Qtn0, m.Qtn0 );

  qtn2dc( m.Qtn0, m.A0 );
  matrix_trans( 3, 3, m.A0, m.A0 );

  // base linear and angular velocity
  matrix_scale( 3, 1, step, m.vd0, k3_v0 );
  matrix_scale( 3, 1, step, m.wd0, k3_w0 );

  matrix_add( 3, 1, v0_ori, k3_v0, m.v0 );
  matrix_add( 3, 1, w0_ori, k3_w0, m.w0 );
  
  // copy the acceleration0
  matrix_add( 3, 1, vd0_tmp, m.vd0, vd0_tmp );
  matrix_add( 3, 1, wd0_tmp, m.wd0, wd0_tmp );
 
  // ================================================
  m.q[0] = 0; // invalid, to avoid computing problems
  m.qd[0] = 0;
  m.qdd[0] = 0;

  // 4th step
  w2dQtn( m.w0, m.Qtn0, m.dQtn0 );
  f_dyn_ffg( m, Gravity, m.vd0, m.wd0, m.qdd );

  // each joint 
  if( m.LINKNUM != 1){
    matrix_scale( m.LINKNUM, 1, step, m.qd, k4_q );
    matrix_scale( m.LINKNUM, 1, step, m.qdd, k4_qd );
    
    // copy the acceleration
    matrix_add( m.LINKNUM, 1, qdd_tmp, m.qdd, qdd_tmp );
  }

  // base position and orientation (quartanion)
  matrix_scale( 3, 1, step, m.v0, k4_POS0 );
  matrix_scale( 4, 1, step, m.dQtn0, k4_Qtn0 );

  qtn2dc( m.Qtn0, m.A0 );
  matrix_trans( 3, 3, m.A0, m.A0 );

  // base linear and angular velocity
  matrix_scale( 3, 1, step, m.vd0, k4_v0 );
  matrix_scale( 3, 1, step, m.wd0, k4_w0 );

  // copy the acceleration0
  matrix_add( 3, 1, vd0_tmp, m.vd0, vd0_tmp);
  matrix_add( 3, 1, wd0_tmp, m.wd0, wd0_tmp);
 
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
  // ---- final solution ----
  if( m.LINKNUM != 1){
    // ==== joint angle ====
    matrix_scale( m.LINKNUM, 1, 2, k2_q, k2_q );
    matrix_scale( m.LINKNUM, 1, 2, k3_q, k3_q );
    
    matrix_add( m.LINKNUM, 1, k1_q, k2_q, k1_q );
    matrix_add( m.LINKNUM, 1, k1_q, k3_q, k1_q );
    matrix_add( m.LINKNUM, 1, k1_q, k4_q, k1_q );
    
    matrix_scale( m.LINKNUM, 1, recip, k1_q, k1_q );
    matrix_add( m.LINKNUM, 1, q_ori, k1_q, m.q ); // joint angle = q
    
 
    // ==== joint angular velocity ====
    matrix_scale( m.LINKNUM, 1, 2, k2_qd, k2_qd );
    matrix_scale( m.LINKNUM, 1, 2, k3_qd, k3_qd );
    
    matrix_add( m.LINKNUM, 1, k1_qd, k2_qd, k1_qd );
    matrix_add( m.LINKNUM, 1, k1_qd, k3_qd, k1_qd );
    matrix_add( m.LINKNUM, 1, k1_qd, k4_qd, k1_qd );
    
    matrix_scale( m.LINKNUM, 1, recip, k1_qd, k1_qd );
    matrix_add( m.LINKNUM, 1, qd_ori, k1_qd, m.qd ); // joint angular velocity = qd
    
    } 
     
    // ==== base position ====
    matrix_scale( 3, 1, 2, k2_POS0, k2_POS0 );
    matrix_scale( 3, 1, 2, k3_POS0, k3_POS0 );
    
    matrix_add( 3, 1, k1_POS0, k2_POS0, k1_POS0 );
    matrix_add( 3, 1, k1_POS0, k3_POS0, k1_POS0 );
    matrix_add( 3, 1, k1_POS0, k4_POS0, k1_POS0 );
    
    matrix_scale( 3, 1, recip, k1_POS0, k1_POS0 );
    matrix_add( 3, 1, POS0_ori, k1_POS0, m.POS0 ); // base position = POS0
    
    // ==== base orientation (quartanion) ====
    matrix_scale( 4, 1, 2, k2_Qtn0, k2_Qtn0 );
    matrix_scale( 4, 1, 2, k3_Qtn0, k3_Qtn0 );
    
    matrix_add( 4, 1, k1_Qtn0, k2_Qtn0, k1_Qtn0 );
    matrix_add( 4, 1, k1_Qtn0, k3_Qtn0, k1_Qtn0 );
    matrix_add( 4, 1, k1_Qtn0, k4_Qtn0, k1_Qtn0 );
    
    matrix_scale( 4, 1, recip, k1_Qtn0, k1_Qtn0 );
    matrix_add( 4, 1, Qtn0_ori, k1_Qtn0, m.Qtn0 ); // base orientation = Qtn0
    
    qtn2dc( m.Qtn0, m.A0 );
    matrix_trans( 3, 3, m.A0, m.A0 ); // base orientation = A0
    
    // ==== base linear velocity ====
    matrix_scale( 3, 1, 2, k2_v0, k2_v0 );
    matrix_scale( 3, 1, 2, k3_v0, k3_v0 );
    
    matrix_add( 3, 1, k1_v0, k2_v0, k1_v0 );
    matrix_add( 3, 1, k1_v0, k3_v0, k1_v0 );
    matrix_add( 3, 1, k1_v0, k4_v0, k1_v0 );
    
    matrix_scale( 3, 1, recip, k1_v0, k1_v0 );
    matrix_add( 3, 1, v0_ori, k1_v0, m.v0 ); // base linear velocity = v0 

    // ==== base angular velocity ====
    matrix_scale( 3, 1, 2, k2_w0, k2_w0 );
    matrix_scale( 3, 1, 2, k3_w0, k3_w0 );
    
    matrix_add( 3, 1, k1_w0, k2_w0, k1_w0 );
    matrix_add( 3, 1, k1_w0, k3_w0, k1_w0 );
    matrix_add( 3, 1, k1_w0, k4_w0, k1_w0 );
    
    matrix_scale( 3, 1, recip, k1_w0, k1_w0 );
    matrix_add( 3, 1, w0_ori, k1_w0, m.w0 ); // base angular velocity = w0
    
    // accelearation
    matrix_scale( 3, 1, recip4, vd0_tmp, m.vd0 );
    matrix_scale( 3, 1, recip4, wd0_tmp, m.wd0 );
    matrix_scale( m.LINKNUM, 1, recip4, qdd_tmp, m.qdd );
    
    
    // ---- memory release ----
    delete [] k1_q;
    delete [] k1_qd;
    delete [] k1_POS0;
    delete [] k1_Qtn0;
    delete [] k1_v0;
    delete [] k1_w0;
    delete [] k2_q;
    delete [] k2_qd;
    delete [] k2_POS0;
    delete [] k2_Qtn0;
    delete [] k2_v0;
    delete [] k2_w0;
    delete [] k1_q_tmp;
    delete [] k1_qd_tmp;
    delete [] k1_POS0_tmp;
    delete [] k1_Qtn0_tmp;
    delete [] k1_v0_tmp;
    delete [] k1_w0_tmp;
    delete [] k2_q_tmp;
    delete [] k2_qd_tmp;
    delete [] k2_POS0_tmp;
    delete [] k2_Qtn0_tmp;
    delete [] k2_v0_tmp;
    delete [] k2_w0_tmp;

    delete [] k3_q;
    delete [] k3_qd;
    delete [] k3_POS0;
    delete [] k3_Qtn0;
    delete [] k3_v0;
    delete [] k3_w0;
    delete [] k4_q;
    delete [] k4_qd;
    delete [] k4_POS0;
    delete [] k4_Qtn0;
    delete [] k4_v0;
    delete [] k4_w0;
  
    delete [] qdd_ori;
    delete [] qd_ori;
    delete [] q_ori; 
    delete [] vd0_ori;
    delete [] wd0_ori;
    delete [] v0_ori;
    delete [] w0_ori;
    delete [] POS0_ori;
    delete [] Qtn0_ori;
    delete [] A0_ori;
  
    delete [] qdd_tmp;
    delete [] vd0_tmp;
    delete [] wd0_tmp;
  }

// === EOF ===
