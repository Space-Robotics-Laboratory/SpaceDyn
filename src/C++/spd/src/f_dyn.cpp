//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
//
// Function : f_dyn( MODEL, Gravity, vd0, wd0, qdd )
//            Forward dynamics
//            calculate the acceleration of the base ( vd0, wd0 ) and the arm ( qdd ) 
//
// Output : vd0  3x1
//          wd0  3x1
//          qdd  nx1
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
void f_dyn(MODEL &m, double *Gravity, double *vd0, double *wd0, double *qdd){

  int i, j;
  Vector6 *v, *a;
  Matrix6 *IA;
  
  double *u;
  double *d;
  
  Vector6 *fe, *h, *p, *c;
 
  Matrix6 Xup_T;
  Vector6 S_T, h_T, F0, Fe, f0;
  double d_recip = 1;

  double *tmp, *tmp1, *tmp2, *tmp3, *tmp4, *tmp5, *tmp6;

  Matrix6 Eye6;
  Vector3 v0_tmp, w0_tmp, vd0_tmp;

 
  IA = new Matrix6[m.LINKNUM];
  fe = new Vector6[m.LINKNUM];
  h  = new Vector6[m.LINKNUM];
  p  = new Vector6[m.LINKNUM];
  c  = new Vector6[m.LINKNUM];
  u  = new double[m.LINKNUM];
  d  = new double[m.LINKNUM];

  v  = new Vector6[m.LINKNUM];
  a  = new Vector6[m.LINKNUM];

  for(i=0;i<m.LINKNUM;i++){
    matrix_Z( 6, 1, v[i]);
    matrix_Z( 6, 1, a[i]);
    
    matrix_Z( 6, 6, IA[i]);
    matrix_Z( 6, 1, fe[i]);
    matrix_Z( 6, 1, h[i]);
    matrix_Z( 6, 1, p[i]);
    matrix_Z( 6, 1, c[i]);

    u[i] = 0;
    d[i] = 0;
  }

  tmp  = new Vector6;
  tmp1 = new Matrix6;
  tmp2 = new Vector6;
  tmp3 = new double;
  tmp4 = new double;
  tmp5 = new Vector6;
  tmp6 = new double;

  matrix_Z( 6, 1, F0 );
  matrix_I( 6, Eye6 );

  for(i=0;i<3;i++){
    tmp[i] = m.v0[i];
    tmp[i+3] = m.w0[i];
  }

  matrix_trans( 6, 6, m.Xup[0], Xup_T );
  matrix_mult( 6, 6, 1, Xup_T, tmp, v[0] );

  matrix_cpy( 6, 6, m.Isp[0], IA[0] );

  crossF( v[0], tmp1 );
  matrix_mult( 6, 6, 1, m.Isp[0], v[0], tmp );
  matrix_mult( 6, 6, 1, tmp1, tmp, p[0] );

  matrix_Z( 6, 1, c[0]);
  for(i=0; i<m.LINKNUM;i++)
    matrix_Z( 6, 1, fe[i]);

  if( m.LINKNUM != 1 ){
    for(i=1;i<m.LINKNUM;i++){

      matrix_trans( 6, 6, m.Xup[i], Xup_T );
      matrix_mult( 6, 6, 1, Xup_T, v[m.BB[i]], tmp );
      matrix_scale( 6, 1, m.qd[i], m.S[i], tmp2 );
      matrix_add( 6, 1, tmp, tmp2, v[i] );
      crossM( v[i], tmp1 );

      matrix_mult( 6, 6, 1, tmp1, tmp2, c[i] );

      matrix_cpy( 6, 6, m.Isp[i], IA[i] );
      crossF( v[i], tmp1 );

      matrix_mult( 6, 6, 1, m.Isp[i], v[i], tmp );
      matrix_mult( 6, 6, 1, tmp1, tmp, p[i] );

      if(m.EE[i] != 0){
	for(j=0;j<3;j++){
	  Fe[j] = m.Fe[i][j];
	  Fe[j+3] = m.Te[i][j];
	}
	matrix_mult( 6, 6, 1, m.jXe[i], Fe, fe[i] );
      }
    }

    for(i=m.LINKNUM-1;i>0;i--){
      
      matrix_trans( 6, 6, m.Xup[i], Xup_T );
      
      matrix_mult( 6, 6, 1, IA[i], m.S[i], h[i] );
      matrix_trans( 6, 1, m.S[i], S_T );
      matrix_mult( 1, 6, 1, S_T, h[i], tmp3);
      d[i] = *tmp3;

      matrix_trans( 6, 1, h[i], h_T );
      matrix_mult( 1, 6, 1, h_T, c[i], tmp3 );
      matrix_mult( 1, 6, 1, S_T, p[i], tmp4 );
      u[i] = m.tau[i] - *tmp3 - *tmp4;
      
      d_recip = 1/(d[i]);

      matrix_mult( 6, 1, 6, h[i], h_T, tmp1 );
      matrix_scale( 6, 6, d_recip, tmp1, tmp1 );

      matrix_sub( 6, 6, IA[i], tmp1, tmp1 );
      matrix_mult( 6, 6, 6, tmp1, Xup_T, tmp1 );
      matrix_mult( 6, 6, 6, m.Xup[i], tmp1, tmp1 );

      matrix_add( 6, 6, IA[m.BB[i]], tmp1, IA[m.BB[i]] );

      matrix_scale( 6, 1, u[i], h[i], tmp );
      matrix_scale( 6, 1, d_recip, tmp, tmp );

      matrix_mult( 6, 6, 1, IA[i], c[i], tmp5 );
      matrix_add( 6, 1, tmp, tmp5, tmp );
      matrix_add( 6, 1, tmp, p[i], tmp );
      matrix_mult( 6, 6, 1, m.Xup[i], tmp, tmp );

      matrix_add( 6, 1, p[m.BB[i]], tmp, p[m.BB[i]] );


      matrix_mult( 6, 1, 6, h[i], S_T, tmp1 );
      matrix_scale( 6, 6, d_recip, tmp1, tmp1 );
      matrix_sub( 6, 6, Eye6, tmp1, tmp1 );
      matrix_mult( 6, 6, 1, tmp1, fe[i], tmp );
      matrix_mult( 6, 6, 1, m.Xup[i], tmp, tmp );
      matrix_add( 6, 1, fe[m.BB[i]], tmp, fe[m.BB[i]] );
    }
  }

  for(i=0;i<3;i++){
    F0[i] = m.F0[i];
    F0[i+3] = m.T0[i];
  }

  matrix_trans( 6, 6, m.Xup[0], Xup_T );
  matrix_mult( 6, 6, 1, Xup_T, F0, tmp );
  matrix_add( 6, 1, tmp, fe[0], f0 );
  
  matrix_sub( 6, 1, f0, p[0], a[0] );
  matrix_inv( 6, IA[0], tmp1 );
  matrix_mult( 6, 6, 1, tmp1, a[0], a[0] );
  matrix_mult( 6, 6, 1, Xup_T, Gravity, tmp );
  matrix_add( 6, 1, a[0], tmp, a[0] );

  if(m.LINKNUM != 1){
    for(i=1;i<m.LINKNUM;i++){
      
      matrix_trans( 6, 6, m.Xup[i], Xup_T );
      matrix_mult( 6, 6, 1, Xup_T, a[m.BB[i]], a[i] );
      // add
//      matrix_mult( 6, 6, 1, Xup_T, Gravity, tmp );
//      matrix_scale( 6, 1, tmp[i], m.JtoC[i], tmp );
//      matrix_add( 6, 1, a[i], tmp, a[i] );
      //
      d_recip = 1/(d[i]);
      matrix_trans( 6, 1, h[i], h_T);
      matrix_trans( 6, 1, m.S[i], S_T);
      
      matrix_mult( 1, 6, 1, h_T, a[i], tmp3 );
      matrix_mult( 1, 6, 1, S_T, fe[i], tmp4 );
      qdd[i] = ( u[i] - (*tmp3) + (*tmp4) ) * (d_recip);

      matrix_scale( 6, 1, m.qdd[i], m.S[i], tmp );
      matrix_add( 6, 1, a[i], c[i], a[i] );
      matrix_add( 6, 1, a[i], tmp, a[i] );
    
    }
  }

  matrix_mult(6, 6, 1, m.Xup[0], a[0], tmp );

  for(i=0;i<3;i++){
    v0_tmp[i] = v[0][i];
    w0_tmp[i] = v[0][i+3];
    vd0_tmp[i] = a[0][i];
  }

  vector_cross3( v0_tmp, w0_tmp, tmp5 );
  matrix_sub( 3, 1, vd0_tmp, tmp5, vd0_tmp);
  
  for(i=0;i<3;i++){
    a[0][i] = vd0_tmp[i];
  }


  matrix_mult( 6, 6, 1, m.Xup[0], a[0], a[0] );

  for(i=0;i<3;i++){
    vd0[i] = a[0][i];
    wd0[i] = a[0][i+3];
  } 

  delete [] v;
  delete [] a;

  delete [] tmp;
  delete [] tmp1;
  delete [] tmp2;
  delete tmp3;
  delete tmp4;
  delete [] tmp5;
  delete tmp6;
  
  delete [] IA;
  delete [] fe;
  delete [] h;
  delete [] p;
  delete [] c;
  
  delete [] u;
  delete [] d;
  
}

// === EOF ===
