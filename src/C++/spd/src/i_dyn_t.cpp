//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
//
// Function : i_dyn( MODEL, Gravity, Force )
//            Inverse Dynamics
//            calculate the force on the system
//
// Output : Force ( 6 + (LINKNUM-1) ) x 1
//
// s.abiko [2007.5]
// tknrk
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/rot.h"
#include "../include/spn.h"
#include "../include/spd.h"


// ----  Inverse Dynamics ----
void i_dyn_t(MODEL &m, double *Gravity, double *Force){

    int i, j;
    Vector6 *f, *v, *a;

    Matrix6 Xup_T;
    Vector6 S_T, Fe, a0;

    double *tmp, *tmp1, *tmp2, *tmp3, *tmp4, *tmp5, *tmp6;

    f  = new Vector6[m.LINKNUM];
    v  = new Vector6[m.LINKNUM];
    a  = new Vector6[m.LINKNUM];

    tmp  = new Vector6;
    tmp1 = new Vector3;
    tmp2 = new Matrix6;
    tmp3 = new Matrix6;
    tmp4 = new Vector6;
    tmp5 = new Vector6;
    tmp6 = new double;


    matrix_Z( 6, 1, v[0] );
    matrix_Z( 6, 1, a[0] );
    matrix_Z( 6, 1, f[0] );

    if( m.LINKNUM != 1 ){
      for( i=1;i<m.LINKNUM;i++){

        matrix_trans( 6, 6, m.Xup[i], Xup_T );

        matrix_mult( 6, 6, 1, Xup_T, v[m.BB[i]], tmp );
        matrix_scale( 6, 1, m.qd[i], m.S[i], tmp4 );
        matrix_add( 6, 1, tmp, tmp4, v[i] );

        matrix_mult( 6, 6, 1, Xup_T, a[m.BB[i]], tmp );
        crossM( v[i], tmp3 );
        matrix_mult( 6, 6, 1, tmp3, tmp4, tmp4 );

        matrix_add( 6, 1, tmp, tmp4, a[i] );

        matrix_mult( 6, 6, 1, m.Isp[i], a[i], tmp );
        crossF( v[i], tmp3 );
        matrix_mult( 6, 6, 6, tmp3, m.Isp[i], tmp2 );
        matrix_mult( 6, 6, 1, tmp2, v[i], tmp4 );

        matrix_add( 6, 1, tmp, tmp4, f[i] );
      }

      for( i=m.LINKNUM-1; i>0; i-- ){

        matrix_trans( 6, 1, m.S[i], S_T );
        //      matrix_mult( 1, 6, 1, S_T, f[i], Force[i+5] );
        // double* ‚É•ÏŠ·‚Å‚«‚È‚¢•ª‚©‚ç‚È‚¢
        matrix_mult( 1, 6, 1, S_T, f[i], tmp6 );
        Force[i+5] = *tmp6;


        matrix_mult( 6, 6, 1, m.Xup[i], f[i], tmp );
        matrix_add( 6, 1, f[m.BB[i]], tmp, f[m.BB[i]] );
      }
    }

//    matrix_mult( 6, 6, 1, m.Xup[0], f[0], f[0]);

//    for(i=0;i<3;i++){

//      Force[i+3] = f[0][i+3];
//    }

    delete [] f;
    delete [] v;
    delete [] a;

    delete [] tmp;
    delete [] tmp1;
    delete [] tmp2;
    delete [] tmp3;
    delete [] tmp4;
    delete [] tmp5;
    delete tmp6;

}

// === EOF ===

