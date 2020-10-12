//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
//
// Function : model_init( MODEL, Gravity )
//            Initialze State Values of Model and Envirounment
//
// s.abiko [2007.5]
// t.ikuta [2011.10]
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/spd.h"
#include "../include/rot.h"
#include "../include/spn.h"

void model_init( MODEL &m )
{  
  int i, j;
  
  for( i=0; i<m.LINKNUM;i++){
	matrix_Z( 6, 6, m.Xup[i] );
  	matrix_Z( 6, 6, m.Xup_I[i] );
  }

  m.Qtn0[0] = 1;
  m.Qtn0[1] = 0;
  m.Qtn0[2] = 0;
  m.Qtn0[3] = 0;

  m.dQtn0[0] = 0;
  m.dQtn0[1] = 0;
  m.dQtn0[2] = 0;
  m.dQtn0[3] = 0;

  m.A0[0] = 1.0;
  m.A0[1] = 0.0;
  m.A0[2] = 0.0;
  m.A0[3] = 0.0;
  m.A0[4] = 1.0;
  m.A0[5] = 0.0;
  m.A0[6] = 0.0;
  m.A0[7] = 0.0;
  m.A0[8] = 1.0;

  for( i=0; i<3; i++){
    m.POS0[i] = 0.0;
    
    m.v0[i] = 0;
    m.w0[i] = 0;

    m.vd0[i]= 0;
    m.wd0[i]= 0;

    m.F0[i] = 0;
    m.T0[i] = 0;
  }
  
//  if( m.LINKNUM != 1){
    for(i=0;i<m.LINKNUM;i++){
      m.q[i] = 0;
      m.qd[i] = 0;
      m.qdd[i] = 0;

      m.qm[i] = 0;
      m.qmd[i] = 0;
      m.qmdd[i] = 0;

      m.dq[i] = 0;
      m.dqd[i] = 0;
      m.dqdd[i] = 0;

      m.tau[i] = 0;
      m.tauM[i] = 0;

	for(j=0;j<3;j++){
	  m.Fe[i][j] = 0;
	  m.Te[i][j] = 0;
      }
    }
//  }
  
  for( i=0; i<(6+m.LINKNUM-1);i++)
    m.Force[i] = 0;

  cout << endl << " /_/_/_/ State Values are initialized, OK !! /_/_/_/ " << endl << endl;
  
}

// === EOF ===
