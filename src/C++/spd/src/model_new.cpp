//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
//
// Function : model_new( string, MODEL )
//            Set a model file from new connectivity
//
//
// s.abiko [2007.5]
//
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/rot.h"
#include "../include/spn.h"
#include "../include/spd.h"

void model_new( MODEL &m_o, MODEL &m )
{
  int i;
  int LINKNUM_New;
  int *BB;

  newnumberofjoints(&n_max, &n_new, &B, Bcon, Cnew, mass, inertia, rcmvec, lvec, uvec, mass_new, inertia_new, rcmvec_new, lvec_new, uvec_new, theta);

  // to check the file w/o fortran 
  LINKNUM_New = 5;
  BB = new int[LINKNUM_New];
  for(i=1;i<LINKNUM_New; i++){
    BB[i] = i-1;
  }

  cout << "link_number = " << LINKNUM_New << endl;
  m.LINKNUM = LINKNUM_New;

  /* create arrays of MODEL class */
  m.constructor();

  m.BB = Bnew;
  m.J_Type = J_Type_new;
  m.EE = EE_new;
  m.construct_ee( E_NUM );
  
  m.Lflag = Lflag_new;
 
  for(j=0; j<3; j++){
      m.Qi[i][j] = Qi_new;
    }
  
    for(j=0;j<3;j++){
      m.CtoJ[i][j] = CtoJ;
    }
    for(j=0;j<3;j++){
      m.JtoC[i][j] = JtoC;
    }

    for(j=0; j<3; j++){
      m.CtoE[i][j] = CtoE;
      m.Qe[i][j] = Qe;
    }

    m.Link_M[i] = Mass_new;
    m.Link_I[i][j] = inertia _new;

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // convert all parameters to Spatial Notation
  calc_SPN( m );
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/


}

// === EOF ===
