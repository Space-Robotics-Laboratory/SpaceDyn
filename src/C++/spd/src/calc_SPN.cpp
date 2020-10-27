//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
//
// Function : calc_SPN( MODEL )
//            calculate Spatial Notation based on Model Parameters
//
// s.abiko [2007.5]
// takano bug fixed
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/rot.h"
#include "../include/spn.h"
#include "../include/spd.h"

#define _LINESKIP ifs >> str;

#include <QDebug>
/* /_/_/_/_/_/_/_/_/_/_/ calc_SPN  /_/_/_/_/_/_/_/_/_/_/_/_/
//   calculate Spatial Notation       */
void calc_SPN( MODEL &m )
{
  int i;
  
  double *tmp;
  double *zero3;
  double *LL;
  double *A_BB_i, *A_BB_e, *trans_tmp;

  // ---- initialization ----
  tmp = matrix_get(3,1);
  zero3 = matrix_get(3,1);
  LL = matrix_get(3,1);
  A_BB_i = matrix_get(6,6);
  A_BB_e = matrix_get(6,6);
  trans_tmp = matrix_get(6,6);

  // ---- base ----
  matrix_I(6,m.Xsp[0]);
  matrix_Z(3,1,zero3);

  spI( m.link_M[0], zero3, m.link_I[0], m.Isp[0] );
  
 
  // --- link ---- if multibody system ----
  if( m.LINKNUM != 1){
    for( i=1; i<m.LINKNUM; i++){
      // each link's parameters
      XrotBB( m.Qi[i], A_BB_i );
      if( m.BB[i] != 0 )
	matrix_add( 3, 1, m.JtoC[m.BB[i]], m.CtoJ[i], LL );
      // LL = m.JtoC[BB[i]] + m.CtoJ[i];  link length
      else // connected with base
	// ---- spatial transformation matrix i+1 -> i(base) ---
	matrix_add( 3, 1, m.JtoC[m.BB[i]], m.CtoJ[i], LL ); // modified by abiko to adjust varying structure
 	//matrix_cpy( 3, 1, m.CtoJ[i], LL);
      // ---- spatial transformation matrix i+1 -> i ---
      Xtrans(LL, trans_tmp);
      matrix_mult( 6, 6, 6, trans_tmp, A_BB_i, m.Xsp[i]);
      // m.Xsp[i] = Xtrans(LL)*rotE;
      
      // ---- spatial transformation matrix CoM i -> Joint i ----
      Xtrans(m.JtoC[i], m.jXc[i]);
      
      if(m.EE[i] != 0){
	XrotBB( m.Qe[m.EE[i]-1], A_BB_e );
	matrix_add( 3, 1, m.JtoC[i], m.CtoE[m.EE[i]-1], LL);
	//   LL = LP.ce(:,i) - LP.cc(:,i,i); link length for last link
	
	// ---- spatial transformation matrix EE -> Joint i ---
	Xtrans(LL, trans_tmp);
	matrix_mult( 6, 6, 6, trans_tmp, A_BB_e, m.jXe[i]);
	// m.jXe[i] = Xtrans(LL)*rotE;
      }
      else{
	matrix_I(6, m.jXe[i]);
      }
      
      // rigid body inertia MF6
      spI( m.link_M[i], m.JtoC[i], m.link_I[i], m.Isp[i] );  
      //matrix_print(6,6,m.Isp[i]);
      
     // joint type matrix
      if(m.J_type[i] != 1){  // revolute joint
        m.S[i][0] = 0.0;
        m.S[i][1] = 0.0;
        m.S[i][2] = 0.0;
        m.S[i][3] = 0.0;
        m.S[i][4] = 0.0;
        m.S[i][5] = 1.0;
      }else{ // prismatic joint
        m.S[i][0] = 0.0;
        m.S[i][1] = 0.0;
        m.S[i][2] = 1.0;
        m.S[i][3] = 0.0;
        m.S[i][4] = 0.0;
        m.S[i][5] = 0.0;
      }//matrix_print(6,1,m.S[i-1]);
    }
  }
/*      // to check the conversion
       cout << "\n" << endl;
       cout << "Isp - Regid Body Inertia :\n";
       for( i=0; i<m.LINKNUM; i++)
         matrix_print(6,6,m.Isp[i]);
       cout << "\n" << endl;
  
       cout << "Xsp - Spatial Transformation Matrix :\n";
       for( i=0; i<m.LINKNUM; i++)
         {      matrix_print(6,6,m.Xsp[i]);    
   	cout << "\n" << endl;
         }
  
       cout << "S - Joint Type Vector :\n";
       for( i=1; i<m.LINKNUM; i++)
         matrix_print(6,1,m.S[i]);
       cout << "\n" << endl;
  
       cout << "jXc - Spatial Transformation Matrix c_i to joint i:\n";
       for( i=1; i<m.LINKNUM; i++)
         matrix_print(6,6,m.jXc[i]);
       cout << "\n" << endl;
 
       cout << "jXe - Spatial Transformation Matrix EE to joint n:\n";
       for( i=1; i<m.LINKNUM; i++)
         matrix_print(6,6,m.jXe[i]);
       cout << "\n" << endl;
*/    
  delete [] tmp;
  delete [] zero3;
  delete [] LL;
  delete [] A_BB_i;
  delete [] A_BB_e;
  delete [] trans_tmp;
 
  cout << " /_/_/_/ model: Spatial Notation Convert finished /_/_/_/ " << endl;
}

