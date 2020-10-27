//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
//
// Function : i_dyn_CL( MODEL, Gravity, Force )
//            Inverse Dynamics
//            calculate the force on the system for closed loops
//
//
// Output : Force ( 6 + (LINKNUM-1)-(closed loop constraints) ) x 1
//
// s.abiko [2007.8]
//
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
#include <iostream>
#include <cmath>
using namespace std;

#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/rot.h"
#include "../include/spn.h"
#include "../include/spd.h"

void i_dyn_CL( MODEL &m, double *Gravity, double *Force_g ){
  
  //_/_/_/_/_/_/_/ calculate P, sensitivity matrix /_/_/_/_/_/_/_/_/_/_/_//
  int Nj, Nf, L_num;
  int i, j, rank;  

  double norm_dPOS;
  
  double **Je_Li, *Je_a, *Je_b, *JC;
  double *tmp, *identity3;
  double *s, *U, *V;
  double *dPOS_e, *J_de, *til_dPOS;
  
  L_num = 0;
  Nj = m.LINKNUM-1;

  Je_a = matrix_get( 6, Nj );
  Je_b = matrix_get( 6, Nj );
  
  dPOS_e = matrix_get( 3, 1 );
  til_dPOS = matrix_get( 3, 3 );
  J_de = matrix_get( 6, 6 );
  matrix_Z( 6, 6, J_de );

  identity3 = matrix_get(3,3);
  matrix_I( 3, identity3 );
  

  matrix_cpy_sub( 6, 6, 1, 3, 1, 3, identity3, J_de );
  matrix_cpy_sub( 6, 6, 4, 6, 4, 6, identity3, J_de );

  // count the number of closed loops
  for( i=1; i<m.LINKNUM; i++){
    if( m.Lflag[i] != 0){
      L_num++;
    }
  }
  // initialize the Jacobian Jc and Je_Li
  JC = matrix_get( 6*L_num, Nj );
  Je_Li = new double*[L_num];
  for( i=0 ; i<Nj; i++ )
 //   Je_Li[i] = new double[ 6*Nj ];
    Je_Li[i] = matrix_get(6, Nj);
 
  L_num = 0;
  // i-th closed loop constraint Jacobian Je_Li // case for connected to the link
  for( i=1; i<m.LINKNUM; i++ ){
    cout << m.EE[i] << ' ' << m.Lflag[i] << endl ;
    if(( m.Lflag[i] != 0 ) && ( m.Lflag[i] != -1 )){
    
      calc_Je( m, m.EE[i], Je_b );
      calc_Je( m, m.EE[m.Lflag[i]], Je_a );
      
      
      f_kin_e(m, m.EE[i]); // position and orientation of EE -> m.POS_e, m.ORI_e
      f_kin_e(m, m.EE[m.Lflag[i]]);
      
      matrix_sub( 3, 1, m.POS_e[m.EE[m.Lflag[i]]], m.POS_e[m.EE[i]], dPOS_e );
      norm_dPOS = sqrt( dPOS_e[0]*dPOS_e[0] + dPOS_e[1]*dPOS_e[1] + dPOS_e[2]*dPOS_e[2] );
      if( norm_dPOS >= 1e-12 ){ // gap between the EEs, but the distance to be fixed
      	tilde( 3, dPOS_e, til_dPOS );
	matrix_cpy_sub( 6, 6, 1, 3, 4, 6, til_dPOS, J_de );
	matrix_mult( 6, 6, Nj, J_de, Je_b, Je_b );
	matrix_sub( 6, Nj, Je_a, Je_b, Je_Li[L_num] );
      }
      else{//no gap between the EEs
	matrix_sub( 6, Nj, Je_a, Je_b, Je_Li[L_num] );
      }
    
       L_num++;
    }
    else if(( m.Lflag[i] == 0 ) && ( m.Lflag[i] != -1 )){// case for connected to the base

      calc_Je( m, m.EE[i], Je_b );
      
      f_kin_e(m, m.EE[i]); // position and orientation of EE -> m.POS_e, m.ORI_e
      
      matrix_sub( 3, 1, m.POS0, m.POS_e[m.EE[i]], dPOS_e );
      norm_dPOS = sqrt( dPOS_e[0]*dPOS_e[0] + dPOS_e[1]*dPOS_e[1] + dPOS_e[2]*dPOS_e[2] );
      if( norm_dPOS >= 1e-12 ){ // gap between the EEs, but the distance to be fixed
      	tilde( 3, dPOS_e, til_dPOS );
	matrix_cpy_sub( 6, 6, 1, 3, 4, 6, til_dPOS, J_de );
	matrix_mult( 6, 6, Nj, J_de, Je_b, Je_b );
	matrix_scale( 6, Nj, -1, Je_b, Je_Li[L_num] );
      }
      else{//no gap between the EEs
	matrix_scale( 6, Nj, -1, Je_b, Je_Li[L_num] );
      }
    
       L_num++;
    } 
  }
  
  // augumented closed loop constraints Jacobian Jc
  for( i=0; i<L_num; i++ ){
    matrix_cpy_sub( 6*L_num, Nj, (i+1)*6-5, (i+1)*6, 1, Nj, Je_Li[i], JC );
  }

  cout << "\n" << endl;
  cout << "Nj = " << Nj << endl;
  cout << "Loop Number = " << L_num << endl; 
  cout << "\n" << endl;
  
  matrix_print(6*L_num, Nj, JC );

  tmp = matrix_get( Nj, 6*L_num );
  // check the rank of Jc 
  if( 6*L_num > Nj ){
    s = matrix_get( 6*L_num, 1 );
    U = matrix_get( 6*L_num, 6*L_num );
    V = matrix_get( Nj, Nj );

    matrix_Z( 6*L_num, 1, s );
    matrix_Z( 6*L_num, 6*L_num, U );
    matrix_Z( Nj, Nj, V );
  }
  else{
    s = matrix_get( Nj, 1 );
    U = matrix_get( Nj, Nj );
    V = matrix_get( 6*L_num, 6*L_num );

    matrix_Z( Nj, 1, s );
    matrix_Z( Nj, Nj, U );
    matrix_Z( 6*L_num, 6*L_num, V );
  }

  rank = 0;

  if( 6*L_num > Nj ){
    matrix_svd( 6*L_num, Nj, JC,  U, s, V );
    matrix_print( 6*L_num, 1, s );
    for( i=0; i<6*L_num; i++){
      if( s[i] >= 1e-12)
	rank++;
    }
  }
  else{
    matrix_trans( 6*L_num, Nj, JC, tmp );
    matrix_svd( Nj, 6*L_num, tmp, U, s, V );
    matrix_print( Nj, 1, s );
    for( i=0; i<Nj; i++){
      if( s[i] >= 1e-12)
	rank++;
    }
  }
  cout << "rank = " << rank << endl;
  
  // convert to reduced row echelon form
  rref( 6*L_num, Nj, JC, JC );
  cout << "Jc = " << endl;
  matrix_print( 6*L_num, Nj, JC );
  

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
  double *H;
  double *Jcm, *Js, *Jg, *iJs;
  double *col, *row;
  int ck = 0, ck_Nf = 0;
  int num_I;
  int q_S[rank], q_G[Nj-rank];
  
  for(i=0; i<rank; i++){
  	q_S[i] = 0;
  }
  for(i=0; i<Nj-rank; i++){
  	q_G[i] = 0;
  }

  Nf = Nj - rank;
  H = matrix_get( rank, Nf );
  
  Jcm = matrix_get( rank, Nj );
  matrix_Z( rank, Nj, Jcm );
  matrix_ext_sub( 6*L_num, Nj, 1, rank, 1, Nj, JC, Jcm );
  cout << "Jcm = " << endl;
  matrix_print( rank, Nj, Jcm  ); // extract the independent constraints from Jc
  
  Js = matrix_get( rank, rank );
  Jg = matrix_get( rank, Nf );
  col = matrix_get( rank, 1 );
  row = matrix_get( 1, Nf );
  iJs = matrix_get( rank, rank );

  matrix_Z( rank, rank, Js );
  matrix_Z( rank, rank, iJs ); // for inverse of Js
  matrix_Z( rank, Nf, Jg );
  matrix_Z( rank, 1, col );

  for( i=0; i<rank; i++){
    num_I = 0;
    for( j=0; j<Nj; j++ ){
    	if( ( Jcm[i*Nj + j] == 1.0 ) && ( num_I != 1 ) ){
		q_S[ck] = j+1; // index for dependent joints
		matrix_ext_col(rank, Nj, q_S[ck], Jcm, col);
		matrix_cpy_col(rank, rank, ck+1, col, Js);
		num_I = 1;
		ck++;
	}
    }
  }
  cout << "Js = " << endl;
  matrix_print(rank, rank, Js);
  
  for( i=0; i<Nj; i++ ){
	num_I = 0;
  	for( j=0; j<rank; j++){
		if( q_S[j] == i+1 )
			num_I = 1;
	}
	if( num_I != 1 ){
		q_G[ck_Nf] = i+1;
		matrix_ext_col(rank, Nj, q_G[ck_Nf], Jcm, col);
		matrix_cpy_col(rank, Nf, ck_Nf+1, col, Jg);
		ck_Nf++;
		}
  }  
  cout << "Jg = " << endl;
  matrix_print(rank, Nf, Jg);
  
  matrix_inv( rank, Js, iJs );	
  matrix_mult( rank, rank, Nf, iJs, Jg, H );
  matrix_scale( rank, Nf, -1.0, H, H );
 
  cout << "H = " << endl;
  matrix_print( rank, Nf, H );

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
  double *Force;
  double *P_T, *P;
  
  Force = matrix_get( 6+Nj, 1 );
  matrix_Z( 6+Nj, 1, Force );
  
  P = matrix_get( 6+Nj, 6+Nf );
  matrix_Z( 6+Nj, 6+Nf, P );

  P_T = matrix_get( 6+Nf, 6+Nj );
  matrix_Z( 6+Nf, 6+Nj, P_T );


  matrix_cpy_sub( 6+Nj, 6+Nf, 1, 3, 1, 3, identity3, P );
  matrix_cpy_sub( 6+Nj, 6+Nf, 4, 6, 4, 6, identity3, P );
  
  for( i=0; i<rank; i++){
  	matrix_ext_row( rank, Nf, i+1, H, row );
  	matrix_cpy_sub( 6+Nj, 6+Nf, 6+q_S[i], 6+q_S[i], 7, 6+Nf, row, P ); 
  }
  for( i=0; i<Nf; i++){
  	matrix_Z( 1, Nf, row );
	row[i] = 1;
  	matrix_cpy_sub( 6+Nj, 6+Nf, 6+q_G[i], 6+q_G[i], 7, 6+Nf, row, P ); 
  }
  
  cout << " P = " << endl;
  matrix_print( 6+Nj, 6+Nf, P );
  matrix_trans( 6+Nj, 6+Nf, P, P_T );

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
  // calculate force for the corresponding open chain structure 
  i_dyn( m, Gravity, Force );
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
  // calculate force for the closed chain structure  
  matrix_mult( 6+Nf, 6+Nj, 1, P_T, Force, Force_g );
  cout << " force for closed loop = " << endl;
  matrix_print( 6+Nf, 1, Force_g );
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
  

  // ---- clear the memory ---- I need to learn precisely since some parts do not work properly
  delete [] H;
  delete [] Jcm;
  delete [] Js;
  delete [] Jg;
  delete [] iJs;
  delete [] col;
  delete [] row;
  
  delete [] Force;
  delete [] P_T;
  delete [] P;
  
  delete [] JC;
  delete [] Je_a;
  delete [] Je_b;
  delete [] Je_Li;

  delete [] tmp;
  delete [] identity3;

  delete [] dPOS_e;
  delete [] til_dPOS;
  delete [] J_de;
  delete [] s;
  delete [] U;
  delete [] V;
}

// === EOF ===
