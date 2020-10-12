//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
//
// Function : calc_CLC2( MODEL &m, double *W, double *S )
//            calculate the matrix for closed loop constraints (CLC).
//            
//            W = No x Nf, W = Na x Nf
//            
// s.abiko [2008.3]
//
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
#include <iostream>
#include <cmath>
using namespace std;

#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/rot.h"
#include "../include/spn.h"
#include "../include/spd.h"

//void calc_CLC( MODEL &m, double *W,  double *S ){
void calc_CLC2( MODEL &m ){
  
  int Nj, Nf,  L_num;
//  int No, Na;
  int i, j, rank;  

  double *Jc, **Je_Li, *Je_a, *Je_b;
  double *tmp, *tmp2;
  double *s, *U, *V;
  
  L_num = 0;
  Nj = m.LINKNUM-1;

  Je_a = matrix_get( 6, Nj );
  Je_b = matrix_get( 6, Nj );

  matrix_Z( 6, Nj, Je_a );
  matrix_Z( 6, Nj, Je_b );
  
  calc_SP( m );
 
  // count the number of closed loops
  for( i=1; i<m.LINKNUM; i++){
    if( m.Lflag[i] != -1){
      L_num++;
    }
  }
  
  // initialize the Jacobian Jc and Je_Li
  Jc = matrix_get( 6*L_num, Nj );
  matrix_Z( 6*L_num, Nj, Jc );
  Je_Li = new double*[L_num];
  for( i=0 ; i<L_num; i++ ){
    Je_Li[i] = matrix_get(6, Nj);
    matrix_Z( 6, Nj, Je_Li[i] );
  }

  L_num = 0;
 // i-th closed loop constraint Jacobian Je_Li
  for( i=1; i<m.LINKNUM; i++ ){
    if( m.Lflag[i] != -1 ){
       cout << m.EE[i] << ' ' << m.Lflag[i] << ' ' << m.EE[m.Lflag[i]] << endl;
      
       cout << " i = " << i << endl;
        
       calc_Je( m, m.EE[i], Je_b );
//       matrix_print( 6, m.LINKNUM-1, Je_b );
       calc_Je( m, m.EE[m.Lflag[i]], Je_a );
//       matrix_print( 6, m.LINKNUM-1, Je_a );
       
       matrix_sub( 6, Nj, Je_a, Je_b, Je_Li[L_num] );
       L_num++;
       }
  }
 
  // augumented closed loop constraints Jacobian Jc
  for( i=0; i<L_num; i++ ){
    matrix_cpy_sub( 6*L_num, Nj, (i+1)*6-5, (i+1)*6, 1, Nj, Je_Li[i], Jc );
  }

  cout << Nj << endl;
  cout << L_num << endl; 

  matrix_print(6*L_num, Nj, Jc);

  tmp = matrix_get( Nj, 6*L_num );
  // check the rank of Jc 
  if( 6*L_num > Nj ){
    s = matrix_get( 6*L_num, 1 );
    U = matrix_get( 6*L_num, 6*L_num );
    V = matrix_get( Nj, Nj );
  }
  else{
    s = matrix_get( Nj, 1 );
    U = matrix_get( Nj, Nj );
    V = matrix_get( 6*L_num, 6*L_num );
  }

  rank = 0;

  if( 6*L_num > Nj ){
    matrix_svd( 6*L_num, Nj, Jc, U, s, V );
    // matrix_print( 6*L_num, 1, s );
    for( i=0; i<Nj; i++ ){
      if( s[i] >= 1e-12 )
	rank++;
    }
  }
  else{
    matrix_trans( 6*L_num, Nj, Jc, tmp );
    matrix_svd( Nj, 6*L_num, tmp, U, s, V );
    // matrix_print( 6*L_num, 1, s );
    for( i=0; i<6*L_num; i++ ){
      if( s[i] >= 1e-12 )
	rank++;
    }
  }
  cout << "rank = " << rank << endl;
  
  // convert to reduced row echelon form
  rref( 6*L_num, Nj, Jc, Jc );
  cout << "Jc = " << endl;
  matrix_print( 6*L_num, Nj, Jc );

  double *H, *JH;// *W, *S,
  double *Js, *Jg, *iJs;
  double *col;
  int ck = 0, ck_Nf = 0;
  int num_z, num_I;
  int q_S[rank], q_G[Nj-rank];//, q_O[Nj-L_num];

/*  for(i=1; i<m.LINKNUM; i++){
    if( m.Lflag[i] == -1)
      q_O[i-1] = i;  // there is bug here
      cout << "q_O = " << q_O[i-1] << endl; 
  }
*/
  Nf = Nj - rank; // degrees-of-freedom (mobility)
  H = matrix_get( rank, Nf );
  tmp2 = matrix_get( 1, Nf );
  
  Js = matrix_get( rank, rank );
  Jg = matrix_get( rank, Nf );
  col = matrix_get( 1, 6*L_num );

  iJs = matrix_get( rank, rank );

  matrix_Z( rank, rank, Js );
  matrix_Z( rank, rank, iJs ); // for inverse of Js
  matrix_Z( rank, Nf, Jg );
  matrix_Z( 1, 6*L_num, col );
  matrix_Z( 1, Nf, tmp2 );

  for( i=0; i<Nj; i++ ){
    
    num_z = 0;
    num_I = 0;
    
    matrix_ext_col( 6*L_num, Nj, (i+1), Jc, col );
    
    for( j=0; j<6*L_num; j++){
      if( col[j] == 0 ) num_z++;
      if( col[j] == 1 ) num_I++;
    }
    
    if( (num_z == (6*L_num-1)) && (num_I == 1) ){
      q_S[ck] = i+1; // index of joints --- dependent joints (passive)
      // cout << "q_S[" << ck << "] = " << i+1 << endl;

      ck++;
      matrix_cpy_col( rank, rank, ck, col, Js );
    }
    else{
      q_G[ck_Nf] = i+1; // index of joints --- independent joints (active)
      cout << "q_G[" << ck_Nf << "] = " << i+1 << endl;

      ck_Nf++;
      matrix_cpy_col( rank, Nf, ck_Nf, col, Jg );

    }
  }

  cout << "Js = " << endl;
  matrix_print( rank, rank, Js );
  
  cout << "Jg = " << endl;
  matrix_print( rank, Nf, Jg );

  matrix_inv( rank, Js, iJs );
  matrix_mult( rank, rank, Nf, iJs, Jg, H );
  matrix_scale( rank, Nf, -1.0, H, H );
 
  cout << "H = " << endl;
  matrix_print( rank, Nf, H );

  JH = matrix_get( Nj, Nf );
  matrix_Z( Nj, Nf, JH );
  
  // reallocate each row for each corresponding joint - JH
  for(i=1;i<m.LINKNUM;i++){
  	for( j=0; j<rank; j++){
      		if( i == q_S[j] ){
			matrix_ext_row( rank, Nf, (j+1), H, tmp2 );
			matrix_cpy_row( Nj, Nf, i, tmp2, JH );
      		}
    	}
	for( j=0; j<Nf; j++){
      		if( i == q_G[j]){
			matrix_Z( 1, Nf, tmp2 );
			tmp2[j] = 1;
			//matrix_print( 1, Nf, tmp2 );
			matrix_cpy_row( Nj, Nf, i, tmp2, JH );
      		}
    	}	  
  }
  cout << " JH = " << endl;
  matrix_print( Nj, Nf, JH );

  /*
  // size(W) = No x Nf, size(S) = Na x Nf
  No = Nj - L_num;
  Na = Nj;

  W = matrix_get( No, Nf );
  S = matrix_get( Na, Nf );

  matrix_Z( No, Nf, W );
  matrix_Z( Na, Nf, S );

  // calculate W
  for( i=0; i<No; i++){
    for( j=0; j<rank; j++){
      if( q_O[i] == q_S[j] ){
	matrix_ext_row( rank, Nf, (j+1), H, tmp2 );
	matrix_cpy_row( No, Nf, (i+1), tmp2, W );
      }
    }
    for( j=0; j<Nj-rank; j++){
      if( q_O[i] == q_G[j]){
	matrix_Z( 1, Nf, tmp2 );
	tmp[j] = 1;
	matrix_cpy_row( No, Nf, (i+1), tmp2, W );
      }
    }
  }
  // calculate S
  for( i=0; i<Na; i++){
    for( j=0; j<rank; j++){
      if( q_O[i] == q_S[j] ){
	matrix_ext_row( rank, Nf, (j+1), H, tmp2 );
	matrix_cpy_row( No, Nf, (i+1), tmp2, S );
      }
    }
    for( j=0; j<Nj-rank; j++){
      if( q_O[i] == q_G[j]){
	matrix_Z( 1, Nf, tmp2 );
	tmp[j] = 1;
	matrix_cpy_row( No, Nf, (i+1), tmp2, S );
      }
    }
  }
 
  cout << " W = " << endl;
  matrix_print( No, Nf, W );
  cout << " S = " << endl;
  matrix_print( Na, Nf, S );
*/
  // ---- clear the memory ----
  delete [] Jc;
  delete [] Je_a;
  delete [] Je_b;
  delete [] tmp;
  delete [] tmp2;
  delete [] s;
  delete [] U;
  delete [] V;
  for( i=0; i<L_num; i++ ){
    delete [] Je_Li[i];
  }
  delete [] Je_Li;
  
  delete [] H;
  delete [] JH;
  //delete [] W;
  //delete [] S;
  delete [] Js;
  delete [] Jg;
  delete [] iJs;
  delete [] col;

}
// === EOF ===
