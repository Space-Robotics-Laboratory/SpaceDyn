//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
//
// Function : calc_hh( MODEL, HH )
//            calculate inertia matrix of free-floating robot
//
// Output: HH ( 6+LINKNUM-1 ) x ( 6+LINKNUM-1 ) matrix 
//
// s.abiko [2007.5]
//
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/rot.h"
#include "../include/spn.h"
#include "../include/spd.h"

void calc_hh( MODEL &m, double *HH )
{

  matrix_Z( (6+m.LINKNUM-1), (6+m.LINKNUM-1), HH );      
  
  int i, j, k, p;
  Matrix6 *IC;
  Matrix6 *f;
  double ***fm; // don't confuse
  
  double *tmp, *tmp1;
  double *tmp2;
  double *Xup_T;
  double *IC_i;
  double *S_T;

  IC = new Matrix6[m.LINKNUM];
  f  = new Matrix6[1];

  fm = new double**[m.LINKNUM];
  for( i=0; i<m.LINKNUM; i++ )
      fm[i] = new double*[(m.LINKNUM)];
  for( i=0; i<m.LINKNUM; i++ )
    for( j=0; j<m.LINKNUM; j++ ){     
      fm[i][j] = new Vector6;
      matrix_Z(6, 1, fm[i][j]);
    }

  // initialization
  IC_i = new Matrix6;
  matrix_Z(6,6, IC_i);

  Xup_T = new Matrix6;
  matrix_Z(6,6,Xup_T);

  tmp = new Vector6;
  tmp1 = new Vector6;
  tmp2 = new double;

  matrix_Z( 6, 1, tmp );
  matrix_Z( 6, 1, tmp1 );
 
  S_T = matrix_get( 1, 6 );
  matrix_Z( 1, 6, S_T );

  double *SS;
  SS = matrix_get(m.LINKNUM, m.LINKNUM); // children links
  calc_SS(m, SS);
  
  // calculate the Inertia Matrix in joint space
  for(i=0;i<m.LINKNUM;i++){
    matrix_cpy( 6, 6, m.Isp[i], IC[i] );
  }

  if( m.LINKNUM != 1 ){
    for( i=m.LINKNUM-1; i>0; i-- ){
      matrix_trans( 6, 6, m.Xup[i], Xup_T );
      matrix_mult( 6, 6, 6, IC[i], Xup_T, IC_i );
      matrix_mult( 6, 6, 6, m.Xup[i], IC_i, IC_i );
      matrix_add( 6, 6, IC[m.BB[i]], IC_i, IC[m.BB[i]] );
      
      matrix_mult( 6, 6, 1, IC[i], m.S[i], fm[i][i] );
    }
  }

  matrix_trans( 6, 6, m.Xup[0], Xup_T );
  matrix_mult( 6, 6, 6, IC[0], Xup_T, IC_i );
  matrix_mult( 6, 6, 6, m.Xup[0], IC_i, f[0] );

  // ---- multi body system ----
  if( m.LINKNUM != 1 ){
    j = 0;
    p = j;
    
    for( i=1; i<m.LINKNUM; i++ ){
      j = m.BB[i];
      while( j > 0 ){
	for( k=1; k<m.LINKNUM; k++ ){
	  if( SS[j*m.LINKNUM+k] == 1 ){//child link 	    
	    matrix_mult( 6, 6, 1, m.Xup[k], fm[i][k], tmp1 );
	    matrix_add( 6, 1, fm[i][j], tmp1, fm[i][j] );
	  }
	}
	p = j;
	j = m.BB[j];
      }

      if(m.BB[i] != 0)
	matrix_mult( 6, 6, 1, m.Xup[p],fm[i][p], fm[i][0] );
      else
	matrix_mult( 6, 6, 1, m.Xup[i],fm[i][i], fm[i][0] );   
      
      matrix_mult(6,6,1,m.Xup[0],fm[i][0],fm[i][0]);
    }

    // ---- allocation of the inertia matrix HH ----
    // ---- matrix in terms of the arms ----
    for( i=1; i<m.LINKNUM; i++ ){

      j = i;
      while( j > 0 ){		
	matrix_trans( 6, 1, m.S[i], S_T );
        matrix_mult( 1, 6, 1, S_T, fm[i][j], tmp2 );

	HH[ 6*( (m.LINKNUM-1)+6 )-1 + i*6 + (i-1)*(m.LINKNUM-1) + j ] = *tmp2;
	HH[ 6*( (m.LINKNUM-1)+6 )-1 + j*6 + (j-1)*(m.LINKNUM-1) + i ] = *tmp2;

	j = m.BB[j];
      }
      // coupling matrix between the base and the arm
      for(k=0;k<6;k++){
	HH[(6-1) + i + k*( 6+(m.LINKNUM-1) ) ] = fm[i][0][k];
	HH[ 6*( (m.LINKNUM-1)+6 ) + k + (i-1)*( (m.LINKNUM-1) + 6 ) ] = fm[i][0][k];
      }
      
    }
    // matrix in terms of the base
    for(i=0;i<6;i++){
      for(j=0;j<6;j++)
	HH[i*(6+m.LINKNUM-1)+j] = f[0][6*i+j];
    }     
  }
  // ---- single body system ----
  else{
    for(i=0;i<6;i++){
      for(j=0;j<6;j++)
	HH[ i*(6+m.LINKNUM-1)+j] = f[0][6*i+j];
    }
  }
  //matrix_print((6+m.LINKNUM-1),(6+m.LINKNUM-1),HH);

  delete [] IC;
  delete [] f;
  delete [] tmp;
  delete [] tmp1;
  delete tmp2;
  delete [] Xup_T;
  delete [] IC_i;
  delete [] S_T;
  delete [] SS;
  
  for( i=0; i<m.LINKNUM; i++ ){
    for( j=0; j<m.LINKNUM; j++ )
      	delete [] fm[i][j];
  }
  for(i=0; i<m.LINKNUM; i++ )
	delete [] fm[i];
  delete [] fm;

}


// === EOF ===
