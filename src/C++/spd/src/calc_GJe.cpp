//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//
// Function : calc_GJe( MODEL, e_num, ans )
//            calculate generalized jacobian matrix for selected end-effector
//
//            ans = (6 x LINKNUM-1) matrix
//
// s.abiko [2007.8]
//
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/rot.h"
#include "../include/spn.h"
#include "../include/spd.h"

void calc_GJe(MODEL &m, int e_num, double *ans){

	double *HH;
	double *Hb, *Hbm, *Jb, *Jm;
	
	HH = matrix_get( (6+m.LINKNUM-1), (6+m.LINKNUM-1) );
	Hb = matrix_get( 6, 6 );
	Hbm = matrix_get( 6, m.LINKNUM-1 );
	
	Jb = matrix_get( 6, 6 );
	Jm = matrix_get( 6, m.LINKNUM-1 );
	
	calc_Jb( m, e_num, Jb );	
	calc_Je( m, e_num, Jm );

	calc_hh( m, HH );
	matrix_ext_sub( (6+m.LINKNUM-1), (6+m.LINKNUM-1), 1, 6, 1, 6, HH, Hb );
	matrix_ext_sub( (6+m.LINKNUM-1), (6+m.LINKNUM-1), 1, 6, 7, 6+m.LINKNUM-1, HH, Hbm );
		
	matrix_inv( 6, Hb, Hb );
	matrix_mult( 6, 6, m.LINKNUM-1, Hb, Hbm, Hbm );
	matrix_mult( 6, 6, m.LINKNUM-1, Jb, Hbm, Hbm );
	
	matrix_sub( 6, m.LINKNUM-1, Jm, Hbm, ans );
	
	delete [] HH;
	delete [] Hb;
	delete [] Hbm;
	delete [] Jb;
	delete [] Jm;
	
}

// === EOF ===
