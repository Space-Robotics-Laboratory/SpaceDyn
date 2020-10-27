//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//
// Function : calc_Jb( MODEL, e_num, ans )
//            calculate jacobian matrix related to the base for selected end-effector
//
//            ans = (6 x 6) matrix
//
// s.abiko [2007.8]
//
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/rot.h"
#include "../include/spn.h"
#include "../include/spd.h"


void calc_Jb( MODEL &m, int e_num, double *ans){

	double *tmp, *t_tmp;
	
	tmp = matrix_get( 3, 1 );
	t_tmp = matrix_get( 3, 3 );
	
	matrix_Z( 3, 1, tmp );
	matrix_Z( 3, 3, t_tmp );
	
	matrix_I( 6, ans );
	
	f_kin_e(  m, e_num );
	matrix_sub( 3, 1, m.POS_e[e_num], m.POS0, tmp );
	
	tilde( 3, tmp, t_tmp );
	matrix_scale( 3, 3, -1, t_tmp, t_tmp );
	
	matrix_cpy_sub( 6, 6, 1, 3, 4, 6, t_tmp, ans );
	
	delete [] tmp;
	delete [] t_tmp;
	
}

// === EOF ===
