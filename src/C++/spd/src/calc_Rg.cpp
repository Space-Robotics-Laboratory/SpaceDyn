//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
//
// Function : calc_Rg( MODEL, Rg )
//            calculate the center of mass of the system
//
// Output : (3 x 1) Rg           
//
// s.abiko [2007.5]
//
// modified by s.abiko [2008.10] 
// to adjust the case of changing the structure
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/rot.h"
#include "../include/spn.h"
#include "../include/spd.h"

void calc_Rg( MODEL &m, double* Rg ){

	int i;
	double *CoM0, *MoM, *tmp2, *tmp3, *RR_i, *Xc_I, total_M;
	
	total_M = 0.0;
	for(i=0;i<m.LINKNUM;i++)
		total_M += m.link_M[i];
		
	
	CoM0 = matrix_get( 3, 1 );
	MoM  = matrix_get( 3, 1 );
	tmp2 = matrix_get( 3, 3 );
	tmp3 = matrix_get( 3, 3 );
	RR_i = matrix_get( 3, 1 );
	
	Xc_I = matrix_get( 6, 6 );
	
//	matrix_scale( 3, 1, m.link_M[0], m.POS0, MoM );
	matrix_add( 3, 1, m.POS0, m.JtoC[0], CoM0);
	matrix_scale( 3, 1, m.link_M[0], CoM0, MoM );

	calc_Xup_I( m );
	
	for( i=1;i<m.LINKNUM;i++){
		matrix_mult( 6, 6, 6, m.Xup_I[i], m.jXc[i], Xc_I );
		
		matrix_ext_sub( 6, 6, 1, 3, 1, 3, Xc_I, tmp2 ); // Xc_I(1:3, 1:3)
		matrix_trans( 3, 3, tmp2, tmp2 ); // Xc_I(1:3, 1:3)'
		
		matrix_ext_sub( 6, 6, 4, 6, 1, 3, Xc_I, tmp3 ); // Xc_I(4:6, 1:3)
		
		matrix_mult( 3, 3, 3, tmp3, tmp2, tmp2 ); //  Xc_I(4:6, 1:3)*Xc_I(1:3, 1:3)'
		i_tilde( tmp2, RR_i ); // RR_i = Centroid of each link
		
		matrix_scale( 3, 1, m.link_M[i], RR_i, RR_i );
		matrix_add( 3, 1, MoM, RR_i, MoM );		
	}
	
	total_M = 1/total_M;

	matrix_scale( 3, 1, total_M, MoM, Rg );
	
	delete [] CoM0;
 	delete [] MoM;
	delete [] tmp2;
	delete [] tmp3;
	delete [] RR_i;
	delete [] Xc_I;
	
}

// === EOF ===
