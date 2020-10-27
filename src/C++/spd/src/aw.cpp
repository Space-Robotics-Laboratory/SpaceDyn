//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
//
// Function : aw( w0, d_time, E0 )
//            calculate a 3x3 transformation representing 
//            a rotation about the vector w0.
//
// Output : E0 3 x 3
//
// s.abiko [2007.8]
//
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
#include <cmath>

#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/rot.h"
#include "../include/spn.h"
#include "../include/spd.h"

void aw(double *w0, double d_time, double *E0){

	int i;
	
	double th, norm_w0;
	double *w;
	
	th = 0.0;
	w = matrix_get(3,1);

	matrix_I(3, E0);
	norm_w0 = sqrt( w0[0]*w0[0] + w0[1]*w0[1] + w0[2]*w0[2] );
	
	if( norm_w0 != 0 ){
		th = norm_w0 * d_time;
		for(i=0; i<3; i++)
			w[i] = w0[i] / norm_w0;
		
		E0[0] = cos(th) + w[0]*w[0]*(1-cos(th));
		E0[1] = w[0]*w[1]*(1-cos(th)) - w[2]*sin(th);
		E0[2] = w[2]*w[0]*(1-cos(th)) + w[1]*sin(th);
		
		E0[3] = w[0]*w[1]*(1-cos(th)) + w[2]*sin(th);
		E0[4] = cos(th) + w[1]*w[1]*(1-cos(th));
		E0[5] = w[2]*w[1]*(1-cos(th)) - w[0]*sin(th);
		
		E0[6] = w[2]*w[0]*(1-cos(th)) - w[1]*sin(th);
		E0[7] = w[2]*w[1]*(1-cos(th)) + w[0]*sin(th);
		E0[8] = cos(th) + w[2]*w[2]*(1-cos(th));
	}
	
	delete [] w;
}

// --- EOF ---
