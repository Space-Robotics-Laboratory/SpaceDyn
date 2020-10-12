//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
//
// Function : w2dQtn( w, dQtn )
//            convert from angular velocity to time-derivative of the quotation 
//
// Output : dQtn 4 x 1 ( scalar 1x1, vector 3x1 )
//
// s.abiko [2007.11]
//
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/rot.h"
#include "../include/spn.h"
#include "../include/spd.h"

void w2dQtn( double *w, double *Qtn, double *dQtn){
  
  //matrix_Z( 4, 1, dQtn ); // test

  dQtn[0] = - ( Qtn[1]*w[0] + Qtn[2]*w[1] + Qtn[3]*w[2] ) / 2;
  
  dQtn[1] = ( Qtn[0]*w[0] + Qtn[3]*w[1] - Qtn[2]*w[2] ) / 2;
  dQtn[2] = (-Qtn[3]*w[0] + Qtn[0]*w[1] + Qtn[1]*w[2] ) / 2;
  dQtn[3] = ( Qtn[2]*w[0] - Qtn[1]*w[1] + Qtn[0]*w[2] ) / 2;

}

// === EOF === //
