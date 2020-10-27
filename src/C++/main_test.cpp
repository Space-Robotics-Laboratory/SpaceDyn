//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
//
// Main Function : dynamic simulator 
//
//
// s.abiko [2007.5]
//
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
#include <iostream>
#include <ctime>
using namespace std;

#include "./spd/include/spd_m.h"
#include "./spd/include/rot.h"
#include "./spd/include/spn.h"
#include "./spd/matrix/matrix.h"
#include "./spd/matrix/vector.h"

static MODEL m;
static double *Gravity;

int main(void)
{
  
  // set the gravity
  Gravity = matrix_get(6, 1);
  matrix_Z(6, 1, Gravity);

  // ------------------------------------
  // call the model parameters
  // ---- example 1 LINKNUM 4 = base + 3 links
  model("./spd/model_sample/TECSAS.def", m );
  model_init( m );

  calc_SP( m );

  f_kin_e( m, 1 );

 matrix_print( 3, 1, m.POS_e[1] );

 matrix_print( 3, 3, m.ORI_e[1] );

  //  --- memory clear ---  
  delete [] Gravity;
    
  return 0;

}

// --- EOF ---

