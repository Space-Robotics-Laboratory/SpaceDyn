//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
//
// Function : gravity( g_set, Gravity )
//            set the gravity -- default is No Gravity [ 0 0 0 ]            
//
// Input  : g_set 3 x 1
// Output : Gravity 6 x 1 : to adapt with spatial notation
//
// s.abiko [2007.5]
//
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/rot.h"
#include "../include/spn.h"
#include "../include/spd.h"

void gravity( double g_set, double *Gravity ){
  
  matrix_Z(6, 1, Gravity);
//   if(g_set){
//     for(i=0;i<3;i++)
//       Gravity[i] = g_set[i];
//   }
}

// --- EOF ---
