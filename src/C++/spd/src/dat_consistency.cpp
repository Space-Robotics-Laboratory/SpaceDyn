//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
//
// Function : dat_consistency( MODEL(original), MDOEL(new model) );
//            store the state values in the original model storage
//
//
// s.abiko [2008.10]
//
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/rot.h"
#include "../include/spn.h"
#include "../include/spd.h"

void dat_consistency( MODEL &m, int *DD, MODEL &new_m, int flag ){
	
	int i, j;
	
	// copy the data ( original model m -> new model new_m )
	if( flag == 0 ){
	
		model_init( new_m );
		 
		matrix_cpy( 4, 1, m.Qtn0,  new_m.Qtn0 );
		matrix_cpy( 4, 1, m.dQtn0, new_m.dQtn0 );
		matrix_cpy( 4, 1, m.A0, new_m.A0 );
		
		matrix_cpy( 3, 1, m.POS0, new_m.POS0 );
		
		matrix_cpy( 3, 1, m.v0, new_m.v0 );
		matrix_cpy( 3, 1, m.w0, new_m.w0 );
		matrix_cpy( 3, 1, m.vd0, new_m.vd0 );
		matrix_cpy( 3, 1, m.wd0, new_m.wd0 );

		matrix_cpy( 3, 1, m.F0, new_m.F0 );
		matrix_cpy( 3, 1, m.T0, new_m.T0 );
	
		if( new_m.LINKNUM != 1 ){
			int num_link[new_m.LINKNUM]; // num_link[0] invalid;
			// initialize 
			for(i=0;i<new_m.LINKNUM;i++)
				num_link[i] = -1;
	
			j = 0;
			for(i=1;i<m.LINKNUM;i++){
				if(DD[i] != 1){
					j++;
					num_link[j] = i;
				}
			}		
			for(i=1;i<new_m.LINKNUM;i++){
				if(num_link[i] != -1){
					new_m.q[i] = m.q[num_link[i]];
					new_m.qd[i] = m.qd[num_link[i]];
					new_m.qdd[i] = m.qdd[num_link[i]];
				
					new_m.tau[i] = m.tau[num_link[i]];
				
					cout << "i= " << i << ", num_link[i]=" << num_link[i] << " q = " << new_m.q[i] << endl; 

		  			matrix_cpy( 3, 1, m.Fe[num_link[i]], new_m.Fe[i] );
	  				matrix_cpy( 3, 1, m.Te[num_link[i]], new_m.Te[i] );
	  			}
			}
		}
	}
	// copy the data ( new model new_m -> original model m )
	else{
		matrix_cpy( 4, 1, new_m.Qtn0, m.Qtn0 );
		matrix_cpy( 4, 1, new_m.dQtn0, m.dQtn0 );
		matrix_cpy( 4, 1, new_m.A0, m.A0 );
		
		matrix_cpy( 3, 1, new_m.POS0, m.POS0 );
		
		matrix_cpy( 3, 1, new_m.v0, m.v0 );
		matrix_cpy( 3, 1, new_m.w0, m.w0 );
		matrix_cpy( 3, 1, new_m.vd0, m.vd0 );
		matrix_cpy( 3, 1, new_m.wd0, m.wd0 );

		matrix_cpy( 3, 1, new_m.F0, m.F0 );
		matrix_cpy( 3, 1, new_m.T0, m.T0 );
	
		if( m.LINKNUM != 1 ){
			int num_link2[m.LINKNUM]; // num_link[0] invalid;
			// initialize 
			for(i=0;i<m.LINKNUM;i++)
				num_link2[i] = -1;
	
			j = 0;
			for(i=1;i<m.LINKNUM;i++){
				if(DD[i] == 0){
					j++;
					num_link2[i] = j;
				}
			}		
			for(i=1;i<m.LINKNUM;i++){
				if(num_link2[i] != -1){
					m.q[i] = new_m.q[num_link2[i]];
					m.qd[i] = new_m.qd[num_link2[i]];
					m.qdd[i] = new_m.qdd[num_link2[i]];
					
					m.tau[i] = new_m.tau[num_link2[i]];

		  			matrix_cpy( 3, 1, new_m.Fe[num_link2[i]], m.Fe[i] );
	  				matrix_cpy( 3, 1, new_m.Te[num_link2[i]], m.Te[i] );
				}
				
				//cout << "i= " << i << ", num_link2[i]=" << num_link2[i] << " m.q = " << m.q[i] << endl; 

			}
		}
	}
	
}

// === EOF ===
