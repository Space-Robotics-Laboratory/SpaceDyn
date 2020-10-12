//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
//
// Function : ch_conf( MODEL, int, MODEL )
//            change the model for the case the joints are locked
//
// Output : Model new_m
//
// s.abiko [2008.4]
//
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/rot.h"
#include "../include/spn.h"
#include "../include/spd.h"

void ch_conf( MODEL &m, int *DD, MODEL &new_m ){

	int i,j,k,l;
	int ck_DD=0, e_num, new_e_num;
			
	double inv_mass;
	
	double *mass_vec;
	Matrix3 *inertia_vec, *inertia_MM;
	double *a_vec, *b_vec;
	double *a_mat, *a_mat_T, *b_mat, *b_mat_T;
	double *link_BBm, *inertia_BBi, *link_Im, *link_I, *inertia_i;
	
	mass_vec = matrix_get( m.LINKNUM, 1 );
	inertia_vec = new Matrix3[m.LINKNUM];
	inertia_MM = new Matrix3[m.LINKNUM];
	
	for( i=0; i<m.LINKNUM; i++)
		matrix_Z( 3, 3, inertia_MM[i]);
		
	a_vec = matrix_get( 3, 1 );
	b_vec = matrix_get( 3, 1 );
	a_mat = matrix_get( 3, 3 );
	a_mat_T = matrix_get( 3, 3 );
	b_mat = matrix_get( 3, 3 );
	b_mat_T = matrix_get( 3, 3 );
	
	link_BBm = matrix_get( 3, 3 );
	inertia_BBi = matrix_get( 3, 3 );
	link_Im = matrix_get( 3, 3 );
	link_I = matrix_get( 3, 3 );
	inertia_i = matrix_get( 3, 3 );
	

	double *AA, *AA_tmp, *DC, *AA_i, *AA_i_T;
	double *tmp, *tmp1, *rg_i, *rj_i, *Rc_i;
	double *MM1_i, *MM1_BB_i, *MM1;
	double *SS;
	
	AA = matrix_get( 3, 3 );
	AA_tmp = matrix_get( 3, 3 );
	DC = matrix_get( 3, 3 );
	AA_i = matrix_get( 3, 3 );
	AA_i_T = matrix_get( 3, 3 );
	
	matrix_I( 3, AA );
	matrix_I( 3, AA_tmp );
	matrix_I( 3, DC );
	matrix_I( 3, AA_i );
	matrix_I( 3, AA_i_T );
	
	tmp  = matrix_get( 3, 1 );
	tmp1  = matrix_get( 3, 1 );
	rg_i = matrix_get( 3, 1 );
	rj_i = matrix_get( 3, 1 );
	Rc_i = matrix_get( 3, 1 );

	MM1_i  = matrix_get( 3, 1 );
	MM1_BB_i = matrix_get( 3, 1 );
	MM1 = matrix_get( 3, 1 );
	
	SS = matrix_get( m.LINKNUM, m.LINKNUM );
	
	Vector3 *CtoJ, *JtoC, *CtoE;
	CtoJ = new Vector3[m.LINKNUM];
	JtoC = new Vector3[m.LINKNUM];
	CtoE = new Vector3[m.E_NUM];
	
	// !!! DD[0] invalid !!!
	new_m.LINKNUM = 1; // base body 
	for(i=1;i<m.LINKNUM;i++ ){
		if(DD[i] == 0 )
			new_m.LINKNUM++;
	}
	new_m.constructor();
	
	int num_link[new_m.LINKNUM]; // num_link[0] invalid;	
	j = 0;
	for(i=1;i<m.LINKNUM;i++){
		if(DD[i] == 0){
			j++;
			num_link[j] = i;
			//cout << "num_link[j], j = " << j << ", i = " << i << endl; 
		}
	}		
	
	e_num = 0;
	for( i=1; i<m.LINKNUM;i++){
		if( m.EE[i] >= e_num ){
			e_num = m.EE[i];
		} 
	}
 	int EE_link[e_num+1];
	for( i=1; i<m.LINKNUM;i++){
		if( m.EE[i] != 0 ){
			j = m.EE[i];
			EE_link[j] = i; // EE_link[0] invalid
		} 
	}
	
 	// ==== initialitation ====
	for(j=0;j<new_m.LINKNUM;j++){
		new_m.BB[j] = 0;		
		new_m.EE[j] = 0;		
	}
	//_/_/_/_/_/ BB /_/_/_/_/_/_/
	for(j=1;j<new_m.LINKNUM;j++){
		if( m.BB[num_link[j]] != 0){	
			k = m.BB[num_link[j]];
			ck_DD = DD[k];
			while( ck_DD != 0 ){
				k = m.BB[k];
				ck_DD = DD[k];
			}
			if( k != 0 ){
				for( l=1;l<new_m.LINKNUM;l++){
					if( num_link[l] == k ){
						new_m.BB[j] = l;	
					}
				}
			}
			else
				new_m.BB[j] = 0;
		}
		else
		new_m.BB[j] = 0;
	}	
	
	//_/_/_/_/_/ EE /_/_/_/_/_/_/
	int new_EE_link[new_m.LINKNUM];
	for( j=0; j<new_m.LINKNUM; j++ ){
		new_EE_link[j] = 0;
	}

	new_e_num = 1;
	for(i=1;i<e_num+1;i++){
		
		k = EE_link[i];
		while( ( DD[k] > 0 ) && ( m.BB[k] > 0 )){
			k = m.BB[k];
		}
		
		if( k != 0 ){
			for( j=1;j<new_m.LINKNUM;j++ ){
				if( num_link[j] == k ){
					new_EE_link[j] = EE_link[i];	
					cout << " new_EE_link[" << j << "]= "<< new_EE_link[j] << endl;
					new_m.EE[j] = 1;	
					new_e_num++;
				}
			}	
		}	
	}
	new_m.E_NUM = new_e_num-1;
	
	new_e_num = 0;
	for( j=1; j<new_m.LINKNUM; j++){
		new_e_num = new_e_num + new_m.EE[j];
		if( new_m.EE[j] != 0 ){
			new_m.EE[j] = new_e_num;
		}
	}
	new_m.construct_ee( new_m.E_NUM );
	
	//_/_/_/_/_/ J_type /_/_/_/_/_/_/
	for( j=1;j<new_m.LINKNUM;j++ ){
		new_m.J_type[j] = m.J_type[ num_link[j] ];
	}
	
	//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
	// to check the input for new model
	cout << " new LINKNUM = " << new_m.LINKNUM << endl;
	cout << " new_m_BB = [ " ;
	for(j=1;j<new_m.LINKNUM;j++){
		cout << new_m.BB[j] << " " ;
	}
	cout <<"] \n" <<  endl;
	
	cout << " new e_num = " << new_m.E_NUM << endl;
	cout << " new_m_EE = [ " ;
	for(j=1;j<new_m.LINKNUM;j++){
		cout << new_m.EE[j] << " " ;
	}
	cout <<"] \n" <<  endl;


	cout << " new_m_J_type = [ " ;
	for(j=1;j<new_m.LINKNUM;j++){
		cout << new_m.J_type[j] << " " ;
	}
	cout <<"] \n" <<  endl;
	
	//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/	
	for( j=1; j<new_m.LINKNUM; j++){
		matrix_Z( 3, 1, new_m.Qi[j] ); // Qi[0] is invalid
	}
	for( j=0; j<new_m.E_NUM; j++){
		matrix_Z( 3, 1, new_m.Qe[j] );
	}
	//_/_/_/_/_/ Qi /_/_/_/_/_/
	for( j=1; j<new_m.LINKNUM; j++){
		
		matrix_I( 3, AA );
		matrix_I( 3, AA_tmp );
		matrix_I( 3, DC );

		k = num_link[j];
			
		rpy2dc( m.Qi[k], AA_tmp );
		matrix_trans( 3, 3, AA_tmp, AA_tmp);			
		matrix_mult( 3, 3, 3, AA_tmp, AA, AA );
		
		while( DD[m.BB[k]] > 0){
			//cout << "k_inside = " << k << endl;
			
			k = m.BB[k];

			matrix_ext_sub( 6, 6, 1, 3, 1, 3, m.Xup[k], AA_tmp); // includes joint angle	
//			rpy2dc( m.Qi[k], AA_tmp );
//			matrix_trans( 3, 3, AA_tmp, AA_tmp);			
			matrix_mult( 3, 3, 3, AA_tmp, AA, AA );
			 
		}
		
		matrix_trans( 3, 3, AA, DC );
		dc2rpy( DC, new_m.Qi[j] );
		
	}
	
	for( i=1 ; i<new_m.LINKNUM ; i++ )
	{
	  cout << "rpy"<< i << " = [ " <<  new_m.Qi[i][0] << " " << new_m.Qi[i][1] << " " << new_m.Qi[i][2] << " ]"<< endl;
	}
	cout << endl;;  
		
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
	//_/_/_/_/_/ Qe /_/_/_/_/_/	
	for( i=1; i<m.E_NUM+1; i++ ){
			
		matrix_I( 3, AA );
		matrix_I( 3, AA_tmp );
		matrix_I( 3, DC );
	
		k = EE_link[i];

		//cout << "k = " << k << endl;
		//matrix_print( 3, 1, m.Qe[m.EE[k]-1]);
		rpy2dc( m.Qe[m.EE[k]-1], AA_tmp );
		matrix_trans( 3, 3, AA_tmp, AA_tmp);			
		matrix_mult( 3, 3, 3, AA_tmp, AA, AA );

		while( DD[k] > 0 ){
								
			matrix_ext_sub( 6, 6, 1, 3, 1, 3, m.Xup[k], AA_tmp); // includes joint angle		
//			rpy2dc( m.Qi[k], AA_tmp );
//			matrix_trans( 3, 3, AA_tmp, AA_tmp);			
			matrix_mult( 3, 3, 3, AA_tmp, AA, AA );

			k = m.BB[k];
		}
				
		//cout << "k = " << k << endl;
		for( j=1; j<new_m.LINKNUM; j++ ){
			if( num_link[j] == k ){
				matrix_trans( 3, 3, AA, DC );
				dc2rpy( DC, new_m.Qe[new_m.EE[j]-1] );
				cout << " new_m.EE[j] = " << new_m.EE[j] << endl;  
				//matrix_print( 3, 1, new_m.Qe[new_m.EE[j]-1] );
			}
		}		
	}

	for( i=0 ; i<new_m.E_NUM; i++ ){
	  cout << " rpy_e"<< i+1 << " = [ " <<  new_m.Qe[i][0] << " " << new_m.Qe[i][1] << " " << new_m.Qe[i][2] << " ]"<< endl;
	}
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

	//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/	
	// mass, inertia copy 
	for( i=0; i<m.LINKNUM;i++ ){
		mass_vec[i] = m.link_M[i];
		matrix_cpy( 3, 3, m.link_I[i], inertia_vec[i] );
	}
	for( i=m.LINKNUM-1;i>0;i-- ){
		if( DD[i] != 0 ) 
			mass_vec[m.BB[i]] = mass_vec[m.BB[i]] + mass_vec[i];
	
	}
	//_/_/_/_/_/ mass /_/_/_/_/_/_/
	new_m.link_M[0] = mass_vec[0];
	for(i=1;i<m.LINKNUM;i++ ){
		for( j=1; j<new_m.LINKNUM;j++){
			if( num_link[j] == i ) 
				new_m.link_M[j] = mass_vec[i];	
		}
	}

	cout << " new_m_link_M = [ " ;
	for(j=0;j<new_m.LINKNUM;j++){
		cout << new_m.link_M[j] << " " ;
	}
	cout <<"] \n" <<  endl;	


	//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/	
	//_/_/_/_/_/_/ JtoC, CtoJ, CtoE /_/_/_/_/_/_/	
	matrix_Z( 3, 1, m.JtoC[0] );
	matrix_Z( 3, 1, m.CtoJ[0] );
	matrix_Z( 3, 1, JtoC[0] );
	matrix_Z( 3, 1, CtoJ[0] );
	matrix_Z( 3, 1, m.Qi[0] ); // basically it is not used for definition
	
	for( i=1; i<m.LINKNUM; i++ ){
		matrix_cpy( 3, 1, m.CtoJ[i], CtoJ[i] );
		matrix_cpy( 3, 1, m.JtoC[i], JtoC[i] );	
		
		for( l=1; l<new_m.LINKNUM; l++){
			if( num_link[l] == i )
				// XXXXX strange XXXX
				matrix_sub( 3, 1, CtoJ[i], new_m.JtoC[new_m.BB[l]], new_m.CtoJ[l] );
		}
	}
	cout << " E_NUM = " << m.E_NUM << endl;
	for( i=0; i<m.E_NUM; i++){
		matrix_cpy( 3, 1, m.CtoE[i], CtoE[i] );
	}
		
	//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
	// mass copy 
	for( i=0; i<m.LINKNUM;i++ )
		mass_vec[i] = m.link_M[i];
	
	//_/_/_/_/_/ Center of Mass, JtoC /_/_/_/_/_/_/_/
	for( i=m.LINKNUM-1; i>=0; i--){		
		matrix_Z( 3, 1, MM1 );		
		
		if( DD[i] != 0 ){
			matrix_scale( 3, 1, mass_vec[m.BB[i]], JtoC[m.BB[i]], MM1_i ); // first order momentum parent link
			
			matrix_ext_sub( 6, 6, 1, 3, 1, 3, m.Xup[i], AA_i);		
			matrix_mult( 3, 3, 1, AA_i, JtoC[i], JtoC[i] ); // w.r.t. parent link frame BB[i]
			
			matrix_add( 3, 1, JtoC[i], m.CtoJ[i], tmp1 );
			cout << " i = " << i << endl;			
			
			matrix_print( 3, 1, m.JtoC[m.BB[i]] );
			matrix_add( 3, 1, tmp1, m.JtoC[m.BB[i]], JtoC[i] );						
			matrix_print( 3, 1, JtoC[i] );
			
			matrix_scale( 3, 1, mass_vec[i], JtoC[i], MM1 ); // first order momentum w.r.t. parent link frame BB[i]
			matrix_add( 3, 1, MM1_i, MM1, MM1 );	
			

			mass_vec[m.BB[i]] = mass_vec[m.BB[i]] + mass_vec[i];
			inv_mass = 1 / mass_vec[m.BB[i]];		
	
			//_/_/_/_/_/_/_/
			matrix_scale( 3, 1, inv_mass, MM1, JtoC[m.BB[i]] ); // Vector from Joint to CoM of composite link BB[i] w.r.t.BB[i]
			//_/_/_/_/_/_/_/
		}
	}
	//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
	for( i=m.LINKNUM-1; i>=0; i--){		
		for( l=1; l<new_m.LINKNUM; l++){
			if( num_link[l] == i ){
				matrix_cpy( 3, 1, JtoC[i], new_m.JtoC[l] );
			}
		}
	}	
	matrix_cpy( 3, 1, JtoC[0], new_m.JtoC[0] );
	//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
	
	
	for( i=m.LINKNUM-1; i > 0; i-- ){	
		
		//_/_/_/_/_/_/_/_/ CtoE /_/_/_/_/_/_/_/
		if( m.EE[i] != 0 ){
		
			j = i;
			matrix_add( 3, 1, m.JtoC[i], m.CtoE[m.EE[i]-1], CtoE[m.EE[i]-1] ); // w.r.t. BB[i] from the origin of i		
		
			while( DD[j] > 0 ){
			
				matrix_ext_sub( 6, 6, 1, 3, 1, 3, m.Xup[j], AA_i);		
							
				//rpy2dc( m.Qi[j], AA_i_T );
				//matrix_trans( 3, 3, AA_i_T, AA_i );
				matrix_mult( 3, 3, 1, AA_i, CtoE[m.EE[i]-1], CtoE[m.EE[i]-1] ); // w.r.t. BB[i] from the origin of i
			
				matrix_add( 3, 1, m.JtoC[m.BB[j]], m.CtoJ[j], tmp1 );
				matrix_add( 3, 1, tmp1, CtoE[m.EE[i]-1], CtoE[m.EE[i]-1] );
				j = m.BB[j];
			}				
		
			for( l=1; l<new_m.LINKNUM; l++ ){
				if( num_link[l] == j ){
					matrix_sub( 3, 1, CtoE[m.EE[i]-1], new_m.JtoC[l], new_m.CtoE[new_m.EE[l]-1] );
				}
			}
		} 		
		
		
		if( DD[m.BB[i]] != 0){
			//_/_/_/_/_/_/_/_/ CtoJ /_/_/_/_/_/_/_/_/
			matrix_add( 3, 1, m.JtoC[m.BB[i]], CtoJ[i], CtoJ[i] ); // w.r.t. BB[i] from the origin of i	
			
			j = i;
				
			while( DD[m.BB[j]] > 0 ){
				j = m.BB[j];
			
				matrix_ext_sub( 6, 6, 1, 3, 1, 3, m.Xup[j], AA_i);
				matrix_mult( 3, 3, 1, AA_i, CtoJ[i], CtoJ[i] ); // w.r.t. BB[i] from the origin of i
			
				matrix_add( 3, 1, m.JtoC[m.BB[j]], CtoJ[j], tmp1 );
				matrix_add( 3, 1, tmp1, CtoJ[i], CtoJ[i] );
				
			}
				for( l=1; l<new_m.LINKNUM; l++){
					if( num_link[l] == i ){
						matrix_sub( 3, 1, CtoJ[i], new_m.JtoC[new_m.BB[l]], new_m.CtoJ[l] );
				}
			}	
		}
		else{
			for( l=1; l<new_m.LINKNUM; l++){
				if( num_link[l] == i ){
					matrix_add( 3, 1, m.JtoC[m.BB[i]], CtoJ[i], CtoJ[i] );
					matrix_sub( 3, 1, CtoJ[i], new_m.JtoC[new_m.BB[l]], new_m.CtoJ[l] );
				}
			}
		}			
	}
	matrix_Z( 3, 1, CtoJ[0] );
	matrix_cpy( 3, 1, CtoJ[0], new_m.CtoJ[0] );
	
	//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/	
	cout << endl;
	cout << " new_m_link_CtoJ & JtoC & CtoE " << endl;
	cout << " CtoJ_0 = [ " <<  new_m.CtoJ[0][0] << " " << new_m.CtoJ[0][1] << " " << new_m.CtoJ[0][2] << " ]"<< endl;
	cout << " JtoC_0 = [ " <<  new_m.JtoC[0][0] << " " << new_m.JtoC[0][1] << " " << new_m.JtoC[0][2] << " ]"<< endl;
	cout << endl;
	for(j=1;j<new_m.LINKNUM;j++){
		  cout << " CtoJ_"<< j << " = [ " <<  new_m.CtoJ[j][0] << " " << new_m.CtoJ[j][1] << " " << new_m.CtoJ[j][2] << " ]"<< endl;
		  cout << " JtoC_"<< j << " = [ " <<  new_m.JtoC[j][0] << " " << new_m.JtoC[j][1] << " " << new_m.JtoC[j][2] << " ]"<< endl;
	  	  cout << endl;
		 if( new_m.EE[j] != 0 ){
			cout << " CtoE_"<< new_m.EE[j] << " = [ " <<  new_m.CtoE[new_m.EE[j]-1][0] << " " <<
			new_m.CtoE[new_m.EE[j]-1][1] << " " << new_m.CtoE[new_m.EE[j]-1][2] << " ]"<< endl << endl;;
		 }
	}
	//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/	
	
	
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/ Moment of Inertia /_/_/_/_/_/_/_/_/

	for( i=m.LINKNUM-1;i>=0;i-- ){
	
		//I = A*Ig*A' + m*tilde(r)T*tilde(r);
		if( DD[i] != 0 ){						
			matrix_ext_sub( 6, 6, 1, 3, 1, 3, m.Xup[i], AA_i);		
			matrix_trans( 3, 3, AA_i, AA_i_T );

			//rpy2dc( m.Qi[i], AA_i_T );
			//matrix_trans( 3, 3, AA_i_T, AA_i );
						
			matrix_mult( 3, 3, 3, AA_i, inertia_vec[i], link_I ); // recursive calculation
			matrix_mult( 3, 3, 3, link_I, AA_i_T, inertia_i );

			matrix_add( 3, 3, inertia_vec[m.BB[i]], inertia_i, inertia_vec[m.BB[i]] ); // composite inertia matrix
		}
				
		j = i;
		matrix_cpy( 3, 1, m.JtoC[j], tmp );
		while( DD[j] > 0 ){
			matrix_ext_sub( 6, 6, 1, 3, 1, 3, m.Xup[j], AA_i);		
			matrix_trans( 3, 3, AA_i, AA_i_T );

			//rpy2dc( m.Qi[j], AA_tmp );
			//matrix_trans( 3, 3, AA_tmp, AA );
			matrix_mult( 3, 3, 1, AA_i, tmp, tmp );
			
			matrix_add( 3, 1, tmp, m.CtoJ[j], tmp );
			matrix_add( 3, 1, tmp, m.JtoC[m.BB[j]], tmp );
		
			j = m.BB[j];
		}
		
		matrix_sub( 3, 1, JtoC[j], tmp, b_vec );
		tilde( 3, b_vec, b_mat );
		matrix_trans( 3, 3, b_mat, b_mat_T );
		matrix_mult( 3, 3, 3, b_mat_T, b_mat, link_BBm );
		matrix_scale( 3, 3, m.link_M[i], link_BBm, link_BBm );
		
		matrix_add( 3, 3, link_BBm, inertia_MM[j], inertia_MM[j] );

	} 
		
	matrix_add( 3, 3, inertia_vec[0], inertia_MM[0], inertia_vec[0] );
	matrix_cpy( 3, 3, inertia_vec[0], new_m.link_I[0]);
	for( i=1;i<m.LINKNUM;i++ ){
		matrix_add( 3, 3, inertia_vec[i], inertia_MM[i], inertia_vec[i] );
		for( j=1; j<new_m.LINKNUM;j++){
			if( num_link[j] == i ){
				matrix_cpy( 3, 3, inertia_vec[i], new_m.link_I[j]);
			}
		}	
	} 
	
	
	for( j=0; j<new_m.LINKNUM;j++){	
		cout << "link_I[" << j << "]" << endl;
		matrix_print( 3, 3, new_m.link_I[j] ); 
	}

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//  calc_SPN( new_m );  
/* /_/_/_/_/_/_/_/_/_/_/ calc_SPN  /_/_/_/_/_/_/_/_/_/_/_/_/
// calculate Spatial Notation       */
// convert all parameters to Spatial Notation
// basically it same as the function "calc_SPN", except for calculation of new_m.Isp[0]
// caution JtoC[0] is the center of mass of new base, new_m.JtoC[0] is the center of mass of the new base
// caution new_m.CtoJ[1] is the vector from new center of mass of the base to the joint 1
  
  double *zero3;
  double *LL;
  double *A_BB_i, *A_BB_e, *trans_tmp;

  // ---- initialization ----
  zero3 = matrix_get(3,1);
  LL = matrix_get(3,1);
  A_BB_i = matrix_get(6,6);
  A_BB_e = matrix_get(6,6);
  trans_tmp = matrix_get(6,6);

  // ---- base ----
  matrix_I(6,new_m.Xsp[0]);
  matrix_Z(3,1,zero3);
  
// spI( new_m.link_M[0], zero3, new_m.link_I[0], new_m.Isp[0] ); // original calc_SPN
  spI( new_m.link_M[0], JtoC[0], new_m.link_I[0], new_m.Isp[0] ); // modification
  
  // --- link ---- if multibody system ----
  if( new_m.LINKNUM != 1){
    for( i=1; i<new_m.LINKNUM; i++){
      // each link's parameters
      XrotBB( new_m.Qi[i], A_BB_i );
	matrix_add( 3, 1, new_m.JtoC[new_m.BB[i]], new_m.CtoJ[i], LL );


      // ---- spatial transformation matrix i+1 -> i ---
      Xtrans(LL, trans_tmp);
      matrix_mult( 6, 6, 6, trans_tmp, A_BB_i, new_m.Xsp[i]);
      // m.Xsp[i] = Xtrans(LL)*rotE;
      
      // ---- spatial transformation matrix CoM i -> Joint i ----
      Xtrans(new_m.JtoC[i], new_m.jXc[i]);
      
      if(new_m.EE[i] != 0){
	XrotBB( new_m.Qe[new_m.EE[i]-1], A_BB_e );
	matrix_add( 3, 1, new_m.JtoC[i], new_m.CtoE[new_m.EE[i]-1], LL);
	//   LL = LP.ce(:,i) - LP.cc(:,i,i); link length for last link
	
	// ---- spatial transformation matrix EE -> Joint i ---
	Xtrans(LL, trans_tmp);
	matrix_mult( 6, 6, 6, trans_tmp, A_BB_e, new_m.jXe[i]);
	// m.jXe[i] = Xtrans(LL)*rotE;
      }
      else{
	matrix_I(6, new_m.jXe[i]);
      }
      
      // rigid body inertia MF6
      spI( new_m.link_M[i], new_m.JtoC[i], new_m.link_I[i], new_m.Isp[i] );  
      
     // joint type matrix
      if(new_m.J_type[i] != 1)  // revolute joint
	new_m.S[i][5] = 1.0;
      else // prismatic joint
	new_m.S[i][2] = 1.0; 
      //matrix_print(6,1,m.S[i-1]);
    }
  }

  delete [] zero3;
  delete [] LL;
  delete [] A_BB_i;
  delete [] A_BB_e;
  delete [] trans_tmp;

 
  cout << " /_/_/_/ model: Spatial Notation Convert finished /_/_/_/ " << endl;
   

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

	// clear memory 
	delete [] AA;
	delete [] AA_tmp;
	delete [] DC;
	delete [] AA_i;
	delete [] AA_i_T;
	
	delete [] tmp;
	delete [] tmp1;
	
	delete [] rg_i;
	delete [] rj_i;
	delete [] Rc_i;
	
	delete [] mass_vec;
	delete [] inertia_vec;
	delete [] inertia_MM;
	delete [] a_vec;
	delete [] b_vec;
	
	delete [] a_mat;
	delete [] a_mat_T;
	delete [] b_mat;
	delete [] b_mat_T;
	
	delete [] link_BBm;
	delete [] inertia_BBi;
	delete [] link_Im;
	delete [] link_I;
	delete [] inertia_i;	
	
	delete [] CtoJ;
	delete [] JtoC;
	delete [] CtoE;
	
	delete [] MM1_i;
	delete [] MM1_BB_i;
	delete [] MM1;
	
	delete [] SS;
	
}

// === EOF ===
