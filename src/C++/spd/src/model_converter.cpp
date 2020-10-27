//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
//
// Function : model_converter( string, string )
//            Convert the Model defined by Rainer to that by Satoko
//
//
// s.abiko [2007.11]
//
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
#include "../matrix/matrix.h"
#include "../matrix/vector.h"
#include "../include/rot.h"
#include "../include/spn.h"
#include "../include/spd.h"

#define _LINESKIP fin >> str;
#define _LINEOUT fout << "\n";

void model_converter( string file_in, string file_out ){

  int i, j;
  
  string str;
  string type[5000];
  double data[5000]; // need to change the size for arbitrary size

  int nLink=0, nEE=0; 
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
  // ---- open the file to read and to write ----
  ifstream fin(file_in.c_str());
  if( !fin )
    {
      cerr << "Cannot open " << file_in << endl;
      exit(1);
    }
  ofstream fout(file_out.c_str());
  if( !fout )
    {
      cerr << "Cannot open " << file_out << endl;
      exit(1);
    }
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
  // ===== read the file_in =====
  i = 0;
  _LINESKIP;
  cout << str <<" ";
  _LINESKIP;
  cout << str << endl;
  _LINESKIP;
  cout << str <<" ";
  _LINESKIP;
  cout << str << endl;
  _LINESKIP;
  cout << str <<" ";
  _LINESKIP;
  cout << str << endl;

  while( !fin.eof() ){
    	fin >> type[i] >> data[i];
	
	if( type[i].compare(0,5,"nLink") == 0 ){
		nLink = (int)data[i];
		cout << " n LINK = " << nLink << endl;
	}
	if( type[i].compare(0,5,"nLoop") == 0 ){
		nEE = (int)data[i]+1;
		cout << " n EE = " << nEE << endl;
	}	
	//	cout << "i = " << i << ", type = " << type[i] << ", data = " << data[i] << endl;
	i++;  
  }
  cout << " #### the def file has been read !! ####" << endl << endl;
  
 int *BB, *J_type, *EE, *EE_Link;
 double *mass, rpy_deg; 
 Vector3 *Link_L, *Link_Le, *CtoJ, *JtoC, *rpyQ;
 //Vector3 *CtoE, *rpyE;
 Vector3 *rpyE;
 
 int pp;
 double *CtoE[nEE];
 for( pp=0; pp<nEE; pp++){
 	CtoE[nEE] = new double[3];
 }
 
 Matrix3 *iAj;
 Matrix3 *iAe;
 Matrix3 *Inertia;
 
 double *tmp, *tmp1, *tmp2, *tmp3;
 tmp = matrix_get(3, 3);
 tmp1 = matrix_get(3, 3);
 tmp2 = matrix_get(3, 1);
 tmp3 = matrix_get(3, 3);
   
 BB = new int[nLink-1];
 J_type = new int[nLink-1];
 EE = new int[nLink];
 EE_Link = new int[nEE];
 
 Link_L = new Vector3[nLink];
 Link_Le = new Vector3[nEE];
 
 CtoJ = new Vector3[nLink];
 JtoC = new Vector3[nLink];
 rpyQ = new Vector3[nLink-1];
 
 //CtoE = new Vector3[nEE];
 rpyE = new Vector3[nEE];
 
 iAj = new Matrix3[nLink];
 iAe = new Matrix3[nEE];
 
 mass = new double[nLink];
 
 Inertia = new Matrix3[nLink];  
 
 int k=0, m=0, p=0;
 int nChar = 0;
 
 string Link, num_EE, num_Child, row, col, xyz, Child;
 int SS[nLink-1][nLink-1], pc[nLink-1][nLink-1];


  // ==== children links connectivity ====
  for( j=0; j<nLink-1; j++){
	  for( k=0; k<nLink-1; k++){
  		SS[j][k] = 0;
		pc[nLink-1][nLink-1] = -1;
	}
  }
  for( j=0; j<nLink-1; j++)
  	SS[j][j] = -1; 
	
  for( j=0; j<nLink; j++)
	EE[j] = 0;
  
  for( j=0; j<i; j++ ){
  	
	// ==== BB, SS ====
	if( type[j].compare(0,10,"Link.pLink") == 0){	
		nChar = 11;	
		Link = "";
		while( type[j].compare( nChar, 1, ")" ) != 0 ){
			Link = Link + type[j][nChar];
			nChar++;
		}
		//cout << " Link = "<< Link << endl;
		if( atoi( Link.c_str() ) != 1 ){
			BB[ atoi( Link.c_str() )-2 ] = (int)data[j]-1; 
			SS[ (int)data[j]-1 ][ atoi( Link.c_str() )-1 ] = 1;						
		//	cout << "######"<<type[j];
		//	cout << "######" << BB[ atoi( Link.c_str() )-2 ] << endl;
		}
	}	

	// ==== EE ====
	if( type[j].compare(0,10,"Loop.pLink") == 0){
		nChar = 11;	
		Link = "";
		while( type[j].compare( nChar, 1, "," ) != 0 ){
			Link = Link + type[j][nChar];
			nChar++;
		}
		nChar++;
		num_EE = "";
		while( type[j].compare( nChar, 1, ")" ) != 0 ){
			num_EE = num_EE + type[j][nChar];
			nChar++;
		}
		EE[ (int)data[j] -1 ] = atoi( num_EE.c_str() );
		EE_Link[ atoi( num_EE.c_str() )-1 ] = (int)data[j];
		//cout << "######"<<type[j] << " " ;
		//cout << (int)data[j] << "######" << EE[ (int)data[j]-1  ] << "######" << EE_Link[ atoi( num_EE.c_str() ) -1] << endl;
	}
	
  }
  
  // ==== cout the number of children connected to link j ====
  for( j=0; j<nLink-1; j++){
  	p = 0;
	for( k=0; k<nLink; k++){
	  	if( SS[j][k] == 1 ){
			pc[j][p] = k;
			cout << "p = " << p << j << ", PC = " << pc[j][p] << endl; 
			p++;
  		}
  	}
  }
  
  
  for( j=0; j<i; j++ ){
  
  	// ==== J_type ====
	if( type[j].compare(0,10,"Joint.type") == 0){
		nChar = 11;	
		Link = "";
		while( type[j].compare( nChar, 1, ")" ) != 0 ){
			Link = Link + type[j][nChar];
			nChar++;
		}
		//cout << " Link ( Joint type ) = " <<Link << endl;
		
		J_type[ atoi( Link.c_str() ) -1 ] = (int)data[j];
		//cout << type[j] <<" = " <<  J_type[atoi( Link.c_str() )-1 ] << endl;
	}
	
	// ==== mass ====
  	if( type[j].compare(0,6,"Link.m") == 0){
		nChar = 7;	
		Link = "";
		while( type[j].compare( nChar, 1, ")" ) != 0 ){
			Link = Link + type[j][nChar];
			nChar++;
		}
		//cout << " Link ( mass ) = " <<Link << endl;

		mass[ atoi( Link.c_str() )-1 ] = data[j];
		//cout << " Link mass " << atoi( Link.c_str() )-1 <<" = " << mass[atoi( Link.c_str() )-1] << endl;
	}
	
	// ==== Inertia ====
	if( type[j].compare(0,6,"Link.J") == 0){
		nChar = 7;	
		Link = "";
		while( type[j].compare( nChar, 1, "," ) != 0 ){
			Link = Link + type[j][nChar];
			nChar++;
		}
		nChar++;
		row = "";
		while( type[j].compare( nChar, 1, "," ) != 0 ){
			row = row + type[j][nChar];
			nChar++;
		}
		nChar++;
		col = "";
		while( type[j].compare( nChar, 1, ")" ) != 0 ){
			col = col + type[j][nChar];
			nChar++;
		}

		//cout << " Link ( Inertia ) = " << Link << " row = " << row << ", col = " << col << endl;

		Inertia[ atoi(Link.c_str())-1 ][ 3*( atoi(col.c_str())-1 ) + ( atoi(row.c_str()) -1 ) ] = data[j];
		//cout << "######"<<type[j];
		//cout << "######" << Inertia[ atoi(Link.c_str())-1 ][ 3*( atoi(row.c_str())-1 ) + ( atoi(col.c_str()) -1 ) ] << endl;							
	}
	
	// ==== JtoC ==== rCOG ====
	if( type[j].compare(0,9,"Link.rCOG") == 0){
		nChar = 10;	
		Link = "";
		while( type[j].compare( nChar, 1, "," ) != 0 ){
			Link = Link + type[j][nChar];
			nChar++;
		}
		nChar++;
		xyz = "";
		while( type[j].compare( nChar, 1, ")" ) != 0 ){
			xyz = xyz + type[j][nChar];
			nChar++;
		}
		// cout << " Link ( rCOG ) = " << Link << " xyz = " << xyz << endl;
		
//		if( atoi( Link.c_str() ) > 1 ){
			JtoC[ atoi( Link.c_str() )-1 ][ atoi(xyz.c_str())-1 ] = data[j];
			//cout << "######"<<type[j];
			//cout << "######" << JtoC[ atoi( Link.c_str() )-2 ][ atoi(xyz.c_str())-1 ]<< endl;		
//		}
	}

	// ==== rpy, CtoJ, rpyEE ==== Link.B ====
	if( type[j].compare(0,6,"Link.B") == 0 ){	
		nChar = 7;	
		Link = "";
		while( type[j].compare( nChar, 1, "," ) != 0 ){
			Link = Link + type[j][nChar];
			nChar++;
		}
		nChar++;	
		Child = "";
		while( type[j].compare( nChar, 1, "," ) != 0 ){
			Child = Child + type[j][nChar];
			nChar++;
		}
		nChar++;	
		row = "";
		while( type[j].compare( nChar, 1, "," ) != 0 ){
			row = row + type[j][nChar];
			nChar++;
		}
		nChar++;	
		col = "";
		while( type[j].compare( nChar, 1, ")" ) != 0 ){
			col = col + type[j][nChar];
			nChar++;
		}

		if( atoi( row.c_str() ) != 4 ){
			if( atoi( Link.c_str()) != 1 ){
				if( EE[ atoi( Link.c_str())-1] == 0 ){ // the link is not a last link
					m = pc[atoi( Link.c_str() )-1][atoi( Child.c_str() )-1]-1;
					if( atoi( col.c_str() ) != 4 ){
						iAj[m][ 3*( atoi(row.c_str())-1 ) + ( atoi(col.c_str()) -1 ) ] = data[j];
					}	
					else{
						Link_L[m][atoi(row.c_str())-1] = data[j];
						cout << "m = " << m <<  " Link_I_J " << Link_L[m][atoi(row.c_str())-1] << endl;
					}
				}
				else{ // the link is a last link, which has an end-effector
					m = EE[atoi( Link.c_str())-1]-1;
					if( atoi( col.c_str() ) != 4 ){
						iAe[m][ 3*( atoi(row.c_str())-1 ) + ( atoi(col.c_str()) -1 ) ] = data[j];
					}	
					else{
						Link_Le[m][atoi(row.c_str())-1] = data[j];
					//cout << Link_Le[m][atoi(row.c_str())-1] << endl;
					}
				}
			}
			else{
				m = 0;
				if( atoi( col.c_str() ) != 4 ){
					iAj[m][ 3*( atoi(row.c_str())-1 ) + ( atoi(col.c_str()) -1 ) ] = data[j];
				}	
				else{
					Link_L[m][atoi(row.c_str())-1] = data[j];
				}
			}
		}
		//cout << " Link ( Link.B ) = " << Link << " Child = " << Child << ", row = " << row << ", col = " << col << endl;
	}
 }
  
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
  // ==== write the model to the file_out ====
  fout << "###Parameters_for_MODEL(___DO_NOT_CHANGE_THIS_ORDER!!___)" << endl;
  _LINEOUT;
  
  fout << "###LINK_NUMBER " << nLink << endl;
  _LINEOUT;
  
  fout << "###LINK_CONNECTIVITY" << endl;
  fout << "BB[ ";
  	for( j=0; j<nLink-1; j++)
		fout << BB[j] << " ";
  fout << "]" << endl;
  _LINEOUT;
  
  fout << "###Joint_Type_[_0=rotational__1=prismatic]" << endl;
  fout << "J_type[ ";
  	for( j=0; j<nLink-1; j++)
		fout << J_type[j] << " ";
  fout << "]" << endl;;
  _LINEOUT;
  
  fout << "###EE" << endl;
  fout << "EE[ ";
  	for( j=1; j< nLink; j++)
		fout << EE[j] << " ";
  fout << "]" <<endl;
  _LINEOUT;
  
  fout << "###Relative_Coordinate_For_Link:Roll_Pitch_Yaw_Angle_of_Each_Bodies{x-y-z}_[deg]" << endl;
  for(j=1;j<nLink;j++){
  	
	matrix_trans( 3, 3, iAj[j-1], tmp );
	dc2rpy( tmp, rpyQ[j] );
	
//	cout << " iAj ( " << j << " )= " << endl;
//	matrix_print( 3, 3, iAj[j] );	
  	fout << "rpy" << j << "[ ";
	for( k=0; k<3; k++){
		rpy_deg = pi2deg( rpyQ[j][k] );
		fout << rpy_deg << " ";
	}
	fout << "]" << endl;
  }
  _LINEOUT;
  
  fout << "###Vector_Of_Link_Length_JtoC_0_0_is_always_[_0_0_0_]" << endl;
  for(j=1;j<nLink;j++){
  	fout << "CtoJ_" << BB[j-1] <<"_" << j <<"[ ";
	
	matrix_scale( 3, 1, -1, JtoC[BB[j-1]], tmp2 );
	matrix_add( 3, 1, Link_L[BB[j-1]], tmp2, tmp2 );
	
	cout << "Link_L (" << BB[j-1] << ", " << j << ") = "<< endl;
	matrix_print( 3, 1, Link_L[ BB[j-1]]);

	//cout << "JtoC (" << BB[j-1] << ", " << j << ") = "<< endl;
	//matrix_print( 3, 1, JtoC[BB[j-1]]);
	
	for( k=0; k<3; k++){
		fout << tmp2[k] << " "; 
	}
	fout << "]" << endl;
	
  	fout << "JtoC_" << j <<"_" << j <<"[ ";
	for( k=0; k<3; k++){
		fout << JtoC[j][k] << " "; 
	}
	fout << "]" << endl;
	_LINEOUT;
  }
  
  fout << "###Vector_To_End-Effector" << endl;
  for(j=0;j<nEE;j++){
  	fout << "CtoE_" << j+1 << "[ ";
	matrix_scale( 3, 1, -1, JtoC[EE_Link[j]-1], tmp2 );
	matrix_add( 3, 1, Link_Le[j], tmp2, tmp2 );
	
	for( k=0; k<3; k++){
		fout << tmp2[k] << " "; 
	}
	fout << "]" << endl;
  }
  _LINEOUT;	
  
  fout << "###Relative_Coordinate_For_End-Effector" << endl;
  for(j=0;j<nEE;j++){
	
	matrix_trans( 3, 3, iAe[j], tmp );
	dc2rpy( tmp, rpyE[j] );
	
  	fout << "rpyE_" << j+1 << "[ ";
	for( k=0; k<3; k++){
		rpy_deg = pi2deg( rpyE[j][k] );
		fout << rpy_deg << " ";
	}
	fout << "]" << endl;	
 }
  _LINEOUT;	
  
  fout << "###Mass_Parameter" << endl;
  fout << "mass[ ";
  for( j=0; j<nLink; j++)
  	fout << mass[j] << " ";
  fout << "]" << endl;
  _LINEOUT;
  	

  fout << "###Inertia_Matrix" << endl;
  fout << "###[_I11_I12_I13]" << endl;
  fout << "###[_I21_I22_I23]" << endl;
  fout << "###[_I31_I32_I33]" << endl;
  _LINEOUT;
  
  for(j=0;j<nLink;j++){

	tilde( 3, JtoC[j], tmp );
	matrix_trans( 3, 3, tmp, tmp1 );
	matrix_mult( 3, 3, 3, tmp1, tmp, tmp3 );
	matrix_scale( 3, 3, -	mass[j], tmp3, tmp ); // calculate the moment of inertia around CoM by using Parallel axis theorem
	
	matrix_add( 3, 3, Inertia[j], tmp, Inertia[j] );

  	fout << "###Link" << j << endl;
	fout << "I11= " << Inertia[j][0]<< endl;
	fout << "I22= " << Inertia[j][4]<< endl;
	fout << "I33= " << Inertia[j][8]<< endl;
	fout << "I12= " << Inertia[j][1]<< endl;
	fout << "I13= " << Inertia[j][2]<< endl;
	fout << "I23= " << Inertia[j][5]<< endl;
	_LINEOUT;
  }
  
  fout << "###EOF 777" << endl; 
  
  cout << " #### The def file is converted into the model file of library ! ####" << endl << endl;
  
  //cout << " ### memody check" << endl;
  delete [] BB;
  delete [] J_type;
  delete [] EE;
  delete [] EE_Link;
  
  delete [] mass;

  delete [] tmp;
  delete [] tmp1;
  delete [] tmp2;
  delete [] tmp3;
 //cout << " ### memody check finished" << endl;

}

// === EOF ===
