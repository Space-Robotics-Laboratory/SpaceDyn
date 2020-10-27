//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//
//
// Function : functions for rotation (attitude transformation)
//            Utilities for RPY, Rotation Transform, 
//            Direct Cosine, Quaternion etc...
//
// s.abiko [2007.5]
//
// original files are in src produced by y.sato
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../include/rot_a.h"
#include "../include/common.h"

#define _R(i,j) R[(i)*3+(j)]
#define _R1(i,j) R1[(i)*3+(j)]
#define _R2(i,j) R2[(i)*3+(j)]

#define dSqrt(x) (sqrt(x))
#define REAL(x) (x)
#define dRecip(x) (1.0/(x))

void rpy2R ( double *rpy, double *R )
{
    double s1,s2,s3,c1,c2,c3;
    
    s1 = sin(rpy[2]);
    s2 = sin(rpy[1]);
    s3 = sin(rpy[0]);
    c1 = cos(rpy[2]);
    c2 = cos(rpy[1]);
    c3 = cos(rpy[0]);
    
    _R(0,0) = c1*c2;
    _R(0,1) = c1*s2*s3 - s1*c3;
    _R(0,2) = c1*s2*c3 + s1*s3;
    _R(1,0) = s1*c2;
    _R(1,1) = s1*s2*s3 + c1*c3;
    _R(1,2) = s1*s2*c3 - c1*s3;
    _R(2,0) = -s2;
    _R(2,1) = c2*s3;
    _R(2,2) = c2*c3;
    /*
    R[0]     = c1*c2;
    R[1]     = c1*s2*s3 - s1*c3;
    R[2]     = c1*s2*c3 + s1*s3;
    R[3*1]   = s1*c2;
    R[3*1+1] = s1*s2*s3 + c1*c3;
    R[3*1+2] = s1*s2*c3 - c1*s3;
    R[3*2]   = -s2;
    R[3*2+1] = c2*s3;
    R[3*2+2] = c2*c3;
    */
}

void R2rpy ( double *R, double *rpy )
{
    
    double co=0.0;
    
    rpy[1] = atan2( -_R(2,0),
		    sqrt( _R(0,0)*_R(0,0) + _R(1,0)*_R(1,0) ));
    
    if( rpy[1] == M_PI/2.0)
    {
	rpy[2] = 0.0;
	rpy[0] = atan2( _R(0,1), _R(1,1) );
    }
    else if( rpy[1] == -M_PI/2.0)
    {
	rpy[2] = 0.0;
	rpy[0] = -atan2( _R(0,1), _R(1,1) );
    }
    else
    {
	co = cos( rpy[1] );
	rpy[2] = atan2( _R(1,0)/co, _R(0,0)/co );
	rpy[0] = atan2( _R(2,1)/co, _R(2,2)/co);
    }
    
    /*
      double s2 = 0.0;
      double c2 = 0.0;
      
      if( ( fabs( _R(1,0) ) < 1e-15 ) && ( fabs( _R(0,0) ) < 1e-15 ) )
      {
      rpy[2] = 0.0;
      rpy[1] = atan2( _R(2,0), _R(0,0) );
      rpy[0] = atan2( _R(1,2), _R(1,1) );
      }
      else{
	
      rpy[2] = atan2( -_R(1,0), _R(0,0) );
      c2 = cos(rpy[2]);
      s2 = sin(rpy[2]);
      rpy[1] = atan2( _R(2,0), c2*_R(0,0)-s2*_R(1,0));
      rpy[0] = atan2( -_R(2,1), _R(2,2) );
      }
    */
}


void R2qtn( double *R, double *qtn )
{
    double tr,s;
    tr = _R(0,0) + _R(1,1) + _R(2,2);
    if (tr >= 0) {
	s = dSqrt (tr + 1);
	qtn[0] = REAL(0.5) * s;
	s = REAL(0.5) * dRecip(s);
	qtn[1] = (_R(2,1) - _R(1,2)) * s;
	qtn[2] = (_R(0,2) - _R(2,0)) * s;
	qtn[3] = (_R(1,0) - _R(0,1)) * s;
    }
    else {
	// find the largest diagonal element and jump to the appropriate case
	if (_R(1,1) > _R(0,0)) {
	    if (_R(2,2) > _R(1,1)) goto case_2;
	    goto case_1;
	}
	if (_R(2,2) > _R(0,0)) goto case_2;
	goto case_0;
	
    case_0:
	s = dSqrt((_R(0,0) - (_R(1,1) + _R(2,2))) + 1);
	qtn[1] = REAL(0.5) * s;
	s = REAL(0.5) * dRecip(s);
	qtn[2] = (_R(0,1) + _R(1,0)) * s;
	qtn[3] = (_R(2,0) + _R(0,2)) * s;
	qtn[0] = (_R(2,1) - _R(1,2)) * s;
	return;
	
    case_1:
	s = dSqrt((_R(1,1) - (_R(2,2) + _R(0,0))) + 1);
	qtn[2] = REAL(0.5) * s;
	s = REAL(0.5) * dRecip(s);
	qtn[3] = (_R(1,2) + _R(2,1)) * s;
	qtn[1] = (_R(0,1) + _R(1,0)) * s;
	qtn[0] = (_R(0,2) - _R(2,0)) * s;
	return;
	
    case_2:
	s = dSqrt((_R(2,2) - (_R(0,0) + _R(1,1))) + 1);
	qtn[3] = REAL(0.5) * s;
	s = REAL(0.5) * dRecip(s);
	qtn[1] = (_R(2,0) + _R(0,2)) * s;
	qtn[2] = (_R(1,2) + _R(2,1)) * s;
	qtn[0] = (_R(1,0) - _R(0,1)) * s;
	return;
    }
    
}

void rpy2qtn( double *rpy, double *qtn )
{
    double cr, cp, cy;
    double sr, sp, sy;
    
    cr = cos( rpy[2]/2 );
    sr = sin( rpy[2]/2 );
    cp = cos( rpy[1]/2 );
    sp = sin( rpy[1]/2 );
    cy = cos( rpy[0]/2 );
    sy = sin( rpy[0]/2 );
    
    qtn[0] = cr*cp*cy - sr*sp*sy;
    qtn[1] = cr*sy*sy + sr*cp*cy;
    qtn[2] = cr*sp*cy + sr*cp*sy;
    qtn[3] = cr*cp*sy - sr*sp*cy;
    
}


double deg2pi( double deg )
{
    return (deg/180.0)*M_PI;
}

double pi2deg( double pi )
{
    return (pi/M_PI)*180.0;
}

double pi_conv( char *data )
{
    
    int sign = 1,flag;
    int i = 0;
    double fact = 1.0;
    
    switch( *data )
    {
    case '+' : sign =  1; i = 3; break; 
    case '-' : sign = -1; i = 3; break; 
    case 'p' : sign =  1; i = 2; break;
    case '0' : sign =  0; i = 1; break;
    }
    
    switch( *(data + i) )
    {
    case '*' : flag = 1; break; 
    case '/' : flag = 2; break; 
    default  : return( sign * M_PI ); break;
    }
    
    if ( sscanf( data + i + 1, " %lf ", &fact ) != 1 )
    {
	printf("error : sscanf\n");
	exit(1);
	
    }
    
    switch( flag )
    {
    case 1 : return( sign * M_PI * fact );
    case 2 : return( sign * M_PI / fact );
    default: return 0;
    }
}

void QMultiply ( double *qb, double *qc, double *qa )
{
    qa[0] = qb[0]*qc[0] - qb[1]*qc[1] - qb[2]*qc[2] - qb[3]*qc[3];
    qa[1] = qb[0]*qc[1] + qb[1]*qc[0] + qb[2]*qc[3] - qb[3]*qc[2];
    qa[2] = qb[0]*qc[2] + qb[2]*qc[0] + qb[3]*qc[1] - qb[1]*qc[3];
    qa[3] = qb[0]*qc[3] + qb[3]*qc[0] + qb[1]*qc[2] - qb[2]*qc[1];
}
