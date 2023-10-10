#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
using namespace std;

//#include "../include/model.h"
#include "../include/define.h"
#include "../include/common.h"

void model(string, MODEL &);
void model_CL(string, MODEL &); // need to debug
void model_new(string, MODEL &); // not finished
void model_init(MODEL &);
void model_converter( string, string ); // Rainer def to Satoko def.

void calc_SPN(MODEL &);
void ch_conf(MODEL &, int *, MODEL &);
void dat_consistency(MODEL &, int *, MODEL &, int);

void calc_SP (MODEL &);
void calc_hh (MODEL &, double *);
void calc_Rg (MODEL &, double *);
void calc_JJ (MODEL &, double **);
void calc_Je (MODEL &, int, double * );
void calc_GJe (MODEL &, int, double * );
void calc_Jb (MODEL &, int, double * );
void calc_Lg( MODEL &, double *);
void f_kin_e (MODEL &, int );
void f_kin_j (MODEL &);
void i_dyn (MODEL &, double *, double *);
void f_dyn (MODEL &, double *, double *, double *, double *);

void calc_C (MODEL &, double *, double * );
//void f_dyn_rk (MODEL &m, double *, double, double *, double *, double *, double *, double *, double *, double *); // fixed step runge-kutta
void f_dyn_rk(MODEL &m, double *Gravity, double step);

void calc_Xup_I(MODEL &);
void calc_SS (MODEL &, double *);
void j_num(MODEL &m, int, int *, int, int );
void aw(double *, double, double *);

// simple explicit euler-method integrator
void int_eu(MODEL &, double *, double *, double *, double *, double *, double *);

// linear algebra
void rref( int, int, double *, double *);
void aw( double*, double, double* );

//void calc_CLC( MODEL &, double *, double *);
//void calc_CLC( MODEL & );
void i_dyn_CL( MODEL &, double *, double *);

void w2dQtn( double *, double*, double *);
// --- EOF ---
