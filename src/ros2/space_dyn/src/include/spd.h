#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include "../include/model.h"
#include "../include/define.h"
#include "../include/common.h"

void model(string, MODEL &);
void model_test( string );
void model_CL(string, MODEL &);
void model_new(string, MODEL &);
void model_init(MODEL &);
void model_converter( string, string ); // Rainer def to Satoko def.

void calc_SPN(MODEL &);
void ch_conf(MODEL &, int *, MODEL &);
void dat_consistency(MODEL &, int *, MODEL &, int );

void calc_SP (MODEL &);
void calc_hh (MODEL &, double *);
void calc_Rg (MODEL &, double *);
void calc_JJ (MODEL &, double **);
void calc_Je (MODEL &, int, double * );
void calc_GJb (MODEL &, double * );
void calc_GJe (MODEL &, int, double * );
void calc_GJe_r (MODEL &, int, double * );
void calc_Jb (MODEL &, int, double * );
void calc_Lg( MODEL &, double *);
void f_kin_e (MODEL &, int );
void f_kin_j (MODEL &);
void i_dyn (MODEL &, double *, double *);
void i_dyn_fix(MODEL &, double *, double *);
void i_dyn_t (MODEL &, double *, double *);

void f_dyn (MODEL &, double *, double *, double *, double *);
void f_dyn_ffg(MODEL &, double *, double *, double *, double *);
void f_dyn_fix(MODEL &, double *, double *);
//void f_dyn (MODEL &, double * );
void inner_force (MODEL &, double *, double **);

void calc_C(MODEL &, double *, double * );
//void f_dyn_rk (MODEL &m, double *, double, double *, double *, double *, double *, double *, double *, double *); // fixed step runge-kutta
void f_dyn_rk(MODEL &m, double *Gravity, double step); // fixed step runge-kutta
void f_dyn_rk1(MODEL &m, double *Gravity, double step); // fixed step runge-kutta
void f_dyn_rk_ffg(MODEL &m, double *Gravity, double step); // fixed step runge-kutta
void int_eu( MODEL &m, double step );

void calc_Xup_I(MODEL &);
void calc_SS (MODEL &, double *);
void j_num(MODEL &m, int, int *, int, int );
//void aw(double *, double, double *);

// simple explicit euler-method integrator
void int_eu(MODEL &, double *, double *, double *, double *, double *, double *);

// linear algebra
void rref( int, int, double *, double *);
void aw( double*, double, double* );

//void calc_CLC( MODEL &, double *, double *);
//void calc_CLC( MODEL & );
void calc_CLC2( MODEL & );
void i_dyn_CL( MODEL &, double *, double *);

void w2dQtn( double *, double*, double *);



// --- EOF ---
