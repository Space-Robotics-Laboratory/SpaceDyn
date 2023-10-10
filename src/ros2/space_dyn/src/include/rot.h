
//_/_/_/_/ header file /_/_/_/_//

//_SD_//
void aw ( double *, double *);

void cx ( double *, double *);
void cy ( double *, double *);
void cz ( double *, double *);

void rpy2dc ( double *, double *);
void dc2rpy ( double *, double *);
void dc2qtn ( double *, double *);
void qtn2dc ( double *, double *);
void w2dqtn ( double *, double *, double *);

//_Sato_you need to modify again//
void rpy2R ( double *, double * );
void R2rpy ( double *, double * );
void R2qtn ( double *, double * );
void rpy2qtn ( double *, double * );

void tilde ( int, double *, double *);

double pi_conv( char *data );
double pi2deg( double pi );
double deg2pi( double deg );

void QMultiply( double *, double *, double * );
