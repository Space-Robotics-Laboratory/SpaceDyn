//_/_/_/_/ vector library _/_/_/_//

double *vector_get( int );
void vector_cpy( int, double*, double* );
void vector_print( int, double* );
void vector_print_int( int, int* );
void vector_add( int, double*, double*, double* );
void vector_sub( int, double*, double*, double* );
void vector_cross3( double*, double*, double* );
void vector_tilde3( double*, double* );
void vector_z( int, double* );
void vector_i( int, double* );
void vector_scale( int, double, double*, double* );
double vector_inner( int, double*, double* );
