#ifndef CBLAS1_H
// ============================================================================

/*
// x <= a*x
void dscal( int n, double alpha, double *x, int incx );


// y <= x
void dcopy( int n, const double *x, int incx, double *y, int incy );


// y <= a*x+y
void daxpy( int n, double alpha, const double *x, int incx, double *y,
	    int incy );


// dot <= x^T*y
double ddot( int n, const double *x, int incx, const double *y, int incy );
*/

// x <= a*x
extern "C"
void dscal_( const int *n, const double *alpha, double *x, const int *incx );


// y <= x
extern "C"
void dcopy_( const int *n, const double *x, const int *incx, double *y,
	     const int *incy );


// y <= a*x+y
extern "C"
void daxpy_( const int *n, const double *alpha, const double *x,
	     const int *incx, double *y, const int *incy );


// dot <= x^T*y
extern "C"
double ddot_( const int *n, const double *x, const int *incx, const double *y,
	      const int *incy );



inline
void dscal( int n, double alpha, double *x, int incx ) {
  dscal_(&n,&alpha,x,&incx);
}

inline
void dcopy( int n, const double *x, int incx, double *y, int incy ) {
  dcopy_(&n,x,&incx,y,&incy);
}

inline
void daxpy( int n, double alpha, const double *x, int incx, double *y,
	    int incy ) {
  daxpy_(&n,&alpha,x,&incx,y,&incy);
}

inline
double ddot( int n, const double *x, int incx, const double *y, int incy ) {
  return ddot_(&n,x,&incx,y,&incy);
}
