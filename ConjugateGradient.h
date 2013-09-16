#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H

#include <vector>

extern "C"
void dscal_( const int *n, const double *alpha, double *x, const int *incx );

extern "C"
void dcopy_( const int *n, const double *x, const int *incx, double *y, const int *incy );

extern "C"
void daxpy_( const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy );

extern "C"
double ddot_( const int *n, const double *x, const int *incx, const double *y, const int *incy );

inline
void dscal( double alpha, std::vector<double> & x) {
	int n = static_cast<int>(x.size());
	int incx=1;
	dscal_(&n,&alpha,&x.front(),&incx);
}

inline
void dcopy( const std::vector<double> & x, std::vector<double> & y ) {
	int n = static_cast<int>(x.size());
	int incx = 1, incy = 1;
	dcopy_(&n,&x.front(),&incx,&y.front(),&incy);
}

inline
void daxpy(  double alpha, const std::vector<double> & x, std::vector<double> & y ) {
  	int n = static_cast<int>(x.size());
	int incx = 1, incy = 1;
	daxpy_(&n,&alpha,&x.front(),&incx,&y.front(),&incy);
}

inline
double ddot( const std::vector<double> & x, const std::vector<double> & y ) {
  	int n = static_cast<int>(x.size());
	int incx = 1, incy = 1;
	return ddot_(&n,&x.front(),&incx,&y.front(),&incy);
}

class Matrix  
{
public:
	Matrix() {}
	~Matrix() {}

	virtual void MultiplyVector(const std::vector<double> & from, std::vector<double> & to) const = 0;

	bool ConjugateGradientSolve( const std::vector<double> & b, std::vector<double> & x, double eps, int maxIterations = 100 ) const;

};


#endif 
