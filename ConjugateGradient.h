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
	int n = x.size();
	int incx=1;
	dscal_(&n,&alpha,&x.front(),&incx);
}

inline
void dcopy( const std::vector<double> & x, std::vector<double> & y ) {
	int n = x.size();
	int incx = 1, incy = 1;
	dcopy_(&n,&x.front(),&incx,&y.front(),&incy);
}

inline
void daxpy(  double alpha, const std::vector<double> & x, std::vector<double> & y ) {
  	int n = x.size();
	int incx = 1, incy = 1;
	daxpy_(&n,&alpha,&x.front(),&incx,&y.front(),&incy);
}

inline
double ddot( const std::vector<double> & x, const std::vector<double> & y ) {
  	int n = x.size();
	int incx = 1, incy = 1;
	return ddot_(&n,&x.front(),&incx,&y.front(),&incy);
}

class Matrix  
{
public:
	Matrix() {}
	~Matrix() {}

	virtual void MultiplyVector(const std::vector<double> & from, std::vector<double> & to) const = 0;

	bool ConjugateGradientSolve( const std::vector<double> & b, std::vector<double> & x, double eps, int maxIterations = 100 ) const
	{
		std::vector<double> g(b.size());
		std::vector<double> r(b.size());
		std::vector<double> p(b.size());
		double t, tau, sig, rho, gam;
		double err=eps*eps*ddot(b,b);

		MultiplyVector(x,g);
		daxpy(-1.0,b,g);
		dscal(-1.0,g);
		dcopy(g,r);
		int iterations = 0;
		double error;
		while ( (error=ddot(g,g)) > err ) 
		{
			double total = 0.0;
			for(int i=0;i<x.size();i++)
			{
				total += x[i];
			}

			MultiplyVector(r,p);
			rho = ddot(p,p);
			sig=ddot(r,p);
			tau=ddot(g,r);
			t=tau/sig;
			daxpy(t,r,x);

			daxpy(-t,p,g);
			gam=(t*t*rho-tau)/tau;
			dscal(gam,r);
			daxpy(1.0,g,r);
			iterations++;
			if( iterations > maxIterations )
				return false;
		}
		return true;
	}
};


#endif 
