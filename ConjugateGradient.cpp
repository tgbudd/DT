#include "ConjugateGradient.h"

bool Matrix::ConjugateGradientSolve( const std::vector<double> & b, std::vector<double> & x, double eps, int maxIterations ) const
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
		for(int i=0;i<static_cast<int>(x.size());i++)
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