#include <cmath>
#include <fstream>
#include <complex>

#include "boost/random/mersenne_twister.hpp"
#include "boost/random/uniform_real.hpp"


#include "LinearAlgebra.h"
#include "utilities.h"

#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include "lapacke.h"

void linearalgebra::Scale( double factor, std::vector<double> & x )
{
	int n = static_cast<int>(x.size());
	int incx=1;
	dscal_(&n,&factor,&x.front(),&incx);
}

void linearalgebra::Copy( const std::vector<double> & from, std::vector<double> & to )
{
	int n = static_cast<int>(from.size());
	int incx = 1, incy = 1;
	dcopy_(&n,&from.front(),&incx,&to.front(),&incy);
}

void linearalgebra::AddTo( double factor, const std::vector<double> & v, std::vector<double> & to )
{
  	int n = static_cast<int>(v.size());
	int incx = 1, incy = 1;
	daxpy_(&n,&factor,&v.front(),&incx,&to.front(),&incy);
}

void linearalgebra::Swap( std::vector<double> & v1, std::vector<double> & v2 )
{
 	int n = static_cast<int>(v1.size());
	int incx = 1, incy = 1;
	dswap_(&n,&v1.front(),&incx,&v2.front(),&incy);
}

double linearalgebra::DotProduct( const std::vector<double> & v1, const std::vector<double> & v2 )
{
	int n = static_cast<int>(v1.size());
	int incx = 1, incy = 1;
	return ddot_(&n,&v1.front(),&incx,&v2.front(),&incy);
}

double linearalgebra::NormSquared( const std::vector<double> & v )
{
	int n = static_cast<int>(v.size());
	int incx = 1,incy=1;
	return ddot_(&n,&v.front(),&incx,&v.front(),&incy);
}

double linearalgebra::Total( const std::vector<double> & v )
{
	double total = 0.0;
	for(std::vector<double>::const_iterator it = v.begin(); it != v.end(); ++it )
	{
		total += *it;
	}
	return total;
}

bool linearalgebra::Matrix::ConjugateGradientSolve( const std::vector<double> & b, std::vector<double> & x, double eps, int maxIterations ) const
{
	std::vector<double> g(b.size());
	std::vector<double> r(b.size());
	std::vector<double> p(b.size());
	double t, tau, sig, rho, gam;
	double err = eps * eps * linearalgebra::NormSquared(b);

	MultiplyVector(x,g);
	linearalgebra::AddTo(-1.0,b,g);
	linearalgebra::Scale(-1.0,g);
	linearalgebra::Copy(g,r);
	int iterations = 0;
	while ( linearalgebra::NormSquared(g) > err ) 
	{
		MultiplyVector(r,p);
		rho = linearalgebra::NormSquared(p);
		sig = linearalgebra::DotProduct(r,p);
		tau = linearalgebra::DotProduct(g,r);
		t = tau/sig;
		linearalgebra::AddTo(t,r,x);

		linearalgebra::AddTo(-t,p,g);
		gam = (t*t*rho-tau)/tau;
		linearalgebra::Scale(gam,r);
		linearalgebra::AddTo(1.0,g,r);
		iterations++;
		if( iterations > maxIterations )
			return false;
	}
	return true;
}

void linearalgebra::Matrix::ComputeEigenvalues()
{
	eigenvalues_.clear();
	eigenvalues_.reserve(dimension_);
	eigenvectors_.resize(dimension_,std::vector<double>(dimension_,0.0));
	for(int i=0;i<dimension_;i++)
	{
		ComputeNextEigenvalue();
	}
}

void linearalgebra::Matrix::ComputeNextEigenvalue()
{
	int MaxIterations = 200;
	double accuracy = 1e-7;
	double MinNorm = 1e-10;

	static boost::mt19937 rng;

	int index = eigenvalues_.size();
	double invsqrtsize = 1.0/std::sqrt((double)dimension_);
	boost::uniform_real<> distribution(0.0,1.0);
	for(int i=0;i<dimension_;i++)
	{
		eigenvectors_[index][i] = invsqrtsize*distribution(rng);
	}
	for(int i=0;i<index;i++)
	{
		ProjectOut(eigenvectors_[i],eigenvectors_[index]);
	}

	int iterations = 0;
	std::vector<double> result(dimension_);
	eigenvalues_.push_back(0.0);
	while( iterations < MaxIterations )
	{
		MultiplyVector(eigenvectors_[index],result);
		for(int i=0;i<index;i++)
		{
			ProjectOut(eigenvectors_[i],result);
		}

		eigenvalues_[index] = linearalgebra::DotProduct(result,eigenvectors_[index]);
		double normresult = std::sqrt( linearalgebra::NormSquared(result) );

		if( normresult < MinNorm )
		{
			// eigenvalue is roughly zero
			break;
		} else
		{
			linearalgebra::Scale( 1.0/normresult, result );
			linearalgebra::AddTo( -1.0, result, eigenvectors_[index] );
			double difference = linearalgebra::NormSquared( eigenvectors_[index] );
			linearalgebra::Copy( result, eigenvectors_[index] );

			if( difference < accuracy*accuracy )
			{
				break;
			}
		}
		iterations++;
	}
}

void linearalgebra::Matrix::ProjectOut(const std::vector<double> & projector, std::vector<double> & v) const
{
	double factor = linearalgebra::DotProduct(projector,v) / linearalgebra::NormSquared(projector);
	linearalgebra::AddTo(-factor,projector,v);
}

double linearalgebra::Matrix::GetEigenvalue(int i) const
{
	return eigenvalues_[i];
}

double linearalgebra::Matrix::GetEigenvector(int i, int vertex) const
{
	return eigenvectors_[i][vertex];
}

linearalgebra::DenseMatrix::DenseMatrix(int n) :
	Matrix(n),
	size_(n)
{
	matrix_.resize(n*n,0.0);
}


void linearalgebra::DenseMatrix::MultiplyVector(const std::vector<double> & from, std::vector<double> & to) const
{
	for(int i=0;i<size_;i++)
	{
		to[i] = 0.0;
		for(int j=0;j<size_;j++)
		{
			to[i] += matrix_[size_*i + j] * from[j];
		}
	}

}

void linearalgebra::DenseMatrix::Set(int x, int y, double value)
{
	matrix_[size_*x + y] = value;
}

double linearalgebra::DenseMatrix::Get(int x, int y) const
{
	return (x>y?matrix_[size_*x+y]:matrix_[size_*y+x]);
}

linearalgebra::PositiveDefiniteDenseMatrix::PositiveDefiniteDenseMatrix(int n) :
	DenseMatrix(n)
{
	ipiv_.resize(n);
	workspace_.resize(n*n);
}



bool linearalgebra::PositiveDefiniteDenseMatrix::ComputeInverse()
{
	std::cout << "LU - ";
	int info=LAPACKE_dpotrf(LAPACK_COL_MAJOR,'U',size_,&matrix_[0],size_);
	if( info != 0)
	{
		return false;
	}
	std::cout << "invert - ";
	info = LAPACKE_dpotri(LAPACK_COL_MAJOR,'U',size_,&matrix_[0],size_);
	std::cout << "done.\n";
	return info == 0;
}