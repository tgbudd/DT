#ifndef DT_LINEAR_ALGEBRA_H
#define DT_LINEAR_ALGEBRA_H

#include <vector>

extern "C"
void dscal_( const int *n, const double *alpha, double *x, const int *incx );

extern "C"
void dcopy_( const int *n, const double *x, const int *incx, double *y, const int *incy );

extern "C"
void dswap_( const int *n, const double *x, const int *incx, double *y, const int *incy );

extern "C"
void daxpy_( const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy );

extern "C"
double ddot_( const int *n, const double *x, const int *incx, const double *y, const int *incy );

extern "C"
double dnrm2_( const int *n, const double *x, const int *incx );

namespace linearalgebra {

	void Scale( double factor, std::vector<double> & vector );
	void Copy( const std::vector<double> & from, std::vector<double> & to );
	void AddTo( double factor, const std::vector<double> & v, std::vector<double> & to );
	void Swap( std::vector<double> & v1, std::vector<double> & v2 );

	double DotProduct( const std::vector<double> & v1, const std::vector<double> & v2 );
	double NormSquared( const std::vector<double> & v );

	double Total( const std::vector<double> & v );

	
	class Matrix  
	{
	public:
		Matrix(int dimension) : dimension_(dimension) {}
		~Matrix() {}

		virtual void MultiplyVector(const std::vector<double> & from, std::vector<double> & to) const = 0;

		bool ConjugateGradientSolve( const std::vector<double> & b, std::vector<double> & x, double eps, int maxIterations = 100 ) const;
	
		void ComputeEigenvalues();
		double GetEigenvalue(int i) const;
		double GetEigenvector(int i, int vertex) const;
	private:
		void ComputeNextEigenvalue();
		void ProjectOut(const std::vector<double> & projector, std::vector<double> & v) const;

		const int dimension_;
		std::vector<double> eigenvalues_;
		std::vector<std::vector<double> > eigenvectors_;
	};

	class DenseMatrix : public Matrix
	{
	public:
		DenseMatrix(int n);
		void MultiplyVector(const std::vector<double> & from, std::vector<double> & to) const;
		void Set(int x, int y, double value);
		double Get(int x, int y) const;
	protected:
		const int size_;
		std::vector<double> matrix_;
	};

	class PositiveDefiniteDenseMatrix : public DenseMatrix
	{
	public:
		PositiveDefiniteDenseMatrix(int n);
		bool ComputeInverse();
	private:
		std::vector<int> ipiv_;
		std::vector<double> workspace_;
	};
}

#endif