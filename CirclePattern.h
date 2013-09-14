#pragma once

#include <map>

#include "boost/array.hpp"

#include "embedding.h"
#include "ThetaModel.h"
#include "CohomologyBasis.h"
#include "ConjugateGradient.h"

double ClausensIntegral(double x);
double Lobachevsky(double x);
double ImLi2plusinv(double x, double theta);

class CirclePatternHessian : public Matrix
{
public:
	CirclePatternHessian(const Triangulation * const triangulation, const ThetaModel * const thetamodel, const std::vector<double> & logradius); 
	void MultiplyVector(const std::vector<double> & from, std::vector<double> & to) const;
private:
	std::vector<std::map<int,double> > matrix_rules_;
};

class CirclePattern :
	public Embedding
{
public:
	CirclePattern(Triangulation * const triangulation, const CohomologyBasis * const cohomologybasis, const ThetaModel * const thetamodel);
	~CirclePattern() {}

	bool FindEdgeMeasure();

private:
	void EuclideanGradient(const std::vector<double> & rho, std::vector<double> & grad);
	double EuclideanFunctional(const std::vector<double> & rho);
	bool FindRadii();
	void RadiiToAngles();

	const Triangulation * const triangulation_;
	const ThetaModel * const thetamodel_;

	std::vector<double> logradius_;
	std::vector<boost::array<double,3> > angles_;

	int max_newton_iterations_;
};

