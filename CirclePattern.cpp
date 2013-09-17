#include "CirclePattern.h"


double ClausensIntegral(double x)
{
	double sign = 1.0;
	if( x < 0.0 )
	{
		x=2*PI-fmod(-x,2*PI);
	}else if( x > 2*PI )
	{
		x=fmod(x,2*PI);
	}
	if( x > PI )
	{
		sign = -1.0;
		x = 2*PI-x;
	}
	if( x == 0.0 )
		return 0.0;

	const double coeff[7] = {1.0, -0.027777761035, -0.000277828908, -0.000004691463, -0.000000100652, -0.000000000773, -0.000000000104};
	double xsquare = x*x;
	double xpow = x;
	double claus = - x * log(2.0*sin(x/2.0));
	for(int i=0;i<7;i++)
	{
		claus += xpow * coeff[i];
		xpow *= xsquare;
	}
	return sign*claus;
}

double Lobachevsky(double x)
{
	return ClausensIntegral(2.0*x)/2.0;
}

double ImLi2plusinv(double x, double theta)
{
	// this is Im[ Li_2[ e^(x+i theta)]] + Im[ Li_2[ e^(-x+i theta)]]
	double p = 2.0 * atan(tanh(x/2.0)*tan((PI-theta)/2.0));
	return p * x + ClausensIntegral(p + PI - theta) + ClausensIntegral(-p+PI-theta)-ClausensIntegral(2.0*(PI-theta));
}

CirclePatternHessian::CirclePatternHessian(const Triangulation * const triangulation, const ThetaModel * const thetamodel, const std::vector<double> & logradius )
{
	matrix_rules_.resize(triangulation->NumberOfTriangles());
	for(int i=0;i<triangulation->NumberOfTriangles();i++)
	{
		Triangle * triangle = triangulation->getTriangle(i);
		for(int j=0;j<3;j++)
		{
			int nbr = triangle->getEdge(j)->getAdjacent()->getParent()->getId();
			double matrixElement = sin(thetamodel->getRealTheta(i,j))/(cosh(logradius[i]-logradius[nbr])-cos(thetamodel->getRealTheta(i,j)));
			std::pair<std::map<int,double>::iterator,bool> insertion;
			insertion = matrix_rules_[i].insert(std::pair<int,double>(nbr,matrixElement));
			if( !insertion.second )
			{
				insertion.first->second += matrixElement;
			}
			insertion = matrix_rules_[i].insert(std::pair<int,double>(i,-matrixElement));
			if( !insertion.second )
			{
				insertion.first->second -= matrixElement;
			}
		}
	}
}

void CirclePatternHessian::MultiplyVector(const std::vector<double> & from, std::vector<double> & to) const
{
	for(int i=0;i<static_cast<int>(from.size());i++)
	{
		to[i] = 0.0;
		for(std::map<int,double>::const_iterator it = matrix_rules_[i].begin(); it != matrix_rules_[i].end(); it++)
		{
			to[i] += it->second * from[it->first];
		}
	}
}

CirclePattern::CirclePattern(Triangulation * const triangulation, CohomologyBasis * const cohomologybasis, const ThetaModel * const thetamodel)
	: triangulation_(triangulation), thetamodel_(thetamodel), Embedding(triangulation,cohomologybasis)
{
	max_newton_iterations_ = 50;
}

bool CirclePattern::FindEdgeMeasure()
{
	if( !FindRadii() )
	{
		return false;
	}
	RadiiToAngles();

	for(int i=0, end=triangulation_->NumberOfTriangles();i<end;i++)
	{
		for(int j=0;j<3;j++)
		{
			Edge * adjEdge = triangulation_->getTriangle(i)->getEdge(j)->getAdjacent();
			double measure = 0.5 / tan( angles_[i][j] ) + 0.5 / tan( angles_[adjEdge->getParent()->getId()][adjEdge->getId()]);
			setEdgeMeasure(i,j,measure);
		}
	}

	return true;
}

bool CirclePattern::FindRadii()
{
	// given the theta, determine the radii of the circles in the circle pattern

	double alpha = 0.25;
	double beta = 0.3;

	std::vector<double> change(triangulation_->NumberOfTriangles(),0.0);
	std::vector<double> grad(triangulation_->NumberOfTriangles(),0.0);

	if( static_cast<int>(logradius_.size()) != triangulation_->NumberOfTriangles() )
	{
		logradius_.resize(triangulation_->NumberOfTriangles(),0.0);
	}

	int steps=0;
	std::vector<double> logradius2(triangulation_->NumberOfTriangles(),0.0);
	std::fill(logradius_.begin(),logradius_.end(),0.0);
	while(true)
	{
		steps++;
		EuclideanGradient(logradius_,grad);

		double maxgrad=0.0;
		for(int i=0;i<static_cast<int>(grad.size());i++)
		{
			if( fabs(grad[i]) > maxgrad )
				maxgrad = fabs(grad[i]);
		}
		if( maxgrad < 1e-7 )
			break;

		CirclePatternHessian hessian(triangulation_,thetamodel_,logradius_);
		hessian.ConjugateGradientSolve(grad,change,1e-8);

		double lambda=0.0;
		for(int i=0;i<static_cast<int>(grad.size());i++)
		{
			lambda += -grad[i]*change[i];
		}
		if( lambda*lambda/2.0 < 1e-8 )
		{
			break;
		}

		double t=1.0;
		double functional = EuclideanFunctional(logradius_);
		for(int i=0;i<static_cast<int>(logradius_.size());i++)
		{
			logradius2[i] = logradius_[i] + change[i];
		}
		double functional2 = EuclideanFunctional(logradius2);
		while( functional2 > functional - alpha * t * lambda )
		{
			t *= beta;
			for(int i=0;i<static_cast<int>(logradius2.size());i++)
			{
				logradius2[i] = logradius_[i] + t*change[i];
			}
			functional2 = EuclideanFunctional(logradius2);
		}
		for(int i=0;i<static_cast<int>(logradius_.size());i++)
		{
			logradius_[i] = logradius2[i];
		}

		if( steps > max_newton_iterations_ )
		{
			return false;
		}
	}

	std::cout << steps << "\n";

	double totlograd=0.0;
	for(int i=0;i<static_cast<int>(logradius_.size());i++)
	{
		totlograd += logradius_[i];
	}
	for(int i=0;i<static_cast<int>(logradius_.size());i++)
	{
		logradius_[i] -= totlograd/static_cast<int>(logradius_.size());
	}
	return true;
}

void CirclePattern::RadiiToAngles()
{
	if( static_cast<int>(angles_.size()) != triangulation_->NumberOfTriangles() )
		angles_.resize(triangulation_->NumberOfTriangles());

	for(int i=0;i<static_cast<int>(angles_.size());i++)
	{
		for(int j=0;j<3;j++)
		{
			int nbr = triangulation_->getTriangle(i)->getEdge(j)->getAdjacent()->getParent()->getId();
			angles_[i][j] = atan2(sin(thetamodel_->getRealTheta(i,j) ),exp(logradius_[i]-logradius_[nbr])-cos(thetamodel_->getRealTheta(i,j)));
		}
	}
}

void CirclePattern::EuclideanGradient(const std::vector<double> & rho, std::vector<double> & grad)
{
	if( static_cast<int>(grad.size()) != triangulation_->NumberOfTriangles() )
		grad.resize(triangulation_->NumberOfTriangles());
	for(int i=0;i<triangulation_->NumberOfTriangles();i++)
	{
		grad[i]=-PI;
		for(int j=0;j<3;j++)
		{
			int nbr = triangulation_->getTriangle(i)->getEdge(j)->getAdjacent()->getParent()->getId();
			double sine = sin(thetamodel_->getRealTheta(i,j) );
			double cosine = cos(thetamodel_->getRealTheta(i,j) );
			grad[i] += thetamodel_->getRealTheta(i,j) + atan2(sine,exp(rho[nbr]-rho[i])-cosine) - atan2(sine,exp(rho[i]-rho[nbr])-cosine);
		}
	}
}

double CirclePattern::EuclideanFunctional(const std::vector<double> & rho)
{
	double action=0.0;
	for(int i=0;i<static_cast<int>(rho.size());i++)
	{
		action += 2.0 * PI * rho[i];
		for(int j=0;j<3;j++)
		{
			int nbr = triangulation_->getTriangle(i)->getEdge(j)->getAdjacent()->getParent()->getId();
			action += 0.5*(ImLi2plusinv(rho[i]-rho[nbr],thetamodel_->getRealTheta(i,j)) - (PI - thetamodel_->getRealTheta(i,j))*(rho[i]+rho[nbr]));
		}
	}
	return action;
}