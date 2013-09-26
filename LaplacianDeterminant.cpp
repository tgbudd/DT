#include "LaplacianDeterminant.h"
#include <cmath>

void DualLaplacianMatrix::MultiplyVector(const std::vector<double> & from, std::vector<double> & to) const
{
	for(int i=0;i<triangulation_->NumberOfTriangles();i++)
	{
		Triangle * triangle = triangulation_->getTriangle(i);
		to[i] = 3.0*from[i];
		for(int j=0;j<3;j++)
		{
			to[i] -= from[triangle->getEdge(j)->getAdjacent()->getParent()->getId()];
		}
	}
}

double LaplacianDeterminant::BoltzmannChangeUnderFlipMove(const Edge * const edge) const
{
	boost::array<Triangle *,4> triangles;
	triangles[0] = edge->getParent();
	triangles[1] = edge->getAdjacent()->getParent();
	triangles[2] = edge->getPrevious()->getAdjacent()->getParent();
	triangles[3] = edge->getAdjacent()->getPrevious()->getAdjacent()->getParent();

	boost::array<boost::array<double,4>,2> y;

	if( triangles[2] == triangles[3] )
	{
		return 1.0;
	}

	for(int k=0;k<2;k++)
	{
		righthandside_[triangles[2*(1-k)]->getId()] = 1.0;
		righthandside_[triangles[2*(1-k)+1]->getId()] -= 1.0;

		std::fill(solution_.begin(),solution_.end(),0.0);

		duallaplacianmatrix_.ConjugateGradientSolve(righthandside_,solution_,1.0e-4,100);

		righthandside_[triangles[2*(1-k)]->getId()] = 0.0;
		righthandside_[triangles[2*(1-k)+1]->getId()] = 0.0;
	
		double total=0.0;
		for(int i=0;i<static_cast<int>(solution_.size());i++)
		{
			total+=solution_[i];
		}
		total = total / static_cast<int>(solution_.size());
		for(int i=0;i<4;i++)
		{
			y[k][i] = solution_[triangles[i]->getId()] - total;
		}
	}
	double ddet = 1 + y[0][0] - y[0][1] + y[1][2] - y[1][3] - y[0][2]*y[1][0] + y[0][3]*y[1][0]
					+ y[0][2]*y[1][1] - y[0][3]*y[1][1] + y[0][0]*y[1][2] 
					- y[0][1]*y[1][2] - y[0][0]*y[1][3] + y[0][1]*y[1][3];
	
	return pow(ddet,-0.5*centralcharge_);
}
