#include "CirclePacking.h"
#include "TriangulationProperties.h"

CirclePacking::CirclePacking(Triangulation * const triangulation, CohomologyBasis * const cohomologybasis)
	: triangulation_(triangulation), Embedding(triangulation,cohomologybasis)
{
	max_iterations_ = 50;
	epsilon_ = 1.0e-7;
	delta_ = 0.01;
}

double CirclePacking::AngleSum(int vertexId)
{
	Vertex * vertex = triangulation_->getVertex(vertexId);
	Edge * edge = vertex->getParent()->getPrevious();
	double angle = 0.0;
	do {
		double nbrRadius1 = radius_[edge->getPrevious()->getOpposite()->getId()];
		double nbrRadius2 = radius_[edge->getOpposite()->getId()];
		angle += 2.0 * std::asin(std::sqrt( nbrRadius1 / (nbrRadius1 + radius_[vertexId]) * nbrRadius2 / (nbrRadius2 + radius_[vertexId]) ) );
		edge = edge->getPrevious()->getAdjacent();
	} while( edge->getNext() != vertex->getParent() );
	return angle;
}

bool CirclePacking::FindRadii()
{
	// given the theta, determine the radii of the circles in the circle pattern

	if( static_cast<int>(radius_.size()) != triangulation_->NumberOfVertices() )
	{
		radius_.resize(triangulation_->NumberOfVertices());
	}
	std::fill(radius_.begin(),radius_.end(),1.0);
	std::vector<double> radius2(triangulation_->NumberOfVertices(),1.0);

	properties::DegreeList(triangulation_,degree_);

	double c = epsilon_ + 1;
	double lambda = -1.0;
	bool flag = false;
	int steps = 0;
	while(c > epsilon_ )
	{
		double c0 = c;
		double lambda0 = lambda;
		bool flag0 = flag;

		for(int i=0,endi=radius_.size();i<endi;i++)
		{
			double theta = AngleSum(i);
			double sintheta = std::sin(theta/(2*degree_[i]));
			double sintarget = std::sin(PI/degree_[i]);
			radius2[i] = (1.0 - sintarget)/sintarget * sintheta / (1.0 - sintheta) * radius_[i];
			c += (theta - 2.0 * PI)*(theta - 2.0 * PI);
		}

		c = std::sqrt(c);
		lambda = c / c0;
		flag = true;
		if( flag0 && lambda < 1.0 )
		{
			c *= lambda;
			if( std::abs( lambda - lambda0 ) < delta_ )
			{
				lambda = lambda / (1.0 - lambda);
			}
			double maxlambda = 1000.0;
			for(int i=0,endi=radius_.size();i<endi;i++)
			{
				if( radius2[i] + maxlambda*(radius2[i]-radius_[i]) < 0.0 )
				{
					maxlambda = - radius2[i] / (radius2[i]-radius_[i]);
				}
			}
			lambda = std::min(lambda, 0.5*maxlambda );
			for(int i=0,endi=radius_.size();i<endi;i++)
			{
				radius_[i] = radius2[i] + lambda * (radius2[i]-radius_[i]);
			}
			flag = false;
		}
		steps++;
		if( steps > max_iterations_ )
		{
			return false;
		}
	}

	return true;
}