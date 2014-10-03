#include <queue>

#include "CirclePacking.h"
#include "TriangulationProperties.h"
#include "HarmonicEmbedding.h"
#include "utilities.h"

CirclePacking::CirclePacking(Triangulation * const triangulation, CohomologyBasis * const cohomologybasis)
	: triangulation_(triangulation), Embedding(triangulation,cohomologybasis), babyuniversedetector_(triangulation)
	, cohomologybasis_(cohomologybasis)
{
	max_iterations_ = 5000;
	epsilon_ = 1.0e-5;
	delta_ = 0.01;
}

double CirclePacking::Angle(Edge * edge)
{
	double radius = radius_[edge->getOpposite()->getId()];
	double nbrRadius1 = radius_[edge->getPrevious()->getOpposite()->getId()];
	double nbrRadius2 = radius_[edge->getNext()->getOpposite()->getId()];
	double x = nbrRadius1 / (nbrRadius1 + radius) * nbrRadius2 / (nbrRadius2 + radius);
	return 2.0 * std::asin(std::sqrt( x ) );
}

double CirclePacking::AngleSum(int vertexId)
{
	Vertex * vertex = triangulation_->getVertex(vertexId);
	Edge * edge = vertex->getParent();
	double angle = 0.0;
	do {
		angle += Angle(edge);
		edge = edge->getNext()->getAdjacent()->getNext();
	} while( edge != vertex->getParent() );
	return angle;
}

bool CirclePacking::FindRadii()
{
	if( static_cast<int>(radius_.size()) != triangulation_->NumberOfVertices() )
	{
		radius_.resize(triangulation_->NumberOfVertices());
	}

	double startRadius = std::sqrt(1.0/triangulation_->NumberOfVertices());
	std::fill(radius_.begin(),radius_.end(),startRadius);
	std::vector<double> radius2(triangulation_->NumberOfVertices(),startRadius);

	properties::DegreeList(triangulation_,degree_);

	double c = epsilon_ + 1;
	double lambda = -1.0;
	bool flag = false;
	int steps = 0;
	while(c > epsilon_ )
	{
		double c0 = c;
		c = 0.0;
		double lambda0 = lambda;
		bool flag0 = flag;

		//double totradiussq = 0.0;
		for(int i=0,endi=radius_.size();i<endi;i++)
		{
			double theta = AngleSum(i);
			double sintheta = std::sin(theta/(2*degree_[i]));
			double sintarget = std::sin(PI/degree_[i]);
			radius2[i] = (1.0 - sintarget)/sintarget * sintheta / (1.0 - sintheta) * radius_[i];
			//totradiussq += radius2[i]*radius2[i];
			c += (theta - 2.0 * PI)*(theta - 2.0 * PI);
		}
		//double normalize = 1.0/std::sqrt(totradiussq);
		//for(int i=0,endi=radius_.size();i<endi;i++)
		//{
		//	radius2[i] *= normalize;
		//}
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
		} else
		{
			std::copy(radius2.begin(),radius2.end(),radius_.begin());
		}
		steps++;
		if( steps > max_iterations_ )
		{
			return false;
		}
	}

	double area=0.0;
	for(int i=0,endi=radius_.size();i<endi;i++)
	{
		area += PI*radius_[i]*radius_[i];
	}
	double norm = 1.0/std::sqrt(area);
	for(int i=0,endi=radius_.size();i<endi;i++)
	{
		radius_[i] *= norm;
	}
	return true;
}

bool CirclePacking::FindEdgeMeasure()
{
	if( !FindRadii() )
	{
		return false;
	}

	for(int i=0,end=triangulation_->NumberOfTriangles();i<end;i++)
	{
		for(int j=0;j<3;j++)
		{
			Edge * edge = triangulation_->getTriangle(i)->getEdge(j);
			Edge * adjEdge = edge->getAdjacent();
			double measure = 0.5 / std::tan( Angle(edge) ) + 0.5 / std::tan( Angle(adjEdge) );
			setEdgeMeasure(i,j,measure);
		}
	}

	return true;
}

void CirclePacking::CalculateEdgeMeasure(std::vector<boost::array<double,3> > & measure)
{
	if( static_cast<int>(measure.size()) !=  triangulation_->NumberOfTriangles() )
	{
		measure.resize(triangulation_->NumberOfTriangles());
	}
	for(int i=0,end=triangulation_->NumberOfTriangles();i<end;i++)
	{
		for(int j=0;j<3;j++)
		{
			Edge * edge = triangulation_->getTriangle(i)->getEdge(j);
			Edge * adjEdge = edge->getAdjacent();
			double angle1 = Angle(edge), angle2 = Angle(adjEdge);
			if( angle1 > 1.0e-7 && angle2 > 1.0e-7 )
			{
				measure[i][j] = 0.5 / std::tan( Angle(edge) ) + 0.5 / std::tan( Angle(adjEdge) );
			} else
			{
				measure[i][j] = 1.0;
			}
		}
	}
}

bool CirclePacking::GetRadii(std::vector<double> & radii)
{
	radii.clear();
	radii.reserve(triangulation_->NumberOfVertices());

	radii = radius_;

	return true;
}
void CirclePacking::GetAbsoluteRadii(std::vector<double>& radii)
{
	// return the radii of the circles when they are embedded in a torus of area 1
	std::pair<double,double> modulus = CalculateModuli();
	
	const Edge * edge = triangulation_->getTriangle(0)->getEdge(0);
	double edgeLength = std::sqrt(NormSquaredTransformedByModulus(getForm(edge),modulus));
	double fraction = edgeLength / (radius_[edge->getEnd()->getId()]+radius_[edge->getStart()->getId()]);
	
	radii.resize(radius_.size());
	std::transform(radius_.begin(),radius_.end(),radii.begin(),std::bind1st(std::multiplies<double>(),fraction));
}
