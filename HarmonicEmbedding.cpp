#include "HarmonicEmbedding.h"
#include <vector>
#include <algorithm>
#include <queue>
#include <cmath>

HarmonicEmbedding::HarmonicEmbedding(const Triangulation * const triangulation, CohomologyBasis * const cohomologybasis) 
	: triangulation_(triangulation), Embedding(triangulation,cohomologybasis), radiusdefinition_(SQUARE_ROOT_OF_AREA)
{

}

bool HarmonicEmbedding::FindEdgeMeasure()
{
	for(int i=0, end=triangulation_->NumberOfTriangles();i<end;i++)
	{
		for(int j=0;j<3;j++)
		{
			setEdgeMeasure(i,j,1.0);
		}
	}
	return true;
}

void HarmonicEmbedding::SetRadiusDefintion( RadiusDefinition radiusdefinition )
{
	radiusdefinition_ = radiusdefinition;
}

bool HarmonicEmbedding::GetRadii(std::vector<double> & radii)
{
	radii.clear();
	radii.reserve(triangulation_->NumberOfTriangles());

	std::pair<double,double> modulus = CalculateModuli();

	double totalarea = 0.0;
	for(int i=0,endi=triangulation_->NumberOfTriangles();i<endi;i++)
	{
		boost::array<Vector2D,3> side;
		for(int j=0;j<3;j++)
		{
			side[j] = ScaleVector2D(TransformByModulus(getForm(i,j),modulus),1.0/std::sqrt(modulus.second));
		}
		double radius;

		if( radiusdefinition_ == SQUARE_ROOT_OF_AREA )
		{
			double area = std::fabs(0.5 * ( - side[0][0] * side[1][1] + side[0][1] * side[1][0] ));
			double area2 = std::fabs(0.5 * ( - getForm(i,0)[0] * getForm(i,1)[1] + getForm(i,0)[1] * getForm(i,1)[0]) );
			BOOST_ASSERT( std::fabs(area - area2) < 1e-7 );
			totalarea += area;
			radius = std::sqrt(area);
		} else if( radiusdefinition_ == CIRCUMSCRIBED_CIRCLE )
		{
			double area = 0.5 * ( - side[0][0] * side[1][1] + side[0][1] * side[1][0] );
			radius = 0.25 * Norm2D(side[0]) * Norm2D(side[1]) * Norm2D(side[2]) / area;
		} else if( radiusdefinition_ == INSCRIBED_CIRCLE )
		{
			double a = Norm2D(side[0]), b = Norm2D(side[1]), c = Norm2D(side[2]);
			double s = 0.5*(a+b+c);
			radius = std::sqrt( (s-a) * (s-b) * (s-c) / s );
		} else if( radiusdefinition_ == AVERAGE_EDGE_LENGTH )
		{
			radius = (Norm2D(side[0]) + Norm2D(side[1]) + Norm2D(side[2]))/3.0;
		}
		radii.push_back( radius );
	}

	//BOOST_ASSERT( std::fabs(1.0 - totalarea) < 0.001 );

	return true;
}
