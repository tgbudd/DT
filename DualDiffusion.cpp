#include <queue>

#include "utilities.h"
#include "Triangle.h"
#include "Edge.h"
#include "DualDiffusion.h"

void DualDiffusionMatrix::MultiplyVector(const std::vector<double> & from, std::vector<double> & to) const
{
	for(int i=0;i<triangulation_->NumberOfTriangles();i++)
	{
		Triangle * triangle = triangulation_->getTriangle(i);
		to[i] = (1.0-walkprobability_)*from[i];
		for(int j=0;j<3;j++)
		{
			to[i] += walkprobability_ /3.0 * from[triangle->getEdge(j)->getAdjacent()->getParent()->getId()];
		}
	}
}


DualDiffusion::DualDiffusion( Triangulation * const triangulation) : triangulation_(triangulation),
	samples_(20), measurements_(0), diffusion_steps_(500), walkprobability_(0.8), max_distance_(50)
{
	distance_.resize( triangulation_->NumberOfTriangles(),0 );
	probability_.resize(triangulation_->NumberOfTriangles(),0.0 );
	next_probability_.resize(triangulation_->NumberOfTriangles(),0.0 );
	distribution_.resize(diffusion_steps_,std::vector<double>(max_distance_+1,0.0));
}


DualDiffusion::~DualDiffusion(void)
{
}

void DualDiffusion::Measure()
{
	DualDiffusionMatrix diffusionmatrix( triangulation_, walkprobability_);

	for(int i=0;i<samples_;i++)
	{
		Triangle * startTriangle = triangulation_->getRandomTriangle();
		DetermineDistance( startTriangle );
		DoDiffusion( startTriangle, &diffusionmatrix );
		measurements_++;
	}
}

void DualDiffusion::DetermineDistance(Triangle * startTriangle)
{
	std::fill(distance_.begin(),distance_.end(),-1);

	std::queue<Triangle *> q;
	q.push( startTriangle );
	distance_[startTriangle->getId()] = 0;

	while( !q.empty() )
	{
		Triangle * triangle = q.front();
		q.pop();

		for(int i=0;i<3;i++)
		{
			Triangle * nbr = triangle->getEdge(i)->getAdjacent()->getParent();
			if( distance_[nbr->getId()] == -1 )
			{
				distance_[nbr->getId()] = distance_[triangle->getId()] + 1;
				q.push( nbr );
			}
		}
	}
}

void DualDiffusion::DoDiffusion(Triangle * startTriangle, const linearalgebra::Matrix * const matrix )
{
	std::fill(probability_.begin(),probability_.end(),0.0);
	probability_[startTriangle->getId()] = 1.0;

	for(int i=0;i<diffusion_steps_;i++)
	{
		matrix->MultiplyVector(probability_,next_probability_);
		linearalgebra::Copy(next_probability_,probability_);
		DoMeasurementOnDistribution(i+1);
	}
}

void DualDiffusion::DoMeasurementOnDistribution(int time)
{
	for(int j=0,end=distance_.size();j<end;j++)
	{
		if( distance_[j] <= max_distance_ )
			distribution_[time-1][distance_[j]] += probability_[j];
	}
}

std::string DualDiffusion::OutputData() const
{
	std::ostringstream stream;
	stream << std::fixed << "dualdiffusion -> {measurements -> " << measurements_ << ", walkprobability -> " << walkprobability_;
	stream << ", diffusionsteps -> " << diffusion_steps_ << ", samples -> " << samples_;
	stream << ", distribution -> ";
	PrintToStream2D(stream,distribution_.begin(),distribution_.end());
	stream << "/" << measurements_ << "}";
	return stream.str();
}