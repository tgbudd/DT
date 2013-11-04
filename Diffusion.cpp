#include <queue>

#include "utilities.h"
#include "Triangle.h"
#include "Edge.h"
#include "Diffusion.h"
#include "TriangulationProperties.h"

DiffusionMatrix::DiffusionMatrix(const Triangulation * const triangulation, const std::vector<int> & degree)
	: Matrix(triangulation->NumberOfVertices())
{
	laplacianRules_.resize(triangulation->NumberOfVertices());
	for(int i=0;i<triangulation->NumberOfVertices();i++)
	{
		Vertex * vertex = triangulation->getVertex(i);
		
		std::pair<std::map<int,double>::iterator,bool> insertion;
		Edge * edge = vertex->getParent()->getPrevious();
		do {
			int nbrId = edge->getPrevious()->getOpposite()->getId();

			insertion = laplacianRules_[i].insert(std::pair<int,double>(nbrId,1.0/degree[nbrId]));
			if( !insertion.second )
			{
				insertion.first->second += 1.0/degree[nbrId];
			}
			edge = edge->getPrevious()->getAdjacent();
		} while( edge->getNext() != vertex->getParent() );
	}
}

void DiffusionMatrix::MultiplyVector(const std::vector<double> & from, std::vector<double> & to) const
{
	for(int i=0;i<static_cast<int>(from.size());i++)
	{
		to[i] = 0.0;
		for(std::map<int,double>::const_iterator it = laplacianRules_[i].begin(); it != laplacianRules_[i].end(); it++)
		{
			to[i] += it->second * from[it->first];
		}
	}
}


Diffusion::Diffusion( Triangulation * const triangulation ) : triangulation_(triangulation),
	samples_(20), measurements_(0), max_distance_(50), distance_only_(false)
{
	int size = triangulation_->NumberOfVertices();
	distance_.resize( size, 0 );
	probability_.resize( size, 0.0 );
	next_probability_.resize( size, 0.0 );
	degree_.resize( size, 0 );
	for(int i=1;i<=100;i++)
	{
		diffusion_times_.push_back(i);
	}
	for(int i=11;2*i*i<triangulation_->NumberOfVertices();i++)
	{
		diffusion_times_.push_back(i*i);
	}
	distribution_.resize(diffusion_times_.size(),std::vector<double>(max_distance_+1,0.0));

}


Diffusion::~Diffusion(void)
{
}

void Diffusion::Measure()
{
	DetermineDegrees();
	DiffusionMatrix diffusionmatrix( triangulation_, degree_);
	
	for(int i=0;i<samples_;i++)
	{
		Vertex * startVertex = triangulation_->getRandomVertex();
		DetermineDistance( startVertex );
		MeasureDistanceDistribution();
		if( !distance_only_ )
		{
			DoDiffusion( startVertex, &diffusionmatrix );
		}
		measurements_++;
	}
}

void Diffusion::DetermineDistance(Vertex * startVertex)
{
	properties::VertexDistanceList(triangulation_,startVertex,distance_);
}

void Diffusion::DetermineDegrees()
{
	properties::DegreeList(triangulation_,degree_);
}

void Diffusion::DoDiffusion(Vertex * startVertex, const linearalgebra::Matrix * const matrix )
{
	std::fill(probability_.begin(),probability_.end(),0.0);
	probability_[startVertex->getId()] = 1.0;
	int previousTime=0;
	for(int i=0,end=diffusion_times_.size();i<end;i++)
	{
		for(int j=previousTime;j<diffusion_times_[i];j++)
		{
			matrix->MultiplyVector(probability_,next_probability_);
			linearalgebra::Copy(next_probability_,probability_);
		}
		previousTime = diffusion_times_[i];
		DoMeasurementOnDistribution(i);
	}
}

void Diffusion::DoMeasurementOnDistribution(int time)
{
	for(int j=0,end=distance_.size();j<end;j++)
	{
		if( distance_[j] <= max_distance_ )
			distribution_[time][distance_[j]] += probability_[j];
	}
}

void Diffusion::MeasureDistanceDistribution()
{
	for(int i=0,end=distance_.size();i<end;i++)
	{
		ResizeAndAdd(distance_distribution_,distance_[i],1);
	}
}

std::string Diffusion::OutputData() const
{
	std::ostringstream stream;
	stream << std::fixed << "diffusion -> {measurements -> " << measurements_ << ", walkprobability -> " << 1.0;
	stream << ", samples -> " << samples_;
	if( !distance_only_ )
	{
		stream << ", diffusiontimes -> ";
		PrintToStream(stream,diffusion_times_.begin(),diffusion_times_.end());
		stream << ", distribution -> ";
		PrintToStream2D(stream,distribution_.begin(),distribution_.end());
		stream << "/" << measurements_;
	}
	stream << ", distancedistribution -> ";
	PrintToStream(stream,distance_distribution_.begin(),distance_distribution_.end());
	stream << "/" << measurements_ << "}";
	return stream.str();
}

void Diffusion::setMeasureDistanceOnly(bool distanceonly)
{
	distance_only_ = distanceonly;
}