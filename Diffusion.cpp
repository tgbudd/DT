#include <queue>

#include "utilities.h"
#include "Triangle.h"
#include "Edge.h"
#include "Diffusion.h"

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
	samples_(20), measurements_(0), diffusion_steps_(100), max_distance_(50)
{
	int size = triangulation_->NumberOfVertices();
	distance_.resize( size, 0 );
	probability_.resize( size, 0.0 );
	next_probability_.resize( size, 0.0 );
	degree_.resize( size, 0 );
	distribution_.resize(diffusion_steps_,std::vector<double>(max_distance_+1,0.0));

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
		DoDiffusion( startVertex, &diffusionmatrix );
		measurements_++;
	}
}

void Diffusion::DetermineDistance(Vertex * startVertex)
{
	std::fill(distance_.begin(),distance_.end(),-1);

	std::queue<Vertex *> q;
	q.push( startVertex );
	distance_[startVertex->getId()] = 0;

	while( !q.empty() )
	{
		Vertex * vertex = q.front();
		q.pop();

		Edge * edge = vertex->getParent()->getPrevious();
		do {
			Vertex * nbr = edge->getPrevious()->getOpposite();

			if( distance_[nbr->getId()] == -1 )
			{
				distance_[nbr->getId()] = distance_[vertex->getId()] + 1;
				q.push(nbr);
			}
			edge = edge->getPrevious()->getAdjacent();
		} while( edge->getNext() != vertex->getParent() );
	}
}

void Diffusion::DetermineDegrees()
{
	for(int i=0,end=triangulation_->NumberOfVertices();i<end;i++)
	{
		Vertex * vertex = triangulation_->getVertex(i);
		Edge * edge = vertex->getParent()->getPrevious();
		degree_[i]=0;
		do {
			edge = edge->getPrevious()->getAdjacent();
			++degree_[i];
		} while( edge->getNext() != vertex->getParent() );
	}
}

void Diffusion::DoDiffusion(Vertex * startVertex, const linearalgebra::Matrix * const matrix )
{
	std::fill(probability_.begin(),probability_.end(),0.0);
	probability_[startVertex->getId()] = 1.0;

	for(int i=0;i<diffusion_steps_;i++)
	{
		BOOST_ASSERT( std::fabs( linearalgebra::Total(probability_) - 1.0 ) < 1e-6 );
		matrix->MultiplyVector(probability_,next_probability_);
		linearalgebra::Copy(next_probability_,probability_);
		DoMeasurementOnDistribution(i+1);
	}
}

void Diffusion::DoMeasurementOnDistribution(int time)
{
	for(int j=0,end=distance_.size();j<end;j++)
	{
		if( distance_[j] <= max_distance_ )
			distribution_[time-1][distance_[j]] += probability_[j];
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
	stream << ", diffusionsteps -> " << diffusion_steps_ << ", samples -> " << samples_;
	stream << ", distribution -> ";
	PrintToStream2D(stream,distribution_.begin(),distribution_.end());
	stream << "/" << measurements_ << ", distancedistribution -> ";
	PrintToStream(stream,distance_distribution_.begin(),distance_distribution_.end());
	stream << "/" << measurements_ << "}";
	return stream.str();
}