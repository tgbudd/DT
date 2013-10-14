#include <queue>

#include "DistanceDistribution.h"


DistanceDistribution::DistanceDistribution(const Triangulation * const triangulation)
	: triangulation_(triangulation), max_distance_(100), in_universal_cover_(false)
	, genus_(triangulation->CalculateGenus()), samples_(10), measurements_(0)
{
	distance_distribution_.resize(max_distance_+1,0);
	distance_.resize(triangulation_->NumberOfVertices());
}


DistanceDistribution::~DistanceDistribution(void)
{
}

void DistanceDistribution::UseUniversalCover(CohomologyBasis * cohomologybasis)
{
	in_universal_cover_ = true;
	cohomologybasis_ = cohomologybasis;
}

void DistanceDistribution::Measure()
{
	if( genus_ == 1 && in_universal_cover_ )
	{
		for(int i=0;i<samples_;i++)
		{
			MeasureDistanceInUniversalCover();
			measurements_++;
		}
	}
}

std::string DistanceDistribution::OutputData() const
{
	std::ostringstream stream;
	stream << std::fixed << "distance -> {measurements -> " << measurements_;
	stream << ", samples -> " << samples_;
	stream << ", distribution -> ";
	PrintToStream(stream,distance_distribution_.begin(),distance_distribution_.end());
	stream << "/" << measurements_ << "}";
	return stream.str();
}

void DistanceDistribution::MeasureDistanceInUniversalCover()
{
	Vertex * startVertex = triangulation_->getRandomVertex();
	std::queue<std::pair<Vertex *,IntForm2D> > queue;
	IntForm2D ZeroForm = {0,0};
	queue.push(std::pair<Vertex*,IntForm2D>(startVertex,ZeroForm));

	while( !queue.empty() )
	{
		std::pair<Vertex *,IntForm2D> vertexnode = queue.front();
		queue.pop();
		int distance = distance_[vertexnode.first->getId()][vertexnode.second];

		distance_distribution_[distance]++;

		if( distance == max_distance_ )
		{
			continue;
		}

		Edge * edge = vertexnode.first->getParent()->getPrevious();
		do {
			std::pair<Vertex *,IntForm2D> nbr( edge->getPrevious()->getOpposite(), AddForms( vertexnode.second, cohomologybasis_->getOmega(edge) ) );
			std::pair< std::map<IntForm2D,int>::iterator, bool> insertresult = distance_[nbr.first->getId()].insert( std::pair<IntForm2D,int>(nbr.second,distance+1) ); 
			if( insertresult.second )
			{
				queue.push( nbr );
			}
			edge = edge->getAdjacent()->getNext();
		} while( edge != vertexnode.first->getParent()->getPrevious() );
	}
	for(int i=0;i<triangulation_->NumberOfVertices();i++)
	{
		distance_[i].clear();
	}
}