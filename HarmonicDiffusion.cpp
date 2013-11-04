#include <queue>

#include "HarmonicDiffusion.h"
#include "TriangulationProperties.h"

HarmonicDiffusion::HarmonicDiffusion(const Triangulation * const triangulation, Embedding * const embedding, CohomologyBasis * const cohomologybasis)
	: triangulation_(triangulation), 
	  embedding_(embedding),
	  cohomologybasis_(cohomologybasis),
	  samples_(200),
	  measurements_(0),
	  walk_time_(triangulation->NumberOfTriangles()/10),
	  max_graph_distance_(35),
	  graph_distance_samples_(5)
{
	average_distance_.resize(walk_time_,0.0);
	square_distance_.resize(walk_time_,0.0);
	graph_distance_distribution_.resize(max_graph_distance_+1,0);
	average_distance_per_graph_distance_.resize(max_graph_distance_+1,0.0);
	square_distance_per_graph_distance_.resize(max_graph_distance_+1,0.0);
	graph_distance_.resize(triangulation_->NumberOfVertices());
}


HarmonicDiffusion::~HarmonicDiffusion(void)
{
}

void HarmonicDiffusion::Measure()
{
	if( !embedding_->IsUpToDate() )
	{
		embedding_->MakeUpToDate();
	}
	moduli_ = embedding_->CalculateModuli();

	for(int i=0;i<samples_;i++)
	{
		DoRandomWalk();
		measurements_++;
	}
	for(int i=0;i<graph_distance_samples_;i++)
	{
		DoGraphDistanceMeasurement();
	}
}

void HarmonicDiffusion::DoRandomWalk()
{
	Vertex * vertex = triangulation_->getRandomVertex();
	Vector2D position = {0.0,0.0};

	for(int t=0;t<walk_time_;t++)
	{
		Edge * edge = RandomEdge(vertex);
		
		position = AddVectors2D(position,embedding_->getForm(edge));
		Vector2D absPosition = {(position[0] + moduli_.first * position[1])/std::sqrt(moduli_.second), position[1] * std::sqrt(moduli_.second)};
		double squareDistance = absPosition[0]*absPosition[0] + absPosition[1]*absPosition[1];
		double distance = std::sqrt(squareDistance);

		average_distance_[t] += distance;
		square_distance_[t] += squareDistance;

		vertex = edge->getPrevious()->getOpposite();
	}
}

Edge * HarmonicDiffusion::RandomEdge(Vertex * v)
{
	Edge * edge = v->getParent()->getPrevious();
	int degree = 0;
	do {
		degree++;
		edge = edge->getAdjacent()->getNext();
	} while( edge != v->getParent()->getPrevious() );
	int randomedge = triangulation_->RandomInteger(0,degree-1);
	degree = 0;
	do {
		if( degree == randomedge )
		{
			break;
		}
		degree++;
		edge = edge->getAdjacent()->getNext();
	} while( edge != v->getParent()->getPrevious() );
	return edge;
}

std::string HarmonicDiffusion::OutputData() const
{
	std::ostringstream stream;
	stream << std::fixed << "harmonicdiffusion -> {measurements -> " << measurements_;
	stream << ", averagedistance -> ";
	PrintToStream(stream,average_distance_.begin(),average_distance_.end());
	stream << ", squaredistance -> ";
	PrintToStream(stream,square_distance_.begin(),square_distance_.end());
	stream << ", averagedistancepergraphdistance -> ";
	PrintToStream(stream,average_distance_per_graph_distance_.begin(),average_distance_per_graph_distance_.end());
	stream << ", squaredistancepergraphdistance -> ";
	PrintToStream(stream,square_distance_per_graph_distance_.begin(),square_distance_per_graph_distance_.end());
	stream << ", graphdistancedistribution -> ";
	PrintToStream(stream,graph_distance_distribution_.begin(),graph_distance_distribution_.end());
	stream << "}";
	return stream.str();
}

void HarmonicDiffusion::DoGraphDistanceMeasurement()
{
	for(int i=0;i<triangulation_->NumberOfVertices();i++)
	{
		graph_distance_[i].clear();
	}

	Vertex * startVertex = triangulation_->getRandomVertex();

	std::queue<std::pair<Vertex *,std::pair<IntForm2D,Vector2D> > > queue;
	IntForm2D ZeroForm = {0,0};
	Vector2D ZeroVector = {0.0,0.0};
	queue.push(std::pair<Vertex*,std::pair<IntForm2D,Vector2D>>(startVertex,std::pair<IntForm2D,Vector2D>(ZeroForm,ZeroVector)));

	while( !queue.empty() )
	{
		std::pair<Vertex *,std::pair<IntForm2D,Vector2D> > vertexnode = queue.front();
		queue.pop();
		int graphDistance = graph_distance_[vertexnode.first->getId()][vertexnode.second.first];

		Vector2D absPosition = {(vertexnode.second.second[0] + moduli_.first * vertexnode.second.second[1])/std::sqrt(moduli_.second), vertexnode.second.second[1] * std::sqrt(moduli_.second)};
		double squareDistance = absPosition[0]*absPosition[0] + absPosition[1]*absPosition[1];
		double distance = std::sqrt(squareDistance);


		graph_distance_distribution_[graphDistance]++;
		average_distance_per_graph_distance_[graphDistance] += distance;
		square_distance_per_graph_distance_[graphDistance] += squareDistance;

		if( graphDistance == max_graph_distance_ )
		{
			continue;
		}


		Edge * edge = vertexnode.first->getParent()->getPrevious();
		do {
			std::pair<Vertex *,std::pair<IntForm2D,Vector2D> > nbr( edge->getPrevious()->getOpposite(), 
				std::pair<IntForm2D,Vector2D>( AddForms( vertexnode.second.first, cohomologybasis_->getOmega(edge) ) ,
				AddVectors2D( vertexnode.second.second, embedding_->getForm(edge) ) ) );
			std::pair< std::map<IntForm2D,int>::iterator, bool> insertresult = graph_distance_[nbr.first->getId()].insert( std::pair<IntForm2D,int>(nbr.second.first,graphDistance+1) ); 
			if( insertresult.second )
			{
				queue.push( nbr );
			}
			edge = edge->getAdjacent()->getNext();
		} while( edge != vertexnode.first->getParent()->getPrevious() );
	}
}