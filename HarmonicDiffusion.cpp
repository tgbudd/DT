#include <queue>

#include "HarmonicDiffusion.h"
#include "TriangulationProperties.h"

HarmonicDiffusion::HarmonicDiffusion(const Triangulation * const triangulation, Embedding * const embedding, CohomologyBasis * const cohomologybasis)
	: triangulation_(triangulation), 
	  embedding_(embedding),
	  cohomologybasis_(cohomologybasis),
	  samples_(3),
	  measurements_(0),
	  walk_time_(triangulation->NumberOfTriangles()/10),
	  max_graph_distance_(35),
	  graph_distance_samples_(5),
	  max_trajectories_(5),
	  max_wiener_trajectories_(10000),
	  time_step_(0.0005),
	  num_wiener_steps_(200)
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
	for(int i=0;i<samples_;i++)
	{
		DoRandomWienerWalk();
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

	bool trajectory = ( max_trajectories_ > static_cast<int>(wiener_time_.size()) );

	if( trajectory )
	{
		wiener_time_.push_back(std::vector<double>());
		wiener_time_.back().reserve(walk_time_);
		euclidean_distance_.push_back(std::vector<double>());
		euclidean_distance_.back().reserve(walk_time_);
	}

	double wienerTime = 0.0;
	for(int t=0;t<walk_time_;t++)
	{
		Edge * edge = RandomEdge(vertex);
		
		position = AddVectors2D(position,embedding_->getForm(edge));
		Vector2D absPosition = {(position[0] + moduli_.first * position[1])/std::sqrt(moduli_.second), position[1] * std::sqrt(moduli_.second)};
		double squareDistance = absPosition[0]*absPosition[0] + absPosition[1]*absPosition[1];
		double distance = std::sqrt(squareDistance);

		average_distance_[t] += distance;
		square_distance_[t] += squareDistance;

		if( trajectory )
		{
			euclidean_distance_.back().push_back(distance);
			Vector2D v0 = embedding_->getForm(edge);
			Vector2D v1 = NegateVector2D(embedding_->getForm(edge->getPrevious()));
			Vector2D v2 = embedding_->getForm(edge->getAdjacent()->getPrevious());
			double area = std::fabs( 0.5 * ( v0[0] * (v1[1]+v2[1]) - v0[1] * (v1[0]+v2[0]) ) );
			wienerTime += area;
			wiener_time_.back().push_back( wienerTime );
		}

		vertex = edge->getPrevious()->getOpposite();
	}
}

void HarmonicDiffusion::DoRandomWienerWalk()
{
	Vertex * vertex = triangulation_->getRandomVertex();
	Vector2D position = {0.0,0.0};

	bool trajectory = ( max_wiener_trajectories_ > static_cast<int>(walking_time_.size()) );

	if( trajectory )
	{
		walking_time_.push_back(std::vector<int>());
		walking_time_.back().reserve(num_wiener_steps_);
		wiener_euclidean_distance_.push_back(std::vector<double>());
		wiener_euclidean_distance_.back().reserve(num_wiener_steps_);
	}

	double wienerTime = 0.0;
	int time=0;
	for(int t=0;t<num_wiener_steps_;t++)
	{
		while( wienerTime < (t+1)*time_step_ )
		{
			Edge * edge = RandomEdge(vertex);
		
			position = AddVectors2D(position,embedding_->getForm(edge));

			Vector2D v0 = embedding_->getForm(edge);
			Vector2D v1 = NegateVector2D(embedding_->getForm(edge->getPrevious()));
			Vector2D v2 = embedding_->getForm(edge->getAdjacent()->getPrevious());
			double area = std::fabs( 0.5 * ( v0[0] * (v1[1]+v2[1]) - v0[1] * (v1[0]+v2[0]) ) );
			wienerTime += area;
			time++;
			vertex = edge->getPrevious()->getOpposite();
		}

		Vector2D absPosition = {(position[0] + moduli_.first * position[1])/std::sqrt(moduli_.second), position[1] * std::sqrt(moduli_.second)};
		double squareDistance = absPosition[0]*absPosition[0] + absPosition[1]*absPosition[1];
		double distance = std::sqrt(squareDistance);

		if( trajectory )
		{
			wiener_euclidean_distance_.back().push_back(distance);
			walking_time_.back().push_back( time );
		}
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
	stream << ", wienertime -> ";
	PrintToStream2D(stream,wiener_time_.begin(),wiener_time_.end());
	stream << ", euclideandistance -> ";
	PrintToStream2D(stream,euclidean_distance_.begin(),euclidean_distance_.end());
	stream << ", walkingtime -> ";
	PrintToStream2D(stream,walking_time_.begin(),walking_time_.end());
	stream << ", wienereuclideandistance -> ";
	PrintToStream2D(stream,wiener_euclidean_distance_.begin(),wiener_euclidean_distance_.end());
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
	queue.push(std::pair<Vertex*,std::pair<IntForm2D,Vector2D> >(startVertex,std::pair<IntForm2D,Vector2D>(ZeroForm,ZeroVector)));

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
