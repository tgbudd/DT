#include <iostream>
#include <sstream>
#include <iomanip>
#include <queue>
#include <set>

#include "Simulation.h"
#include "Triangulation.h"
#include "utilities.h"
#include "CMinusTwoBuilder.h"
#include "SimplexAttribute.h"
#include "TriangulationProperties.h"
#include "Histogram.h"

class TriangulationExploration : public Observable {
public:
	TriangulationExploration(const Triangulation * triangulation, int maxlength, int lengthbins, int maxedges, int edgebins)
		: triangulation_(triangulation), dual_distance_(triangulation), distance_(triangulation),
		explored_(triangulation), edge_length_(triangulation),
		maxlength_(maxlength), lengthbins_(lengthbins), maxedges_(maxedges), edgebins_(edgebins),
		max_logs_(20), log_count_(0)
	{
		histogram_.resize(lengthbins,Histogram<int>(0,maxedges,edgebins));
	}
	void Measure();
	std::string OutputData() const;
private:
	const Triangulation * triangulation_;
	TriangleAttribute<int> dual_distance_;
	VertexAttribute<std::pair<double,int> > distance_;
	EdgeAttribute<double> edge_length_;
	EdgeAttribute<bool> explored_;

	int maxlength_, lengthbins_;
	int maxedges_, edgebins_;
	std::vector<Histogram<int> > histogram_;

	void setRandomEdgeLengths();
	std::pair<int,int> getFrontierSize(double time);
	void InspectExploration();

	int MinDistAlongFrontier(const EdgeAttribute<bool> & explored, const Edge * startEdge) const;
	int LengthFrontier(const EdgeAttribute<bool> & explored, const Edge * startEdge) const;
	int IngoingEdgesFrontier(const EdgeAttribute<bool> & explored, const Edge * startEdge) const;
	void SetComponentToTrue(EdgeAttribute<bool> & explored, const Edge * startEdge) const;

	const Vertex * initial_vertex_;
	const Triangle * final_triangle_;

	int max_logs_;
	int log_count_;
	std::ostringstream logstream_;
};

void TriangulationExploration::setRandomEdgeLengths()
{
	const double logLengthMinimum = 1.0e-10;
	for(int i=0,endi=triangulation_->NumberOfTriangles();i<endi;++i)
	{
		const Triangle * triangle = triangulation_->getTriangle(i);
		for(int j=0;j<3;j++)
		{
			const Edge * edge = triangle->getEdge(j);
			double explength = -std::log(std::max(logLengthMinimum,triangulation_->RandomReal(0.0,1.0)));
			edge_length_[edge] = explength;
			edge_length_[edge->getAdjacent()] = explength;
		}
	}
}

void TriangulationExploration::Measure()
{
	initial_vertex_ = triangulation_->getRandomVertex();
	final_triangle_ = triangulation_->getRandomTriangle();
	setRandomEdgeLengths();
	properties::VertexWeightedDistanceList(triangulation_,initial_vertex_,edge_length_,distance_);
	properties::TriangleDistanceList(triangulation_,final_triangle_, dual_distance_ );

	double dist = 1.0e10;
	for(int i=0;i<3;i++)
	{
		Vertex * vertex = final_triangle_->getEdge(i)->getOpposite();
		if( distance_[vertex].first < dist )
		{
			dist = distance_[vertex].first;
		}
	}

	if( dist > 1.0 )
	{
		InspectExploration();
		
		std::pair<int,int> size = getFrontierSize(0.5*dist);

		if( size.first < maxlength_ )
		{
			histogram_[size.first*lengthbins_/maxlength_].Insert(size.second);
		}
	}
}

std::pair<int,int> TriangulationExploration::getFrontierSize(double time)
{
	for(EdgeAttribute<bool>::iterator edgeIt = explored_.begin(); edgeIt != explored_.end(); edgeIt++)
	{
		const Edge * edge = explored_.getSubject(edgeIt);
		double dist1 = distance_[edge->getNext()->getOpposite()].first;
		double dist2 = distance_[edge->getPrevious()->getOpposite()].first;
		*edgeIt = 0.5*(dist1 + dist2 + edge_length_[edge]) < time;
	}

	int length = 0;
	int frontieredges = 0;

	TriangleAttribute<bool> visited(triangulation_,false);

	std::queue<const Triangle *> queue;
	queue.push(final_triangle_);
	visited[final_triangle_]=true;
	while( !queue.empty() )
	{
		const Triangle * triangle = queue.front();
		queue.pop();

		for(int i=0;i<3;i++)
		{
			const Edge * edge = triangle->getEdge(i);

			if( explored_[edge] )
			{
				length++;
			} else
			{
				double distanceStartVertex = distance_[edge->getNext()->getOpposite()].first;

				if( distanceStartVertex  < time )
				{
					frontieredges++;
				}

				const Triangle * nbr = edge->getAdjacent()->getParent();
				if( !visited[nbr] )
				{
					visited[nbr] = true;
					queue.push(nbr);
				}
			}
		}
	}
	return std::pair<int,int>(length,frontieredges);
}

template< typename A, typename B >
class OrderBy {
public:
	OrderBy( const B & b) : b_(b) {}
	bool operator()(const A & a1, const A & a2) {
		return b_[a1] < b_[a2];
	}
private:
	const B & b_;
};

void TriangulationExploration::InspectExploration()
{
	EdgeAttribute<double> edgeDistance(triangulation_);
	std::vector<const Edge *> edges;
	for(EdgeAttribute<double>::iterator edgeIt = edgeDistance.begin(); edgeIt != edgeDistance.end(); ++edgeIt)
	{
		const Edge * edge = edgeDistance.getSubject(edgeIt);
		double dist1 = distance_[edge->getEnd()].first;
		double dist2 = distance_[edge->getStart()].first;
		*edgeIt = 0.5*(dist1 + dist2 + edge_length_[edge]);
		if( edge->getParent()->getId() < edge->getAdjacent()->getParent()->getId() ||
			(edge->getParent()->getId() == edge->getAdjacent()->getParent()->getId() &&
			edge->getId() < edge->getAdjacent()->getId() ) )
		{
			edges.push_back(edge);
		}
	}

	std::sort(edges.begin(),edges.end(),OrderBy<const Edge *,EdgeAttribute<double> >(edgeDistance));

	bool logthis = false;
	if( log_count_ < max_logs_ )
	{
		logthis = true;
		log_count_++;
	}
	if( logthis )
	{
		logstream_ << (log_count_ > 1 ?",":"") << "{";
	}

	VertexAttribute<bool> vertexExplored(triangulation_,false);
	EdgeAttribute<bool> edgeExplored(triangulation_,false);
	vertexExplored[initial_vertex_] = true;
	int length = 0;
	int ingoing = initial_vertex_->getDegree();
	if( logthis )
	{
		logstream_ << "{" << length << "," << ingoing << "}";
	}
	for(std::vector<const Edge *>::iterator edgeIt = edges.begin();edgeIt!= edges.end();++edgeIt)
	{
		if( (*edgeIt)->getParent() == final_triangle_ || (*edgeIt)->getAdjacent()->getParent() == final_triangle_ )
		{
			break;
		}
		if( edgeExplored[*edgeIt] )
		{
			continue;
		}
		edgeExplored[*edgeIt] = true;
		edgeExplored[(*edgeIt)->getAdjacent()] = true;
		const Edge * frontierEdge = *edgeIt;
		if( vertexExplored[(*edgeIt)->getEnd()] && vertexExplored[(*edgeIt)->getStart()] )
		{
			// Need to determine which side of the edge is the component containing final_triangle_.
			// To do this, look for the triangle adjacent to an explored edge that has minimal distance to final_triangle_.
			if( MinDistAlongFrontier(edgeExplored,*edgeIt) > MinDistAlongFrontier(edgeExplored,(*edgeIt)->getAdjacent()) )
			{
				frontierEdge = (*edgeIt)->getAdjacent();
			}
			length = LengthFrontier(edgeExplored,frontierEdge);
			ingoing = IngoingEdgesFrontier(edgeExplored,frontierEdge);
			SetComponentToTrue(edgeExplored,frontierEdge->getAdjacent());
		} else
		{
			length += 2;
			if( vertexExplored[(*edgeIt)->getEnd()] )
			{
				ingoing += (*edgeIt)->getStart()->getDegree()-2;
				vertexExplored[(*edgeIt)->getStart()] = true;
			} else
			{
				ingoing += (*edgeIt)->getEnd()->getDegree()-2;
				vertexExplored[(*edgeIt)->getEnd()] = true;
			}
		}
		/*
		int lengthshouldbe = LengthFrontier(edgeExplored,frontierEdge);
		int ingoingshouldbe = IngoingEdgesFrontier(edgeExplored,frontierEdge);
		BOOST_ASSERT( length == lengthshouldbe );
		BOOST_ASSERT( ingoing == ingoingshouldbe );
		*/
		if( logthis )
		{
			logstream_ << ",{" << length << "," << ingoing << "}";
		}
	}
	if( logthis )
	{
		logstream_ << "}";
	}
	int hy=0;
}

int TriangulationExploration::MinDistAlongFrontier(const EdgeAttribute<bool> & explored, const Edge * startEdge) const
{
	int mindist = 1000000;
	const Edge * edge = startEdge;
	do
	{
		int dist = dual_distance_[edge->getParent()];
		if( dist < mindist )
		{
			mindist = dist;
		}

		edge = edge->getNext();
		while( !explored[edge] )
		{
			edge = edge->getAdjacent()->getNext();
		}
	} while(edge != startEdge);
	return mindist;
}

int TriangulationExploration::LengthFrontier(const EdgeAttribute<bool> & explored, const Edge * startEdge) const
{
	int length=0;
	const Edge * edge = startEdge;
	do
	{
		edge = edge->getNext();
		while( !explored[edge]  )
		{
			edge = edge->getAdjacent()->getNext();
		}
		length++;
	} while(edge != startEdge);
	return length;
}

int TriangulationExploration::IngoingEdgesFrontier(const EdgeAttribute<bool> & explored, const Edge * startEdge) const
{
	int ingoing = 0;
	const Edge * edge = startEdge;
	do
	{
		int dist = dual_distance_[edge->getParent()];
		edge = edge->getNext();
		while( !explored[edge]  )
		{
			edge = edge->getAdjacent()->getNext();
			ingoing++;
		}
	} while(edge != startEdge);
	return ingoing;
}

void TriangulationExploration::SetComponentToTrue(EdgeAttribute<bool> & explored, const Edge * startEdge) const
{
	std::queue<const Triangle *> queue;
	std::set<const Triangle *> set;
	queue.push(startEdge->getParent());
	set.insert(startEdge->getParent());
	while( !queue.empty() )
	{
		const Triangle * triangle = queue.front();
		queue.pop();

		for(int i=0;i<3;i++)
		{
			const Edge * edge = triangle->getEdge(i);
			if( !explored[edge] )
			{
				explored[edge] = true;
				explored[edge->getAdjacent()] = true;
				if( set.insert(edge->getAdjacent()->getParent()).second )
				{
					queue.push(edge->getAdjacent()->getParent());
				}
			}
		}
	}
}

std::string TriangulationExploration::OutputData() const
{
	std::ostringstream os;
	os << "exploration -> {";
	os << "maxlength -> " << maxlength_;
	os << ", lengthbins -> " << lengthbins_;
	os << ", histograms -> {";
	for(int i=0,endi=histogram_.size();i<endi;++i)
	{
		os << (i>0?",":"");
		histogram_[i].PrintTo(os);
	}
	os << "}, log -> {" << logstream_.str() << "}}";
	return os.str();
}

int main(int argc, char* argv[])
{
	ParameterStream param(argc,argv);
	
	Triangulation triangulation;
	int matter = param.Read<int>("matter (0=none,1=minus2)");
	int thermalizationSweeps=0;
	int measurementSweeps=1;
	int n = param.Read<int>("triangles");
	if( matter == 0 )
	{
		thermalizationSweeps = param.Read<int>("thermalization sweeps");
		measurementSweeps = param.Read<int>("measurement sweeps");
	}
	int maxlength = param.Read<int>("max length");
	int lengthbins = param.Read<int>("length bins");
	int maxedges = param.Read<int>("max edges");
	int edgebins = param.Read<int>("edge bins");
	int secondsperoutput = param.Read<int>("seconds per output");
	bool output = param.UserInput();

	CMinusTwoBuilder builder(&triangulation,0,n);
	triangulation.setDominantMatter( &builder );
	triangulation.DoSweep();
	if( matter == 0 )
	{
		triangulation.clearDominantMatter();
	}

	Simulation simulation( &triangulation, thermalizationSweeps, secondsperoutput, output );

	TriangulationExploration exploration( &triangulation, maxlength, lengthbins,maxedges,edgebins );
	
	simulation.AddObservable( &exploration, measurementSweeps );
	simulation.SetDirectory("./output/");
	simulation.Run();
	return 0;
}