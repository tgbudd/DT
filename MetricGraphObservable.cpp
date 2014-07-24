#include <queue>
#include <algorithm>

#include "MetricGraphObservable.h"
#include "TriangulationProperties.h"

MetricGraphObservable::MetricGraphObservable(const Triangulation * triangulation, double MaxLength, int LengthBins, int samples) :
	triangulation_(triangulation),
	dual_edge_length_(triangulation),
	distance_(triangulation),
	num_edges_(triangulation),
	samples_(samples),
	length_max_(MaxLength),
	length_bins_(LengthBins),
	length_min_(0.0),
	distance_is_zero_(0)
{
	int maxEdges = 2*static_cast<int>(length_max_);
	edges_per_length_.resize(length_bins_,Histogram<int>(0,maxEdges,maxEdges));
	dist_per_length_.resize(length_bins_,Histogram<int>(0,maxEdges,maxEdges));
}

MetricGraphObservable::~MetricGraphObservable()
{ 
}

void MetricGraphObservable::Measure()
{
	setRandomEdgeLengths();

	for(int i=0;i<samples_;i++)
	{
		MeasureTwoPointDistance();
	}
}

void MetricGraphObservable::setRandomEdgeLengths()
{
	const double logLengthMinimum = 1.0e-10;
	for(int i=0,endi=triangulation_->NumberOfTriangles();i<endi;++i)
	{
		const Triangle * triangle = triangulation_->getTriangle(i);
		for(int j=0;j<3;j++)
		{
			const Edge * edge = triangle->getEdge(j);
			double explength = -std::log(std::max(logLengthMinimum,triangulation_->RandomReal(0.0,1.0)));
			dual_edge_length_[edge] = explength;
			dual_edge_length_[edge->getAdjacent()] = explength;
		}
	}
}

void MetricGraphObservable::MeasureTwoPointDistance()
{
	const Triangle * t1 = triangulation_->getRandomTriangle();
	const Triangle * t2 = triangulation_->getRandomTriangle();

	if( t1 == t2 )
	{
		distance_is_zero_++;
	} else
	{
		std::pair<double,int> distance = Distance(t1,t2);

		int discrete = properties::TriangleDistance(triangulation_,t1,t2);

		int lengthbin = static_cast<int>((distance.first-length_min_)/(length_max_-length_min_)*length_bins_);
		if( lengthbin >= 0 && lengthbin < length_bins_ )
		{
			edges_per_length_[lengthbin].Insert(distance.second);
			dist_per_length_[lengthbin].Insert(discrete);
		}
	}
}


std::pair<double, int> MetricGraphObservable::Distance(const Triangle * t1, const Triangle * t2)
{
	std::fill(distance_.begin(),distance_.end(),100.0 * length_max_);
	std::fill(num_edges_.begin(),num_edges_.end(),0);

	typedef std::pair<const Triangle *,double> TriangleDist;

	std::priority_queue<TriangleDist, std::vector<TriangleDist>, TriangleDistanceComparator> queue;
	distance_[t1] = 0.0;
	queue.push(TriangleDist(t1,0.0));
	while(!queue.empty())
	{
		TriangleDist triangle = queue.top();
		queue.pop();
		if( triangle.first == t2 )
		{
			return std::pair<double,int>(distance_[t2],num_edges_[t2]);
		}

		for(int i=0;i<3;i++)
		{
			TriangleDist neighbour;
			neighbour.first = triangle.first->getEdge(i)->getAdjacent()->getParent();
			neighbour.second = triangle.second + dual_edge_length_[triangle.first->getEdge(i)];
			if( neighbour.second < distance_[neighbour.first] )
			{
				distance_[neighbour.first] = neighbour.second;
				num_edges_[neighbour.first] = num_edges_[triangle.first] + 1;
				queue.push(neighbour);
			}
		}
	}
	return std::pair<double,int>(-1.0,-1);
}

std::string MetricGraphObservable::OutputData() const
{
	std::ostringstream stream;
	stream << std::fixed << "metricgraphobservable -> {";
	stream << "samples -> " << samples_;
	stream << ", lengthmin -> " << length_min_;
	stream << ", lengthmax -> " << length_max_;
	stream << ", lengthbins -> " << length_bins_;
	stream << ", distanceiszero -> " << distance_is_zero_;
	stream << ", edgesperlength -> {";
	for(int i=0,endi=edges_per_length_.size();i<endi;++i)
	{
		stream << (i>0?",":"");
		edges_per_length_[i].PrintTo(stream);
	}
	stream << "}, distperlength -> {";
	for(int i=0,endi=dist_per_length_.size();i<endi;++i)
	{
		stream << (i>0?",":"");
		dist_per_length_[i].PrintTo(stream);
	}
	stream << "}}";
	return stream.str();
}