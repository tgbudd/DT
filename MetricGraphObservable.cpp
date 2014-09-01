#include <queue>
#include <algorithm>

#include "MetricGraphObservable.h"
#include "TriangulationProperties.h"

MetricGraphObservable::MetricGraphObservable(const Triangulation * triangulation, double MaxLength, int LengthBins, int samples) :
	triangulation_(triangulation),
	dual_edge_length_(triangulation),
	distance_(triangulation),
	num_edges_(triangulation),
	v_distance_(triangulation),
	v_num_edges_(triangulation),
	samples_(samples),
	length_max_(MaxLength),
	length_bins_(LengthBins),
	length_min_(0.0),
	distance_is_zero_(0),
	vertex_degree_(0,50,50),
	vertex_degree_geodesic_(0,50,50)
{
	int maxEdges = 2*static_cast<int>(length_max_);
	edges_per_length_.resize(length_bins_,Histogram<int>(0,maxEdges,maxEdges));
	dist_per_length_.resize(length_bins_,Histogram<int>(0,maxEdges,maxEdges));
	tri_dist_per_length_.resize(length_bins_,Histogram<int>(0,maxEdges,maxEdges));
	tri_dist_per_dist_.resize(maxEdges,Histogram<int>(0,maxEdges,maxEdges));
	tri_time_per_length_.resize(maxEdges,Histogram<double>(0,0.25*maxEdges,maxEdges));
	tri_hop_per_length_.resize(maxEdges,Histogram<int>(0,maxEdges,maxEdges));
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
		//MeasureDegreeOfVerticesAlongGeodesic(triangulation_->getRandomVertex(),triangulation_->getRandomVertex());
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
	
		const Vertex * v1 = t1->getEdge(triangulation_->RandomInteger(0,2))->getOpposite();
		const Vertex * v2 = t2->getEdge(triangulation_->RandomInteger(0,2))->getOpposite();
		int tri_discrete = properties::VertexDistance(triangulation_,v1,v2);
		std::pair<double,int> tri_distance = TriangulationDistance(v1,v2);

		int lengthbin = static_cast<int>((distance.first-length_min_)/(length_max_-length_min_)*length_bins_);
		if( lengthbin >= 0 && lengthbin < length_bins_ )
		{
			edges_per_length_[lengthbin].Insert(distance.second);
			dist_per_length_[lengthbin].Insert(discrete);
			tri_dist_per_length_[lengthbin].Insert(tri_discrete);
			tri_time_per_length_[lengthbin].Insert(tri_distance.first);
			tri_hop_per_length_[lengthbin].Insert(tri_distance.second);
		}

		if( discrete < static_cast<int>(tri_dist_per_dist_.size()) )
		{
			tri_dist_per_dist_[discrete].Insert(tri_discrete);
		}
	}
}


std::pair<double, int> MetricGraphObservable::Distance(const Triangle * t1, const Triangle * t2)
{
	return properties::TriangleWeightedDistance(triangulation_,t1,t2,dual_edge_length_);
}

std::pair<double, int> MetricGraphObservable::TriangulationDistance(const Vertex * v1, const Vertex * v2)
{
	return properties::VertexWeightedDistance(triangulation_,v1,v2,dual_edge_length_);
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
	stream << "}, tridistperlength -> {";
	for(int i=0,endi=tri_dist_per_length_.size();i<endi;++i)
	{
		stream << (i>0?",":"");
		tri_dist_per_length_[i].PrintTo(stream);
	}
	stream << "}, tridistperdist -> {";
	for(int i=0,endi=tri_dist_per_dist_.size();i<endi;++i)
	{
		stream << (i>0?",":"");
		tri_dist_per_dist_[i].PrintTo(stream);
	}
	stream << "}, tritimeperlength -> {";
	for(int i=0,endi=tri_time_per_length_.size();i<endi;++i)
	{
		stream << (i>0?",":"");
		tri_time_per_length_[i].PrintTo(stream);
	}
	stream << "}, trihopperlength -> {";
	for(int i=0,endi=tri_hop_per_length_.size();i<endi;++i)
	{
		stream << (i>0?",":"");
		tri_hop_per_length_[i].PrintTo(stream);
	}
	stream << "}, vertexdegree -> ";
	vertex_degree_.PrintTo(stream);
	stream << ", vertexdegreegeodesic -> ";
	vertex_degree_geodesic_.PrintTo(stream);
	stream << "}";
	return stream.str();
}

void MetricGraphObservable::MeasureDegreeOfVerticesAlongGeodesic(const Vertex * v1, const Vertex * v2)
{
	boost::array<std::vector<int>,2 > distances;
	properties::VertexDistanceList(triangulation_,v1,distances[0]);
	properties::VertexDistanceList(triangulation_,v2,distances[1]);
	
	std::vector<int> degree;
	properties::DegreeList(triangulation_,degree);

	for(int i=0,endi=triangulation_->NumberOfVertices();i<endi;++i)
	{
		if( distances[0][i] > 0 && distances[1][i] > 0 &&
			distances[0][i] + distances[1][i] == distances[0][v2->getId()] )
		{
			vertex_degree_geodesic_.Insert(degree[i]);
		}
		vertex_degree_.Insert(degree[i]);
	}
}