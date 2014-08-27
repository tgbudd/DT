#ifndef METRIC_GRAPH_OBSERVABLE_H
#define METRIC_GRAPH_OBSERVABLE_H

#include "Observable.h"
#include "SimplexAttribute.h"
#include "Histogram.h"

class MetricGraphObservable : public Observable {
public:
	MetricGraphObservable(const Triangulation * triangulation, double MaxLength, int LengthBins, int samples);
	~MetricGraphObservable();

	void Measure();
	std::string OutputData() const;
private:
	void setRandomEdgeLengths();
	void MeasureTwoPointDistance();
	void MeasureDegreeOfVerticesAlongGeodesic(const Vertex * v1, const Vertex * v2);
	std::pair<double, int> Distance(const Triangle * t1, const Triangle * t2);
	std::pair<double, int> TriangulationDistance(const Vertex * v1, const Vertex * v2);

	EdgeAttribute<double> dual_edge_length_;
	TriangleAttribute<double> distance_;
	TriangleAttribute<int> num_edges_;
	VertexAttribute<double> v_distance_;
	VertexAttribute<int> v_num_edges_;
	const Triangulation * triangulation_;

	int samples_;
	double length_min_, length_max_;
	int length_bins_;
	int distance_is_zero_;
	std::vector<Histogram<int> > edges_per_length_;
	std::vector<Histogram<int> > tri_dist_per_length_;
	std::vector<Histogram<int> > dist_per_length_;
	std::vector<Histogram<int> > tri_dist_per_dist_;
	std::vector<Histogram<double> > tri_time_per_length_;
	std::vector<Histogram<int> > tri_hop_per_length_;
	Histogram<int> vertex_degree_geodesic_;
	Histogram<int> vertex_degree_;

	template<typename T>
	class SecondComparator {
	public:
		bool operator()(const std::pair<T,double> & t1, const std::pair<T,double> & t2 )
		{
			return t1.second > t2.second;
		}
	};
};

#endif