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
	std::pair<double, int> Distance(const Triangle * t1, const Triangle * t2);

	EdgeAttribute<double> dual_edge_length_;
	TriangleAttribute<double> distance_;
	TriangleAttribute<int> num_edges_;
	const Triangulation * triangulation_;

	int samples_;
	double length_min_, length_max_;
	int length_bins_;
	int distance_is_zero_;
	std::vector<Histogram<int> > edges_per_length_;
	std::vector<Histogram<int> > dist_per_length_;

	class TriangleDistanceComparator {
	public:
		bool operator()(const std::pair<const Triangle *,double> & t1, const std::pair<const Triangle *,double> & t2 )
		{
			return t1.second > t2.second;
		}
	};
};

#endif