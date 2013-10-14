#ifndef DISTANCE_DISTRIBUTION_H
#define DISTANCE_DISTRIBUTION_H

#include <map>

#include "Triangulation.h"
#include "observable.h"
#include "CohomologyBasis.h"

class DistanceDistribution :
	public Observable
{
public:
	DistanceDistribution(const Triangulation * const triangulation);
	~DistanceDistribution();

	void UseUniversalCover(CohomologyBasis * cohomologybasis);

	void Measure();
	std::string OutputData() const;

private:
	void MeasureDistanceInUniversalCover();

	const Triangulation * const triangulation_;
	CohomologyBasis * cohomologybasis_;

	std::vector<std::map<IntForm2D,int> > distance_;

	std::vector<int> distance_distribution_;
	int samples_;
	int measurements_;
	int max_distance_;
	bool in_universal_cover_;
	int genus_;
};

#endif