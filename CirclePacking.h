#ifndef CIRCLE_PACKING_H
#define CIRCLE_PACKING_H

#include "boost/array.hpp"

#include "Embedding.h"
#include "CohomologyBasis.h"
#include "LinearAlgebra.h"

class CirclePacking :
	public Embedding
{
public:
	CirclePacking(Triangulation * const triangulation, CohomologyBasis * const cohomologybasis);
	~CirclePacking() {}

	bool FindEdgeMeasure();
	bool FindRadii();
private:
//	void RadiiToAngles();
	double AngleSum(int vertexId);
	const Triangulation * const triangulation_;
	
	std::vector<int> degree_;
	std::vector<double> radius_;
	std::vector<boost::array<double,3> > angles_;

	int max_iterations_;
	double epsilon_, delta_;
};

#endif
