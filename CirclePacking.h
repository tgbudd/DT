#ifndef CIRCLE_PACKING_H
#define CIRCLE_PACKING_H

#include "boost/array.hpp"

#include "Embedding.h"
#include "CohomologyBasis.h"
#include "LinearAlgebra.h"
#include "BabyUniverseDetector.h"
#include "LaplacianMatrix.h"

class CirclePacking :
	public Embedding
{
public:
	CirclePacking(Triangulation * const triangulation, CohomologyBasis * const cohomologybasis);
	~CirclePacking() {}

	bool FindEdgeMeasure();
	void CalculateEdgeMeasure(std::vector<boost::array<double,3> > & measure);
	bool FindRadii();
	void RadiiToAngles();

	bool GetRadii(std::vector<double> & radii);
private:
//	void RadiiToAngles();
	double AngleSum(int vertexId);
	double Angle(Edge * edge);
	const Triangulation * const triangulation_;

	BabyUniverseDetector babyuniversedetector_;
	CohomologyBasis * cohomologybasis_;

	std::vector<int> degree_;
	std::vector<double> radius_;
	std::vector<boost::array<double,3> > angles_;
	

	int max_iterations_;
	double epsilon_, delta_;
};

#endif
