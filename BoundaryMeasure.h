#ifndef BOUNDARY_MEASURE_H
#define BOUNDARY_MEASURE_H

#include <Eigen/Sparse>

#include "Triangulation.h"
#include "Observable.h"
#include "BabyUniverseDetector.h"
#include "Triangle.h"
#include "Edge.h"

class BoundaryMeasure : public Observable {
public:
	BoundaryMeasure( Triangulation * triangulation );
	void Measure();
	std::string OutputData() const;

	virtual void DetermineBoundary() = 0;
protected:
	std::list<const Edge *> boundary_;
private:
	Triangulation * triangulation_;
	BabyUniverseDetector babyuniversedetector_;

	std::vector<Triangle *> triangles_;

	Eigen::MatrixXd invLaplacian_;
	bool DetermineInverseLaplacian();
	void UpdateHistogram();

	int samples_;
	int max_delta_;
	double minlogmeasure_, maxlogmeasure_;
	int measure_bins_;
	std::vector<std::vector<int> > measure_histogram_;
	std::vector<int> measurements_;
};

#endif