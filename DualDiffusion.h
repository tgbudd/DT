#ifndef DT_DUAL_DIFFUSION_H
#define DT_DUAL_DIFFUSION_H

#include <map>
#include <vector>

#include "Observable.h"
#include "LinearAlgebra.h"
#include "Triangulation.h"

class DualDiffusionMatrix : public linearalgebra::Matrix
{
public:
	DualDiffusionMatrix(const Triangulation * const triangulation, double walkprobability) : Matrix(triangulation->NumberOfTriangles()), triangulation_(triangulation), walkprobability_(walkprobability) {}

	void MultiplyVector(const std::vector<double> & from, std::vector<double> & to) const;
private:
	const Triangulation * const triangulation_;
	double walkprobability_;
};

class DualDiffusion :
	public Observable
{
public:
	DualDiffusion( Triangulation * const triangulation );
	~DualDiffusion(void);

	void Measure();
	std::string OutputData() const;
private:
	void DetermineDistance(Triangle * startTriangle);
	void DoDiffusion(Triangle * startTriangle, const linearalgebra::Matrix * const matrix );
	void DoMeasurementOnDistribution(int time);
	
	Triangulation * const triangulation_;

	std::vector<int> distance_;
	std::vector<double> probability_;
	std::vector<double> next_probability_;

	int measurements_;
	int samples_;
	int diffusion_steps_;
	double walkprobability_;

	int max_distance_;
	std::vector<std::vector<double> > distribution_;  // distribution[t][r] is fraction of random walks of length t ending at distance r
};

#endif
