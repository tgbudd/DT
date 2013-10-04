#ifndef DT_DIFFUSION_H
#define DT_DIFFUSION_H

#include <map>
#include <vector>

#include "Observable.h"
#include "LinearAlgebra.h"
#include "Triangulation.h"

class DiffusionMatrix : public linearalgebra::Matrix
{
public:
	DiffusionMatrix(const Triangulation * const triangulation, const std::vector<int> & degree);

	void MultiplyVector(const std::vector<double> & from, std::vector<double> & to) const;
private:
	std::vector<std::map<int,double> > laplacianRules_;
};

class Diffusion :
	public Observable
{
public:
	Diffusion( Triangulation * const triangulation );
	~Diffusion(void);

	void Measure();
	std::string OutputData() const;
private:
	void DetermineDistance(Vertex * startVertex);
	void DetermineDegrees();
	void DoDiffusion(Vertex * startVertex, const linearalgebra::Matrix * const matrix );
	void DoMeasurementOnDistribution(int time);
	void MeasureDistanceDistribution();
	
	Triangulation * const triangulation_;

	std::vector<int> degree_;
	std::vector<int> distance_;
	std::vector<double> probability_;
	std::vector<double> next_probability_;

	int measurements_;
	int samples_;
	std::vector<int> diffusion_times_;

	int max_distance_;
	std::vector<std::vector<double> > distribution_;  // distribution[t][r] is fraction of random walks of length t ending at distance r
	std::vector<int> distance_distribution_;
};

#endif
