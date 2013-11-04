#ifndef HARMONIC_DIFFUSION_H
#define HARMONIC_DIFFUSION_H

#include "observable.h"
#include "Triangulation.h"
#include "Edge.h"
#include "Embedding.h"
#include "CohomologyBasis.h"

class HarmonicDiffusion :
	public Observable
{
public:
	HarmonicDiffusion(const Triangulation * const triangulation, Embedding * const embedding, CohomologyBasis * const cohomologybasis);
	~HarmonicDiffusion(void);

	void Measure();
	std::string OutputData() const;
private:
	void DoRandomWalk();
	Edge * RandomEdge(Vertex * v);
	void DoGraphDistanceMeasurement();

	const Triangulation * const triangulation_;
	Embedding * const embedding_;
	CohomologyBasis * const cohomologybasis_;
	std::pair<double,double> moduli_;

	std::vector<double> average_distance_;
	std::vector<double> square_distance_;

	std::vector<std::map<IntForm2D,int> > graph_distance_;
	std::vector<int> graph_distance_distribution_;
	std::vector<double> average_distance_per_graph_distance_;
	std::vector<double> square_distance_per_graph_distance_;

	int samples_;
	int graph_distance_samples_;
	int measurements_;
	int walk_time_;
	int max_graph_distance_;
};

#endif