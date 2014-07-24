#ifndef BALL_SIZE_DISTRIBUTION_H
#define BALL_SIZE_DISTRIBUTION_H


#include "Observable.h"
#include "Embedding.h"
#include "utilities.h"
#include "Histogram.h"

class BallSizeDistribution :
	public Observable
{
public:
	BallSizeDistribution(Embedding * const embedding, 
		const Triangulation * const triangulation,
		double absmaxlogradius, int logradiusbins, int samples);
	~BallSizeDistribution(void);

	void AddReferenceRadius(double epsilon);
	void AddVolumeFraction(double delta);
	void InitializeHistograms();
	void Measure();
	std::string OutputData() const;
private:
	Embedding * const embedding_;
	const Triangulation * const triangulation_;
	void ComputeDistanceHistogram(const std::vector<Vector2D> & x, const Vector2D & c, std::vector<std::vector<Histogram<int> > > & histograms);
	double DistanceSquared(Vector2D x) const;
	Histogram<double> log_dist_;
	std::vector<std::vector<Histogram<int> > > ball_size_;
	std::vector<std::vector<Histogram<int> > > unrooted_ball_size_;
	std::vector<double> ref_epsilons_;
	std::vector<double> deltas_;
	int samples_;

	std::pair<double,double> modulus_;
	std::vector<Vector2D> points_;
};

#endif