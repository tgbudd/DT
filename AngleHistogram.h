#ifndef ANGLE_HISTOGRAM_H
#define ANGLE_HISTOGRAM_H

#include "Observable.h"
#include "Embedding.h"

class AngleHistogram :
	public Observable
{
public:
	AngleHistogram(const Embedding * const embedding, int bins=180, double min=0.0, double max=PI);
	~AngleHistogram(void);

	void Measure();
	std::string OutputData() const;
private:
	double CalculateAngle(int i, int j) const;
	const Embedding * const embedding_;
	std::vector<int> histogram_;
	double min_;
	double max_;
	int bins_;
	int measurements_;
};

#endif
