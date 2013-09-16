#ifndef THETA_HISTOGRAM_H
#define THETA_HISTOGRAM_H

#include "Observable.h"
#include "ThetaModel.h"

class ThetaHistogram :
	public Observable
{
public:
	ThetaHistogram(const ThetaModel * const thetamodel, int bins=180, double min=0.0, double max=PI);
	~ThetaHistogram(void);

	void Measure();
	std::string OutputData() const;
private:
	const ThetaModel * const thetamodel_;
	std::vector<int> histogram_;
	double min_;
	double max_;
	int bins_;
	int measurements_;
};

#endif
