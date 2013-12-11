#ifndef CONFORMAL_DISTRIBUTION_H
#define CONFORMAL_DISTRIBUTION_H

#include "Observable.h"
#include "Embedding.h"

class ConformalDistribution :
	public Observable
{
public:
	ConformalDistribution(Embedding * const embedding);
	~ConformalDistribution(void);

	void Measure();
	std::string OutputData() const;
private:
	int measurements_;
	Embedding * const embedding_;

	std::vector<int> histogram_;
	double minlograd_, maxlograd_;
	int bins_;
};

#endif