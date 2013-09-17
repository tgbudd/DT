#ifndef MODULI_OBSERVABLE_H
#define MODULI_OBSERVABLE_H

#include "Observable.h"
#include "Embedding.h"

void ToFundamentalDomain( std::pair< double, double > & x);

class ModuliObservable :
	public Observable
{
public:
	ModuliObservable(Embedding * const embedding);
	~ModuliObservable(void);

	void Measure();
	std::string OutputData() const;
private:
	std::vector<std::vector<int> > histogram_;
	double tau2_min_;
	double tau2_max_;
	int tau2_bins_;
	double tau1_min_;
	double tau1_max_;
	int tau1_bins_;
	int measurements_;
	Embedding * const embedding_;
};

#endif