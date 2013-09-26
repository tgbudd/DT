#pragma once

#include "Triangulation.h"
#include "observable.h"

class LaplacianSpectrum :
	public Observable
{
public:
	LaplacianSpectrum(const Triangulation * const triangulation);
	~LaplacianSpectrum(void);

	void Measure();
	std::string OutputData() const;
private:
	const Triangulation * const triangulation_;

	std::vector<int> eigenvalue_histogram_;
	int ev_bins_;
	double ev_max_;
	double ev_min_;

	int measurements_;
};

