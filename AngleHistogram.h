#ifndef ANGLE_HISTOGRAM_H
#define ANGLE_HISTOGRAM_H

#include "Triangulation.h"
#include "Observable.h"
#include "Embedding.h"

class AngleHistogram :
	public Observable
{
public:
	AngleHistogram(Embedding * const embedding, Triangulation * triangulation = NULL, int bins=180, double min=0.0, double max=PI);
	~AngleHistogram(void);

	void Measure();
	std::string OutputData() const;
private:
	const Triangulation * triangulation_;
	double CalculateAngle(int i, int j) const;
	double CalculateShear(const Edge * edge) const;
	double CalculateLogArea(int triangle) const;
	double CalculateLogLength(int triangle, int edge) const;
	Embedding * const embedding_;
	std::vector<int> histogram_;
	double min_;
	double max_;
	int bins_;
	int measurements_;

	bool measure_other_;
	std::vector<int> theta_histogram_;
	std::vector<int> shear_histogram_;
	double min_shear_;
	double max_shear_;
	int shear_bins_;
	std::vector<int> log_length_histogram_;
	std::vector<int> log_area_histogram_;
	double min_log_length_;
	double max_log_length_;
	int log_length_bins_;
	double min_log_area_;
	double max_log_area_;
	int log_area_bins_;
	std::vector<int> degree_histogram_;

	std::pair<double,double> modulus_;
};

#endif
