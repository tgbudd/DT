#ifndef HYPERBOLIC_STRUCTURE_H
#define HYPERBOLIC_STRUCTURE_H

#include "Observable.h"
#include "CirclePattern.h"
#include "ThetaModel.h"
#include "utilities.h"

class HyperbolicStructure : public Observable
{
public:
	HyperbolicStructure(const Triangulation * const triangulation, ThetaModel * const thetamodel, CirclePattern * const circlepattern);
	~HyperbolicStructure();

	void Measure();
	std::string OutputData() const;
	void setMaxTheta(double max);
private:
	const Triangulation * const triangulation_;
	ThetaModel * const thetamodel_;
	CirclePattern * const circlepattern_;
	double HyperbolicLength(const std::list<const Edge*> & path) const;
	int NumberOfLegs(const std::list<const Edge*> & path) const;
	void MeasureLegShears(const std::list<const Edge*> & path);

	double maxtheta_;

	double sqlength_max_;
	int sqlength_bins_;
	std::vector<int> sqlength_histogram_;
	int measurements_;
	std::vector<int> legs_histogram_;
	double angle_surplus_max_;
	int angle_surplus_bins_;
	std::vector<int> angle_surplus_histogram_;
	double leg_shear_min_;
	double leg_shear_max_;
	int leg_shear_bins_;
	std::vector<int> leg_shear_histogram_;

	int num_log_;
	bool use_log_;
	bool log_fresh_;
	std::ostringstream log_;
};

#endif