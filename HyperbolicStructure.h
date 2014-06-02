#ifndef HYPERBOLIC_STRUCTURE_H
#define HYPERBOLIC_STRUCTURE_H

#include <Eigen/dense>
#include <complex>

#include "Observable.h"
#include "CirclePattern.h"
#include "ThetaModel.h"
#include "utilities.h"
#include "Histogram.h"


class HyperbolicStructure : public Observable
{
public:
	HyperbolicStructure(const Triangulation * const triangulation, ThetaModel * const thetamodel, CirclePattern * const circlepattern, DualCohomologyBasis * const dualcohomologybasis);
	~HyperbolicStructure();

	void Measure();
	std::string OutputData() const;
	void setMaxTheta(double max);
	void setMaxLength(double length);
	int FindShortHyperbolicCurves();

	struct Curve {
		double length;
		std::list<boost::array<Vector2D,2> > segments;
	};

	void RetrieveCurves(std::list<Curve> & curves) const;
private:
	typedef Eigen::Matrix2d SL2Mat;
	typedef std::complex<double> Complex;

	const Triangulation * const triangulation_;
	ThetaModel * const thetamodel_;
	CirclePattern * const circlepattern_;
	DualCohomologyBasis * const dualcohomologybasis_;

	void MeasureCuspSize(Complex modulus);
	double CuspSize(const Vertex * vertex, Complex modulus) const;

	SL2Mat Holonomy(const std::list<const Edge*> & path, bool startWithTurnRight = false) const;
	void GeodesicCoordinates(const std::list<const Edge*> & path, std::list<boost::array<Vector2D,2> > & coor, Complex modulus) const;
	double HyperbolicLength(const std::list<const Edge*> & path) const;
	int NumberOfLegs(const std::list<const Edge*> & path) const;
	void MeasureLegShears(const std::list<const Edge*> & path, double length);
	static boost::array<double,2> FixedPoints(const SL2Mat & mat); 
	static Complex ToComplex(Vector2D v);
	static Vector2D ToVector2D(Complex x);
	static Complex MapToTriangle(Complex x, const boost::array<Complex,3> & v);
	static Complex PoincareDiskToKleinDisk(Complex x, const boost::array<Complex,3> & v);

	struct PathInfo {
		SL2Mat holonomy;
		IntForm2D omega;
		bool has_turned_left, has_turned_right;
	};
	void SearchPath(std::list<const Edge *> & path, PathInfo & info, double maxcoshlength);
	SL2Mat ShearMatrix(const Edge * edge) const;
	ReusableFlag visited_;
	std::list<std::pair<double,std::list<const Edge *> > > paths_;
	SL2Mat turn_right_mat_;

	double maxtheta_;
	double maxhyplength_;

	Histogram<double> sqlength_;
	ExtendableHistogram<int> legs_;
	Histogram<double> angle_surplus_;
	Histogram<double> leg_shear_;
	Histogram<double> log_sqlength_shear_ratio_;
	Histogram<double> cusp_size_;

	int num_log_;
	bool use_log_;
	bool log_fresh_;
	std::ostringstream log_;
};

#endif