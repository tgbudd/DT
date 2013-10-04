#ifndef DUAL_SCALAR_FIELD_H
#define DUAL_SCALAR_FIELD_H

#include "Matter.h"
#include "Edge.h"

class DualScalarField :
	public Matter
{
public:
	DualScalarField() : triangulation_(NULL) {}
	DualScalarField(Triangulation * const triangulation);
	DualScalarField(Triangulation * const triangulation, double massSquared);
	~DualScalarField();

	void Initialize();
	double BoltzmannChangeUnderFlipMove(const Edge * const ) const;
	double BoltzmannChangeUnderGeneralMove(const std::vector<boost::array<Triangle *,2> > & toBeDeleted, const std::vector<boost::array<Triangle *,2> > & toBeAdded ) const;
	void DoSweep();
	double getField(Triangle * triangle) const;

	double CentralCharge() const;
	std::string ConfigurationData() const;
private:
	double getField(int triangle) const;
	void setField(Triangle * triangle, const double & field);
	void setField(int triangle, const double & field);	void DoMove();
	Triangulation * triangulation_;

	bool massive_;
	double mass_squared_;
	std::vector<double> field_;
};

#endif
