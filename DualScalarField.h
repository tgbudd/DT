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
	void DoSweep();
	double getField(Triangle * triangle) const;
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
