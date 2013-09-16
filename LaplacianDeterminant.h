#ifndef LAPLACIAN_DETERMINANT_H
#define LAPLACIAN_DETERMINANT_H

#include "Matter.h"
#include "Triangulation.h"
#include "ConjugateGradient.h"

class DualLaplacianMatrix : public Matrix
{
public:
	DualLaplacianMatrix(const Triangulation * const triangulation) : triangulation_(triangulation) {}

	void MultiplyVector(const std::vector<double> & from, std::vector<double> & to) const;
private:
	const Triangulation * const triangulation_;

};

class LaplacianDeterminant :
	public Matter
{
public:
	LaplacianDeterminant(const Triangulation * const triangulation, double centralcharge) 
		: triangulation_(triangulation), solution_(triangulation->NumberOfTriangles(),0.0),
		righthandside_(triangulation->NumberOfTriangles(),0.0), duallaplacianmatrix_(triangulation), centralcharge_(centralcharge) {}
	~LaplacianDeterminant(){}

	void DoSweep() {}
	void UpdateAfterFlipmove(const Edge * const edge) {}
	void Initialize() {}
	//bool IsFlipMoveAllowed(const Edge * const edge) { return true;}
	double BoltzmannChangeUnderFlipMove(const Edge * const ) const;
	void setCentralCharge(double c)
	{
		centralcharge_ = c;
	}
private:
	const Triangulation * const triangulation_;
	mutable std::vector<double> solution_;			// Make these mutable such that they can be called by the 
	mutable std::vector<double> righthandside_;		// const functions. They are only stored for speed purposes.
	const DualLaplacianMatrix duallaplacianmatrix_;
	double centralcharge_;
};

#endif
