#pragma once

#include <vector>

#include "triangulation.h"
#include "matter.h"
#include "CohomologyBasis.h"

class ThetaModel :
	public Matter
{
public:
	//ThetaModel() : triangulation_(NULL) {}
	ThetaModel(Triangulation * const triangulation, const CohomologyBasis * const cohomologybasis);
	~ThetaModel(void);

	bool ReplacesFlipMove() {
		return true;
	}

	void Initialize();
	void UpdateAfterFlipMove(const Edge * const) {}
	void DoSweep();

	void setTheta(Edge * const & edge, double theta) {
		theta_[edge->getParent()->getId()][edge->getId()] = theta;
		theta_[edge->getOpposite()->getParent()->getId()][edge->getOpposite()->getId()] = theta;
	}

	double getTheta(Edge * const & edge) {
		return theta_[edge->getParent()->getId()][edge->getId()]; 
	}

private:
	std::vector<boost::array<double,3> > theta_;
		
	Triangulation * const triangulation_;
	const CohomologyBasis * const cohomologybasis_;


	bool TestCutCondition();

};

