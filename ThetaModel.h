#pragma once

#include <vector>

#include "triangulation.h"
#include "DominantMatter.h"
#include "DualCohomologyBasis.h"

class ThetaModel :
	public DominantMatter
{
public:
	//ThetaModel() : triangulation_(NULL) {}
	ThetaModel(Triangulation * const triangulation, const DualCohomologyBasis * const dualcohomologybasis, int PiInUnits = 1800);
	~ThetaModel(void);

	void Initialize();
	void DoSweep();

	void setTheta(Edge * const & edge, int theta) {
		BOOST_ASSERT( theta > 0 && theta < 2 * pi_in_units_ );
		theta_[edge->getParent()->getId()][edge->getId()] = theta;
		theta_[edge->getAdjacent()->getParent()->getId()][edge->getAdjacent()->getId()] = theta;
	}
	void addToTheta(Edge * const & edge, int theta) {
		BOOST_ASSERT( theta_[edge->getParent()->getId()][edge->getId()] == theta_[edge->getAdjacent()->getParent()->getId()][edge->getAdjacent()->getId()] );
		int firsttheta = (theta_[edge->getParent()->getId()][edge->getId()] += theta);
		BOOST_ASSERT( firsttheta > 0 && firsttheta < 2 * pi_in_units_ );
		firsttheta = (theta_[edge->getAdjacent()->getParent()->getId()][edge->getAdjacent()->getId()] += theta);
		BOOST_ASSERT( firsttheta > 0 && firsttheta < 2 * pi_in_units_ );
	}
	int getTheta(Edge * const & edge) const {
		return theta_[edge->getParent()->getId()][edge->getId()]; 
	}
	int getTheta(int triangleId, int edgeId) const {
		return theta_[triangleId][edgeId]; 
	}
	double getRealTheta(Edge * const & edge ) const {
		return PI * (double)(getTheta(edge))/pi_in_units_;
	}
	double getRealTheta(int triangleId, int edgeId) const {
		return PI * (double)(getTheta(triangleId,edgeId))/pi_in_units_;
	}
private:
	std::vector<boost::array<int,3> > theta_;
	const int pi_in_units_;	// use integer angles where pi_in_units corresponds to an angle of pi

	Triangulation * const triangulation_;
	const DualCohomologyBasis * const dualcohomologybasis_;

	bool TryThetaMove();
	bool TryThetaMove(Edge * edge);

	bool TestCutCondition(Edge *, Edge *, const IntForm2D &, int theta) const;

	Edge * previousEdge(Edge * edge);  // functions of which pointers are used in the thetamove
	Edge * nextEdge(Edge * edge);

	bool TestVertexSum(Vertex *);
};

