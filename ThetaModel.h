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
	ThetaModel(Triangulation * const triangulation, const DualCohomologyBasis * const dualcohomologybasis);
	~ThetaModel(void);

	void Initialize();
	void DoSweep();

	void setTheta(Edge * const & edge, double theta) {
		BOOST_ASSERT( theta > 0.0 && theta < 2.0 * PI );
		theta_[edge->getParent()->getId()][edge->getId()] = theta;
		theta_[edge->getAdjacent()->getParent()->getId()][edge->getAdjacent()->getId()] = theta;
	}
	void addToTheta(Edge * const & edge, double theta) {
		BOOST_ASSERT( 1.0e-7 > fabs(theta_[edge->getParent()->getId()][edge->getId()] - theta_[edge->getAdjacent()->getParent()->getId()][edge->getAdjacent()->getId()]) );
		double firsttheta = (theta_[edge->getParent()->getId()][edge->getId()] += theta);
		BOOST_ASSERT( firsttheta > 0.0 && firsttheta < 2.0 * PI );
		firsttheta = (theta_[edge->getAdjacent()->getParent()->getId()][edge->getAdjacent()->getId()] += theta);
		BOOST_ASSERT( firsttheta > 0.0 && firsttheta < 2.0 * PI );
	}
	double getTheta(Edge * const & edge) const {
		return theta_[edge->getParent()->getId()][edge->getId()]; 
	}

private:
	std::vector<boost::array<double,3> > theta_;
		
	Triangulation * const triangulation_;
	const DualCohomologyBasis * const dualcohomologybasis_;

	bool TryThetaMove();
	bool TryThetaMove(Edge * edge);

	bool TestCutCondition(Edge *, Edge *, const DualCohomologyBasis::IntForm2D &, double theta) const;

	Edge * previousEdge(Edge * edge);  // functions of which pointers are used in the thetamove
	Edge * nextEdge(Edge * edge);

	bool TestVertexSum(Vertex *);
};

