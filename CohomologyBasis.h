#pragma once

#include <vector>

#include "decoration.h"
#include "triangulation.h"

class CohomologyBasis :
	public Decoration
{
public:
	CohomologyBasis() : triangulation_(NULL) {}
	CohomologyBasis(Triangulation * const triangulation);
	~CohomologyBasis(void);

	void setOmega(const Edge * const & edge, int x, int omega) {
		omega_[edge->getParent()->getId()][edge->getId()][x] = omega;
	}
	void setOmegaToMinusAdjacent(const Edge * const & edge) {
		for(int x=0;x<2;x++)
		{
			setOmega(edge,x,getOmega(edge->getAdjacent(),x));
		}
	}
	int getOmega(const Edge * const & edge, int x) const {
		return omega_[edge->getParent()->getId()][edge->getId()][x];
	}

	void Initialize();
	void Initialize(int width, int height);

	void UpdateAfterFlipMove(const Edge * const);
private:
	std::vector<boost::array<boost::array<int,2>,3> > omega_;		// two closed one-forms on the triangulation
	Triangulation * triangulation_;

};

