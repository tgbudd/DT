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

	typedef boost::array<int,2> IntForm2D;

	void setOmega(const Edge * const & edge, int x, int omega) {
		omega_[edge->getParent()->getId()][edge->getId()][x] = omega;
	}
	void setOmega(const Edge * const & edge, const IntForm2D & omega) {
		omega_[edge->getParent()->getId()][edge->getId()] = omega;
	}
	void addToOmega(const Edge * const & edge, int x, int omega) {
		omega_[edge->getParent()->getId()][edge->getId()][x] += omega;
	}
	void addToOmega(const Edge * const & edge, const IntForm2D & omega) {
		omega_[edge->getParent()->getId()][edge->getId()][0] += omega[0];
		omega_[edge->getParent()->getId()][edge->getId()][1] += omega[1];
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
	const boost::array<int,2> & getOmega(const Edge * const & edge) const {
		return omega_[edge->getParent()->getId()][edge->getId()];
	}

	void Initialize();
	void Initialize(int width, int height);

	void UpdateAfterFlipMove(const Edge * const);
protected:
	std::vector<boost::array<IntForm2D,3> > omega_;		// two closed one-forms on the triangulation
	Triangulation * triangulation_;

};

