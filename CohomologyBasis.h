#ifndef COHOMOLOGY_BASIS_H
#define COHOMOLOGY_BASIS_H

#include <vector>

#include "Decoration.h"
#include "Triangulation.h"

typedef boost::array<int,2> IntForm2D;

inline IntForm2D AddForms( const IntForm2D & f1, const IntForm2D & f2 )
{
	IntForm2D sum(f1);
	sum[0] += f2[0];
	sum[1] += f2[1];
	return sum;
}

inline IntForm2D SubtractForms( const IntForm2D & f1, const IntForm2D & f2 )
{
	IntForm2D diff(f1);
	diff[0] -= f2[0];
	diff[1] -= f2[1];
	return diff;
}
inline IntForm2D NegateForm( const IntForm2D & f )
{
	IntForm2D neg;
	neg[0] = -f[0];
	neg[1] = -f[1];
	return neg;
}
inline bool FormsIndependent( const IntForm2D & f1, const IntForm2D & f2 )
{
	return f1[0] * f2[1] - f1[1] * f2[0] != 0;
}
inline bool FormIsZero( const IntForm2D & f )
{
	return f[0] == 0 && f[1] == 0;
}

class DualCohomologyBasis;

class CohomologyBasis :
	public Decoration
{
public:
	CohomologyBasis() : triangulation_(NULL) {}
	CohomologyBasis(const Triangulation * const triangulation);
	CohomologyBasis(const CohomologyBasis & cohomologybasis) : omega_(cohomologybasis.omega_), triangulation_(cohomologybasis.triangulation_) {}
	CohomologyBasis(const Triangulation * const triangulation, const DualCohomologyBasis & dualcohomologybasis);
	~CohomologyBasis(void);



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
			setOmega(edge,x,-getOmega(edge->getAdjacent(),x));
		}
	}
	int getOmega(int triangle, int edge, int x) const {
		return omega_[triangle][edge][x];
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

	bool CheckClosedness() const;

	IntForm2D Integrate(const std::list<Edge *> & path) const;

	void Simplify(bool StayInSameClass = false);
	void SetToDualOf(const DualCohomologyBasis & dualOmega);	// change the DualCohomologyBasis to an equivalent CohomologyBasis

protected:
	std::vector<boost::array<IntForm2D,3> > omega_;		// two closed one-forms on the triangulation
	void ClearOmega();
	const Triangulation * triangulation_;
private:

	
};

#endif
