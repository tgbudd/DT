#ifndef DUAL_COHOMOLOGY_BASIS_H
#define DUAL_COHOMOLOGY_BASIS_H

#include "CohomologyBasis.h"
#include "Triangulation.h"


class DualCohomologyBasis :
	public CohomologyBasis
{
public:
	DualCohomologyBasis() {}
	DualCohomologyBasis(const CohomologyBasis & cohomologybasis);
	DualCohomologyBasis(const Triangulation * const triangulation) : CohomologyBasis(triangulation) {}
	DualCohomologyBasis(const Triangulation * const triangulation, const std::vector<std::list<Edge*> > & generators, const std::vector<IntForm2D> & integrals );
	~DualCohomologyBasis(void);

	void Initialize(int width, int height);
	void Initialize();

	void UpdateAfterFlipMove(const Edge * const);

	bool CheckClosedness() const;
	IntForm2D IntegrateToParent(Edge * edge) const;

	void SetToDualOf(const CohomologyBasis & omega);		// change the CohomologyBasis into a DualCohomologyBasis
	void SetAccordingToGenerators( const std::vector<std::list<Edge*> > & generators, const std::vector<IntForm2D> & integrals);

	void Simplify(bool StayInSameClass=false);
};

#endif
