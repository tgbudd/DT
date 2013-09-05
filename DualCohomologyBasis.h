#pragma once

#include "CohomologyBasis.h"
#include "triangulation.h"


class DualCohomologyBasis :
	public CohomologyBasis
{
public:
	DualCohomologyBasis() {}
	DualCohomologyBasis(const CohomologyBasis & cohomologybasis);
	DualCohomologyBasis(Triangulation * const triangulation) : CohomologyBasis(triangulation) {}
	DualCohomologyBasis(Triangulation * const triangulation, const std::vector<std::list<Edge*> > & generators, const std::vector<IntForm2D> & integrals );
	~DualCohomologyBasis(void);

	void Initialize(int width, int height);

	void UpdateAfterFlipMove(const Edge * const);

	bool CheckClosedness() const;
	IntForm2D IntegrateToParent(Edge * edge) const;

	void SetToDualOf(const CohomologyBasis & omega);		// change the CohomologyBasis into a DualCohomologyBasis
};

