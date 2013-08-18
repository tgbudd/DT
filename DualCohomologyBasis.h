#pragma once

#include "CohomologyBasis.h"
#include "triangulation.h"

class DualCohomologyBasis :
	public CohomologyBasis
{
public:
	DualCohomologyBasis() {}
	DualCohomologyBasis(Triangulation * const triangulation) : CohomologyBasis(triangulation) {}
	~DualCohomologyBasis(void);

	void Initialize(int width, int height);

	void UpdateAfterFlipMove(const Edge * const);

private:
	void Dualize();		// change the CohomologyBasis into a DualCohomologyBasis
};

