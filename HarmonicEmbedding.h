#pragma once

#include "triangulation.h"
#include "embedding.h"

class HarmonicEmbedding :
	public Embedding
{
public:
	HarmonicEmbedding(Triangulation * const triangulation, const CohomologyBasis * const cohomologybasis);
	~HarmonicEmbedding() {}
	bool FindEdgeMeasure();
private:
	const Triangulation * const triangulation_;
};

