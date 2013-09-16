#ifndef HARMONIC_EMBEDDING_H
#define HARMONIC_EMBEDDING_H

#include "Triangulation.h"
#include "Embedding.h"

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

#endif
