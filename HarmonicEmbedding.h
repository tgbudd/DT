#ifndef HARMONIC_EMBEDDING_H
#define HARMONIC_EMBEDDING_H

#include "Triangulation.h"
#include "Embedding.h"

class HarmonicEmbedding :
	public Embedding
{
public:
	HarmonicEmbedding(const Triangulation * const triangulation, CohomologyBasis * const cohomologybasis);
	~HarmonicEmbedding() {}
	bool FindEdgeMeasure();
	bool GetRadii(std::vector<double> & radii);

	enum RadiusDefinition {
		SQUARE_ROOT_OF_AREA,
		CIRCUMSCRIBED_CIRCLE,
		INSCRIBED_CIRCLE,
		AVERAGE_EDGE_LENGTH
	};
	void SetRadiusDefintion( RadiusDefinition radiusdefinition );
private:
	const Triangulation * const triangulation_;
	RadiusDefinition radiusdefinition_;
};

#endif
