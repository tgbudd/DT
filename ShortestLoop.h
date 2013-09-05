#pragma once

#include <list>

#include "observable.h"
#include "triangulation.h"
#include "CohomologyBasis.h"


class ShortestLoop :
	public Observable
{
public:
	ShortestLoop(const Triangulation * const triangulation, const CohomologyBasis * const cohomologybasis) : triangulation_(triangulation), cohomologybasis_(cohomologybasis) {}
	~ShortestLoop(void)  {}

	void FindShortestLoop(int maxLength = -1);
	void FindGenerators();
	std::list<Edge*> FindShortestLoop(Vertex * startVertex, int maxLength = -1 );
	std::list<Edge*> FindShortestLoop(const std::list<Vertex *> & startVertices, int maxLength = -1 );
	std::list<Edge*> FindShortestLoop(const std::list<Vertex *> & startVertices, const IntForm2D & notMultipleOf, int maxLength = -1 );
	std::list<Edge*> FindShortestLoop(Vertex * startVertex, const IntForm2D & notMultipleOf, int maxLength = -1 );

	bool CheckPathIsLoop(const std::list<Edge*> & path) const;

	void Measure() {}

	const std::vector<std::list<Edge*> > & getGenerators() const
	{
		return generators_;
	}
	const std::vector<IntForm2D> & getGeneratorIntegrals() const
	{
		return generator_integrals_;
	}
	const std::list<Edge*> & getShortestLoop() const
	{
		return shortestloop_;
	}
private:
	const Triangulation * const triangulation_;
	const CohomologyBasis * const cohomologybasis_;

	std::vector<std::list<Edge*> > generators_;
	std::vector<IntForm2D > generator_integrals_;

	std::list<Edge*> shortestloop_;

	struct VertexNode {
		int distance;
		Edge * parent;
		IntForm2D integral;
	};

	std::vector<VertexNode> visit;

	void TrimLoop( std::list<Edge*> & path );
};

