#ifndef DT_TRIANGULATION_H
#define DT_TRIANGULATION_H

#include <vector>
#include <list>
#include <time.h>

#include "boost/random/mersenne_twister.hpp"
#include "boost/foreach.hpp"

#include "utilities.h"

class Triangle;
class Vertex;
class Matter;
class DominantMatter;
class Decoration;
class Edge;

class Triangulation
{
public:
	Triangulation(void);
	~Triangulation(void);

	void LoadRegularLattice(int width, int height);
	void LoadFromAdjacencyList(const std::vector<boost::array<std::pair<int,int>,3 > > & adj);

	int NumberOfTriangles() const
	{
		return n_triangles_;
	}
	
	int NumberOfVertices() const
	{
		return n_vertices_;
	}

	Triangle * const & getTriangle(int id) const;
	
	Vertex * const & getVertex(int id) const;

	Triangle * const & getRandomTriangle();
	Edge * const & getRandomEdge();


	int RandomInteger(int min, int max);
	bool SucceedWithProbability(double probability);
	double RandomReal();
	double RandomReal(double min, double max);
	void setDominantMatter(DominantMatter * const & dominantmatter);
	void AddMatter(Matter * matter);
	void AddDecoration(Decoration * decoration);
	
	void DoSweep();
	void DoSweep(int);

	bool TryFlipMove();
	bool TryFlipMove(Edge * edge);

	void SeedRandom(unsigned int seed);

	std::string OutputData() const;

private:
	void DetermineVertices();

	std::vector<Triangle*> triangles_;
	std::vector<Vertex* > vertices_;
	int n_triangles_, n_vertices_;

	boost::mt19937 rng_;

	DominantMatter* dominantmatter_;   // matter may be present in the system that replaces the standard flip move
	std::list<Decoration* > decoration_;   // contains all additional structure that needs updating after a flip move
	std::list<Matter* > matter_;	// subset of decoration_ that has its own moves and boltzmann weights

	bool use_flipmove_;  // Defaults to true. May be set to false because some matter fields require a custom flip move.
};

#endif
