#pragma once

#include <vector>
#include <list>
#include <time.h>

#include "boost/random/mersenne_twister.hpp"
#include "boost/random/uniform_int_distribution.hpp"
#include "boost/random/uniform_real_distribution.hpp"
#include "boost/random/discrete_distribution.hpp"
#include "boost/foreach.hpp"

class Triangulation;
class Triangle;
class Vertex;
class Matter;

#include "utilities.h"
#include "Triangle.h"
#include "Vertex.h"
#include "Matter.h"


class Triangulation
{
public:
	Triangulation(void);
	~Triangulation(void);

	void LoadRegularLattice(int width, int height);

	int NumberOfTriangles() const
	{
		return n_triangles_;
	}
	
	int NumberOfVertices() const
	{
		return n_vertices_;
	}

	Triangle * const & getTriangle(int id) const
	{
		return triangles_[id];
	}
	
	Vertex * const & getVertex(int id) const
	{
		return vertices_[id];
	}
	Triangle * const & getRandomTriangle()
	{
		return triangles_[RandomInteger(0,n_triangles_-1)];
	}
	Edge * const & getRandomEdge()
	{
		return triangles_[RandomInteger(0,n_triangles_-1)]->getEdge(RandomInteger(0,2));
	}


	int RandomInteger(int min, int max)
	{
		boost::random::uniform_int_distribution<> distribution(min, max);
		return distribution(rng_);
	}
	bool SucceedWithProbability(double probability)
	{
		if( probability < 0.0 )
			return false;
		if( probability > 1.0 )
			return true;
		double probabilities[] = {probability,1.0-probability};
		boost::random::discrete_distribution<> distribution (probabilities);
		return distribution(rng_) == 0;
	}

	void AddMatter(Matter * matter)
	{
		decoration_.push_back(matter);
		matter_.push_back(matter);
		if( matter->ReplacesFlipMove() )
		{
			use_flipmove_ = false;
		}
	}
	void AddDecoration(Decoration * decoration)
	{
		decoration_.push_back(decoration);
	}
	
	void DoSweep();
	void DoSweep(int);

	bool TryFlipMove();

private:
	void DetermineVertices();

	std::vector<Triangle*> triangles_;
	std::vector<Vertex* > vertices_;
	int n_triangles_, n_vertices_;

	boost::random::mt19937 rng_;

	std::list<Decoration* > decoration_;   // contains all additional structure that needs updating after a flip move
	std::list<Matter* > matter_;	// subset of decoration_ that has its own moves and boltzmann weights

	bool use_flipmove_;  // Defaults to true. May be set to false because some matter fields require a custom flip move.
};

