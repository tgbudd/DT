#pragma once

#include <vector>
#include <list>
#include <time.h>

#include "boost/random/mersenne_twister.hpp"
#include "boost/random/uniform_int_distribution.hpp"
#include "boost/random/uniform_real_distribution.hpp"
#include "boost/random/discrete_distribution.hpp"

class Triangulation;
class Triangle;
class Vertex;
class Matter;

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
		matter_.push_back(matter);
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

	std::list<Matter* > matter_;
};

