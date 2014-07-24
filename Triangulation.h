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

typedef unsigned int TriangulationState;

class Triangulation
{
public:
	Triangulation(void);
	~Triangulation(void);

	void LoadRegularLattice(int width, int height);
	void LoadFromAdjacencyList(const std::vector<boost::array<std::pair<int,int>,3 > > & adj);
	void LoadSphericalBySubdivision(int NumberOfTriangles);

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

	Triangle * const & getRandomTriangle() const;
	Edge * const & getRandomEdge() const;
	Vertex * const & getRandomVertex() const;


	int RandomInteger(int min, int max) const;
	bool SucceedWithProbability(double probability) const;
	double RandomReal() const;
	double RandomReal(double min, double max) const;
	double RandomNormal(double mean, double sigma) const;
	void RandomSample(int min, int max, int n, std::vector<int> & sample) const;
	void setDominantMatter(DominantMatter * const & dominantmatter);
	void clearDominantMatter();
	void AddMatter(Matter * matter);
	void AddDecoration(Decoration * decoration);
	
	void DoSweep();
	void DoSweep(int);

	bool TryFlipMove();
	bool TryFlipMove(Edge * edge);

	bool TryCutMove(const boost::array< Edge *, 2> & edges, double combinatorialBoltzmann=1.0);
	void DoCutMove(const boost::array< Edge *, 2> & edges);

	void SeedRandom(unsigned int seed);

	std::string OutputData() const;

	bool IsState(const TriangulationState & state) const
	{
		return state == state_;
	}
	TriangulationState getState() const
	{
		return state_;
	}

	Triangle * NewTriangle();
	Triangle * NewTriangle(const Triangle & triangle);
	void DeleteTriangle(Triangle * triangle);

	void Clear();
	void DetermineVertices();

	double TotalCentralCharge() const;

	int CalculateGenus() const;

	void SetCustomSweepSize(int n);

	void RemoveBabyUniverses();
private:
	void LoadTetrahedron();
	void Subdivide(Triangle * triangle);
	void IncreaseState();
	bool CheckVertexNeighbourhood(const Vertex * const vertex) const;

	std::vector<Triangle*> triangles_;
	std::vector<Vertex* > vertices_;
	int n_triangles_, n_vertices_;

	mutable boost::mt19937 rng_;

	DominantMatter* dominantmatter_;   // matter may be present in the system that replaces the standard flip move
	std::list<Decoration* > decoration_;   // contains all additional structure that needs updating after a flip move
	std::list<Matter* > matter_;	// subset of decoration_ that has its own moves and boltzmann weights

	bool use_flipmove_;  // Defaults to true. May be set to false because some matter fields require a custom flip move.

	TriangulationState state_;

	bool use_custom_sweep_size_;	// Only used in special circumstances where we want to do 
	int custom_sweep_size_;			// more frequent observations (like when a video is produced).
};

#endif
