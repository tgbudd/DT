#ifndef SPANNING_TREE_H
#define SPANNING_TREE_H

#include "Triangulation.h"
#include "DominantMatter.h"
#include "Edge.h"

class SpanningTree :
	public DominantMatter
{
public:
	SpanningTree(Triangulation * const triangulation);
	~SpanningTree();

	void Initialize();

	void DoSweep();
	bool InSpanningTree(int triangle,int edge) const;
	bool InSpanningTree(Edge * edge) const;
	void setInSpanningTree(int triangle, int edge, bool ist);
	void setInSpanningTree(Edge * edge, bool ist);

	double CentralCharge() const;
	std::string ConfigurationData() const;
private:
	int SpanningDegree(const Vertex * vertex) const;
	void UpdateSpanningDegree(const Vertex * vertex);
	void UpdateSpanningDegree();
	Edge * RandomOtherSpanningEdge(Edge * edge);

	bool TryFlipMove(Edge * edge);
	bool TrySpanningMove(Edge * edge);

	Triangulation * const triangulation_;
	std::vector<boost::array<bool,3> > in_spanning_tree_;
	std::vector<int> spanning_degree_;  // number of spanning tree edges incident at vertex
};

#endif