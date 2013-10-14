#ifndef CMINUSTWOBUILDER_H
#define CMINUSTWOBUILDER_H

#include "DominantMatter.h"
#include "Triangulation.h"

class CMinusTwoBuilder :
	public DominantMatter
{
public:
	CMinusTwoBuilder(Triangulation * const triangulation, int Genus, int NumberOfTriangles);
	~CMinusTwoBuilder();

	void DoSweep();
	
	double CentralCharge() const;
	std::string ConfigurationData() const;
	void getSpanningTree(std::vector<boost::array<bool,3> > & intree) const;
private:
	void RandomDiskTriangulation();
	void RandomBoundaryMatching();
	void MatchingToGenusOneMatching();
	void ApplyBoundaryMatching(const std::vector<std::pair<int,int> > & matching);

	Triangulation * const triangulation_;
	int n_triangles_;
	int genus_;

	std::vector<Edge *> boundary_;
	std::vector<int> tree_list_;
	std::vector<std::pair<int,int> > matching_;
};

#endif
