#ifndef CMINUSTWOBUILDER_H
#define CMINUSTWOBUILDER_H

#include "DominantMatter.h"
#include "Triangulation.h"

class LogFactorialTable 
{
public:
	LogFactorialTable(int max);
	double LogFactorial(int n) const;
	double LogBinomial(int n, int k) const;
	double LogBinomialDifference(int n, int k1, int k2) const;
	double LogNumberWalks(int u, int d) const;
	double LogNumberWalksMark(int u, int d, int q) const;
private:
	std::vector<double> logfactorial_;
	int max_;
	double veryNegative_;
};

class CMinusTwoBuilder :
	public DominantMatter
{
public:
	CMinusTwoBuilder(Triangulation * const triangulation, int Genus, int NumberOfTriangles, int boundaryLength = 0);
	~CMinusTwoBuilder();

	void DoSweep();
	
	double CentralCharge() const;
	std::string ConfigurationData() const;
	void getSpanningTree(std::vector<boost::array<bool,3> > & intree) const;
	void setRemoveBabyUniverses(bool remove);

	void getDiskBoundary(std::list<const Edge *> & boundary) const;
	void getPath(const Vertex * v0, const Vertex * v1, std::list<const Edge *> & boundary) const;
private:
	void RandomDiskTriangulation();
	void RandomBoundaryMatching();
	void RandomBoundaryMatchingDisk();
	void GenusZeroToGenusOneMatching();
	void ApplyBoundaryMatching(const std::vector<std::pair<int,int> > & matching);
	void ApplyBoundaryMatchingDisk();
	void InsertDiskAtBoundary();
	void getDiskBoundary(std::list<Edge *> & boundary) const;

	int ChooseRandomWithLogWeights(const std::vector<double> & w) const;

	Triangulation * const triangulation_;
	int n_triangles_;
	int genus_;
	bool remove_baby_universes_;

	bool make_disk_;
	int boundary_length_;
	std::vector<int> boundary_positions_;
	std::list<const Edge*> reverse_boundary_;

	std::vector<Edge *> boundary_;
	std::vector<int> tree_list_;
	std::vector<std::pair<int,int> > matching_;

	LogFactorialTable table_;
};

#endif
