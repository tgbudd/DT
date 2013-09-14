#pragma once

#include "observable.h"
#include "CohomologyBasis.h"

class BabyUniverseDistribution :
	public Observable
{
public:
	BabyUniverseDistribution(const Triangulation * const triangulation, const CohomologyBasis * const cohomologybasis, int MinbuNeckSize) 
		: triangulation_(triangulation), cohomologybasis_(cohomologybasis), minbu_necksize_(MinbuNeckSize), measurements_(0) {}
	~BabyUniverseDistribution() {}

	void FindMinbuNecksOfLength3(std::list<std::list<const Edge *> > & paths);
	void FindMinbuNecksOfLength3(const Edge * edge, std::list<std::list<const Edge *> > & paths);

	void FindMinbuSizes(int necklength, std::vector<int> & sizes);
	void FindMinbuNecks(int necklength, std::list<std::list<const Edge *> > & paths);

	void Measure() 
	{
		FindMinbuSizes(minbu_necksize_,sizes_);
		measurements_++;
	}

	std::string OutputData() const;
private:
	struct VertexNode {
		Vertex * vertex;
		IntForm2D integral;
		bool operator<(const VertexNode & v) const {
			return vertex < v.vertex || (vertex == v.vertex && (integral[0] < v.integral[0] || (integral[0] == v.integral[0] && integral[1] < v.integral[1]) ) );
		}
		bool operator>(const VertexNode & v) const {
			return vertex > v.vertex || (vertex == v.vertex && (integral[0] > v.integral[0] || (integral[0] == v.integral[0] && integral[1] > v.integral[1]) ) );
		}
	};

	int VolumeEnclosed(const std::list<const Edge *> & boundary) const;

	const Triangulation * const triangulation_;
	const CohomologyBasis * const cohomologybasis_;
	std::vector<int> sizes_;
	int minbu_necksize_;
	int measurements_;
};

