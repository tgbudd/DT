#ifndef BABY_UNIVERSE_DISTRIBUTION_H
#define BABY_UNIVERSE_DISTRIBUTION_H

#include "Observable.h"
#include "CohomologyBasis.h"
#include "DualCohomologyBasis.h"

class BabyUniverseDistribution :
	public Observable
{
public:
	BabyUniverseDistribution(const Triangulation * const triangulation, const CohomologyBasis * const cohomologybasis, int MinbuNeckSize) 
		: triangulation_(triangulation), cohomologybasis_(cohomologybasis), dualcohomologybasis_(NULL), minbu_necksize_(MinbuNeckSize), measurements_(0) {}
	BabyUniverseDistribution(const Triangulation * const triangulation, const DualCohomologyBasis * const dualcohomologybasis, int MinbuNeckSize) 
		: triangulation_(triangulation), cohomologybasis_(NULL), dualcohomologybasis_(dualcohomologybasis), minbu_necksize_(MinbuNeckSize), measurements_(0) {}
	~BabyUniverseDistribution() {}

	void FindMinbuNecksOfLength3(std::list<std::list<const Edge *> > & paths);
	void FindMinbuNecksOfLength3(const Edge * edge, std::list<std::list<const Edge *> > & paths);

	void FindMinbuSizes(int necklength, std::vector<int> & sizes);
	void FindMinbuNecks(int necklength, std::list<std::list<const Edge *> > & paths);

	void Measure() 
	{
		if( cohomologybasis_ == NULL )
		{
			CohomologyBasis cohom(triangulation_,*dualcohomologybasis_);
			cohomologybasis_ = &cohom;
			FindMinbuSizes(minbu_necksize_,sizes_);
			cohomologybasis_ = NULL;
		} else
		{
			FindMinbuSizes(minbu_necksize_,sizes_);
		}
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
	const CohomologyBasis * cohomologybasis_;
	const DualCohomologyBasis * dualcohomologybasis_;
	std::vector<int> sizes_;
	int minbu_necksize_;
	int measurements_;
};

#endif

