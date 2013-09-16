#ifndef CONNECTIVITY_RESTRICTOR_H
#define CONNECTIVITY_RESTRICTOR_H

#include "Triangulation.h"
#include "Matter.h"
#include "CohomologyBasis.h"

class ConnectivityRestrictor :
	public Matter
{
public:
	enum Restriction {
		NO_DOUBLE_EDGES,
		NO_CONTRACTIBLE_DOUBLE_EDGES,
		NO_SELF_LOOPS,
		NO_CONTRACTIBLE_SELF_LOOPS
	};

	ConnectivityRestrictor(Triangulation * const & triangulation, Restriction restriction)
		: triangulation_(triangulation), cohomologybasis_(NULL), restriction_(restriction)  {
	}
	ConnectivityRestrictor(Triangulation * const & triangulation, const CohomologyBasis * const cohomologybasis, Restriction restriction)
		: triangulation_(triangulation), cohomologybasis_(cohomologybasis), restriction_(restriction) {
	}
	~ConnectivityRestrictor() {}
	void DoSweep() {}
	void UpdateAfterFlipmove(const Edge * const edge) {}
	void Initialize() {}
	bool IsFlipMoveAllowed(const Edge * const edge);

private:
	Triangulation * const triangulation_;
	const CohomologyBasis * const cohomologybasis_;
	Restriction restriction_;
};

#endif
