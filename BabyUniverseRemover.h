#ifndef BABY_UNIVERSE_REMOVER_H
#define BABY_UNIVERSE_REMOVER_H

#include "Triangulation.h"
#include "CohomologyBasis.h"

class BabyUniverseRemover
{
public:
	BabyUniverseRemover(Triangulation * triangulation);
	BabyUniverseRemover(Triangulation * triangulation, CohomologyBasis * cohomologybasis);

	void RemoveBabyUniverses();
private:
	void RemoveBabyUniversesSpherical();
	void RemoveBabyUniversesHigherGenus();
	Triangulation * triangulation_;
	CohomologyBasis * cohomologybasis_;

	int genus_;
};

#endif