#ifndef DUAL_RANDOM_WALK_H
#define DUAL_RANDOM_WALK_H

#include "Triangulation.h"
#include "observable.h"

class DualRandomWalk :
	public Observable
{
public:
	DualRandomWalk(const Triangulation * const triangulation);
	~DualRandomWalk();
};

#endif
