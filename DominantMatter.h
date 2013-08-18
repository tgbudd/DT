#pragma once

#include "triangulation.h"


class DominantMatter
{
public:
	DominantMatter() : triangulation_(NULL) {}
	~DominantMatter() {}
	virtual void DoSweep() = 0;
	std::string PrintState() const { return ""; }
private:
	Triangulation * const triangulation_;
};

