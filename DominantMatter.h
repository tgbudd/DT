#ifndef DOMINANT_MATTER_H
#define DOMINANT_MATTER_H

#include "Triangulation.h"


class DominantMatter
{
public:
	DominantMatter() : triangulation_(NULL) {}
	~DominantMatter() {}
	virtual void DoSweep() = 0;
	virtual std::string ExportState() const { return ""; }
private:
	Triangulation * const triangulation_;
};

#endif
