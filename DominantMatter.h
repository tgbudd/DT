#ifndef DOMINANT_MATTER_H
#define DOMINANT_MATTER_H

#include <string>

class DominantMatter
{
public:
	DominantMatter() {}
	~DominantMatter() {}
	virtual void DoSweep() = 0;
	virtual std::string ExportState() const { return ""; }
};

#endif
