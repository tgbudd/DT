#ifndef DECORATION_H
#define DECORATION_H

#include <string>

#include "Triangulation.h"
#include "Edge.h"

class Decoration
{
public:
	Decoration() {}
	Decoration(const Triangulation * const triangulation) : triangulation_(triangulation) {}
	~Decoration(void) {}

	virtual void Initialize() = 0;

	virtual void UpdateAfterFlipMove(const Edge * const) {}
	virtual std::string ExportState() const { return ""; }

private:
	const Triangulation * triangulation_;

};

#endif

