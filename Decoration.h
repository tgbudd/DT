#pragma once

#include <string>

#include "triangulation.h"
#include "Edge.h"

class Decoration
{
public:
	Decoration() {}
	Decoration(const Triangulation * const triangulation) : triangulation_(triangulation) {}
	~Decoration(void) {}

	virtual void Initialize() = 0;

	virtual void UpdateAfterFlipMove(const Edge * const) {}
	virtual std::string PrintState() const { return ""; }

private:
	const Triangulation * triangulation_;

};
