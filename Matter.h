#pragma once

#include <vector>
#include <string>

#include "triangulation.h"
#include "Edge.h"

class Matter
{
public:
	Matter() {}
	Matter(const Triangulation * const triangulation) : triangulation_(triangulation) {}
	~Matter(void) {}
	bool IsFlipMoveAllowed(const Edge * const) {
		return true;
	}
	virtual void Initialize() = 0;
	virtual double BoltzmannChangeUnderFlipMove(const Edge * const ) const = 0;
	virtual void UpdateAfterFlipMove(const Edge * const) = 0;
	virtual void DoSweep() = 0;
	std::string PrintState() const { return ""; }

private:
	const Triangulation * triangulation_;

};

