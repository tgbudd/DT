#pragma once

#include <vector>
#include <string>

#include "triangulation.h"
#include "Decoration.h"
#include "Edge.h"

class Matter : public Decoration
{
public:
	Matter() {}
	~Matter(void) {}
	virtual bool IsFlipMoveAllowed(const Edge * const) {
		return true;
	}
	virtual void Initialize() = 0;
	virtual double BoltzmannChangeUnderFlipMove(const Edge * const ) const { 
		return 1.0; 
	}
	virtual void UpdateAfterFlipMove(const Edge * const) {}
	virtual void DoSweep() = 0;

private:
	//const Triangulation * triangulation_;

};

