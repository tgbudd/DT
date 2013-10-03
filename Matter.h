#ifndef DT_MATTER_H
#define DT_MATTER_H

#include <vector>
#include <string>

#include "Triangulation.h"
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
	virtual double BoltzmannChangeUnderGeneralMove(const std::vector<boost::array<Triangle *,2> > & toBeDeleted, const std::vector<boost::array<Triangle *,2> > & toBeAdded ) const { 
		return 1.0; 
	}
	virtual void UpdateAfterFlipMove(const Edge * const) {}
	virtual void DoSweep() = 0;

private:
	//const Triangulation * triangulation_;

};

#endif
