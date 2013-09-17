#ifndef DECORATION_H
#define DECORATION_H

#include <string>

#include "Triangulation.h"
#include "Edge.h"

class Decoration
{
public:
	Decoration() {}
	Decoration(const Triangulation * const triangulation) : triangulation_(triangulation), last_state_(0) {}
	~Decoration(void) {}

	virtual void Initialize() = 0;

	virtual void UpdateAfterFlipMove(const Edge * const) {}
	virtual std::string ExportState() const { return ""; }

	bool IsUpToDate() const { 
		return triangulation_->IsState(last_state_);
	}
	virtual bool MakeUpToDate() { return true; }
	void SetUpToDate() {
		last_state_ = triangulation_->getState();
	}
private:
	const Triangulation * triangulation_;
	TriangulationState last_state_;
};

#endif

