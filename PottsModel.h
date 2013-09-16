#ifndef POTTS_MODEL_H
#define POTTS_MODEL_H

#include "Matter.h"
#include "Triangulation.h"


class PottsModel :
	public Matter
{
public:
	PottsModel() : triangulation_(NULL) , states_(2) {}
	PottsModel(Triangulation * const triangulation);
	PottsModel(Triangulation * const triangulation, int states);
	~PottsModel(void) {}

	void Initialize();

	void setStates(const int & states ) {
		states_ = states;
	}
	void setInteraction(const double & interaction ) {
		interaction_ = interaction;
		edgeweight_ = exp( - 2.0 * interaction_ );
	}

	double BoltzmannChangeUnderFlipMove(const Edge * const ) const;
	void UpdateAfterFlipMove(const Edge * const);
	void DoSweep();

	int getSpin(int i) const {
		return spin_[i];
	}

	std::string ExportState() const;


private:
	int DoWolffMove();

	Triangulation * triangulation_;
	int states_;	// 2 states = Ising model
	double interaction_;
	double edgeweight_;

	std::vector<int> spin_;
};

#endif
