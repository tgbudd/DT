#ifndef PEELING_PROCEDURE_H
#define PEELING_PROCEDURE_H

#include "Observable.h"
#include "Triangulation.h"
#include "Edge.h"
#include "BabyUniverseDetector.h"
#include "CohomologyBasis.h"

class Snapshot;

class PeelingProcedure :
	public Observable
{
public:
	///// temporary /////
	friend class Snapshot;
	friend class PeelingMeasurement;
	/////////////////////

	PeelingProcedure(const Triangulation * const triangulation);
	PeelingProcedure(const Triangulation * const triangulation, CohomologyBasis * const cohomologybasis);
	~PeelingProcedure();
	void Measure();
	std::string OutputData() const;
private:
	void DoPeeling();
	void DoPeelingOnTorus(Triangle * startTriangle);
	bool DoPeelingStepOnTorus();
	bool DoPeelingStepOnSphere();
	void PreparePeelingOnSphere(Triangle * startTriangle);
	void InitializePeeling();
	void InitializePeeling(Triangle * startTriangle);
	bool PeelingStep();
	bool FirstLoopContainsFinalVertex(std::list<const Edge*>::const_iterator loopBegin) const;
	bool FirstLoopIsDisk(std::list<const Edge*>::const_iterator loopBegin) const;
	std::list<const Edge*>::iterator FindSecondLoopBegin();
	void ProcessBabyUniverse(std::list<const Edge*>::iterator secondLoopBegin, bool FirstLoopContainsBabyUniverse);
	Vertex * ChooseFinalVertex();
	void RandomWalk(bool StayInMotherUniverse);
	void DualRandomWalk(bool StayInMotherUniverse);
	Triangle * RandomStartTriangle();

	bool FrontierIsSimpleClosedCurve();
	bool FinalVertexBeyondFrontier();
	bool FinalVertexInLoop(const std::list<const Edge *> & boundary);
	Vertex * RandomNeighbour(Vertex * v);
	
	void DoMeasurementOnBabyUniverse(const std::list<const Edge*>::const_iterator & begin, const std::list<const Edge*>::const_iterator & end, int volume=-1);
	int DoBabyWalk(Vertex * startVertex, int maxTime);

	const Triangulation * const triangulation_;
	BabyUniverseDetector babyuniversedetector_;
	CohomologyBasis * const cohomologybasis_;

	Triangle * last_triangle_;
	Triangle * start_triangle_;
	Vertex * start_vertex_;
	Vertex * final_vertex_;
	std::vector<bool> vertex_in_mother_universe_;
	std::vector<bool> in_mother_universe_;  // for each triangle
	std::vector<bool> in_frontier_;			// for each vertex
	std::list<const Edge*> frontier_;
	int frontier_size_;
	int volume_within_frontier_;
	std::vector<int> distance_;
	std::vector<int> distance_to_final_;
	std::vector<std::vector<int> > frontier_size_trajectory_;
	std::vector<std::vector<int> > volume_trajectory_;
	std::vector<std::vector<int> > distance_trajectory_;
	int max_trajectories_;

	bool random_walk_;
	int random_walk_measurements_;
	int random_walk_measurements_mother_;
	int walk_time_;
	std::vector<int> waiting_times_;
	std::vector<int> distance_square_;
	std::vector<int> distance_square_mother_;
	std::vector<int> distance_square_mother_effective_;
	
	std::vector<Triangle *> triangle_order_;

	bool baby_walk_;
	int baby_walk_samples_;
	std::vector<std::vector<std::vector<int> > > baby_return_time_;
	std::vector<std::vector<int> > baby_universes_;
};

#endif
