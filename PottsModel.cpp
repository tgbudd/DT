#include <stack>
#include <algorithm>
#include <sstream>

#include "PottsModel.h"

PottsModel::PottsModel(Triangulation * const triangulation) : triangulation_(triangulation), states_(2) {
	
}

PottsModel::PottsModel(Triangulation * const triangulation, int states) : triangulation_(triangulation), states_(states) {
	
}

void PottsModel::Initialize()
{
	spin_.resize(triangulation_->NumberOfTriangles());
	std::fill(spin_.begin(),spin_.end(),0);
}

double PottsModel::BoltzmannChangeUnderFlipMove(const Edge * const edge) const {
	int SpinOfFirst =			 spin_[edge->getParent()->getId()];
	int SpinOfSecond =			 spin_[edge->getAdjacent()->getParent()->getId()];
	int SpinOfNeigbourOfFirst =	 spin_[edge->getPrevious()->getAdjacent()->getParent()->getId()];
	int SpinOfNeigbourOfSecond = spin_[edge->getAdjacent()->getPrevious()->getAdjacent()->getParent()->getId()];
	
	int ChangeInNumberOfFrustratedEdges = 0;
	if( SpinOfNeigbourOfFirst == SpinOfSecond ) ChangeInNumberOfFrustratedEdges++;
	if( SpinOfNeigbourOfSecond == SpinOfFirst ) ChangeInNumberOfFrustratedEdges++;
	if( SpinOfNeigbourOfFirst == SpinOfFirst ) ChangeInNumberOfFrustratedEdges--;
	if( SpinOfNeigbourOfSecond == SpinOfSecond ) ChangeInNumberOfFrustratedEdges--;

	return pow(edgeweight_,ChangeInNumberOfFrustratedEdges);
}

void PottsModel::UpdateAfterFlipMove(const Edge * const) {
	// nothing needs to be done
}

void PottsModel::DoSweep() {
	// perform a sweep of Wolff moves

	int NumberOfFlips = 0;

	while( NumberOfFlips < triangulation_->NumberOfTriangles() )
	{
		NumberOfFlips += DoWolffMove();
	}
}

int PottsModel::DoWolffMove() {
	Triangle * triangle = triangulation_->getRandomTriangle();
	
	int oldSpin = spin_[triangle->getId()];
	int changeSpin = triangulation_->RandomInteger(1,states_-1);
	int newSpin = (oldSpin + changeSpin)%states_;

	spin_[triangle->getId()] = newSpin;

	std::stack<Triangle *> cluster;
	cluster.push(triangle);

	double bondProbability = 1.0 - edgeweight_;

	int NumberOfFlips = 0;

	while( !cluster.empty() )
	{
		Triangle * currentTriangle = cluster.top();
		cluster.pop();
		for(int i=0;i<3;i++)
		{
			Triangle * neighbour = currentTriangle->getEdge(i)->getAdjacent()->getParent();
			if( spin_[neighbour->getId()] == oldSpin && triangulation_->SucceedWithProbability(bondProbability) )
			{
				cluster.push(neighbour);
				spin_[neighbour->getId()] = newSpin;
				NumberOfFlips++;
			}
		}
	}
	return NumberOfFlips;
}

std::string PottsModel::PrintState() const {
	std::ostringstream string;
	string << "{";
	for(int i=0;i<(int)spin_.size();i++)
	{
		string << (i>0?",":"") << spin_[i];
	}
	string << "}";
	return string.str();
}