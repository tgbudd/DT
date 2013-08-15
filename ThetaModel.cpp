
#include "ThetaModel.h"


ThetaModel::ThetaModel(Triangulation * const triangulation, const CohomologyBasis * const cohomologybasis) : triangulation_(triangulation), cohomologybasis_(cohomologybasis)
{
}


ThetaModel::~ThetaModel(void)
{
}

void ThetaModel::Initialize()
{
	// For now we have to assume that the lattice is regular since we have 
	// no algorithm to determine an admissable set of theta's for an arbitrary triangulation.
	for(int i=0;i<triangulation_->NumberOfVertices();i++)
	{
		BOOST_ASSERT(triangulation_->getVertex(i)->getDegree() == 6);
	}

	boost::array<double,3> equilateral = {PI/3.0,PI/3.0,PI/3.0};
	theta_.assign(triangulation_->NumberOfTriangles(),equilateral);
	
}
