#include "DualCohomologyBasis.h"


DualCohomologyBasis::~DualCohomologyBasis(void)
{
}

void DualCohomologyBasis::Initialize(int width, int height)
{
	CohomologyBasis::Initialize(width,height);

	Dualize();
}

void DualCohomologyBasis::Dualize()
{
	// Take the dual form from triangles t to t' to be the form integrated from the 
	// vertex at position 0 of t to the vertex at position 0 of t'.

	std::vector<boost::array<IntForm2D,3> > oldOmega(omega_);
	for(int triangleId1=0;triangleId1<triangulation_->NumberOfTriangles();triangleId1++)
	{
		Triangle * triangle = triangulation_->getTriangle(triangleId1);
		for(int direction1=0;direction1<3;direction1++)
		{
			Edge * edge = triangle->getEdge(direction1);
			int triangleId2 = edge->getAdjacent()->getParent()->getId();
			int direction2 = edge->getAdjacent()->getId();
			for(int x=0;x<2;x++)
			{
				setOmega(edge,x, (direction1==0 ? - oldOmega[triangleId1][2][x] : (direction1 == 1 ? oldOmega[triangleId1][1][x] : 0) )
							    +(direction2==0 ? - oldOmega[triangleId2][1][x] : (direction2 == 1 ? 0 : oldOmega[triangleId2][2][x]) ) );
			}
		}
	}
}

void DualCohomologyBasis::UpdateAfterFlipMove(const Edge * const edge)
{
	Edge * otherEdge = edge->getPrevious()->getAdjacent()->getNext();
	
	addToOmega(otherEdge->getNext()->getAdjacent(),getOmega(otherEdge));
	setOmegaToMinusAdjacent(otherEdge->getNext());
	
	addToOmega(edge->getAdjacent(),getOmega(otherEdge));
	setOmegaToMinusAdjacent(edge);

	setOmegaToMinusAdjacent(otherEdge);
	IntForm2D zero = {0,0};
	setOmega(edge->getPrevious(),zero);
	setOmega(otherEdge->getPrevious(),zero);
}
