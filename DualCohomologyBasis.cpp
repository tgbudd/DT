#include "DualCohomologyBasis.h"


DualCohomologyBasis::~DualCohomologyBasis(void)
{
}

void DualCohomologyBasis::Initialize(int width, int height)
{
	CohomologyBasis cohom(triangulation_);
	cohom.Initialize(width,height);
	SetToDualOf(cohom);
}

DualCohomologyBasis::DualCohomologyBasis(const CohomologyBasis & cohomologybasis)
	: CohomologyBasis(cohomologybasis)
{
	SetToDualOf(cohomologybasis);
}

DualCohomologyBasis::DualCohomologyBasis(Triangulation * const triangulation, const std::vector<std::list<Edge*> > & generators, const std::vector<IntForm2D> & integrals )
	: CohomologyBasis(triangulation)
{
	omega_.resize(triangulation_->NumberOfTriangles());
	ClearOmega();
	
	// construct a basis of dual closed forms such that the integral along a path homotopic to generators[i] gives integrals[i]
	BOOST_ASSERT(generators.size() == 2);
	boost::array<IntForm2D,2> forms = {NegateForm(integrals[1]),integrals[0]};
	for(int i=0;i<2;i++)
	{
		for( std::list<Edge*>::const_iterator edgeIt = generators[i].begin(); edgeIt != generators[i].end(); edgeIt++)
		{
			addToOmega(*edgeIt,forms[i]);
			setOmegaToMinusAdjacent((*edgeIt)->getAdjacent());
		}
	}

	// The pair of generators might not be properly oriented. Therefore we have to check whether
	// the overall sign is correct.


}

void DualCohomologyBasis::SetToDualOf(const CohomologyBasis & cohom)
{
	omega_.resize(triangulation_->NumberOfTriangles());

	// Take the dual form from triangles t to t' to be the form integrated from the 
	// vertex at position 0 of t to the vertex at position 0 of t'.
	
	for(int triangleId1=0;triangleId1<triangulation_->NumberOfTriangles();triangleId1++)
	{
		Triangle * triangle1 = triangulation_->getTriangle(triangleId1);
		for(int direction1=0;direction1<3;direction1++)
		{
			Edge * edge = triangle1->getEdge(direction1);
			Triangle * triangle2 = edge->getAdjacent()->getParent();
			int direction2 = edge->getAdjacent()->getId();
			for(int x=0;x<2;x++)
			{
				setOmega(edge,x, (direction1==0 ? cohom.getOmega(triangle1->getEdge(2),x) : (direction1 == 1 ? - cohom.getOmega(triangle1->getEdge(1),x) : 0) )
							    +(direction2==0 ? cohom.getOmega(triangle2->getEdge(1),x) : (direction2 == 1 ? 0 : - cohom.getOmega(triangle2->getEdge(2),x)) ) );
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

bool DualCohomologyBasis::CheckClosedness() const
{
	for(int i=0;i<triangulation_->NumberOfVertices();i++)
	{
		Vertex * vertex = triangulation_->getVertex(i);
		for(int j=0;j<2;j++)
		{
			int total = 0;
			Edge * startEdge = vertex->getParent()->getNext();
			Edge * edge = startEdge;
			do
			{
				total += getOmega(edge,j);
				edge = edge->getAdjacent()->getPrevious();
			} while( edge != startEdge );
			if( total != 0 )
				return false;
		}
	}
	return true;
}

IntForm2D DualCohomologyBasis::IntegrateToParent(Edge * edge) const
{
	// Consider the path starting at edge->getParent() and ending at the parent of the end-point of edge in anti-clockwise direction.
	// Return the integral of dualOmega along this path.

	IntForm2D integral = {0,0};
	while( edge->getPrevious() != edge->getPrevious()->getOpposite()->getParent() )
	{
		integral = AddForms(integral,getOmega(edge));
		edge = edge->getAdjacent()->getPrevious();
	}
	return integral;
}