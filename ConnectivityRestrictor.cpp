#include "ConnectivityRestrictor.h"

bool ConnectivityRestrictor::IsFlipMoveAllowed(const Edge * const edge) {
	Vertex * otherVertex = edge->getAdjacent()->getOpposite();
	if( edge->getOpposite() == otherVertex )
	{
		if( restriction_ == NO_SELF_LOOPS || restriction_ == NO_DOUBLE_EDGES )
		{
			return false;
		}else if( cohomologybasis_->getOmega(edge->getNext()) == cohomologybasis_->getOmega(edge->getAdjacent()->getPrevious()->getAdjacent()) )
		{
			return false;
		}
	}
	if( restriction_ == NO_DOUBLE_EDGES || restriction_ == NO_CONTRACTIBLE_DOUBLE_EDGES )
	{
		Vertex * otherVertex = edge->getAdjacent()->getOpposite();
		Edge * startEdge = edge->getPrevious();
		Edge * curEdge = startEdge;
		do {
			if( curEdge->getPrevious()->getOpposite() == otherVertex )
			{
				if( restriction_ == NO_DOUBLE_EDGES )
				{
					return false;
				}
				IntForm2D integral = AddForms(cohomologybasis_->getOmega(edge->getNext()),cohomologybasis_->getOmega(edge->getAdjacent()->getPrevious()));
				if( FormIsZero(AddForms(integral,cohomologybasis_->getOmega(curEdge))) )
				{
					return false;
				}
			}
			curEdge = curEdge->getPrevious()->getAdjacent();
		} while(curEdge != startEdge);
	}
	return true;
}