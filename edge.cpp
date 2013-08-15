#include "Edge.h"


Edge::~Edge(void)
{
}


bool Edge::IsFlipMovePossible() const
{
	// Disallow if the triangle is adjacent to itself along this edge.
	// Also, if diagonals are connected a move will not change the triangulation, so don't allow it.
	return (getParent() != getAdjacent()->getParent())	
		&& (getNext()->getAdjacent() != getAdjacent()->getNext())
		&& (getPrevious()->getAdjacent() != getAdjacent()->getPrevious()) ;
}

void Edge::DoFlipMove()
{
	Edge * adjacent = adjacent_;

	// make sure the vertices keep a correct parent
	getNext()->getOpposite()->setParent(adjacent->getPrevious());
	adjacent->getNext()->getOpposite()->setParent(getPrevious());

	// update the adjacency of the triangles getParent() and adjacent->getParent()
	bindAdjacent(adjacent->getPrevious()->getAdjacent());
	adjacent->bindAdjacent(getPrevious()->getAdjacent());
	getPrevious()->bindAdjacent(adjacent->getPrevious());

	// update the vertices
	getNext()->setOpposite( adjacent->getOpposite() );
	adjacent->getNext()->setOpposite( getOpposite() );

	BOOST_ASSERT(getNext()->getOpposite() == getNext()->getOpposite()->getParent()->getOpposite());
}
