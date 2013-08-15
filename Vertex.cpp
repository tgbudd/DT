#include "Vertex.h"



Vertex::~Vertex(void)
{
}

int Vertex::getDegree() const
{
	int degree = 0;
	Edge * edge = getParent();
	do
	{
		edge = edge->getNext()->getAdjacent()->getNext();
		degree++;
	} while( edge != getParent() );
	return degree;
}