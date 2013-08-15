#include <boost/array.hpp>

#include "Triangle.h"
#include "Edge.h"

Triangle::Triangle()
{
	for(int i=0;i<3;i++)
	{
		edges_[i] = new Edge(this,i);
	}
	for(int i=0;i<3;i++)
	{
		edges_[i]->setNext(edges_[(i+1)%3]);
		edges_[i]->setPrevious(edges_[(i+2)%3]);
	}
}

Triangle::~Triangle()
{
	for(int i=0;i<3;i++)
	{
		delete edges_[i];
	}
}
