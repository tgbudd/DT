#include <vector>
#include <queue>
#include <algorithm>

#include "TriangulationProperties.h"
#include "Vertex.h"
#include "Edge.h"

void properties::DegreeList(const Triangulation * const triangulation, std::vector<int> & degree)
{
	if( degree.size() != triangulation->NumberOfVertices() )
	{
		degree.resize(triangulation->NumberOfVertices());
	}
	for(int i=0,end=triangulation->NumberOfVertices();i<end;i++)
	{
		Vertex * vertex = triangulation->getVertex(i);
		Edge * edge = vertex->getParent()->getPrevious();
		degree[i]=0;
		do {
			edge = edge->getPrevious()->getAdjacent();
			++degree[i];
		} while( edge->getNext() != vertex->getParent() );
	}
}

void properties::VertexDistanceList(const Triangulation * const triangulation, const Vertex * const startVertex, std::vector<int> & distance)
{
	if( distance.size() != triangulation->NumberOfVertices() )
	{
		distance.resize(triangulation->NumberOfVertices());
	}
	std::fill(distance.begin(),distance.end(),-1);

	std::queue<const Vertex *> q;
	q.push( startVertex );
	distance[startVertex->getId()] = 0;

	while( !q.empty() )
	{
		const Vertex * vertex = q.front();
		q.pop();

		Edge * edge = vertex->getParent()->getPrevious();
		do {
			Vertex * nbr = edge->getPrevious()->getOpposite();

			if( distance[nbr->getId()] == -1 )
			{
				distance[nbr->getId()] = distance[vertex->getId()] + 1;
				q.push(nbr);
			}
			edge = edge->getPrevious()->getAdjacent();
		} while( edge->getNext() != vertex->getParent() );
	}
}
