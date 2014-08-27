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

int properties::VertexDistance(const Triangulation * const triangulation, const Vertex * const startVertex, const Vertex * const endVertex)
{
	if( startVertex == endVertex )
	{
		return 0;
	}
	std::vector<int> distance(triangulation->NumberOfVertices(),-1);
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
				if( nbr == endVertex )
				{
					return distance[nbr->getId()];
				}
				q.push(nbr);
			}
			edge = edge->getPrevious()->getAdjacent();
		} while( edge->getNext() != vertex->getParent() );
	}
	return -1;
}

void properties::TriangleDistanceList(const Triangulation * const triangulation, const Triangle * const startTriangle, std::vector<int> & distance)
{
	if( distance.size() != triangulation->NumberOfTriangles() )
	{
		distance.resize(triangulation->NumberOfTriangles());
	}
	std::fill(distance.begin(),distance.end(),-1);

	std::queue<const Triangle *> q;
	q.push( startTriangle );
	distance[startTriangle->getId()] = 0;

	while( !q.empty() )
	{
		const Triangle * triangle = q.front();
		q.pop();

		for(int i=0;i<3;i++)
		{
			Triangle * nbr = triangle->getEdge(i)->getAdjacent()->getParent();
			if( distance[nbr->getId()] == -1 )
			{
				distance[nbr->getId()] = distance[triangle->getId()] + 1;
				q.push(nbr);
			}
		}
	}
}

int properties::TriangleDistance(const Triangulation * const triangulation, const Triangle * const startTriangle, const Triangle * const endTriangle)
{
	if( startTriangle == endTriangle )
	{
		return 0;
	}
	std::vector<int> distance(triangulation->NumberOfTriangles(),-1);
	std::queue<const Triangle *> q;
	q.push( startTriangle );
	distance[startTriangle->getId()] = 0;

	while( !q.empty() )
	{
		const Triangle * triangle = q.front();
		q.pop();

		for(int i=0;i<3;i++)
		{
			Triangle * nbr = triangle->getEdge(i)->getAdjacent()->getParent();
			if( distance[nbr->getId()] == -1 )
			{
				distance[nbr->getId()] = distance[triangle->getId()] + 1;
				if( nbr == endTriangle )
				{
					return distance[nbr->getId()];
				}
				q.push(nbr);
			}
		}
	}
	return -1;
}
