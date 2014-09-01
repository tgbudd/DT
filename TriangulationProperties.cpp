#include <vector>
#include <queue>
#include <algorithm>

#include "TriangulationProperties.h"
#include "Vertex.h"
#include "Edge.h"

namespace properties {

void DegreeList(const Triangulation * const triangulation, std::vector<int> & degree)
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

void VertexDistanceList(const Triangulation * const triangulation, const Vertex * const startVertex, std::vector<int> & distance)
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

void VertexDistanceList(const Triangulation * const triangulation, const Vertex * const startVertex, VertexAttribute<int> & distance)
{
	return VertexDistanceList(triangulation,startVertex,distance.data());
}


int VertexDistance(const Triangulation * const triangulation, const Vertex * const startVertex, const Vertex * const endVertex)
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

void TriangleDistanceList(const Triangulation * const triangulation, const Triangle * const startTriangle, std::vector<int> & distance)
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

void TriangleDistanceList(const Triangulation * const triangulation, const Triangle * const startTriangle, TriangleAttribute<int> & distance)
{
	TriangleDistanceList(triangulation,startTriangle,distance.data());
}


int TriangleDistance(const Triangulation * const triangulation, const Triangle * const startTriangle, const Triangle * const endTriangle)
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

template<typename T>
class SecondGreater {
public:
	bool operator()(const std::pair<T,double> & t1, const std::pair<T,double> & t2 )
	{
		return t1.second > t2.second;
	}
};


void VertexWeightedDistanceAlgorithm(const Triangulation * const triangulation, const Vertex * const startVertex, bool withEndVertex, const Vertex * const endVertex, const EdgeAttribute<double> & weight, VertexAttribute<std::pair<double,int> > & distance)
{
	typedef std::pair<const Vertex *,double> VertexDist;

	std::priority_queue<VertexDist, std::vector<VertexDist>, SecondGreater<const Vertex *> > queue;
	distance[startVertex].first = 0.0;
	distance[startVertex].second = 0;
	queue.push(VertexDist(startVertex,0.0));
	while(!queue.empty())
	{
		VertexDist vertex = queue.top();
		queue.pop();

		if( withEndVertex && vertex.first == endVertex )
		{
			return;
		}

		const Edge * edge = vertex.first->getParent()->getPrevious();
		do{
			VertexDist neighbour;
			neighbour.first = edge->getPrevious()->getOpposite();
			neighbour.second = vertex.second + weight[edge];
			if( distance[neighbour.first].second == -1 || neighbour.second < distance[neighbour.first].first )
			{
				distance[neighbour.first].first = neighbour.second;
				distance[neighbour.first].second = distance[vertex.first].second + 1;
				queue.push(neighbour);
			}
			edge = edge->getAdjacent()->getNext();
		} while( edge != vertex.first->getParent()->getPrevious() );
	}
}

void TriangleWeightedDistanceAlgorithm(const Triangulation * const triangulation, const Triangle * const startTriangle, bool withEndTriangle, const Triangle * const endTriangle, const EdgeAttribute<double> & weight, TriangleAttribute<std::pair<double,int> > & distance)
{
	typedef std::pair<const Triangle *,double> TriangleDist;

	std::priority_queue<TriangleDist, std::vector<TriangleDist>, SecondGreater<const Triangle*> > queue;
	distance[startTriangle].first = 0.0;
	distance[startTriangle].second = 0;
	queue.push(TriangleDist(startTriangle,0.0));
	while(!queue.empty())
	{
		TriangleDist triangle = queue.top();
		queue.pop();
		if( withEndTriangle && triangle.first == endTriangle )
		{
			return;
		}

		for(int i=0;i<3;i++)
		{
			TriangleDist neighbour;
			neighbour.first = triangle.first->getEdge(i)->getAdjacent()->getParent();
			neighbour.second = triangle.second + weight[triangle.first->getEdge(i)];
			if( distance[neighbour.first].second == -1 || neighbour.second < distance[neighbour.first].first )
			{
				distance[neighbour.first].first = neighbour.second;
				distance[neighbour.first].second = distance[triangle.first].second + 1;
				queue.push(neighbour);
			}
		}
	}
}

std::pair<double,int> VertexWeightedDistance(const Triangulation * const triangulation, const Vertex * const startVertex, const Vertex * const endVertex, const EdgeAttribute<double> & weight)
{
	VertexAttribute<std::pair<double,int> > distance(triangulation,std::pair<double,int>(0.0,-1));

	VertexWeightedDistanceAlgorithm(triangulation,startVertex,true,endVertex,weight,distance);
	return distance[endVertex];
}

void VertexWeightedDistanceList(const Triangulation * const triangulation, const Vertex * const startVertex, const EdgeAttribute<double> & weight, VertexAttribute<std::pair<double,int> > & distance)
{
	std::fill(distance.begin(),distance.end(),std::pair<double,int>(0.0,-1));
	VertexWeightedDistanceAlgorithm(triangulation,startVertex,false,startVertex,weight,distance);
}

std::pair<double,int> TriangleWeightedDistance(const Triangulation * const triangulation, const Triangle * const startTriangle, const Triangle * const endTriangle, const EdgeAttribute<double> & weight)
{
	TriangleAttribute<std::pair<double,int> > distance(triangulation,std::pair<double,int>(0.0,-1));

	TriangleWeightedDistanceAlgorithm(triangulation,startTriangle,true,endTriangle,weight,distance);
	return distance[endTriangle];
}

void TriangleWeightedDistanceList(const Triangulation * const triangulation, const Triangle * const startTriangle, const EdgeAttribute<double> & weight, TriangleAttribute<std::pair<double,int> > & distance)
{
	std::fill(distance.begin(),distance.end(),std::pair<double,int>(0.0,-1));
	TriangleWeightedDistanceAlgorithm(triangulation,startTriangle,false,startTriangle,weight,distance);
}

}