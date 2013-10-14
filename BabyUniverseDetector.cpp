#include <list>
#include <queue> 
#include "boost/array.hpp"

#include "Edge.h"
#include "BabyUniverseDetector.h"


BabyUniverseDetector::BabyUniverseDetector(const Triangulation * triangulation)
	: triangulation_(triangulation), visited_(triangulation->NumberOfTriangles()), genus_(triangulation->CalculateGenus())
{
}


BabyUniverseDetector::~BabyUniverseDetector()
{
}

std::pair<int,bool> BabyUniverseDetector::VolumeEnclosed(const std::list<const Edge *> & boundary) const
{
	BOOST_ASSERT( genus_ == 0 );
	return VolumeEnclosedSpherical(boundary);
}

std::pair<int,bool> BabyUniverseDetector::VolumeEnclosedSpherical(const std::list<const Edge *> & boundary) const
{
	// To find the volume enclosed by the boundary: perform a breadth-first 
	// search on both sides simultaneously.
	// The volume is returned together with a bool saying whether the edges
	// in the boundary are on the smallest side.
		
	boost::array<std::queue<Triangle *>,2> queues;
	visited_.Reset();

	queues[0].push(boundary.front()->getParent());
	queues[1].push(boundary.front()->getAdjacent()->getParent());
	visited_.Set(queues[0].front()->getId());
	visited_.Set(queues[1].front()->getId());

	int numTriangles = 0;
	while( !queues[0].empty() && !queues[1].empty() )
	{
		for(int i=0;i<2;i++)
		{
			Triangle * triangle = queues[i].front();
			queues[i].pop();
			for(int j=0;j<3;j++)
			{
				Edge * edge = triangle->getEdge(j);
				Triangle * nbrTriangle = edge->getAdjacent()->getParent();
				if( !visited_.isSet(nbrTriangle->getId()) &&
					std::find(boundary.begin(),boundary.end(), (i==0? edge : edge->getAdjacent())) == boundary.end() )
				{
					queues[i].push(nbrTriangle);
					visited_.Set(nbrTriangle->getId());
				}
			}
		}
		numTriangles++;
	}
	return std::pair<int,bool>(numTriangles,queues[0].empty());
}

void BabyUniverseDetector::EnclosedTriangles(const std::list<const Edge *> & boundary, std::vector<Triangle *> & triangles, bool thisSide) const
{
	std::queue<Triangle *> queue;
	visited_.Reset();
	if( thisSide )
	{
		queue.push(boundary.front()->getParent());
	} else
	{
		queue.push(boundary.front()->getAdjacent()->getParent());
	}
	visited_.Set(queue.front()->getId());
	while( !queue.empty() )
	{
		Triangle * triangle = queue.front();
		queue.pop();
		triangles.push_back(triangle);
		for(int j=0;j<3;j++)
		{
			Edge * edge = triangle->getEdge(j);
			Triangle * nbrTriangle = edge->getAdjacent()->getParent();
			if( !visited_.isSet(nbrTriangle->getId()) &&
				std::find(boundary.begin(),boundary.end(),(thisSide? edge : edge->getAdjacent())) == boundary.end() )
			{
				queue.push(nbrTriangle);
				visited_.Set(nbrTriangle->getId());
			}
		}
	}
}
