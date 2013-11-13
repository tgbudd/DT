#include <list>
#include <queue> 
#include <set>
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

std::pair<int,bool> BabyUniverseDetector::VolumeEnclosed(std::list<const Edge *>::const_iterator begin, std::list<const Edge *>::const_iterator end) const
{
	BOOST_ASSERT( genus_ == 0 || genus_ == 1 );
	if( genus_ == 0 )
	{
		return VolumeEnclosedSpherical(begin,end);
	} else
	{
		return VolumeEnclosedTorus(begin,end);
	}
}


std::pair<int,bool> BabyUniverseDetector::VolumeEnclosed(const std::list<const Edge *> & boundary) const
{
	return VolumeEnclosed(boundary.begin(),boundary.end());
}

std::pair<int,bool> BabyUniverseDetector::VolumeEnclosedSpherical(std::list<const Edge *>::const_iterator begin, std::list<const Edge *>::const_iterator end) const
{
	// To find the volume enclosed by the boundary: perform a breadth-first 
	// search on both sides simultaneously.
	// The volume is returned together with a bool saying whether the edges
	// in the boundary are on the smallest side.
		
	boost::array<std::queue<Triangle *>,2> queues;
	visited_.Reset();

	queues[0].push((*begin)->getParent());
	queues[1].push((*begin)->getAdjacent()->getParent());
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
					std::find(begin,end, (i==0? edge : edge->getAdjacent())) == end )
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

std::pair<int,bool> BabyUniverseDetector::VolumeEnclosedTorus(std::list<const Edge *>::const_iterator begin, std::list<const Edge *>::const_iterator end) const
{
	// To find the volum enclosed by the boundary: perform a breadth-first 
	// search on both sides simultaneously. When one of them is finished
	// we can determine the volume on either side. By also counting the number of
	// vertices one can determine via the Euler formula which is the disk.

	boost::array<std::queue<Triangle *>,2> queues;
	visited_.Reset();

	queues[0].push((*begin)->getParent());
	queues[1].push((*begin)->getAdjacent()->getParent());
	visited_.Set(queues[0].front()->getId());
	visited_.Set(queues[1].front()->getId());

	int numTriangles = 0;
	boost::array<std::set<Vertex*>,2> setOfVertices;
	
	while( !queues[0].empty() && !queues[1].empty() )
	{
		for(int i=0;i<2;i++)
		{
			Triangle * triangle = queues[i].front();
			queues[i].pop();
			for(int j=0;j<3;j++)
			{
				Edge * edge = triangle->getEdge(j);
				setOfVertices[i].insert(edge->getOpposite());
				Triangle * nbrTriangle = edge->getAdjacent()->getParent();
				if( !visited_.isSet(nbrTriangle->getId()) &&
					std::find(begin,end, (i==0? edge : edge->getAdjacent())) == end )
				{
					queues[i].push(nbrTriangle);
					visited_.Set(nbrTriangle->getId());
				}
			}
		}
		numTriangles++;
	}
	int smallestSide = ( queues[0].empty() ? 0 : 1 );

	// Euler's formula for the disk
	if( numTriangles + static_cast<int>(std::distance(begin,end)) + 2 == 2 * setOfVertices[smallestSide].size() )
	{
		return std::pair<int,bool>(numTriangles,smallestSide == 0);
	}
	// Take the complement
	return std::pair<int,bool>(triangulation_->NumberOfTriangles() - numTriangles,smallestSide == 1);
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

int BabyUniverseDetector::VolumeEnclosed(std::list<const Edge *>::const_iterator begin, std::list<const Edge *>::const_iterator end, bool thisSide) const
{
	int triangles=0;
	std::queue<Triangle *> queue;
	visited_.Reset();
	if( thisSide )
	{
		queue.push((*begin)->getParent());
	} else
	{
		queue.push((*begin)->getAdjacent()->getParent());
	}
	visited_.Set(queue.front()->getId());
	while( !queue.empty() )
	{
		Triangle * triangle = queue.front();
		queue.pop();
		triangles++;
		for(int j=0;j<3;j++)
		{
			Edge * edge = triangle->getEdge(j);
			Triangle * nbrTriangle = edge->getAdjacent()->getParent();
			if( !visited_.isSet(nbrTriangle->getId()) &&
				std::find(begin,end,(thisSide? edge : edge->getAdjacent())) == end )
			{
				queue.push(nbrTriangle);
				visited_.Set(nbrTriangle->getId());
			}
		}
	}
	return triangles;
}

int BabyUniverseDetector::VolumeEnclosed(const std::list<const Edge *> & boundary, bool thisSide) const
{
	return VolumeEnclosed(boundary.begin(),boundary.end(),thisSide);
}
