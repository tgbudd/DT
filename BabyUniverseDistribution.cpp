#include <map>
#include <queue>
#include <set>
#include <sstream>

#include "boost/array.hpp"

#include "BabyUniverseDistribution.h"


void BabyUniverseDistribution::FindMinbuNecksOfLength3(std::list<std::list<const Edge *> > & paths)
{
	for(int i=0;i<triangulation_->NumberOfTriangles();i++)
	{
		Triangle * triangle = triangulation_->getTriangle(i);
		for(int j=0;j<3;j++)
		{
			Edge * edge = triangle->getEdge(j);
			// make sure that we do every only once
			if( edge->getAdjacent()->getParent()->getId() > i || (edge->getAdjacent()->getParent() == triangle && edge->getAdjacent()->getId() > j ) )
			{
				FindMinbuNecksOfLength3(triangle->getEdge(j),paths);
			}
		}
	}
}

void BabyUniverseDistribution::FindMinbuSizes(int necklength, std::vector<int> & sizes)
{
	std::list<std::list<const Edge *> > paths;
	FindMinbuNecks(necklength,paths);
	BOOST_FOREACH(std::list<const Edge *> & path, paths)
	{
		int volume = VolumeEnclosed(path);
		BOOST_ASSERT( volume % 2 == 1 && volume >= 3 );
		ResizeAndAdd(sizes,(volume-3)/2,1);
	}
}

void BabyUniverseDistribution::FindMinbuNecks(int necklength, std::list<std::list<const Edge *> > & paths)
{
	if( necklength == 3 )
	{
		FindMinbuNecksOfLength3(paths);
	}
}

void BabyUniverseDistribution::FindMinbuNecksOfLength3(const Edge * edge, std::list<std::list<const Edge *> > & paths)
{
	std::map<VertexNode,Edge *> targets;

	Vertex * startVertex = edge->getNext()->getOpposite();
	Vertex * endVertex = edge->getPrevious()->getOpposite();

	IntForm2D edgeForm = cohomologybasis_->getOmega(edge);
	Edge * stopEdge = edge->getAdjacent()->getPrevious()->getAdjacent();
	for( Edge * currentEdge = edge->getNext()->getAdjacent()->getNext(); currentEdge != stopEdge; currentEdge = currentEdge->getAdjacent()->getNext() )
	{
		VertexNode node;
		node.vertex = currentEdge->getPrevious()->getOpposite();
		node.integral = AddForms(edgeForm,cohomologybasis_->getOmega(currentEdge));
		targets.insert(std::pair<VertexNode,Edge *>(node,currentEdge));
	}

	stopEdge = edge->getPrevious()->getAdjacent();
	for( Edge * currentEdge = edge->getAdjacent()->getNext()->getAdjacent()->getNext(); currentEdge != stopEdge; currentEdge = currentEdge->getAdjacent()->getNext() )
	{
		VertexNode node;
		node.vertex = currentEdge->getPrevious()->getOpposite();
		node.integral = cohomologybasis_->getOmega(currentEdge);
		std::map<VertexNode,Edge*>::iterator it=targets.find(node);
		if( it != targets.end() )
		{
			if( node.vertex->getId() < startVertex->getId() && node.vertex->getId() < endVertex->getId() )
			{
				paths.push_back(std::list<const Edge *>());
				paths.back().push_back(edge);
				paths.back().push_back(it->second);
				paths.back().push_back(currentEdge->getAdjacent());
			}
		}
	}
}

int BabyUniverseDistribution::VolumeEnclosed(const std::list<const Edge *> & boundary) const
{
	// To find the volum enclosed by the boundary: perform a breadth-first 
	// search on both sides simultaneously. When one of them is finished
	// we can determine the volume on either side. By also counting the number of
	// vertices one can determine via the Euler formula which is the disk.

	boost::array<std::queue<Triangle *>,2> queues;
	std::vector<bool> visit(triangulation_->NumberOfTriangles(),false);

	queues[0].push(boundary.front()->getParent());
	queues[1].push(boundary.front()->getAdjacent()->getParent());
	visit[queues[0].front()->getId()] = true;
	visit[queues[1].front()->getId()] = true;

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
				if( !visit[nbrTriangle->getId()] &&
					std::find(boundary.begin(),boundary.end(), (i==0? edge : edge->getAdjacent())) == boundary.end() )
				{
					queues[i].push(nbrTriangle);
					visit[nbrTriangle->getId()] = true;
				}
			}
		}
		numTriangles++;
	}
	int smallestSide = ( queues[0].empty() ? 0 : 1 );

	// Euler's formula for the disk
	if( numTriangles + static_cast<int>(boundary.size()) + 2 == 2 * setOfVertices[smallestSide].size() )
	{
		return numTriangles;
	}
	// Take the complement
	return triangulation_->NumberOfTriangles() - numTriangles;
}

std::string BabyUniverseDistribution::OutputData() const
{
	std::ostringstream output;
	output << "babyuniverses -> { measurements -> " << measurements_;
	output << ", necksize -> " << minbu_necksize_;
	output << ", distribution -> ";
	PrintToStream(output,sizes_.begin(),sizes_.end());
	output << "}";

	return output.str();
}