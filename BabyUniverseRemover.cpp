#include <map>
#include <queue>

#include "BabyUniverseRemover.h"
#include "Vertex.h"
#include "Edge.h"

BabyUniverseRemover::BabyUniverseRemover(Triangulation * triangulation)
	: triangulation_(triangulation), 
	cohomologybasis_(NULL),
	genus_(triangulation_->CalculateGenus())
{

}

BabyUniverseRemover::BabyUniverseRemover(Triangulation * triangulation, CohomologyBasis * cohomologybasis)
	: triangulation_(triangulation), 
	cohomologybasis_(cohomologybasis),
	genus_(triangulation_->CalculateGenus())
{

}

void BabyUniverseRemover::RemoveBabyUniverses()
{
	if( genus_ == 0 )
	{
		RemoveBabyUniversesSpherical();
	} else
	{
		RemoveBabyUniversesHigherGenus();
	}
}

void BabyUniverseRemover::RemoveBabyUniversesSpherical()
{
	// Identify all pairs of edges that form a neck of length 2
	std::vector<boost::array<Edge *,3> > newNbr(triangulation_->NumberOfTriangles());
	for(int i=0;i<triangulation_->NumberOfVertices();i++)
	{
		Vertex * vertex = triangulation_->getVertex(i);
		std::map<Vertex*,Edge*> edgeToVertex;
		Edge * edge = vertex->getParent()->getPrevious();
		bool firstround = true;
		while(true)
		{
			if( edge->getPrevious()->getOpposite() == vertex )
			{
				newNbr[edge->getParent()->getId()][edge->getId()] = edge->getAdjacent();
			} else
			{
				std::pair<std::map<Vertex*,Edge*>::iterator,bool> ret = edgeToVertex.insert(std::pair<Vertex*,Edge*>(edge->getPrevious()->getOpposite(),edge));
				if( !ret.second )
				{
					newNbr[ret.first->second->getAdjacent()->getParent()->getId()][ret.first->second->getAdjacent()->getId()] = edge;
					newNbr[edge->getParent()->getId()][edge->getId()] = ret.first->second->getAdjacent();
					ret.first->second = edge;
				}
			}

			edge = edge->getAdjacent()->getNext();
			if( edge == vertex->getParent()->getPrevious() )
			{
				if( firstround )
				{
					firstround = false;
				} else
				{
					break;
				}
			}
		}
	}

	// Determine largest component
	std::vector<int> flag(triangulation_->NumberOfTriangles(),-1);
	int largest = 0;
	int largestFlag = 0;
	std::queue<Triangle *> q;
	for(int i=0;i<triangulation_->NumberOfTriangles();i++)
	{
		if( flag[i] != -1 )
			continue;

		int size = 0;
		q.push(triangulation_->getTriangle(i));
		flag[i] = i;
		while( !q.empty() )
		{
			Triangle * t = q.front();
			q.pop();
			size++;
			for(int j=0;j<3;j++)
			{
				Triangle * nbr = newNbr[t->getId()][j]->getParent();
				BOOST_ASSERT( flag[nbr->getId()] == -1 || flag[nbr->getId()] == i );
				if( flag[nbr->getId()] == -1 )
				{
					flag[nbr->getId()] = i;
					q.push(nbr);
				}
			}
		}
		if( size > largest )
		{
			largest = size;
			largestFlag = i;
		}
	}

	// Update neighbour info
	for(int i=0;i<triangulation_->NumberOfTriangles();i++)
	{
		Triangle * triangle = triangulation_->getTriangle(i);
		if( flag[i] == largestFlag )
		{
			for(int j=0;j<3;j++)
			{
				triangle->getEdge(j)->setAdjacent(newNbr[i][j]);
			}
		}
	}

	// Delete all triangles that do not have flag=largestFlag
	std::vector<Triangle *> tmpTriangles;
	for(int i=0;i<triangulation_->NumberOfTriangles();i++)
	{
		tmpTriangles.push_back(triangulation_->getTriangle(i));
	}
	for(int i=0,endi=tmpTriangles.size();i<endi;i++)
	{
		if( flag[i] != largestFlag )
		{
			triangulation_->DeleteTriangle(tmpTriangles[i]);
		}
	}

	triangulation_->DetermineVertices();
}

void BabyUniverseRemover::RemoveBabyUniversesHigherGenus()
{
	// Identify all pairs of edges that form a neck of length 2
	std::vector<boost::array<Edge *,3> > newNbr(triangulation_->NumberOfTriangles());
	for(int i=0;i<triangulation_->NumberOfVertices();i++)
	{
		Vertex * vertex = triangulation_->getVertex(i);
		std::map<std::pair<IntForm2D,Vertex*>,Edge*> edgeToVertex;
		Edge * edge = vertex->getParent()->getPrevious();
		bool firstround = true;
		while(true)
		{
			if( edge->getPrevious()->getOpposite() == vertex )
			{
				newNbr[edge->getParent()->getId()][edge->getId()] = edge->getAdjacent();
			} else
			{
				IntForm2D form = cohomologybasis_->getOmega(edge);
				std::pair<std::map<std::pair<IntForm2D,Vertex*>,Edge*>::iterator,bool> ret = edgeToVertex.insert(std::pair<std::pair<IntForm2D,Vertex*>,Edge*>(std::pair<IntForm2D,Vertex*>(form,edge->getPrevious()->getOpposite()),edge));
				if( !ret.second )
				{
					newNbr[ret.first->second->getAdjacent()->getParent()->getId()][ret.first->second->getAdjacent()->getId()] = edge;
					newNbr[edge->getParent()->getId()][edge->getId()] = ret.first->second->getAdjacent();
					ret.first->second = edge;
				}
			}

			edge = edge->getAdjacent()->getNext();
			if( edge == vertex->getParent()->getPrevious() )
			{
				if( firstround )
				{
					firstround = false;
				} else
				{
					break;
				}
			}
		}
	}

	// Update neighbour info
	for(int i=0;i<triangulation_->NumberOfTriangles();i++)
	{
		Triangle * triangle = triangulation_->getTriangle(i);
		for(int j=0;j<3;j++)
		{
			triangle->getEdge(j)->setAdjacent(newNbr[i][j]);
		}
	}
	triangulation_->DetermineVertices();

	// Determine which component is the torus
	std::vector<int> flag(triangulation_->NumberOfTriangles(),-1);
	std::vector<int> visited(triangulation_->NumberOfVertices(),false);
	int torusFlag = 0;
	std::queue<Triangle *> q;
	for(int i=0;i<triangulation_->NumberOfTriangles();i++)
	{
		if( flag[i] != -1 )
			continue;

		int nTriangles = 0;
		int nVertices = 0;
		q.push(triangulation_->getTriangle(i));
		flag[i] = i;
		while( !q.empty() )
		{
			Triangle * t = q.front();
			q.pop();
			nTriangles++;
			for(int j=0;j<3;j++)
			{
				Vertex * vertex = t->getEdge(j)->getOpposite();
				if( !visited[vertex->getId()] )
				{
					visited[vertex->getId()] = true;
					nVertices++;
				}

				Triangle * nbr = t->getEdge(j)->getAdjacent()->getParent();
				BOOST_ASSERT( flag[nbr->getId()] == -1 || flag[nbr->getId()] == i );
				if( flag[nbr->getId()] == -1 )
				{
					flag[nbr->getId()] = i;
					q.push(nbr);
				}
			}
		}
		BOOST_ASSERT( nTriangles - 2 * nVertices == 0 || nTriangles - 2 * nVertices == -4 );
		if( nTriangles - 2 * nVertices == 0 )
		{
			torusFlag = i;
			break;
		}
	}

	// Delete all triangles that do not have flag=torusFlag
	std::vector<Triangle *> tmpTriangles;
	for(int i=0;i<triangulation_->NumberOfTriangles();i++)
	{
		tmpTriangles.push_back(triangulation_->getTriangle(i));
	}
	for(int i=0,endi=tmpTriangles.size();i<endi;i++)
	{
		if( flag[i] != torusFlag )
		{
			triangulation_->DeleteTriangle(tmpTriangles[i]);
		}
	}

	triangulation_->DetermineVertices();
}