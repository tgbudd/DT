#include <queue>
#include <set>
#include <list>

#include "ShortestLoop.h"

void ShortestLoop::FindShortestLoop(int maxLength)
{
	BOOST_ASSERT( !generators_.empty() );
	std::set<Vertex *> setVertices;
	for(std::vector<std::list<Edge*> >::iterator pathIt = generators_.begin(); pathIt != generators_.end(); pathIt++)
	{
		for(std::list<Edge*>::iterator edgeIt = pathIt->begin(); edgeIt != pathIt->end(); edgeIt++)
		{
			setVertices.insert((*edgeIt)->getPrevious()->getOpposite());
		}
		if( maxLength == -1 || (pathIt->size() < maxLength) )
		{
			maxLength = pathIt->size();
		}
	}
	std::list<Vertex *> startVertices(setVertices.begin(),setVertices.end());

	shortestloop_ = FindShortestLoop(startVertices, maxLength);
}

void ShortestLoop::FindGenerators()
{
	generators_.clear();
	generator_integrals_.clear();
	generators_.reserve(2);

	Vertex * startVertex = triangulation_->getVertex(0);

	generators_.push_back(FindShortestLoop(startVertex));
	BOOST_ASSERT(CheckPathIsLoop(generators_.back()));
	generator_integrals_.push_back(cohomologybasis_->Integrate(generators_.back()));
	generators_.push_back(FindShortestLoop(startVertex,generator_integrals_.back()));
	BOOST_ASSERT(CheckPathIsLoop(generators_.back()));
	generator_integrals_.push_back(cohomologybasis_->Integrate(generators_.back()));

	// Check that the orientation is ok. The determinant of generator_integrals_ should be 1.
	int determinant = generator_integrals_[0][0]*generator_integrals_[1][1]-generator_integrals_[0][1]*generator_integrals_[1][0];
	BOOST_ASSERT( abs(determinant) == 1 );
	if( determinant == -1 )
	{
		std::swap(generators_[0],generators_[1]);
		std::swap(generator_integrals_[0],generator_integrals_[1]);
	}
}

std::list<Edge*> ShortestLoop::FindShortestLoop(Vertex * startVertex, int maxLength )
{
	IntForm2D zeroform = {0,0};
	return FindShortestLoop(startVertex,zeroform,maxLength);
}

std::list<Edge*> ShortestLoop::FindShortestLoop(const std::list<Vertex *> & startVertices, int maxLength  )
{
	IntForm2D zeroform = {0,0};
	return FindShortestLoop(startVertices,zeroform,maxLength);
}

std::list<Edge*> ShortestLoop::FindShortestLoop(const std::list<Vertex *> & startVertices, const IntForm2D & notMultipleOf, int maxLength  )
{
	int maxlength = maxLength;
	std::list<Edge*> shortestpath; 
	for(std::list<Vertex *>::const_iterator vert = startVertices.begin(); vert != startVertices.end(); vert++)
	{
		const std::list<Edge*> path = FindShortestLoop(*vert,notMultipleOf,maxlength);
		if( !path.empty() && (maxlength == -1 || path.size() <= maxlength) )
		{
			maxlength = ((int)path.size())-1;
			shortestpath = path;
		}
	}
	return shortestpath;
}
	
std::list<Edge*> ShortestLoop::FindShortestLoop(Vertex * startVertex, const IntForm2D & notMultipleOf, int maxLength )
{
	if( visit.size() != triangulation_->NumberOfVertices() )
	{
		visit.resize(triangulation_->NumberOfVertices());
	}
	for(int i=0;i<visit.size();i++)
	{
		// set all vertices to unvisited
		visit[i].distance = -1;
	}

	std::list<Edge*> ShortestPath;

	std::queue<Vertex *> q;
	q.push(startVertex);
	visit[startVertex->getId()].distance = 0;
	visit[startVertex->getId()].parent = NULL;
	visit[startVertex->getId()].integral = IntForm2D();

	while( !q.empty() )
	{
		Vertex * vertex = q.front();
		int v = vertex->getId();
		q.pop();
	
		// if we have arrived at len2/2 stop
		if( maxLength != -1 && 2*visit[v].distance > maxLength ) 
			break;

		// scan through its neighbours
		Edge * firstEdge = vertex->getParent()->getPrevious();
		Edge * edge = firstEdge;
		do {
			Vertex * nbrVertex = edge->getPrevious()->getOpposite();
			int v2 = nbrVertex->getId();
			IntForm2D dx = cohomologybasis_->getOmega(edge);
			IntForm2D totalForm = AddForms(visit[v].integral,dx);
			if( visit[v2].distance == -1 )
			{
				visit[v2].distance = visit[v].distance + 1;
				visit[v2].parent = edge;
				visit[v2].integral = totalForm;
				q.push(nbrVertex);
			}else if( !FormIsZero(SubtractForms(totalForm, visit[v2].integral)) 
				&& ( FormIsZero(notMultipleOf) || FormsIndependent(SubtractForms(totalForm, visit[v2].integral),notMultipleOf) ) )
			{
				// non-trivial loop
				if( maxLength == -1 || visit[v].distance + visit[v2].distance + 1 <= maxLength )
				{
					maxLength = visit[v].distance + visit[v2].distance;

					// store path
					ShortestPath.clear();
					Vertex * current = vertex;
					while( current != startVertex )
					{
						ShortestPath.push_back(visit[current->getId()].parent);
						current = visit[current->getId()].parent->getNext()->getOpposite();
					}
					std::reverse(ShortestPath.begin(),ShortestPath.end());
					ShortestPath.push_back(edge);
					current = nbrVertex;
					while( current != startVertex )
					{
						ShortestPath.push_back(visit[current->getId()].parent->getAdjacent());
						current = visit[current->getId()].parent->getNext()->getOpposite();
					}
				}
				if( visit[v2].distance <= visit[v].distance || maxLength == 0 )
				{
					return ShortestPath;	
				}
			}
			edge = edge->getPrevious()->getAdjacent();
		} while( edge != firstEdge );
	}
	return ShortestPath;
}

bool ShortestLoop::CheckPathIsLoop(const std::list<Edge*> & path) const
{
	if( path.empty() )
	{
		return false;
	}
	Vertex * previousVertex = path.back()->getPrevious()->getOpposite();
	for(std::list<Edge *>::const_iterator edgeIt = path.begin(); edgeIt!=path.end(); edgeIt++)
	{
		if( (*edgeIt)->getNext()->getOpposite() != previousVertex )
			return false;
		previousVertex = (*edgeIt)->getPrevious()->getOpposite();
	}
	return true;
}