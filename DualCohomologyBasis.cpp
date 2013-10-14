#include <queue>

#include "DualCohomologyBasis.h"
#include "ShortestLoop.h"

DualCohomologyBasis::~DualCohomologyBasis(void)
{
}


DualCohomologyBasis::DualCohomologyBasis(const CohomologyBasis & cohomologybasis)
	:  CohomologyBasis(cohomologybasis)
{
	SetToDualOf(cohomologybasis);
}

DualCohomologyBasis::DualCohomologyBasis(const Triangulation * const triangulation, const std::vector<std::list<Edge*> > & generators, const std::vector<IntForm2D> & integrals )
	: CohomologyBasis(triangulation)
{
	omega_.resize(triangulation_->NumberOfTriangles());
	
	SetAccordingToGenerators(generators,integrals);
}

void DualCohomologyBasis::Initialize(int width, int height)
{
	CohomologyBasis cohom(triangulation_);
	cohom.Initialize(width,height);
	SetToDualOf(cohom);
}

void DualCohomologyBasis::Initialize()
{
	// Construct a unicellular map dual to a dual spanning tree.
	boost::array<bool,3> alltrue = {true,true,true};
	std::vector<boost::array<bool,3> > inMap(triangulation_->NumberOfTriangles(),alltrue);
	Triangle * startTriangle = triangulation_->getRandomTriangle();

	std::vector<bool> visited(triangulation_->NumberOfTriangles(),false);
	std::vector<Triangle *> queue;

	visited[startTriangle->getId()] = true;
	queue.push_back( startTriangle );

	while( !queue.empty() )
	{
		int entry = triangulation_->RandomInteger(0,queue.size()-1);
		std::swap( queue[entry], queue.back() );
		Triangle * triangle = queue.back();
		queue.pop_back();

		for(int i=0;i<3;i++)
		{
			Edge * edge = triangle->getEdge(i);
			int nbrid = edge->getAdjacent()->getParent()->getId();
			if( !visited[nbrid] )
			{
				inMap[edge->getParent()->getId()][edge->getId()] = false;
				inMap[nbrid][edge->getAdjacent()->getId()] = false;
				visited[nbrid] = true;
				queue.push_back(edge->getAdjacent()->getParent());
			}
		}
	}

	// Traverse the unicellular map and point vertices towards the root.
	std::vector<Edge *> pointToRoot(triangulation_->NumberOfVertices(),NULL);
	visited.resize(triangulation_->NumberOfVertices());
	std::fill(visited.begin(),visited.end(),false);

	Vertex * startVertex = triangulation_->getRandomVertex();
	std::queue<Vertex *> vertexqueue;
	vertexqueue.push(startVertex);
	visited[startVertex->getId()] = true;

	std::vector<Edge *> cycleEdges;
	while( !vertexqueue.empty() && cycleEdges.size() != 2 )
	{
		Vertex * vertex = vertexqueue.front();
		vertexqueue.pop();

		Edge * edge = vertex->getParent()->getPrevious();
		do{
			Vertex * nbr = edge->getPrevious()->getOpposite();
			if( edge != pointToRoot[vertex->getId()] && inMap[edge->getParent()->getId()][edge->getId()] )
			{
				if( visited[nbr->getId()] )
				{
					if( cycleEdges.empty() || cycleEdges[0] != edge->getAdjacent() )
					{
						cycleEdges.push_back(edge);
						if( cycleEdges.size() == 2 )
						{
							break;
						}
					}
				} else
				{
					visited[nbr->getId()] = true;
					pointToRoot[nbr->getId()] = edge->getAdjacent();
					vertexqueue.push( nbr );
				}
			}
			edge = edge->getAdjacent()->getNext();
		} while( edge != vertex->getParent()->getPrevious());
	}

	std::vector<std::list<Edge*> > generators(2);
	for(int i=0;i<2;i++)
	{
		Edge * currentEdge = cycleEdges[i]->getAdjacent();
		while( currentEdge != NULL )
		{
			generators[i].push_back(currentEdge->getAdjacent());
			currentEdge = pointToRoot[currentEdge->getPrevious()->getOpposite()->getId()];
		}
		std::reverse(generators[i].begin(),generators[i].end());
		currentEdge = pointToRoot[cycleEdges[i]->getPrevious()->getOpposite()->getId()];
		while( currentEdge != NULL )
		{
			generators[i].push_back(currentEdge);
			currentEdge = pointToRoot[currentEdge->getPrevious()->getOpposite()->getId()];
		}
	}
	std::vector<IntForm2D> integrals(2);
	integrals[0][0] = 1;
	integrals[0][1] = 0;
	integrals[1][0] = 0;
	integrals[1][1] = 1;

	omega_.resize(triangulation_->NumberOfTriangles());
	ClearOmega();
	SetAccordingToGenerators(generators,integrals);

	SetUpToDate();
}

void DualCohomologyBasis::SetAccordingToGenerators(const std::vector<std::list<Edge*> > & generators, const std::vector<IntForm2D> & integrals )
{
	// construct a basis of dual closed forms such that the integral along a path homotopic to generators[i] gives integrals[i]
	ClearOmega();
	BOOST_ASSERT(generators.size() == 2);
	boost::array<IntForm2D,2> forms = {NegateForm(integrals[1]),integrals[0]};
	for(int i=0;i<2;i++)
	{
		for( std::list<Edge*>::const_iterator edgeIt = generators[i].begin(); edgeIt != generators[i].end(); edgeIt++)
		{
			addToOmega(*edgeIt,forms[i]);
			setOmegaToMinusAdjacent((*edgeIt)->getAdjacent());
		}
	}

}

void DualCohomologyBasis::SetToDualOf(const CohomologyBasis & cohom)
{
	omega_.resize(triangulation_->NumberOfTriangles());

	// Take the dual form from triangles t to t' to be the form integrated from the 
	// vertex at position 0 of t to the vertex at position 0 of t'.
	
	for(int triangleId1=0;triangleId1<triangulation_->NumberOfTriangles();triangleId1++)
	{
		Triangle * triangle1 = triangulation_->getTriangle(triangleId1);
		for(int direction1=0;direction1<3;direction1++)
		{
			Edge * edge = triangle1->getEdge(direction1);
			Triangle * triangle2 = edge->getAdjacent()->getParent();
			int direction2 = edge->getAdjacent()->getId();
			for(int x=0;x<2;x++)
			{
				setOmega(edge,x, (direction1==0 ? cohom.getOmega(triangle1->getEdge(2),x) : (direction1 == 1 ? - cohom.getOmega(triangle1->getEdge(1),x) : 0) )
							    +(direction2==0 ? cohom.getOmega(triangle2->getEdge(1),x) : (direction2 == 1 ? 0 : - cohom.getOmega(triangle2->getEdge(2),x)) ) );
			}
		}
	}
	if( cohom.IsUpToDate() )
	{
		SetUpToDate();
	}
}

void DualCohomologyBasis::UpdateAfterFlipMove(const Edge * const edge)
{
	Edge * otherEdge = edge->getPrevious()->getAdjacent()->getNext();
	
	addToOmega(otherEdge->getNext()->getAdjacent(),getOmega(otherEdge));
	setOmegaToMinusAdjacent(otherEdge->getNext());
	
	addToOmega(edge->getAdjacent(),getOmega(otherEdge));
	setOmegaToMinusAdjacent(edge);

	setOmegaToMinusAdjacent(otherEdge);
	IntForm2D zero = {0,0};
	setOmega(edge->getPrevious(),zero);
	setOmega(otherEdge->getPrevious(),zero);

	SetUpToDate();
}

bool DualCohomologyBasis::CheckClosedness() const
{
	for(int i=0;i<triangulation_->NumberOfVertices();i++)
	{
		Vertex * vertex = triangulation_->getVertex(i);
		for(int j=0;j<2;j++)
		{
			int total = 0;
			Edge * startEdge = vertex->getParent()->getNext();
			Edge * edge = startEdge;
			do
			{
				total += getOmega(edge,j);
				edge = edge->getAdjacent()->getPrevious();
			} while( edge != startEdge );
			if( total != 0 )
				return false;
		}
	}
	return true;
}

IntForm2D DualCohomologyBasis::IntegrateToParent(Edge * edge) const
{
	// Consider the path starting at edge->getParent() and ending at the parent of the end-point of edge in anti-clockwise direction.
	// Return the integral of dualOmega along this path.

	IntForm2D integral = {0,0};
	while( edge->getPrevious() != edge->getPrevious()->getOpposite()->getParent() )
	{
		integral = AddForms(integral,getOmega(edge));
		edge = edge->getAdjacent()->getPrevious();
	}
	return integral;
}

void DualCohomologyBasis::Simplify(bool StayInSameClass)
{
	CohomologyBasis cohomologybasis(triangulation_,*this);

	ShortestLoop shortestloop(triangulation_, &cohomologybasis );
	shortestloop.FindGenerators();
	std::vector<std::list<Edge*> > generators = shortestloop.getGenerators();
	std::vector<IntForm2D> integrals;
	if( StayInSameClass )
	{
		integrals = shortestloop.getGeneratorIntegrals();
	} else
	{
		integrals.resize(2);
		integrals[0][0] = 1;
		integrals[0][1] = 0;
		integrals[1][0] = 0;
		integrals[1][1] = 1;
	}

	SetAccordingToGenerators(generators,integrals);
}