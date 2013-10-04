#include "SpanningTree.h"


SpanningTree::SpanningTree(Triangulation * const triangulation)
	: triangulation_(triangulation)
{
}


SpanningTree::~SpanningTree()
{
}

void SpanningTree::Initialize()
{
	// Grow a spanning tree on the dual graph.
	boost::array<bool,3> alltrue = {true,true,true};
	in_spanning_tree_.clear();
	in_spanning_tree_.resize(triangulation_->NumberOfTriangles(),alltrue);

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
			if( !visited[edge->getAdjacent()->getParent()->getId()] )
			{
				setInSpanningTree(edge, false);
				setInSpanningTree(edge->getAdjacent(),false);
				visited[edge->getAdjacent()->getParent()->getId()] = true;
				queue.push_back(edge->getAdjacent()->getParent());
			}
		}
	}

	UpdateSpanningDegree();
}

void SpanningTree::DoSweep()
{
	int SuccesfulMoves = 0;
	int MovesPerSweep = triangulation_->NumberOfTriangles();
	while( SuccesfulMoves < MovesPerSweep)
	{
		Edge * edge = triangulation_->getRandomEdge();
		if( InSpanningTree( edge ) )
		{
			if( TrySpanningMove( edge ) )
				SuccesfulMoves++;
		} else
		{
			if( TryFlipMove( edge ) )
				SuccesfulMoves++;
		}
	}
}

bool SpanningTree::InSpanningTree(int triangle,int edge) const
{
	return in_spanning_tree_[triangle][edge];
}

bool SpanningTree::InSpanningTree(Edge * edge) const
{
	return InSpanningTree(edge->getParent()->getId(),edge->getId());
}

void SpanningTree::setInSpanningTree(int triangle, int edge, bool ist)
{
	in_spanning_tree_[triangle][edge] = ist;
}

void SpanningTree::setInSpanningTree(Edge * edge, bool ist)
{
	setInSpanningTree(edge->getParent()->getId(),edge->getId(),ist);
}

int SpanningTree::SpanningDegree(const Vertex * vertex) const
{
	Edge * edge = vertex->getParent()->getPrevious();
	int degree = 0;
	do {
		if( InSpanningTree(edge) )
		{
			degree++;
		}
		edge = edge->getPrevious()->getAdjacent();
	} while( edge->getNext() != vertex->getParent() );
	return degree;
}

void SpanningTree::UpdateSpanningDegree(const Vertex * vertex)
{
	spanning_degree_[vertex->getId()] = SpanningDegree(vertex);
}

void SpanningTree::UpdateSpanningDegree()
{
	if( spanning_degree_.size() != triangulation_->NumberOfVertices() )
	{
		spanning_degree_.resize( triangulation_->NumberOfVertices() );
	}
	for(int i=0,end=triangulation_->NumberOfVertices();i<end;i++)
	{
		UpdateSpanningDegree(triangulation_->getVertex(i));
	}
}

bool SpanningTree::TryFlipMove(Edge * edge)
{
	BOOST_ASSERT( !InSpanningTree( edge ) );
	
	Edge * adjEdge = edge->getAdjacent();
	if( !triangulation_->TryFlipMove( edge ) )
	{
		return false;
	}

	setInSpanningTree( edge, InSpanningTree( adjEdge->getPrevious() ) );
	setInSpanningTree( adjEdge, InSpanningTree( edge->getPrevious() ) );
	setInSpanningTree( edge->getPrevious(), false );
	setInSpanningTree( adjEdge->getPrevious(), false );

	return true;
}

bool SpanningTree::TrySpanningMove(Edge * edge)
{
	Vertex * vertex = edge->getNext()->getOpposite();
	int degree = spanning_degree_[vertex->getId()];
	if( degree < 2 )
	{
		return false;
	}

	Edge * otherEdge = RandomOtherSpanningEdge( edge );

	boost::array<Edge*, 2> edges = {edge,otherEdge};
	if( edges[0]->getPrevious()->getOpposite() == edges[1]->getPrevious()->getOpposite() || edges[0]->getPrevious()->getOpposite() == vertex )
	{
		return false;
	}

	boost::array<int,2> degrees;
	for(int i=0;i<2;i++)
	{
		degrees[i] = spanning_degree_[ edges[i]->getPrevious()->getOpposite()->getId() ];
	}

	// After the move the degree of the vertex will be newdegree
	int newdegree = degrees[0] + degrees[1];

	if( !triangulation_->TryCutMove( edges, static_cast<double>(degree-1)/static_cast<double>(newdegree-1) ) )
	{
		return false;
	}

	UpdateSpanningDegree(edges[0]->getPrevious()->getOpposite());
	UpdateSpanningDegree(edges[0]->getNext()->getOpposite());
	UpdateSpanningDegree(edges[1]->getNext()->getOpposite());

	return true;
}

Edge * SpanningTree::RandomOtherSpanningEdge(Edge * edge)
{
	int direction = triangulation_->RandomInteger(1,spanning_degree_[edge->getNext()->getOpposite()->getId()]-1);
	Edge * otherEdge = edge;
	for(int degree = 0;degree < direction; degree++)
	{
		do {
			otherEdge = otherEdge->getPrevious()->getAdjacent();
		} while( !InSpanningTree(otherEdge) );
	}
	return otherEdge;
}

double SpanningTree::CentralCharge() const
{
	return -2.0;
}

std::string SpanningTree::ConfigurationData() const
{
	std::ostringstream stream;
	stream << std::fixed << "{type -> \"spanningtree\", centralcharge -> -2 }";
	return stream.str();
}