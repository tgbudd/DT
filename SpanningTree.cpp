#include <queue>

#include "SpanningTree.h"


SpanningTree::SpanningTree(Triangulation * const triangulation)
	: triangulation_(triangulation)
{
}


SpanningTree::~SpanningTree()
{
}

void SpanningTree::Initialize(const std::vector<boost::array<bool,3> > & intree)
{
	in_spanning_tree_ = intree;
	UpdateSpanningDegree();
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
	int MovesPerSweep = (UsesCustomSweepSize() ? CustomSweepSize() : triangulation_->NumberOfTriangles() );
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
	if( !UsesCustomSweepSize() && triangulation_->SucceedWithProbability(0.02) )
	{
		InitializeWithLoopErasedRandomWalk();
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

void SpanningTree::UpdateSpanningDegree(const Vertex * vertex, int degree)
{
	spanning_degree_[vertex->getId()] = degree;
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
	Vertex * midVertex = edge->getNext()->getOpposite();
	int degree = spanning_degree_[midVertex->getId()];
	if( degree < 2 )
	{
		return false;
	}
	
	std::pair<Edge *,int> otherEdge = RandomOtherSpanningEdge( edge );

	Vertex * firstVertex = edge->getPrevious()->getOpposite();
	Vertex * secondVertex = otherEdge.first->getPrevious()->getOpposite();
	if( firstVertex == secondVertex || firstVertex == midVertex )
	{
		return false;
	}

	// After the move the degree of the vertex will be newdegree
	int newdegree = spanning_degree_[ firstVertex->getId() ] + spanning_degree_[ secondVertex->getId() ];

	boost::array<Edge*, 2> edges = {edge,otherEdge.first};
	if( !triangulation_->TryCutMove( edges, static_cast<double>(degree-1)/static_cast<double>(newdegree-1) ) )
	{
		return false;
	}

	UpdateSpanningDegree( firstVertex, newdegree );
	UpdateSpanningDegree( secondVertex, otherEdge.second );
	UpdateSpanningDegree( midVertex, degree - otherEdge.second );
	
	return true;
}

std::pair<Edge *,int> SpanningTree::RandomOtherSpanningEdge(Edge * edge)
{
	int direction = triangulation_->RandomInteger(1,spanning_degree_[edge->getNext()->getOpposite()->getId()]-1);
	Edge * otherEdge = edge;
	for(int degree = 0;degree < direction; degree++)
	{
		do {
			otherEdge = otherEdge->getPrevious()->getAdjacent();
		} while( !InSpanningTree(otherEdge) );
	}
	return std::pair<Edge *, int>(otherEdge,direction);
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

void SpanningTree::InitializeWithLoopErasedRandomWalk()
{
	// This is Wilson's algorithm.

	boost::array<bool,3> alltrue = {true,true,true};
	in_spanning_tree_.clear();
	in_spanning_tree_.resize(triangulation_->NumberOfTriangles(),alltrue);

	const int IN_TREE = 3;
	std::vector<int> status(triangulation_->NumberOfTriangles(),0);

	std::vector<Triangle *> randomNotInTree;
	std::vector<int> positionInRandom;
	randomNotInTree.reserve(triangulation_->NumberOfTriangles());
	positionInRandom.reserve(triangulation_->NumberOfTriangles());
	for(int i=0,end=triangulation_->NumberOfTriangles();i<end;i++)
	{
		randomNotInTree.push_back(triangulation_->getTriangle(i));
		positionInRandom.push_back(i);
	}
	
	int rootindex = triangulation_->RandomInteger(0,randomNotInTree.size()-1);
	Triangle * root = randomNotInTree[rootindex];
	status[root->getId()] = IN_TREE;
	positionInRandom[randomNotInTree.back()->getId()] = rootindex;
	std::swap(randomNotInTree[rootindex],randomNotInTree.back());
	randomNotInTree.pop_back();

	while( !randomNotInTree.empty() )
	{
		int startindex = triangulation_->RandomInteger(0,randomNotInTree.size()-1);

		// Perform random walk until it hits the tree.
		Triangle * currentTriangle = randomNotInTree[startindex];
		while( status[currentTriangle->getId()] != IN_TREE )
		{
			int direction = triangulation_->RandomInteger(0,2);
			status[currentTriangle->getId()] = direction;
			currentTriangle = currentTriangle->getEdge(direction)->getAdjacent()->getParent();
		}

		// Add the loop-erased walk to the tree.
		currentTriangle = randomNotInTree[startindex];
		while( status[currentTriangle->getId()] != IN_TREE )
		{
			positionInRandom[randomNotInTree.back()->getId()] = positionInRandom[currentTriangle->getId()];
			std::swap(randomNotInTree[positionInRandom[currentTriangle->getId()]],randomNotInTree.back());
			randomNotInTree.pop_back();
			
			int direction = status[currentTriangle->getId()];
			status[currentTriangle->getId()] = IN_TREE;
			setInSpanningTree(currentTriangle->getEdge(direction),false);
			setInSpanningTree(currentTriangle->getEdge(direction)->getAdjacent(),false);
			currentTriangle = currentTriangle->getEdge(direction)->getAdjacent()->getParent();
		}
	}
	UpdateSpanningDegree();

}

bool SpanningTree::CheckTree()
{
	// Check number of edges not in spanning tree (=number of edges in dual spanning tree)
	int num=0;
	for(int i=0;i<triangulation_->NumberOfTriangles();i++)
	{
		for(int j=0;j<3;j++)
		{
			if( !InSpanningTree(i,j) )
			{
				num++;
				if( InSpanningTree(triangulation_->getTriangle(i)->getEdge(j)->getAdjacent()) )
				{
					return false;
				}
			}
		}
	}
	if( num != triangulation_->NumberOfTriangles()-1 )
	{
		return false;
	}

	// Check it is connected.
	std::vector<bool> visited(triangulation_->NumberOfTriangles(),false);
	std::queue<Triangle *> q;
	q.push(triangulation_->getTriangle(0));
	visited[0]=true;
	num=1;
	while( !q.empty() )
	{
		Triangle * t = q.front();
		q.pop();
		for(int i=0;i<3;i++)
		{
			Triangle * nbr = t->getEdge(i)->getAdjacent()->getParent();
			if( !InSpanningTree(t->getEdge(i)) && !visited[ nbr->getId() ] )
			{
				visited[ nbr->getId() ] = true;
				num++;
				q.push(nbr);
			}
		}
	}
	if( num != triangulation_->NumberOfTriangles() )
	{
		return false;
	}
	return true;
}
