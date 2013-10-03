#include "boost/random/uniform_int.hpp"
#include "boost/random/uniform_real.hpp"
#include "boost/random/normal_distribution.hpp"

#include "Triangulation.h"
#include "Edge.h"
#include "Decoration.h"
#include "Matter.h"
#include "DominantMatter.h"
#include <time.h>

Triangulation::Triangulation(void) : use_flipmove_(true) , dominantmatter_(NULL)
{
	rng_.seed(static_cast<unsigned int>(time(NULL)));
	state_ = 1;
}

void Triangulation::SeedRandom(unsigned int seed)
{
	rng_.seed(seed);
}

Triangulation::~Triangulation(void)
{
	for(int i=0;i<(int)triangles_.size();i++)
	{
		delete triangles_[i];
	}
}

Triangle * const & Triangulation::getTriangle(int id) const
{
	return triangles_[id];
}

Vertex * const & Triangulation::getVertex(int id) const
{
	return vertices_[id];
}
	
Triangle * const & Triangulation::getRandomTriangle()
{
	return triangles_[RandomInteger(0,n_triangles_-1)];
}

Edge * const & Triangulation::getRandomEdge()
{
	return triangles_[RandomInteger(0,n_triangles_-1)]->getEdge(RandomInteger(0,2));
}

Vertex * const & Triangulation::getRandomVertex()
{
	return vertices_[RandomInteger(0,n_vertices_-1)];
}

int Triangulation::RandomInteger(int min, int max) const
{
	boost::uniform_int<> distribution(min, max);
	return distribution(rng_);
}

bool Triangulation::SucceedWithProbability(double probability) const
{
	if( probability < 1.0e-8 )
		return false;
	if( 1.0 - probability < 1.0e-8 )
		return true;
	return probability > RandomReal(0.0,1.0); 
}

double Triangulation::RandomReal() const
{
	return RandomReal(0.0,1.0);
}
double Triangulation::RandomReal(double min, double max) const
{
	boost::uniform_real<> distribution(min,max);
	return distribution(rng_);
}
double Triangulation::RandomNormal(double mean, double sigma) const
{
	boost::normal_distribution<> distribution(mean,sigma);
	return distribution(rng_);
}

void Triangulation::setDominantMatter(DominantMatter * const & dominantmatter)
{
	dominantmatter_ = dominantmatter;
}

void Triangulation::clearDominantMatter()
{
	dominantmatter_ = NULL;
}

void Triangulation::AddMatter(Matter * matter)
{
	decoration_.push_back(matter);
	matter_.push_back(matter);
}

void Triangulation::AddDecoration(Decoration * decoration)
{
	decoration_.push_back(decoration);
}

void Triangulation::LoadRegularLattice(int width, int height)
{
	n_triangles_ = width * height * 2;

	triangles_.reserve(n_triangles_);
	for(int i=0;i<n_triangles_;i++)
	{
		triangles_.push_back(new Triangle());
		triangles_.back()->setId(i);
	}

	for(int y=0;y<height;y++)
	{
		for(int x=0;x<width;x++)
		{
			int UpwardPointingTriangle = x + 2*y*width;
			// fix adjacency between UpwardPointingTriangle and the appropriate downward pointing triangles
			triangles_[UpwardPointingTriangle]->getEdge(0)->bindAdjacent(triangles_[x + (2*((y+height-1)%height)+1)*width]->getEdge(0));
			triangles_[UpwardPointingTriangle]->getEdge(1)->bindAdjacent(triangles_[(x+1)%width + (2*y+1)*width]->getEdge(1));
			triangles_[UpwardPointingTriangle]->getEdge(2)->bindAdjacent(triangles_[x + (2*y+1)*width]->getEdge(2));
		}
	}

	DetermineVertices();
	IncreaseState();
}

void Triangulation::DoSweep()
{
	if( dominantmatter_ != NULL )
	{
		dominantmatter_->DoSweep();
	} else
	{
		// Perform a sweep of triangle flips
		int SuccesfulMoves = 0;
		while( SuccesfulMoves < n_triangles_ )
		{
			if( TryFlipMove() )
				SuccesfulMoves++;
		}
	}

	// Then perform a sweep for each matter field
	for(std::list<Matter *>::iterator matter = matter_.begin(); matter != matter_.end(); matter++ )
	{
		(*matter)->DoSweep();
	}
}

void Triangulation::DoSweep(int NumberOfSweeps)
{
	for(int i=0;i<NumberOfSweeps;i++)
	{
		DoSweep();
	}
}

bool Triangulation::TryFlipMove()
{
	Edge * randomEdge = getRandomEdge();
	return TryFlipMove(randomEdge);
}

/// <image url="$(SolutionDir)\images\FlipMove.png" scale="1.2" />

bool Triangulation::TryFlipMove(Edge * edge)
{
	if( !edge->IsFlipMovePossible() )
		return false;

	for(std::list<Matter *>::iterator matter = matter_.begin(); matter != matter_.end(); matter++ )
	{
		if( !(*matter)->IsFlipMoveAllowed(edge) )
			return false;
	}

	double BoltzmannChange = 1.0;
	for(std::list<Matter *>::iterator matter = matter_.begin(); matter != matter_.end(); matter++ )
	{
		BoltzmannChange *= (*matter)->BoltzmannChangeUnderFlipMove(edge);
	}

	if( !SucceedWithProbability(BoltzmannChange) )
	{
		return false;
	}

	edge->DoFlipMove();
	IncreaseState();

	for(std::list<Decoration *>::iterator decoration = decoration_.begin(); decoration != decoration_.end(); decoration++ )
	{
		(*decoration)->UpdateAfterFlipMove(edge);
	}

	return true;
}

bool Triangulation::TryCutMove(const boost::array< Edge *, 2> & edges, double combinatorialBoltzmann)
{
	if( edges[0]->getNext()->getOpposite() != edges[1]->getNext()->getOpposite() ||
		edges[0]->getPrevious()->getOpposite() == edges[1]->getPrevious()->getOpposite() ||
		edges[0]->getPrevious()->getOpposite() == edges[0]->getNext()->getOpposite() ||
		edges[1]->getPrevious()->getOpposite() == edges[1]->getNext()->getOpposite() )
	{
		return false;
	}

	// Make explicit the dual bonds that are deleted and created
	std::vector<boost::array<Triangle *,2> > toBeDeleted(2);
	for(int i=0;i<2;i++)
	{
		toBeDeleted[i][0] = edges[i]->getParent();
		toBeDeleted[i][1] = edges[i]->getAdjacent()->getParent();
	}
	std::vector<boost::array<Triangle *,2> > toBeAdded(toBeDeleted);
	std::swap( toBeAdded[0][1], toBeAdded[1][1] );

	double BoltzmannChange = combinatorialBoltzmann;
	for(std::list<Matter *>::iterator matter = matter_.begin(); matter != matter_.end(); matter++ )
	{
		BoltzmannChange *= (*matter)->BoltzmannChangeUnderGeneralMove(toBeDeleted,toBeAdded);
	}

	if( !SucceedWithProbability(BoltzmannChange) )
	{
		return false;
	}

	DoCutMove(edges);
	IncreaseState();

	for(std::list<Decoration *>::iterator decoration = decoration_.begin(); decoration != decoration_.end(); decoration++ )
	{
		(*decoration)->UpdateAfterCutMove(edges);
	}


	return true;

}

void Triangulation::DoCutMove(const boost::array< Edge *, 2> & edges)
{
	Vertex * midVertex = edges[0]->getNext()->getOpposite();
	Vertex * firstVertex = edges[0]->getPrevious()->getOpposite();
	Vertex * secondVertex = edges[1]->getPrevious()->getOpposite();

	// Update vertices
	secondVertex->setParent(edges[0]->getNext());
	midVertex->setParent(edges[1]->getNext());

	Edge * oppositeEdge = edges[0]->getNext();
	while( oppositeEdge->getPrevious() != edges[1] )
	{
		oppositeEdge->setOpposite( secondVertex );
		oppositeEdge = oppositeEdge->getNext()->getAdjacent()->getNext();
	}
	oppositeEdge = edges[1]->getAdjacent()->getNext();
	do {
		oppositeEdge->setOpposite( firstVertex );
		oppositeEdge = oppositeEdge->getNext()->getAdjacent()->getNext();
	} while( oppositeEdge->getPrevious()->getAdjacent() != edges[1] );

	// Update edge adjacency
	Edge * adjFirstEdge = edges[0]->getAdjacent();
	edges[0]->bindAdjacent(edges[1]->getAdjacent());
	edges[1]->bindAdjacent(adjFirstEdge);

	BOOST_ASSERT( midVertex->getParent()->getOpposite() == midVertex );
	BOOST_ASSERT( firstVertex->getParent()->getOpposite() == firstVertex );
	BOOST_ASSERT( secondVertex->getParent()->getOpposite() == secondVertex );
	BOOST_ASSERT( CheckVertexNeighbourhood(midVertex) );
	BOOST_ASSERT( CheckVertexNeighbourhood(firstVertex) );
	BOOST_ASSERT( CheckVertexNeighbourhood(secondVertex) );
}

void Triangulation::DetermineVertices()
{
	if( !vertices_.empty() )
	{
		for(std::vector<Vertex *>::iterator vert = vertices_.begin();vert != vertices_.end();vert++)
		{
			delete *vert;
		}
	}

	vertices_.clear();

	for(std::vector<Triangle *>::iterator tri = triangles_.begin(); tri != triangles_.end(); tri++)
	{
		for(int i=0;i<3;i++)
		{
			(*tri)->getEdge(i)->setOpposite(NULL);
		}
	}

	for(std::vector<Triangle *>::iterator tri = triangles_.begin(); tri != triangles_.end(); tri++)
	{
		for(int i=0;i<3;i++)
		{
			Edge * edge = (*tri)->getEdge(i);
			if( edge->getOpposite() == NULL )
			{
				// Create new vertex and make this edge its parent
				vertices_.push_back(new Vertex(static_cast<int>(vertices_.size()),edge));

				// Walk around the new vertex and make the appropriate edge point to the vertex
				do {
					edge = edge->getNext()->getAdjacent()->getNext();
					edge->setOpposite(vertices_.back());
				} while( edge != (*tri)->getEdge(i) );
			}
		}
	}
	n_vertices_ = static_cast<int>(vertices_.size());
}

void Triangulation::LoadFromAdjacencyList(const std::vector<boost::array<std::pair<int,int>,3 > > & adj)
{
	n_triangles_ = static_cast<int>(adj.size());

	triangles_.reserve(n_triangles_);
	for(int i=0;i<n_triangles_;i++)
	{
		triangles_.push_back(new Triangle());
		triangles_.back()->setId(i);
	}

	for(int i=0, end=static_cast<int>(adj.size());i<end;i++)
	{
		for(int j=0;j<3;j++)
		{
			triangles_[i]->getEdge(j)->setAdjacent(triangles_[adj[i][j].first]->getEdge(adj[i][j].second));
		}
	}

	DetermineVertices();
}

std::string Triangulation::OutputData() const
{
	std::ostringstream stream;
	stream << "triangulation -> {numberoftriangles -> " << NumberOfTriangles() << "}";
	return stream.str();
}

void Triangulation::IncreaseState()
{
	state_++;
}

Triangle * Triangulation::NewTriangle()
{
	triangles_.push_back(new Triangle());
	triangles_.back()->setId(static_cast<int>(triangles_.size())-1);
	n_triangles_++;
	return triangles_.back();
}

void Triangulation::Clear()
{
	for(std::vector<Triangle *>::iterator it = triangles_.begin(); it != triangles_.end(); ++it )
	{
		delete *it;
	}
	triangles_.clear();
	n_triangles_ = 0;
	for(std::vector<Vertex *>::iterator it = vertices_.begin(); it != vertices_.end(); ++it )
	{
		delete *it;
	}
	vertices_.clear();
	n_vertices_ = 0;
	IncreaseState();
}

bool Triangulation::CheckVertexNeighbourhood(const Vertex * const vertex) const
{
	Edge * edge = vertex->getParent();
	do {
		if( edge->getOpposite() != vertex )
		{
			return false;
		}
		edge = edge->getNext()->getAdjacent()->getNext();
	} while( edge != vertex->getParent() );

	return true;
}