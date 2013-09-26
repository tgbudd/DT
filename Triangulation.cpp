#include "boost/random/uniform_int.hpp"
#include "boost/random/uniform_real.hpp"

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

int Triangulation::RandomInteger(int min, int max)
{
	boost::uniform_int<> distribution(min, max);
	return distribution(rng_);
}

bool Triangulation::SucceedWithProbability(double probability)
{
	if( probability < 1.0e-8 )
		return false;
	if( 1.0 - probability < 1.0e-8 )
		return true;
	return probability > RandomReal(0.0,1.0); 
}

double Triangulation::RandomReal()
{
	return RandomReal(0.0,1.0);
}
double Triangulation::RandomReal(double min, double max)
{
	boost::uniform_real<> distribution(min,max);
	return distribution(rng_);
}

void Triangulation::setDominantMatter(DominantMatter * const & dominantmatter)
{
	dominantmatter_ = dominantmatter;
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