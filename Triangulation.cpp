#include "triangulation.h"
#include "Edge.h"
#include <time.h>

Triangulation::Triangulation(void) : use_flipmove_(true)
{
	rng_.seed(static_cast<unsigned int>(time(NULL)));
}


Triangulation::~Triangulation(void)
{
	for(int i=0;i<(int)triangles_.size();i++)
	{
		delete triangles_[i];
	}
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
}

void Triangulation::DoSweep()
{
	// First perform a sweep of triangle flips
	int SuccesfulMoves = 0;
	while( SuccesfulMoves < n_triangles_ )
	{
		if( TryFlipMove() )
			SuccesfulMoves++;
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

	if( !randomEdge->IsFlipMovePossible() )
		return false;

	for(std::list<Matter *>::iterator matter = matter_.begin(); matter != matter_.end(); matter++ )
	{
		if( !(*matter)->IsFlipMoveAllowed(randomEdge) )
			return false;
	}

	double BoltzmannChange = 1.0;
	for(std::list<Matter *>::iterator matter = matter_.begin(); matter != matter_.end(); matter++ )
	{
		BoltzmannChange *= (*matter)->BoltzmannChangeUnderFlipMove(randomEdge);
	}

	if( !SucceedWithProbability(BoltzmannChange) )
	{
		return false;
	}

	randomEdge->DoFlipMove();
	for(std::list<Decoration *>::iterator decoration = decoration_.begin(); decoration != decoration_.end(); decoration++ )
	{
		(*decoration)->UpdateAfterFlipMove(randomEdge);
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
				vertices_.push_back(new Vertex(vertices_.size(),edge));

				// Walk around the new vertex and make the appropriate edge point to the vertex
				do {
					edge = edge->getNext()->getAdjacent()->getNext();
					edge->setOpposite(vertices_.back());
				} while( edge != (*tri)->getEdge(i) );
			}
		}
	}
}