#include "CohomologyBasis.h"


CohomologyBasis::CohomologyBasis(Triangulation * const triangulation) : triangulation_(triangulation)
{
}


CohomologyBasis::~CohomologyBasis(void)
{
}

void CohomologyBasis::Initialize()
{
	IntForm2D zeroForm = {0,0};
	boost::array<IntForm2D,3> zeroTriangle = {zeroForm,zeroForm,zeroForm};
	omega_.assign(triangulation_->NumberOfTriangles(),zeroTriangle);
}

void CohomologyBasis::Initialize(int width, int height)
{
	Initialize();
		
	BOOST_ASSERT( triangulation_->NumberOfTriangles() == 2 * width * height );

	for(int i=0;i<triangulation_->NumberOfVertices();i++)
	{
		BOOST_ASSERT(triangulation_->getVertex(i)->getDegree() == 6);
	}
	for(int y=0;y<height;y++)
	{
		omega_[2*y*width][0][0] = 1;
		omega_[2*y*width][2][0] = -1;
		omega_[(2*y+1)*width][0][0] = -1;
		omega_[(2*y+1)*width][2][0] = 1;
	}
	for(int x=0;x<width;x++)
	{
		omega_[x][1][1] = 1;
		omega_[x][2][1] = -1;
		omega_[x+width][1][1] = -1;
		omega_[x+width][2][1] = 1;
	}
}

void CohomologyBasis::UpdateAfterFlipMove(const Edge * const edge)
{
	setOmegaToMinusAdjacent(edge);
	setOmegaToMinusAdjacent(edge->getPrevious()->getAdjacent()->getNext());
	for(int i=0;i<2;i++)
	{
		int newomega = - getOmega(edge,i) - getOmega(edge->getNext(),i);
		BOOST_ASSERT( abs(newomega) < 100000000 );
		setOmega(edge->getPrevious(),i,newomega);
		setOmega(edge->getPrevious()->getAdjacent(),i,-newomega);
	}
}

