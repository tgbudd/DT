#include "CohomologyBasis.h"
#include "ShortestLoop.h"
#include "DualCohomologyBasis.h"

CohomologyBasis::CohomologyBasis(Triangulation * const triangulation) : triangulation_(triangulation)
{
}


CohomologyBasis::~CohomologyBasis(void)
{
}

CohomologyBasis::CohomologyBasis(Triangulation * const triangulation, const DualCohomologyBasis & dualcohomologybasis) : triangulation_(triangulation)
{
	omega_.resize(triangulation_->NumberOfTriangles());
	SetToDualOf(dualcohomologybasis);
}

void CohomologyBasis::Initialize()
{
	omega_.resize(triangulation_->NumberOfTriangles());
	ClearOmega();
}

void CohomologyBasis::ClearOmega()
{
	IntForm2D zeroForm = {0,0};
	boost::array<IntForm2D,3> zeroTriangle = {zeroForm,zeroForm,zeroForm};
	std::fill(omega_.begin(),omega_.end(),zeroTriangle);
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

	BOOST_ASSERT(CheckClosedness());
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

	/*// test
	for(int i=0;i<2;i++)
	{
		int tot=getOmega(edge,i);
		tot += getOmega(edge->getNext(),i);
		tot += getOmega(edge->getPrevious(),i);
		BOOST_ASSERT( tot == 0 );
	}*/
}

bool CohomologyBasis::CheckClosedness() const
{
	for(int i=0;i<triangulation_->NumberOfTriangles();i++)
	{
		for(int j=0;j<2;j++)
		{
			int total = 0;
			for(int k=0;k<3;k++)
			{
				total += omega_[i][k][j];
			}
			if( total != 0 )
				return false;
		}
	}
	return true;
}

IntForm2D CohomologyBasis::Integrate(const std::list<Edge *> & path) const
{
	IntForm2D integral = {0,0};
	for(std::list<Edge *>::const_iterator edgeIt = path.begin(); edgeIt != path.end(); edgeIt++ )
	{
		integral = AddForms(integral,getOmega(*edgeIt));
	}
	return integral;
}

void CohomologyBasis::Simplify()
{
	ShortestLoop shortestloop(triangulation_,this);
	shortestloop.FindGenerators();
	std::vector<std::list<Edge*> > generators = shortestloop.getGenerators();
	std::vector<IntForm2D> integrals = shortestloop.getGeneratorIntegrals();

	DualCohomologyBasis dualbasis(triangulation_,generators,integrals);

	SetToDualOf(dualbasis);
}



void CohomologyBasis::SetToDualOf(const DualCohomologyBasis & dualOmega)
{
	// Take the form from vertex v to vertex v' to be the dual form integrated from the parent of v to the parent of v'
	for(int triangleId=0;triangleId<triangulation_->NumberOfTriangles();triangleId++)
	{
		Triangle * triangle = triangulation_->getTriangle(triangleId);
		for(int direction=0;direction<3;direction++)
		{
			Edge * edge = triangle->getEdge(direction);
			IntForm2D toInitial = dualOmega.IntegrateToParent(edge->getPrevious());
			IntForm2D toFinal = dualOmega.IntegrateToParent(edge);
			setOmega(edge,SubtractForms(toFinal,toInitial));
		}
	}
}