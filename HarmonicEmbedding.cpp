#include "HarmonicEmbedding.h"
#include <vector>
#include <algorithm>

HarmonicEmbedding::HarmonicEmbedding(Triangulation * const triangulation, const CohomologyBasis * const cohomologybasis) 
	: triangulation_(triangulation), cohomologybasis_(cohomologybasis)
{
	accuracy_ = 1.0e-8;
	maxiterations_ = 120;
}

bool HarmonicEmbedding::FindEmbedding()
{
	setSize(triangulation_->NumberOfTriangles(),triangulation_->NumberOfVertices());

	LaplacianMatrix laplacian(triangulation_);

	std::vector<double> deltaOmega(triangulation_->NumberOfVertices(),0.0);
	std::vector<double> coordinate(triangulation_->NumberOfVertices(),0.0);
	for(int i=0;i<2;i++)
	{
		if( i > 0 )
		{
			std::fill(deltaOmega.begin(),deltaOmega.end(),0.0);
		}

		for(int j=0;j<triangulation_->NumberOfTriangles();j++)
		{
			Triangle * triangle = triangulation_->getTriangle(j);
			for(int k=0;k<3;k++)
			{
				Edge * edge = triangle->getEdge(k);
				deltaOmega[edge->getPrevious()->getOpposite()->getId()] += (double)(cohomologybasis_->getOmega(edge)[i]);
			}
		}

		// it can happen that deltaOmega is identically zero; in that case x[i]=0 is the correct solution
		bool iszero = true;
		for(int j=0;j<triangulation_->NumberOfVertices();j++)
		{
			if( fabs(deltaOmega[j]) > 0.01 )
			{
				iszero = false;
				break;
			}
		}


		if( iszero )
		{
			for(int j=0;j<triangulation_->NumberOfVertices();j++)
			{
				coordinate[j] = 0.0;
			}
		}else
		{
			for(int j=0;j<triangulation_->NumberOfVertices();j++)
			{
				coordinate[j] = getCoordinate(j)[i];
			}
			if( !laplacian.ConjugateGradientSolve(deltaOmega,coordinate,accuracy_,maxiterations_) )
			{
				return false;
			}
		}

		for(int j=0;j<triangulation_->NumberOfTriangles();j++)
		{
			Triangle * triangle = triangulation_->getTriangle(j);
			for(int k=0;k<3;k++)
			{
				Edge * edge = triangle->getEdge(k);
				setForm(edge,i,(double)(cohomologybasis_->getOmega(edge)[i]) + coordinate[edge->getNext()->getOpposite()->getId()] - coordinate[edge->getPrevious()->getOpposite()->getId()]); 
			}
		}

		for(int j=0;j<triangulation_->NumberOfVertices();j++)
		{
			setCoordinate(j,i,properfmod(coordinate[j],1.0));
		}
	}
	return true;
}

std::pair< double, double > HarmonicEmbedding::CalculateModuli()
{
	boost::array<boost::array<double,2>, 2> inproducts;
	for(int i=0;i<2;i++)
	{
		for(int j=0;j<2;j++)
		{
			inproducts[i][j] = 0.0;
			for(int k=0;k<triangulation_->NumberOfTriangles();k++)
			{
				for(int l=0;l<3;l++)
				{
					inproducts[i][j] += 0.5*getForm(k,l)[i]*getForm(k,l)[j];
				}
			}
		}
	}
	return std::pair<double,double>( -inproducts[0][1] / inproducts[1][1], sqrt(inproducts[0][0] * inproducts[1][1] - inproducts[0][1] * inproducts[1][0])/inproducts[1][1] );
}
