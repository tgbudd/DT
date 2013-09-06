#pragma once

#include <map>

#include "triangulation.h"
#include "CohomologyBasis.h"
#include "embedding.h"
#include "ConjugateGradient.h"

class LaplacianMatrix : public Matrix
{
public:
	LaplacianMatrix (Triangulation * const triangulation_) {
		laplacianRules_.resize(triangulation_->NumberOfVertices());
		for(int i=0;i<triangulation_->NumberOfTriangles();i++)
		{
			Triangle * triangle = triangulation_->getTriangle(i);
			for(int j=0;j<3;j++)
			{
				Edge * edge = triangle->getEdge(j);
				int start = edge->getNext()->getOpposite()->getId();
				int end = edge->getPrevious()->getOpposite()->getId();

				std::pair<std::map<int,int>::iterator,bool> insertion;
				insertion = laplacianRules_[start].insert(std::pair<int,int>(end,-1));
				if( !insertion.second )
				{
					insertion.first->second -= 1;
				}
				insertion = laplacianRules_[start].insert(std::pair<int,int>(start,1));
				if( !insertion.second )
				{
					insertion.first->second += 1;
				}
			}
		}
	}

	void MultiplyVector(const std::vector<double> & from, std::vector<double> & to) const
	{
		for(int i=0;i<from.size();i++)
		{
			to[i] = 0.0;
			for(std::map<int,int>::const_iterator it = laplacianRules_[i].begin(); it != laplacianRules_[i].end(); it++)
			{
				to[i] += it->second * from[it->first];
			}
		}
	}
private:
	std::vector<std::map<int,int> > laplacianRules_;
};

class HarmonicEmbedding :
	public Embedding
{
public:
	HarmonicEmbedding(Triangulation * const triangulation, const CohomologyBasis * const cohomologybasis);
	~HarmonicEmbedding() {}

	bool FindEmbedding();

	std::pair< double, double > CalculateModuli();

	void setMaxIterations(int iterations)
	{
		maxiterations_ = iterations;
	}
private:
	void LoadInitialCoordinates( std::vector<double> & coordinates, int i, Vertex * startVertex ) const;

	Triangulation * triangulation_;
	const CohomologyBasis * const cohomologybasis_;

	double accuracy_;
	int maxiterations_;
};

