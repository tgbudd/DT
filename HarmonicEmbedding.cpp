#include "HarmonicEmbedding.h"
#include <vector>
#include <algorithm>
#include <queue>
#include <math.h>

HarmonicEmbedding::HarmonicEmbedding(Triangulation * const triangulation, CohomologyBasis * const cohomologybasis) 
	: triangulation_(triangulation), Embedding(triangulation,cohomologybasis)
{

}

bool HarmonicEmbedding::FindEdgeMeasure()
{
	for(int i=0, end=triangulation_->NumberOfTriangles();i<end;i++)
	{
		for(int j=0;j<3;j++)
		{
			setEdgeMeasure(i,j,1.0);
		}
	}
	return true;
}