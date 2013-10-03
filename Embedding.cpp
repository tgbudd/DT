#include <queue>
#include <vector>
#include <map>

#include "Embedding.h"

LaplacianMatrix::LaplacianMatrix(const Triangulation * const triangulation) : Matrix(triangulation->NumberOfVertices() ) {
	boost::array<double,3> unitmeasure = {1.0,1.0,1.0};
	std::vector<boost::array<double,3> > edgemeasure(triangulation->NumberOfTriangles(),unitmeasure);
	Initialize(triangulation,edgemeasure);
}

LaplacianMatrix::LaplacianMatrix(const Triangulation * const triangulation, const std::vector<boost::array<double,3> > & edge_measure) : Matrix(edge_measure.size()) {
	Initialize(triangulation,edge_measure);
}

void LaplacianMatrix::Initialize(const Triangulation * const triangulation, const std::vector<boost::array<double,3> > & edge_measure)
{
	laplacianRules_.resize(triangulation->NumberOfVertices());
	for(int i=0;i<triangulation->NumberOfTriangles();i++)
	{
		Triangle * triangle = triangulation->getTriangle(i);
		for(int j=0;j<3;j++)
		{
			Edge * edge = triangle->getEdge(j);
			int start = edge->getNext()->getOpposite()->getId();
			int end = edge->getPrevious()->getOpposite()->getId();

			std::pair<std::map<int,double>::iterator,bool> insertion;
			insertion = laplacianRules_[start].insert(std::pair<int,double>(end,-edge_measure[i][j]));
			if( !insertion.second )
			{
				insertion.first->second -= edge_measure[i][j];
			}
			insertion = laplacianRules_[start].insert(std::pair<int,double>(start,edge_measure[i][j]));
			if( !insertion.second )
			{
				insertion.first->second += edge_measure[i][j];
			}
		}
	}
}

void LaplacianMatrix::MultiplyVector(const std::vector<double> & from, std::vector<double> & to) const
{
	for(int i=0;i<static_cast<int>(from.size());i++)
	{
		to[i] = 0.0;
		for(std::map<int,double>::const_iterator it = laplacianRules_[i].begin(); it != laplacianRules_[i].end(); it++)
		{
			to[i] += it->second * from[it->first];
		}
	}
}


Embedding::Embedding(const Triangulation * const triangulation, CohomologyBasis * const cohomologybasis) : Decoration(triangulation), triangulation_(triangulation), cohomologybasis_(cohomologybasis) {
	accuracy_ = 1.0e-6;
	maxiterations_ = 2000;
	edge_measure_.resize(triangulation->NumberOfTriangles());
}

void Embedding::UpdateAfterFlipMove(const Edge * const edge) 
{
	setFormToMinusAdjacent(edge);
	setFormToMinusAdjacent(edge->getPrevious()->getAdjacent()->getNext());
	Vector2D newform = NegateVector2D(AddVectors2D(getForm(edge),getForm(edge->getNext())));
	setForm(edge->getPrevious(),newform);
	setForm(edge->getPrevious()->getAdjacent(),NegateVector2D(newform));
}

void Embedding::Coderivative( const std::vector<boost::array<double,3> > & oneform, std::vector<double> & result )
{
	result.resize( triangulation_->NumberOfVertices() );
	std::fill(result.begin(),result.end(),0.0);
	for(int j=0;j<triangulation_->NumberOfTriangles();j++)
	{
		Triangle * triangle = triangulation_->getTriangle(j);
		for(int k=0;k<3;k++)
		{
			Edge * edge = triangle->getEdge(k);
			result[edge->getPrevious()->getOpposite()->getId()] += edge_measure_[j][k] * oneform[j][k];
		}
	}
}


void Embedding::LoadInitialCoordinates( std::vector<double> & coordinates, int i, Vertex * startVertex ) const
{
	std::queue<Vertex *> q;
	q.push(startVertex);
	std::vector<bool> visited(triangulation_->NumberOfVertices(),false);
	visited[startVertex->getId()] = true;
	coordinates[startVertex->getId()] = getCoordinate(startVertex)[i];
	while( !q.empty() )
	{
		Vertex * v = q.front();
		q.pop();
		Edge * startEdge = v->getParent()->getPrevious();
		Edge * edge = startEdge;
		do {
			Vertex * nbrVertex = edge->getPrevious()->getOpposite();
			if( !visited[nbrVertex->getId()] )
			{
				visited[nbrVertex->getId()] = true;
				q.push(nbrVertex);
				coordinates[nbrVertex->getId()] = coordinates[v->getId()] + getForm(edge)[i] - cohomologybasis_->getOmega(edge,i);
			}
			edge = edge->getAdjacent()->getNext();
		} while( edge != startEdge );
	}
}


bool Embedding::FindEmbedding()
{
	if( !cohomologybasis_->IsUpToDate() )
	{
		if( !cohomologybasis_->MakeUpToDate() )
		{
			return false;
		}
	}
	cohomologybasis_->Simplify();

	if( !FindEdgeMeasure() )
	{
		return false;
	}

	setSize(triangulation_->NumberOfTriangles(),triangulation_->NumberOfVertices());

	LaplacianMatrix laplacian(triangulation_,edge_measure_);

	std::vector<double> minDeltaOmega(triangulation_->NumberOfVertices(),0.0);
	std::vector<double> coordinate(triangulation_->NumberOfVertices(),0.0);
	std::vector<boost::array<double,3> > minomega(triangulation_->NumberOfTriangles());
	for(int i=0;i<2;i++)
	{
		for(int j=0,end=minomega.size();j<end;j++)
		{
			for(int k=0;k<3;k++)
			{
				minomega[j][k] = - (double)(cohomologybasis_->getOmega(j,k,i));
			}
		}
		Coderivative(minomega,minDeltaOmega);

		// it can happen that deltaOmega is identically zero; in that case x[i]=0 is the correct solution
		bool iszero = true;
		for(int j=0;j<triangulation_->NumberOfVertices();j++)
		{
			if( fabs(minDeltaOmega[j]) > 0.01 )
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
			LoadInitialCoordinates(coordinate,i,triangulation_->getVertex(0));
			if( !laplacian.ConjugateGradientSolve(minDeltaOmega,coordinate,accuracy_,maxiterations_) )
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
				setForm(edge,i,(double)(cohomologybasis_->getOmega(edge)[i]) - coordinate[edge->getNext()->getOpposite()->getId()] + coordinate[edge->getPrevious()->getOpposite()->getId()]); 
			}
		}

		for(int j=0;j<triangulation_->NumberOfVertices();j++)
		{
			setCoordinate(j,i,properfmod(coordinate[j],1.0));
		}
	}

	SetUpToDate();
	return true;
}

std::pair< double, double > Embedding::CalculateModuli()
{
	if( !IsUpToDate() )
	{
		FindEmbedding();
	}

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
					inproducts[i][j] += 0.5*getForm(k,l)[i]*getForm(k,l)[j]*edge_measure_[k][l];
				}
			}
		}
	}
	return std::pair<double,double>( -inproducts[0][1] / inproducts[1][1], sqrt(inproducts[0][0] * inproducts[1][1] - inproducts[0][1] * inproducts[1][0])/inproducts[1][1] );
}

bool Embedding::MakeUpToDate()
{
	if( IsUpToDate() )
	{
		return true;
	}

	return FindEmbedding();
}