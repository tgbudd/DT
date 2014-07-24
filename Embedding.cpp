#include <queue>
#include <vector>
#include <map>
#include <fstream>

#include "Embedding.h"

Embedding::Embedding(const Triangulation * const triangulation, CohomologyBasis * const cohomologybasis) 
	: Decoration(triangulation), triangulation_(triangulation), cohomologybasis_(cohomologybasis) 
{
	work_with_harmonic_forms_ = false;
	accuracy_ = 1.0e-6;
	maxiterations_ = 10000;
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


void Embedding::UpdateAfterCutMove(const boost::array< Edge *, 2> & edges)
{
	Vector2D integral = SubtractVectors2D(getForm(edges[0]),getForm(edges[1]));
	Edge * edge = edges[1];
	while( edge != edges[0] )
	{
		setForm(edge,AddVectors2D(getForm(edge),integral));
		setForm(edge->getNext(),AddVectors2D(getForm(edge->getNext()),NegateVector2D(integral)));
		edge = edge->getNext()->getAdjacent();
	}
	setCoordinate(edges[0]->getNext()->getOpposite(),getCoordinate(edges[1]->getNext()->getOpposite()));
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
	//cohomologybasis_->Simplify();

	setSize(triangulation_->NumberOfTriangles(),triangulation_->NumberOfVertices());

	if( !FindEdgeMeasure() )
	{
		return false;
	}


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
				if( work_with_harmonic_forms_ )
				{
					minomega[j][k] = - getForm(j,k)[i];
				} else
				{
					minomega[j][k] = - (double)(cohomologybasis_->getOmega(j,k,i));
				}
			}
		}
		Coderivative(minomega,minDeltaOmega);

		// it can happen that deltaOmega is identically zero; in that case x[i]=0 is the correct solution
		bool iszero = true;
		for(int j=0;j<triangulation_->NumberOfVertices();j++)
		{
			if( fabs(minDeltaOmega[j]) > 0.00001 )
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
			if( work_with_harmonic_forms_ )
			{
				std::fill(coordinate.begin(),coordinate.end(),0.0);
			} else
			{
				LoadInitialCoordinates(coordinate,i,triangulation_->getVertex(0));
			}
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
				setForm(edge,i,-minomega[j][k] - coordinate[edge->getNext()->getOpposite()->getId()] + coordinate[edge->getPrevious()->getOpposite()->getId()]); 
			}
		}

		for(int j=0;j<triangulation_->NumberOfVertices();j++)
		{
			if( work_with_harmonic_forms_ )
			{
				setCoordinate(j,i,properfmod(getCoordinate(j)[i]+coordinate[j],1.0));
			} else
			{
				setCoordinate(j,i,properfmod(coordinate[j],1.0));
			}
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

bool Embedding::CheckClosedness()
{
	for(int i=0,endi=form_.size();i<endi;i++)
	{
		Vector2D tot = {0.0,0.0};
		for(int j=0;j<3;j++)
		{
			tot = AddVectors2D(tot,form_[i][j]);
		}
		if( std::fabs(tot[0]) > 1e-6 || std::fabs(tot[1]) > 1e-6 )
		{
			return false;
		}
	}
	return true;
}

bool Embedding::GetRadii(std::vector<double> & radii)
{
	return false;
}

Vector2D Embedding::GetCentroid(const Triangle * triangle)
{
	Vector2D x = getCoordinate(triangle->getEdge(0)->getOpposite());
	x = AddVectors2D(x, AddScaledVectors2D(1/3.0,getForm(triangle->getEdge(1)),-1/3.0,getForm(triangle->getEdge(2))));
	return x;
}

void Embedding::SaveEmbedding(std::string filename)
{
	std::ofstream file(filename.c_str());

	std::pair<double,double> tau = CalculateModuli();
	file << std::fixed << "{ modulus -> " << tau.first << " + I " << tau.second << ", points -> {";
	for(int i=0,endi=triangulation_->NumberOfTriangles();i<endi;i++)
	{
		Triangle * triangle = triangulation_->getTriangle(i);
		Vector2D x = coordinate_[triangle->getEdge(0)->getOpposite()->getId()];
		x = AddVectors2D(x,ScaleVector2D(getForm(i,2),1/3.0));
		x = AddVectors2D(x,ScaleVector2D(getForm(i,1),-1/3.0));
		x[0] = properfmod(x[0],1.0);
		x[1] = properfmod(x[1],1.0);
		file << (i>0?",":"") << "{" << x[0] << ", " << x[1] << "}";
	}
	
	/*for(int i=0,endi=coordinate_.size();i<endi;i++)
	{
		file << (i>0?",":"") << "{" << coordinate_[i][0] << ", " << coordinate_[i][1] << "}";
	}*/
	file << "}}\n";
}