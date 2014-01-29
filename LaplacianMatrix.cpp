#include <map>

#include "LaplacianMatrix.h"

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


DomainLaplacianMatrix::DomainLaplacianMatrix(const Triangulation * const triangulation, const std::vector<const Vertex*> & vertices, const std::vector<const Vertex*> & fixed) 
	: Matrix(vertices.size()),
	vertices_(vertices)
{
	boost::array<double,3> unitmeasure = {1.0,1.0,1.0};
	std::vector<boost::array<double,3> > edgemeasure(triangulation->NumberOfTriangles(),unitmeasure);
	Initialize(triangulation,edgemeasure,vertices,fixed);
}

DomainLaplacianMatrix::DomainLaplacianMatrix(const Triangulation * const triangulation, const std::vector<boost::array<double,3> > & edge_measure, const std::vector<const Vertex*> & vertices, const std::vector<const Vertex*> & fixed) 
	: Matrix(vertices.size() ),
	vertices_(vertices)
{
	Initialize(triangulation,edge_measure,vertices,fixed);
}

void DomainLaplacianMatrix::Initialize(const Triangulation * const triangulation, const std::vector<boost::array<double,3> > & edge_measure, const std::vector<const Vertex*> & vertices, const std::vector<const Vertex*> & fixed) 
{
	int NotInSet = vertices.size()+1;
	std::vector<int> type(triangulation->NumberOfVertices(),NotInSet);
	for(int i=0,endi=vertices.size();i<endi;i++)
	{
		type[vertices[i]->getId()] = i;
	}
	for(int i=0,endi=fixed.size();i<endi;i++)
	{
		type[fixed[i]->getId()] = -i-1;
	}

	laplacianRules_.resize(vertices.size());
	targetRules_.resize(vertices.size());
	for(int i=0;i<triangulation->NumberOfTriangles();i++)
	{
		Triangle * triangle = triangulation->getTriangle(i);
		for(int j=0;j<3;j++)
		{
			Edge * edge = triangle->getEdge(j);
			int start = edge->getNext()->getOpposite()->getId();
			int end = edge->getPrevious()->getOpposite()->getId();

			if( type[start] >= 0 && type[start] != NotInSet )
			{
				std::pair<std::map<int,double>::iterator,bool> insertion;
				insertion = laplacianRules_[type[start]].insert(std::pair<int,double>(type[start],edge_measure[i][j]));
				if( !insertion.second )
				{
					insertion.first->second += edge_measure[i][j];
				}
				if( type[end] >= 0 && type[end] != NotInSet )
				{
					insertion = laplacianRules_[type[start]].insert(std::pair<int,double>(type[end],-edge_measure[i][j]));
					if( !insertion.second )
					{
						insertion.first->second -= edge_measure[i][j];
					}
				} else
				{
					insertion = targetRules_[type[start]].insert(std::pair<int,double>(fixed[-type[end]-1]->getId(),edge_measure[i][j]));
					if( !insertion.second )
					{
						insertion.first->second += edge_measure[i][j];
					}
				}
			}
		}
	}
}

void DomainLaplacianMatrix::MultiplyVector(const std::vector<double> & from, std::vector<double> & to) const
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

bool DomainLaplacianMatrix::FindHarmonic(std::vector<double> & x, double eps, int maxIterations )
{
	std::vector<double> target(targetRules_.size());
	GetTarget(x,target);
	std::vector<double> xsub(vertices_.size());
	for(int i=0,endi=vertices_.size();i<endi;i++)
	{
		xsub[i] = x[vertices_[i]->getId()];
	}
	return ConjugateGradientSolve(target,xsub,eps,maxIterations);
}

void DomainLaplacianMatrix::GetTarget(const std::vector<double> & x, std::vector<double> & target)
{
	for(int i=0;i<static_cast<int>(targetRules_.size());i++)
	{
		target[i] = 0.0;
		for(std::map<int,double>::const_iterator it = targetRules_[i].begin(); it != targetRules_[i].end(); it++)
		{
			target[i] += it->second * x[it->first];
		}
	}
}