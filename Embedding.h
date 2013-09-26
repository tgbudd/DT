#ifndef DT_EMBEDDING_H
#define DT_EMBEDDING_H

#include <vector>
#include <map>

#include "boost/array.hpp"

#include "Triangulation.h"
#include "Decoration.h"
#include "CohomologyBasis.h"
#include "DualCohomologyBasis.h"
#include "LinearAlgebra.h"
#include "Vertex.h"
#include "Edge.h"
#include "utilities.h"


class LaplacianMatrix : public linearalgebra::Matrix
{
public:
	LaplacianMatrix (const Triangulation * const triangulation);
	LaplacianMatrix (const Triangulation * const triangulation, const std::vector<boost::array<double,3> > & edge_measure);
	void MultiplyVector(const std::vector<double> & from, std::vector<double> & to) const;
private:
	void Initialize(const Triangulation * const triangulation, const std::vector<boost::array<double,3> > & edge_measure);
	std::vector<std::map<int,double> > laplacianRules_;
};

class Embedding : public Decoration
{
public:
	Embedding(const Triangulation * const triangulation, CohomologyBasis * const cohomologybasis);

	~Embedding() {}

	void Initialize() {}
	void UpdateAfterFlipMove(const Edge * const edge); 
	bool FindEmbedding();
	virtual bool FindEdgeMeasure() = 0;
	std::pair< double, double > CalculateModuli();

	bool MakeUpToDate();

	const Vector2D & getCoordinate(Vertex * const & vertex) const
	{
		return coordinate_[vertex->getId()];
	}
	const Vector2D & getCoordinate(int i) const
	{
		return coordinate_[i];
	}
	void setCoordinate(Vertex * const & vertex, Vector2D coor) 
	{
		coordinate_[vertex->getId()] = coor;
	}
	void setCoordinate(Vertex * const & vertex, int i, double coor) 
	{
		coordinate_[vertex->getId()][i] = coor;
	}
	void setCoordinate(int vertex, int i, double coor) 
	{
		coordinate_[vertex][i] = coor;
	}
	const Vector2D & getForm(const Edge * const & edge) const
	{
		return form_[edge->getParent()->getId()][edge->getId()];
	}
	const Vector2D & getForm(int triangle, int edge) const
	{
		return form_[triangle][edge];
	}
	void setForm(const Edge * const & edge, const Vector2D & coor )
	{
		 form_[edge->getParent()->getId()][edge->getId()] = coor;
	}
	void setForm(const Edge * const & edge, int i, double coor )
	{
		 form_[edge->getParent()->getId()][edge->getId()][i] = coor;
	}
	void setForm(int triangle, int edge, int i, double coor )
	{
		 form_[triangle][edge][i] = coor;
	}
	void setFormToMinusAdjacent(const Edge * const & edge)
	{
		Edge * adj = edge->getAdjacent();
		form_[edge->getParent()->getId()][edge->getId()] = NegateVector2D(form_[adj->getParent()->getId()][adj->getId()] );
	}
	void setSize(int NumberOfTriangles,int NumberOfVertices)
	{
		if( static_cast<int>(coordinate_.size()) != NumberOfVertices )
		{
			coordinate_.resize(NumberOfVertices);
		}
		if( static_cast<int>(form_.size()) != NumberOfTriangles )
		{
			form_.resize(NumberOfTriangles );
		}
	}
	void setMaxIterations(int iterations)
	{
		maxiterations_ = iterations;
	}

	void setEdgeMeasure(int triangle, int edge, double measure)
	{
		edge_measure_[triangle][edge] = measure;
	}
	int getSize() const
	{
		return static_cast<int>(form_.size());
	}
private:
	void Coderivative( const std::vector<boost::array<double,3> > & oneform, std::vector<double> & result );
	void LoadInitialCoordinates( std::vector<double> & coordinates, int i, Vertex * startVertex ) const;

	const Triangulation * const triangulation_;
	CohomologyBasis * const cohomologybasis_;
	std::vector<Vector2D> coordinate_;
	std::vector<boost::array<Vector2D,3> > form_;
	std::vector<boost::array<double,3> > edge_measure_;
	double accuracy_;
	int maxiterations_;
};

#endif
