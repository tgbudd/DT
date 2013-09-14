#pragma once

#include <vector>
#include <map>

#include "boost/array.hpp"

#include "triangulation.h"
#include "Decoration.h"
#include "CohomologyBasis.h"
#include "ConjugateGradient.h"
#include "Vertex.h"
#include "Edge.h"
#include "utilities.h"


class LaplacianMatrix : public Matrix
{
public:
	LaplacianMatrix (const Triangulation * const triangulation, const std::vector<boost::array<double,3> > & edge_measure);
	void MultiplyVector(const std::vector<double> & from, std::vector<double> & to) const;
private:
	std::vector<std::map<int,double> > laplacianRules_;
};

class Embedding : public Decoration
{
public:
	Embedding(const Triangulation * const triangulation, const CohomologyBasis * const cohomologybasis) : triangulation_(triangulation), cohomologybasis_(cohomologybasis) {
		accuracy_ = 1.0e-6;
		maxiterations_ = 2000;
		edge_measure_.resize(triangulation->NumberOfTriangles());
	}
	~Embedding() {}

	void Initialize() {}
	void UpdateAfterFlipMove(const Edge * const edge); 
	bool FindEmbedding();
	virtual bool FindEdgeMeasure() = 0;
	std::pair< double, double > CalculateModuli() const;


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
		if( coordinate_.size() != NumberOfVertices )
		{
			coordinate_.resize(NumberOfVertices);
		}
		if( form_.size() != NumberOfTriangles )
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
private:
	void Coderivative( const std::vector<boost::array<double,3> > & oneform, std::vector<double> & result );
	void LoadInitialCoordinates( std::vector<double> & coordinates, int i, Vertex * startVertex ) const;

	const Triangulation * const triangulation_;
	const CohomologyBasis * const cohomologybasis_;
	std::vector<Vector2D> coordinate_;
	std::vector<boost::array<Vector2D,3> > form_;
	std::vector<boost::array<double,3> > edge_measure_;
	double accuracy_;
	int maxiterations_;
};
