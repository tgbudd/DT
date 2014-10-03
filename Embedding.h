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
#include "LaplacianMatrix.h"

class Embedding : public Decoration
{
public:
	Embedding(const Triangulation * const triangulation, CohomologyBasis * const cohomologybasis);

	~Embedding() {}

	void Initialize() {}
	void UpdateAfterFlipMove(const Edge * const edge); 
	void UpdateAfterCutMove(const boost::array< Edge *, 2> & edges);
	bool FindEmbedding();
	virtual bool FindEdgeMeasure() = 0;
	std::pair< double, double > CalculateModuli();

	bool MakeUpToDate();

	virtual bool GetRadii(std::vector<double> & radii);

	const Vector2D & getCoordinate(const Vertex * const & vertex) const
	{
		return coordinate_[vertex->getId()];
	}
	const Vector2D & getCoordinate(int i) const
	{
		return coordinate_[i];
	}
	void setCoordinate(const Vertex * const & vertex, Vector2D coor) 
	{
		coordinate_[vertex->getId()] = coor;
	}
	void setCoordinate(const Vertex * const & vertex, int i, double coor) 
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
			/*for(int i=0;i<NumberOfVertices;i++)
			{
				coordinate_[i][0] = triangulation_->RandomReal(0.0,1.0);
				coordinate_[i][1] = triangulation_->RandomReal(0.0,1.0);
			}*/
		}
		if( static_cast<int>(form_.size()) != NumberOfTriangles )
		{
			form_.resize(NumberOfTriangles );
			/*for(int i=0;i<NumberOfVertices;i++)
			{
				for(int j=0;j<3;j++)
				{
					form_[i][j][0] = triangulation_->RandomReal(0.0,1.0);
					form_[i][j][1] = triangulation_->RandomReal(0.0,1.0);
				}
			}*/
		}
		if( static_cast<int>(edge_measure_.size()) != NumberOfTriangles )
		{
			edge_measure_.resize(NumberOfTriangles);
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
	void setWorkWithHarmonicForms(bool value)
	{
		work_with_harmonic_forms_ = value;
	}
	void SetAccuracy(double accuracy) 
	{
		accuracy_=accuracy;
	}
	void SaveEmbedding(std::string filename);
	Vector2D GetCentroid(const Triangle * triangle);
private:
	void Coderivative( const std::vector<boost::array<double,3> > & oneform, std::vector<double> & result );
	void LoadInitialCoordinates( std::vector<double> & coordinates, int i, Vertex * startVertex ) const;
	bool CheckClosedness();


	const Triangulation * const triangulation_;
	CohomologyBasis * const cohomologybasis_;
	std::vector<Vector2D> coordinate_;
	std::vector<boost::array<Vector2D,3> > form_;
	std::vector<boost::array<double,3> > edge_measure_;
	double accuracy_;
	int maxiterations_;
	bool work_with_harmonic_forms_;
};

#endif
