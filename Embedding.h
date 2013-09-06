#pragma once

#include <vector>

#include "triangulation.h"
#include "Decoration.h"
#include "Vertex.h"
#include "Edge.h"
#include "utilities.h"


class Embedding : public Decoration
{
public:
	Embedding() {}
	~Embedding() {}

	void Initialize() {}

	void UpdateAfterFlipMove(const Edge * const edge) 
	{
		setFormToMinusAdjacent(edge);
		setFormToMinusAdjacent(edge->getPrevious()->getAdjacent()->getNext());
		Vector2D newform = NegateVector2D(AddVectors2D(getForm(edge),getForm(edge->getNext())));
		setForm(edge->getPrevious(),newform);
		setForm(edge->getPrevious()->getAdjacent(),NegateVector2D(newform));
	}

	virtual bool FindEmbedding() = 0;

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
private:
	std::vector<Vector2D> coordinate_;
	std::vector<boost::array<Vector2D,3> > form_;
};
